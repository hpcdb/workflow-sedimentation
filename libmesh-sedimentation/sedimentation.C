/*
 * File:   sedimentation.cpp
 * Author: camata
 *
 * Created on December 9, 2014, 3:05 PM
 */
// C++ include files that we need
#include <iostream>
#include <chrono>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <string>
#include <unistd.h>
#include <vector>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/vtk_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/perf_log.h"
#include "libmesh/boundary_info.h"
#include "libmesh/utility.h"

#include "libmesh/getpot.h"
#include "libmesh/mesh_refinement.h"

#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

#define XDMF_

#include "xdmf.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;
using namespace std;

#include "sedimentation_flow.h"
#include "sedimentation_transport.h"
#include "sedimentation_deposition.h"
#include "mesh_moviment.h"
#include "provenance.h"
#include "performance.h"
#include "FEAdaptor.h"

const int textArraySize = 64;
const int jsonArraySize = 512;

double ramp(double t) {
    double x[3];
    double y[3];

    x[0] = 0.0;
    y[0] = 0.01;
    x[1] = 1.0;
    y[1] = 1.0;
    x[2] = 1.0E04;
    y[2] = 1.0;

    double L0 = ((t - x[1])*(t - x[2])) / ((x[0] - x[1])*(x[0] - x[2]));
    double L1 = ((t - x[0])*(t - x[2])) / ((x[1] - x[0])*(x[1] - x[2]));
    double L2 = ((t - x[0])*(t - x[1])) / ((x[2] - x[0])*(x[2] - x[1]));

    return (y[0] * L0 + y[1] * L1 + y[2] * L2);
}

void WriteRestartFile(EquationSystems &es, int t_step, std::string rname) {
    const std::string mesh_restart = rname + "_mesh_restart.xda";
    const std::string solution_restart = rname + "_solution_restart.xda";

    Real time = es.parameters.get<Real> ("time");
    Real dt = es.parameters.get<Real> ("dt");
}

bool is_file_exist(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

//#define MESH_MOVIMENT

// We can now begin the main program.  Note that this
// example will fail if you are using complex numbers
// since it was designed to be run only with real numbers.

int main(int argc, char** argv) {
    Performance solverPerf;

#ifdef PERFORMANCE
    if (libMesh::global_processor_id() == 0) {
        solverPerf.start();
    }
#endif

    int simulationID = 1;

    // Initialize libMesh.
    LibMeshInit init(argc, argv);

    PerfLog perf_log("Sedimentation Solver");

    // This example requires Adaptive Mesh Refinement support - although
    // it only refines uniformly, the refinement code used is the same
    // underneath
    GetPot infile("sedimentation.in");

    int init_tstep = 0;

    int dim = infile("dim", 2);
    int ncellx = infile("ncellx", 10);
    int ncelly = infile("ncelly", 10);
    int ncellz = infile("ncellz", 10);
    double xmin = infile("xmin", 0.0);
    double ymin = infile("ymin", 0.0);
    double zmin = infile("zmin", 0.0);
    double xmax = infile("xmax", 1.0);
    double ymax = infile("ymax", 1.0);
    double zmax = infile("zmax", 1.0);
    int ref_interval = infile("r_interval", 1);

    Provenance prov;

#ifdef PROV
    // Mesh Generation
    prov.inputMeshGeneration(simulationID, dim, ncellx, ncelly, ncellz, xmin, ymin, zmin, xmax, ymax, zmax, ref_interval);
#endif
    // Create a mesh object, with dimension to be overridden later,
    // distributed across the default MPI communicator.
    Mesh mesh(init.comm());

    double r_fraction = infile("r_fraction", 0.70);
    double c_fraction = infile("c_fraction", 0.10);
    double max_h_level = infile("max_h_level", 1);
    const unsigned int hlevels = infile("hlevels", 0);

    MeshRefinement refinement(mesh);

    refinement.refine_fraction() = r_fraction;
    refinement.coarsen_fraction() = c_fraction;
    refinement.max_h_level() = max_h_level;

    bool first_step_refinement = true;

    // Create an equation systems object.
    EquationSystems equation_systems(mesh);

    SedimentationFlow sediment_flow(equation_systems);
    SedimentationTransport sediment_transport(equation_systems);
    SedimentationDeposition sediment_deposition(equation_systems);

#ifdef MESH_MOVIMENT
    MeshMoviment moving_mesh(equation_systems);
#endif

    // INPUT: Physical parameters
    Real Reynolds = infile("Reynolds", 1.0E03);
    Real Gr = infile("Grashof", 1.0E05);
    Real Sc = infile("Sc", 0.71);
    Real Us = infile("Us", 0.0);
    Real xlock = infile("xlock", 1.0);
    Real alfa = infile("alfa", 1.0);
    Real theta = infile("theta", 0.5);
    Real ex = infile("ex", 0.0);
    Real ey = infile("ey", 0.0);
    Real ez = infile("ez", -1.0);
    Real c_factor = infile("c_factor", 1.0);
    Real fopc = infile("fopc", 0.2);

    if (Reynolds == 0.0)
        Reynolds = std::sqrt(Gr);

    Real Diffusivity = 1.0 / (Sc * Reynolds);

    cout << "Parameters: " << endl;
    cout << "  Reynolds: " << Reynolds << endl;
    cout << "  Grashof : " << Gr << endl;
    cout << "  Sc      : " << Sc << endl;

    equation_systems.parameters.set<Real> ("Reynolds") = Reynolds;
    equation_systems.parameters.set<Real> ("Diffusivity") = Diffusivity;
    equation_systems.parameters.set<Real> ("Us") = Us;

    /////// 
    equation_systems.parameters.set<Real> ("xlock") = xlock;
    equation_systems.parameters.set<Real> ("theta") = theta;
    equation_systems.parameters.set<Real> ("ex") = ex;
    equation_systems.parameters.set<Real> ("ey") = ey;
    equation_systems.parameters.set<Real> ("ez") = ez;
    equation_systems.parameters.set<Real> ("c_factor") = c_factor;
    equation_systems.parameters.set<Real> ("fopc") = fopc;

    ////

    // LOOP 
    Real init_time = 0.0;

    // INPUT: TIME INTEGRATION
    Real dt = infile("deltat", 0.005);
    Real tmax = infile("tmax", 0.1);
    const unsigned int n_time_steps = infile("n_time_steps", 10);

    //INPUT: NONLINEAR SOLVER
    const unsigned int n_nonlinear_steps = infile("n_nonlinear_steps", 10);
    double nonlinear_tolerance = infile("nonlinear_tolerance", 1.0E-03);
    int max_linear_iters = infile("max_linear_iterations", 2000);
    int max_r_steps = infile("max_r_steps", 1);

    const unsigned int write_interval = infile("write_interval", 10);
    std::string rname = infile("output", "out");

#ifdef XDMF_
    XDMFWriter xdmf_writer(mesh);
#endif

    if (!is_file_exist("restart.in") || true) {

        const string mesh_file = infile("mesh_file", "0");

        std::cout << "Opening file: " << mesh_file << endl;

        if (dim == 2) {
            MeshTools::Generation::build_square(mesh, ncellx, ncelly,
                    xmin, xmax,
                    ymin, ymax, QUAD4);
        } else {
            MeshTools::Generation::build_cube(mesh, ncellx, ncelly, ncellz,
                    xmin, xmax,
                    ymin, ymax,
                    zmin, zmax, HEX8);
        }

        refinement.uniformly_refine(hlevels);

        sediment_flow.setup();
        sediment_transport.setup();
        sediment_deposition.setup();

#ifdef MESH_MOVIMENT
        moving_mesh.setup();
#endif

#ifdef PROV
        // Mesh Refinement
        prov.outputMeshGeneration(simulationID, r_fraction, c_fraction, max_h_level, hlevels);
#endif

        // Initialize the data structures for the equation system.
        equation_systems.init();

#ifdef PROV
        // Generate solver parameters
        prov.outputCreateEquationSystems(simulationID, Reynolds, Gr, Sc, Us, Diffusivity, xlock, alfa, theta, ex, ey, ez, c_factor);
#endif

    } else {
        GetPot restart("restart.in");
        const string mesh_restart = restart("mesh_restart", "0");
        const string solution_restart = restart("solution_restart", "0");
        init_time = restart("init_time", 0.0);
        dt = restart("deltat", 0.0);

#ifdef LIBMESH_HAVE_HDF5
        int xdmf_file_id = restart("xdmf_file_id", 0);
        xdmf_writer.set_file_id(xdmf_file_id);
#endif

        mesh.read(mesh_restart);

#ifdef PROV
        // Mesh Refinement
        prov.outputMeshGeneration(simulationID, r_fraction, c_fraction, max_h_level, hlevels);
#endif

        equation_systems.read(solution_restart, READ);

        // Get a reference to the Convection-Diffusion system object.
        TransientLinearImplicitSystem & transport_system =
                equation_systems.get_system<TransientLinearImplicitSystem> ("sediment");

        //transport_system.add_vector("volume");

        transport_system.update();

        // Get a reference to the Convection-Diffusion system object.
        TransientLinearImplicitSystem & flow_system =
                equation_systems.get_system<TransientLinearImplicitSystem> ("flow");

        flow_system.update();

        ExplicitSystem & deposition_system = equation_systems.get_system<ExplicitSystem>("deposition");
        deposition_system.add_vector("deposition_rate");

        deposition_system.update();

#ifdef PROV
        // Generate solver parameters
        prov.outputCreateEquationSystems(simulationID, Reynolds, Gr, Sc, Us, Diffusivity, xlock, alfa, theta, ex, ey, ez, c_factor);
#endif

    }

    // Print information about the mesh to the screen.
    mesh.print_info();

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & transport_system =
            equation_systems.get_system<TransientLinearImplicitSystem> ("sediment");

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & flow_system =
            equation_systems.get_system<TransientLinearImplicitSystem> ("flow");

#ifdef MESH_MOVIMENT
    LinearImplicitSystem & mesh_system = equation_systems.get_system<LinearImplicitSystem> ("mesh_moviment");
#endif
    // Prints information about the system to the screen.
    equation_systems.print_info();

    double time = init_time;
    transport_system.time = time;
    flow_system.time = time;

    unsigned int t_step = 0;
    unsigned int n_linear_iterations_flow = 0;
    unsigned int n_nonlinear_iterations_flow = 0;
    unsigned int n_nonlinear_iterations_transport = 0;
    unsigned int n_linear_iterations_transport = 0;
    bool redo_nl;

    string* current_files;

#ifdef LIBMESH_HAVE_HDF5
    current_files = xdmf_writer.write_time_step(equation_systems, time);
    cout << "[WRITE] " + current_files[0] + " - " + current_files[1] << endl;
#else
    int exodus_step = 0;
    {
        std::ostringstream out;
        out << rname << "_" << std::setw(5) << std::setfill('0') << exodus_step << ".e";
        ExodusII_IO(mesh).write_equation_systems(out.str(), equation_systems);
        exodus_step++;
    }
    //std::string exodus_filename = "output.e";
#endif

#ifdef PROV
    // Generate loop iterations
    prov.outputGetMaximumIterations(simulationID, dt, tmax, n_time_steps, n_nonlinear_steps, nonlinear_tolerance, max_linear_iters, max_r_steps, write_interval, current_files[1]);
#endif

    Performance perf;

    if (dim == 2) {
        // 2D analysis
        char firstFilename[textArraySize];
        sprintf(firstFilename, "init_ext_plane_%d.csv", t_step);
        char finalFilename[textArraySize];
        sprintf(finalFilename, "ext_plane_%d.csv", t_step);
#ifdef PROV
        // Mesh Writer
        prov.inputInitDataExtraction(simulationID, "initdataextraction");
#endif

#ifdef PERFORMANCE
        if (libMesh::global_processor_id() == 0) {
            perf.start();
        }
#endif

#ifdef USE_CATALYST
        FEAdaptor::Initialize(argc, argv);
        FEAdaptor::CoProcess(argc, argv, equation_systems, 0.0, t_step, false, false);
        if (libMesh::global_processor_id() == 0) {
            char commandLine[jsonArraySize];
            sprintf(commandLine, "python clean-csv.py %s %s;rm %s", firstFilename, finalFilename, firstFilename);
            system(strdup(commandLine));
        }
#endif  

#ifdef PERFORMANCE
        if (libMesh::global_processor_id() == 0) {
            perf.end();
            double elapsedTime = perf.elapsedTime();
            char buffer[textArraySize];
            sprintf(buffer, "Data Extraction Cost: %.5f", elapsedTime);
            cout << buffer << endl;
            prov.storeDataExtractionCost(elapsedTime);
        }
#endif  

#ifdef PROV
        // Mesh Writer
        prov.outputInitDataExtraction(simulationID, "initdataextraction", "oinitdataextraction", 0, current_files[1], finalFilename, dim, "irde");
#endif
    } else if (dim == 3) {
        // 3D analysis
        for (int ik = 0; ik <= 3; ik++) {
            char firstFilename[textArraySize];
            sprintf(firstFilename, "init_ext_line_%d_%d.csv", ik, t_step);
            char finalFilename[textArraySize];
            sprintf(finalFilename, "ext_line_%d_%d.csv", ik, t_step);

#ifdef PROV
            // Mesh Writer
            char argument1[textArraySize];
            sprintf(argument1, "iline%dextraction", ik);
            prov.inputInitDataExtraction(simulationID, argument1);
#endif

#ifdef PERFORMANCE
            if (libMesh::global_processor_id() == 0) {
                perf.start();
            }
#endif

#ifdef USE_CATALYST
            if (ik == 0) {
                FEAdaptor::Initialize(argc, argv);
                FEAdaptor::CoProcess(argc, argv, equation_systems, 0.0, t_step, false, false);
            }
            if (libMesh::global_processor_id() == 0) {
                char commandLine[jsonArraySize];
                sprintf(commandLine, "python clean-csv.py %s %s;rm %s", firstFilename, finalFilename, firstFilename);
                system(strdup(commandLine));
            }
#endif  

#ifdef PERFORMANCE
            if (libMesh::global_processor_id() == 0) {
                perf.end();
                double elapsedTime = perf.elapsedTime();
                char buffer[textArraySize];
                sprintf(buffer, "Data Extraction Cost: %.5f", elapsedTime);
                cout << buffer << endl;
                prov.storeDataExtractionCost(elapsedTime);
            }
#endif  

#ifdef PROV
            // Mesh Writer
            sprintf(argument1, "iline%dextraction", ik);
            char argument2[textArraySize];
            sprintf(argument2, "oline%diextraction", ik);
            char* extractorName = (char*) malloc(jsonArraySize);
            sprintf(extractorName, "iline%d", ik);
            prov.outputInitDataExtraction(simulationID, argument1, argument2, 0, current_files[1], finalFilename, dim, extractorName);
            free(extractorName);
#endif
        }
    }


    // STEP LOOP
    // Loop in time steps
    int taskID = 0;
    int numberIterationsFluid = 0;
    int numberIterationsSediments = 0;
    int numberIterationsMeshRefinements = 0;
    int numberOfWrites = 0;
    vector<string> meshDependencies;

    for (t_step = init_tstep; (t_step < n_time_steps)&&(time < tmax); t_step++) {
        taskID++;
        if (is_file_exist("abort.run")) break;

        // Increment the time counter, set the time and the
        // time step size as parameters in the EquationSystem.
        time += dt;
        flow_system.time += dt;
        transport_system.time += dt;

        Real tmp = Reynolds;
        equation_systems.parameters.set<Real> ("Reynolds") = tmp;
        equation_systems.parameters.set<Real> ("Diffusivity") = 1.0 / (Sc * tmp);
        equation_systems.parameters.set<Real> ("time") = time;
        equation_systems.parameters.set<Real> ("dt") = dt;

        // A pretty update message
        std::cout << std::setw(55)
                << std::setfill('=')
                << "\n";
        std::cout << " Time step ";
        {
            std::ostringstream out;

            out << std::setw(2)
                    << std::right
                    << t_step
                    << "  simulation time = "
                    << std::fixed
                    << std::setw(6)
                    << std::setprecision(3)
                    << std::setfill('0')
                    << std::left
                    << transport_system.time
                    << "...";

            std::cout << out.str() << std::endl;
        }
        std::cout << std::setw(55)
                << std::setfill('=')
                << "\n";

        // At the beginning of each solve, reset the linear solver tolerance
        // to a "reasonable" starting value.
        const Real initial_linear_solver_tol = 1.e-6;
        equation_systems.parameters.set<Real> ("linear solver tolerance") = initial_linear_solver_tol;
        equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = max_linear_iters;

        // AMR/C Loop
        redo_nl = true;
        *flow_system.old_local_solution = *flow_system.current_local_solution;
        *transport_system.old_local_solution = *transport_system.current_local_solution;

        // Loop in linear steps
        for (unsigned int r = 0; r < max_r_steps; r++) {
            if (!redo_nl) break;
            std::cout << " Solving Navier-Stokes equation..." << std::endl;
            {
                std::ostringstream out;

                // We write the file in the ExodusII format.
                out << std::setw(55)
                        << std::setfill('-')
                        << "\n"
                        << std::setfill(' ')
                        << std::setw(5)
                        << std::left
                        << "STEP"
                        << std::setw(5)
                        << "NLI"
                        << std::setw(15)
                        << "|b-AX|"
                        << std::setw(15)
                        << "|du|"
                        << std::setw(15)
                        << "|du|/|u|"
                        << "\n"
                        << std::right
                        << std::setw(55)
                        << std::setfill('-')
                        << "\n";
                std::cout << out.str() << std::flush;
            }

            // determine if we can exit the nonlinear loop.
            UniquePtr<NumericVector<Number> >
                    flow_last_nonlinear_soln(flow_system.solution->clone());

            // Fluid
            // FLOW NONLINEAR LOOP
            for (unsigned int l = 0; l < n_nonlinear_steps; ++l) {
                numberIterationsFluid++;
#ifdef PROV
                // Fluids
                prov.inputSolverSimulationFluid(taskID, simulationID, numberIterationsFluid);
#endif

                // Update the nonlinear solution.
                flow_last_nonlinear_soln->zero();
                flow_last_nonlinear_soln->add(*flow_system.solution);

                // Assemble & solve the linear system.
                flow_system.solve();

                // Compute the difference between this solution and the last
                // nonlinear iterate.
                flow_last_nonlinear_soln->add(-1., *flow_system.solution);

                // Close the vector before computing its norm
                flow_last_nonlinear_soln->close();

                // Compute the l2 norm of the difference
                const Real norm_delta = flow_last_nonlinear_soln->l2_norm();
                const Real u_norm = flow_system.solution->l2_norm();

                // How many iterations were required to solve the linear system?
                const unsigned int n_linear_iterations = flow_system.n_linear_iterations();
                //Total number of linear iterations (so far)
                n_linear_iterations_flow = n_linear_iterations_flow + n_linear_iterations;

                // What was the final residual of the linear system?
                const Real final_linear_residual = flow_system.final_linear_residual();

                {
                    std::ostringstream out;

                    // We write the file in the ExodusII format.
                    out << std::setw(5)
                            << std::left
                            << l
                            << std::setw(5)
                            << n_linear_iterations
                            << std::setw(15)
                            << final_linear_residual
                            << std::setw(15)
                            << norm_delta
                            << std::setw(15)
                            << norm_delta / u_norm
                            << "\n";

                    std::cout << out.str() << std::flush;
                }

                bool converged = (norm_delta < nonlinear_tolerance);

                //Total number of non-linear iterations (so far)
                n_nonlinear_iterations_flow++;

#ifdef PROV
                // Fluids
                prov.outputSolverSimulationFluid(taskID, simulationID, numberIterationsFluid, t_step, transport_system.time, r, l, n_linear_iterations, final_linear_residual, norm_delta, norm_delta / u_norm, converged);
#endif  

                // Terminate the solution iteration if the difference between
                // this nonlinear iterate and the last is sufficiently small, AND
                // if the most recent linear system was solved to a sufficient tolerance.
                if ((norm_delta < nonlinear_tolerance)
                        // && (flow_system.final_linear_residual() < nonlinear_tolerance)
                        ) {
                    std::ostringstream out;
                    // We write the file in the ExodusII format.
                    out << std::setw(55)
                            << std::setfill('-')
                            << "\n";
                    std::cout << out.str()
                            << " Nonlinear solver converged at step "
                            << l
                            << std::endl;
                    break;
                }

                // Otherwise, decrease the linear system tolerance.  For the inexact Newton
                // method, the linear solver tolerance needs to decrease as we get closer to
                // the solution to ensure quadratic convergence.  The new linear solver tolerance
                // is chosen (heuristically) as the square of the previous linear system residual norm.
                //Real flr2 = final_linear_residual*final_linear_residual;
                equation_systems.parameters.set<Real> ("linear solver tolerance") =
                        std::min(Utility::pow<2>(final_linear_residual), initial_linear_solver_tol);

            } // end nonlinear loop

            std::cout << " Solving sedimentation equation..." << std::endl;
            {
                std::ostringstream out;

                // We write the file in the ExodusII format.
                out << std::setw(55)
                        << std::setfill('-')
                        << "\n"
                        << std::setfill(' ')
                        << std::setw(5)
                        << std::left
                        << "STEP"
                        << std::setw(5)
                        << "NLI"
                        << std::setw(15)
                        << "|b-AX|"
                        << std::setw(15)
                        << "|du|"
                        << std::setw(15)
                        << "|du|/|u|"
                        << "\n"
                        << std::right
                        << std::setw(55)
                        << std::setfill('-')
                        << "\n";
                std::cout << out.str() << std::flush;
            }

            equation_systems.parameters.set<Real> ("linear solver tolerance") = initial_linear_solver_tol;

            // determine if we can exit the nonlinear loop.
            UniquePtr<NumericVector<Number> > sed_last_nonlinear_soln(transport_system.solution->clone());

            // Sediments
            // FLOW NON-LINEAR LOOP
            for (unsigned int l = 0; l < n_nonlinear_steps; ++l) {
                numberIterationsSediments++;

#ifdef PROV
                // Fluids
                prov.inputSolverSimulationSediments(taskID, simulationID, numberIterationsSediments);
#endif

                // Update the nonlinear solution.
                sed_last_nonlinear_soln->zero();
                sed_last_nonlinear_soln->add(*transport_system.solution);

                // Assemble & solve the linear system.
                transport_system.solve();

                // Compute the difference between this solution and the last
                // nonlinear iterate.
                sed_last_nonlinear_soln->add(-1., *transport_system.solution);

                // Close the vector before computing its norm
                sed_last_nonlinear_soln->close();

                // Compute the l2 norm of the difference
                const Real norm_delta = sed_last_nonlinear_soln->l2_norm();
                const Real u_norm = transport_system.solution->l2_norm();

                // How many iterations were required to solve the linear system?
                const unsigned int n_linear_iterations = transport_system.n_linear_iterations();
                //Total number of linear iterations (so far)
                n_linear_iterations_transport += n_linear_iterations;

                // What was the final residual of the linear system?
                const Real final_linear_residual = transport_system.final_linear_residual();

                {
                    std::ostringstream out;

                    // We write the file in the ExodusII format.
                    out << std::setw(5)
                            << std::left
                            << l
                            << std::setw(5)
                            << n_linear_iterations
                            << std::setw(15)
                            << final_linear_residual
                            << std::setw(15)
                            << norm_delta
                            << std::setw(15)
                            << norm_delta / u_norm
                            << "\n";

                    std::cout << out.str() << std::flush;
                }

                bool converged = (norm_delta < nonlinear_tolerance) &&
                        (transport_system.final_linear_residual() < nonlinear_tolerance);

                //Total number of non-linear iterations (so far)
                n_nonlinear_iterations_transport++;

#ifdef PROV
                // Sediments
                prov.outputSolverSimulationSediments(taskID, simulationID, numberIterationsSediments, t_step, transport_system.time, r, l, n_linear_iterations, final_linear_residual, norm_delta, norm_delta / u_norm, converged);
#endif

                // Terminate the solution iteration if the difference between
                // this nonlinear iterate and the last is sufficiently small, AND
                // if the most recent linear system was solved to a sufficient tolerance.
                if ((norm_delta < nonlinear_tolerance) &&
                        (transport_system.final_linear_residual() < nonlinear_tolerance)) {
                    // std::ostringstream out;
                    //  // We write the file in the ExodusII format.
                    //  out << std::setw(55)
                    //      << std::setfill('-')
                    //      << "\n";
                    //  std::cout << out.str()
                    //           << " Nonlinear solver converged at step "
                    //           << l
                    //           << std::endl;
                    break;
                }

                // Otherwise, decrease the linear system tolerance.  For the inexact Newton
                // method, the linear solver tolerance needs to decrease as we get closer to
                // the solution to ensure quadratic convergence.  The new linear solver tolerance
                // is chosen (heuristically) as the square of the previous linear system residual norm.
                //Real flr2 = final_linear_residual*final_linear_residual;
                equation_systems.parameters.set<Real> ("linear solver tolerance") =
                        std::min(Utility::pow<2>(final_linear_residual), initial_linear_solver_tol);

            } // end nonlinear loop

            redo_nl = false;

            if (first_step_refinement || (((r + 1) != max_r_steps) && (t_step + 1) % ref_interval == 0)) {
                std::cout << "\n****************** Mesh Refinement ********************  " << std::endl;
                numberIterationsMeshRefinements++;
                int beforeNActiveElem = mesh.n_active_elem();
                std::cout << "Number of elements before AMR step: " << mesh.n_active_elem() << std::endl;
                Real H1norm = transport_system.calculate_norm(*transport_system.solution, SystemNorm(H1));
                ErrorVector error;
                KellyErrorEstimator error_estimator;
                error_estimator.estimate_error(transport_system, error);
                refinement.flag_elements_by_error_fraction(error);
                refinement.refine_and_coarsen_elements();
                equation_systems.reinit();
                redo_nl = true;
                std::cout << "Number of elements after AMR step: " << mesh.n_active_elem() << std::endl;

#ifdef PROV
                // Mesh Refinement
                prov.outputMeshRefinement(taskID, simulationID, numberIterationsMeshRefinements, first_step_refinement, t_step, beforeNActiveElem, mesh.n_active_elem());
#endif

                first_step_refinement = false;
            }

            sediment_deposition.ComputeDeposition();
            sediment_deposition.print();

#ifdef MESH_MOVIMENT
            std::cout << "Solving Mesh moviment..." << std::endl;

            equation_systems.parameters.set<Real> ("linear solver tolerance") = initial_linear_solver_tol;
            mesh_system.solve();

            // How many iterations were required to solve the linear system?
            const unsigned int n_linear_iterations = flow_system.n_linear_iterations();
            std::cout << " Number of iterations: " << n_linear_iterations << std::endl;


            // What was the final residual of the linear system?
            const Real final_linear_residual = flow_system.final_linear_residual();
            std::cout << " final residual " << final_linear_residual << std::endl;

            moving_mesh.updateMesh();
#endif

            // Output every 10 timesteps to file.
            if ((t_step + 1) % write_interval == 0) {
                numberOfWrites++;
                const std::string mesh_restart = rname + "_mesh_restart.xda";
                const std::string solution_restart = rname + "_solution_restart.xda";

                mesh.write(mesh_restart);
                equation_systems.write(solution_restart, WRITE);

                std::ofstream fout;
                fout.open("restart.in");
                fout << "#Restart file: " << std::endl;
                fout << "time = " << time << std::endl;
                fout << "dt   = " << dt << std::endl;
                fout << "init_tstep = " << t_step << std::endl;
                fout << "mesh_restart = " << mesh_restart << std::endl;
                fout << "solution_restart = " << solution_restart << std::endl;

#ifdef LIBMESH_HAVE_HDF5
                fout << "xdmf_file_id = " << xdmf_writer.get_file_id() << std::endl;
#endif

                fout.close();

#ifdef PROV
                // Mesh Writer
                prov.inputMeshWriter(taskID, simulationID, numberOfWrites);
#endif

#ifdef LIBMESH_HAVE_HDF5
                current_files = xdmf_writer.write_time_step(equation_systems, time);
                cout << "[WRITE] " + current_files[0] + " - " + current_files[1] << endl;
#else
                //ExodusII_IO exo(mesh);
                //exo.append(true);
                //exo.write_time_step (exodus_filename, equation_systems, t_step+1, flow_system.time);

                {
                    std::ostringstream out;
                    out << rname << "_" << std::setw(5) << std::setfill('0') << exodus_step << ".e";
                    ExodusII_IO(mesh).write_equation_systems(out.str(), equation_systems);
                    exodus_step++;
                }
#endif

#ifdef PROV
                prov.outputMeshWriter(taskID, simulationID, numberOfWrites, t_step, current_files[1]);
#endif

                int step = t_step + 1;
                if (dim == 2) {
#ifdef PROV
                    prov.inputDataExtraction(taskID, simulationID, numberOfWrites, "dataextraction");
#endif

#ifdef PERFORMANCE
                    if (libMesh::global_processor_id() == 0) {
                        perf.start();
                    }
#endif

                    char firstFilename[textArraySize];
                    sprintf(firstFilename, "init_ext_plane_%d.csv", step);
                    char finalFilename[textArraySize];
                    sprintf(finalFilename, "ext_plane_%d.csv", step);

#ifdef USE_CATALYST
                    FEAdaptor::CoProcess(argc, argv, equation_systems, transport_system.time, step, false, false);
                    if (libMesh::global_processor_id() == 0) {
                        char commandLine[jsonArraySize];
                        sprintf(commandLine, "python clean-csv.py %s %s;rm %s", firstFilename, finalFilename, firstFilename);
                        system(strdup(commandLine));
                    }
#endif

#ifdef PERFORMANCE
                    if (libMesh::global_processor_id() == 0) {
                        perf.end();
                        double elapsedTime = perf.elapsedTime();
                        char buffer[textArraySize];
                        sprintf(buffer, "Data Extraction Cost: %.5f", elapsedTime);
                        cout << buffer << endl;
                        prov.storeDataExtractionCost(elapsedTime);
                    }
#endif

#ifdef PROV
                    char* extractorName = (char*) malloc(jsonArraySize);
                    sprintf(extractorName, "rde%d", numberOfWrites);
                    prov.outputDataExtraction(taskID, simulationID, numberOfWrites, "dataextraction", "odataextraction", step, current_files[1], finalFilename, dim, extractorName);
#endif
                } else if (dim == 3) {
                    // 3D analysis
                    for (int ik = 0; ik <= 3; ik++) {
                        char firstFilename[textArraySize];
                        sprintf(firstFilename, "init_ext_line_%d_%d.csv", ik, step);
                        char finalFilename[textArraySize];
                        sprintf(finalFilename, "ext_line_%d_%d.csv", ik, step);

#ifdef PROV
                        // Mesh Writer
                        char argument1[textArraySize];
                        sprintf(argument1, "line%dextraction", ik);
                        prov.inputDataExtraction(taskID, simulationID, numberOfWrites, argument1);
#endif

#ifdef PERFORMANCE
                        if (libMesh::global_processor_id() == 0) {
                            perf.start();
                        }
#endif

#ifdef USE_CATALYST
                        if (ik == 0) {
                            FEAdaptor::CoProcess(argc, argv, equation_systems, transport_system.time, step, false, false);
                        }
                        if (libMesh::global_processor_id() == 0) {
                            char commandLine[jsonArraySize];
                            sprintf(commandLine, "python clean-csv.py %s %s;rm %s", firstFilename, finalFilename, firstFilename);
                            system(strdup(commandLine));
                        }
#endif  

#ifdef PERFORMANCE
                        if (libMesh::global_processor_id() == 0) {
                            perf.end();
                            double elapsedTime = perf.elapsedTime();
                            char buffer[textArraySize];
                            sprintf(buffer, "Data Extraction Cost: %.5f", elapsedTime);
                            cout << buffer << endl;
                            prov.storeDataExtractionCost(elapsedTime);
                        }
#endif  

#ifdef PROV
                        // Mesh Writer
                        sprintf(argument1, "line%dextraction", ik);
                        char* extractorName = (char*) malloc(jsonArraySize);
                        sprintf(extractorName, "line%d%d", ik, numberOfWrites);
                        char argument2[textArraySize];
                        sprintf(argument2, "oline%dextraction", ik);
                        prov.outputDataExtraction(taskID, simulationID, numberOfWrites, argument1, argument2, 0, current_files[1], finalFilename, dim, extractorName);
                        free(extractorName);
#endif
                    }
                }

                char* depID = (char*) malloc(jsonArraySize);
                sprintf(depID, "%d", taskID);
                meshDependencies.push_back(depID);
                free(depID);
            }
        }
    }

    if ((t_step + 1) % write_interval != 0) {
        numberOfWrites++;
#ifdef PROV
        // Mesh Writer
        prov.inputMeshWriter(taskID, simulationID, numberOfWrites);
#endif

#ifdef LIBMESH_HAVE_HDF5
        current_files = xdmf_writer.write_time_step(equation_systems, time);
        cout << "[WRITE] " + current_files[0] + " - " + current_files[1] << endl;
#else
        //        ExodusII_IO exo(mesh);
        //        exo.append(true);
        //        exo.write_time_step (exodus_filename, equation_systems, t_step+1, flow_system.time);
        {
            std::ostringstream out;
            out << rname << "_" << std::setw(5) << std::setfill('0') << exodus_step << ".e";
            ExodusII_IO(mesh).write_equation_systems(out.str(), equation_systems);
            exodus_step++;
        }
#endif

#ifdef PROV
        // Mesh Writer
        prov.outputMeshWriter(taskID, simulationID, numberOfWrites, t_step, current_files[1]);
#endif

        int step = t_step + 1;
        if (dim == 2) {
#ifdef PROV
            prov.inputDataExtraction(taskID, simulationID, numberOfWrites, "dataextraction");
#endif

#ifdef PERFORMANCE
            if (libMesh::global_processor_id() == 0) {
                perf.start();
            }
#endif

            char firstFilename[textArraySize];
            sprintf(firstFilename, "init_ext_plane_%d.csv", step);
            char finalFilename[textArraySize];
            sprintf(finalFilename, "ext_plane_%d.csv", step);

#ifdef USE_CATALYST
            FEAdaptor::CoProcess(argc, argv, equation_systems, transport_system.time, step, true, false);
            if (libMesh::global_processor_id() == 0) {
                char commandLine[jsonArraySize];
                sprintf(commandLine, "python clean-csv.py %s %s;rm %s", firstFilename, finalFilename, firstFilename);
                system(strdup(commandLine));
            }
#endif

#ifdef PERFORMANCE
            if (libMesh::global_processor_id() == 0) {
                perf.end();
                double elapsedTime = perf.elapsedTime();
                char buffer[textArraySize];
                sprintf(buffer, "Data Extraction Cost: %.5f", elapsedTime);
                cout << buffer << endl;
                prov.storeDataExtractionCost(elapsedTime);
            }
#endif

#ifdef PROV
            char* extractorName = (char*) malloc(jsonArraySize);
            sprintf(extractorName, "rde%d", numberOfWrites);
            prov.outputDataExtraction(taskID, simulationID, numberOfWrites, "dataextraction", "odataextraction", step, current_files[1], finalFilename, dim, extractorName);
#endif
        } else if (dim == 3) {
            // 3D analysis
            for (int ik = 0; ik <= 3; ik++) {
                char firstFilename[textArraySize];
                sprintf(firstFilename, "init_ext_line_%d_%d.csv", ik, step);
                char finalFilename[textArraySize];
                sprintf(finalFilename, "ext_line_%d_%d.csv", ik, step);

#ifdef PROV
                // Mesh Writer
                char argument1[textArraySize];
                sprintf(argument1, "line%dextraction", ik);
                prov.inputDataExtraction(taskID, simulationID, numberOfWrites, argument1);
#endif

#ifdef PERFORMANCE
                if (libMesh::global_processor_id() == 0) {
                    perf.start();
                }
#endif

#ifdef USE_CATALYST
                if (ik == 0) {
                    FEAdaptor::CoProcess(argc, argv, equation_systems, transport_system.time, step, true, false);
                }
                if (libMesh::global_processor_id() == 0) {
                    char commandLine[jsonArraySize];
                    sprintf(commandLine, "python clean-csv.py %s %s;rm %s", firstFilename, finalFilename, firstFilename);
                    system(strdup(commandLine));
                }
#endif  

#ifdef PERFORMANCE
                if (libMesh::global_processor_id() == 0) {
                    perf.end();
                    double elapsedTime = perf.elapsedTime();
                    char buffer[textArraySize];
                    sprintf(buffer, "Data Extraction Cost: %.5f", elapsedTime);
                    cout << buffer << endl;
                    prov.storeDataExtractionCost(elapsedTime);
                }
#endif  

#ifdef PROV
                // Mesh Writer
                sprintf(argument1, "line%dextraction", ik);
                char argument2[textArraySize];
                sprintf(argument2, "oline%diextraction", ik);
                char* extractorName = (char*) malloc(jsonArraySize);
                sprintf(extractorName, "line%d%d", ik, numberOfWrites);
                prov.outputDataExtraction(taskID, simulationID, numberOfWrites, argument1, argument2, 0, current_files[1], finalFilename, dim, extractorName);
                free(extractorName);
#endif
            }
        }

        char* depID = (char*) malloc(jsonArraySize);
        sprintf(depID, "%d", taskID);
        meshDependencies.push_back(depID);
        free(depID);
    }

    std::cout << "FLOW SOLVER - TOTAL LINEAR ITERATIONS : " << n_linear_iterations_flow << std::endl;
    std::cout << "TRANSPORT SOLVER - TOTAL LINEAR ITERATIONS : " << n_linear_iterations_transport << std::endl;

#ifdef PROV
    // Mesh Aggregator
    char out_filename[256];
    sprintf(out_filename, "%s_%d.xmf", rname.c_str(), libMesh::global_n_processors());
    prov.meshAggregator(simulationID, out_filename, libMesh::global_n_processors(), meshDependencies);
    prov.finishDataIngestor();
#endif

#ifdef PERFORMANCE
    if (libMesh::global_processor_id() == 0) {
        solverPerf.end();
        double elapsedTime = solverPerf.elapsedTime();
        prov.storeSolverCost(elapsedTime);
    }
#endif

    // All done.
    return 0;
}
