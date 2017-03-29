/*
 * File:   sedimentation.cpp
 * Author: camata, vitor
 *
 * Created on December 9, 2014, 3:05 PM
 */
// C++ include files that we need
#include <iostream>
#include <ctime>
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
#include "libmesh/fe_interface.h"

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

#include "FEAdaptor.h"

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
#include "extractor.h"

void PrintStats(EquationSystems& es);

void ComputeMassAndEnergy(EquationSystems& es, std::ostream &out, double &mass, double &mass_dep);

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

void WriteRestartFile(EquationSystems &es, std::string rname) {
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
    // Initialize libMesh.
    LibMeshInit init(argc, argv);

    PerfLog perf_log("Sedimentation Solver");

    // This example requires Adaptive Mesh Refinement support - although
    // it only refines uniformly, the refinement code used is the same
    // underneath
    GetPot command_line(argc, argv);
    std::string input;

    if (command_line.search(1, "-i"))
        input = command_line.next(input);
    else {
        std::cout << "Usage: " << argv[0] << " -i [input file].in -m [gmsh file].msh -o [output file prefix] -d [output file path] -e [extraction script].py -v [visualization script].py" << std::endl;
        libmesh_error_msg("You need specify a input file!");
    }

    GetPot infile(input);

    char meshDependenciesList[256];
    int indexerID = 0;

#ifdef PROV
    Performance solverPerformance;
    solverPerformance.begin();
    Provenance prov(libMesh::global_processor_id());
#endif

    cout << "PROCESSOR ID:" << endl;
    cout << libMesh::global_n_processors() << endl;
    cout << libMesh::global_processor_id() << endl;

    int init_tstep = 1;

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

#ifdef PROV
    prov.inputInputMesh(dim);
#endif
    // Create a mesh object, with dimension to be overridden later,
    // distributed across the default MPI communicator.
    Mesh mesh(init.comm());

    double r_fraction = infile("r_fraction", 0.70);
    double c_fraction = infile("c_fraction", 0.10);
    double max_h_level = infile("max_h_level", 1);
    const unsigned int hlevels = infile("hlevels", 0);
    bool first_step_refinement = infile("first_step_refinement", false);
    bool amrc_flow_transp = infile("amrc_flow_transp", false);

    int ref_interval = infile("r_interval", 1);
    int max_r_steps = infile("max_r_steps", 1);

    MeshRefinement refinement(mesh);
    refinement.refine_fraction() = r_fraction;
    refinement.coarsen_fraction() = c_fraction;
    refinement.max_h_level() = max_h_level;

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
    Real hlock = infile("hlock", 1.0);
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
    equation_systems.parameters.set<Real> ("xlock") = xlock;
    equation_systems.parameters.set<Real> ("hlock") = hlock;
    equation_systems.parameters.set<Real> ("alfa") = alfa;
    equation_systems.parameters.set<Real> ("theta") = theta;
    equation_systems.parameters.set<Real> ("ex") = ex;
    equation_systems.parameters.set<Real> ("ey") = ey;
    equation_systems.parameters.set<Real> ("ez") = ez;
    equation_systems.parameters.set<Real> ("c_factor") = c_factor;
    equation_systems.parameters.set<Real> ("fopc") = fopc;
    equation_systems.parameters.set<Real> ("dt_stab") = infile("dt_stab", 0.1);
    equation_systems.parameters.set<Real> ("erosion/Rp") = infile("erosion/Rp", 0.0);

    // LOOP 
    Real init_time = 0.0;

    //#ifdef PROV
    //    prov.createIndexDirectory();
    //#endif

    // INPUT: TIME INTEGRATION
    Real dt = infile("deltat", 0.005);
    Real tmax = infile("tmax", 0.1);
    unsigned int n_time_steps = infile("n_time_steps", 10);

    //INPUT: NONLINEAR SOLVER
    unsigned int n_nonlinear_steps = infile("n_nonlinear_steps", 10);
    double nonlinear_tolerance = infile("nonlinear_tolerance", 1.0E-03);
    int max_linear_iters = infile("max_linear_iterations", 2000);

    unsigned int write_interval = infile("write_interval", 10);
    std::string rname = "out";
    std::string dpath = "output";
    int numberOfScripts = 0;
    std::string extractionScript = "";
    std::string visualizationScript = "";

    if (command_line.search(1, "-o"))
        rname = command_line.next(rname);

    if (command_line.search(1, "-d"))
        dpath = command_line.next(dpath);

    if (command_line.search(1, "-e"))
        extractionScript = command_line.next(extractionScript);
    numberOfScripts++;

    if (command_line.search(1, "-v"))
        visualizationScript = command_line.next(visualizationScript);
    numberOfScripts++;

    std::cout << "File name: " << rname << endl;
    std::cout << "  path: " << dpath << endl;
    std::cout << "Number of scripts: " << numberOfScripts << endl;
    std::cout << "  Extraction script: " << extractionScript << endl;
    std::cout << "  Visualization script: " << visualizationScript << endl;


    XDMFWriter xdmf_writer(mesh);
    xdmf_writer.set_file_name(rname);
    xdmf_writer.set_dir_path(dpath);

    if (!is_file_exist("restart.run")) {

        string mesh_file;
        if (command_line.search(1, "-m"))
            mesh_file = command_line.next(mesh_file);
        else {
            std::cout << "Usage: " << argv[0] << " -i [input file].in -m [gmsh file].msh" << std::endl;
            libmesh_error_msg("You need specify a mesh file!");
        }

        std::cout << "Opening file: " << mesh_file << std::endl;

        mesh.read(mesh_file);

        refinement.uniformly_refine(hlevels);

        equation_systems.parameters.set<int> ("dim") = mesh.mesh_dimension();

        sediment_flow.init();
        sediment_transport.init();
        sediment_deposition.init();

        sediment_flow.setup(infile);
        sediment_transport.setup(infile);
        sediment_deposition.setup(infile);

#ifdef MESH_MOVIMENT  
        moving_mesh.init();
        moving_mesh.setup(infile);
#endif

#ifdef PROV
        prov.outputInputMesh(r_fraction, c_fraction, max_h_level, hlevels);
#endif
        // Initialize the data structures for the equation system.
        equation_systems.init();

#ifdef PROV
        prov.outputCreateEquationSystems(Reynolds, Gr, Sc, Us, Diffusivity, xlock, alfa, theta, ex, ey, ez, c_factor);
#endif

    } else {

        GetPot restart("restart.in");
        const string mesh_restart = restart("mesh_restart", "0");
        const string solution_restart = restart("solution_restart", "0");
        init_time = restart("time", 0.0);
        dt = restart("dt", 0.0);
        init_tstep = restart("init_tstep", 0);

#ifdef PROV
        indexerID = restart("indexerID", 0);
#endif

        sediment_transport.init_mass = restart("initial_mass", 0.0);
        sediment_transport.mass_dep = restart("mass_dep", 0.0);


        int xdmf_file_id = restart("xdmf_file_id", 0);
        xdmf_writer.set_file_id(xdmf_file_id);


        mesh.read(mesh_restart);


        equation_systems.read(solution_restart, READ);
        equation_systems.parameters.set<int> ("dim") = mesh.mesh_dimension();

        sediment_flow.setup(infile);
        sediment_transport.setup(infile);
        sediment_deposition.setup(infile);

#ifdef MESH_MOVIMENT
        moving_mesh.setup(infile);
#endif

        // Get a reference to the Convection-Diffusion system object.
        TransientLinearImplicitSystem & flow_system =
                equation_systems.get_system<TransientLinearImplicitSystem> ("flow");

        // Get a reference to the Convection-Diffusion system object.
        TransientLinearImplicitSystem & transport_system =
                equation_systems.get_system<TransientLinearImplicitSystem> ("sediment");
        //transport_system.add_vector("volume");


        ExplicitSystem & deposition_system = equation_systems.get_system<ExplicitSystem>("deposition");


        flow_system.update();
        transport_system.update();
        deposition_system.update();
        equation_systems.update();


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
    perf_log.start_event("XDMF:Write");
    current_files = xdmf_writer.write_time_step(equation_systems, time);
    perf_log.stop_event("XDMF:Write");
    cout << "[WRITE] " + current_files[0] + " - " + current_files[1] << endl;

#ifdef PROV
    prov.outputGetMaximumIterations(dt, tmax, n_time_steps, n_nonlinear_steps, nonlinear_tolerance, max_linear_iters, max_r_steps, write_interval, current_files[1]);
    Performance performance;
#endif

#ifdef USE_CATALYST    
    if (dim == 2) {
        // 2D analysis        
#ifdef PROV
        prov.inputInitDataExtraction(-1);
        performance.begin();
#endif
        perf_log.start_event("CATALYST", "Init");
        FEAdaptor::Initialize(numberOfScripts, extractionScript, visualizationScript);
        perf_log.stop_event("CATALYST", "Init");
        perf_log.start_event("CATALYST", "CoProcess");
        FEAdaptor::CoProcess(numberOfScripts, extractionScript, visualizationScript, equation_systems, 0.0, t_step, write_interval, false, false);
        perf_log.stop_event("CATALYST", "CoProcess");
#ifdef PROV
        performance.end();
        prov.storeCatalystCost(0, 1, performance.getElapsedTime());
        extractor::invoke2DRawDataExtractor(libMesh::global_processor_id(), t_step);
        indexerID++;
        prov.outputInitDataExtraction(-1, current_files[1], dim, indexerID);
#endif        
    } else if (dim == 3) {
        // 3D analysis
        for (int lineID = 0; lineID <= 3; lineID++) {
#ifdef PROV
            prov.inputInitVisualization(lineID);
            prov.inputInitDataExtraction(lineID);
            performance.begin();
#endif
            if (lineID == 0) {
                perf_log.start_event("CATALYST", "Init");
                FEAdaptor::Initialize(numberOfScripts, extractionScript, visualizationScript);
                perf_log.stop_event("CATALYST", "Init");
                perf_log.start_event("CATALYST", "CoProcess");
                FEAdaptor::CoProcess(numberOfScripts, extractionScript, visualizationScript, equation_systems, 0.0, t_step, write_interval, false, false);
                perf_log.stop_event("CATALYST", "CoProcess");
            }
#ifdef PROV
            performance.end();
            extractor::invoke3DRawDataExtractor(libMesh::global_processor_id(), t_step, lineID);
            prov.storeCatalystCost(0, 1, performance.getElapsedTime());
            indexerID++;
            prov.outputInitVisualization(lineID, t_step);
            prov.outputInitDataExtraction(lineID, current_files[1], dim, indexerID);
#endif
        }
    }
#endif

    // STEP LOOP
    // Loop in time steps
    // TODO: PROV 
    int taskID = 0;
    int numberIterationsFluid = 0;
    int numberIterationsSediments = 0;
    int numberIterationsMeshRefinements = 0;
    int subTaskID = 0;
    vector<string> meshDependencies;

    equation_systems.parameters.set<PerfLog*> ("PerfLog") = &perf_log;

    for (t_step = init_tstep; (t_step < n_time_steps)&&(time < tmax); t_step++) {
        taskID++;
        if (is_file_exist("abort.run")) break;

        if (is_file_exist("reset.run")) {

            GetPot reset(input);
            dt = reset("deltat", 0.005);
            tmax = reset("tmax", 0.1);
            n_time_steps = reset("n_time_steps", 10);
            n_nonlinear_steps = reset("n_nonlinear_steps", 10);
            write_interval = reset("write_interval", 10);
            nonlinear_tolerance = reset("nonlinear_tolerance", 1.0E-03);
            max_linear_iters = reset("max_linear_iterations", 2000);

        }

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
                // TODO: PROV
                numberIterationsFluid++;
#ifdef PROV
                prov.inputSolverSimulationFluid(taskID, numberIterationsFluid);
#endif
                // Update the nonlinear solution.
                flow_last_nonlinear_soln->zero();
                flow_last_nonlinear_soln->add(*flow_system.solution);

                // Assemble & solve the linear system.
                perf_log.start_event("Solver", "Flow");
                flow_system.solve();
                perf_log.stop_event("Solver", "Flow");

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
                prov.outputSolverSimulationFluid(taskID, numberIterationsFluid, t_step, transport_system.time, r, l, n_linear_iterations, final_linear_residual, norm_delta, norm_delta / u_norm, converged);
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
                //TODO: PROV
                numberIterationsSediments++;
#ifdef PROV
                prov.inputSolverSimulationSediments(taskID, numberIterationsSediments);
#endif
                // Update the nonlinear solution.
                sed_last_nonlinear_soln->zero();
                sed_last_nonlinear_soln->add(*transport_system.solution);

                // Assemble & solve the linear system.
                perf_log.start_event("Solver", "Transport");
                transport_system.solve();
                perf_log.stop_event("Solver", "Transport");
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
                prov.outputSolverSimulationSediments(taskID, numberIterationsSediments, t_step, transport_system.time, r, l, n_linear_iterations, final_linear_residual, norm_delta, norm_delta / u_norm, converged);
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
                numberIterationsMeshRefinements++;
                std::cout << "\n****************** Mesh Refinement ********************  " << std::endl;
                std::cout << " Considering Transport" << ((amrc_flow_transp && !first_step_refinement) ? " & Flow Variables\n" : " Variable\n");
                std::cout << "Number of elements before AMR step: " << mesh.n_active_elem() << std::endl;

                int beforeNActiveElem = mesh.n_active_elem();
                ErrorVector error_flow, error;
                KellyErrorEstimator error_estimator_flow;
                KellyErrorEstimator error_estimator_transp;

                perf_log.start_event("estimate_error", "AMR");
                // First compute error for transport only
                error_estimator_transp.estimate_error(transport_system, error);

                // Now compute the error also for weighted flow variables but only if not "first_step_refinement"
                if (amrc_flow_transp && !first_step_refinement) {
                    // Weights to compute error for the Flow variables (to exclude pressure)
                    std::vector<Real> weights(flow_system.n_vars(), 1.0);
                    weights[flow_system.n_vars() - 1] = 0.0; // excluding pressure for while
                    error_estimator_flow.error_norm = SystemNorm(std::vector<FEMNormType>(flow_system.n_vars(), error_estimator_flow.error_norm.type(0)), weights);
                    error_estimator_flow.estimate_error(flow_system, error_flow);

                    libmesh_assert(error_flow == error);
                    for (unsigned int i = 0; i < error.size(); i++)
                        error[i] += error_flow[i];
                }

                perf_log.stop_event("estimate_error", "AMR");

                perf_log.start_event("flag_elements", "AMR");
                refinement.flag_elements_by_error_fraction(error);
                perf_log.stop_event("flag_elements", "AMR");

                perf_log.start_event("refine_and_coarse", "AMR");
                refinement.refine_and_coarsen_elements();
                perf_log.stop_event("refine_and_coarse", "AMR");

                //equation_systems.update();
                perf_log.start_event("reinit systems", "AMR");
                equation_systems.reinit();
                perf_log.stop_event("reinit systems", "AMR");

                redo_nl = true;
                first_step_refinement = false;
#ifdef USE_CATALYST
                FEAdaptor::mark_to_rebuild_grid();
#endif 
                std::cout << "Number of elements after AMR step: " << mesh.n_active_elem() << std::endl;
#ifdef PROV
                prov.outputMeshRefinement(taskID, numberIterationsMeshRefinements, first_step_refinement, t_step, beforeNActiveElem, mesh.n_active_elem());

#endif

                xdmf_writer.mesh_changed_on();
            }

            sediment_deposition.ComputeDeposition();
            //sediment_deposition.print();

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

            // Output every write_interval timesteps to file.
            if ((t_step + 1) % write_interval == 0) {
                //TODO: PROV
                subTaskID++;

                const std::string mesh_restart = rname + "_mesh_restart.xdr";
                const std::string solution_restart = rname + "_solution_restart.xdr";

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
                fout << "xdmf_file_id = " << xdmf_writer.get_file_id() << std::endl;
#ifdef PROV
                fout << "indexerID = " << indexerID << std::endl;

#endif         

                fout.close();
#ifdef PROV
                prov.inputMeshWriter(taskID, subTaskID);
#endif

                perf_log.start_event("Write", "XDMF_IO");
                current_files = xdmf_writer.write_time_step(equation_systems, time);
                perf_log.stop_event("Write", "XDMF_IO");
                cout << "[WRITE] " + current_files[0] + " - " + current_files[1] << endl;

#ifdef PROV
                prov.outputMeshWriter(taskID, subTaskID, t_step, current_files[1]);
#endif

                int step = t_step + 1;

                if (!redo_nl) {
#ifdef USE_CATALYST
                    if (dim == 2) {
#ifdef PROV
                        prov.inputDataExtraction(taskID, subTaskID, -1);
                        performance.begin();
#endif
                        perf_log.start_event("CATALYST:CoProcess");
                        FEAdaptor::CoProcess(numberOfScripts, extractionScript, visualizationScript, equation_systems, transport_system.time, step, write_interval, false, false);
                        perf_log.stop_event("CATALYST:CoProcess");
#ifdef PROV
                        performance.end();
                        prov.storeCatalystCost(taskID, subTaskID, performance.getElapsedTime());
                        extractor::invoke2DRawDataExtractor(libMesh::global_processor_id(), step);
                        indexerID++;
                        prov.outputDataExtraction(taskID, subTaskID, -1, step, current_files[1], dim, indexerID);
#endif
                    } else if (dim == 3) {
                        // 3D analysis
                        for (int lineID = 0; lineID <= 3; lineID++) {
#ifdef PROV
                            prov.inputVisualization(lineID, taskID);
                            prov.inputDataExtraction(taskID, subTaskID, lineID);
                            performance.begin();
#endif
                            if (lineID == 0) {
                                perf_log.start_event("CATALYST:CoProcess");
                                FEAdaptor::CoProcess(numberOfScripts, extractionScript, visualizationScript, equation_systems, transport_system.time, step, write_interval, false, false);
                                perf_log.stop_event("CATALYST:CoProcess");
                            }
#ifdef PROV
                            performance.end();
                            extractor::invoke3DRawDataExtractor(libMesh::global_processor_id(), step, lineID);
                            prov.storeCatalystCost(taskID, subTaskID, performance.getElapsedTime());
                            indexerID++;
                            prov.outputVisualization(lineID, taskID, step);
                            prov.outputDataExtraction(taskID, subTaskID, lineID, step, current_files[1], dim, indexerID);
#endif
                        }
                        sprintf(meshDependenciesList, "%d", taskID);
                        meshDependencies.push_back(meshDependenciesList);
                    }
#endif
                }
            }
        }

        PrintStats(equation_systems);

    }

    if ((t_step + 1) % write_interval != 0) {
        //TODO: Prov
        subTaskID++;
#ifdef PROV
        prov.inputMeshWriter(taskID, subTaskID);
#endif


        perf_log.start_event("Write", "XDMF_IO");
        current_files = xdmf_writer.write_time_step(equation_systems, time);
        perf_log.stop_event("Write", "XDMF_IO");
        cout << "[WRITE] " + current_files[0] + " - " + current_files[1] << endl;


#ifdef PROV
        prov.outputMeshWriter(taskID, subTaskID, t_step, current_files[1]);
#endif

        int step = t_step + 1;
#ifdef USE_CATALYST
        if (dim == 2) {
#ifdef PROV
            prov.inputDataExtraction(taskID, subTaskID, -1);
            performance.begin();
#endif
            perf_log.start_event("CATALYST:CoProcess");
            FEAdaptor::CoProcess(numberOfScripts, extractionScript, visualizationScript, equation_systems, transport_system.time, step, write_interval, true, false);
            perf_log.stop_event("CATALYST:CoProcess");
#ifdef PROV
            performance.end();
            prov.storeCatalystCost(taskID, subTaskID, performance.getElapsedTime());

            extractor::invoke2DRawDataExtractor(libMesh::global_processor_id(), step);
            indexerID++;
            prov.outputDataExtraction(taskID, subTaskID, -1, step, current_files[1], dim, indexerID);
#endif
        } else if (dim == 3) {
            // 3D analysis
            for (int lineID = 0; lineID <= 3; lineID++) {
#ifdef PROV
                prov.inputVisualization(lineID, taskID);
                prov.inputDataExtraction(taskID, subTaskID, lineID);
                performance.begin();
#endif
                if (lineID == 0) {
                    perf_log.start_event("CATALYST:CoProcess");
                    FEAdaptor::CoProcess(numberOfScripts, extractionScript, visualizationScript, equation_systems, transport_system.time, step, write_interval, true, false);
                    perf_log.stop_event("CATALYST:CoProcess");
                }
#ifdef PROV
                performance.end();
                extractor::invoke3DRawDataExtractor(libMesh::global_processor_id(), step, lineID);
                prov.storeCatalystCost(taskID, subTaskID, performance.getElapsedTime());
                indexerID++;
                prov.outputVisualization(lineID, taskID, step);
                prov.outputDataExtraction(taskID, subTaskID, lineID, step, current_files[1], dim, indexerID);
#endif
            }
        }
#endif
        sprintf(meshDependenciesList, "%d", taskID);
        meshDependencies.push_back(meshDependenciesList);
    }

    std::cout << "FLOW SOLVER - TOTAL LINEAR ITERATIONS : " << n_linear_iterations_flow << std::endl;
    std::cout << "TRANSPORT SOLVER - TOTAL LINEAR ITERATIONS : " << n_linear_iterations_transport << std::endl;

#ifdef PROV
    //    #ifdef USE_CATALYST
    //        sprintf(memalloc, "rm video.mp4;cat image_*.png | ffmpeg -i - -r 30 video.mp4");
    //        cout << memalloc << endl;
    //        system(memalloc);
    //    #endif

    char out_filename[256];
    sprintf(out_filename, "%s_%d.xmf", rname.c_str(), libMesh::global_n_processors());
    sort(meshDependencies.begin(), meshDependencies.end());
    meshDependencies.erase(unique(meshDependencies.begin(), meshDependencies.end()), meshDependencies.end());
    prov.meshAggregator(out_filename, libMesh::global_n_processors(), meshDependencies);
    prov.finishDataIngestor();
#endif

#ifdef PROV
    solverPerformance.end();
    prov.storeSolverCost(solverPerformance.getElapsedTime());
#endif

    // All done.
    cout << "All done!" << endl;
    return 0;
}

void PrintStats(EquationSystems& es) {
    std::vector<string> vname;
    std::vector<Real> local_min;
    std::vector<Real> local_max;
    std::vector<Real> local_l2norm;
    std::vector<unsigned int> vars;
    std::vector<unsigned int> sys;


    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();


    const System & flow_system = es.get_system("flow");
    const System & sediment_system = es.get_system("sediment");

    vars.push_back(flow_system.variable_number("u"));
    sys.push_back(flow_system.number());
    vname.push_back("u");

    vars.push_back(flow_system.variable_number("v"));
    sys.push_back(flow_system.number());
    vname.push_back("v");

    if (dim > 2) {
        vars.push_back(flow_system.variable_number("w"));
        sys.push_back(flow_system.number());
        vname.push_back("w");
    }

    vars.push_back(flow_system.variable_number("p"));
    sys.push_back(flow_system.number());
    vname.push_back("p");

    vars.push_back(sediment_system.variable_number("s"));
    sys.push_back(sediment_system.number());
    vname.push_back("s");

    const int nvars = vname.size();

    local_min.resize(nvars);
    local_max.resize(nvars);
    local_l2norm.resize(nvars);

    for (int v = 0; v < nvars; v++) {
        local_min[v] = 1.0E10;
        local_max[v] = -1.0E10;
        local_l2norm[v] = 0.0;
    }


    const DofMap & dof_map_flow = flow_system.get_dof_map();
    const DofMap & dof_map_sediment = sediment_system.get_dof_map();

    /*
    for(int v = 0; v < nvars-1; v++) {
        std::vector<unsigned int> dof_indices;
        dof_map_flow.local_variable_indices(dof_indices, mesh, vars[v]);
        for(int idof; idof  < dof_indices.size(); idof++)
        {
            Real value = flow_system.solution->el(idof);
            local_min[v]     = std::min(local_min[v], value);
            local_max[v]     = std::max(local_min[v], value);
            local_l2norm[v]  += value*value;
        }
        
    }
    
    {
         std::vector<unsigned int> dof_indices;
        dof_map_sediment.local_variable_indices(dof_indices, mesh, vars[nvars-1]);
        for(int idof; idof  < dof_indices.size(); idof++)
        {
            Real value = flow_system.solution->el(idof);
            local_min[nvars-1]     = std::min(local_min[nvars-1], value);
            local_max[nvars-1]     = std::max(local_min[nvars-1], value);
            local_l2norm[nvars-1]  += value*value;
        }
    }
     */

    MeshBase::const_node_iterator iter = mesh.local_nodes_begin();
    for (; iter != mesh.local_nodes_end(); iter++) {

        const Node* node = (*iter);
        for (int v = 0; v < nvars; v++) {
            int dof = node->dof_number(sys[v], vars[v], 0);
            const System & s = es.get_system(sys[v]);
            if (dof >= s.solution->first_local_index() && dof < s.solution->last_local_index()) {
                Real value = s.solution->el(dof);
                local_min[v] = std::min(local_min[v], value);
                local_max[v] = std::max(local_max[v], value);
                local_l2norm[v] += value*value;
            }
        }


    }


    std::vector<Real> global_min(vars.size());
    std::vector<Real> global_max(vars.size());
    std::vector<Real> global_l2norm(vars.size());

    MPI_Reduce(&local_min[0], &global_min[0], nvars, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_max[0], &global_max[0], nvars, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_l2norm[0], &global_l2norm[0], nvars, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    //if(libMesh::global_processor_id() == 0) 
    // {

    //           123456789012345678901234567890123456789012345678901234567890123456789
    std::cout << "______________________________________________________" << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << "   S I M U L A T I O N   I N F O                      " << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    std::cout << " Vars      Min             Max             L2-Norm    " << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;
    for (int v = 0; v < vars.size(); v++) {
        std::ostringstream out;
        out.precision(5);
        out << std::left << std::setfill(' ') << std::setw(10) << vname[v]
                << std::setw(15) << std::scientific << global_min[v]
                << std::setw(15) << global_max[v]
                << std::setw(15) << std::sqrt(global_l2norm[v]) << "\n";
        std::cout << out.str() << std::flush;
    }
    std::cout << "______________________________________________________" << std::endl;


}