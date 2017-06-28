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

#include "xdmf.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;
using namespace std;

#ifndef PROVENANCE_H
#define PROVENANCE_H
#include "include/provenance/provenance.h"
#endif

#include "performance.h"
#include "contrib/dfanalyzer/extractor.h"
#include "viscosity_flow.h"
#include "include/catalyst/FEAdaptor.h"

#include "mesh_moviment.h"
#include "sedimentation_transport.h"
#include "sedimentation_deposition.h"
#include "sedimentation_flow.h"

#include "timeStepControlBase.h"
#include "timeStepControlPID.h"
#include "timeStepControlResidual.h"
#include "timeStepControlPC11.h"

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

bool is_file_exist(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

void InSituCatalystInitCoprocessing(int dim, int numberOfScripts,
        string &extractionScript, string& visualizationScript,
        int t_step, int write_interval, string* current_files,
        EquationSystems & equation_systems,
        PerfLog& perf_log, Provenance &provenance);



void InSituCatalystCoprocessing(int dim, int numberOfScripts,
        string &extractionScript, string& visualizationScript,
        int t_step, double time, int write_interval, string* current_files,
        EquationSystems & equation_systems,
        PerfLog& perf_log, Provenance &provenance);



// We can now begin the main program.  Note that this
// example will fail if you are using complex numbers
// since it was designed to be run only with real numbers.

int main(int argc, char** argv) {
    // Initialize libMesh.
    LibMeshInit init(argc, argv);

    //    Restart
    bool restartControl = is_file_exist("restart.run");

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

    string mesh_file;
    if (command_line.search(1, "-m"))
        mesh_file = command_line.next(mesh_file);
    else {
        std::cout << "Usage: " << argv[0] << " -i [input file].in -m [gmsh file].msh" << std::endl;
        libmesh_error_msg("You need specify a mesh file!");
    }

    GetPot infile(input);
    double v = infile("version", 1.0);
    if (v == 1.0) {
        libmesh_error_msg("This input file is not compatible with this libMesh-Sedimentation version !");
    }

    perf_log.start_event("Init", "Provenance");
    Provenance provenance(libMesh::global_processor_id());
    perf_log.stop_event("Init", "Provenance");

#ifdef PROVENANCE
    Performance solverPerformance;
    solverPerformance.begin();
    provenance.SetUp();
    provenance.resetIndexerID();
#endif

    //Define performance parameters
    unsigned int n_flow_nonlinear_iterations_total = 0;
    unsigned int n_flow_linear_iterations_total = 0;

    unsigned int n_transport_nonlinear_iterations_total = 0;
    unsigned int n_transport_linear_iterations_total = 0;

    unsigned int n_rejected_flow_nonlinear_iterations_total = 0;
    unsigned int n_rejected_flow_linear_iterations_total = 0;

    unsigned int n_rejected_transport_nonlinear_iterations_total = 0;
    unsigned int n_rejected_transport_linear_iterations_total = 0;

    int dim = infile("dim", 3);

#ifdef PROVENANCE
    perf_log.start_event("InputMesh", "Provenance");
    provenance.inputInputMesh();
    perf_log.stop_event("InputMesh", "Provenance");
#endif

    // Create a mesh object, with dimension to be overridden later,
    // distributed across the default MPI communicator.
    Mesh mesh(init.comm());

#ifdef PROVENANCE
    perf_log.start_event("InputMesh", "Provenance");
    provenance.outputInputMesh(dim, mesh_file, restartControl);
    perf_log.stop_event("InputMesh", "Provenance");
#endif

    // Getting mesh refinement parameters
    double r_fraction = infile("amr/r_fraction", 0.70);
    double c_fraction = infile("amr/c_fraction", 0.10);
    double max_h_level = infile("amr/max_h_level", 1);
    const unsigned int initial_unif_ref_mesh = infile("amr/initial_unif_ref_mesh", 0);
    bool first_step_refinement = infile("amr/first_step_refinement", false);
    bool amrc_flow_transp = infile("amr/amrc_flow_transp", false);
    int ref_interval = infile("amr/r_interval", 1);
    int max_r_steps = infile("amr/max_r_steps", 1);
    bool use_amr = (max_r_steps > 1);

    MeshRefinement refinement(mesh);
    refinement.refine_fraction() = r_fraction;
    refinement.coarsen_fraction() = c_fraction;
    refinement.max_h_level() = max_h_level;

#ifdef PROVENANCE
    perf_log.start_event("AMRConfig", "Provenance");
    provenance.outputAMRConfig(r_fraction, c_fraction, max_h_level, initial_unif_ref_mesh, first_step_refinement,
            amrc_flow_transp, ref_interval, max_r_steps);
    perf_log.stop_event("AMRConfig", "Provenance");
#endif

    // Create an equation systems object.
    EquationSystems equation_systems(mesh);

    SedimentationFlow sediment_flow(equation_systems);
    SedimentationTransport sediment_transport(equation_systems);
    SedimentationDeposition sediment_deposition(equation_systems);

#ifdef MESH_MOVIMENT
    MeshMoviment moving_mesh(equation_systems);
#endif

    //Getting physical parameters
    Real Re = infile("Reynolds", 1.0E03);
    Real Gr = infile("Grashof", 1.0E05);
    Real Sc = infile("Schmidt", 0.71);
    Real Ri = infile("Richardson", 1.0);
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

    if (Re == 0.0) {
        Re = std::sqrt(Gr);
        Ri = 1.0;
    }

    Real Df = 1.0 / (Sc * Re);

    cout << "\nDimensionless Physical Parameters" << endl;
    cout << "---------------------------------" << endl;
    cout << " Reynolds    : " << Re << endl;
    cout << " Grashof     : " << Gr << endl;
    cout << " Schmidt     : " << Sc << endl;
    cout << " Diffusivity : " << Df << endl;
    cout << " Richardson  : " << Ri << endl;

    equation_systems.parameters.set<Real> ("Reynolds") = Re;
    equation_systems.parameters.set<Real> ("Diffusivity") = Df;
    equation_systems.parameters.set<Real> ("Richardson") = Ri;
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
    equation_systems.parameters.set<Real> ("erosion/Rp") = infile("erosion/Rp", 0.0);
    equation_systems.parameters.set<Real> ("Cs") = infile("Cs", 0.0);

#ifdef PROVENANCE
    perf_log.start_event("CreateEquationSystems", "Provenance");
    provenance.outputCreateEquationSystems(Re, Gr, Sc, Us, Df, xlock, alfa, theta, ex, ey, ez, c_factor);
    perf_log.stop_event("CreateEquationSystems", "Provenance");
#endif

    // LOOP 
    Real init_time = 0.0;
    int init_tstep = 0;
    unsigned int t_step = 0;

    // INPUT: TIME INTEGRATION
    Real dt_init = infile("time/deltat", 0.005);
    unsigned int n_time_steps = infile("time/n_time_steps", 1000);
    Real tmax = infile("time/tmax", 0.1);
    equation_systems.parameters.set<Real> ("tmax") = tmax;

    // Time step control input data
    //-----------------------------
    std::string ts_control_model_name = infile("ts_control/model_name", "PC11");

    double dt_min = infile("ts_control/dt_min", 0.001);
    double dt_max = infile("ts_control/dt_max", 0.1);
    double tol_u = infile("ts_control/tol_u", 0.1);
    double tol_s = infile("ts_control/tol_s", 0.1);
    double kp = infile("ts_control/kp", 0.075);
    double ki = infile("ts_control/ki", 0.175);
    double kd = infile("ts_control/kd", 0.01);
    unsigned int nsa_max = infile("ts_control/nsa_max", 10);
    unsigned int nsa_target_flow = infile("ts_control/nsa_target_flow", 4);
    unsigned int nsa_target_transport = infile("ts_control/nsa_target_transport", 4);
    unsigned int nsa_limit_flow = infile("ts_control/nsa_limit_flow", 8);
    unsigned int nsa_limit_transport = infile("ts_control/nsa_limit_transport", 8);
    double mult_factor_max = infile("ts_control/mult_factor_max", 2.0);
    double mult_factor_min = infile("ts_control/mult_factor_min", 0.7);
    double pc11_theta = infile("ts_control/pc11_theta", 1.0);
    double alpha = infile("ts_control/alpha", 0.8);
    double k_exp = infile("ts_control/k_exp", 2.0);
    double s_min = infile("ts_control/s_min", 0.1);
    double s_max = infile("ts_control/s_max", 2.0);
    double reduct_factor = infile("ts_control/reduct_factor", 0.5);
    bool complete_flow_norm = infile("ts_control/complete_flow_norm", false);

    // Stabilization Parameters
    equation_systems.parameters.set<Real> ("dt_stab") = infile("stabilization/dt_stab", 0.1);
    equation_systems.parameters.set<Real> ("s_ref_bar_yzBeta") = infile("stabilization/s_ref_bar_yzBeta", 1.0);
    equation_systems.parameters.set<Real> ("delta_transient_factor") = infile("stabilization/delta_transient_factor", 0.5);
    std::string fem_model = infile("stabilization/fem_model", "");
    equation_systems.parameters.set<std::string> ("fem_model") = fem_model;
    bool yzBeta = infile("stabilization/yzBeta", true);
    equation_systems.parameters.set<bool> ("yzBeta") = yzBeta;

#ifdef PROVENANCE
    provenance.outputTSControlConfig(ts_control_model_name, dt_min, dt_max, tol_u, tol_s,
            kp, ki, kd, nsa_max, nsa_target_flow, nsa_target_transport,
            nsa_limit_flow, nsa_limit_transport, mult_factor_max, mult_factor_min,
            pc11_theta, alpha, k_exp, s_min, s_max, reduct_factor, complete_flow_norm);
#endif

    // Linear and non-linear solver parameters 
    int flow_n_nonlinear_steps = infile("flow_n_nonlinear_steps", 10);
    int transport_n_nonlinear_steps = infile("transport_n_nonlinear_steps", 10);
    double flow_nonlinear_tolerance = infile("flow_nonlinear_tolerance", 1.0E-03);
    double transport_nonlinear_tolerance = infile("transport_nonlinear_tolerance", 1.0E-03);
    int max_linear_iter = infile("max_linear_iterations", 2000);
    double flow_initial_linear_solver_tol = infile("flow_initial_linear_solver_tolerance", 1.0E-3);
    double transport_initial_linear_solver_tol = infile("transport_initial_linear_solver_tolerance", 1.0E-4);
    double minimum_linear_solver_tol = infile("minimum_linear_solver_tolerance", 1.0E-8);
    double linear_tolerance_power = infile("linear_tolerance_power", 2.0);

    sediment_flow.non_linear_tolerance() = flow_nonlinear_tolerance;
    sediment_flow.initial_linear_tolerance() = flow_initial_linear_solver_tol;
    sediment_flow.max_nonlinear_iteractions() = flow_n_nonlinear_steps;
    sediment_flow.linear_tolerance_power() = linear_tolerance_power;

    sediment_transport.nonlinear_tolerance() = transport_nonlinear_tolerance;
    sediment_transport.max_nonlinear_iteractions() = transport_n_nonlinear_steps;
    sediment_transport.initial_linear_tolerance() = transport_initial_linear_solver_tol;
    sediment_transport.linear_tolerance_power() = linear_tolerance_power;

    unsigned int n_flow_sstate = infile("n_flow_sstate", 50);
    unsigned int n_transport_sstate = infile("n_transport_sstate", 50);

    // Output parameters
    unsigned int write_interval = infile("write_interval", 10);
    unsigned int catalyst_interval = infile("catalyst_interval", 10);
    bool write_restart = infile("write_restart", false);

    std::string rname = "out";
    std::string dpath = "output";

    int numberOfScripts = 0;
    std::string extractionScript = "not defined";
    std::string visualizationScript = "not defined";

    if (command_line.search(1, "-o"))
        rname = command_line.next(rname);

    if (command_line.search(1, "-d"))
        dpath = command_line.next(dpath);

    if (command_line.search(1, "-e")) {
        extractionScript = command_line.next(extractionScript);
        numberOfScripts++;
    }

    if (command_line.search(1, "-v")) {
        visualizationScript = command_line.next(visualizationScript);
        numberOfScripts++;
    }

    std::cout << "\nOutput data" << endl;
    cout << "-----------" << endl;
    std::cout << " File name           : " << rname << endl;
    std::cout << " Path                : " << dpath << endl;
    std::cout << " Write interval      : " << to_string(write_interval) << endl;
    std::cout << " Catalyst interval   : " << to_string(catalyst_interval) << endl;
    std::cout << " Number of scripts   : " << numberOfScripts << endl;
    std::cout << " Extraction script   : " << extractionScript << endl;
    std::cout << " Visualization script: " << visualizationScript << endl;

#ifdef PROVENANCE
    provenance.outputIOConfig(dpath, rname, write_interval, catalyst_interval, write_restart);
#endif

    XDMFWriter xdmf_writer(mesh);
    xdmf_writer.set_file_name(rname);
    xdmf_writer.set_dir_path(dpath);

    if (!restartControl) {

        // Starting simulation from the beginning...
        std::cout << "\nOpening mesh file: " << mesh_file << std::endl;

        mesh.read(mesh_file);

        refinement.uniformly_refine(initial_unif_ref_mesh);

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

        // Initialize the data structures for the equation system.
        equation_systems.init();

    } else {

        GetPot restart("restart.in");
        const string mesh_restart = restart("mesh_restart", "0");
        const string solution_restart = restart("solution_restart", "0");
        init_time = restart("time", 0.0);
        dt_init = restart("dt", 0.0);
        init_tstep = restart("init_tstep", 0);

#ifdef PROVENANCE
        provenance.setIndexerID(restart("indexerID", 0));
#endif

        sediment_transport.init_mass = restart("initial_mass", 0.0);
        sediment_transport.mass_dep = restart("mass_dep", 0.0);

        int xdmf_file_id = restart("xdmf_file_id", 0);
        xdmf_writer.set_file_id(xdmf_file_id);

        n_flow_nonlinear_iterations_total = restart("n_flow_nonlinear_iterations_total", 0);
        n_transport_nonlinear_iterations_total = restart("n_transport_nonlinear_iterations_total", 0);
        n_flow_linear_iterations_total = restart("n_flow_linear_iterations_total", 0);
        n_transport_linear_iterations_total = restart("n_transport_linear_iterations_total", 0);
        n_rejected_flow_nonlinear_iterations_total = restart("n_rejected_flow_nonlinear_iterations_total", 0);
        n_rejected_transport_nonlinear_iterations_total = restart("n_rejected_transport_nonlinear_iterations_total", 0);
        n_rejected_flow_linear_iterations_total = restart("n_rejected_flow_linear_iterations_total", 0);
        n_rejected_transport_linear_iterations_total = restart("n_rejected_transport_linear_iterations_total", 0);


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
                equation_systems.get_system<TransientLinearImplicitSystem> ("transport");
        //transport_system.add_vector("volume");

        ExplicitSystem & deposition_system = equation_systems.get_system<ExplicitSystem>("deposition");

        flow_system.update();
        transport_system.update();
        deposition_system.update();
        equation_systems.update();
    }

    // Print information about the mesh to the screen.
    mesh.print_info();
    
    if (initial_unif_ref_mesh)
        cout<<"  Applying "<<initial_unif_ref_mesh<<" level(s) of initial refinement throughout the mesh\n\n";

    // Get a reference to the Convection-Diffusion system object.
    TransientLinearImplicitSystem & transport_system =
            equation_systems.get_system<TransientLinearImplicitSystem> ("transport");

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

    // defining an instance of Time-step Control
    timeStepControlBase* ts_control = 0;

    bool control_ts = false;
    double dt;
    bool accepted = true;

    std::ofstream foutDt, foutMass;
    std::string out_dat_name;
    if (ts_control_model_name == "PID") {
        ts_control = new timeStepControlPID(dt_init, dt_min, dt_max, nsa_max, tol_u, tol_s, kp, ki, kd);
        dt = dt_init;
        control_ts = true;
        if (mesh.processor_id() == 0) {
            out_dat_name = rname + "_dtXtime_PID.dat";
            foutDt.open(out_dat_name);
            out_dat_name = rname + "_massXtime_PID.dat";
            foutMass.open(out_dat_name);
        }
    }/*
        else if (ts_control_model_name == "RES") {
            ts_control = new timeStepControlResidual(dt_init, dt_min, dt_max, nsa_target_flow, nsa_target_transport, nsa_max, nsa_limit_flow, nsa_limit_transport, mult_factor_max, mult_factor_min, reduct_factor);
            // to ensure not to use number NLI ratio for the first time-step control
            equation_systems.parameters.set<unsigned int> ("old_n_non_linear_iter_flow")      = flow_n_nonlinear_steps + 1;
            equation_systems.parameters.set<unsigned int> ("old_n_non_linear_iter_transport") = transport_n_nonlinear_steps + 1;
            dt = dt_init;
            control_ts = true;
            if (mesh.processor_id() == 0) {
                out_dat_name = rname + "_dtXtime_RES.dat";
                foutDt.open(out_dat_name);
                out_dat_name = rname + "_massXtime_RES.dat";
                foutMass.open(out_dat_name);
            }
       
        }*/
    else if (ts_control_model_name == "PC11") {
        ts_control = new timeStepControlPC11(dt_init, dt_min, dt_max, nsa_max, tol_u, tol_s, pc11_theta, alpha, k_exp, s_min, s_max, complete_flow_norm);
        dt = dt_init;
        control_ts = true;
        if (mesh.processor_id() == 0) {
            out_dat_name = rname + "_dtXtime_PC11.dat";
            foutDt.open(out_dat_name);
            out_dat_name = rname + "_massXtime_PC11.dat";
            foutMass.open(out_dat_name);
        }
    } else {
        cout << "\n Fixed time-step adopted!\n";
        dt = dt_init;
        delete ts_control;
        if (mesh.processor_id() == 0) {
            out_dat_name = rname + "_dtXtime_CTE.dat";
            foutDt.open(out_dat_name);
            out_dat_name = rname + "_massXtime_CTE.dat";
            foutMass.open(out_dat_name);
        }
    }

    equation_systems.parameters.set<Real> ("dt") = dt;

    // assert for time-step control the maximum number of successive approximation (nsa_max) parameter
    if (ts_control) {
        if (nsa_max > flow_n_nonlinear_steps || nsa_max > transport_n_nonlinear_steps) {
            std::cout << "\nError: Maximum number of successive approximation greater than number NLI for Flow and/or Transport problem! ";
            std::cout << "Simulation will be abort.";
            return 0;
        }
    }

    // Steady-state control parameters
    unsigned int flow_sstate_count = 0;
    unsigned int transport_sstate_count = 0;

    // non-linear iteration counter
    unsigned int flow_nli_counter = 0;
    unsigned int transport_nli_counter = 0;
    bool diverged_flow = false, diverged_transport = false;

    string* current_files;
    perf_log.start_event("Write", "XDMF");
    current_files = xdmf_writer.write_time_step(equation_systems, time);
    perf_log.stop_event("Write", "XDMF");
    //cout << "[WRITE] " + current_files[0] + " - " + current_files[1] << endl;

#ifdef PROVENANCE
    perf_log.start_event("GetMaximumIterationsToFlow", "Provenance");
    provenance.outputGetMaximumIterationsToFlow(dt, tmax, n_time_steps, n_flow_nonlinear_iterations_total, flow_nonlinear_tolerance, n_flow_linear_iterations_total, current_files[1]);
    perf_log.stop_event("GetMaximumIterationsToFlow", "Provenance");
    perf_log.start_event("GetMaximumIterationsToTransport", "Provenance");
    provenance.outputGetMaximumIterationsToTransport(dt, tmax, n_time_steps, n_transport_nonlinear_iterations_total, transport_nonlinear_tolerance, n_transport_linear_iterations_total, current_files[1]);
    perf_log.stop_event("GetMaximumIterationsToTransport", "Provenance");
    // prov.outputFlowSolverConfig();
#endif

#ifdef USE_CATALYST   
    InSituCatalystInitCoprocessing(dim, numberOfScripts, extractionScript, visualizationScript,
            t_step, write_interval, current_files,
            equation_systems, perf_log, provenance);
#endif

    // STEP LOOP
    // Loop in time steps
#ifdef PROVENANCE
    provenance.resetTaskID();
    provenance.resetSubTaskID();
    provenance.resetIterationsFlow();
    provenance.resetIterationsTransport();
    provenance.resetIterationsMeshRefinement();
    provenance.resetIterationsComputeSolutionChange();
    provenance.resetIterationsComputeTimeStep();
    provenance.resetMeshDependencies();
    sediment_flow.attach_provenance(&provenance);
    sediment_transport.attach_provenance(&provenance);
#endif

    equation_systems.parameters.set<PerfLog*> ("PerfLog") = &perf_log;
    equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = max_linear_iter;
    equation_systems.parameters.set<double>("minimum_linear_solver_tolerance") = minimum_linear_solver_tol;
    equation_systems.parameters.set<unsigned int>("write_interval") = write_interval;

    cout<<"\nAdopting "<<((fem_model=="SUPG/PSPG")? "SUPG/PSPG": "RbVMS")<< " to solve Flow and Transport problems\n"<<endl;

    // Writing into a file the initial time-step value at the beginning of the simulation
    if (mesh.processor_id() == 0)
        foutDt << 0 << "  " << 0.0 << "  " << dt << endl;

    // Writing into a file the suspended mass normalized by the initial one at the simulation beginning
    sediment_transport.PrintMass(foutMass, t_step);

    for (t_step = init_tstep; (t_step < n_time_steps) && abs(time - tmax) > 1.0e-08 && (!diverged_flow && !diverged_transport); ++t_step) {

#ifdef PROVENANCE       
        provenance.incrementTaskID();
#endif

        if (is_file_exist("abort.run")) break;

        if (is_file_exist("reset.run")) {
            GetPot reset(input);
            dt = reset("time/deltat", 0.005);
            tmax = reset("time/tmax", 0.1);
            n_time_steps = reset("time/n_time_steps", 10);

            flow_n_nonlinear_steps = reset("flow_n_nonlinear_steps", 10);
            transport_n_nonlinear_steps = reset("transport_n_nonlinear_steps", 10);
            flow_nonlinear_tolerance = reset("flow_nonlinear_tolerance", 1.0E-03);
            transport_nonlinear_tolerance = reset("transport_nonlinear_tolerance", 1.0E-03);
            write_interval = reset("write_interval", 10);
            max_linear_iter = reset("max_linear_iterations", 2000);

            // AKI DEVO INCLUIR TODOS OS PARÂMETROS LIDOS NO RESET
            equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = max_linear_iter;
            equation_systems.parameters.set<unsigned int>("write_interval") = write_interval;
            equation_systems.parameters.set<Real> ("tmax") = tmax;
        }

        // storing a copy of the solution vectors to be used in case of rejecting current time-step
        UniquePtr<NumericVector<Number> > flow_soln_saved = flow_system.solution->clone();
        UniquePtr<NumericVector<Number> > transp_soln_saved = transport_system.solution->clone();

        // storing the last accepted solutions as "old solutions"
        *flow_system.old_local_solution = *flow_system.current_local_solution;
        *transport_system.old_local_solution = *transport_system.current_local_solution;

        bool TimeStepAccepted = true;
        do {
            time += dt;
            flow_system.time = time;
            transport_system.time = time;

            equation_systems.parameters.set<Real> ("time") = time;
            equation_systems.parameters.set<Real> ("dt") = dt;

            // A pretty update message
            std::cout << "\n" << std::setw(55)
                    << std::setfill('=')
                    << "\n";
            std::cout << " Time step ";
            {
                std::ostringstream out;
                out << std::setw(2)
                        << std::right
                        << t_step + 1
                        << "  simulation time = "
                        << std::fixed
                        << std::setw(6)
                        << std::setprecision(5)
                        << std::setfill('0')
                        << std::left
                        << time
                        << "...";
                std::cout << out.str() << std::endl;
            }
            std::cout << std::setw(55)
                    << std::setfill('=')
                    << "\n";

            // SolverSimulationFluid
            sediment_flow.solve(t_step, dt, time, 0, diverged_flow);
            // SolverSimulationTransport
            sediment_transport.solve(t_step, dt, time, 0, diverged_transport);

            n_flow_nonlinear_iterations_total += sediment_flow.nonlinear_iteractions();
            n_transport_nonlinear_iterations_total += sediment_transport.nonlinear_iteractions();
            n_flow_linear_iterations_total += sediment_flow.linear_iteractions();
            n_transport_linear_iterations_total += sediment_transport.linear_iteractions();

            if (ts_control) {
#ifdef PROVENANCE
                perf_log.start_event("ComputeSolutionChange", "Provenance");
                provenance.inputComputeSolutionChange();
                perf_log.stop_event("ComputeSolutionChange", "Provenance");
#endif

                perf_log.start_event("computeSolutionChange", "Time-Step Control");
                ts_control->computeSolutionChangeInTime(equation_systems);
                perf_log.stop_event("computeSolutionChange", "Time-Step Control");

                if (t_step >= ts_control->getStartTimeStepControl() - 1)
                    ts_control->checkTimeStepAcceptance(dt, sediment_flow.nonlinear_iteractions(), sediment_transport.nonlinear_iteractions(), TimeStepAccepted);
                else
                    ts_control->storeSolutionChangeinTime();

#ifdef PROVENANCE
                perf_log.start_event("ComputeSolutionChange", "Provenance");
                provenance.outputComputeSolutionChange(t_step, time, dt,
                        n_flow_linear_iterations_total, n_flow_nonlinear_iterations_total,
                        n_transport_linear_iterations_total, n_transport_nonlinear_iterations_total,
                        TimeStepAccepted, ts_control->getError());
                perf_log.stop_event("ComputeSolutionChange", "Provenance");
#endif
            }

            if (!TimeStepAccepted) {
                cout << "\nTimeStep was rejected! " << endl;
                time -= dt;
                flow_system.time = time;
                transport_system.time = time;

                cout << "Recovering previous solution " << endl;

                perf_log.start_event("RestoreSolution", "Time-Step Control");
                flow_system.solution->zero();
                flow_system.solution->add(*flow_soln_saved);
                flow_system.solution->close();

                transport_system.solution->zero();
                transport_system.solution->add(*transp_soln_saved);
                transport_system.solution->close();

                *flow_system.current_local_solution = *flow_system.old_local_solution;
                *transport_system.current_local_solution = *transport_system.old_local_solution;
                perf_log.stop_event("RestoreSolution", "Time-Step Control");

                // updating iteration counters
                n_rejected_flow_nonlinear_iterations_total += sediment_flow.nonlinear_iteractions();
                n_rejected_transport_nonlinear_iterations_total += sediment_transport.nonlinear_iteractions();
                n_rejected_flow_linear_iterations_total += sediment_flow.linear_iteractions();
                n_rejected_transport_linear_iterations_total += sediment_transport.linear_iteractions();
            }

            if (ts_control) {
#ifdef PROVENANCE
                perf_log.start_event("ComputeTimeStep", "Provenance");
                provenance.inputComputeTimeStep();
                perf_log.stop_event("ComputeTimeStep", "Provenance");
#endif
                
                perf_log.start_event("computeTimeStep", "Time-Step Control");
                if (t_step >= ts_control->getStartTimeStepControl() - 1 && (abs(equation_systems.parameters.get<Real> ("time") - tmax) > 1.0e-08 || !TimeStepAccepted))
                    ts_control->computeTimeStep(TimeStepAccepted, time, tmax, dt);
                perf_log.start_event("computeTimeStep", "Time-Step Control");

#ifdef PROVENANCE
                perf_log.start_event("ComputeTimeStep", "Provenance");
                provenance.outputComputeTimeStep(t_step, time, dt, TimeStepAccepted);
                perf_log.stop_event("ComputeTimeStep", "Provenance");
#endif
            }

        } while (!TimeStepAccepted);

        if (use_amr || first_step_refinement) {
            if ((t_step + 1) % ref_interval == 0 || first_step_refinement) {
#ifdef PROVENANCE
                perf_log.start_event("MeshRefinement", "Provenance");
                provenance.inputMeshRefinement();
                perf_log.stop_event("MeshRefinement", "Provenance");
#endif

                int beforeNActiveElem = mesh.n_active_elem();
                std::cout << "\n****************** Mesh Refinement *******************  " << std::endl;
                std::cout << " Considering Transport" << ((amrc_flow_transp && !first_step_refinement) ? " & Flow Variables\n" : " Variable\n");
                std::cout << "Number of elements before AMR step: " << beforeNActiveElem << std::endl;

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

                // After the initial refinement, refinements are going to be done only if max
                // number of refinements per time step (max_r_step) is greater than 1
                first_step_refinement = false;

#ifdef USE_CATALYST
                if (numberOfScripts > 0) {
                    FEAdaptor::mark_to_rebuild_grid();
                }
#endif 

                int AfterNActiveElem = mesh.n_active_elem();
                std::cout << "Number of elements after AMR step: " << AfterNActiveElem << std::endl;
                double var_nelem = double((AfterNActiveElem - beforeNActiveElem)) / beforeNActiveElem * 100;
                std::cout << "Nelem variation = " << var_nelem << " %." << endl;

                std::cout << "******************************************************\n" << std::endl;
#ifdef PROVENANCE
                perf_log.start_event("MeshRefinement", "Provenance");
                provenance.outputMeshRefinement(first_step_refinement, t_step, beforeNActiveElem, AfterNActiveElem);
                perf_log.stop_event("MeshRefinement", "Provenance");
#endif
                xdmf_writer.mesh_changed_on();

                // Solve Flow and Transport to new mesh
                sediment_flow.solve(t_step, dt, time, 1, diverged_flow);

                sediment_transport.solve(t_step, dt, time, 1, diverged_transport);

            }
        }

        sediment_deposition.ComputeDeposition();

#ifdef MESH_MOVIMENT
        std::cout << "Solving Mesh moviment..." << std::endl;

        //equation_systems.parameters.set<Real> ("linear solver tolerance") = initial_linear_solver_tol;
        mesh_system.solve();

        // How many iterations were required to solve the linear system?
        const unsigned int n_linear_iterations = flow_system.n_linear_iterations();
        std::cout << " Number of iterations: " << n_linear_iterations << std::endl;


        // What was the final residual of the linear system?
        const Real final_linear_residual = flow_system.final_linear_residual();
        std::cout << " final residual " << final_linear_residual << std::endl;

        moving_mesh.updateMesh();
#endif

        // Computing the suspended mass normalized by the initial one at current simulation time
        sediment_transport.PrintMass(foutMass, t_step + 1); // Output every write_interval time steps to file.       

        // Output every write_interval timesteps to file.
        if ((t_step + 1) % write_interval == 0) {
            const std::string mesh_restart = rname + "_mesh_restart.xdr";
            const std::string solution_restart = rname + "_solution_restart.xdr";

            mesh.write(mesh_restart);
            equation_systems.write(solution_restart, WRITE);

            std::ofstream f_restart;
            f_restart.open("restart.in");
            f_restart << "#Restart file: " << std::endl;
            f_restart << "time = " << time << std::endl;
            f_restart << "dt   = " << dt << std::endl;
            f_restart << "init_tstep = " << t_step << std::endl;
            f_restart << "mesh_restart = " << mesh_restart << std::endl;
            f_restart << "solution_restart = " << solution_restart << std::endl;
            f_restart << "xdmf_file_id = " << xdmf_writer.get_file_id() << std::endl;
            f_restart << "first_step_refinement = " << false << std::endl;
            f_restart << "mesh_restart = " << mesh_restart << std::endl;
            f_restart << "solution_restart = " << solution_restart << std::endl;
            f_restart << "mass_dep = " << sediment_transport.mass_dep << std::endl;
            f_restart << "initial_mass = " << sediment_transport.init_mass << std::endl;
            f_restart << "init_transport = " << false << std::endl;
            f_restart << "n_flow_nonlinear_iterations_total = " << n_flow_nonlinear_iterations_total << std::endl;
            f_restart << "n_transport_nonlinear_iterations_total = " << n_transport_nonlinear_iterations_total << std::endl;
            f_restart << "n_flow_linear_iterations_total    = " << n_flow_linear_iterations_total << std::endl;
            f_restart << "n_transport_linear_iterations_total    =  " << n_transport_linear_iterations_total << std::endl;
            f_restart << "n_rejected_flow_nonlinear_iterations_total = " << n_rejected_flow_nonlinear_iterations_total << std::endl;
            f_restart << "n_rejected_transport_nonlinear_iterations_total = " << n_rejected_transport_nonlinear_iterations_total << std::endl;
            f_restart << "n_rejected_flow_linear_iterations_total    =  " << n_rejected_flow_linear_iterations_total << std::endl;
            f_restart << "n_rejected_transport_linear_iterations_total    =  " << n_rejected_transport_linear_iterations_total << std::endl;

#ifdef PROVENANCE
            provenance.incrementSubTaskID();
            f_restart << "indexerID = " << provenance.getIndexerID() << std::endl;
            perf_log.start_event("MeshWriter", "Provenance");
            provenance.inputMeshWriter();
            perf_log.stop_event("MeshWriter", "Provenance");
#endif
            f_restart.close();

            perf_log.start_event("Write", "XDMF");
            current_files = xdmf_writer.write_time_step(equation_systems, time);
            perf_log.stop_event("Write", "XDMF");
            //cout << "[WRITE] " + current_files[0] + " - " + current_files[1] << endl;

#ifdef PROVENANCE
            perf_log.start_event("MeshWriter", "Provenance");
            provenance.outputMeshWriter(t_step, current_files[1]);
            perf_log.stop_event("MeshWriter", "Provenance");
#endif

            // Writing into a file the time-step size at current simulation time
            if (control_ts && abs(time - tmax) > 1.0e-08) {
                if (mesh.processor_id() == 0)
                    foutDt << t_step + 1 << "  " << time << "  " << ts_control->getLastAcceptedTS() << endl;
            } else {
                if (mesh.processor_id() == 0)
                    foutDt << t_step + 1 << "  " << time << "  " << dt << endl;
            }


        }

        // Prints solution variables range and L2 vector norm
        PrintStats(equation_systems);


#ifdef USE_CATALYST
        if (((t_step + 1) % catalyst_interval == 0)) {

            InSituCatalystCoprocessing(dim, numberOfScripts, extractionScript, visualizationScript,
                    t_step, time, write_interval, current_files,
                    equation_systems, perf_log, provenance);
        }
#endif

    } // end time step loop

    if (!diverged_flow && !diverged_transport) {
        if (t_step % write_interval != 0) {

#ifdef PROVENANCE
            provenance.incrementSubTaskID();
            perf_log.start_event("MeshWriter", "Provenance");
            provenance.inputMeshWriter();
            perf_log.stop_event("MeshWriter", "Provenance");
#endif

            perf_log.start_event("Write", "XDMF");
            current_files = xdmf_writer.write_time_step(equation_systems, time);
            perf_log.stop_event("Write", "XDMF");
            //cout << "[WRITE] " + current_files[0] + " - " + current_files[1] << endl;

#ifdef PROVENANCE
            perf_log.start_event("MeshWriter", "Provenance");
            provenance.outputMeshWriter(t_step, current_files[1]);
            perf_log.stop_event("MeshWriter", "Provenance");
#endif

            // Writing into a file time-step at current simulation time
            if (mesh.processor_id() == 0)
                foutDt << t_step << "  " << equation_systems.parameters.get<Real> ("time") << "  " << dt << endl;

        }

#ifdef USE_CATALYST
        if ((t_step) % catalyst_interval != 0) {
            InSituCatalystCoprocessing(dim, numberOfScripts, extractionScript, visualizationScript,
                    t_step, time, write_interval, current_files,
                    equation_systems, perf_log, provenance);
            FEAdaptor::Finalize();
        }
#endif

#ifdef PROVENANCE
        char out_filename[256];
        sprintf(out_filename, "%s_%d.xmf", rname.c_str(), libMesh::global_n_processors());

        perf_log.start_event("MeshAggregator", "Provenance");
        provenance.meshAggregator(out_filename, libMesh::global_n_processors());
        perf_log.stop_event("MeshAggregator", "Provenance");
        provenance.finishDataIngestor();
        solverPerformance.end();
#endif

        // Write time-step control performance
        if (control_ts)
            std::cout << *ts_control;

        // Write performance parameters to output file
        std::cout << "\n End of Simulation: " << ((flow_sstate_count > n_flow_sstate && transport_sstate_count > n_transport_sstate) ? "Steady-State reached!\n                   ---------------------" : "Final Simulation time reached!\n                   ----------------------------")
                << "\n Final number of active elements = " << mesh.n_active_elem() << endl;
        std::cout << "\nPERFORMANCE PARAMETERS"
                << "\n Total number of non-linear iterations: "
                << n_flow_nonlinear_iterations_total + n_transport_nonlinear_iterations_total
                << "\n Total number of non-linear iterations for Flow : "
                << n_flow_nonlinear_iterations_total
                << "\n Total number of non-linear iterations for Transport : "
                << n_transport_nonlinear_iterations_total
                << "\n Total number of linear iterations: "
                << n_flow_linear_iterations_total + n_transport_linear_iterations_total
                << "\n Total number of linear iterations for flow: "
                << n_flow_linear_iterations_total
                << "\n Total number of linear iterations for transport: "
                << n_transport_linear_iterations_total
                << "\n Total number of rejected non-linear iterations: "
                << n_rejected_flow_nonlinear_iterations_total + n_rejected_transport_nonlinear_iterations_total
                << "\n Total number of rejected non-linear iterations for flow: "
                << n_rejected_flow_nonlinear_iterations_total
                << "\n Total number of rejected non-linear iterations for transport: "
                << n_rejected_transport_nonlinear_iterations_total
                << "\n Total number of rejected linear iterations: "
                << n_rejected_flow_linear_iterations_total + n_rejected_transport_linear_iterations_total
                << "\n Total number of rejected linear iterations for flow: "
                << n_rejected_flow_linear_iterations_total
                << "\n Total number of rejected linear iterations for transport: "
                << n_rejected_transport_linear_iterations_total
                << "\n Rate Rejected/Total non-linear iterations: "
                << 100.0 * (n_rejected_flow_nonlinear_iterations_total + n_rejected_transport_nonlinear_iterations_total) / (n_flow_nonlinear_iterations_total + n_transport_nonlinear_iterations_total) << "%"
                << "\n Rate Rejected/Total non-linear iterations for Flow: "
                << 100.0 * n_rejected_flow_nonlinear_iterations_total / (n_flow_nonlinear_iterations_total) << "%"
                << "\n Rate Rejected/Total non-linear iterations for Transport: "
                << 100.0 * n_rejected_transport_nonlinear_iterations_total / (n_transport_nonlinear_iterations_total) << "%"
                << "\n Rate Rejected/Total linear iterations: "
                << 100.0 * (n_rejected_flow_linear_iterations_total + n_rejected_transport_linear_iterations_total) / (n_flow_linear_iterations_total + n_transport_linear_iterations_total) << "%"
                << "\n Rate Rejected/Total linear iterations for Flow: "
                << 100.0 * n_rejected_flow_linear_iterations_total / (n_flow_linear_iterations_total) << "%"
                << "\n Rate Rejected/Total linear iterations for Transport: "
                << 100.0 * n_rejected_transport_linear_iterations_total / (n_transport_linear_iterations_total) << "%"
                << std::endl;

        // All done.
        cout << "\nAll done!" << endl;
    } else {
        cout << "\nSimulation diverged. Aborting!\n";
    }

    // Closing "data against time" files
    foutDt.close();
    foutMass.close();

    return 0;
}

void PrintStats(EquationSystems & es) {

    std::vector<string> vname;
    std::vector<Real> local_min;
    std::vector<Real> local_max;
    std::vector<Real> local_l2norm;
    std::vector<unsigned int> vars;
    std::vector<unsigned int> sys;


    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();


    const System & flow_system = es.get_system("flow");
    const System & sediment_system = es.get_system("transport");

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


#ifdef USE_CATALYST

inline void InSituCatalystInitCoprocessing(int dim, int numberOfScripts,
        string &extractionScript, string& visualizationScript,
        int t_step, int write_interval, string* current_files,
        EquationSystems & equation_systems,
        PerfLog& perf_log, Provenance &provenance) {
    if (dim == 2 && numberOfScripts > 0 && !extractionScript.empty()) {
        // 2D analysis       
#ifdef PROVENANCE
        perf_log.start_event("InitDataExtraction", "Provenance");
        provenance.inputInitDataExtraction(-1);
        perf_log.stop_event("InitDataExtraction", "Provenance");
#endif
        perf_log.start_event("Init", "Catalyst");
        FEAdaptor::Initialize(numberOfScripts, extractionScript, visualizationScript);
        perf_log.stop_event("Init", "Catalyst");
        perf_log.start_event("CoProcess", "Catalyst");
        FEAdaptor::CoProcess(equation_systems, 0.0, t_step, write_interval, false, false);
        perf_log.stop_event("CoProcess", "Catalyst");
#ifdef PROVENANCE
        extractor::invoke2DRawDataExtractor(libMesh::global_processor_id(), t_step);
        provenance.incrementIndexerID();
        perf_log.start_event("InitDataExtraction", "Provenance");
        provenance.outputInitDataExtraction(-1, current_files[1], dim);
        perf_log.stop_event("InitDataExtraction", "Provenance");
#endif  
    } else if (dim == 3 && numberOfScripts > 0) {
        // 3D analysis
        for (int lineID = 0; lineID <= 3; lineID++) {
#ifdef PROVENANCE
            if (!visualizationScript.empty()) {
                perf_log.start_event("InitVisualization", "Provenance");
                provenance.inputInitVisualization(lineID);
                perf_log.stop_event("InitVisualization", "Provenance");
            }

            if (!extractionScript.empty()) {
                perf_log.start_event("InitDataExtraction", "Provenance");
                provenance.inputInitDataExtraction(lineID);
                perf_log.stop_event("InitDataExtraction", "Provenance");
            }
#endif
            if (lineID == 0) {
                perf_log.start_event("Init", "Catalyst");
                FEAdaptor::Initialize(numberOfScripts, extractionScript, visualizationScript);
                perf_log.stop_event("Init", "Catalyst");
                perf_log.start_event("CoProcess", "Catalyst");
                FEAdaptor::CoProcess(equation_systems, 0.0, t_step, write_interval, false, false);
                perf_log.stop_event("CoProcess", "Catalyst");
            }
#ifdef PROVENANCE
            if (!extractionScript.empty()) {
                extractor::invoke3DRawDataExtractor(libMesh::global_processor_id(), t_step, lineID);
            }

            provenance.incrementIndexerID();

            if (!visualizationScript.empty()) {
                perf_log.start_event("InitVisualization", "Provenance");
                provenance.outputInitVisualization(lineID, t_step);
                perf_log.stop_event("InitVisualization", "Provenance");
            }

            if (!extractionScript.empty()) {
                perf_log.start_event("InitDataExtraction", "Provenance");
                provenance.outputInitDataExtraction(lineID, current_files[1], dim);
                perf_log.stop_event("InitDataExtraction", "Provenance");
            }
#endif
        }
    }
}

inline void InSituCatalystCoprocessing(int dim, int numberOfScripts,
        string &extractionScript, string& visualizationScript,
        int t_step, double time, int write_interval, string* current_files,
        EquationSystems & equation_systems,
        PerfLog& perf_log, Provenance &provenance) {
    if (dim == 2 && numberOfScripts > 0 && !extractionScript.empty()) {
#ifdef PROVENANCE
        perf_log.start_event("DataExtraction", "Provenance");
        provenance.inputDataExtraction(-1);
        perf_log.stop_event("DataExtraction", "Provenance");
#endif
        perf_log.start_event("CoProcess", "Catalyst");
        FEAdaptor::CoProcess(equation_systems, time, t_step, write_interval, false, false);
        perf_log.stop_event("CoProcess", "Catalyst");
#ifdef PROVENANCE
        extractor::invoke2DRawDataExtractor(libMesh::global_processor_id(), t_step);
        provenance.incrementIndexerID();

        perf_log.start_event("DataExtraction", "Provenance");
        provenance.outputDataExtraction(-1, t_step, current_files[1], dim);
        perf_log.stop_event("DataExtraction", "Provenance");
#endif
    } else if (dim == 3 && numberOfScripts > 0) {
        // 3D analysis
        for (int lineID = 0; lineID <= 3; lineID++) {
#ifdef PROVENANCE
            if (!visualizationScript.empty()) {
                perf_log.start_event("Visualization", "Provenance");
                provenance.inputVisualization(lineID);
                perf_log.stop_event("Visualization", "Provenance");
            }

            if (!extractionScript.empty()) {
                perf_log.start_event("DataExtraction", "Provenance");
                provenance.inputDataExtraction(lineID);
                perf_log.stop_event("DataExtraction", "Provenance");
            }
#endif
            if (lineID == 0) {
                perf_log.start_event("CoProcess", "Catalyst");
                FEAdaptor::CoProcess(equation_systems, time, t_step, write_interval, false, false);
                perf_log.stop_event("CoProcess", "Catalyst");
            }
#ifdef PROVENANCE
            if (!extractionScript.empty()) {
                extractor::invoke3DRawDataExtractor(libMesh::global_processor_id(), t_step, lineID);
            }

            provenance.incrementIndexerID();

            if (!visualizationScript.empty()) {
                perf_log.start_event("Visualization", "Provenance");
                provenance.outputVisualization(lineID, t_step);
                perf_log.stop_event("Visualization", "Provenance");
            }

            if (!extractionScript.empty()) {
                perf_log.start_event("DataExtraction", "Provenance");
                provenance.outputDataExtraction(lineID, t_step, current_files[1], dim);
                perf_log.stop_event("DataExtraction", "Provenance");
            }
#endif
        }
    }
}
#endif