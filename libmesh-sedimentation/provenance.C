/* 
 * File:   provenance.C
 * Author: vitor
 *
 * Created on May 18, 2016, 09:42 AM
 */

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

#include "provenance.h"

using namespace std;
using namespace libMesh;

string space = "      ";
string directory = "/Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/example/prov";
string pgCommandLine = "java -jar /Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/dfa/PG-1.0.jar ";

void Provenance::inputMeshGeneration(int simulationID, int dim, int ncellx, int ncelly, int ncellz, 
			double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, int ref_interval)
{
	clock_t begin = clock();
	
	// run PG
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation meshGeneration -id " 
		+ to_string(simulationID) 
		+ " -workspace /Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/example/prov -status FINISHED";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-element -dataflow sedimentation -transformation meshGeneration -id "
	+  to_string(simulationID) 
	+ " -set imeshgeneration -element [{'" 
	+ to_string(simulationID) + ";"
	+ to_string(dim) + ";"
	+ to_string(ncellx) + ";"
	+ to_string(ncelly) + ";"
	+ to_string(ncellz) + ";"
	+ to_string(xmin) + ";"
	+ to_string(ymin) + ";"
	+ to_string(zmin) + ";"
	+ to_string(xmax) + ";"
	+ to_string(ymax) + ";"
	+ to_string(zmax) + ";"
	+ to_string(ref_interval)
	+ "'}]";
	system(strdup(str.c_str()));

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);

	ofstream file;
	file.open(directory + "/1.mesh-generation.prov", ios_base::app);
	file << "PROV:MeshGeneration:Input" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "dim(" + to_string(dim) + ")" << endl <<
	    space << "ncellx(" + to_string(ncellx) + ")" << endl <<
	    space << "ncelly(" + to_string(ncelly) + ")" << endl <<
	    space << "ncellz(" + to_string(ncellz) + ")" << endl <<
	    space << "xmin(" + to_string(xmin) + ")" << endl <<
	    space << "ymin(" + to_string(ymin) + ")" << endl <<
	    space << "zmin(" + to_string(zmin) + ")" << endl <<
	    space << "xmax(" + to_string(xmax) + ")" << endl <<
	    space << "ymax(" + to_string(ymax) + ")" << endl <<
	    space << "zmax(" + to_string(zmax) + ")" << endl <<
	    space << "ref_interval(" + to_string(ref_interval) + ")" << endl <<
	    space << "xmax(" + to_string(xmax) + ")" << endl;

	file << space << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
	file.close();
}

void Provenance::outputMeshGeneration(int simulationID, double r_fraction,double c_fraction,double max_h_level,unsigned int hlevels){
	clock_t begin = clock();

	// run PG
	// mesh generation
	string str = pgCommandLine + "-element -dataflow sedimentation -transformation meshGeneration -id "
	+  to_string(simulationID) 
	+ " -set omeshgeneration -element [{'" 
	+ to_string(simulationID) + ";"
	+ to_string(r_fraction) + ";"
	+ to_string(c_fraction) + ";"
	+ to_string(max_h_level) + ";"
	+ to_string(hlevels)
	+ "'}]";
	system(strdup(str.c_str()));

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation meshGeneration " + to_string(simulationID);
	system(strdup(str.c_str()));

	// mesh refinement
	str = pgCommandLine + "-task -dataflow sedimentation -transformation meshRefinement -id " 
		+ to_string(simulationID) 
		+ " -workspace /Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/example/prov -status FINISHED"
		+ " -dependencies [{meshGeneration},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);

	ofstream file;
	file.open(directory + "/2.mesh-refinement.prov", ios_base::app);
	file << "PROV:MeshRefinement:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "r_fraction(" + to_string(r_fraction) + ")" << endl <<
	    space << "c_fraction(" + to_string(c_fraction) + ")" << endl <<
	    space << "max_h_level(" + to_string(max_h_level) + ")" << endl <<
	    space << "hlevels(" + to_string(hlevels) + ")" << endl;

	file << space << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
	file.close();
}

void Provenance::outputMeshRefinement(int simulationID, bool first_step_refinement){
	clock_t begin = clock();

	// run PG
	// mesh refinement
	string str = pgCommandLine + "-element -dataflow sedimentation -transformation meshRefinement -id "
	+  to_string(simulationID) 
	+ " -set omeshrefinement -element [{'" 
	+ to_string(simulationID) + ";"
	+ to_string(first_step_refinement)
	+ "'}]";
	system(strdup(str.c_str()));

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation meshRefinement " + to_string(simulationID);
	system(strdup(str.c_str()));

	// create equation systems
	str = pgCommandLine + "-task -dataflow sedimentation -transformation createEquationSystems -id " 
		+ to_string(simulationID) 
		+ " -workspace /Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/example/prov -status FINISHED"
		+ " -dependencies [{meshRefinement},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);

	ofstream file;
	file.open(directory + "/1.mesh-refinement.prov", ios_base::app);
	file << "PROV:MeshRefinement:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "first_step_refinement(" + to_string(first_step_refinement) + ")" << endl;

	file << space << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
	file.close();
}

void Provenance::outputCreateEquationSystems(int simulationID, Real Reynolds,Real Gr,Real Sc,Real Us,Real Diffusivity,Real xlock,Real alfa,
	Real theta,Real ex,Real ey,Real ez,Real c_factor){
	clock_t begin = clock();

	// run PG
	// create equation systems
	string str = pgCommandLine + "-element -dataflow sedimentation -transformation createEquationSystems -id "
	+  to_string(simulationID) 
	+ " -set ocreateequationsystems -element [{'" 
	+ to_string(simulationID) + ";"
	+ to_string(Reynolds) + ";"
	+ to_string(Gr) + ";"
	+ to_string(Sc) + ";"
	+ to_string(Us) + ";"
	+ to_string(Diffusivity) + ";"
	+ to_string(xlock) + ";"
	+ to_string(alfa) + ";"
	+ to_string(theta) + ";"
	+ to_string(ex) + ";"
	+ to_string(ey) + ";"
	+ to_string(ez) + ";"
	+ to_string(c_factor)
	+ "'}]";
	system(strdup(str.c_str()));

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation createEquationSystems " + to_string(simulationID);
	system(strdup(str.c_str()));

	// get maximum iterations
	str = pgCommandLine + "-task -dataflow sedimentation -transformation getMaximumIterations -id " 
		+ to_string(simulationID) 
		+ " -workspace /Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/example/prov -status FINISHED"
		+ " -dependencies [{createEquationSystems},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);

	ofstream file;
	file.open(directory + "/3.create-equation-systems.prov", ios_base::app);
	file << "PROV:CreateEquationSystems:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "Reynolds(" + to_string(Reynolds) + ")" << endl <<
	    space << "Gr(" + to_string(Gr) + ")" << endl <<
	    space << "Sc(" + to_string(Sc) + ")" << endl <<
	    space << "Us(" + to_string(Us) + ")" << endl <<
	    space << "Diffusivity(" + to_string(Diffusivity) + ")" << endl <<
	    space << "xlock(" + to_string(xlock) + ")" << endl <<
	    space << "alfa(" + to_string(alfa) + ")" << endl <<
	    space << "theta(" + to_string(theta) + ")" << endl <<
	    space << "ex(" + to_string(ex) + ")" << endl <<
	    space << "ey(" + to_string(ey) + ")" << endl <<
	    space << "ez(" + to_string(ez) + ")" << endl <<
	    space << "c_factor(" + to_string(c_factor) + ")" << endl;

	file << space << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
	file.close();
}

void Provenance::outputGetMaximumIterations(int simulationID, Real dt, Real tmax, unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance, 
	int max_linear_iters, int max_r_steps, unsigned int write_interval){
	clock_t begin = clock();

	// pg
	// get maximum iterations
	string str = pgCommandLine + "-element -dataflow sedimentation -transformation getMaximumIterations -id "
	+  to_string(simulationID) 
	+ " -set ogetmaximumiterations -element [{'" 
	+ to_string(simulationID) + ";"
	+ to_string(dt) + ";"
	+ to_string(tmax) + ";"
	+ to_string(n_time_steps) + ";"
	+ to_string(n_nonlinear_steps) + ";"
	+ to_string(nonlinear_tolerance) + ";"
	+ to_string(max_linear_iters) + ";"
	+ to_string(max_r_steps) + ";"
	+ to_string(write_interval)
	+ "'}]";
	system(strdup(str.c_str()));

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation getMaximumIterations " + to_string(simulationID);
	system(strdup(str.c_str()));

	// solver simulation to the fluid
	// str = pgCommandLine + "-task -dataflow sedimentation -transformation solverSimulationFluid -id " 
	// 	+ to_string(simulationID) + " -subid " + to_string(n_linear_step)
	// 	+ " -workspace /Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/example/prov -invocation SolverSimulationFluid -status FINISHED";
	// system(strdup(str.c_str()));

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);

	ofstream file;
	file.open(directory + "/4.get-maximum-iterations.prov", ios_base::app);
	file << "PROV:GetMaximumIterations:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "dt(" + to_string(dt) + ")" << endl <<
	    space << "tmax(" + to_string(tmax) + ")" << endl <<
	    space << "n_time_steps(" + to_string(n_time_steps) + ")" << endl <<
	    space << "n_nonlinear_steps(" + to_string(n_nonlinear_steps) + ")" << endl <<
	    space << "nonlinear_tolerance(" + to_string(nonlinear_tolerance) + ")" << endl <<
	    space << "max_linear_iters(" + to_string(max_linear_iters) + ")" << endl <<
	    space << "max_r_steps(" + to_string(max_r_steps) + ")" << endl <<
	    space << "write_interval(" + to_string(write_interval) + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
	file.close();
}

void Provenance::outputSolverSimulationFluid(int simulationID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, 
	Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged){
	clock_t begin = clock();

	// pg
	// solver simulation to the fluid
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation solverSimulationFluid -id " 
		+ to_string(simulationID) + " -status FINISHED"
		+ " -workspace /Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/example/prov -invocation SolverSimulationFluid"
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{getMaximumIterations},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-element -dataflow sedimentation -transformation solverSimulationFluid -id "
	+  to_string(simulationID) 
	+ " -subid " + to_string(subTaskID)
	+ " -set osolversimulationfluid -element [{'" 
	+ to_string(simulationID) + ";"
	+ to_string(time_step) + ";"
	+ to_string(time) + ";"
	+ to_string(linear_step) + ";"
	+ to_string(n_linear_step) + ";"
	+ to_string(n_linear_iterations) + ";"
	+ to_string(linear_residual) + ";"
	+ to_string(norm_delta) + ";"
	+ to_string(norm_delta_u) + ";"
	+ to_string(converged)
	+ "'}]";
	system(strdup(str.c_str()));

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation solverSimulationFluid " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

	// solver simulation to the sediments
	// str = pgCommandLine + "-task -dataflow sedimentation -transformation solverSimulationSediments -id " 
	// 	+ to_string(simulationID) + " -subid " + to_string(n_linear_step)
	// 	+ " -workspace /Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/example/prov -invocation SolverSimulationSediments -status FINISHED";
	// system(strdup(str.c_str()));

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);

  	ofstream file;
	file.open(directory + "/5.solver-simulation-fluid.prov", ios_base::app);
	file << "PROV:SolverSimulationFluid:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "time_step(" + to_string(time_step) + ")" << endl <<
	    space << "time(" + to_string(time) + ")" << endl <<
	    space << "linear_step(" + to_string(linear_step) + ")" << endl <<
	    space << "n_linear_step(" + to_string(n_linear_step) + ")" << endl <<
	    space << "n_linear_iterations(" + to_string(n_linear_iterations) + ")" << endl <<
	    space << "linear_residual(" + to_string(linear_residual) + ")" << endl <<
	    space << "norm_delta(" + to_string(norm_delta) + ")" << endl <<
	    space << "norm_delta_u(" + to_string(norm_delta_u) + ")" << endl <<
	    space << "converged(" + to_string(converged) + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
  	file.close();

 //  	file.open(directory + "/debugging.prov", ios_base::app);
	// file << to_string(simulationID) << ";" 
	// 	<< to_string(time_step) << ";" 
	// 	<< to_string(time) << ";" 
	// 	<< to_string(linear_step) << ";" 
	// 	<< to_string(n_linear_step) << ";" 
	// 	<< to_string(n_linear_iterations) << ";" 
	//     << to_string(linear_residual) << ";" 
	//     << to_string(norm_delta) << ";" 
	//     << to_string(norm_delta_u) << ";" 
	//     << to_string(converged) << endl;
 //    file.close();
}

void Provenance::outputSolverSimulationSediments(int simulationID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, 
	Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged){
	clock_t begin = clock();

	// pg
	// solver simulation to the sediments
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation solverSimulationSediments -id " 
		+ to_string(simulationID) + " -status FINISHED"
		+ " -workspace /Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/example/prov -invocation SolverSimulationSediments"
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{solverSimulationFluid},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-element -dataflow sedimentation -transformation solverSimulationSediments -id "
	+  to_string(simulationID) 
	+ " -subid " + to_string(subTaskID)
	+ " -set osolversimulationsediments -element [{'" 
	+ to_string(simulationID) + ";"
	+ to_string(time_step) + ";"
	+ to_string(time) + ";"
	+ to_string(linear_step) + ";"
	+ to_string(n_linear_step) + ";"
	+ to_string(n_linear_iterations) + ";"
	+ to_string(linear_residual) + ";"
	+ to_string(norm_delta) + ";"
	+ to_string(norm_delta_u) + ";"
	+ to_string(converged)
	+ "'}]";
	system(strdup(str.c_str()));

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation solverSimulationSediments " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);

	ofstream file;
	file.open(directory + "/6.solver-simulation-sediments.prov", ios_base::app);
	file << "PROV:SolverSimulationSediments:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "time_step(" + to_string(time_step) + ")" << endl <<
	    space << "time(" + to_string(time) + ")" << endl <<
	    space << "linear_step(" + to_string(linear_step) + ")" << endl <<
	    space << "n_linear_step(" + to_string(n_linear_step) + ")" << endl <<
	    space << "n_linear_iterations(" + to_string(n_linear_iterations) + ")" << endl <<
	    space << "linear_residual(" + to_string(linear_residual) + ")" << endl <<
	    space << "norm_delta(" + to_string(norm_delta) + ")" << endl <<
	    space << "norm_delta_u(" + to_string(norm_delta_u) + ")" << endl <<
	    space << "converged(" + to_string(converged) + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
  	file.close();
}




