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

#include "libmesh/getpot.h"

#include "provenance.h"

using namespace std;
using namespace libMesh;

string space = "      ";
string directory = "";
string pgCommandLine = "";

Provenance::Provenance(){
	GetPot infile("provenance.in");
	directory = infile("directory", "/Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/example/prov");
	string pgFilePath = infile("pgFilePath", "/Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/dfa/PG-1.0.jar");
	pgCommandLine = "java -jar " + pgFilePath + " ";
	processor_id = libMesh::global_processor_id();
}

void Provenance::inputMeshGeneration(int simulationID, int dim, int ncellx, int ncelly, int ncellz, 
			double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, int ref_interval)
{
	if(processor_id != 0) return; 

	clock_t begin = clock();
	
	// run PG
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation meshGeneration -id " 
		+ to_string(simulationID) 
		+ " -workspace " + directory + " -status FINISHED";
	system(strdup(str.c_str()));
	str = pgCommandLine + "-performance -starttime -dataflow sedimentation -transformation meshGeneration -task "
	+  to_string(simulationID) 
	+ " -computation libMeshSedimentation::MeshGeneration";
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
	file.open("prov/log/1.mesh-generation.prov", ios_base::app);
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
	    space << "ref_interval(" + to_string(ref_interval) + ")" << endl;

	file << space << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
	file.close();
}

void Provenance::outputMeshGeneration(int simulationID, double r_fraction,double c_fraction,double max_h_level,unsigned int hlevels){
	if(processor_id != 0) return; 
	
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

	str = pgCommandLine + "-performance -endtime -dataflow sedimentation -transformation meshGeneration -task "
	+  to_string(simulationID) ;
	system(strdup(str.c_str()));

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation meshGeneration " + to_string(simulationID);
	system(strdup(str.c_str()));

	// create equation systems
	str = pgCommandLine + "-task -dataflow sedimentation -transformation createEquationSystems -id " 
		+ to_string(simulationID) 
		+ " -workspace " + directory + " -status FINISHED"
		+ " -dependencies [{meshGeneration},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));
	
	str = pgCommandLine + "-performance -starttime -dataflow sedimentation -transformation createEquationSystems -task "
	+  to_string(simulationID) 
	+ " -computation libMeshSedimentation::CreateEquationSystems";
	system(strdup(str.c_str()));

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);

	ofstream file;
	file.open("prov/log/1.mesh-generation.prov", ios_base::app);
	file << "PROV:MeshGeneration:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "r_fraction(" + to_string(r_fraction) + ")" << endl <<
	    space << "c_fraction(" + to_string(c_fraction) + ")" << endl <<
	    space << "max_h_level(" + to_string(max_h_level) + ")" << endl <<
	    space << "hlevels(" + to_string(hlevels) + ")" << endl;

	file << space << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
	file.close();
}

void Provenance::outputCreateEquationSystems(int simulationID, Real Reynolds,Real Gr,Real Sc,Real Us,Real Diffusivity,Real xlock,Real alfa,
	Real theta,Real ex,Real ey,Real ez,Real c_factor){
	if(processor_id != 0) return; 
	
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

	str = pgCommandLine + "-performance -endtime -dataflow sedimentation -transformation createEquationSystems -task "
	+  to_string(simulationID) ;
	system(strdup(str.c_str()));

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation createEquationSystems " + to_string(simulationID);
	system(strdup(str.c_str()));

	// get maximum iterations
	str = pgCommandLine + "-task -dataflow sedimentation -transformation getMaximumIterations -id " 
		+ to_string(simulationID) 
		+ " -workspace " + directory + " -status FINISHED"
		+ " -dependencies [{createEquationSystems},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -starttime -dataflow sedimentation -transformation getMaximumIterations -task "
	+  to_string(simulationID) 
	+ " -computation libMeshSedimentation::GetMaximumIterations";
	system(strdup(str.c_str()));

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);

	ofstream file;
	file.open("prov/log/2.create-equation-systems.prov", ios_base::app);
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
	int max_linear_iters, int max_r_steps, unsigned int write_interval, string hdf5, string xdmf){
	if(processor_id != 0) return; 
	
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
	+ to_string(write_interval) + ";"
	+ hdf5 + ";"
	+ xdmf
	+ "'}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -endtime -dataflow sedimentation -transformation getMaximumIterations -task "
	+  to_string(simulationID) ;
	system(strdup(str.c_str()));

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation getMaximumIterations " + to_string(simulationID);
	system(strdup(str.c_str()));

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);

	ofstream file;
	file.open("prov/log/3.get-maximum-iterations.prov", ios_base::app);
	file << "PROV:GetMaximumIterations:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "dt(" + to_string(dt) + ")" << endl <<
	    space << "tmax(" + to_string(tmax) + ")" << endl <<
	    space << "n_time_steps(" + to_string(n_time_steps) + ")" << endl <<
	    space << "n_nonlinear_steps(" + to_string(n_nonlinear_steps) + ")" << endl <<
	    space << "nonlinear_tolerance(" + to_string(nonlinear_tolerance) + ")" << endl <<
	    space << "max_linear_iters(" + to_string(max_linear_iters) + ")" << endl <<
	    space << "max_r_steps(" + to_string(max_r_steps) + ")" << endl <<
	    space << "write_interval(" + to_string(write_interval) + ")" << endl <<
	    space << "hdf5(" + hdf5 + ")" << endl <<
	    space << "xdmf(" + xdmf + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
	file.close();
}

void Provenance::inputSolverSimulationFluid(int simulationID, int subTaskID){
	if(processor_id != 0) return; 
	
	clock_t begin = clock();

	// pg
	// solver simulation to the fluid
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation solverSimulationFluid -id " 
		+ to_string(simulationID) + " -status RUNNING"
		+ " -workspace " + directory + " -invocation SolverSimulationFluid"
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{getMaximumIterations},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -starttime -dataflow sedimentation -transformation solverSimulationFluid -task "
	+  to_string(simulationID) + " -subid " + to_string(subTaskID)
	+ " -computation libMeshSedimentation::SolverSimulationFluid-" + to_string(simulationID) + "-" + to_string(subTaskID);
	system(strdup(str.c_str()));

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation solverSimulationFluid " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);

  	ofstream file;
	file.open("prov/log/4.solver-simulation-fluid.prov", ios_base::app);
	file << "PROV:SolverSimulationFluid:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "subTaskID(" + to_string(subTaskID) + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
  	file.close();
}

void Provenance::outputSolverSimulationFluid(int simulationID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, 
	Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged){
	if(processor_id != 0) return; 
	
	clock_t begin = clock();

	// pg
	// solver simulation to the fluid
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation solverSimulationFluid -id " 
		+ to_string(simulationID) + " -status FINISHED"
		+ " -subid " + to_string(subTaskID);
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
	file.open("prov/log/4.solver-simulation-fluid.prov", ios_base::app);
	file << "PROV:SolverSimulationFluid:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "subTaskID(" + to_string(subTaskID) + ")" << endl <<
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

void Provenance::inputSolverSimulationSediments(int simulationID, int subTaskID){
	if(processor_id != 0) return; 
	
	clock_t begin = clock();

	// pg
	// solver simulation to the fluid
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation solverSimulationSediments -id " 
		+ to_string(simulationID) + " -status RUNNING"
		+ " -workspace " + directory + " -invocation SolverSimulationSediments"
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{solverSimulationFluid},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -starttime -dataflow sedimentation -transformation solverSimulationSediments -task "
	+  to_string(simulationID) + " -subid " + to_string(subTaskID)
	+ " -computation libMeshSedimentation::SolverSimulationSediments-" + to_string(simulationID) + "-" + to_string(subTaskID);
	system(strdup(str.c_str()));

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation solverSimulationSediments " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);

  	ofstream file;
	file.open("prov/log/5.solver-simulation-sediments.prov", ios_base::app);
	file << "PROV:SolverSimulationSediments:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "subTaskID(" + to_string(subTaskID) + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
  	file.close();
}

void Provenance::outputSolverSimulationSediments(int simulationID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, 
	Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged){
	if(processor_id != 0) return; 
	
	clock_t begin = clock();

	// pg
	// solver simulation to the sediments
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation solverSimulationSediments -id " 
		+ to_string(simulationID) + " -status FINISHED"
		+ " -subid " + to_string(subTaskID);
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
	file.open("prov/log/5.solver-simulation-sediments.prov", ios_base::app);
	file << "PROV:SolverSimulationSediments:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "subTaskID(" + to_string(subTaskID) + ")" << endl <<
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

void Provenance::outputMeshRefinement(int simulationID, int subTaskID, bool first_step_refinement, int before_n_active_elem, int after_n_active_elem){
	if(processor_id != 0) return; 
	
	clock_t begin = clock();

	// run PG
	// mesh refinement
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation meshRefinement -id " 
		+ to_string(simulationID) 
		+ " -workspace " + directory + " -status FINISHED"
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{solverSimulationSediments},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	// input element
	str = pgCommandLine + "-element -dataflow sedimentation -transformation meshRefinement -id "
	+  to_string(simulationID) 
	+ " -subid " + to_string(subTaskID)
	+ " -set omeshrefinement -element [{'" 
	+ to_string(simulationID) + ";"
	+ to_string(first_step_refinement) + ";"
	+ to_string(before_n_active_elem) + ";"
	+ to_string(after_n_active_elem)
	+ "'}]";
	system(strdup(str.c_str()));

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation meshRefinement " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);

	ofstream file;
	file.open("prov/log/6.mesh-refinement.prov", ios_base::app);
	file << "PROV:MeshRefinement:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "subTaskID(" + to_string(subTaskID) + ")" << endl <<
	    space << "first_step_refinement(" + to_string(first_step_refinement) + ")" << endl << 
	    space << "before_n_active_elem(" + to_string(before_n_active_elem) + ")" << endl <<
	    space << "after_n_active_elem(" + to_string(after_n_active_elem) + ")" << endl;

	file << space << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
	file.close();
}

void Provenance::finishDataIngestor(){
	if(processor_id != 0) return; 
	
	string str = "cp ../dfa/finish.token prov/di/sedimentation";
	system(strdup(str.c_str()));

	cout << "[Provenance] Finish Data Ingestor" << endl;
}



