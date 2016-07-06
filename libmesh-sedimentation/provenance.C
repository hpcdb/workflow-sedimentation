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
#include <chrono>

#include "libmesh/libmesh.h"

#include "libmesh/getpot.h"

#include "provenance.h"
#include "performance.h"

using namespace std;
using namespace libMesh;

string space = "      ";
string directory = "";
string pgCommandLine = "";

Provenance::Provenance(){
	GetPot infile("provenance.in");
	directory = infile("directory", "/Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/sedimentation");
	string pgFilePath = infile("pgFilePath", "/Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/dfa/PG-1.0.jar");
	pgCommandLine = "java -jar " + pgFilePath + " ";
	processor_id = libMesh::global_processor_id();
}

void Provenance::inputMeshGeneration(int simulationID, int dim, int ncellx, int ncelly, int ncellz, 
			double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, int ref_interval)
{
	if(processor_id != 0) return; 

	Performance perf;
	perf.start();
        
    // run PG
    char buffer[4096];
	sprintf(buffer, "%s-task -dataflow sedimentation -transformation meshGeneration -id %d -workspace %s -status FINISHED",pgCommandLine.c_str(),simulationID,directory.c_str());
    system(strdup(buffer));

    sprintf(buffer, "%s-performance -starttime -dataflow sedimentation -transformation meshGeneration -task %d -computation libMeshSedimentation::MeshGeneration",pgCommandLine.c_str(),simulationID);
	system(strdup(buffer));

	sprintf(buffer, "%s-element -dataflow sedimentation -transformation meshGeneration -id %d -set imeshgeneration -element [{'%d;%d;%d;%d;%d;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%d'}]",pgCommandLine.c_str(),simulationID,simulationID,dim,ncellx,ncelly,ncellz,xmin,ymin,zmin,xmax,ymax,zmax,ref_interval);
	system(strdup(buffer));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	ofstream file;
	file.open("prov/log/mesh-generation.prov", ios_base::app);
	file << "PROV:MeshGeneration:Input" << endl;
	sprintf(buffer, "%ssimulationID(%d)", space.c_str(), simulationID);
	file << buffer << endl;
	sprintf(buffer, "%sdim(%d)", space.c_str(), dim);
	file << buffer << endl;
	sprintf(buffer, "%smcellx(%d)", space.c_str(), ncellx);
	file << buffer << endl;
	sprintf(buffer, "%sncelly(%d)", space.c_str(), ncelly);
	file << buffer << endl;
	sprintf(buffer, "%sncellz(%d)", space.c_str(), ncellz);
	file << buffer << endl;
	sprintf(buffer, "%sxmin(%.2f)", space.c_str(), xmin);
	file << buffer << endl;
	sprintf(buffer, "%symin(%.2f)", space.c_str(), ymin);
	file << buffer << endl;
	sprintf(buffer, "%szmin(%.2f)", space.c_str(), zmin);
	file << buffer << endl;
	sprintf(buffer, "%sxmax(%.2f)", space.c_str(), xmax);
	file << buffer << endl;
	sprintf(buffer, "%symax(%.2f)", space.c_str(), ymax);
	file << buffer << endl;
	sprintf(buffer, "%szmax(%.2f)", space.c_str(), zmax);
	file << buffer << endl;
	sprintf(buffer, "%sref_interval(%d)", space.c_str(), ref_interval);
	file << buffer << endl;
	sprintf(buffer, "%selapsed-time: %.2f seconds.", space.c_str(), elapsedTime);
	file << buffer << endl;
	file.close();
}

void Provenance::outputMeshGeneration(int simulationID, double r_fraction,double c_fraction,double max_h_level,unsigned int hlevels){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

	// run PG
	// mesh generation
	char buffer[4096];
	printf(buffer,"%s-element -dataflow sedimentation -transformation meshGeneration -id %d -set omeshgeneration -element [{'%d;%.2f;%.2f;%.2f;%d'}]",pgCommandLine.c_str(),simulationID,simulationID,r_fraction,c_fraction,max_h_level,hlevels);
	system(strdup(buffer));

	string str = pgCommandLine + "-performance -endtime -dataflow sedimentation -transformation meshGeneration -task "
	+  to_string(simulationID) 
	+ " -computation libMeshSedimentation::MeshGeneration";
	system(strdup(str.c_str()));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation meshGeneration " + to_string(simulationID);
	system(strdup(str.c_str()));

	perf.start();

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

	perf.end();
  	elapsedTime += perf.elapsedTime();

	ofstream file;
	file.open("prov/log/mesh-generation.prov", ios_base::app);
	file << "PROV:MeshGeneration:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "r_fraction(" + to_string(r_fraction) + ")" << endl <<
	    space << "c_fraction(" + to_string(c_fraction) + ")" << endl <<
	    space << "max_h_level(" + to_string(max_h_level) + ")" << endl <<
	    space << "hlevels(" + to_string(hlevels) + ")" << endl;

	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
	file.close();
}

void Provenance::outputCreateEquationSystems(int simulationID, Real Reynolds,Real Gr,Real Sc,Real Us,Real Diffusivity,Real xlock,Real fopc,
	Real theta,Real ex,Real ey,Real ez,Real c_factor){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

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
	+ to_string(fopc) + ";"
	+ to_string(theta) + ";"
	+ to_string(ex) + ";"
	+ to_string(ey) + ";"
	+ to_string(ez) + ";"
	+ to_string(c_factor)
	+ "'}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -endtime -dataflow sedimentation -transformation createEquationSystems -task "
	+  to_string(simulationID) 
	+ " -computation libMeshSedimentation::createEquationSystems"
	+ " -dependencies [{meshGeneration},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation createEquationSystems " + to_string(simulationID);
	system(strdup(str.c_str()));

	perf.start();

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

	perf.end();
  	elapsedTime += perf.elapsedTime();

	ofstream file;
	file.open("prov/log/create-equation-systems.prov", ios_base::app);
	file << "PROV:CreateEquationSystems:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "Reynolds(" + to_string(Reynolds) + ")" << endl <<
	    space << "Gr(" + to_string(Gr) + ")" << endl <<
	    space << "Sc(" + to_string(Sc) + ")" << endl <<
	    space << "Us(" + to_string(Us) + ")" << endl <<
	    space << "Diffusivity(" + to_string(Diffusivity) + ")" << endl <<
	    space << "xlock(" + to_string(xlock) + ")" << endl <<
	    space << "fopc(" + to_string(fopc) + ")" << endl <<
	    space << "theta(" + to_string(theta) + ")" << endl <<
	    space << "ex(" + to_string(ex) + ")" << endl <<
	    space << "ey(" + to_string(ey) + ")" << endl <<
	    space << "ez(" + to_string(ez) + ")" << endl <<
	    space << "c_factor(" + to_string(c_factor) + ")" << endl;

	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
	file.close();
}

void Provenance::outputGetMaximumIterations(int simulationID, Real dt, Real tmax, unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance, 
	int max_linear_iters, int max_r_steps, unsigned int write_interval, string xdmf){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

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
	+ directory + "/" + xdmf
	+ "'}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-file -dataflow sedimentation -transformation getMaximumIterations -id " +  to_string(simulationID) 
	+ " -name \"" + xdmf + "\" -path " + directory;
	system(strdup(str.c_str()));
	
	str = pgCommandLine + "-performance -endtime -dataflow sedimentation -transformation getMaximumIterations -task "
	+  to_string(simulationID) 
	+ " -computation libMeshSedimentation::GetMaximumIterations";
	system(strdup(str.c_str()));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation getMaximumIterations " + to_string(simulationID);
	system(strdup(str.c_str()));

	ofstream file;
	file.open("prov/log/get-maximum-iterations.prov", ios_base::app);
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
	    space << "xdmf(" + xdmf + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
	file.close();
}

void Provenance::inputInitDataExtraction(int simulationID, string transformation, string extractionFileName){
	if(processor_id != 0) return; 
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the sediments
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation " + transformation + " -id " 
		+ to_string(simulationID) + " -status RUNNING"
		+ " -workspace " + directory
		+ " -dependencies [{getMaximumIterations},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -starttime -dataflow sedimentation -transformation " + transformation + " -task "
	+  to_string(simulationID)
	+ " -computation libMeshSedimentation::" + transformation + "-" + to_string(simulationID);
	system(strdup(str.c_str()));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation " + transformation + " " + to_string(simulationID);
	system(strdup(str.c_str()));

	ofstream file;
	file.open("prov/log/" + extractionFileName + ".prov", ios_base::app);
	file << "PROV:" + transformation + ":Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
  	file.close();
}

void Provenance::outputInitDataExtraction(int simulationID, string transformation, string extractionFileName, string outDataSet, int time_step, string xdmf, string rawDataFile){
	if(processor_id != 0) return; 
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the sediments
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation " + transformation + " -id " 
		+ to_string(simulationID) + " -status FINISHED"
		+ " -workspace " + directory
		+ " -dependencies [{getMaximumIterations},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-file -dataflow sedimentation -transformation " + transformation + " -id " +  to_string(simulationID) 
	+ " -name \"" + xdmf + "\" -path " + directory;
	system(strdup(str.c_str()));

	str = pgCommandLine + "-file -dataflow sedimentation -transformation " + transformation + " -id " +  to_string(simulationID) 
	+ " -name \"" + rawDataFile + "\" -path " + directory;
	system(strdup(str.c_str()));

	str = pgCommandLine + "-element -dataflow sedimentation -transformation " + transformation + " -id "
	+  to_string(simulationID) 
	+ " -set " + outDataSet + " -element [{'" 
	+ to_string(simulationID) + ";"
	+ to_string(time_step) + ";"
	+ directory + "/" + xdmf + ";"
	+ directory + "/" + rawDataFile
	+ "'}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -endtime -dataflow sedimentation -transformation " + transformation + " -task "
	+  to_string(simulationID)
	+ " -computation libMeshSedimentation::" + transformation + "-" + to_string(simulationID);
	system(strdup(str.c_str()));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation " + transformation + " " + to_string(simulationID);
	system(strdup(str.c_str()));

	ofstream file;
	file.open("prov/log/" + extractionFileName + ".prov", ios_base::app);
	file << "PROV:" + transformation + ":Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "time_step(" + to_string(time_step) + ")" << endl <<
	    space << "xdmf(" + directory + "/" + xdmf + ")" << endl <<
	    space << "rawDataFile(" + directory + "/" + rawDataFile + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
  	file.close();
}

void Provenance::inputSolverSimulationFluid(int simulationID, int subTaskID){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the fluid
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation solverSimulationFluid -id " 
		+ to_string(simulationID) + " -status RUNNING"
		+ " -workspace " + directory + " -invocation SolverSimulationFluid"
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{getMaximumIterations},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -starttime -dataflow sedimentation -transformation solverSimulationFluid -task "
	+  to_string(simulationID) + " -subtask " + to_string(subTaskID)
	+ " -computation libMeshSedimentation::SolverSimulationFluid-" + to_string(simulationID) + "-" + to_string(subTaskID);
	system(strdup(str.c_str()));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation solverSimulationFluid " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

  	ofstream file;
	file.open("prov/log/solver-simulation-fluid.prov", ios_base::app);
	file << "PROV:SolverSimulationFluid:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "subTaskID(" + to_string(subTaskID) + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
  	file.close();
}

void Provenance::outputSolverSimulationFluid(int simulationID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, 
	Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the fluid
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation solverSimulationFluid -id " 
		+ to_string(simulationID) + " -status FINISHED"
		+ " -workspace " + directory
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{getMaximumIterations},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -endtime -dataflow sedimentation -transformation solverSimulationFluid -task "
	+  to_string(simulationID) + " -subtask " + to_string(subTaskID)
	+ " -computation libMeshSedimentation::SolverSimulationFluid-" + to_string(simulationID) + "-" + to_string(subTaskID);
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

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation solverSimulationFluid " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

  	ofstream file;
	file.open("prov/log/solver-simulation-fluid.prov", ios_base::app);
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

  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
  	file.close();
}

void Provenance::inputSolverSimulationSediments(int simulationID, int subTaskID){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the fluid
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation solverSimulationSediments -id " 
		+ to_string(simulationID) + " -status RUNNING"
		+ " -workspace " + directory + " -invocation SolverSimulationSediments"
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{solverSimulationFluid},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -starttime -dataflow sedimentation -transformation solverSimulationSediments -task "
	+  to_string(simulationID) + " -subtask " + to_string(subTaskID)
	+ " -computation libMeshSedimentation::SolverSimulationSediments-" + to_string(simulationID) + "-" + to_string(subTaskID);
	system(strdup(str.c_str()));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation solverSimulationSediments " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

  	ofstream file;
	file.open("prov/log/solver-simulation-sediments.prov", ios_base::app);
	file << "PROV:SolverSimulationSediments:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "subTaskID(" + to_string(subTaskID) + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
  	file.close();
}

void Provenance::outputSolverSimulationSediments(int simulationID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, 
	Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the sediments
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation solverSimulationSediments -id " 
		+ to_string(simulationID) + " -status FINISHED"
		+ " -workspace " + directory
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{solverSimulationFluid},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -endtime -dataflow sedimentation -transformation solverSimulationSediments -task "
	+  to_string(simulationID) + " -subtask " + to_string(subTaskID)
	+ " -computation libMeshSedimentation::SolverSimulationSediments-" + to_string(simulationID) + "-" + to_string(subTaskID);
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

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation solverSimulationSediments " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

	ofstream file;
	file.open("prov/log/solver-simulation-sediments.prov", ios_base::app);
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

  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
  	file.close();
}

void Provenance::outputMeshRefinement(int simulationID, int subTaskID, bool first_step_refinement, int time_step, int before_n_active_elem, int after_n_active_elem){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

	// run PG
	// mesh refinement
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation meshRefinement -id " 
		+ to_string(simulationID) 
		+ " -workspace " + directory + " -status FINISHED"
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{solverSimulationSediments},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -starttime -dataflow sedimentation -transformation meshRefinement -task "
	+  to_string(simulationID) + " -subtask " + to_string(subTaskID)
	+ " -computation libMeshSedimentation::MeshRefinement-" + to_string(simulationID) + "-" + to_string(subTaskID);
	system(strdup(str.c_str()));

	// input element
	str = pgCommandLine + "-element -dataflow sedimentation -transformation meshRefinement -id "
	+  to_string(simulationID) 
	+ " -subid " + to_string(subTaskID)
	+ " -set omeshrefinement -element [{'" 
	+ to_string(simulationID) + ";"
	+ to_string(first_step_refinement) + ";"
	+ to_string(time_step) + ";"
	+ to_string(before_n_active_elem) + ";"
	+ to_string(after_n_active_elem)
	+ "'}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -endtime -dataflow sedimentation -transformation meshRefinement -task "
	+  to_string(simulationID) + " -subtask " + to_string(subTaskID)
	+ " -computation libMeshSedimentation::MeshRefinement-" + to_string(simulationID) + "-" + to_string(subTaskID);
	system(strdup(str.c_str()));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation meshRefinement " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

	ofstream file;
	file.open("prov/log/mesh-refinement.prov", ios_base::app);
	file << "PROV:MeshRefinement:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "subTaskID(" + to_string(subTaskID) + ")" << endl <<
	    space << "first_step_refinement(" + to_string(first_step_refinement) + ")" << endl << 
	    space << "time_step(" + to_string(time_step) + ")" << endl << 
	    space << "before_n_active_elem(" + to_string(before_n_active_elem) + ")" << endl <<
	    space << "after_n_active_elem(" + to_string(after_n_active_elem) + ")" << endl;

	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
	file.close();
}

void Provenance::inputMeshWriter(int simulationID, int subTaskID){
	if(processor_id != 0) return; 
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the fluid
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation meshWriter -id " 
		+ to_string(simulationID) + " -status RUNNING"
		+ " -workspace " + directory + " -invocation MeshWriter"
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{solverSimulationSediments},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -starttime -dataflow sedimentation -transformation meshWriter -task "
	+  to_string(simulationID) + " -subtask " + to_string(subTaskID)
	+ " -computation libMeshSedimentation::MeshWriter-" + to_string(simulationID) + "-" + to_string(subTaskID);
	system(strdup(str.c_str()));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation meshWriter " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

  	ofstream file;
	file.open("prov/log/mesh-writer.prov", ios_base::app);
	file << "PROV:MeshWriter:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "subTaskID(" + to_string(subTaskID) + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
  	file.close();
  	return;
}

void Provenance::outputMeshWriter(int simulationID, int subTaskID, int time_step, string xdmf){
	if(processor_id != 0) return; 
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the sediments
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation meshWriter -id " 
		+ to_string(simulationID) + " -status FINISHED"
		+ " -workspace " + directory
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{solverSimulationSediments},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -endtime -dataflow sedimentation -transformation meshWriter -task "
	+  to_string(simulationID) + " -subtask " + to_string(subTaskID)
	+ " -computation libMeshSedimentation::MeshWriter-" + to_string(simulationID) + "-" + to_string(subTaskID);
	system(strdup(str.c_str()));

	str = pgCommandLine + "-file -dataflow sedimentation -transformation meshWriter -id " +  to_string(simulationID) 
	+ " -subid " + to_string(subTaskID)
	+ " -name \"" + xdmf + "\" -path " + directory;
	system(strdup(str.c_str()));

	str = pgCommandLine + "-element -dataflow sedimentation -transformation meshWriter -id "
	+  to_string(simulationID) 
	+ " -subid " + to_string(subTaskID)
	+ " -set omeshwriter -element [{'" 
	+ to_string(simulationID) + ";"
	+ to_string(time_step) + ";"
	+ directory + "/" + xdmf
	+ "'}]";
	system(strdup(str.c_str()));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation meshWriter " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

	ofstream file;
	file.open("prov/log/mesh-writer.prov", ios_base::app);
	file << "PROV:MeshWriter:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "subTaskID(" + to_string(subTaskID) + ")" << endl <<
	    space << "time_step(" + to_string(time_step) + ")" << endl <<
	    space << "xdmf(" + directory + "/" + xdmf + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
  	file.close();
}

void Provenance::inputDataExtraction(int simulationID, int subTaskID, string transformation, string extractionFileName){
	if(processor_id != 0) return; 
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the sediments
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation " + transformation + " -id " 
		+ to_string(simulationID) + " -status RUNNING"
		+ " -workspace " + directory
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{meshWriter},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -starttime -dataflow sedimentation -transformation " + transformation + " -task "
	+  to_string(simulationID) + " -subtask " + to_string(subTaskID)
	+ " -computation libMeshSedimentation::" + transformation + "-" + to_string(simulationID);
	system(strdup(str.c_str()));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation " + transformation + " " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

	ofstream file;
	file.open("prov/log/" + extractionFileName + ".prov", ios_base::app);
	file << "PROV:" + transformation + ":Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl << 
	    space << "subTaskID(" + to_string(subTaskID) + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
  	file.close();
}

void Provenance::outputDataExtraction(int simulationID, int subTaskID, string transformation, string extractionFileName, string outDataSet, int time_step, string xdmf, string rawDataFile){
	if(processor_id != 0) return; 
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the sediments
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation " + transformation + " -id " 
		+ to_string(simulationID) + " -status FINISHED"
		+ " -workspace " + directory
		+ " -subid " + to_string(subTaskID)
		+ " -dependencies [{meshWriter},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -endtime -dataflow sedimentation -transformation " + transformation + " -task "
	+  to_string(simulationID) + " -subtask " + to_string(subTaskID)
	+ " -computation libMeshSedimentation::" + transformation + "-" + to_string(simulationID);
	system(strdup(str.c_str()));

	str = pgCommandLine + "-file -dataflow sedimentation -transformation " + transformation + " -id " +  to_string(simulationID) 
	+ " -subid " + to_string(subTaskID)
	+ " -name \"" + xdmf + "\" -path " + directory;
	system(strdup(str.c_str()));

	str = pgCommandLine + "-file -dataflow sedimentation -transformation " + transformation + " -id " +  to_string(simulationID) 
	+ " -subid " + to_string(subTaskID)
	+ " -name \"" + rawDataFile + "\" -path " + directory;
	system(strdup(str.c_str()));

	str = pgCommandLine + "-element -dataflow sedimentation -transformation " + transformation + " -id "
	+  to_string(simulationID) 
	+ " -subid " + to_string(subTaskID)
	+ " -set " + outDataSet + " -element [{'" 
	+ to_string(simulationID) + ";"
	+ to_string(time_step) + ";"
	+ directory + "/" + xdmf + ";"
	+ directory + "/" + rawDataFile
	+ "'}]";
	system(strdup(str.c_str()));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation " + transformation + " " + to_string(simulationID) + " " + to_string(subTaskID);
	system(strdup(str.c_str()));

	ofstream file;
	file.open("prov/log/" + extractionFileName + ".prov", ios_base::app);
	file << "PROV:" + transformation + ":Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "subTaskID(" + to_string(subTaskID) + ")" << endl <<
	    space << "time_step(" + to_string(time_step) + ")" << endl <<
	    space << "xdmf(" + directory + "/" + xdmf + ")" << endl <<
	    space << "rawDataFile(" + directory + "/" + rawDataFile + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
  	file.close();
}

void Provenance::meshAggregator(int simulationID, string xdmf, int n_processors){
	if(processor_id != 0) return; 
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the sediments
	string str = pgCommandLine + "-task -dataflow sedimentation -transformation meshAggregator -id " 
		+ to_string(simulationID) + " -status FINISHED"
		+ " -workspace " + directory
		+ " -dependencies [{meshWriter},{" + to_string(simulationID) + "}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -starttime -dataflow sedimentation -transformation meshAggregator -task "
	+  to_string(simulationID)
	+ " -computation libMeshSedimentation::MeshAggregator-" + to_string(simulationID);
	system(strdup(str.c_str()));

	str = pgCommandLine + "-file -dataflow sedimentation -transformation meshAggregator -id " +  to_string(simulationID) 
	+ " -name \"" + xdmf + "\" -path " + directory;
	system(strdup(str.c_str()));

	str = pgCommandLine + "-element -dataflow sedimentation -transformation meshAggregator -id "
	+  to_string(simulationID) 
	+ " -set omeshaggregator -element [{'" 
	+ to_string(simulationID) + ";"
	+ directory + "/" + xdmf + ";"
	+ to_string(n_processors)
	+ "'}]";
	system(strdup(str.c_str()));

	str = pgCommandLine + "-performance -endtime -dataflow sedimentation -transformation meshAggregator -task "
	+  to_string(simulationID)
	+ " -computation libMeshSedimentation::MeshAggregator-" + to_string(simulationID);
	system(strdup(str.c_str()));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	str = pgCommandLine + "-ingest -task sedimentation meshAggregator " + to_string(simulationID);
	system(strdup(str.c_str()));

	ofstream file;
	file.open("prov/log/mesh-aggregator.prov", ios_base::app);
	file << "PROV:MeshAggregator:Output" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "xdmf(" + directory + "/" + xdmf + ")" << endl <<
	    space << "n_processors(" + to_string(n_processors) + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
  	file.close();
}

void Provenance::storeDataExtractionCost(int simulationID, int subTaskID, int time_step, string xdmf, string rawDataFile, double elapsedTime){
	if(processor_id != 0) return; 
	ofstream file;
	file.open("prov/rde/data-extraction.prov", ios_base::app);
	file << "RDE:DataExtraction:Process" << endl << 
	    space << "simulationID(" + to_string(simulationID) + ")" << endl <<
	    space << "subTaskID(" + to_string(subTaskID) + ")" << endl <<
	    space << "time_step(" + to_string(time_step) + ")" << endl <<
	    space << "xdmf(" + directory + "/" + xdmf + ")" << endl <<
	    space << "rawDataFile(" + directory + "/" + rawDataFile + ")" << endl;

  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
  	file.close();
}

void Provenance::storeSolverCost(double elapsedTime){
	if(processor_id != 0) return; 
	ofstream file;
	file.open("prov/solver/time.prov", ios_base::app);
	file << "Solver:Time:Process" << endl;
  	file << space << "elapsed-time: " << to_string(elapsedTime) << " seconds." << endl;
  	file.close();
}

void Provenance::finishDataIngestor(){
	if(processor_id != 0) return; 
	
	string str = "cp ../dfa/finish.token prov/di/sedimentation";
	system(strdup(str.c_str()));

	cout << "[Provenance] Finish Data Ingestor" << endl;
}



