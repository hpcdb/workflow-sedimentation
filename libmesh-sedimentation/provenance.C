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
}

void Provenance::outputMeshGeneration(int simulationID, double r_fraction,double c_fraction,double max_h_level,unsigned int hlevels){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

	// run PG
	// mesh generation
	char buffer[4096];
	sprintf(buffer,"%s-element -dataflow sedimentation -transformation meshGeneration -id %d -set omeshgeneration -element [{'%d;%.2f;%.2f;%.2f;%d'}]",pgCommandLine.c_str(),simulationID,simulationID,r_fraction,c_fraction,max_h_level,hlevels);
	system(strdup(buffer));

	sprintf(buffer,"%s-performance -endtime -dataflow sedimentation -transformation meshGeneration -task %d -computation libMeshSedimentation::MeshGeneration",pgCommandLine.c_str(),simulationID);
	system(strdup(buffer));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	sprintf(buffer,"%s-ingest -task sedimentation meshGeneration %d",pgCommandLine.c_str(),simulationID);
	system(strdup(buffer));

	perf.start();

	// create equation systems
	sprintf(buffer,"%s-task -dataflow sedimentation -transformation createEquationSystems -id %d -workspace %s -status FINISHED -dependencies [{meshGeneration},{%d}]",pgCommandLine.c_str(),simulationID,directory.c_str(),simulationID);
	system(strdup(buffer));
	
	sprintf(buffer,"%s-performance -starttime -dataflow sedimentation -transformation createEquationSystems -task %d -computation libMeshSedimentation::CreateEquationSystems",pgCommandLine.c_str(),simulationID);
	system(strdup(buffer));

	perf.end();
  	elapsedTime += perf.elapsedTime();
}

void Provenance::outputCreateEquationSystems(int simulationID, Real Reynolds,Real Gr,Real Sc,Real Us,Real Diffusivity,Real xlock,Real fopc,
	Real theta,Real ex,Real ey,Real ez,Real c_factor){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

	// run PG
	// create equation systems
	char buffer[4096];
	sprintf(buffer,"%s-element -dataflow sedimentation -transformation createEquationSystems -id %d -set ocreateequationsystems -element [{'%d;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f'}]",pgCommandLine.c_str(),simulationID,simulationID,Reynolds,Gr,Sc,Us,Diffusivity,xlock,fopc,theta,ex,ey,ez,c_factor);
	system(strdup(buffer));

	sprintf(buffer,"%s-performance -endtime -dataflow sedimentation -transformation createEquationSystems -task %d -computation libMeshSedimentation::createEquationSystems -dependencies [{meshGeneration},{%d}]",pgCommandLine.c_str(),simulationID,simulationID);
	system(strdup(buffer));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	sprintf(buffer,"%s-ingest -task sedimentation createEquationSystems %d",pgCommandLine.c_str(),simulationID);
	system(strdup(buffer));

	perf.start();

	// get maximum iterations
	sprintf(buffer,"%s-task -dataflow sedimentation -transformation getMaximumIterations -id %d -workspace %s -status FINISHED -dependencies [{createEquationSystems},{%d}]",pgCommandLine.c_str(),simulationID,directory.c_str(),simulationID);
	system(strdup(buffer));

	sprintf(buffer,"%s-performance -starttime -dataflow sedimentation -transformation getMaximumIterations -task %d -computation libMeshSedimentation::GetMaximumIterations",pgCommandLine.c_str(),simulationID);
	system(strdup(buffer));

	perf.end();
  	elapsedTime += perf.elapsedTime();
}

void Provenance::outputGetMaximumIterations(int simulationID, Real dt, Real tmax, unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance, 
	int max_linear_iters, int max_r_steps, unsigned int write_interval, string xdmf){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

	// pg
	// get maximum iterations
	char buffer[4096];
	sprintf(buffer,"%s-element -dataflow sedimentation -transformation getMaximumIterations -id %d -set ogetmaximumiterations -element [{'%d;%.2f;%.2f;%d;%d;%.2f;%d;%d;%d;%s/%s'}]",pgCommandLine.c_str(),simulationID,simulationID,dt,tmax,n_time_steps,n_nonlinear_steps,nonlinear_tolerance,max_linear_iters,max_r_steps,write_interval,directory.c_str(),xdmf.c_str());
	system(strdup(buffer));

	sprintf(buffer,"%s-file -dataflow sedimentation -transformation getMaximumIterations -id %d -name \"%s\" -path %s",pgCommandLine.c_str(),simulationID,xdmf.c_str(),directory.c_str());
	system(strdup(buffer));
	
	sprintf(buffer,"%s-performance -endtime -dataflow sedimentation -transformation getMaximumIterations -task %d -computation libMeshSedimentation::GetMaximumIterations",pgCommandLine.c_str(),simulationID);
	system(strdup(buffer));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	sprintf(buffer,"%s-ingest -task sedimentation getMaximumIterations %d",pgCommandLine.c_str(),simulationID);
	system(strdup(buffer));
}

void Provenance::inputInitDataExtraction(int simulationID, string transformation, string extractionFileName){
	if(processor_id != 0) return; 
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the sediments
	char buffer[4096];
	sprintf(buffer,"%s-task -dataflow sedimentation -transformation %s -id %d -status RUNNING -workspace %s -dependencies [{getMaximumIterations},{%d}]",pgCommandLine.c_str(),transformation.c_str(),simulationID,directory.c_str(),simulationID);
	system(strdup(buffer));

	sprintf(buffer,"%s-performance -starttime -dataflow sedimentation -transformation %s -task %d -computation libMeshSedimentation::%s-%d",pgCommandLine.c_str(),transformation.c_str(),simulationID,transformation.c_str(),simulationID);
	system(strdup(buffer));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	sprintf(buffer,"%s-ingest -task sedimentation %s %d",pgCommandLine.c_str(),transformation.c_str(),simulationID);
	system(strdup(buffer));
}

void Provenance::outputInitDataExtraction(int simulationID, string transformation, string extractionFileName, string outDataSet, int time_step, string xdmf, string rawDataFile){
	if(processor_id != 0) return; 
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the sediments
	char buffer[4096];
	sprintf(buffer,"%s-task -dataflow sedimentation -transformation %s -id %d -status FINISHED -workspace %s -dependencies [{getMaximumIterations},{%d}]",pgCommandLine.c_str(),transformation.c_str(),simulationID,directory.c_str(),simulationID);
	system(strdup(buffer));

	sprintf(buffer,"%s-file -dataflow sedimentation -transformation %s -id %d -name \"%s\" -path %s",pgCommandLine.c_str(),transformation.c_str(),simulationID,xdmf.c_str(),directory.c_str());
	system(strdup(buffer));

	sprintf(buffer,"%s-file -dataflow sedimentation -transformation %s -id %d -name \"%s\" -path %s",pgCommandLine.c_str(),transformation.c_str(),simulationID,rawDataFile.c_str(),directory.c_str());
	system(strdup(buffer));

	sprintf(buffer,"%s-element -dataflow sedimentation -transformation %s -id %d -set %s -element [{'%d;%d;%s/%s;%s/%s'}]",pgCommandLine.c_str(),transformation.c_str(),simulationID,outDataSet.c_str(),simulationID,time_step,directory.c_str(),xdmf.c_str(),directory.c_str(),rawDataFile.c_str());
	system(strdup(buffer));

	sprintf(buffer,"%s-performance -endtime -dataflow sedimentation -transformation %s -task %d -computation libMeshSedimentation::%s-%d",pgCommandLine.c_str(),transformation.c_str(),simulationID,transformation.c_str(),simulationID);
	system(strdup(buffer));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	sprintf(buffer,"%s-ingest -task sedimentation %s %d",pgCommandLine.c_str(),transformation.c_str(),simulationID);
	system(strdup(buffer));
}

void Provenance::inputSolverSimulationFluid(int simulationID, int subTaskID){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the fluid
	char buffer[4096];
	sprintf(buffer,"%s-task -dataflow sedimentation -transformation solverSimulationFluid -id %d -status RUNNING -workspace %s -invocation SolverSimulationFluid -subid %d -dependencies [{getMaximumIterations},{%d}]",pgCommandLine.c_str(),simulationID,directory.c_str(),subTaskID,simulationID);
	system(strdup(buffer));

	sprintf(buffer,"%s-performance -starttime -dataflow sedimentation -transformation solverSimulationFluid -task %d -subtask %d -computation libMeshSedimentation::SolverSimulationFluid-%d-%d",pgCommandLine.c_str(),simulationID,subTaskID,simulationID,subTaskID);
	system(strdup(buffer));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	sprintf(buffer,"%s-ingest -task sedimentation solverSimulationFluid %d %d",pgCommandLine.c_str(),simulationID,subTaskID);
	system(strdup(buffer));
}

void Provenance::outputSolverSimulationFluid(int simulationID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, 
	Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the fluid
	char buffer[4096];
	sprintf(buffer,"%s-task -dataflow sedimentation -transformation solverSimulationFluid -id %d -status FINISHED -workspace %s -subid %d -dependencies [{getMaximumIterations},{%d}]",pgCommandLine.c_str(),simulationID,directory.c_str(),subTaskID,simulationID);
	system(strdup(buffer));

	sprintf(buffer,"%s-performance -endtime -dataflow sedimentation -transformation solverSimulationFluid -task %d -subtask %d -computation libMeshSedimentation::SolverSimulationFluid-%d-%d",pgCommandLine.c_str(),simulationID,subTaskID,simulationID,subTaskID);
	system(strdup(buffer));

	sprintf(buffer,"%s-element -dataflow sedimentation -transformation solverSimulationFluid -id %d -subid %d -set osolversimulationfluid -element [{'%d;%d;%.2f;%d;%d;%d;%.2f;%.2f;%.2f;%s'}]",pgCommandLine.c_str(),simulationID,subTaskID,simulationID,time_step,time,linear_step,n_linear_step,n_linear_iterations,linear_residual,norm_delta,norm_delta_u,converged ? "true" : "false");
	system(strdup(buffer));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	sprintf(buffer,"%s-ingest -task sedimentation solverSimulationFluid %d %d",pgCommandLine.c_str(),simulationID,subTaskID);
	system(strdup(buffer));
}

void Provenance::inputSolverSimulationSediments(int simulationID, int subTaskID){
	if(processor_id != 0) return; 
	
	Performance perf;
	perf.start();

	// pg
	// solver simulation to the fluid
	char buffer[4096];
	sprintf(buffer,"%s-task -dataflow sedimentation -transformation solverSimulationSediments -id %d -status RUNNING -workspace %s -invocation SolverSimulationSediments -subid %d -dependencies [{solverSimulationFluid},{%d}]",pgCommandLine.c_str(),simulationID,directory.c_str(),subTaskID,simulationID);
	system(strdup(buffer));

	sprintf(buffer,"%s-performance -starttime -dataflow sedimentation -transformation solverSimulationSediments -task %d -subtask %d -computation libMeshSedimentation::SolverSimulationSediments-%d-%d",pgCommandLine.c_str(),simulationID,subTaskID,simulationID,subTaskID);
	system(strdup(buffer));

	perf.end();
  	double elapsedTime = perf.elapsedTime();

	// ingest
	sprintf(buffer,"%s-ingest -task sedimentation solverSimulationSediments %d %d",pgCommandLine.c_str(),simulationID,subTaskID);
	system(strdup(buffer));
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



