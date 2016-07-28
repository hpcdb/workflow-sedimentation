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

#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"

#include "dfanalyzer/task.h"

#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"

#include "provenance.h"
#include "performance.h"

using namespace std;
using namespace libMesh;

string space = "      ";
string directory = "";
string pgCommandLine = "";
string dataflow = "sedimentation";
string jsonDirectory = "";

Provenance::Provenance() {
    GetPot infile("provenance.in");
    directory = infile("directory", "/Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/sedimentation");
    string pgFilePath = infile("pgFilePath", "/Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/dfa/PG-1.0.jar");
    pgCommandLine = "java -jar " + pgFilePath + " ";
    jsonDirectory = directory + "/prov/di/" + dataflow + "/";
    processor_id = libMesh::global_processor_id();
}

void Provenance::inputMeshGeneration(int simulationID, int dim, int ncellx, int ncelly, int ncellz,
        double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, int ref_interval) {
    if (processor_id != 0) return;
    
    Performance perf;
    perf.start();

    string transformation = "meshgeneration";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    
    PerformanceMetric p;
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    t.addPerformance(p);
//    falta starttime e endtime
    
    char* element = (char*) malloc(4096);
    sprintf(element, "%d;%d;%d;%d;%d;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%d", simulationID, dim, ncellx, ncelly, ncellz, xmin, ymin, zmin, xmax, ymax, zmax, ref_interval);
    cout << element << endl;
    vector<string> e = {element};
    t.addSet("i" + transformation,e);
    free(element);
    
    char* buffer = (char*) malloc(4096);
    sprintf(buffer, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    cout << buffer << endl;
    t.writeJSON(buffer);
    free(buffer);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    char* bPointer = (char*) malloc(512);
    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    file << space << *bPointer << endl;
    file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
    file.close();
    free(bPointer);
}

void Provenance::outputMeshGeneration(int simulationID, double r_fraction, double c_fraction, double max_h_level, unsigned int hlevels) {
    if (processor_id != 0) return;

    Performance perf;
    perf.start();

    // run PG
    // mesh generation
    char buffer[4096];
    sprintf(buffer, "%s-element -dataflow sedimentation -transformation meshGeneration -id %d -set omeshgeneration -element [{'%d;%.2f;%.2f;%.2f;%d'}]", pgCommandLine.c_str(), simulationID, simulationID, r_fraction, c_fraction, max_h_level, hlevels);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -endtime -dataflow sedimentation -transformation meshGeneration -task %d -computation libMeshSedimentation::MeshGeneration", pgCommandLine.c_str(), simulationID);
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation meshGeneration %d", pgCommandLine.c_str(), simulationID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/MeshGeneration.prov", ios_base::app);
    file << "PROV:MeshGeneration:Output" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();

    perf.start();

    // create equation systems
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation createEquationSystems -id %d -workspace %s -status FINISHED -dependencies [{meshGeneration},{%d}]", pgCommandLine.c_str(), simulationID, directory.c_str(), simulationID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -starttime -dataflow sedimentation -transformation createEquationSystems -task %d -computation libMeshSedimentation::CreateEquationSystems", pgCommandLine.c_str(), simulationID);
    //cout << buffer << endl;

    perf.end();
    elapsedTime += perf.elapsedTime();

    file.open("prov/log/CreateEquationSystems.prov", ios_base::app);
    file << "PROV:CreateEquationSystems:Input" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::outputCreateEquationSystems(int simulationID, Real Reynolds, Real Gr, Real Sc, Real Us, Real Diffusivity, Real xlock, Real fopc,
        Real theta, Real ex, Real ey, Real ez, Real c_factor) {
    if (processor_id != 0) return;

    Performance perf;
    perf.start();

    // run PG
    // create equation systems
    char buffer[4096];
    sprintf(buffer, "%s-element -dataflow sedimentation -transformation createEquationSystems -id %d -set ocreateequationsystems -element [{'%d;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f'}]", pgCommandLine.c_str(), simulationID, simulationID, Reynolds, Gr, Sc, Us, Diffusivity, xlock, fopc, theta, ex, ey, ez, c_factor);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -endtime -dataflow sedimentation -transformation createEquationSystems -task %d -computation libMeshSedimentation::createEquationSystems -dependencies [{meshGeneration},{%d}]", pgCommandLine.c_str(), simulationID, simulationID);
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation createEquationSystems %d", pgCommandLine.c_str(), simulationID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/CreateEquationSystems.prov", ios_base::app);
    file << "PROV:CreateEquationSystems:Output" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();

    perf.start();

    // get maximum iterations
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation getMaximumIterations -id %d -workspace %s -status FINISHED -dependencies [{createEquationSystems},{%d}]", pgCommandLine.c_str(), simulationID, directory.c_str(), simulationID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -starttime -dataflow sedimentation -transformation getMaximumIterations -task %d -computation libMeshSedimentation::GetMaximumIterations", pgCommandLine.c_str(), simulationID);
    //cout << buffer << endl;

    perf.end();
    elapsedTime += perf.elapsedTime();

    file.open("prov/log/GetMaximumIterations.prov", ios_base::app);
    file << "PROV:GetMaximumIterations:Input" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::outputGetMaximumIterations(int simulationID, Real dt, Real tmax, unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance,
        int max_linear_iters, int max_r_steps, unsigned int write_interval, string xdmf) {
    if (processor_id != 0) return;

    Performance perf;
    perf.start();

    // pg
    // get maximum iterations
    char buffer[4096];
    sprintf(buffer, "%s-element -dataflow sedimentation -transformation getMaximumIterations -id %d -set ogetmaximumiterations -element [{'%d;%.2f;%.2f;%d;%d;%.2f;%d;%d;%d;%s/%s'}]", pgCommandLine.c_str(), simulationID, simulationID, dt, tmax, n_time_steps, n_nonlinear_steps, nonlinear_tolerance, max_linear_iters, max_r_steps, write_interval, directory.c_str(), xdmf.c_str());
    //cout << buffer << endl;

    sprintf(buffer, "%s-file -dataflow sedimentation -transformation getMaximumIterations -id %d -name \"%s\" -path %s", pgCommandLine.c_str(), simulationID, xdmf.c_str(), directory.c_str());
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -endtime -dataflow sedimentation -transformation getMaximumIterations -task %d -computation libMeshSedimentation::GetMaximumIterations", pgCommandLine.c_str(), simulationID);
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation getMaximumIterations %d", pgCommandLine.c_str(), simulationID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/GetMaximumIterations.prov", ios_base::app);
    file << "PROV:GetMaximumIterations:Output" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::inputInitDataExtraction(int simulationID, string transformation, string extractionFileName) {
    if (processor_id != 0) return;
    Performance perf;
    perf.start();

    // pg
    // solver simulation to the sediments
    char buffer[4096];
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation %s -id %d -status RUNNING -workspace %s -dependencies [{getMaximumIterations},{%d}]", pgCommandLine.c_str(), transformation.c_str(), simulationID, directory.c_str(), simulationID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -starttime -dataflow sedimentation -transformation %s -task %d -computation libMeshSedimentation::%s-%d", pgCommandLine.c_str(), transformation.c_str(), simulationID, transformation.c_str(), simulationID);
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation %s %d", pgCommandLine.c_str(), transformation.c_str(), simulationID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/InitDataExtraction.prov", ios_base::app);
    file << "PROV:InitDataExtraction:Input" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::outputInitDataExtraction(int simulationID, string transformation, string extractionFileName, string outDataSet, int time_step, string xdmf, string rawDataFile) {
    if (processor_id != 0) return;
    Performance perf;
    perf.start();

    // pg
    // solver simulation to the sediments
    char buffer[4096];
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation %s -id %d -status FINISHED -workspace %s -dependencies [{getMaximumIterations},{%d}]", pgCommandLine.c_str(), transformation.c_str(), simulationID, directory.c_str(), simulationID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-file -dataflow sedimentation -transformation %s -id %d -name \"%s\" -path %s", pgCommandLine.c_str(), transformation.c_str(), simulationID, xdmf.c_str(), directory.c_str());
    //cout << buffer << endl;

    sprintf(buffer, "%s-file -dataflow sedimentation -transformation %s -id %d -name \"%s\" -path %s", pgCommandLine.c_str(), transformation.c_str(), simulationID, rawDataFile.c_str(), directory.c_str());
    //cout << buffer << endl;

    sprintf(buffer, "%s-element -dataflow sedimentation -transformation %s -id %d -set %s -element [{'%d;%d;%s/%s;%s/%s'}]", pgCommandLine.c_str(), transformation.c_str(), simulationID, outDataSet.c_str(), simulationID, time_step, directory.c_str(), xdmf.c_str(), directory.c_str(), rawDataFile.c_str());
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -endtime -dataflow sedimentation -transformation %s -task %d -computation libMeshSedimentation::%s-%d", pgCommandLine.c_str(), transformation.c_str(), simulationID, transformation.c_str(), simulationID);
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation %s %d", pgCommandLine.c_str(), transformation.c_str(), simulationID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/InitDataExtraction.prov", ios_base::app);
    file << "PROV:InitDataExtraction:Output" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::inputSolverSimulationFluid(int taskID, int simulationID, int subTaskID) {
    if (processor_id != 0) return;

    Performance perf;
    perf.start();

    // pg
    // solver simulation to the fluid
    char buffer[4096];
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation solverSimulationFluid -id %d -status RUNNING -workspace %s -invocation SolverSimulationFluid -subid %d -dependencies [{getMaximumIterations},{%d}]", pgCommandLine.c_str(), taskID, directory.c_str(), subTaskID, simulationID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -starttime -dataflow sedimentation -transformation solverSimulationFluid -task %d -subtask %d -computation libMeshSedimentation::SolverSimulationFluid-%d-%d", pgCommandLine.c_str(), taskID, subTaskID, simulationID, subTaskID);
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation solverSimulationFluid %d %d", pgCommandLine.c_str(), taskID, subTaskID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/SolverSimulationFluid.prov", ios_base::app);
    file << "PROV:SolverSimulationFluid:Input" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::outputSolverSimulationFluid(int taskID, int simulationID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations,
        Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged) {
    if (processor_id != 0) return;

    Performance perf;
    perf.start();

    // pg
    // solver simulation to the fluid
    char buffer[4096];
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation solverSimulationFluid -id %d -status FINISHED -workspace %s -subid %d -dependencies [{getMaximumIterations},{%d}]",
            pgCommandLine.c_str(), taskID, directory.c_str(), subTaskID, simulationID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -endtime -dataflow sedimentation -transformation solverSimulationFluid -task %d -subtask %d -computation libMeshSedimentation::SolverSimulationFluid-%d-%d",
            pgCommandLine.c_str(), taskID, subTaskID, simulationID, subTaskID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-element -dataflow sedimentation -transformation solverSimulationFluid -id %d -subid %d -set osolversimulationfluid -element [{'%d;%d;%.2f;%d;%d;%d;%.2f;%.2f;%.2f;%s'}]",
            pgCommandLine.c_str(), taskID, subTaskID, simulationID, time_step, time, linear_step, n_linear_step, n_linear_iterations, linear_residual, norm_delta, norm_delta_u, converged ? "true" : "false");
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation solverSimulationFluid %d %d", pgCommandLine.c_str(), taskID, subTaskID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/SolverSimulationFluid.prov", ios_base::app);
    file << "PROV:SolverSimulationFluid:Output" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::inputSolverSimulationSediments(int taskID, int simulationID, int subTaskID) {
    if (processor_id != 0) return;

    Performance perf;
    perf.start();

    // pg
    // solver simulation to the fluid
    char buffer[4096];
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation solverSimulationSediments -id %d -status RUNNING -workspace %s -invocation SolverSimulationSediments -subid %d -dependencies [{solverSimulationFluid},{%d}]",
            pgCommandLine.c_str(), taskID, directory.c_str(), subTaskID, taskID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -starttime -dataflow sedimentation -transformation solverSimulationSediments -task %d -subtask %d -computation libMeshSedimentation::SolverSimulationSediments-%d-%d",
            pgCommandLine.c_str(), taskID, subTaskID, simulationID, subTaskID);
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation solverSimulationSediments %d %d", pgCommandLine.c_str(), taskID, subTaskID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/SolverSimulationSediments.prov", ios_base::app);
    file << "PROV:SolverSimulationSediments:Input" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::outputSolverSimulationSediments(int taskID, int simulationID, int subTaskID, int time_step,
        Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations,
        Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged) {
    if (processor_id != 0) return;

    Performance perf;
    perf.start();

    // pg
    // solver simulation to the sediments
    char buffer[4096];
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation solverSimulationSediments -id %d -status FINISHED -workspace %s -subid %d -dependencies [{solverSimulationFluid},{%d}]",
            pgCommandLine.c_str(), taskID, directory.c_str(), subTaskID, taskID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -endtime -dataflow sedimentation -transformation solverSimulationSediments -task %d -subtask %d -computation libMeshSedimentation::SolverSimulationSediments-%d-%d",
            pgCommandLine.c_str(), taskID, subTaskID, simulationID, subTaskID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-element -dataflow sedimentation -transformation solverSimulationSediments -id %d -subid %d -set osolversimulationsediments -element [{'%d;%d;%.2f;%d;%d;%d;%.2f;%.2f;%.2f;%s'}]",
            pgCommandLine.c_str(), taskID, subTaskID, simulationID, time_step, time, linear_step,
            n_linear_step, n_linear_iterations, linear_residual, norm_delta, norm_delta_u, converged ? "true" : "false");
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation solverSimulationSediments %d %d",
            pgCommandLine.c_str(), taskID, subTaskID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/SolverSimulationSediments.prov", ios_base::app);
    file << "PROV:SolverSimulationSediments:Output" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::outputMeshRefinement(int taskID, int simulationID, int subTaskID, bool first_step_refinement, int time_step, int before_n_active_elem, int after_n_active_elem) {
    if (processor_id != 0) return;

    Performance perf;
    perf.start();

    // run PG
    // mesh refinement
    char buffer[4096];
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation meshRefinement -id %d -workspace %s -status FINISHED -subid %d -dependencies [{solverSimulationSediments},{%d}]",
            pgCommandLine.c_str(), taskID, directory.c_str(), subTaskID, taskID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -starttime -dataflow sedimentation -transformation meshRefinement -task %d -subtask %d -computation libMeshSedimentation::MeshRefinement-%d-%d",
            pgCommandLine.c_str(), taskID, subTaskID, simulationID, subTaskID);
    //cout << buffer << endl;

    // input element
    sprintf(buffer, "%s-element -dataflow sedimentation -transformation meshRefinement -id %d -subid %d -set omeshrefinement -element [{'%d;%s;%d;%d;%d'}]",
            pgCommandLine.c_str(), taskID, subTaskID, simulationID, first_step_refinement ? "true" : "false", time_step, before_n_active_elem, after_n_active_elem);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -endtime -dataflow sedimentation -transformation meshRefinement -task %d -subtask %d -computation libMeshSedimentation::MeshRefinement-%d-%d",
            pgCommandLine.c_str(), taskID, subTaskID, simulationID, subTaskID);
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation meshRefinement %d %d",
            pgCommandLine.c_str(), taskID, subTaskID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/MeshRefinement.prov", ios_base::app);
    file << "PROV:MeshRefinement:Input" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::inputMeshWriter(int taskID, int simulationID, int subTaskID) {
    if (processor_id != 0) return;
    Performance perf;
    perf.start();

    // pg
    // solver simulation to the fluid
    char buffer[4096];
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation meshWriter -id %d -status RUNNING -workspace %s -invocation MeshWriter -subid %d -dependencies [{solverSimulationSediments},{%d}]",
            pgCommandLine.c_str(), taskID, directory.c_str(), subTaskID, taskID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -starttime -dataflow sedimentation -transformation meshWriter -task %d -subtask %d -computation libMeshSedimentation::MeshWriter-%d-%d",
            pgCommandLine.c_str(), taskID, subTaskID, simulationID, subTaskID);
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation meshWriter %d %d",
            pgCommandLine.c_str(), taskID, subTaskID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/MeshWriter.prov", ios_base::app);
    file << "PROV:MeshWriter:Input" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();

    return;
}

void Provenance::outputMeshWriter(int taskID, int simulationID, int subTaskID, int time_step, string xdmf) {
    if (processor_id != 0) return;
    Performance perf;
    perf.start();

    // pg
    // solver simulation to the sediments
    char buffer[4096];
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation meshWriter -id %d -status FINISHED -workspace %s -subid %d -dependencies [{solverSimulationSediments},{%d}]",
            pgCommandLine.c_str(), taskID, directory.c_str(), subTaskID, taskID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -endtime -dataflow sedimentation -transformation meshWriter -task %d -subtask %d -computation libMeshSedimentation::MeshWriter-%d-%d",
            pgCommandLine.c_str(), taskID, subTaskID, simulationID, subTaskID);
    //cout << buffer << endl;

    sprintf(buffer, "%s -file -dataflow sedimentation -transformation meshWriter -id %d -subid %d -name \"%s\" -path %s",
            pgCommandLine.c_str(), taskID, subTaskID, xdmf.c_str(), directory.c_str());
    //cout << buffer << endl;

    sprintf(buffer, "%s-element -dataflow sedimentation -transformation meshWriter -id %d -subid %d -set omeshwriter -element [{'%d;%d;%s/%s'}]",
            pgCommandLine.c_str(), taskID, subTaskID, simulationID, time_step, directory.c_str(), xdmf.c_str());
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation meshWriter %d %d",
            pgCommandLine.c_str(), taskID, subTaskID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/MeshWriter.prov", ios_base::app);
    file << "PROV:MeshWriter:Output" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::inputDataExtraction(int taskID, int simulationID, int subTaskID, string transformation, string extractionFileName) {
    if (processor_id != 0) return;
    Performance perf;
    perf.start();

    // pg
    // solver simulation to the sediments
    char buffer[4096];
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation %s -id %d -status RUNNING -workspace %s -subid %d -dependencies [{meshWriter},{%d}]",
            pgCommandLine.c_str(), transformation.c_str(), taskID, directory.c_str(), subTaskID, taskID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -starttime -dataflow sedimentation -transformation %s -task %d -subtask %d -computation libMeshSedimentation::%s-%d",
            pgCommandLine.c_str(), transformation.c_str(), taskID, subTaskID, transformation.c_str(), simulationID);
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation %s %d %d",
            pgCommandLine.c_str(), transformation.c_str(), taskID, subTaskID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/DataExtraction.prov", ios_base::app);
    file << "PROV:DataExtraction:Input" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::outputDataExtraction(int taskID, int simulationID, int subTaskID, string transformation, string extractionFileName, string outDataSet, int time_step, string xdmf, string rawDataFile) {
    if (processor_id != 0) return;
    Performance perf;
    perf.start();

    // pg
    // solver simulation to the sediments
    char buffer[4096];
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation %s -id %d -status FINISHED -workspace %s -subid %d -dependencies [{meshWriter},{%d}]",
            pgCommandLine.c_str(), transformation.c_str(), taskID, directory.c_str(), subTaskID, taskID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -endtime -dataflow sedimentation -transformation %s -task %d -subtask %d -computation libMeshSedimentation::%s-%d",
            pgCommandLine.c_str(), transformation.c_str(), taskID, subTaskID, transformation.c_str(), simulationID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-file -dataflow sedimentation -transformation %s -id %d -subid %d -name \"%s\" -path %s",
            pgCommandLine.c_str(), transformation.c_str(), taskID, subTaskID, xdmf.c_str(), directory.c_str());
    //cout << buffer << endl;

    sprintf(buffer, "%s-file -dataflow sedimentation -transformation %s -id %d -subid %d -name \"%s\" -path %s",
            pgCommandLine.c_str(), transformation.c_str(), taskID, subTaskID, rawDataFile.c_str(), directory.c_str());
    //cout << buffer << endl;

    sprintf(buffer, "%s-element -dataflow sedimentation -transformation %s -id %d -subid %d -set %s -element [{'%d;%d;%s/%s;%s/%s'}]",
            pgCommandLine.c_str(), transformation.c_str(), taskID, subTaskID, outDataSet.c_str(), simulationID, time_step, directory.c_str(), xdmf.c_str(), directory.c_str(), rawDataFile.c_str());
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation %s %d %d",
            pgCommandLine.c_str(), transformation.c_str(), taskID, subTaskID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/DataExtraction.prov", ios_base::app);
    file << "PROV:DataExtraction:Output" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::meshAggregator(int simulationID, string xdmf, int n_processors, string meshDependencies) {
    if (processor_id != 0) return;
    Performance perf;
    perf.start();

    // pg
    // solver simulation to the sediments
    char buffer[4096];
    sprintf(buffer, "%s-task -dataflow sedimentation -transformation meshAggregator -id %d -status FINISHED -workspace %s -dependencies \"[{meshWriter},{%s}]\"",
            pgCommandLine.c_str(), simulationID, directory.c_str(), meshDependencies.c_str());
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -starttime -dataflow sedimentation -transformation meshAggregator -task %d -computation libMeshSedimentation::MeshAggregator-%d",
            pgCommandLine.c_str(), simulationID, simulationID);
    //cout << buffer << endl;

    sprintf(buffer, "%s-file -dataflow sedimentation -transformation meshAggregator -id %d -name \"%s\" -path %s",
            pgCommandLine.c_str(), simulationID, xdmf.c_str(), directory.c_str());
    //cout << buffer << endl;

    sprintf(buffer, "%s-element -dataflow sedimentation -transformation meshAggregator -id %d -set omeshaggregator -element [{'%d;%s/%s;%d'}]",
            pgCommandLine.c_str(), simulationID, simulationID, directory.c_str(), xdmf.c_str(), n_processors);
    //cout << buffer << endl;

    sprintf(buffer, "%s-performance -endtime -dataflow sedimentation -transformation meshAggregator -task %d -computation libMeshSedimentation::MeshAggregator-%d",
            pgCommandLine.c_str(), simulationID, simulationID);
    //cout << buffer << endl;

    perf.end();
    double elapsedTime = perf.elapsedTime();

    // ingest
    sprintf(buffer, "%s-ingest -task sedimentation meshAggregator %d",
            pgCommandLine.c_str(), simulationID);
    //cout << buffer << endl;

    ofstream file;
    file.open("prov/log/MeshAggregator.prov", ios_base::app);
    file << "PROV:MeshAggregator:Output" << endl;
    file << space << buffer << endl;
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::storeDataExtractionCost(int simulationID, int subTaskID, int time_step, string xdmf, string rawDataFile, double elapsedTime) {
    if (processor_id != 0) return;
    ofstream file;
    file.open("prov/rde/data-extraction.prov", ios_base::app);
    file << "RDE:DataExtraction:Process" << endl;
    char buffer[1024];
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::storeSolverCost(double elapsedTime) {
    if (processor_id != 0) return;
    ofstream file;
    file.open("prov/solver/time.prov", ios_base::app);
    file << "Solver:Time:Process" << endl;
    char buffer[1024];
    sprintf(buffer, "%.2f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::finishDataIngestor() {
    if (processor_id != 0) return;

    string str = "cp ../dfa/finish.token prov/di/sedimentation";
    system(strdup(str.c_str()));

    cout << "[Provenance] Finish Data Ingestor" << endl;
}



