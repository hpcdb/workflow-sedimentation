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
#include "dfanalyzer/extractor.h"

#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"

#include "provenance.h"
#include "performance.h"

using namespace std;
using namespace libMesh;

Provenance::Provenance() {
    GetPot infile("provenance.in");
    directory = infile("directory", "/Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/sedimentation");
    outputDirectory = infile("outputDirectory", "/Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/sedimentation");
    string pgFilePath = infile("pgFilePath", "/Users/vitor/Documents/Repository/Thesis/Workflow-Sedimentation/dfa/PG-1.0.jar");
    string rdeFilePath = infile("rdeFilePath", "/Users/vitor/Documents/Repository/Thesis/Workflow-Sedimentation/dfa/RDE-1.0.jar");
    rawDataAccess = infile("access", "EXTRACTION");
    cartridge = infile("cartridge", "CSV");
    pgCommandLine = "java -jar " + pgFilePath + " ";
    rdeCommandLine = "java -jar " + rdeFilePath + " ";
    jsonDirectory = directory + "/prov/di/" + dataflow + "/";
    processor_id = libMesh::global_processor_id();
}

void Provenance::inputMeshGeneration(int simulationID, int dim, int ncellx, int ncelly, int ncellz,
        double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, int ref_interval) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Mesh Generation" << endl;
#endif    

    Performance perf;
    perf.start();

    string transformation = "meshgeneration";
    PerformanceMetric p;
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(simulationID);
    t.addPerformanceMetric(p);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("RUNNING");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d;%d;%d;%d;%d;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%d", simulationID, dim, ncellx, ncelly, ncellz, xmin, ymin, zmin, xmax, ymax, zmax, ref_interval);
    vector<string> e = {memalloc};
    t.addSet("i" + transformation, e);

    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputMeshGeneration(int simulationID, double r_fraction, double c_fraction,
        double max_h_level, unsigned int hlevels) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Mesh Generation" << endl;
#endif

    Performance perf;
    perf.start();

    string transformation = "meshgeneration";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d;%.2f;%.2f;%.2f;%d", simulationID, r_fraction, c_fraction, max_h_level, hlevels);
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();

    //    Create Equation Systems
    perf.start();

    transformation = "createequationsystems";
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t2(simulationID);
    t2.addPerformanceMetric(p);
    t2.setDataflow(dataflow);
    t2.setTransformation(transformation);
    t2.setWorkspace(directory);
    t2.setStatus("RUNNING");
    t2.addDtDependency("meshgeneration");

    sprintf(memalloc, "%d", simulationID);
    t2.addIdDependency(memalloc);

    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);

    perf.end();
    elapsedTime = perf.elapsedTime();

    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputCreateEquationSystems(int simulationID, Real Reynolds, Real Gr,
        Real Sc, Real Us, Real Diffusivity, Real xlock, Real fopc,
        Real theta, Real ex, Real ey, Real ez, Real c_factor) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Create Equation Systems" << endl;
#endif

    Performance perf;
    //        Create Equation Systems
    perf.start();

    string transformation = "createequationsystems";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f",
            simulationID, Reynolds, Gr, Sc, Us, Diffusivity, xlock, fopc, theta, ex, ey, ez, c_factor);
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();

    perf.start();

    transformation = "getmaximumiterations";
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t2(simulationID);
    t2.addPerformanceMetric(p);
    t2.setDataflow(dataflow);
    t2.setTransformation(transformation);
    t2.setWorkspace(directory);
    t2.setStatus("RUNNING");
    t2.addDtDependency("createequationsystems");

    sprintf(memalloc, "%d", simulationID);
    t2.addIdDependency(memalloc);

    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);

    perf.end();
    elapsedTime = perf.elapsedTime();

    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputGetMaximumIterations(int simulationID, Real dt, Real tmax,
        unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance,
        int max_linear_iters, int max_r_steps, unsigned int write_interval, string xdmf) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Get Maximum Iterations" << endl;
#endif

    Performance perf;
    perf.start();

    string transformation = "getmaximumiterations";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");

    char memalloc[1000];
    sprintf(memalloc, "%d;%.2f;%.2f;%d;%d;%.2f;%d;%d;%d;%s/%s",
            simulationID, dt, tmax, n_time_steps, n_nonlinear_steps,
            nonlinear_tolerance, max_linear_iters, max_r_steps,
            write_interval, directory.c_str(), xdmf.c_str());
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    File f(directory, xdmf);
    t.addFile(f);

    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::inputInitDataExtraction(int simulationID, string transformation) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Init Data Extraction" << endl;
#endif

    Performance perf;
    perf.start();

    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d", transformation.c_str(), simulationID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(simulationID);
    t.addPerformanceMetric(p);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("RUNNING");
    t.addDtDependency("getmaximumiterations");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputInitDataExtraction(int simulationID, string transformation, string dataSet,
        int time_step, string xdmf, string rawDataFile, int dimension, string extractorName) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Init Data Extraction" << endl;
#endif
    Performance perf;
    perf.start();

    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("getmaximumiterations");

    char memalloc[4096];
    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();

    Performance rdePerf;
    rdePerf.start();

    string extension = "data";
    if (rawDataAccess == "INDEXING") {
        extension = "index";
        Extractor ext(rdeCommandLine, rawDataAccess, cartridge, extractorName);
        ext.addAttribute("u", "numeric", false);
        ext.addAttribute("v", "numeric", false);
        if (dimension == 3) {
            ext.addAttribute("w", "numeric", false);
        }
        ext.addAttribute("p", "numeric", false);
        ext.addAttribute("s", "numeric", false);
        ext.addAttribute("d", "numeric", false);
        if (dimension == 3) {
            ext.addAttribute("vtkvalidpointmask", "numeric", false);
            ext.addAttribute("arc_length", "numeric", false);
        }
        ext.addAttribute("points0", "numeric", false);
        ext.addAttribute("points1", "numeric", false);
        ext.addAttribute("points2", "numeric", false);
        ext.extract(directory, rawDataFile);
    }

    rdePerf.end();
    storeRDEComponentCost(rdePerf.elapsedTime());

    perf.start();

    char extractedFileName[jsonArraySize];
    if (rawDataAccess == "INDEXING") {
        sprintf(extractedFileName, "%s.%s", extractorName.c_str(), extension.c_str());
    } else {
        sprintf(extractedFileName, "%s", rawDataFile.c_str());
    }

    sprintf(memalloc, "%d;%d;%s/%s;%s/%s",
            simulationID, time_step, directory.c_str(), xdmf.c_str(),
            directory.c_str(), extractedFileName);
    vector<string> e = {memalloc};
    t.addSet(dataSet, e);

    PerformanceMetric p;
    sprintf(memalloc, "libMeshSedimentation::%s-%d", transformation.c_str(), simulationID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    File f1(directory, xdmf);
    t.addFile(f1);

    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);

    perf.end();
    elapsedTime = perf.elapsedTime();

    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::inputSolverSimulationFluid(int taskID, int simulationID, int subTaskID) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Solver Simulation Fluid" << endl;
#endif

    Performance perf;
    perf.start();

    string transformation = "solversimulationfluid";
    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), simulationID, subTaskID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(taskID);
    t.addPerformanceMetric(p);
    t.setSubID(subTaskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("RUNNING");
    t.addDtDependency("getmaximumiterations");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputSolverSimulationFluid(int taskID, int simulationID, int subTaskID,
        int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations,
        Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Solver Simulation Fluid" << endl;
#endif

    Performance perf;
    perf.start();

    string transformation = "solversimulationfluid";
    Task t(taskID);
    t.setSubID(subTaskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("getmaximumiterations");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%d;%.2f;%d;%d;%d;%.2f;%.2f;%.2f;%s",
            simulationID, time_step, time, linear_step, n_linear_step,
            n_linear_iterations, linear_residual, norm_delta,
            norm_delta_u, converged ? "true" : "false");
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), simulationID, subTaskID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::inputSolverSimulationSediments(int taskID, int simulationID, int subTaskID) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Solver Simulation Sediments" << endl;
#endif
    Performance perf;
    perf.start();

    string transformation = "solversimulationsediments";
    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), simulationID, subTaskID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(taskID);
    t.addPerformanceMetric(p);
    t.setSubID(subTaskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("RUNNING");
    t.addDtDependency("solversimulationfluid");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputSolverSimulationSediments(int taskID, int simulationID, int subTaskID, int time_step,
        Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations,
        Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Solver Simulation Sediments" << endl;
#endif

    Performance perf;
    perf.start();

    string transformation = "solversimulationsediments";
    Task t(taskID);
    t.setSubID(subTaskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("solversimulationfluid");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%d;%.2f;%d;%d;%d;%.2f;%.2f;%.2f;%s",
            simulationID, time_step, time, linear_step, n_linear_step,
            n_linear_iterations, linear_residual, norm_delta, norm_delta_u,
            converged ? "true" : "false");
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), simulationID, subTaskID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputMeshRefinement(int taskID, int simulationID, int subTaskID,
        bool first_step_refinement, int time_step, int before_n_active_elem, int after_n_active_elem) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Mesh Refinement" << endl;
#endif
    Performance perf;
    perf.start();

    string transformation = "meshrefinement";
    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), simulationID, subTaskID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(taskID);
    t.setSubID(subTaskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("solversimulationsediments");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%s;%d;%d;%d",
            simulationID, first_step_refinement ? "true" : "false", time_step,
            before_n_active_elem, after_n_active_elem);
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::inputMeshWriter(int taskID, int simulationID, int subTaskID) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Mesh Writer" << endl;
#endif

    Performance perf;
    perf.start();

    string transformation = "meshwriter";
    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), simulationID, subTaskID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(taskID);
    t.addPerformanceMetric(p);
    t.setSubID(subTaskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("RUNNING");
    t.addDtDependency("solversimulationsediments");

    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputMeshWriter(int taskID, int simulationID, int subTaskID, int time_step, string xdmf) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Mesh Writer" << endl;
#endif

    Performance perf;
    perf.start();

    string transformation = "meshwriter";
    Task t(taskID);
    t.setSubID(subTaskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("solversimulationsediments");

    char memalloc[4096];
    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

    File f1(directory, xdmf);
    t.addFile(f1);

    sprintf(memalloc, "%d;%d;%s/%s",
            simulationID, time_step, directory.c_str(), xdmf.c_str());
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), simulationID, subTaskID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::inputDataExtraction(int taskID, int simulationID, int subTaskID,
        string transformation) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Data Extraction" << endl;
#endif

    Performance perf;
    perf.start();

    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), simulationID, subTaskID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(taskID);
    t.addPerformanceMetric(p);
    t.setSubID(subTaskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("RUNNING");
    t.addDtDependency("meshwriter");

    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputDataExtraction(int taskID, int simulationID, int subTaskID,
        string transformation, string dataSet, int time_step,
        string xdmf, string rawDataFile, int dimension, string extractorName) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Data Extraction" << endl;
#endif

    Performance perf;
    perf.start();

    Task t(taskID);
    t.setSubID(subTaskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("meshwriter");

    char memalloc[4096];
    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();

    Performance rdePerf;
    rdePerf.start();

    string extension = "data";
    if (rawDataAccess == "INDEXING") {
        extension = "index";
        Extractor ext(rdeCommandLine, rawDataAccess, cartridge, extractorName);
        ext.addAttribute("u", "numeric", false);
        ext.addAttribute("v", "numeric", false);
        if (dimension == 3) {
            ext.addAttribute("w", "numeric", false);
        }
        ext.addAttribute("p", "numeric", false);
        ext.addAttribute("s", "numeric", false);
        ext.addAttribute("d", "numeric", false);
        if (dimension == 3) {
            ext.addAttribute("vtkvalidpointmask", "numeric", false);
            ext.addAttribute("arc_length", "numeric", false);
        }
        ext.addAttribute("points0", "numeric", false);
        ext.addAttribute("points1", "numeric", false);
        ext.addAttribute("points2", "numeric", false);
        ext.extract(directory, rawDataFile);
    }

    rdePerf.end();
    storeRDEComponentCost(rdePerf.elapsedTime());

    perf.start();

    char extractedFileName[jsonArraySize];
    if (rawDataAccess == "INDEXING") {
        sprintf(extractedFileName, "%s.%s", extractorName.c_str(), extension.c_str());
    } else {
        sprintf(extractedFileName, "%s", rawDataFile.c_str());
    }

    sprintf(memalloc, "%d;%d;%s/%s;%s/%s",
            simulationID, time_step, directory.c_str(), xdmf.c_str(),
            directory.c_str(), extractedFileName);
    vector<string> e = {memalloc};
    t.addSet(dataSet, e);

    File f1(directory, xdmf);
    t.addFile(f1);

    PerformanceMetric p;
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), simulationID, subTaskID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);

    perf.end();
    elapsedTime = perf.elapsedTime();

    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::meshAggregator(int simulationID, string xdmf, int n_processors, vector<string> meshDependencies) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Mesh Writer" << endl;
#endif

    Performance perf;
    perf.start();
    
    PerformanceMetric ptemp;
    ptemp.IdentifyStartTime();

    string transformation = "meshaggregator";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("meshwriter");
    for (string dep : meshDependencies) {
        t.addIdDependency(dep);
    }

    File f1(directory, xdmf);
    t.addFile(f1);

    char memalloc[4096];
    sprintf(memalloc, "%d;%s/%s;%d",
            simulationID, directory.c_str(), xdmf.c_str(), n_processors);
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    sprintf(memalloc, "libMeshSedimentation::%s-%d",
            transformation.c_str(), simulationID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.SetStartTime(ptemp.GetStartTime());
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
    
//    if (processor_id != 0) return;
//    Performance perf;
//    perf.start();
//
//    string transformation = "meshaggregator";
//    PerformanceMetric p;
//    char memalloc[jsonArraySize];
//    sprintf(memalloc, "libMeshSedimentation::%s-%d",
//            transformation.c_str(), simulationID);
//    p.SetDescription(memalloc);
//    p.SetMethod("COMPUTATION");
//    p.IdentifyStartTime();
//
//    Task t(simulationID);
//    t.setDataflow(dataflow);
//    t.setTransformation(transformation);
//    t.setWorkspace(directory);
//    t.setStatus("FINISHED");
    
//
//    File f1(directory, xdmf);
//    t.addFile(f1);
//
//    sprintf(memalloc, "%d;%s/%s;%d",
//            simulationID, directory.c_str(), xdmf.c_str(), n_processors);
//    vector<string> e = {memalloc};
//    t.addSet("o" + transformation, e);
//
//    p.IdentifyEndTime();
//    t.addPerformanceMetric(p);
//
//    sprintf(memalloc, "%s%s-%d-F.json",
//            jsonDirectory.c_str(), transformation.c_str(), simulationID);
//    t.writeJSON(memalloc);
//
//    perf.end();
//    double elapsedTime = perf.elapsedTime();
//
//    ofstream file;
//    file.open("prov/log/" + transformation + ".prov", ios_base::app);
//    file << "PROV:" + transformation + ":Output" << endl;
//    sprintf(memalloc, "%.5f", elapsedTime);
//    file << space << memalloc << endl;
//    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
//    file.close();
}

void Provenance::storeDataExtractionCost(double elapsedTime) {
    if (processor_id != 0) return;
    ofstream file;
    file.open("prov/rde/data-extraction.prov", ios_base::app);
    file << "RDE:DataExtraction:Process" << endl;
    char buffer[jsonArraySize];
    sprintf(buffer, "%.5f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::storeRDEComponentCost(double elapsedTime) {
    if (processor_id != 0) return;
    ofstream file;
    file.open("prov/indexing/rde-component.prov", ios_base::app);
    file << "RDEComponent:DataExtraction:Process" << endl;
    char buffer[jsonArraySize];
    sprintf(buffer, "%.5f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::storeSolverCost(double elapsedTime) {
    if (processor_id != 0) return;
    ofstream file;
    file.open("prov/solver/time.prov", ios_base::app);
    file << "Solver:Time:Process" << endl;
    char buffer[jsonArraySize];
    sprintf(buffer, "%.5f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();
}

void Provenance::finishDataIngestor() {
    if (processor_id != 0) return;

    string str = "cp ../dfa/finish.token prov/di/sedimentation";
    system(strdup(str.c_str()));

    cout << "[Provenance] Finish Data Ingestor" << endl;
}



