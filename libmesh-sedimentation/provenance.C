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
    t.setStatus("RUNNING");

    char* element = (char*) malloc(4096);
    sprintf(element, "%d;%d;%d;%d;%d;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%d", simulationID, dim, ncellx, ncelly, ncellz, xmin, ymin, zmin, xmax, ymax, zmax, ref_interval);
    vector<string> e = {element};
    t.addSet("i" + transformation, e);
    free(element);

    char* buffer = (char*) malloc(4096);
    sprintf(buffer, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(buffer);
    free(buffer);

    perf.end();
    double elapsedTime = perf.elapsedTime();

    char* bPointer = (char*) malloc(512);
    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(bPointer, "%.2f", elapsedTime);
    file << space << *bPointer << endl;
    file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
    file.close();
    free(bPointer);
}

void Provenance::outputMeshGeneration(int simulationID, double r_fraction, double c_fraction, double max_h_level, unsigned int hlevels) {
    if (processor_id != 0) return;

    Performance perf;
    {
        perf.start();

        string transformation = "meshgeneration";
        Task t(simulationID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("FINISHED");

        char* element = (char*) malloc(4096);
        sprintf(element, "%d;%.2f;%.2f;%.2f;%d", simulationID, r_fraction, c_fraction, max_h_level, hlevels);
        vector<string> e = {element};
        t.addSet("o" + transformation, e);
        free(element);

        PerformanceMetric p;
        p.SetDescription("libMeshSedimentation::" + transformation);
        p.SetMethod("COMPUTATION");
        t.addPerformance(p);
        //    missing starttime and endtime

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Output" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }

    {
        //    Create Equation Systems
        perf.start();

        string transformation = "createequationsystems";
        Task t(simulationID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("RUNNING");
        t.addDtDependency("meshgeneration");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        PerformanceMetric p;
        p.SetDescription("libMeshSedimentation::" + transformation);
        p.SetMethod("COMPUTATION");
        t.addPerformance(p);
        //    missing starttime and endtime

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Input" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::outputCreateEquationSystems(int simulationID, Real Reynolds, Real Gr, Real Sc, Real Us, Real Diffusivity, Real xlock, Real fopc,
        Real theta, Real ex, Real ey, Real ez, Real c_factor) {
    if (processor_id != 0) return;

    Performance perf;
    {
        //        Create Equation Systems
        perf.start();

        string transformation = "createequationsystems";
        Task t(simulationID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("FINISHED");

        char* element = (char*) malloc(4096);
        sprintf(element, "%d;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f",
                simulationID, Reynolds, Gr, Sc, Us, Diffusivity, xlock, fopc, theta, ex, ey, ez, c_factor);
        vector<string> e = {element};
        t.addSet("o" + transformation, e);
        free(element);

        PerformanceMetric p;
        p.SetDescription("libMeshSedimentation::" + transformation);
        p.SetMethod("COMPUTATION");
        t.addPerformance(p);
        //    missing starttime and endtime

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Output" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }

    {
        perf.start();

        string transformation = "getmaximumiterations";
        Task t(simulationID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("RUNNING");
        t.addDtDependency("createequationsystems");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        PerformanceMetric p;
        p.SetDescription("libMeshSedimentation::" + transformation);
        p.SetMethod("COMPUTATION");
        t.addPerformance(p);
        //    missing starttime and endtime

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Input" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::outputGetMaximumIterations(int simulationID, Real dt, Real tmax,
        unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance,
        int max_linear_iters, int max_r_steps, unsigned int write_interval, string xdmf) {
    if (processor_id != 0) return;

    Performance perf;
    {
        perf.start();

        string transformation = "getmaximumiterations";
        Task t(simulationID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("FINISHED");

        char* element = (char*) malloc(4096);
        sprintf(element, "%d;%.2f;%.2f;%d;%d;%.2f;%d;%d;%d;%s/%s",
                simulationID, dt, tmax, n_time_steps, n_nonlinear_steps,
                nonlinear_tolerance, max_linear_iters, max_r_steps,
                write_interval, directory.c_str(), xdmf.c_str());
        vector<string> e = {element};
        t.addSet("o" + transformation, e);
        free(element);

        PerformanceMetric p;
        p.SetDescription("libMeshSedimentation::" + transformation);
        p.SetMethod("COMPUTATION");
        t.addPerformance(p);
        //    missing starttime and endtime

        File f(directory, xdmf);
        t.addFile(f);

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Output" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::inputInitDataExtraction(int simulationID, string transformation, string extractionFileName) {
    if (processor_id != 0) return;
    Performance perf;
    {
        perf.start();

        Task t(simulationID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("RUNNING");
        t.addDtDependency("getmaximumiterations");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Input" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::outputInitDataExtraction(int simulationID, string transformation, string extractionFileName, string outDataSet, int time_step, string xdmf, string rawDataFile) {
    if (processor_id != 0) return;
    Performance perf;
    {
        perf.start();

        Task t(simulationID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("FINISHED");
        t.addDtDependency("getmaximumiterations");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        char* element = (char*) malloc(4096);
        sprintf(element, "%d;%d;%s/%s;%s/%s",
                simulationID, time_step, directory.c_str(), xdmf.c_str(),
                directory.c_str(), rawDataFile.c_str());
        vector<string> e = {element};
        t.addSet("o" + transformation, e);
        free(element);

        PerformanceMetric p;
        char* perfbuffer = (char*) malloc(4096);
        sprintf(perfbuffer, "libMeshSedimentation::%s-%d", transformation.c_str(), simulationID);
        p.SetDescription(perfbuffer);
        p.SetMethod("COMPUTATION");
        t.addPerformance(p);
        //    missing starttime and endtime

        File f1(directory, xdmf);
        t.addFile(f1);
        File f2(directory, rawDataFile);
        t.addFile(f2);

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Output" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::inputSolverSimulationFluid(int taskID, int simulationID, int subTaskID) {
    if (processor_id != 0) return;
    Performance perf;
    {
        perf.start();

        string transformation = "solversimulationfluid";
        Task t(taskID);
        t.setSubID(subTaskID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("RUNNING");
        t.addDtDependency("getmaximumiterations");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Input" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::outputSolverSimulationFluid(int taskID, int simulationID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations,
        Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged) {
    if (processor_id != 0) return;
    Performance perf;
    {
        perf.start();

        string transformation = "solversimulationfluid";
        Task t(taskID);
        t.setSubID(subTaskID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("FINISHED");
        t.addDtDependency("getmaximumiterations");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        char* element = (char*) malloc(4096);
        sprintf(element, "%d;%d;%.2f;%d;%d;%d;%.2f;%.2f;%.2f;%s",
                simulationID, time_step, time, linear_step, n_linear_step,
                n_linear_iterations, linear_residual, norm_delta,
                norm_delta_u, converged ? "true" : "false");
        vector<string> e = {element};
        t.addSet("o" + transformation, e);
        free(element);

        PerformanceMetric p;
        char* perfbuffer = (char*) malloc(4096);
        sprintf(perfbuffer, "libMeshSedimentation::%s-%d-%d",
                transformation.c_str(), simulationID, subTaskID);
        p.SetDescription(perfbuffer);
        p.SetMethod("COMPUTATION");
        t.addPerformance(p);
        //    missing starttime and endtime

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Output" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::inputSolverSimulationSediments(int taskID, int simulationID, int subTaskID) {
    if (processor_id != 0) return;
    Performance perf;
    {
        perf.start();

        string transformation = "solversimulationsediments";
        Task t(taskID);
        t.setSubID(subTaskID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("RUNNING");
        t.addDtDependency("solversimulationfluid");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Input" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::outputSolverSimulationSediments(int taskID, int simulationID, int subTaskID, int time_step,
        Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations,
        Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged) {
    if (processor_id != 0) return;
    Performance perf;
    {
        perf.start();

        string transformation = "solversimulationsediments";
        Task t(taskID);
        t.setSubID(subTaskID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("FINISHED");
        t.addDtDependency("solversimulationfluid");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        char* element = (char*) malloc(4096);
        sprintf(element, "%d;%d;%.2f;%d;%d;%d;%.2f;%.2f;%.2f;%s",
                simulationID, time_step, time, linear_step, n_linear_step,
                n_linear_iterations, linear_residual, norm_delta, norm_delta_u,
                converged ? "true" : "false");
        vector<string> e = {element};
        t.addSet("o" + transformation, e);
        free(element);

        PerformanceMetric p;
        char* perfbuffer = (char*) malloc(4096);
        sprintf(perfbuffer, "libMeshSedimentation::%s-%d-%d",
                transformation.c_str(), simulationID, subTaskID);
        p.SetDescription(perfbuffer);
        p.SetMethod("COMPUTATION");
        t.addPerformance(p);
        //    missing starttime and endtime

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Output" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::outputMeshRefinement(int taskID, int simulationID, int subTaskID,
        bool first_step_refinement, int time_step, int before_n_active_elem, int after_n_active_elem) {
    if (processor_id != 0) return;
    Performance perf;
    {
        perf.start();

        string transformation = "meshrefinement";
        Task t(taskID);
        t.setSubID(subTaskID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("FINISHED");
        t.addDtDependency("solversimulationsediments");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        char* element = (char*) malloc(4096);
        sprintf(element, "%d;%s;%d;%d;%d",
                simulationID, first_step_refinement ? "true" : "false", time_step,
                before_n_active_elem, after_n_active_elem);
        vector<string> e = {element};
        t.addSet("o" + transformation, e);
        free(element);

        PerformanceMetric p;
        char* perfbuffer = (char*) malloc(4096);
        sprintf(perfbuffer, "libMeshSedimentation::%s-%d-%d",
                transformation.c_str(), simulationID, subTaskID);
        p.SetDescription(perfbuffer);
        p.SetMethod("COMPUTATION");
        t.addPerformance(p);
        //    missing starttime and endtime

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Output" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::inputMeshWriter(int taskID, int simulationID, int subTaskID) {
    if (processor_id != 0) return;
    Performance perf;
    {
        perf.start();

        string transformation = "meshwriter";
        Task t(taskID);
        t.setSubID(subTaskID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("RUNNING");
        t.addDtDependency("solversimulationsediments");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Input" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::outputMeshWriter(int taskID, int simulationID, int subTaskID, int time_step, string xdmf) {
    if (processor_id != 0) return;
    Performance perf;
    {
        perf.start();

        string transformation = "meshwriter";
        Task t(taskID);
        t.setSubID(subTaskID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("FINISHED");
        t.addDtDependency("solversimulationsediments");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        File f1(directory, xdmf);
        t.addFile(f1);

        char* element = (char*) malloc(4096);
        sprintf(element, "%d;%d;%s/%s",
                simulationID, time_step, directory.c_str(), xdmf.c_str());
        vector<string> e = {element};
        t.addSet("o" + transformation, e);
        free(element);

        PerformanceMetric p;
        char* perfbuffer = (char*) malloc(4096);
        sprintf(perfbuffer, "libMeshSedimentation::%s-%d-%d",
                transformation.c_str(), simulationID, subTaskID);
        p.SetDescription(perfbuffer);
        p.SetMethod("COMPUTATION");
        t.addPerformance(p);
        //    missing starttime and endtime

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Output" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::inputDataExtraction(int taskID, int simulationID, int subTaskID, string transformation, string extractionFileName) {
    if (processor_id != 0) return;
    Performance perf;
    {
        perf.start();

        Task t(taskID);
        t.setSubID(subTaskID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("RUNNING");
        t.addDtDependency("meshwriter");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Input" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::outputDataExtraction(int taskID, int simulationID, int subTaskID, string transformation, string extractionFileName, string outDataSet, int time_step, string xdmf, string rawDataFile) {
    if (processor_id != 0) return;
    Performance perf;
    {
        perf.start();

        Task t(taskID);
        t.setSubID(subTaskID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("FINISHED");
        t.addDtDependency("meshwriter");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        File f1(directory, xdmf);
        t.addFile(f1);
        File f2(directory, rawDataFile);
        t.addFile(f2);

        char* element = (char*) malloc(4096);
        sprintf(element, "%d;%d;%s/%s;%s/%s",
                simulationID, time_step, directory.c_str(), xdmf.c_str(),
                directory.c_str(), rawDataFile.c_str());
        vector<string> e = {element};
        t.addSet("o" + transformation, e);
        free(element);

        PerformanceMetric p;
        char* perfbuffer = (char*) malloc(4096);
        sprintf(perfbuffer, "libMeshSedimentation::%s-%d-%d",
                transformation.c_str(), simulationID, subTaskID);
        p.SetDescription(perfbuffer);
        p.SetMethod("COMPUTATION");
        t.addPerformance(p);
        //    missing starttime and endtime

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Output" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
}

void Provenance::meshAggregator(int simulationID, string xdmf, int n_processors, string meshDependencies) {
    if (processor_id != 0) return;
    Performance perf;
    {
        perf.start();

        string transformation = "meshaggregator";
        Task t(simulationID);
        t.setDataflow(dataflow);
        t.setTransformation(transformation);
        t.setWorkspace(directory);
        t.setStatus("FINISHED");
        t.addDtDependency("meshwriter");
        char* vs = (char*) malloc(4096);
        sprintf(vs, "%d", simulationID);
        t.addIdDependency(vs);
        free(vs);

        File f1(directory, xdmf);
        t.addFile(f1);

        char* element = (char*) malloc(4096);
        sprintf(element, "%d;%s/%s;%d",
                simulationID, directory.c_str(), xdmf.c_str(), n_processors);
        vector<string> e = {element};
        t.addSet("o" + transformation, e);
        free(element);

        PerformanceMetric p;
        char* perfbuffer = (char*) malloc(4096);
        sprintf(perfbuffer, "libMeshSedimentation::%s-%d",
                transformation.c_str(), simulationID);
        p.SetDescription(perfbuffer);
        p.SetMethod("COMPUTATION");
        t.addPerformance(p);
        //    missing starttime and endtime

        char* buffer = (char*) malloc(4096);
        sprintf(buffer, "%s%s-%d-F.json", 
                jsonDirectory.c_str(), transformation.c_str(), simulationID);
        t.writeJSON(buffer);
        free(buffer);

        perf.end();
        double elapsedTime = perf.elapsedTime();

        char* bPointer = (char*) malloc(512);
        ofstream file;
        file.open("prov/log/" + transformation + ".prov", ios_base::app);
        file << "PROV:" + transformation + ":Output" << endl;
        sprintf(bPointer, "%.2f", elapsedTime);
        file << space << *bPointer << endl;
        file << space << "elapsed-time: " << *bPointer << " seconds." << endl;
        file.close();
        free(bPointer);
    }
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



