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
#include <ctime>

#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"

#include "dfanalyzer/task.h"
#include "dfanalyzer/indexer.h"

#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"

#include "provenance.h"
#include "performance.h"

#define LINUX
#define BACKUP

using namespace std;
using namespace libMesh;

Provenance::Provenance() {
    GetPot infile("provenance.in");
    directory = infile("directory", "/Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/sedimentation");

    string pgFilePath = infile("pgFilePath", "/Users/vitor/Documents/Repository/Thesis/Workflow-Sedimentation/dfa/PG-1.0.jar");
    string rdeFilePath = infile("rdeFilePath", "/Users/vitor/Documents/Repository/Thesis/Workflow-Sedimentation/dfa/RDE-1.0.jar");
    string rdiFilePath = infile("rdiFilePath", "/Users/vitor/Documents/Repository/Thesis/Workflow-Sedimentation/dfa/RDI-1.0.jar");
    bin = infile("bin", "");
    extraArguments = infile("extraArguments", "");
    rawDataAccess = infile("access", "EXTRACTION");
    cartridge = infile("cartridge", "CSV");
    pgCommandLine = "java -jar " + pgFilePath + " ";
    rdeCommandLine = "java -jar " + rdeFilePath + " ";
    rdiCommandLine = "java -jar " + rdiFilePath + " ";

    jsonDirectory = directory + "/prov/di/" + dataflow + "/";
    pgDirectory = directory + "/prov/pg/" + dataflow + "/";
#ifdef LINUX
    directory = directory.substr(0, directory.size() - 1);
    jsonDirectory = directory + "/prov/di/" + dataflow + "/";
    pgDirectory = directory + "/prov/pg/" + dataflow + "/";
#endif

    cout << "###########################" << endl;
    cout << "Provenance Properties" << endl;
    cout << "pgFilePath=";
    cout << pgFilePath.c_str() << endl;
    cout << "rdeFilePath=";
    cout << rdeFilePath.c_str() << endl;
    cout << "rdiFilePath=";
    cout << rdiFilePath.c_str() << endl;
    cout << "rawDataAccess=";
    cout << rawDataAccess.c_str() << endl;
    cout << "cartridge=";
    cout << cartridge.c_str() << endl;
    cout << "bin=";
    cout << bin.c_str() << endl;
    cout << "extraArguments=";
    cout << extraArguments.c_str() << endl;
    cout << "jsonDirectory=";
    cout << jsonDirectory.c_str() << endl;
    cout << "pgDirectory=";
    cout << pgDirectory.c_str() << endl;
    cout << "###########################" << endl;

    processor_id = libMesh::global_processor_id();
    simulationID = 1;
}

void Provenance::inputInputMesh(int dim) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Input Mesh" << endl;
#endif    
    Performance perf;
    perf.begin();

    string transformation = "inputmesh";
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
    sprintf(memalloc, "%d;%d", simulationID, dim);
    vector<string> e = {memalloc};
    t.addSet("i" + transformation, e);

    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputInputMesh(double r_fraction, double c_fraction,
        double max_h_level, unsigned int hlevels) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Input Mesh" << endl;
#endif

    Performance perf;
    perf.begin();

    string transformation = "inputmesh";
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
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();

    //    Create Equation Systems
    perf.begin();

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
    t2.addDtDependency("inputmesh");

    sprintf(memalloc, "%d", simulationID);
    t2.addIdDependency(memalloc);

    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif

    perf.end();
    elapsedTime = perf.getElapsedTime();

    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputCreateEquationSystems(Real Reynolds, Real Gr,
        Real Sc, Real Us, Real Diffusivity, Real xlock, Real fopc,
        Real theta, Real ex, Real ey, Real ez, Real c_factor) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Create Equation Systems" << endl;
#endif

    Performance perf;
    //        Create Equation Systems
    perf.begin();

    char memalloc[jsonArraySize];
    string transformation = "createequationsystems";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("inputmesh");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

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
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();

    perf.begin();

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
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif

    perf.end();
    elapsedTime = perf.getElapsedTime();

    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputGetMaximumIterations(Real dt, Real tmax,
        unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance,
        int max_linear_iters, int max_r_steps, unsigned int write_interval, string xdmf) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Get Maximum Iterations" << endl;
#endif

    Performance perf;
    perf.begin();

    char memalloc[1000];
    string transformation = "getmaximumiterations";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("createequationsystems");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%.7f;%.7f;%d;%d;%.9f;%d;%d;%d;%s",
            simulationID, dt, tmax, n_time_steps, n_nonlinear_steps,
            nonlinear_tolerance, max_linear_iters, max_r_steps,
            write_interval, xdmf.c_str());
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
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::inputInitDataExtraction(int lineID) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Init Data Extraction" << endl;
#endif

    char transformation[arraySize];
    if (lineID == -1) {
        sprintf(transformation, "initdataextraction");
    } else {
        sprintf(transformation, "iline%dextraction", lineID);
    }

    Performance perf;
    perf.begin();

    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d", transformation, simulationID);
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

    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation, simulationID);
    t.writeJSON(memalloc);
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation, simulationID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + string(transformation) + ".prov", ios_base::app);
    file << "PROV:" << transformation << ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputInitDataExtraction(int lineID, string xdmf, int dimension, int indexerID) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Init Data Extraction" << endl;
#endif

    char transformationTag[arraySize];
    char setName[arraySize];
    char extractorName[arraySize];
    char rawDataFile[arraySize];

    if (lineID == -1) {
        sprintf(transformationTag, "initdataextraction");
        sprintf(setName, "oinitdataextraction");
        sprintf(extractorName, "irde");
        sprintf(rawDataFile, "ext_line_0.csv");
    } else {
        sprintf(transformationTag, "iline%dextraction", lineID);
        sprintf(setName, "oline%diextraction", lineID);
        sprintf(extractorName, "iline%d", lineID);
        sprintf(rawDataFile, "ext_line_%d_0.csv", lineID);
    }

    Performance perf;
    perf.begin();

    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformationTag);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("getmaximumiterations");

    char memalloc[4096];
    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + string(transformationTag) + ".prov", ios_base::app);
    file << "PROV:" + string(transformationTag) + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();

    Performance rdiPerf;
    rdiPerf.begin();

    string extension = "data";
    if (rawDataAccess.compare("INDEXING") == 0) {
        extension = "index";
        Indexer idx(rdiCommandLine, rawDataAccess, cartridge, "irde");
        idx.addAttribute("u", "numeric", false);
        idx.addAttribute("v", "numeric", false);
        if (dimension == 3) {
            idx.addAttribute("w", "numeric", false);
        }
        idx.addAttribute("p", "numeric", false);
        idx.addAttribute("s", "numeric", false);
        idx.addAttribute("d", "numeric", false);
        if (dimension == 3) {
            idx.addAttribute("vtkvalidpointmask", "numeric", false);
            idx.addAttribute("arc_length", "numeric", false);
        }
        idx.addAttribute("points0", "numeric", false);
        idx.addAttribute("points1", "numeric", false);
        idx.addAttribute("points2", "numeric", false);

        if (cartridge.compare("FASTBIT") == 0 || cartridge.compare("OPTIMIZED_FASTBIT") == 0) {
            idx.setBin(bin);
            idx.setExtraArguments(extraArguments);
        }

        idx.index(directory, rawDataFile, indexerID);
    }

    rdiPerf.end();
    storeRDIComponentCost(0, 1, rdiPerf.getElapsedTime());

    perf.begin();

    char extractedFileName[jsonArraySize];
    if (rawDataAccess.compare("INDEXING") == 0) {
        sprintf(extractedFileName, "%s.%s", "irde", extension.c_str());
    } else {
        sprintf(extractedFileName, "%s", rawDataFile);
    }

    if (rawDataAccess.compare("INDEXING") == 0 && (cartridge.compare("FASTBIT") == 0 || cartridge.compare("OPTIMIZED_FASTBIT") == 0)) {
        sprintf(memalloc, "%d;%d;%s;%s/index/%d/%s",
                simulationID, 0, xdmf.c_str(),
                directory.c_str(), indexerID, extractedFileName);
    } else {
        sprintf(memalloc, "%d;%d;%s;%s/%s",
                simulationID, 0, xdmf.c_str(),
                directory.c_str(), extractedFileName);
    }
    cout << memalloc << endl;

    vector<string> e = {memalloc};
    t.addSet(setName, e);

    PerformanceMetric p;
    sprintf(memalloc, "libMeshSedimentation::%s-%d", transformationTag, simulationID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    File f1(directory, xdmf);
    t.addFile(f1);

    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformationTag, simulationID);
    t.writeJSON(memalloc);
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformationTag, simulationID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    elapsedTime = perf.getElapsedTime();

    file.open("prov/log/" + string(transformationTag) + ".prov", ios_base::app);
    file << "PROV:" + string(transformationTag) + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::inputInitVisualization(int lineID) {
    if (processor_id != 0 or lineID != 0) return;
#ifdef VERBOSE
    cout << "Input Init Visualization" << endl;
#endif

    string transformation = "ivisualization";

    Performance perf;
    perf.begin();

    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d", transformation, simulationID);
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

    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation, simulationID);
    t.writeJSON(memalloc);
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputInitVisualization(int lineID, int timeStep) {
    if (processor_id != 0 or lineID != 0) return;
#ifdef VERBOSE
    cout << "Output Init Visualization" << endl;
#endif

    char transformation[arraySize];
    char dataSet[arraySize];
    char png[arraySize];
    sprintf(transformation, "ivisualization");
    sprintf(dataSet, "oivisualization");
    sprintf(png, "image_%d.png", timeStep);

    Performance perf;
    perf.begin();

    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("getmaximumiterations");

    char memalloc[4096];
    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    File f1(directory, png);
    t.addFile(f1);

    sprintf(memalloc, "%d;%d;%s%s",
            simulationID, timeStep, pgDirectory.c_str(), png);

    vector<string> e = {memalloc};
    t.addSet(dataSet, e);

    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation, simulationID);
    t.writeJSON(memalloc);
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation, simulationID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + string(transformation) + ".prov", ios_base::app);
    file << "PROV:" + string(transformation) + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::inputVisualization(int lineID, int taskID) {
    if (processor_id != 0 or lineID != 0) return;
#ifdef VERBOSE
    cout << "Input Visualization" << endl;
#endif

    char transformation[arraySize];
    sprintf(transformation, "visualization");

    Performance perf;
    perf.begin();

    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d", transformation, simulationID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(taskID);
    t.addPerformanceMetric(p);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("RUNNING");
    t.addDtDependency("meshwriter");

    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation, taskID);
    t.writeJSON(memalloc);
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation, taskID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + string(transformation) + ".prov", ios_base::app);
    file << "PROV:" + string(transformation) + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputVisualization(int lineID, int taskID, int timeStep) {
    if (processor_id != 0 or lineID != 0) return;
#ifdef VERBOSE
    cout << "Output Visualization" << endl;
#endif

    char transformation[arraySize];
    char dataSet[arraySize];
    char png[arraySize];
    sprintf(transformation, "visualization");
    sprintf(dataSet, "ovisualization");
    sprintf(png, "image_%d.png", timeStep);

    Performance perf;
    perf.begin();

    Task t(taskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("meshwriter");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

    File f1(directory, png);
    t.addFile(f1);

    sprintf(memalloc, "%d;%d;%s%s",
            simulationID, timeStep, pgDirectory.c_str(), png);

    vector<string> e = {memalloc};
    t.addSet(dataSet, e);

    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation, taskID);
    t.writeJSON(memalloc);
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation, taskID);
    t.writeJSON(memalloc);
#endif        

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + string(transformation) + ".prov", ios_base::app);
    file << "PROV:" + string(transformation) + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::inputSolverSimulationFluid(int taskID, int subTaskID) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Solver Simulation Fluid" << endl;
#endif

    Performance perf;
    perf.begin();

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
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputSolverSimulationFluid(int taskID, int subTaskID,
        int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations,
        Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Solver Simulation Fluid" << endl;
#endif

    Performance perf;
    perf.begin();

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

    sprintf(memalloc, "%d;%d;%.7f;%d;%d;%d;%.9f;%.9f;%.9f;%s",
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
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::inputSolverSimulationSediments(int taskID, int subTaskID) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Solver Simulation Sediments" << endl;
#endif
    Performance perf;
    perf.begin();

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

    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputSolverSimulationSediments(int taskID, int subTaskID, int time_step,
        Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations,
        Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Solver Simulation Sediments" << endl;
#endif

    Performance perf;
    perf.begin();

    string transformation = "solversimulationsediments";
    Task t(taskID);
    t.setSubID(subTaskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("solversimulationfluid");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%d;%.7f;%d;%d;%d;%.9f;%.9f;%.9f;%s",
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
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputMeshRefinement(int taskID, int subTaskID,
        bool first_step_refinement, int time_step, int before_n_active_elem, int after_n_active_elem) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Mesh Refinement" << endl;
#endif
    Performance perf;
    perf.begin();

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

    sprintf(memalloc, "%d", taskID);
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
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::inputMeshWriter(int taskID, int subTaskID) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Mesh Writer" << endl;
#endif

    Performance perf;
    perf.begin();

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
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputMeshWriter(int taskID, int subTaskID, int time_step, string xdmf) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Mesh Writer" << endl;
#endif

    Performance perf;
    perf.begin();

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

    sprintf(memalloc, "%d;%d;%s",
            simulationID, time_step, xdmf.c_str());

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
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::inputDataExtraction(int taskID, int subTaskID, int lineID) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Data Extraction" << endl;
#endif

    char transformation[arraySize];
    if (lineID == -1) {
        sprintf(transformation, "dataextraction");
    } else {
        sprintf(transformation, "line%dextraction", lineID);
    }

    Performance perf;
    perf.begin();

    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation, simulationID, subTaskID);
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

    sprintf(memalloc, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation, simulationID, subTaskID);
    t.writeJSON(memalloc);
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-R.json", pgDirectory.c_str(), transformation, simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + string(transformation) + ".prov", ios_base::app);
    file << "PROV:" + string(transformation) + ":Input" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::outputDataExtraction(int taskID, int subTaskID, int lineID, int timeStep,
        string xdmf, int dimension, int indexerID) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Data Extraction" << endl;
#endif

    char transformation[arraySize];
    char dataSet[arraySize];
    char indexerName[arraySize];
    char rawDataFile[arraySize];
    if (lineID == -1) {
        sprintf(transformation, "dataextraction");
        sprintf(transformation, "odataextraction");
        sprintf(transformation, "rde%d", subTaskID);
        sprintf(rawDataFile, "ext_line_%d.csv", timeStep);
    } else {
        sprintf(transformation, "line%dextraction", lineID);
        sprintf(dataSet, "oline%dextraction", lineID);
        sprintf(indexerName, "line%d%d", lineID, subTaskID);
        sprintf(rawDataFile, "ext_line_%d_%d.csv", lineID, timeStep);
    }

    Performance perf;
    perf.begin();

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
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + string(transformation) + ".prov", ios_base::app);
    file << "PROV:" + string(transformation) + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();

    Performance rdiPerf;
    rdiPerf.begin();

    string extension = "data";
    if (rawDataAccess.compare("INDEXING") == 0) {
        extension = "index";
        Indexer idx(rdiCommandLine, rawDataAccess, cartridge, indexerName);
        idx.addAttribute("u", "numeric", false);
        idx.addAttribute("v", "numeric", false);
        if (dimension == 3) {
            idx.addAttribute("w", "numeric", false);
        }
        idx.addAttribute("p", "numeric", false);
        idx.addAttribute("s", "numeric", false);
        idx.addAttribute("d", "numeric", false);
        if (dimension == 3) {
            idx.addAttribute("vtkvalidpointmask", "numeric", false);
            idx.addAttribute("arc_length", "numeric", false);
        }
        idx.addAttribute("points0", "numeric", false);
        idx.addAttribute("points1", "numeric", false);
        idx.addAttribute("points2", "numeric", false);
        if (cartridge.compare("FASTBIT") == 0 || cartridge.compare("OPTIMIZED_FASTBIT") == 0) {
            idx.setBin(bin);
            idx.setExtraArguments(extraArguments);
        }
        idx.index(directory, rawDataFile, indexerID);
    }

    rdiPerf.end();
    storeRDIComponentCost(taskID, subTaskID, rdiPerf.getElapsedTime());

    perf.begin();

    char extractedFileName[jsonArraySize];
    if (rawDataAccess.compare("INDEXING") == 0) {
        sprintf(extractedFileName, "%s.%s", indexerName, extension.c_str());
    } else {
        sprintf(extractedFileName, "%s", rawDataFile);
    }

    if (rawDataAccess.compare("INDEXING") == 0 && (cartridge.compare("FASTBIT") == 0 || cartridge.compare("OPTIMIZED_FASTBIT") == 0)) {
        sprintf(memalloc, "%d;%d;%s;%s/index/%d/%s",
                simulationID, timeStep, xdmf.c_str(),
                directory.c_str(), indexerID, extractedFileName);
    } else {
        sprintf(memalloc, "%d;%d;%s;%s/%s",
                simulationID, timeStep, xdmf.c_str(),
                directory.c_str(), extractedFileName);
    }
    cout << memalloc << endl;

    vector<string> e = {memalloc};
    t.addSet(dataSet, e);

    File f1(directory, xdmf);
    t.addFile(f1);

    PerformanceMetric p;
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation, simulationID, subTaskID);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation, simulationID, subTaskID);
    t.writeJSON(memalloc);
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-F.json", pgDirectory.c_str(), transformation, simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    elapsedTime = perf.getElapsedTime();

    file.open("prov/log/" + string(transformation) + ".prov", ios_base::app);
    file << "PROV:" + string(transformation) + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::meshAggregator(string xdmf, int n_processors, vector<string> meshDependencies) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Mesh Writer" << endl;
#endif

    Performance perf;
    perf.begin();

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

    // with visualization
    //    sprintf(memalloc, "%d;%s%s;%d;%svideo.mp4",
    //            simulationID, pgDirectory.c_str(), xdmf.c_str(), n_processors, pgDirectory.c_str());
    // without visualization
    sprintf(memalloc, "%d;%s%s;%d",
            simulationID, pgDirectory.c_str(), xdmf.c_str(), n_processors);


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
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif

    perf.end();
    double elapsedTime = perf.getElapsedTime();

    ofstream file;
    file.open("prov/log/" + transformation + ".prov", ios_base::app);
    file << "PROV:" + transformation + ":Output" << endl;
    sprintf(memalloc, "%.5f", elapsedTime);
    file << space << memalloc << endl;
    file << space << "elapsed-time: " << memalloc << " seconds." << endl;
    file.close();
}

void Provenance::storeCatalystCost(int taskID, int subTaskID, double elapsedTime) {
    if (processor_id != 0) return;
    ofstream file;
    file.open("prov/paraview/catalyst_" + to_string(taskID) + "_" + to_string(subTaskID) + ".prov", ios_base::app);
    file << "Paraview:Catalyst:Pipeline" << endl;
    char buffer[256];
    sprintf(buffer, "%.5f", elapsedTime);
    file << space << "elapsed-time: " << buffer << " seconds." << endl;
    file.close();

    sprintf(buffer, "Catalyst Cost: %.5f", elapsedTime);
    cout << buffer << endl;
}

void Provenance::storeRDIComponentCost(int taskID, int subTaskID, double elapsedTime) {
    if (processor_id != 0) return;
    ofstream file;
    file.open("prov/rdi/indexing_" + to_string(taskID) + "_" + to_string(subTaskID) + ".prov", ios_base::app);
    file << "RDIComponent:DataExtraction:Process" << endl;
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
    int exitStatus = system(strdup(str.c_str()));

    cout << "[Provenance] Finish Data Ingestor" << endl;
}

void Provenance::createIndexDirectory() {
    char cmd[32];
    sprintf(cmd, "mkdir index");
    int exitStatus = system(cmd);
};

