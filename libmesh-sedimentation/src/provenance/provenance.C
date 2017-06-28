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
#include <unistd.h>
#define GetCurrentDir getcwd

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
#define DATABASE
#define BACKUP

using namespace std;
using namespace libMesh;

Provenance::Provenance(int processorID) {
    processor_id = processorID;
    simulationID = 1;
}

void Provenance::SetUp() {

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
}

std::string GetCurrentWorkingDir(void) {
    char buff[FILENAME_MAX];
    GetCurrentDir(buff, FILENAME_MAX);
    std::string current_working_dir(buff);
    return current_working_dir;
}

void Provenance::inputInputMesh() {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Input Mesh" << endl;
#endif    
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
    sprintf(memalloc, "%d", simulationID);
    vector<string> e = {memalloc};
    t.addSet("i" + transformation, e);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif
}

void Provenance::outputInputMesh(int dim, string mesh_file, bool restartControl) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Input Mesh" << endl;
#endif
    string transformation = "inputmesh";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d;%d;%s/%s;%s",
            simulationID, dim, GetCurrentWorkingDir().c_str(), mesh_file.c_str(), restartControl ? "true" : "false");
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif

    //    AMR Config
    transformation = "amrconfig";
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

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif
}

void Provenance::outputAMRConfig(double r_fraction, double c_fraction, double max_h_level, unsigned int hlevels, bool first_step_refinement,
        bool amrc_flow_transp, int ref_interval, int max_r_steps) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output AMR Config" << endl;
#endif

    string transformation = "amrconfig";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("inputmesh");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%.2f;%.2f;%.2f;%d;%s;%s;%d;%d",
            simulationID, r_fraction, c_fraction, max_h_level, hlevels,
            first_step_refinement ? "true" : "false",
            amrc_flow_transp ? "true" : "false",
            ref_interval, max_r_steps);
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif

    //    Create Equation Systems
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

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif
}

void Provenance::outputCreateEquationSystems(Real Reynolds, Real Gr,
        Real Sc, Real Us, Real Diffusivity, Real xlock, Real fopc,
        Real theta, Real ex, Real ey, Real ez, Real c_factor) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Create Equation Systems" << endl;
#endif

    //        Create Equation Systems
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

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif

    //    Time Step Control Config
    transformation = "timestepcontrolconfig";
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

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif
}

void Provenance::outputTSControlConfig(string ts_control_model_name, double dt_min, double dt_max, double tol_u, double tol_s,
        double kp, double ki, double kd, unsigned int nsa_max, unsigned int nsa_target_flow, unsigned int nsa_target_transport,
        unsigned int nsa_limit_flow, unsigned int nsa_limit_transport, double mult_factor_max, double mult_factor_min,
        double pc11_theta, double alpha, double k_exp, double s_min, double s_max, double reduct_factor, bool complete_flow_norm) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Time Step Control Config" << endl;
#endif

    char memalloc[jsonArraySize];
    string transformation = "timestepcontrolconfig";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("inputmesh");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%s;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%d;%d;%d;%d;%d;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%s",
            simulationID, ts_control_model_name.c_str(), dt_min, dt_max, tol_u, tol_s, kp,
            ki, kd, nsa_max, nsa_target_flow, nsa_target_transport, nsa_limit_flow,
            nsa_limit_transport, mult_factor_max, mult_factor_min, pc11_theta,
            alpha, k_exp, s_min, s_max, reduct_factor, complete_flow_norm ? "true" : "false");

    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif

    //    IO Config
    transformation = "ioconfig";
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

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif
}

void Provenance::outputIOConfig(string dpath, string rname, unsigned int write_interval, unsigned int catalyst_interval, bool write_restart) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output IO Config" << endl;
#endif

    char memalloc[jsonArraySize];
    string transformation = "ioconfig";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("inputmesh");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%s;%s;%d;%d;%s",
            simulationID, dpath.c_str(), rname.c_str(), write_interval, catalyst_interval, write_restart ? "true" : "false");

    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif

    //    Get Maximum Iterations
    transformation = "getmaximumiterationstoflow";
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

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif
}

void Provenance::outputGetMaximumIterationsToFlow(Real dt, Real tmax,
        unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance,
        int max_linear_iters, string xdmf) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Get Maximum Iterations to Flow" << endl;
#endif

    char memalloc[1000];
    string transformation = "getmaximumiterationstoflow";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("createequationsystems");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%.7f;%.7f;%d;%d;%.9f;%d;%s",
            simulationID, dt, tmax, n_time_steps, n_nonlinear_steps,
            nonlinear_tolerance, max_linear_iters, xdmf.c_str());
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    File f(directory, xdmf);
    t.addFile(f);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif

    //    Get Maximum Iterations to Transport
    transformation = "getmaximumiterationstotransport";
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t2(simulationID);
    t2.addPerformanceMetric(p);
    t2.setDataflow(dataflow);
    t2.setTransformation(transformation);
    t2.setWorkspace(directory);
    t2.setStatus("RUNNING");
    t2.addDtDependency("getmaximumiterationstoflow");

    sprintf(memalloc, "%d", simulationID);
    t2.addIdDependency(memalloc);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t2.writeJSON(memalloc);
#endif
}

void Provenance::outputGetMaximumIterationsToTransport(Real dt, Real tmax,
        unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance,
        int max_linear_iters, string xdmf) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Get Maximum Iterations to Transport" << endl;
#endif

    char memalloc[1000];
    string transformation = "getmaximumiterationstotransport";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("getmaximumiterationstoflow");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%.7f;%.7f;%d;%d;%.9f;%d;%s",
            simulationID, dt, tmax, n_time_steps, n_nonlinear_steps,
            nonlinear_tolerance, max_linear_iters, xdmf.c_str());
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    p.SetDescription("libMeshSedimentation::" + transformation);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

    File f(directory, xdmf);
    t.addFile(f);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif
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
    t.addDtDependency("getmaximumiterationstotransport");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation, simulationID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation, simulationID);
    t.writeJSON(memalloc);
#endif
}

void Provenance::outputInitDataExtraction(int lineID, string xdmf, int dimension) {
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

    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformationTag);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("getmaximumiterationstotransport");

    char memalloc[4096];
    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

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

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformationTag, simulationID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformationTag, simulationID);
    t.writeJSON(memalloc);
#endif
}

void Provenance::inputInitVisualization(int lineID) {
    if (processor_id != 0 or lineID != 0) return;
#ifdef VERBOSE
    cout << "Input Init Visualization" << endl;
#endif

    char transformation[arraySize];
    sprintf(transformation, "ivisualization");

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
    t.addDtDependency("getmaximumiterationstotransport");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation, simulationID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation, simulationID);
    t.writeJSON(memalloc);
#endif
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

    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("getmaximumiterationstotransport");

    char memalloc[4096];
    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    File f1(directory, png);
    t.addFile(f1);

    sprintf(memalloc, "%d;%d;%s%s",
            simulationID, timeStep, pgDirectory.c_str(), png);

    vector<string> e = {memalloc};
    t.addSet(dataSet, e);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation, simulationID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation, simulationID);
    t.writeJSON(memalloc);
#endif
}

void Provenance::inputVisualization(int lineID) {
    if (processor_id != 0 or lineID != 0) return;
#ifdef VERBOSE
    cout << "Input Visualization" << endl;
#endif

    char transformation[arraySize];
    sprintf(transformation, "visualization");

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
    t.addDtDependency("solversimulationtransport");

    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-R.json", jsonDirectory.c_str(), transformation, taskID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-R.json", pgDirectory.c_str(), transformation, taskID);
    t.writeJSON(memalloc);
#endif
}

void Provenance::outputVisualization(int lineID, int timeStep) {
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

    Task t(taskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("solversimulationtransport");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

    File f1(directory, png);
    t.addFile(f1);

    sprintf(memalloc, "%d;%d;%s%s",
            simulationID, timeStep, pgDirectory.c_str(), png);

    vector<string> e = {memalloc};
    t.addSet(dataSet, e);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation, taskID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation, taskID);
    t.writeJSON(memalloc);
#endif
}

void Provenance::inputSolverSimulationFlow() {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Solver Simulation Flow" << endl;
#endif

    string transformation = "solversimulationflow";
    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), simulationID, numberIterationsFlow);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(simulationID);
    t.addPerformanceMetric(p);
    t.setSubID(numberIterationsFlow);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("RUNNING");
    t.addDtDependency("getmaximumiterationstotransport");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, numberIterationsFlow);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID, numberIterationsFlow);
    t.writeJSON(memalloc);
#endif
}

void Provenance::outputSolverSimulationFlow(int time_step, Real dt, Real time,
        int linear_step, int n_linear_step, unsigned int n_linear_iterations,
        Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Solver Simulation Flow" << endl;
#endif

    string transformation = "solversimulationflow";
    Task t(simulationID);
    t.setSubID(numberIterationsFlow);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("getmaximumiterationstotransport");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%d;%.7f;%.7f;%d;%d;%d;%.9f;%.9f;%.9f;%s",
            simulationID, time_step, dt, time, linear_step, n_linear_step,
            n_linear_iterations, linear_residual, norm_delta,
            norm_delta_u, converged ? "true" : "false");
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), simulationID, numberIterationsFlow);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, numberIterationsFlow);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID, numberIterationsFlow);
    t.writeJSON(memalloc);
#endif
}

void Provenance::inputSolverSimulationTransport() {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Solver Simulation Transport" << endl;
#endif

    string transformation = "solversimulationtransport";
    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), taskID, numberIterationsTransport);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(taskID);
    t.addPerformanceMetric(p);
    t.setSubID(numberIterationsTransport);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("RUNNING");
    t.addDtDependency("solversimulationflow");

    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), taskID, numberIterationsTransport);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-R.json", pgDirectory.c_str(), transformation.c_str(), taskID, numberIterationsTransport);
    t.writeJSON(memalloc);
#endif
}

void Provenance::outputSolverSimulationTransport(int time_step, Real dt,
        Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations,
        Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Solver Simulation Transport" << endl;
#endif

    string transformation = "solversimulationtransport";
    Task t(taskID);
    t.setSubID(numberIterationsTransport);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("solversimulationflow");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d", simulationID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%d;%.7f;%.7f;%d;%d;%d;%.9f;%.9f;%.9f;%s",
            simulationID, time_step, dt, time, linear_step, n_linear_step,
            n_linear_iterations, linear_residual, norm_delta, norm_delta_u,
            converged ? "true" : "false");
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), taskID, numberIterationsTransport);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), taskID, numberIterationsTransport);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-F.json", pgDirectory.c_str(), transformation.c_str(), taskID, numberIterationsTransport);
    t.writeJSON(memalloc);
#endif
}

void Provenance::inputComputeSolutionChange() {
    incrementIterationsComputeSolutionChange();
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Compute Solution Change" << endl;
#endif

    string transformation = "computesolutionchange";
    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), taskID, numberIterationsComputeSolutionChange);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(taskID);
    t.addPerformanceMetric(p);
    t.setSubID(numberIterationsComputeSolutionChange);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("RUNNING");
    t.addDtDependency("solversimulationtransport");

    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), taskID, numberIterationsComputeSolutionChange);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-R.json", pgDirectory.c_str(), transformation.c_str(), taskID, numberIterationsComputeSolutionChange);
    t.writeJSON(memalloc);
#endif
}

void Provenance::outputComputeSolutionChange(int time_step, Real time, Real dt,
        unsigned int n_flow_linear_iterations_total, unsigned int n_flow_nonlinear_iterations_total,
        unsigned int n_transport_linear_iterations_total, unsigned int n_transport_nonlinear_iterations_total,
        bool timeStepAccepted, double error) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Compute Solution Change" << endl;
#endif

    string transformation = "computesolutionchange";
    Task t(taskID);
    t.setSubID(numberIterationsComputeSolutionChange);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("solversimulationtransport");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%d;%.7f;%.7f;%d;%d;%d;%d;%s;%.7f",
            simulationID, time_step, time, dt,
            n_flow_linear_iterations_total, n_flow_nonlinear_iterations_total,
            n_transport_linear_iterations_total, n_transport_nonlinear_iterations_total,
            timeStepAccepted ? "true" : "false",
            error
            );
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), taskID, numberIterationsComputeSolutionChange);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), taskID, numberIterationsComputeSolutionChange);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-F.json", pgDirectory.c_str(), transformation.c_str(), taskID, numberIterationsComputeSolutionChange);
    t.writeJSON(memalloc);
#endif
}

void Provenance::inputComputeTimeStep() {
    incrementIterationsComputeTimeStep();
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Compute Time Step" << endl;
#endif

    string transformation = "computetimestep";
    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), taskID, numberIterationsComputeTimeStep);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(taskID);
    t.addPerformanceMetric(p);
    t.setSubID(numberIterationsComputeTimeStep);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("RUNNING");
    t.addDtDependency("computesolutionchange");

    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), taskID, numberIterationsComputeTimeStep);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-R.json", pgDirectory.c_str(), transformation.c_str(), taskID, numberIterationsComputeTimeStep);
    t.writeJSON(memalloc);
#endif
}

void Provenance::outputComputeTimeStep(int time_step, Real time, Real dt, bool timeStepAccepted) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Compute Time Step" << endl;
#endif

    string transformation = "computetimestep";
    Task t(taskID);
    t.setSubID(numberIterationsComputeSolutionChange);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("computesolutionchange");

    char memalloc[jsonArraySize];
    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%d;%.7f;%.7f;%s",
            simulationID, time_step, time, dt,
            timeStepAccepted ? "true" : "false"
            );
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    PerformanceMetric p;
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), taskID, numberIterationsComputeTimeStep);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), taskID, numberIterationsComputeTimeStep);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-F.json", pgDirectory.c_str(), transformation.c_str(), taskID, numberIterationsComputeTimeStep);
    t.writeJSON(memalloc);
#endif
}

void Provenance::inputMeshRefinement() {
    incrementIterationsMeshRefinements();

    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Mesh Refinement" << endl;
#endif

    string transformation = "meshrefinement";
    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), simulationID, numberIterationsMeshRefinements);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(taskID);
    t.setSubID(numberIterationsMeshRefinements);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("RUNNING");
    t.addDtDependency("solversimulationtransport");

    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, numberIterationsMeshRefinements);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID, numberIterationsMeshRefinements);
    t.writeJSON(memalloc);
#endif
}

void Provenance::outputMeshRefinement(bool first_step_refinement,
        int time_step, int before_n_active_elem, int after_n_active_elem) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Mesh Refinement" << endl;
#endif

    string transformation = "meshrefinement";
    PerformanceMetric p;
    char memalloc[jsonArraySize];
    sprintf(memalloc, "libMeshSedimentation::%s-%d-%d",
            transformation.c_str(), simulationID, numberIterationsMeshRefinements);
    p.SetDescription(memalloc);
    p.SetMethod("COMPUTATION");
    p.IdentifyStartTime();

    Task t(taskID);
    t.setSubID(numberIterationsMeshRefinements);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("solversimulationtransport");

    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

    sprintf(memalloc, "%d;%s;%d;%d;%d",
            simulationID, first_step_refinement ? "true" : "false", time_step,
            before_n_active_elem, after_n_active_elem);
    vector<string> e = {memalloc};
    t.addSet("o" + transformation, e);

    p.IdentifyEndTime();
    t.addPerformanceMetric(p);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, numberIterationsMeshRefinements);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID, numberIterationsMeshRefinements);
    t.writeJSON(memalloc);
#endif
}

void Provenance::inputMeshWriter() {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Input Mesh Writer" << endl;
#endif

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
    t.addDtDependency("solversimulationtransport");

    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-R.json", pgDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif
}

void Provenance::outputMeshWriter(int time_step, string xdmf) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Mesh Writer" << endl;
#endif

    string transformation = "meshwriter";
    Task t(taskID);
    t.setSubID(subTaskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("solversimulationtransport");

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

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif

    addMeshDependencyToList();
}

void Provenance::inputDataExtraction(int lineID) {
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
    t.addDtDependency("solversimulationtransport");

    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-R.json", jsonDirectory.c_str(), transformation, simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-R.json", pgDirectory.c_str(), transformation, simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif
}

void Provenance::outputDataExtraction(int lineID, int timeStep,
        string xdmf, int dimension) {
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

    Task t(taskID);
    t.setSubID(subTaskID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("solversimulationtransport");

    char memalloc[4096];
    sprintf(memalloc, "%d", taskID);
    t.addIdDependency(memalloc);

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

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-%d-F.json", jsonDirectory.c_str(), transformation, simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-%d-F.json", pgDirectory.c_str(), transformation, simulationID, subTaskID);
    t.writeJSON(memalloc);
#endif
}

void Provenance::meshAggregator(string xdmf, int n_processors) {
    if (processor_id != 0) return;
#ifdef VERBOSE
    cout << "Output Mesh Writer" << endl;
#endif

    PerformanceMetric ptemp;
    ptemp.IdentifyStartTime();

    string transformation = "meshaggregator";
    Task t(simulationID);
    t.setDataflow(dataflow);
    t.setTransformation(transformation);
    t.setWorkspace(directory);
    t.setStatus("FINISHED");
    t.addDtDependency("meshwriter");
    for (int dep : meshDependencies) {
        string depStr = std::to_string(dep);
        t.addIdDependency(depStr);
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

#ifdef DATABASE
    sprintf(memalloc, "%s%s-%d-F.json", jsonDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif
#ifdef BACKUP
    sprintf(memalloc, "%s%s-%d-F.json", pgDirectory.c_str(), transformation.c_str(), simulationID);
    t.writeJSON(memalloc);
#endif
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
}

void Provenance::writeMonitoringDataIntoFile(char* filename, int timeStep, int time, 
        Real initial_norm_delta, Real final_norm_delta, int linear_iteractions) {
    if (processor_id != 0) return;

    Real ratio = Real(final_norm_delta / initial_norm_delta);
    string flag;
    if (ratio == 1) {
        flag.assign("stagnant");
    } else if (ratio > 1) {
        flag.assign("bad_convergence");
    }else{
        flag.assign("ok");
    }

    FILE * monitoringFilePath = fopen(filename, "a");
    fprintf(monitoringFilePath, "%d;%d;%.7f;%.7f;%.7f;%d;%s\n",
            timeStep, time,
            initial_norm_delta, final_norm_delta,
            ratio,
            linear_iteractions,
            flag.c_str());
    fclose(monitoringFilePath);
};