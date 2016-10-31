/* 
 * File:   prov.h
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
#include <vector>

#include "libmesh/libmesh.h"

#define PROV
//#define VERBOSE

using namespace std;
using namespace libMesh;

class Provenance {
public:
    Provenance();
    void inputMeshGeneration(int simulationID, int dim, int ncellx, int ncelly, int ncellz, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, int ref_interval);
    void outputMeshGeneration(int simulationID, double r_fraction, double c_fraction, double max_h_level, unsigned int hlevels);

    void outputCreateEquationSystems(int simulationID, Real Reynolds, Real Gr, Real Sc, Real Us, Real Diffusivity, Real xlock, Real fopc, Real theta, Real ex, Real ey, Real ez, Real c_factor);

    void outputGetMaximumIterations(int simulationID, Real dt, Real tmax, unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance, int max_linear_iters, int max_r_steps, unsigned int write_interval, string xdmf);

    void inputInitDataExtraction(int simulationID, string transformation);
    void outputInitDataExtraction(int simulationID, string transformation, string dataSet, int time_step, string xdmf, string rawDataFile, int dimension, string extractorName);

    void inputSolverSimulationFluid(int taskID, int simulationID, int subTaskID);
    void outputSolverSimulationFluid(int taskID, int simulationID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged);

    void inputSolverSimulationSediments(int taskID, int simulationID, int subTaskID);
    void outputSolverSimulationSediments(int taskID, int simulationID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged);

    void outputMeshRefinement(int taskID, int simulationID, int subTaskID, bool first_step_refinement, int time_step, int before_n_active_elem, int after_n_active_elem);

    void inputMeshWriter(int taskID, int simulationID, int subTaskID);
    void outputMeshWriter(int taskID, int simulationID, int subTaskID, int time_step, string xdmf);

    void inputDataExtraction(int taskID, int simulationID, int subTaskID, string transformation);
    void outputDataExtraction(int taskID, int simulationID, int subTaskID, string transformation, string dataSet, int time_step, string xdmf, string rawDataFile, int dimension, string extractorName);

    void meshAggregator(int simulationID, string xdmf, int n_processors, vector<string> meshDependencies);

    void finishDataIngestor();

    void storeDataExtractionCost(double elapsedTime);
    void storeRDEComponentCost(double elapsedTime);
    void storeSolverCost(double elapsedTime);
private:
    int processor_id;
    const int textArraySize = 64;
    const int jsonArraySize = 128;
    string space = "      ";
    string directory = "";
    string pgCommandLine = "";
    string rdeCommandLine = "";
    string rdiCommandLine = "";
    string bin="";
    string extraArguments = "";
    string dataflow = "sedimentation";
    string jsonDirectory = "";
    string pgDirectory = "";
    string rawDataAccess = "";
    string cartridge = "";
};





