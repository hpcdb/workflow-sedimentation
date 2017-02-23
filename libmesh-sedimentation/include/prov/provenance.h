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
#include <ctime>
#include <vector>

#include "libmesh/libmesh.h"

#define PROV
// #define VERBOSE

using namespace std;
using namespace libMesh;

class Provenance {
public:
    Provenance();
    void inputInputMesh(int dim);
    void outputInputMesh(double r_fraction, double c_fraction, double max_h_level, unsigned int hlevels);

    void outputCreateEquationSystems(Real Reynolds, Real Gr, Real Sc, Real Us, Real Diffusivity, Real xlock, Real fopc, Real theta, Real ex, Real ey, Real ez, Real c_factor);

    void outputGetMaximumIterations(Real dt, Real tmax, unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance, int max_linear_iters, int max_r_steps, unsigned int write_interval, string xdmf);

    void inputInitDataExtraction(int lineID);
    void outputInitDataExtraction(int lineID, string xdmf, int dimension, int indexerID);

    void inputInitVisualization(int lineID);
    void outputInitVisualization(int lineID, int timeStep);

    void inputSolverSimulationFluid(int taskID, int subTaskID);
    void outputSolverSimulationFluid(int taskID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged);

    void inputSolverSimulationSediments(int taskID, int subTaskID);
    void outputSolverSimulationSediments(int taskID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged);

    void outputMeshRefinement(int taskID, int subTaskID, bool first_step_refinement, int time_step, int before_n_active_elem, int after_n_active_elem);

    void inputMeshWriter(int taskID, int subTaskID);
    void outputMeshWriter(int taskID, int subTaskID, int time_step, string xdmf);

    void inputDataExtraction(int taskID, int subTaskID, int lineID);
    void outputDataExtraction(int taskID, int subTaskID, int lineID, int time_step, string xdmf, int dimension, int indexerID);

    void inputVisualization(int lineID, int taskID);
    void outputVisualization(int lineID, int taskID, int time_step);

    void meshAggregator(string xdmf, int n_processors, vector<string> meshDependencies);

    void finishDataIngestor();

    void storeCatalystCost(int taskID, int subTaskID, double elapsedTime);
    void storeRDIComponentCost(int taskID, int subTaskID, double elapsedTime);
    void storeSolverCost(double elapsedTime);

    void createIndexDirectory();

private:
    int processor_id;
    int simulationID;
    const int jsonArraySize = 4096;
    const int arraySize = 256;
    string space = "      ";
    string directory = "";
    string pgCommandLine = "";
    string rdeCommandLine = "";
    string rdiCommandLine = "";
    string bin = "";
    string extraArguments = "";
    string dataflow = "sedimentation";
    string jsonDirectory = "";
    string pgDirectory = "";
    string rawDataAccess = "";
    string cartridge = "";
};





