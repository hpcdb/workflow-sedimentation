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
#include "libmesh/perf_log.h"

// #define VERBOSE

using namespace std;
using namespace libMesh;

class Provenance {
public:

    Provenance(int processorID);

    void inputInputMesh(int dim, string mesh_file);
    void outputInputMesh(double r_fraction, double c_fraction, double max_h_level, unsigned int hlevels);

    void outputCreateEquationSystems(Real Reynolds, Real Gr, Real Sc, Real Us, Real Diffusivity, Real xlock, Real fopc, Real theta, Real ex, Real ey, Real ez, Real c_factor);

    void outputGetMaximumIterations(Real dt, Real tmax, unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance, int max_linear_iters, int max_r_steps, unsigned int write_interval, string xdmf);

    void inputInitDataExtraction(int lineID);
    void outputInitDataExtraction(int lineID, string xdmf, int dimension);

    void inputInitVisualization(int lineID);
    void outputInitVisualization(int lineID, int timeStep);

    void inputSolverSimulationFluid();
    void outputSolverSimulationFluid(int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged);

    void inputSolverSimulationSediments();
    void outputSolverSimulationSediments(int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged);

    void outputMeshRefinement(bool first_step_refinement, int time_step, int before_n_active_elem, int after_n_active_elem);

    void inputMeshWriter();
    void outputMeshWriter(int time_step, string xdmf);

    void inputDataExtraction(int lineID);
    void outputDataExtraction(int lineID, int time_step, string xdmf, int dimension);

    void inputVisualization(int lineID);
    void outputVisualization(int lineID, int time_step);

    void meshAggregator(string xdmf, int n_processors);

    void finishDataIngestor();

    void storeCatalystCost(double elapsedTime);
    void storeRDIComponentCost(double elapsedTime);
    void storeSolverCost(double elapsedTime);

    void createIndexDirectory();

    void incrementTaskID() {
        if (processor_id != 0) return;
        taskID++;
    }

    void incrementSubTaskID() {
        if (processor_id != 0) return;
        subTaskID++;
    }

    void incrementIterationsFluid() {
        if (processor_id != 0) return;
        numberIterationsFluid++;
    }

    void incrementIterationsSediments() {
        if (processor_id != 0) return;
        numberIterationsSediments++;
    }

    void incrementIterationsMeshRefinements() {
        if (processor_id != 0) return;
        numberIterationsMeshRefinements++;
    }

    void setTaskID(int taskID) {
        if (processor_id != 0) return;
        this->taskID = taskID;
    }

    void setSubTaskID(int id) {
        if (processor_id != 0) return;
        this->subTaskID = id;
    }

    void setNumberIterationsFluid(int value) {
        if (processor_id != 0) return;
        this->numberIterationsFluid = value;
    }

    void setNumberIterationsSediments(int value) {
        if (processor_id != 0) return;
        this->numberIterationsSediments = value;
    }

    void setNumberIterationsMeshRefinements(int value) {
        if (processor_id != 0) return;
        this->numberIterationsMeshRefinements = value;
    }

    void addMeshDependencyToList() {
        if (processor_id != 0) return;
        meshDependencies.push_back(taskID);
    }

    void resetTaskID() {
        if (processor_id != 0) return;
        taskID = 0;
    }

    void resetSubTaskID() {
        if (processor_id != 0) return;
        subTaskID = 0;
    }

    void resetIterationsFluid() {
        if (processor_id != 0) return;
        numberIterationsFluid = 0;
    }

    void resetIterationsSediments() {
        if (processor_id != 0) return;
        numberIterationsSediments = 0;
    }

    void resetIterationsMeshRefinement() {
        if (processor_id != 0) return;
        numberIterationsMeshRefinements = 0;
    }
    
    void incrementIndexerID(){
        if (processor_id != 0) return;
        indexerID++;
    }
    
    void setIndexerID(int id){
        if (processor_id != 0) return;
        indexerID = id;
    }
    
    void resetIndexerID(){
        if (processor_id != 0) return;
        indexerID = 0;
    }
    
    void resetMeshDependencies(){
        if (processor_id != 0) return;
        meshDependencies.clear();
    }
    
    int getIndexerID(){
        return indexerID;
    }


private:
    int processor_id;
    int simulationID;

    int indexerID = 0;
    int taskID = 0;
    int subTaskID = 0;
    int numberIterationsFluid = 0;
    int numberIterationsSediments = 0;
    int numberIterationsMeshRefinements = 0;
    vector<int> meshDependencies;

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





