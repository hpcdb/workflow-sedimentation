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
    Provenance(int processorID);
    void inputInputMesh(int dim);
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

    void meshAggregator(string xdmf, int n_processors, vector<int> meshDependencies);

    void finishDataIngestor();

    void storeCatalystCost(double elapsedTime);
    void storeRDIComponentCost(double elapsedTime);
    void storeSolverCost(double elapsedTime);

    void createIndexDirectory();

    int getIndexerID() {
        return this->indexerID;
    }

    void setIndexerID(int indexerID) {
        this->indexerID = indexerID;
    }

    void incrementIndexerID() {
        this->indexerID++;
    }
    
    void addMeshDependencyToList() {
        meshDependencies.push_back(taskID);
    }
    
    vector<int> getMeshDependencies(){
        return meshDependencies;
    }
    
    int getNumberIterationsFluid() const {
        return numberIterationsFluid;
    }

    void setNumberIterationsFluid(int numberIterationsFluid) {
        this->numberIterationsFluid = numberIterationsFluid;
    }

    int getNumberIterationsMeshRefinements() const {
        return numberIterationsMeshRefinements;
    }

    void setNumberIterationsMeshRefinements(int numberIterationsMeshRefinements) {
        this->numberIterationsMeshRefinements = numberIterationsMeshRefinements;
    }

    int getNumberIterationsSediments() const {
        return numberIterationsSediments;
    }

    void setNumberIterationsSediments(int numberIterationsSediments) {
        this->numberIterationsSediments = numberIterationsSediments;
    }

    int getSubTaskID() const {
        return subTaskID;
    }

    void setSubTaskID(int subTaskID) {
        this->subTaskID = subTaskID;
    }

    int getTaskID() const {
        return taskID;
    }

    void setTaskID(int taskID) {
        this->taskID = taskID;
    }
    
    void incrementTaskID(){
        this->taskID++;
    }
    
    void incrementSubTaskID(){
        this->subTaskID++;
    }
    
    void incrementNumberIterationsFluid(){
        this->numberIterationsFluid++;
    }
    
    void incrementNumberIterationsMeshRefinements(){
        this->numberIterationsMeshRefinements++;
    }
    
    void incrementNumberIterationsSediments(){
        this->numberIterationsSediments++;
    }

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

    int indexerID = 0;
    vector<int> meshDependencies;
    
    int taskID = 0;
    int subTaskID = 0;
    int numberIterationsFluid = 0;
    int numberIterationsSediments = 0;
    int numberIterationsMeshRefinements = 0;    
};





