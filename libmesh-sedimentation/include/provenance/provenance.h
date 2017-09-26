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

#include "dfanalyzer/task.h"

// #define VERBOSE

using namespace std;
using namespace libMesh;

class Provenance {
public:

    Provenance(int processorID);

    void SetUp();
    void inputInputMesh();
    void outputInputMesh(int dim, string mesh_file, bool restarControl);

    void outputAMRConfig(double r_fraction, double c_fraction, double max_h_level, unsigned int hlevels, bool first_step_refinement,
            bool amrc_flow_transp, int ref_interval, int max_r_steps);

    void outputCreateEquationSystems(Real Reynolds, Real Gr, Real Sc, Real Us, Real Diffusivity, Real xlock, Real fopc, Real theta, Real ex, Real ey, Real ez, Real c_factor);

    void outputTSControlConfig(string ts_control_model_name, double dt_min, double dt_max, double tol_u, double tol_s,
            double kp, double ki, double kd, unsigned int nsa_max, unsigned int nsa_target_flow, unsigned int nsa_target_transport,
            unsigned int nsa_limit_flow, unsigned int nsa_limit_transport, double mult_factor_max, double mult_factor_min,
            double pc11_theta, double alpha, double k_exp, double s_min, double s_max, double reduct_factor, bool complete_flow_norm);

    void outputIOConfig(string dpath, string rname, unsigned int write_interval, unsigned int catalyst_interval, bool write_restart);

    void outputGetMaximumIterationsToFlow(Real dt, Real tmax, unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance, int max_linear_iters, string xdmf);

    void outputGetMaximumIterationsToTransport(Real dt, Real tmax, unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance, int max_linear_iters, string xdmf);

    void inputInitDataExtraction(int lineID);
    void outputInitDataExtraction(int lineID, string xdmf, int dimension);

    void inputInitVisualization(int lineID);
    void outputInitVisualization(int lineID, int timeStep);

    void inputSolverSimulationFlow();
    Task generateTaskToOutputSolverSimulationFlow();
    Task addElementToOutputSolverSimulationFlow(Task t, int time_step, Real dt, Real time,
        int linear_step, int n_linear_step, unsigned int n_linear_iterations,
        Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged);
    void finishTaskToOutputSolverSimulationFlow(Task t);

    void inputSolverSimulationTransport();
    void outputSolverSimulationTransport(int time_step, Real dt, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged);

    void inputComputeSolutionChange();
    void outputComputeSolutionChange(int time_step, Real time, Real dt,
            unsigned int n_flow_linear_iterations_total, unsigned int n_flow_nonlinear_iterations_total,
            unsigned int n_transport_linear_iterations_total, unsigned int n_transport_nonlinear_iterations_total,
            bool timeStepAccepted, double error);

    void inputComputeTimeStep();
    void outputComputeTimeStep(int time_step, Real time, Real dt, bool timeStepAccepted);

    void inputMeshRefinement();
    void outputMeshRefinement(bool first_step_refinement, int time_step, int before_n_active_elem, int after_n_active_elem);

    void inputMeshWriter();
    void outputMeshWriter(int time_step, string xdmf);

    void inputDataExtraction(int lineID);
    void outputDataExtraction(int lineID, int time_step, string xdmf, int dimension);

    void inputVisualization(int lineID);
    void outputVisualization(int lineID, int time_step);

    void meshAggregator(string xdmf, int n_processors);

    void writeMonitoringDataIntoFile(char* filename, int timeStep, Real time, 
        Real initial_norm_delta, Real final_norm_delta, int linear_iteractions);

    void finishDataIngestor();

    void createIndexDirectory();

    void incrementTaskID() {
        if (processor_id != 0) return;
        taskID++;
    }

    void incrementSubTaskID() {
        if (processor_id != 0) return;
        subTaskID++;
    }

    void incrementIterationsFlow() {
        if (processor_id != 0) return;
        numberIterationsFlow++;
    }

    void incrementIterationsTransport() {
        if (processor_id != 0) return;
        numberIterationsTransport++;
    }

    void incrementIterationsMeshRefinements() {
        if (processor_id != 0) return;
        numberIterationsMeshRefinements++;
    }

    void incrementIterationsComputeSolutionChange() {
        if (processor_id != 0) return;
        numberIterationsComputeSolutionChange++;
    }

    void incrementIterationsComputeTimeStep() {
        if (processor_id != 0) return;
        numberIterationsComputeTimeStep++;
    }

    void setTaskID(int taskID) {
        if (processor_id != 0) return;
        this->taskID = taskID;
    }

    void setSubTaskID(int id) {
        if (processor_id != 0) return;
        this->subTaskID = id;
    }

    void setNumberIterationsFlow(int value) {
        if (processor_id != 0) return;
        this->numberIterationsFlow = value;
    }

    void setNumberIterationsTransport(int value) {
        if (processor_id != 0) return;
        this->numberIterationsTransport = value;
    }

    void setNumberIterationsMeshRefinements(int value) {
        if (processor_id != 0) return;
        this->numberIterationsMeshRefinements = value;
    }

    void setNumberIterationsComputeSolutionChange(int value) {
        if (processor_id != 0) return;
        this->numberIterationsComputeSolutionChange = value;
    }

    void setNumberIterationsComputeTimeStep(int value) {
        if (processor_id != 0) return;
        this->numberIterationsComputeTimeStep = value;
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

    void resetIterationsFlow() {
        if (processor_id != 0) return;
        numberIterationsFlow = 0;
    }

    void resetIterationsTransport() {
        if (processor_id != 0) return;
        numberIterationsTransport = 0;
    }

    void resetIterationsMeshRefinement() {
        if (processor_id != 0) return;
        numberIterationsMeshRefinements = 0;
    }

    void resetIterationsComputeSolutionChange() {
        if (processor_id != 0) return;
        numberIterationsComputeSolutionChange = 0;
    }

    void resetIterationsComputeTimeStep() {
        if (processor_id != 0) return;
        numberIterationsComputeTimeStep = 0;
    }

    void incrementIndexerID() {
        if (processor_id != 0) return;
        indexerID++;
    }

    void setIndexerID(int id) {
        if (processor_id != 0) return;
        indexerID = id;
    }

    void resetIndexerID() {
        if (processor_id != 0) return;
        indexerID = 0;
    }

    void resetMeshDependencies() {
        if (processor_id != 0) return;
        meshDependencies.clear();
    }

    int getIndexerID() {
        return indexerID;
    }

    int getTaskID() {
        return taskID;
    }


private:
    int processor_id;
    int simulationID;

    int indexerID = 0;
    int taskID = 0;
    int subTaskID = 0;
    int numberIterationsFlow = 0;
    int numberIterationsTransport = 0;
    int numberIterationsMeshRefinements = 0;
    int numberIterationsComputeSolutionChange = 0;
    int numberIterationsComputeTimeStep = 0;
    vector<int> meshDependencies;

    //    to control the loop on the flow and the sediments
    int previousTaskID;

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





