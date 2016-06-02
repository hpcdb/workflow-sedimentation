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
#include <math.h>

#include "libmesh/libmesh.h"

#define PROV

using namespace std;
using namespace libMesh;

class Provenance
{
  	public:
  		Provenance();
  		void inputMeshGeneration(int simulationID, int dim, int ncellx, int ncelly, int ncellz, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, int ref_interval);
  		void outputMeshGeneration(int simulationID, double r_fraction,double c_fraction,double max_h_level,unsigned int hlevels);
  		void outputCreateEquationSystems(int simulationID, Real Reynolds,Real Gr,Real Sc,Real Us,Real Diffusivity,Real xlock,Real alfa,Real theta,Real ex,Real ey,Real ez,Real c_factor);
  		void outputGetMaximumIterations(int simulationID, Real dt, Real tmax, unsigned int n_time_steps, unsigned int n_nonlinear_steps, double nonlinear_tolerance, int max_linear_iters, int max_r_steps, unsigned int write_interval);
  		void outputSolverSimulationFluid(int simulationID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged);
  		void outputSolverSimulationSediments(int simulationID, int subTaskID, int time_step, Real time, int linear_step, int n_linear_step, unsigned int n_linear_iterations, Real linear_residual, Real norm_delta, Real norm_delta_u, bool converged);
      void outputMeshRefinement(int simulationID, int subTaskID, bool first_step_refinement, int before_n_active_elem, int after_n_active_elem);
      void finishDataIngestor();
    private:
      int processor_id;
};





