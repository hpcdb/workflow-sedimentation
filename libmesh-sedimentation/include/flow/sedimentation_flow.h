/* 
 * File:   sedimentation_flow.h
 * Author: camata
 *
 * Created on September 16, 2015, 10:12 AM
 */

#ifndef SEDIMENTATION_FLOW_H
#define SEDIMENTATION_FLOW_H

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/equation_systems.h"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/vector_value.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

// Define the PerfLog, a performance logging utility.
// It is useful for timing events in a code and giving
// you an idea where bottlenecks lie.
#include "libmesh/perf_log.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// To impose Dirichlet boundary conditions
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"

#include "timeStepControlBase.h"

#ifndef PROVENANCE_H
#define PROVENANCE_H
#include "../provenance/provenance.h"
#endif

// Bring in everything from the libMesh namespace
using namespace libMesh;

class SedimentationFlow : public System::Assembly {
public:

    SedimentationFlow(EquationSystems &es_in) : es(es_in), Reynolds(1.0) {
    };

    void assemble();
    void init();
    void setup(GetPot &infile, bool restartControl = false);
    void solve(int t_step, Real dt, Real time, int r_step, bool& diverged);
    void restart(GetPot &restart);
    //void attach_time_stepping(timeStepControlBase *ts) {this->tsControl = ts;}

    Real& non_linear_tolerance() {
        return _non_linear_tolerance;
    };

    Real& linear_tolerance() {
        return _linear_tolerance;
    };

    Real& initial_linear_tolerance() {
        return _initial_linear_tolerance;
    };

    Real current_final_linear_residual() {
        return _current_final_linear_residual;
    };

    Real& linear_tolerance_power() {
        return _linear_tolerance_power;
    };

    unsigned int& max_nonlinear_iteractions() {
        return _max_nonlinear_iteractions;
    };

    unsigned int linear_iteractions() {
        return _linear_iteractions;
    };

    unsigned int nonlinear_iteractions() {
        return _nonlinear_iteractions;
    };
    
    unsigned int get_ssteady_count() {
        return flow_ssteady_count;
    };    

    Number strainRateTensorNorm2D(const RealGradient& GradU, const RealGradient& GradV) {
        /* GradU = (dphi_x*U, dphi_y*U)
         * GradV = (dphi_x*V, dphi_y*V)
         * eij = 1/2 * ( dphi_j * v_i + dphi_i * v_j))
         */
        Number exx = GradU(0);
        Number exy = 0.5 * (GradU(1) + GradV(0)); // = eyx (symmetric tensor)
        Number eyy = GradV(1);

        // 2 * (e_ij * e_ij)
        Number t_norm = 2.0 * (exx * exx + 2.0 * exy * exy + eyy * eyy);

        return pow(t_norm, 0.5);
    }

    Number strainRateTensorNorm3D(const RealGradient& GradU, const RealGradient& GradV, const RealGradient& GradW) {
        /* GradU = (dphi_x*U, dphi_y*U, dphi_z*U)
         * GradV = (dphi_x*V, dphi_y*V, dphi_z*V)
         * GradV = (dphi_x*W, dphi_y*W, dphi_z*W)      
         * eij = 1/2 * ( dphi_j * v_i + dphi_i * v_j))
         */
        Number exx = GradU(0);
        Number exy = 0.5 * (GradU(1) + GradV(0)); // = eyx (symmetric tensor)
        Number exz = 0.5 * (GradU(2) + GradW(0)); // = ezx (symmetric tensor)
        Number eyy = GradV(1);
        Number eyz = 0.5 * (GradV(2) + GradW(1)); // = ezy (symmetric tensor)
        Number ezz = GradW(2);

        // 2 * (e_ij * e_ij)
        Number t_norm = 2.0 * (exx * exx + 2.0 * exy * exy + 2.0 * exz * exz + eyy * eyy + 2.0 * eyz * eyz + ezz * ezz);

        return pow(t_norm, 0.5);
    }

#ifdef PROVENANCE
    void attach_provenance(Provenance *provenance) {
        this->prov = provenance;
    }
#endif

private:
    void assemble2D();
    void assemble3D();
    void assembleSUPG2D();
    void assembleSUPG3D();
    void assembleRBVMS2D();
    void assembleRBVMS3D();    

#ifdef PROVENANCE
    Provenance *prov;
#endif

    Real Reynolds;
    Point normal;
    Real gravity;
    Real rho;
    Real viscosity;
    int outflow_id;
    Real _non_linear_tolerance;
    Real _linear_tolerance;
    unsigned int _max_nonlinear_iteractions;
    unsigned int _nonlinear_iteractions;
    unsigned int _current_n_linear_iteractions;
    unsigned int _linear_iteractions;
    unsigned int flow_ssteady_count;
    Real _current_final_linear_residual;
    Real _solution_norm;
    Real _nonlinear_step_norm;
    Real _initial_linear_tolerance;
    Real _linear_tolerance_power;
    int dim;
    EquationSystems &es;
};

#endif







