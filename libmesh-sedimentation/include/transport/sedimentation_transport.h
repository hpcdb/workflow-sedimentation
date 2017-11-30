/* 
 * File:   sedimentation_transport.h
 * Author: camata
 *
 * Created on September 16, 2015, 10:45 AM
 */

#ifndef SEDIMENTATION_TRANSPORT_H
#define SEDIMENTATION_TRANSPORT_H

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

class SedimentationTransport : public System::Assembly {
public:

    SedimentationTransport(EquationSystems &es_in) : es(es_in), Reynolds(1.0) {
    };

    void init();
    void assemble();
    void setup(GetPot &infile, bool restartControl = false);
    void solve(int t_step, Real dt, Real time, int r_step, bool& diverged);
    //void attach_time_stepping(timeStepControlBase *ts) {this->tsControl = ts;}
    void PrintMass(std::ofstream& fmass, unsigned int t_step);

    Real& nonlinear_tolerance() {
        return _nonlinear_tolerance;
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
        return transport_ssteady_count;
    };     

#ifdef PROVENANCE
    void attach_provenance(Provenance *provenance) {
        this->prov = provenance;
    }
#endif

    double mass_dep;
    double init_mass;
    double total_mass;
    double inlet_mass;
    double outlet_volume;


private:
    EquationSystems &es;
    timeStepControlBase *tsControl;

#ifdef PROVENANCE
    Provenance *prov;
#endif

    Real Reynolds;
    Real Grashof;
    int dim;
    void assemble3D();
    void assemble2D();
    void assembleSUPG2D();
    void assembleSUPG3D();
    void assembleRBVMS2D();
    void assembleRBVMS3D();    
    int erosion_bc_id;
    int noflux_bc_id;
    int deposition_id;
    bool apply_bottom_flow;
    int inlet_bc_id;

    Real _nonlinear_tolerance;
    Real _linear_tolerance;
    unsigned int _max_nonlinear_iteractions;
    unsigned int _linear_iteractions;
    unsigned int _nonlinear_iteractions;
    unsigned int _current_n_linear_iteractions;
    unsigned int transport_ssteady_count;
    Real _current_final_linear_residual;
    Real _solution_norm;
    Real _nonlinear_step_norm;
    Real _initial_linear_tolerance;
    Real _linear_tolerance_power;
};


#endif /* SEDIMENTATION_TRANSPORT_H */
