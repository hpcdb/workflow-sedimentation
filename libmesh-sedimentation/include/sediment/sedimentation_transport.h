/* 
 * File:   sedimentation_transport.h
 * Author: camata
 *
 * Created on September 16, 2015, 10:45 AM
 */

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

// Bring in everything from the libMesh namespace
using namespace libMesh;
#ifndef SEDIMENTATION_TRANSPORT_H
#define	SEDIMENTATION_TRANSPORT_H

class SedimentationTransport : public System::Assembly
{
public:
    
  SedimentationTransport (EquationSystems &es_in) : es (es_in), Reynolds(1.0)
  {};

  void init();
  void assemble ();
  void setup(GetPot &infile);
  //void attach_time_stepping(timeStepControlBase *ts) {this->tsControl = ts;}
  void PrintMass(const char* fname);
  
  double   mass_dep;
  double init_mass ;
  double total_mass;
  
  
private:
  EquationSystems &es;
  timeStepControlBase *tsControl;
  
  Real Reynolds;
  Real Grashof;
  int dim;
  void assemble3D();
  void assemble2D();
  int erosion_bc_id;
  int noflux_bc_id;
  int deposition_id;
  
  //unsigned int n_transport_nonlinear_iterations_total;
  //unsigned int n_transport_linear_iterations_total;
  //unsigned int n_rejected_transport_linear_iterations_total;
  
  //unsigned int old_n_non_linear_iter_transport;
  
};


#endif	/* SEDIMENTATION_TRANSPORT_H */
