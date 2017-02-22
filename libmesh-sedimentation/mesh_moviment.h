/*
 * File:   mesh_moviment.h
 * Author: camata
 *
 * Created on September 16, 2015, 10:45 AM
 */

#ifndef MESH_MOVIMENT_H
#define	MESH_MOVIMENT_H

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
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
#include "libmesh/parallel_ghost_sync.h"
#include "define.h"

//#define MESH_MOVIMENT

// Bring in everything from the libMesh namespace
using namespace libMesh;


class MeshMoviment : public System::Assembly
{

  public:

  MeshMoviment (EquationSystems &es_in) : es (es_in), deposition_id(0), fixedwall_id(0)
  {};
  void init();
  void assemble ();
  void setup(GetPot &infile);
  void updateMesh();

   private:

     int deposition_id;
     int fixedwall_id;
     EquationSystems &es;
     
};


#endif	/* SEDIMENTATION_TRANSPORT_H */
