/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   viscosity_flow.h
 * Author: camata
 *
 * Created on 21 de Fevereiro de 2017, 21:31
 */

#ifndef VISCOSITY_FLOW_H
#define VISCOSITY_FLOW_H

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

#include "define.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

class ViscosityFlow : public System::Assembly
{
  public:
    
  ViscosityFlow (EquationSystems &es_in) : es (es_in)
  {};
  

  void assemble ();
  void init();
  //void setup(GetPot &infile);
  //void solve();

private:
  EquationSystems &es;
};

#endif /* VISCOSITY_FLOW_H */

