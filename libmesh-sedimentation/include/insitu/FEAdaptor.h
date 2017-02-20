/* 
 * File:   FEAdaptor.h
 * Author: camata
 *
 * Created on June 2, 2016, 11:46 AM
 */

#ifndef FEADAPTOR_H
#define	FEADAPTOR_H

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"

using namespace libMesh;

namespace FEAdaptor
{
    
  void mark_to_rebuild_grid();
    
  void Initialize(int numScripts, char* scripts[]);

  void Finalize();

  void CoProcess(EquationSystems &eq, double time, unsigned int timeStep, bool lastTimeStep);
  
}

#endif	/* FEADAPTOR_H */

