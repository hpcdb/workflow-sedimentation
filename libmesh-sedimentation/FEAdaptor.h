/* 
 * File:   FEAdaptor.h
 * Author: camata
 *
 * Created on June 2, 2016, 11:46 AM
 */

#ifndef FEADAPTOR_H
#define	FEADAPTOR_H

#define USE_CATALYST

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"

using namespace libMesh;

namespace FEAdaptor
{
  void Initialize(int numScripts, char* scripts[]);

  void Finalize();

  void CoProcess(int numScripts, char* scripts[],EquationSystems &eq, double time, unsigned int timeStep, bool lastTimeStep, bool using_amr);
  
}




#endif	/* FEADAPTOR_H */

