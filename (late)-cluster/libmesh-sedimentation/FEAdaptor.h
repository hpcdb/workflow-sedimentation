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

#include <iostream>
#include <string>

using namespace libMesh;
using namespace std;

namespace FEAdaptor
{
  void Initialize(int numScripts, string extractionScript, string visualizationScript);

  void Finalize();

  void CoProcess(int numScripts, string extractionScript, string visualizationScript, EquationSystems &eq, double time, unsigned int timeStep, bool lastTimeStep, bool using_amr);
  
}




#endif	/* FEADAPTOR_H */

