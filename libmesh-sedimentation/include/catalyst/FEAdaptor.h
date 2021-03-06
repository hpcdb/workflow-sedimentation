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
  void mark_to_rebuild_grid();

  void Initialize(int numScripts, string extractionScript, vector<string> visualizationScripts);

  void Finalize();

  void CoProcess(EquationSystems &eq, double time, unsigned int timeStep, unsigned int analysisInterval, bool lastTimeStep, bool using_amr);
  
}




#endif	/* FEADAPTOR_H */

