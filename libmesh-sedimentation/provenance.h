/* 
 * File:   prov.h
 * Author: vitor
 *
 * Created on May 18, 2016, 09:42 AM
 */

// C++ include files that we need
#include <iostream>
#include <sstream>

#define PROV

using namespace std;

class Provenance
{
  	public:
  		Provenance() {};
  		void step1(int dim, int ncellx, int ncelly, int ncellz, 
			double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, int ref_interval);
};