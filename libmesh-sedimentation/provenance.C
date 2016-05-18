/* 
 * File:   provenance.C
 * Author: vitor
 *
 * Created on May 18, 2016, 09:42 AM
 */

// C++ include files that we need
#include <iostream>

#include "provenance.h"

char space[] = "      ";

void Provenance::step1(int dim, int ncellx, int ncelly, int ncellz, 
			double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, int ref_interval)
{
	clock_t begin = clock();

	cout << "Step-1:Input" << endl;
  	cout << "PROV:" << endl << 
	    space << "dim(" + to_string(dim) + ")" << endl <<
	    space << "ncellx(" + to_string(ncellx) + ")" << endl <<
	    space << "ncelly(" + to_string(ncelly) + ")" << endl <<
	    space << "ncellz(" + to_string(ncellz) + ")" << endl <<
	    space << "xmin(" + to_string(xmin) + ")" << endl <<
	    space << "ymin(" + to_string(ymin) + ")" << endl <<
	    space << "zmin(" + to_string(zmin) + ")" << endl <<
	    space << "xmax(" + to_string(xmax) + ")" << endl <<
	    space << "ymax(" + to_string(ymax) + ")" << endl <<
	    space << "zmax(" + to_string(zmax) + ")" << endl <<
	    space << "ref_interval(" + to_string(ref_interval) + ")" << endl <<
	    space << "xmax(" + to_string(xmax) + ")" << endl;

	clock_t end = clock();
  	double elapsed_secs = (double(end - begin) / CLOCKS_PER_SEC);
  	cout << "elapsed-time: " << to_string(elapsed_secs) << " seconds." << endl;
}