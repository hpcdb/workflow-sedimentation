/* 
 * File:   prov.h
 * Author: vitor
 *
 * Created on May 18, 2016, 09:42 AM
 */

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>

#include "libmesh/libmesh.h"

using namespace std;
using namespace libMesh;

class Performance {
public:

    Performance() {
    };

    void begin() {
        beginTime = clock();
    }

    void end() {
        endTime = clock();
    }

    double getElapsedTime() {
        return float( endTime - beginTime ) /  CLOCKS_PER_SEC;;
    }

private:
    clock_t beginTime;
    clock_t endTime;
};





