/* 
 * File:   prov.h
 * Author: vitor
 *
 * Created on May 18, 2016, 09:42 AM
 */

#ifndef PERFORMANCE_H
#define PERFORMANCE_H

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

class Performance
{
    public:
      Performance(){};

      void begin(){
        clock_gettime(CLOCK_REALTIME, &startTime);
      }

      void end(){
        clock_gettime(CLOCK_REALTIME, &endTime);
      }

      double getElapsedTime(){
          return (endTime.tv_sec - startTime.tv_sec + (endTime.tv_nsec - startTime.tv_nsec) / 1000000000.);
      }

    private:
      timespec startTime;
      timespec endTime;
};

#endif




