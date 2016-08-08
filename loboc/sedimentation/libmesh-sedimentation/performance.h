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
#include <chrono>

#include "libmesh/libmesh.h"

#define PERFORMANCE

using namespace std;
using namespace libMesh;

class Performance
{
  	public:
      Performance(){};

      void start(){
        startTime = std::chrono::system_clock::now();
      }

      void end(){
        endTime = std::chrono::system_clock::now();
      }

      double elapsedTime(){
        return double(double(std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count())/1000000.00);
      }

    private:
      std::chrono::time_point<std::chrono::system_clock> startTime;
      std::chrono::time_point<std::chrono::system_clock> endTime;
};





