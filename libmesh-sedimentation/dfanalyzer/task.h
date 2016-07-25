#include <iostream>

#include "provenance_object.h"

using namespace std;

class Task: public ProvenanceObject
{
  int test;

  public:
  	Task(int initID): ProvenanceObject(initID) {};
};