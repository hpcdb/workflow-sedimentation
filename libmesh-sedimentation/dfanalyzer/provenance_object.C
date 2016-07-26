#include <iostream>

#include "provenance_object.h"

using namespace std;

ProvenanceObject::ProvenanceObject(int initID){ ID = initID; }

int ProvenanceObject::setID(int newID){ ID = newID; }

int ProvenanceObject::getID(){ return ID; }