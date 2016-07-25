#include <iostream>

using namespace std;

class ProvenanceObject {
  protected:
	int ID;

  public:
  	ProvenanceObject(int initID){ ID = initID; };
    int setID(int newID){ ID = newID; }
	int getID(){ return ID; }
};