#include <iostream>

using namespace std;

class ProvenanceObject {
  protected:
	int ID;

  public:
  	ProvenanceObject(int initID);
    int setID(int newID);
	int getID();
};