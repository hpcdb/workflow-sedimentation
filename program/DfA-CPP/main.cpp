/* 
 * File:   main.cpp
 * Author: vitor
 *
 * Created on March 31, 2018, 10:55 AM
 */

#include "dataflow.h"
#include "task.h"

#include <cstdlib>
#include <string>
#include <vector>

/*
 * 
 */
int main(int argc, char** argv) {

    cout << "running main..." << endl;

    //prospective provenance
    Dataflow dataflow = Dataflow("clothing");

    //set ideduplication
    vector<string> attribute_names = {"customerid", "country", "continent", "age", "gender", "children", "status"};
    vector<attribute_type> attribute_types = {NUMERIC, TEXT, TEXT, NUMERIC, TEXT, NUMERIC, TEXT};
    Set& ideduplication = dataflow.add_set("ideduplication", attribute_names, attribute_types);

    //set odeduplication
    Set& odeduplication = dataflow.add_set("odeduplication");
    odeduplication.add_extractor("ext_odeduplication", EXTRACTION, PROGRAM, attribute_names, attribute_types);

    //transformation deduplication
    Transformation& deduplication = dataflow.add_transformation("deduplication", ideduplication, odeduplication);
    
    //set oeurope
    Set& oeurope = dataflow.add_set("oeurope", attribute_names, attribute_types);
    
    //transformation europe
    Transformation& europe = dataflow.add_transformation("europe", odeduplication, oeurope);
    
    dataflow.save();
    
    //retrospective provenance
    

    return 0;
}

