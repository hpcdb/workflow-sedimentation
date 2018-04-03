/* 
 * File:   main.cpp
 * Author: vitor
 *
 * Created on March 31, 2018, 10:55 AM
 */

#include "dfa.h"

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
    Set& odeduplication = dataflow.add_set("odeduplication", attribute_names, attribute_types);
    
//    with raw data extraction
//    Set& odeduplication = dataflow.add_set("odeduplication");
//    odeduplication.add_extractor("ext_odeduplication", EXTRACTION, PROGRAM, attribute_names, attribute_types);

    //transformation deduplication
    Transformation& deduplication = dataflow.add_transformation("deduplication", ideduplication, odeduplication);
    
    //set oeurope
    Set& oeurope = dataflow.add_set("oeurope", attribute_names, attribute_types);
    
    //transformation europe
    Transformation& europe = dataflow.add_transformation("europe", odeduplication, oeurope);
    
    dataflow.save();
    
    //retrospective provenance
    Task task_deduplication = Task(dataflow.get_tag(), deduplication.get_tag(), 1);
    task_deduplication.set_workspace("/path/");
    task_deduplication.set_resource("local");
    
    vector<string> ideduplication_values = {"1","Botswana","Africa","46","female","1","single"};
    Dataset& ds_ideduplication = task_deduplication.add_dataset(ideduplication.get_tag());
    ds_ideduplication.add_element_with_values(ideduplication_values);

    task_deduplication.begin();
    
    ds_ideduplication.clear_elements();
    Dataset& ds_odeduplication = task_deduplication.add_dataset(odeduplication.get_tag());
    ds_odeduplication.add_element_with_values(ideduplication_values);
    task_deduplication.end();
    
    return 0;
}

