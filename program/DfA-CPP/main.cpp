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
    cout << dataflow.get_tag() << endl;

    //set ideduplication
    vector<string> attribute_names = {"customerid", "country", "continent", "age", "gender", "children", "status"};
    vector<attribute_type> attribute_types = {NUMERIC, TEXT, TEXT, NUMERIC, TEXT, NUMERIC, TEXT};
    Set& ideduplication = dataflow.add_set("ideduplication", attribute_names, attribute_types);

    //set odeduplication
    Set& odeduplication = dataflow.add_set("odeduplication", attribute_names, attribute_types);

    //transformation deduplication
    Transformation& deduplication = dataflow.add_transformation("deduplication", ideduplication, odeduplication);
    
    //set oeurope
    Set& oeurope = dataflow.add_set("oeurope");    
    
    //with raw data extraction
    oeurope.add_extractor("ext_oeurope", EXTRACTION, PROGRAM, attribute_names, attribute_types);
    
    //transformation europe
    Transformation& europe = dataflow.add_transformation("europe", odeduplication, oeurope);
    
    dataflow.save();
    
    //retrospective provenance
    int task_id = 1;
    Task task_deduplication = Task(dataflow.get_tag(), deduplication.get_tag(), task_id);
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
    
    task_id++;
    Task task_europe = Task(dataflow.get_tag(), europe.get_tag(), task_id);
    task_europe.set_workspace("/path/");
    task_europe.set_resource("local");
    
    vector<string> transformation_tags = {deduplication.get_tag()};
    task_europe.add_dependent_transformation_tags(transformation_tags);
    task_europe.add_dependent_transformation_id(task_deduplication.get_id());

    task_europe.begin();
    
    Dataset& ds_oeurope = task_europe.add_dataset(oeurope.get_tag());    
    vector<string> oeurope_values = {"/home/vitor/Documents/dev/workflow-sedimentation/program/DfA-CPP/ext.data"};
    ds_oeurope.add_element_with_values(oeurope_values);
    task_europe.end();
    
    return 0;
}

