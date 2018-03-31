/* 
 * File:   task.h
 * Author: vitor
 *
 * Created on March 31, 2018, 9:43 AM
 */

#include "dataset.h"
#include "dependency.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class Task {
protected:
    int ID;
    int subID;
    string dataflow_tag;
    string transformation_tag;
    string workspace;
    vector<Dataset> sets; 
    string resource;
    enum task_status {READY, RUNNING, FINISHED};
    task_status status;
    Dependency dependency;

public:
    Task(string dataflow_tag, string transformation_tag, int ID, int subID=0) {
        this->dataflow_tag = dataflow_tag;
        this->transformation_tag = transformation_tag;
        this->ID = ID;
    };
    
    int begin(); // return a code with the status of this service request
    int end();
    
    void set_workspace(string workspace);
    void set_resource(string resource);
    void add_element_to_dataset(string dataset_tag, int number_of_values, ...);
    void add_element(string dataset_tag, Element element);
    void set_status(task_status status);
    void add_dependent_transformation_tags(int number_of_transformations, ...);
    void add_dependent_task_ids(...);
    
    int get_id();
    string get_dataflow_tag();
    string get_transformation_tag();
    string get_workspace();
    string get_resource();
    Dataset get_dataset(string tag);
    task_status get_status();
    Dependency get_dependency();
};

