/* 
 * File:   task.h
 * Author: vitor
 *
 * Created on March 31, 2018, 9:43 AM
 */

#include "dataset.h"
#include "dependency.h"
#include "task_status_enum.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>

using namespace std;

class Task {
protected:
    int id;
    int sub_id;
    string dataflow_tag;
    string transformation_tag;
    string workspace;
    map<string, Dataset> datasets; 
    string resource;
    task_status status;
    Dependency dependency;

public:
    Task(string dataflow_tag, string transformation_tag, int ID, int sub_id=0) {
        this->dataflow_tag = dataflow_tag;
        this->transformation_tag = transformation_tag;
        this->id = ID;
        this->sub_id = sub_id;
        this->dependency = Dependency();
    };
    
    int begin(); // return a code with the status of this service request
    int end();
    
    void add_dependent_transformation_tag(string transformation_tag);
    void add_dependent_transformation_tags(vector<string> transformation_tags);
    void add_dependent_transformation_id(int task_id);
    void add_dependent_transformation_ids(vector<int> transformation_ids);
    void add_element_to_dataset(string dataset_tag, string value);
    void add_element_to_dataset(string dataset_tag, vector<string> values);
    
    void set_workspace(string workspace);
    void set_resource(string resource);
    void set_status(task_status status);
    
    int get_id();
    int get_sub_id();
    string get_dataflow_tag();
    string get_transformation_tag();
    string get_workspace();
    string get_resource();
    Dataset& get_dataset(string tag);
    string get_status();    
    Dependency& get_dependency();    
    string get_specification();
};

