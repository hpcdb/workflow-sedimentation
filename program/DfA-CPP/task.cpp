#include "task.h"

void Task::set_workspace(string workspace) {
    this->workspace = workspace;
}

void Task::set_resource(string resource) {
    this->resource = resource;
}

void Task::set_status(task_status status) {
    this->status = status;
}

int Task::get_id() {
    return this->id;
}

int Task::get_sub_id() {
    return this->sub_id;
}

string Task::get_dataflow_tag() {
    return this->dataflow_tag;
}

string Task::get_transformation_tag() {
    return this->transformation_tag;
}

string Task::get_workspace() {
    return this->workspace;
}

string Task::get_resource() {
    return this->resource;
}

Dataset& Task::get_dataset(string tag) {
    return this->datasets.find(tag)->second;
}

string Task::get_status() {
    string status_str("unknown");
    switch (this->status) {
        case READY:
        {
            status_str = "READY";
        }
            break;
        case RUNNING:
        {
            status_str = "RUNNING";
        }
            break;
        case FINISHED:
        {
            status_str = "FINISHED";
        }
            break;
    }
    return status_str;
}

Dependency& Task::get_dependency(){
    return this->dependency;
}

