#include "dependency.h"

void Dependency::add_transformation_tag(string transformation_tag){
    this->transformation_tags.push_back(transformation_tag);
}

void Dependency::add_transformation_ids(vector<int> transformation_ids){
    this->transformation_ids.push_back(transformation_ids);
}

vector<string>& Dependency::get_transformation_tags(){
    return this->transformation_tags;
}

vector<vector<int>>& Dependency::get_transformation_ids(){
    return this->transformation_ids;
}