/* 
 * File:   dataflow.h
 * Author: vitor
 *
 * Created on March 31, 2018, 11:03 AM
 */

#include "dfa_configuration.h"
#include "dfa_config.h"
#include "transformation.h"

#include <string>
#include <vector>

using namespace std;

class Dataflow{
protected:
    DfA_Config config;
    string tag;
    vector<Transformation> transformations;
    vector<Set> sets;
    
    string get_post_message();
    
public:
    Dataflow(string tag, string hostname=dfa_hostname, int port=dfa_port){
        this->config.hostname = hostname;
        this->config.port = port;
        this->tag = tag;
    }
    
    void save(); 
    
    void add_transformation(Transformation transformation, vector<Set> input_sets, vector<Set> output_sets);
    void add_transformation(Transformation transformation, Set input_sets, Set output_sets);
    void add_transformation(Transformation transformation, vector<Set> input_sets, Set output_sets);
    void add_transformation(Transformation transformation, Set input_sets, vector<Set> output_sets);
    void add_set(Set set);
};