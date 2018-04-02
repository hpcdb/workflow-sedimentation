
#include "dataflow.h"
#include <curl/curl.h>

void Dataflow::add_set(Set set) {
    this->sets.push_back(set);
}

void Dataflow::add_transformation(Transformation transformation, vector<Set> input_sets, vector<Set> output_sets) {
    transformation.set_input_sets(input_sets);
    transformation.set_output_sets(output_sets);
    this->transformations.push_back(transformation);
}

void Dataflow::add_transformation(Transformation transformation, Set input_set, Set output_set) {
    vector<Set> input_sets = {input_set};
    transformation.set_input_sets(input_sets);
    vector<Set> output_sets = {output_set};
    transformation.set_output_sets(output_sets);
    this->transformations.push_back(transformation);
}

void Dataflow::add_transformation(Transformation transformation, vector<Set> input_sets, Set output_set) {
    transformation.set_input_sets(input_sets);
    vector<Set> output_sets = {output_set};
    transformation.set_output_sets(output_sets);
    this->transformations.push_back(transformation);
}

void Dataflow::add_transformation(Transformation transformation, Set input_set, vector<Set> output_sets) {
    vector<Set> input_sets = {input_set};
    transformation.set_input_sets(input_sets);
    transformation.set_output_sets(output_sets);
    this->transformations.push_back(transformation);
}

string Dataflow::get_post_message() {
    string message = "dataflow(" + this->tag + ")";

    for (Set set : this->sets) {
        message += "\n" + set.get_specification();
    }

    for (Transformation transformation : this->transformations) {
        message += "\n" + transformation.get_specification();
    }
    
    return message;
}

void Dataflow::save() {
    CURL *hnd = curl_easy_init();
    curl_easy_setopt(hnd, CURLOPT_CUSTOMREQUEST, "POST");

    string hostname = "http://" +
            this->config.hostname + ":" + to_string(this->config.port) +
            "/pde/dataflow";
    curl_easy_setopt(hnd, CURLOPT_URL, hostname.c_str());

    struct curl_slist *headers = NULL;
    headers = curl_slist_append(headers, "postman-token: 6afcae02-81cb-821f-379f-f66efb776d94");
    headers = curl_slist_append(headers, "cache-control: no-cache");
    headers = curl_slist_append(headers, "Content-Type: application/text");

    string message = this->get_post_message();
    
    curl_easy_setopt(hnd, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(hnd, CURLOPT_POSTFIELDS, message.c_str());
    curl_easy_setopt(hnd, CURLOPT_VERBOSE, 0L); //0 disable messages

    curl_easy_perform(hnd); //send request
    curl_easy_cleanup(hnd);
    curl_global_cleanup();
}