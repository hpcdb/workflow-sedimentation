
#include "dataflow.h"
#include <curl/curl.h>

string Dataflow::get_post_message(){
    string message = "dataflow(" + this->tag + ")";
    return message;
}

void Dataflow::save(){
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