#include "task.h"
#include <curl/curl.h>

void Task::sendRequest(string dfa_hostname) {
    Document document;
    document.SetObject();

    if (ID != 0) {
        char vs[10];
        int len = sprintf(vs, "%d", ID);
        Value v;
        v.SetString(vs, static_cast<SizeType>(len), document.GetAllocator());
        document.AddMember("id", v, document.GetAllocator());
    }

    if (subID != 0) {
        char vs[10];
        int len = sprintf(vs, "%d", subID);
        Value v;
        v.SetString(vs, static_cast<SizeType>(len), document.GetAllocator());
        document.AddMember("subid", v, document.GetAllocator());
    } else {
        Value v;
        document.AddMember("subid", v, document.GetAllocator());
    }

    if (!dataflow.empty()) {
        Value v;
        v.SetString(dataflow.c_str(), dataflow.size(), document.GetAllocator());
        document.AddMember("dataflow", v, document.GetAllocator());
    }

    if (!transformation.empty()) {
        Value v;
        v.SetString(transformation.c_str(), transformation.size(), document.GetAllocator());
        document.AddMember("transformation", v, document.GetAllocator());
    }

    if (!workspace.empty()) {
        Value v;
        v.SetString(workspace.c_str(), workspace.size(), document.GetAllocator());
        document.AddMember("workspace", v, document.GetAllocator());
    }

    if (!resource.empty()) {
        Value v;
        v.SetString(resource.c_str(), resource.size(), document.GetAllocator());
        document.AddMember("resource", v, document.GetAllocator());
    }

    if (!status.empty()) {
        Value v;
        v.SetString(status.c_str(), status.size(), document.GetAllocator());
        document.AddMember("status", v, document.GetAllocator());
    }

    Value aperfs(kArrayType);
    for (PerformanceMetric pm : this->performances) {
        Value p;
        p.SetObject();

        Value description;
        description.SetString(pm.GetDescription().c_str(), pm.GetDescription().size(), document.GetAllocator());
        p.AddMember("description", description, document.GetAllocator());

        Value method;
        method.SetString(pm.GetMethod().c_str(), pm.GetMethod().size(), document.GetAllocator());
        p.AddMember("method", method, document.GetAllocator());
        
        if(!pm.GetStartTime().empty()){
            Value time;
            time.SetString(pm.GetStartTime().c_str(), pm.GetStartTime().size(), document.GetAllocator());
            p.AddMember("startTime", time, document.GetAllocator());
        }
        
        if(!pm.GetEndTime().empty()){
            Value time;
            time.SetString(pm.GetEndTime().c_str(), pm.GetEndTime().size(), document.GetAllocator());
            p.AddMember("endTime", time, document.GetAllocator());
        }

        aperfs.PushBack(p, document.GetAllocator());
    }
    document.AddMember("performances", aperfs, document.GetAllocator());

    Value asets(kArrayType);
    map<string, vector < string>>::iterator it;
    for (it = this->sets.begin(); it != this->sets.end(); it++) {
        string currentSet = it->first;
        vector<string> currentElements = it->second;

        Value p;
        p.SetObject();

        Value set;
        set.SetString(currentSet.c_str(), currentSet.size(), document.GetAllocator());
        p.AddMember("tag", set, document.GetAllocator());

        Value elements(kArrayType);
        for (string e : currentElements) {
            Value ve;
            // ve.SetString(StringRef(e.c_str()));

            // char vs[1024];
            // int len = sprintf(vs, "%s", e);
            // Value v;
            ve.SetString(e.c_str(), e.size(), document.GetAllocator());

            elements.PushBack(ve, document.GetAllocator());
        }
        p.AddMember("elements", elements, document.GetAllocator());

        asets.PushBack(p, document.GetAllocator());
    }
    document.AddMember("sets", asets, document.GetAllocator());

    Value adeps;
    adeps.SetObject();
    if (this->dtDependencies.size() > 0) {
        Value adts(kArrayType);
        for (string dt : this->dtDependencies) {
            Value cv;
            cv.SetObject();

            Value cvalue;
            cvalue.SetString(dt.c_str(), dt.size(), document.GetAllocator());
            cv.AddMember("tag", cvalue, document.GetAllocator());

            adts.PushBack(cv, document.GetAllocator());
        }
        adeps.AddMember("tags", adts, document.GetAllocator());
    }

    if (this->idDependencies.size() > 0) {
        Value aids(kArrayType);
        for (string id : this->idDependencies) {
            Value cv;
            cv.SetObject();

            Value cvalue;
            cvalue.SetString(id.c_str(), id.size(), document.GetAllocator());
            cv.AddMember("id", cvalue, document.GetAllocator());

            aids.PushBack(cv, document.GetAllocator());
        }
        adeps.AddMember("ids", aids, document.GetAllocator());
    }
    document.AddMember("dependency", adeps, document.GetAllocator());
    
    Value afiles(kArrayType);
    for (File pfile : this->files) {
        Value p;
        p.SetObject();

        Value name;
        name.SetString(pfile.GetName().c_str(), pfile.GetName().size(), document.GetAllocator());
        p.AddMember("name", name, document.GetAllocator());

        Value path;
        path.SetString(pfile.GetPath().c_str(), pfile.GetPath().size(), document.GetAllocator());
        p.AddMember("path", path, document.GetAllocator());
        
        afiles.PushBack(p, document.GetAllocator());
    }
    document.AddMember("files", afiles, document.GetAllocator());

    StringBuffer sb;
    PrettyWriter<StringBuffer> writer2(sb);
    document.Accept(writer2); // Accept() traverses the DOM and generates Handler events.
    //    puts(sb.GetString());

//    ofstream file;
//    file.open(filename, ios_base::out);
//    file << sb.GetString() << endl;
//    file.close();

    CURL *hnd = curl_easy_init();
    curl_easy_setopt(hnd, CURLOPT_CUSTOMREQUEST, "POST");
    string hostname = "http://" + dfa_hostname + ":22000/pde/task/json";
    curl_easy_setopt(hnd, CURLOPT_URL, hostname.c_str());
    struct curl_slist *headers = NULL;
    headers = curl_slist_append(headers, "postman-token: 6afcae02-81cb-821f-379f-f66efb776d94");
    headers = curl_slist_append(headers, "cache-control: no-cache");
    headers = curl_slist_append(headers, "Content-Type: application/text");
    curl_easy_setopt(hnd, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(hnd, CURLOPT_POSTFIELDS, sb.GetString());
    curl_easy_setopt(hnd, CURLOPT_VERBOSE, 0L); //0 disable messages
    curl_easy_perform(hnd); //send request
    curl_easy_cleanup(hnd);
    curl_global_cleanup();
}