#include "task.h"

void Task::writeJSON(string filename) {
//    StringBuffer s;
//    Writer<StringBuffer> writer(s);
//
//    writer.StartObject();
//    writer.Key("ID");
//    writer.Int(ID);
//    writer.EndObject();
//
//    cout << "##############" << endl;
//    cout << s.GetString() << endl;
//    cout << "##############" << endl;
//
//    ofstream file;
//    file.open(filename, ios_base::out);
//    file << s.GetString() << endl;
//    file.close();

    Document document;
    document.SetObject();
    Value v; // Null
    v.SetInt(10);
    document.AddMember("contacts","test",document.GetAllocator());

    printf("\nModified JSON with reformatting:\n");
    StringBuffer sb;
    PrettyWriter<StringBuffer> writer2(sb);
    document.Accept(writer2); // Accept() traverses the DOM and generates Handler events.
    puts(sb.GetString());
    
    ofstream file;
    file.open(filename, ios_base::out);
    file << sb.GetString() << endl;
    file.close();
}