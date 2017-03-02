/* 
 * File:   extractor.h
 * Author: vitor
 *
 * Created on February 23, 2017, 10:23 AM
 */

#ifndef EXTRACTOR_H
#define EXTRACTOR_H

#ifdef __cplusplus
extern "C" {
#endif

    namespace extractor {

        void invoke2DRawDataExtractor(int processorID, int timeStep) {
            if (processorID == 0) {
                char firstFilename[256];
                char finalFilename[256];
                sprintf(firstFilename, "init_ext_line_%d.csv", timeStep);
                sprintf(finalFilename, "ext_line_%d.csv", timeStep);
                
                // char commandLine[256];
                // sprintf(commandLine, "python clean-csv.py %s %s;rm %s", firstFilename, finalFilename, firstFilename);
                // cout << commandLine << endl;
                // system(commandLine);

                ifstream in(firstFilename);
                ofstream out(finalFilename);
                string wordToReplace("\"");
                string wordToReplace2(":");

                string line;
                size_t len = wordToReplace.length();
                size_t len2 = wordToReplace2.length();
                while (getline(in, line))
                {
                    if(line.find("nan") == std::string::npos){
                        size_t pos = line.find(wordToReplace);
                        if (pos != string::npos)
                            line.replace(pos, len, "");

                        size_t pos2 = line.find(wordToReplace2);
                        if (pos2 != string::npos)
                            line.replace(pos2, len2, "");

                        out << line << '\n';
                    }
                }
            }
        }

        void invoke3DRawDataExtractor(int processorID, int timeStep, int lineID) {
            if (processorID == 0) {
                char firstFilename[256];
                char finalFilename[256];
                sprintf(firstFilename, "init_ext_line_%d_%d.csv", lineID, timeStep);
                sprintf(finalFilename, "ext_line_%d_%d.csv", lineID, timeStep);

                // char commandLine[256];
                // sprintf(commandLine, "python clean-csv.py %s %s;rm %s;chmod -R 774 %s", firstFilename, finalFilename, firstFilename, finalFilename);
                // cout << commandLine << endl;
                // system(commandLine);

                ifstream in(firstFilename);
                ofstream out(finalFilename);
                string wordToReplace("\"");
                string wordToReplace2(":");

                string line;
                size_t len = wordToReplace.length();
                size_t len2 = wordToReplace2.length();
                while (getline(in, line))
                {
                    if(line.find("nan") == std::string::npos){
                        size_t pos = line.find(wordToReplace);
                        if (pos != string::npos)
                            line.replace(pos, len, "");

                        size_t pos2 = line.find(wordToReplace2);
                        if (pos2 != string::npos)
                            line.replace(pos2, len2, "");

                        out << line << '\n';
                    }
                }
            }
        }
    }


#ifdef __cplusplus
}
#endif

#endif /* EXTRACTOR_H */

