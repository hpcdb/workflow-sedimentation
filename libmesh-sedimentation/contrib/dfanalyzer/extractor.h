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
            char firstFilename[256];
            char finalFilename[256];
            sprintf(firstFilename, "init_ext_line_%d.csv", timeStep);
            sprintf(finalFilename, "ext_line_%d.csv", timeStep);

            if (libMesh::global_processor_id() == 0) {
                char commandLine[256];
                sprintf(commandLine, "python clean-csv.py %s %s;rm %s", firstFilename, finalFilename, firstFilename);
                cout << commandLine << endl;
                system(commandLine);
            }
        }

        void invoke3DRawDataExtractor(int processorID, int timeStep, int lineID) {
            char firstFilename[256];
            char finalFilename[256];
            sprintf(firstFilename, "init_ext_line_%d_%d.csv", lineID, timeStep);
            sprintf(finalFilename, "ext_line_%d_%d.csv", lineID, timeStep);

            if (libMesh::global_processor_id() == 0) {
                char commandLine[256];
                sprintf(commandLine, "python clean-csv.py %s %s;rm %s;chmod -R 774 %s", firstFilename, finalFilename, firstFilename, finalFilename);
                cout << commandLine << endl;
                system(commandLine);
            }
        }
    }


#ifdef __cplusplus
}
#endif

#endif /* EXTRACTOR_H */

