#!/bin/bash
# solver execution
time ../libmesh-sedimentation/libmesh-sedimentation-opt
# Calculate total elapsed time for provenance gathering
cd prov/log
directory=`pwd`
cd ../..
java -jar ../dfa/LogAnalysis.jar $directory