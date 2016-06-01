#!/bin/bash
clear
# make
cd ../libmesh-sedimentation
make
# solver execution
cd ../sedimentation
time ../libmesh-sedimentation/libmesh-sedimentation-opt
# Data Ingestor - Kill all java processes
killall java
# Calculate total elapsed time for provenance gathering
cd prov/log
directory=`pwd`
cd ../..
java -jar ../dfa/LogAnalysis.jar $directory