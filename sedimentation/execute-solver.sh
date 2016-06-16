#!/bin/bash
# solver execution
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/vitor/Documents/Program_Installations/ParaView-v5.0.1-source/build/CMakeFiles/__macos_install/lib/paraview-5.0
export LD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/vitor/Documents/Program_Installations/ParaView-v5.0.1-source/build/CMakeFiles/__macos_install/lib/paraview-5.0
time mpirun -np 2 ../libmesh-sedimentation/libmesh-sedimentation-opt teste.py
#time ../libmesh-sedimentation/libmesh-sedimentation-opt
# Calculate total elapsed time for provenance gathering
cd prov/log
directory=`pwd`
cd ../..
java -jar ../dfa/LogAnalysis.jar $directory