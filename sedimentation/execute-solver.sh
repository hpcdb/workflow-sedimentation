#!/bin/bash
# solver execution

# mac
# export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/vitor/Documents/Program_Installations/ParaView-v5.0.1-source/build/CMakeFiles/__macos_install/lib/paraview-5.0
# export LD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/vitor/Documents/Program_Installations/ParaView-v5.0.1-source/build/CMakeFiles/__macos_install/lib/paraview-5.0

# virtual box - mint
#export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/home/vitor/programs/paraview/lib/paraview-5.1
#export LD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/home/vitor/programs/paraview/lib/paraview-5.1

time mpirun -np 2 ../libmesh-sedimentation/libmesh-sedimentation-opt teste.py
#time ../libmesh-sedimentation/libmesh-sedimentation-opt
# Calculate total elapsed time for provenance gathering
cd prov/log
directory=`pwd`
cd ../..
java -jar ../dfa/LogAnalysis.jar $directory