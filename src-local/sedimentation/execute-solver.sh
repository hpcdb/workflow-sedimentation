#!/bin/bash
# solver execution

# mac
# export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/vitor/Documents/Program_Installations/ParaView-v5.2.0/build/CMakeFiles/__macos_install/lib/paraview-5.2:/usr/local/opt/gcc/lib/gcc/6/
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/vitor/Documents/Program_Installations/ParaView-v5.2.0/build/lib/paraview-5.2

# virtual box - mint
# export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/home/vitor/src/paraview/lib/paraview-5.1
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vitor/src/paraview/lib/paraview-5.1

# Mint
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/home/vitor/Documents/program/paraview/lib/paraview-5.1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vitor/Documents/program/paraview/lib/paraview-5.1

# stampede
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/03664/silva/programs/paraview/lib/paraview-5.0

#docker
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/programs/paraview/lib/paraview-5.0

rm output-solver.log
rm output-statistics.log

# time mpirun -np 2 ../libmesh-sedimentation/libmesh-sedimentation-opt analysis2D.py | tee -a "output-solver.log"
#time mpirun -np 2 ../libmesh-sedimentation/libmesh-sedimentation-opt analysis3D.py | tee -a "output-solver.log"

time ../libmesh-sedimentation/libmesh-sedimentation-opt -i sedimentation.in -m necker3d.msh -e extraction.py -v visualization.py -o output -d /home/vitor/Documents/dev/workflow-sedimentation/src-local/sedimentation/output/ | tee -a "output-solver.log"

#time ../libmesh-sedimentation/libmesh-sedimentation-opt
# Calculate total elapsed time for provenance gathering
cd prov
directory=`pwd`
cd ..
TMAX=0.015
DELTAT=0.005
WRITE_INTERVAL=2
IMAGES=`echo $TMAX $DELTAT $WRITE_INTERVAL | awk '{print int(($1/$2)/$3)}'`
java -jar ../dfa/ProvenanceAnalysis.jar $directory $IMAGES | tee -a "output-statistics.log"