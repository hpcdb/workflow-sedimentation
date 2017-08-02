#!/bin/bash
# solver execution

# environments:xps-nacad,xps-home
environment="xps-home"
experiment_dir=""
# mpi
mpi=false
processors=2
# solver
SOLVER_IN=lock_necker3D_cte.in
SOLVER=lock_necker3D

if [ "$environment" == "xps-nacad" ]; then
	experiment_dir="/home/vitor/Documents/dev/workflow-sedimentation/src-local/sedimentation"
elif [ "$environment" == "xps-home" ]; then
	experiment_dir="/home/vitor/Documents/dev/workflow-sedimentation/src-local/sedimentation"
fi

export DYLD_LIBRARY_PATH=$PARAVIEW_DIR/lib/paraview-$PARAVIEW_VERSION
export LD_LIBRARY_PATH=$PARAVIEW_DIR/lib/paraview-$PARAVIEW_VERSION

echo $DYLD_LIBRARY_PATH
echo $LD_LIBRARY_PATH

rm output-solver.log
rm output-statistics.log

# without MPI
if [ "$mpi" ]; then
	time mpirun -np $processors ../../libmesh-sedimentation/sediment-opt -i $SOLVER_IN -m $SOLVER.msh -e $SOLVER_extraction.py -v $SOLVER_visualization.py -o output -d $experiment_dir/output -ksp_converged_use_min_initial_residual_norm | tee -a "output-solver.log"
else	
	time ../../libmesh-sedimentation/sediment-opt -i $SOLVER_IN -m $SOLVER.msh -e $SOLVER_extraction.py -v $SOLVER_visualization.py -o output -d $experiment_dir/output -ksp_converged_use_min_initial_residual_norm | tee -a "output-solver.log"
fi

# mac
#export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/vitor/Documents/program/paraview-5.4.0/CMakeFiles/__macos_install/lib/paraview-5.4:/usr/local/opt/gcc/lib/gcc/6/
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/vitor/Documents/program/paraview-5.4.0/CMakeFiles/__macos_install/lib/paraview-5.4

# virtual box - mint
# export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/home/vitor/src/paraview/lib/paraview-5.1
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vitor/src/paraview/lib/paraview-5.1

# Mint
# export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/home/vitor/Documents/program/paraview/lib/paraview-5.1
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vitor/Documents/program/paraview/lib/paraview-5.1

# Ubuntu - VMWare
# export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/home/vitor/programs/paraview-5.3.0/lib/paraview-5.3
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vitor/programs/paraview-5.3.0/lib/paraview-5.3

# XPS - Ubuntu - NACAD
# export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/home/vitor/program/paraview-v5.3.0/lib/paraview-5.3
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vitor/program/paraview-v5.3.0/lib/paraview-5.3

# Mint - NACAD
# export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/home/vitor/program/paraview-5.3.0/lib/paraview-5.3
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vitor/program/paraview-5.3.0/lib/paraview-5.3

# Ubuntu - Inspiron
#export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/home/vitor/Documents/program/paraview-v5.4.0/lib/paraview-5.4
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/vitor/Documents/program/paraview-v5.4.0/lib/paraview-5.4

# stampede
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/03664/silva/programs/paraview/lib/paraview-5.0

#docker
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/programs/paraview/lib/paraview-5.0

# rm output-solver.log
# rm output-statistics.log

# time mpirun -np 2 ../libmesh-sedimentation/libmesh-sedimentation-opt analysis2D.py | tee -a "output-solver.log"
#time mpirun -np 2 ../libmesh-sedimentation/libmesh-sedimentation-opt analysis3D.py | tee -a "output-solver.log"

# without catalyst
# time mpirun -np 2 ../../libmesh-sedimentation/sediment-opt -i sedimentation.in -m necker3d.msh -o output -d /home/vitor/Documents/dev/workflow-sedimentation/src-local/sedimentation/output/ | tee -a "output-solver.log"
# with catalyst
# time mpirun -np 2 ../../libmesh-sedimentation/sediment-opt -i lock_necker3D_pc11.in -m lock_necker3D.msh -e extraction.py -v visualization.py -o output -d /home/vitor/Documents/dev/workflow-sedimentation/src-local/sedimentation/output | tee -a "output-solver.log"

# MacOS
# time mpirun -np 2 ../../libmesh-sedimentation/sediment-opt -i lock_necker3D_pc11.in -m lock_necker3D.msh -e extraction.py -v visualization.py -o output -d /Users/vitor/Documents/repository/workflow-sedimentation/src-local/sedimentation/output -ksp_converged_use_min_initial_residual_norm | tee -a "output-solver.log"

# XPS - Ubuntu - NACAD
# time mpirun -np 2 ../../libmesh-sedimentation/sediment-opt -i lock_necker3D_pc11.in -m lock_necker3D.msh -e extraction.py -v visualization.py -o output -d /home/vitor/Documents/dev/workflow-sedimentation/src-local/sedimentation/output -ksp_converged_use_min_initial_residual_norm | tee -a "output-solver.log"
# time ../../libmesh-sedimentation/sediment-opt -i lock_necker3D_pc11.in -m lock_necker3D.msh -e extraction.py -v visualization.py -o output -d /home/vitor/Documents/dev/workflow-sedimentation/src-local/sedimentation/output -ksp_converged_use_min_initial_residual_norm | tee -a "output-solver.log"

# Ubuntu - Inspiron
#time mpirun -np 2 ../../libmesh-sedimentation/sediment-opt -i lock_necker3D_pc11.in -m lock_necker3D.msh -e extraction.py -v visualization.py -o output -d /media/vitor/data-linux/dev/workflow-sedimentation/src-local/sedimentation/output -ksp_converged_use_min_initial_residual_norm | tee -a "output-solver.log"
