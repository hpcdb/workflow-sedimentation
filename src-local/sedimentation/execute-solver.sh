#!/bin/bash
# solver execution
CASE_STUDY="meiburg2D"
if [ "$CASE_STUDY" == "necker3D" ]; then
	SOLVER="necker3D/lock_container"
	SOLVER_MESH=$SOLVER".msh"
	SOLVER_IN=$SOLVER"_cte.in"
	SOLVER_EXTRACTION=$SOLVER"_extraction.py"
	SOLVER_VISUALIZATION=$SOLVER"_visualization_surface.py,"$SOLVER"_visualization_volume.py,"$SOLVER"_visualization_wireframe.py"
elif [ "$CASE_STUDY" == "meiburg2D" ]; then
	SOLVER="meiburg2D/lock_meiburg2D"
	SOLVER_MESH=$SOLVER".msh"
	SOLVER_IN=$SOLVER"_pc11.in"
	SOLVER_EXTRACTION=$SOLVER"_extraction.py"
	SOLVER_VISUALIZATION=$SOLVER"_visualization.py"
fi
# ENVIRONMENT: macos,xps-nacad,xps-home,inspiron-laptop, macos
ENVIRONMENT="inspiron-laptop"
EXPERIMENT_DIR=""
if [ "$ENVIRONMENT" == "xps-nacad" ]; then
	EXPERIMENT_DIR="/home/vitor/Documents/dev/workflow-sedimentation/src-local/sedimentation"
elif [ "$ENVIRONMENT" == "xps-home" ]; then
	EXPERIMENT_DIR="/home/vitor/Documents/dev/workflow-sedimentation/src-local/sedimentation"	
elif [ "$ENVIRONMENT" == "inspiron-laptop" ]; then
	EXPERIMENT_DIR="/home/vitor/Documents/dev/workflow-sedimentation/src-local/sedimentation"	
elif [ "$ENVIRONMENT" == "macos" ]; then
	EXPERIMENT_DIR="/Users/vitor/Documents/repository/workflow-sedimentation/src-local/sedimentation"	
fi

export DYLD_LIBRARY_PATH=$PARAVIEW_DIR/lib/paraview-$PARAVIEW_VERSION
export LD_LIBRARY_PATH=$PARAVIEW_DIR/lib/paraview-$PARAVIEW_VERSION

echo "DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH
echo "LD_LIBRARY_PATH="$LD_LIBRARY_PATH

# send dataflow specification
./send-dataflow-spec.sh

rm output-solver.log
rm output-statistics.log

echo "time ../../libmesh-sedimentation/sediment-opt -i $SOLVER_IN -m $SOLVER_MESH -e $SOLVER_EXTRACTION -v $SOLVER_VISUALIZATION -o $CASE_STUDY -d $EXPERIMENT_DIR/$CASE_STUDY -dfa localhost -ksp_converged_use_min_initial_residual_norm | tee -a \"output-solver.log\""
time ../../libmesh-sedimentation/sediment-opt -i $SOLVER_IN -m $SOLVER_MESH -e $SOLVER_EXTRACTION -v $SOLVER_VISUALIZATION -o $CASE_STUDY -d $EXPERIMENT_DIR/output_$CASE_STUDY -dfa localhost -ksp_converged_use_min_initial_residual_norm | tee -a "output-solver.log"

killall python