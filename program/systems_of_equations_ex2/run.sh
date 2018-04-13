#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=systems_of_equations_ex2

export PARAVIEW_DIR=/media/vitor/data-linux/program/paraview-5.4
export DFANALYZER_DIR=/home/vitor/Documents/dev/workflow-sedimentation/program/dfa-lib
export FASTBIT_DIR=/media/vitor/data-linux/program/fastbit-2.0.3

export LD_LIBRARY_PATH=$DFANALYZER_DIR/lib:$LD_LIBRARY_PATH

run_example "$example_name" 
