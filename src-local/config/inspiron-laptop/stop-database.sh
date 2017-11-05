#!/bin/bash
SIMULATION_DIR=/media/vitor/data-linux/dev/workflow-sedimentation/src-local
CPATH=$SIMULATION_DIR/sedimentation
DATAPATH=$CPATH/data
$SIMULATION_DIR/bin/database_stopper.sh database.conf $CPATH $DATAPATH