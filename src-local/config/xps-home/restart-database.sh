#!/bin/bash
clear
echo "Starting database system..."
SIMULATION_DIR=/home/vitor/Documents/dev/workflow-sedimentation/src-local
CPATH=$SIMULATION_DIR/sedimentation
DATAPATH=$CPATH/data
$SIMULATION_DIR/bin/database_starter.sh database.conf $CPATH $DATAPATH