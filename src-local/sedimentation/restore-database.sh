#!/bin/bash
clear
# monetdb-start-all
# Restart MonetDB
# cd ../dfa
# ./restart-monetdb.sh
# cd ../sedimentation

echo "Starting database system..."
SIMULATION_DIR=/home/vitor/Documents/dev/workflow-sedimentation/src-local
CPATH=$SIMULATION_DIR/sedimentation
DATAPATH=$CPATH/data
rm -rf data
unzip $SIMULATION_DIR/bin/monetdb-database.zip
$SIMULATION_DIR/bin/database_starter.sh database.conf $CPATH $DATAPATH