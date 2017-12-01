#!/bin/bash
clear
killall monetdbd
killall mserver5
sleep 5
./stop-database.sh
# monetdb-start-all
# Restart MonetDB
# cd ../dfa
# ./restart-monetdb.sh
# cd ../sedimentation

echo "Starting database system..."
SIMULATION_DIR=/Users/vitor/Documents/repository/workflow-sedimentation/src-local
CPATH=$SIMULATION_DIR/sedimentation
DATAPATH=$CPATH/data
rm -rf data
unzip $CPATH/data.zip
$SIMULATION_DIR/bin/database_starter.sh database.conf $CPATH $DATAPATH