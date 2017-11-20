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
SIMULATION_DIR=`pwd`
DATAPATH=$SIMULATION_DIR/data
rm -rf data
unzip $SIMULATION_DIR/data.zip
$SIMULATION_DIR/../bin/database_starter.sh database.conf $SIMULATION_DIR $DATAPATH