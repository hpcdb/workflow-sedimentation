#!/bin/bash
SIMULATION_DIR=`pwd`
DATAPATH=$SIMULATION_DIR/data
$SIMULATION_DIR/../bin/database_stopper.sh database.conf $SIMULATION_DIR $DATAPATH