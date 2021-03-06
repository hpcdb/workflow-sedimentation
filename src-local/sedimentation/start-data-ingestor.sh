#!/bin/bash
clear
echo "Setting up environment variables"
SIMULATION_DIR=`pwd`
DFA_PROPERTIES=$SIMULATION_DIR/DfA.properties
DI_DIR=$SIMULATION_DIR
echo "--------------------------------------------"
echo "Removing data from previous executions"
rm $DFA_PROPERTIES
# organizing provenance directories
rm -rf provenance
./delete.sh
# cleaning up MonetDB directory
killall mserver5
killall monetdbd
sleep 3
rm -rf data
rm prov-db.dump
# configuring computational environment for application run
rm database.conf
rm nodes.txt
# killing previous processes
# killall java
echo "--------------------------------------------"
echo "Configuring DfA.properties file"
echo "di_dir="$DI_DIR >> $DFA_PROPERTIES
echo "dbms=MONETDB" >> $DFA_PROPERTIES
echo "db_server=localhost" >> $DFA_PROPERTIES
echo "db_port=50000" >> $DFA_PROPERTIES
echo "db_name=dataflow_analyzer" >> $DFA_PROPERTIES
echo "db_user=monetdb" >> $DFA_PROPERTIES
echo "db_password=monetdb" >> $DFA_PROPERTIES
echo "--------------------------------------------"
echo "Restoring MonetDB database..."
echo "localhost" >> database.conf
unzip -q data.zip
clear
echo "--------------------------------------------"
echo "Starting database system..."
DATAPATH=$SIMULATION_DIR/data
$SIMULATION_DIR/../dfa/database_starter.sh database.conf $SIMULATION_DIR $DATAPATH
echo "--------------------------------------------"
echo "Starting performance monitor..."
$SIMULATION_DIR/../monitor/start_monitor.sh localhost 
# > out.txt 2> err.txt
echo "--------------------------------------------"
echo "Starting DfA RESTful API"
java -jar $SIMULATION_DIR/../dfa/REST-DfA-1.0.jar
