#!/bin/bash
echo "Setting up environment variables"
db_host=$1
SIMULATION_DIR=$2
DATA_DIR=$3
DFA_PROPERTIES=DfA.properties
DI_DIR=$SIMULATION_DIR/provenance
rm $DFA_PROPERTIES
echo "--------------------------------------------"
echo "Configuring DfA.properties file"
echo "pg_dir="$SIMULATION_DIR >> $DFA_PROPERTIES
echo "di_dir="$DI_DIR >> $DFA_PROPERTIES
echo "dbms=MONETDB" >> $DFA_PROPERTIES
echo "db_server=localhost" >> $DFA_PROPERTIES
echo "db_port=54321" >> $DFA_PROPERTIES
echo "db_name=dataflow_analyzer" >> $DFA_PROPERTIES
echo "db_user=monetdb" >> $DFA_PROPERTIES
echo "db_password=monetdb" >> $DFA_PROPERTIES
echo "dataflow_tag=clothing" >> $DFA_PROPERTIES
echo "directory="$SIMULATION_DIR >> $DFA_PROPERTIES
echo "--------------------------------------------"
cp $DATA_DIR/$DFA_PROPERTIES $SIMULATION_DIR/$DFA_PROPERTIES
cp $SIMULATION_DIR/database.conf $DATA_DIR/database.conf
echo "Restoring MonetDB database..."
cd $DATA_DIR
killall monetdbd
killall mserver5
killall java
sleep 3
rm dfa.log
rm -rf data
#unzip -q backup/data-loboc.zip
unzip /home/users/vitorss/local/all/monetdb/monetdb-database.zip
clear
echo "--------------------------------------------"
echo "Starting DfAnalyzer..."
DATAPATH=$DATA_DIR/data
$SIMULATION_DIR/../dfa/database_starter.sh database.conf $SIMULATION_DIR $DATAPATH
sleep 30
