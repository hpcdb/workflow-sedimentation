#!/bin/bash
clear
# monetdb-start-all

# Restart MonetDB
cd ../dfa
./restart-monetdb.sh

# delete files from the previous execution
cd ../sedimentation
./delete.sh

# Start daemon process to Data Ingestor
cd ../sedimentation
java -jar ../dfa/DI-1.0.jar -daemon start

# commandline
#clear;monetdb-start-all;cd ../dfa;./restart-monetdb.sh;cd ../sedimentation;java -jar ../dfa/DI-1.0.jar -daemon start

# Back up database
# rm prov-db.dump
# mclient -p 50000 -d dataflow_analyzer --dump > prov-db.dump
# Restore database
# mclient -p 50000 -d dataflow_analyzer -lsql prov-db.dump
