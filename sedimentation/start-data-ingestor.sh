#!/bin/bash
clear
# Restart MonetDB
cd ../dfa
./restart-monetdb.sh
cd ../sedimentation
# delete files from the previous execution
./delete.sh
# Start daemon process to Data Ingestor
cd ../sedimentation
java -jar ../dfa/DI-1.0.jar -daemon start


# Back up database
# rm prov-db.dump
# mclient -p 50000 -d dataflow_analyzer --dump > prov-db.dump
# Restore database
# mclient -p 50000 -d dataflow_analyzer -lsql prov-db.dump