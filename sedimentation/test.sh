clear
# Restart MonetDB
cd ../dfa
./restart-monetdb.sh
cd ../sedimentation
# delete files from the previous execution
./delete.sh
# Start daemon process to Data Ingestor
cd ../sedimentation
