clear
# Restart MonetDB
cd ../dfa
./restart-monetdb.sh
cd ../example
# delete files from the previous execution
./delete.sh
mkdir prov/pg/sedimentation
mkdir prov/di/sedimentation
cp dataflow.json prov/pg/sedimentation/
cp dataflow.json prov/di/sedimentation/
# Start daemon process to Data Ingestor
cd ../example
java -jar ../dfa/DI-1.0.jar -daemon start