clear
# make
cd ../libmesh-sedimentation
make
cd ../example
# restart monetdb database
./restart-monetdb.sh
# delete files from the previous execution
./delete.sh
mkdir prov/pg/sedimentation
mkdir prov/di/sedimentation
cp dataflow.json prov/pg/sedimentation/
cp dataflow.json prov/di/sedimentation/
sleep 2
# solver execution
time ../libmesh-sedimentation/libmesh-sedimentation-opt