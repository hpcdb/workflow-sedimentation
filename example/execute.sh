clear
# make
cd ../libmesh-sedimentation
make
cd ../example
# delete files from the previous execution
./delete.sh
cp dataflow.json prov/pg/sedimentation/
# solver execution
time ../libmesh-sedimentation/libmesh-sedimentation-opt