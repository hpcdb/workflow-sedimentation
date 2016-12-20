cd /experiment
rm -rf *
cp -rf /shared/experiment/* /experiment

cd libmesh-sedimentation
make
cp -rf * /experiment/libmesh-sedimentation/

cd ../sedimentation
# ./pg-dataflow.sh
# cp -rf dataflow.json ~/shared/experiment/sedimentation
# ./delete.sh
# ./execute-solver.sh
