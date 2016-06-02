clear
# move data from /shared
cd ..
./copy.sh; 
# compile solver
cd libmesh-sedimentation
# rm .depend
# make clean
make
cd ../sedimentation
