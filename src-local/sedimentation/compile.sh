#!/bin/bash
clear
# compile solver
cd ../../libmesh-sedimentation
# rm .depend
# make clean
make PROVENANCE=1 CATALYST=1
cd ../experiment/sedimentation