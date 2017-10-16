# README - author: vitor
# libmesh-sedimentation folder (source): edit file Makefile (CATALYST?=1,PROVENANCE?=1)
# 
# open a terminal in the following directory: ./src-local/sedimentation
# create prov directory
rm -rf prov
mkdir prov
mkdir prov/pg
mkdir prov/di
# edit the file provenance.in (directory)
# 
# edit and run the script pg-dataflow.sh 
./pg-dataflow.sh
# edit and run the script execute-solver.sh (ENVIRONMENT,DYLD_LIBRARY_PATH,LD_LIBRARY_PATH)
./execute-solver.sh

