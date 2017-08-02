#!/bin/bash
# delete files from the previous execution
rm restart.in
rm -rf out*
rm temp*
rm *.csv
rm -rf prov/pg/*
rm -rf prov/di/*
rm -rf log.*.out
rm -rf log.*.err
# rm -rf *.conf
rm -rf nodes.txt
# rm -rf *.dump
rm -rf data
rm -rf *.data
rm -rf *.index
rm -rf index
rm -rf image*
rm -rf video*
rm abort.run

rm -rf monitoring*.log
echo "time_step;time;initial_norm_delta;final_norm_delta;ratio_norm_delta;linear_iterations;flag" > monitoring-flow.log
echo "time_step;time;initial_norm_delta;final_norm_delta;ratio_norm_delta;linear_iterations;flag" > monitoring-transport.log

# environments: xps-nacad,xps-home,inspiron-laptop
ENVIRONMENT="xps-home"
cp -rf ../config/$ENVIRONMENT/* .

mkdir prov
mkdir prov/di
mkdir prov/pg
mkdir prov/pg/sedimentation
mkdir prov/di/sedimentation
cp dataflow.json prov/pg/sedimentation/
cp dataflow.json prov/di/sedimentation/
mkdir output
mkdir index