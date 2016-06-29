#!/bin/bash
# delete files from the previous execution
rm restart.in
rm out*
rm temp*
rm *.csv
rm -rf prov/log/*
rm -rf prov/rde/*
rm -rf prov/solver/*
rm -rf prov/pg/*
rm -rf prov/di/*

mkdir prov/log
mkdir prov/rde
mkdir prov/solver
mkdir prov/pg/sedimentation
mkdir prov/di/sedimentation
cp dataflow.json prov/pg/sedimentation/
cp dataflow.json prov/di/sedimentation/