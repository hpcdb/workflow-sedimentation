#!/bin/bash
# delete files from the previous execution
rm restart.in
rm out*
rm temp*
rm -rf prov/log/*
rm -rf prov/pg/*
rm -rf prov/di/*

mkdir prov/log
mkdir prov/pg/sedimentation
mkdir prov/di/sedimentation
cp dataflow.json prov/pg/sedimentation/
cp dataflow.json prov/di/sedimentation/