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
rm -rf log.*.out
rm -rf log.*.err
rm -rf *.conf
rm -rf nodes.txt
rm -rf *.dump
rm -rf data
rm -rf *.data
rm -rf *.index

mkdir prov
mkdir prov/di
mkdir prov/pg
mkdir prov/log
mkdir prov/rde
mkdir prov/solver
mkdir prov/pg/sedimentation
mkdir prov/di/sedimentation
cp dataflow.json prov/pg/sedimentation/
cp dataflow.json prov/di/sedimentation/