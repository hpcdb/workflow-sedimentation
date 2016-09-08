#!/bin/bash
# delete files from the previous execution
rm restart.in
rm out*
rm err.txt
rm db-out.txt
rm db-err.txt
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
rm -rf solver*
rm -rf mpirun.sh
rm -rf *.index
rm -rf *.data
rm -rf /mnt/scratch/user4/sc16/*

mkdir prov
mkdir prov/di
mkdir prov/pg
mkdir prov/log
mkdir prov/rde
mkdir prov/solver
mkdir prov/indexing
mkdir prov/pg/sedimentation
mkdir prov/di/sedimentation
cp dataflow.json prov/pg/sedimentation/
cp dataflow.json prov/di/sedimentation/
