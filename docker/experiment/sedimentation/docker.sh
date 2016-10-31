#!/bin/bash
rm -rf /experiment/*
cp -rf /shared/experiment/* /experiment
./compile.sh
./delete.sh
./execute-solver.sh