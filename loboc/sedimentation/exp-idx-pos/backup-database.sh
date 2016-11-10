#!/bin/bash

echo "Starting server..."
monetdbd start /work/03664/silva/simulation/sedimentation/exp-idx-pos/data
monetdbd get all /work/03664/silva/simulation/sedimentation/exp-idx-pos/data

echo "Starting database..."
monetdb start dataflow_analyzer
monetdb status

echo "Backing up provenance database..."
rm prov-db.dump
mclient -p 54321 -d dataflow_analyzer --dump > prov-db.dump

#mclient -u monetdb -p54321 -ddataflow_analyzer