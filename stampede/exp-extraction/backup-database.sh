#!/bin/bash

echo "Starting server..."
monetdbd start /work/03664/silva/simulation/sedimentation/exp-extraction/data
monetdbd get all /work/03664/silva/simulation/sedimentation/exp-extraction/data

echo "Starting database..."
monetdb start dataflow_analyzer
monetdb status

echo "Backing up provenance database..."
rm prov-db.dump
mclient -p 50000 -d dataflow_analyzer --dump > prov-db.dump

#mclient -u monetdb -p50000 -ddataflow_analyzer