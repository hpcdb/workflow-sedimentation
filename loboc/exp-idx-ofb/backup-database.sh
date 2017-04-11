#!/bin/bash

echo "Starting server..."
monetdbd start /home/users/vitorss/simulation/sedimentation/exp-idx-ofb/data
monetdbd get all /home/users/vitorss/simulation/sedimentation/exp-idx-ofb/data

echo "Starting database..."
monetdb start dataflow_analyzer
monetdb status

echo "Backing up provenance database..."
rm prov-db.dump
mclient -p 54321 -d dataflow_analyzer --dump > prov-db.dump

#mclient -u monetdb -p54321 -ddataflow_analyzer