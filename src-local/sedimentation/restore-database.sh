#!/bin/bash
clear
monetdb-start-all
# Restart MonetDB
cd ../dfa
./restart-monetdb.sh
cd ../sedimentation
# delete files from the previous execution
# ./delete.sh