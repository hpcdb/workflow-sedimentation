#!/bin/sh
#set -x
dfa=$1
datapath=$2
rm -rf data
monetdbd stop $datapath
killall mserver5
killall monetdbd
monetdbd create $datapath
monetdbd set port=54321 $datapath
monetdbd start $datapath
monetdb create dataflow_analyzer
monetdb release dataflow_analyzer
mclient -p 54321 -d dataflow_analyzer $MONETDB/scripts/create-schema.sql
mclient -p 54321 -d dataflow_analyzer $dfa/database-script.sql
monetdbd stop $datapath