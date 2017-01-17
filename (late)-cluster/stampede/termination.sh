#!/bin/sh
#set -x

CPATH=$1
DATAFLOW=$2

while [ -d "$CPATH/prov/di/$DATAFLOW" ]; do
   echo "   Waiting to load data into provenance database..." 
   sleep 60
done