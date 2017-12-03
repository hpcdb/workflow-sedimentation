#!/bin/sh
#set -x

CPATH=$1

while [ ! -f "$CPATH/finish.token" ]; do
   echo "   Waiting to load data into provenance database..." 
   sleep 60
done

#mclient -p 54321 -d dataflow_analyzer --dump > prov-db.dump
