#!/bin/sh
#set -x

CPATH=$1
DATAFLOW=$2

while [ -e "$CPATH/prov/di/$DATAFLOW/finish.token" ]; do
   echo "   Waiting to load data into provenance database..." 
   sleep 60
done