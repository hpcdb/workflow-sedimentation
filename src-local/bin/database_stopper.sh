#!/bin/sh
#set -x
conf=$1
cpath=$2
datapath=$3
lines=`cat $conf  | egrep -v "#"`
dir=`pwd`

for i in `echo $lines`; do 
  host=`echo $i`
  echo "cd $cpath;monetdbd stop $datapath\n"
  "cd $cpath;monetdbd stop $datapath"
done