#!/bin/sh
#set -x

conf=$1
cpath=$2
datapath=$3
lines=`cat $conf  | egrep -v "#"`
dir=`pwd`

for i in `echo $lines`; do 
  host=`echo $i`
  echo "cd $cpath;killall monetdb;killall monetdbd;killall mserver5;monetdbd set port=54321 $datapath; monetdbd start $datapath;monetdbd get all $datapath;monetdb start dataflow_analyzer;monetdb status;java -jar /home/user4/sedimentation/dfa/DI-1.0.jar -daemon start"
  ssh $host "cd $cpath;monetdb-stop-all;killall monetdb;killall monetdbd;killall mserver5;monetdbd set port=54321 $datapath; monetdbd start $datapath;monetdbd get all $datapath;monetdb start dataflow_analyzer;monetdb status;java -jar /home/user4/sedimentation/dfa/DI-1.0.jar -daemon start  > db-out.txt 2> db-err.txt" &
done