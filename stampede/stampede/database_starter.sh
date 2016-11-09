#!/bin/sh
#set -x

conf=$1
cpath=$2
datapath=$3
lines=`cat $conf  | egrep -v "#"`
dir=`pwd`

for i in `echo $lines`; do 
  host=`echo $i`
  echo "cd $cpath;monetdbd start $datapath;monetdbd get all $datapath;monetdb start dataflow_analyzer;monetdb status;java -jar /work/03664/silva/experiments/dfa/DI-1.0.jar -daemon start >>solver.output 2>>solver.error"
  ssh $host "cd $cpath;monetdbd start $datapath;monetdbd get all $datapath;monetdb start dataflow_analyzer;monetdb status;java -jar /work/03664/silva/experiments/dfa/DI-1.0.jar -daemon start >>solver.output 2>>solver.error" &
done