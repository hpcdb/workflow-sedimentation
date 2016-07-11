#!/bin/sh
#set -x

conf=$1
cpath=$2
lines=`cat $conf  | egrep -v "#"`
dir=`pwd`

for i in `echo $lines`; do 
  host=`echo $i`
  echo "cd $cpath;monetdb-start-all;java -jar /work/03664/silva/experiments/dfa/DI-1.0.jar -daemon start"
  ssh $host "cd $cpath;monetdb-start-all;java -jar /work/03664/silva/experiments/dfa/DI-1.0.jar -daemon start" &
done