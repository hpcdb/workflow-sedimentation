#!/bin/sh
#set -x
conf=$1
cpath=$2
datapath=$3
lines=`cat $conf  | egrep -v "#"`
dir=`pwd`
SIMULATION_DIR=/media/vitor/data-linux/dev/workflow-sedimentation/src-local

for i in `echo $lines`; do 
  host=`echo $i`
  echo "cd $cpath;monetdbd start $datapath;monetdbd get all $datapath;monetdb start dataflow_analyzer;monetdb status;$JAVA_HOME/bin/java -jar $SIMULATION_DIR/dfa/DI-1.0.jar -daemon start"
  cd $cpath;monetdbd start $datapath;monetdbd get all $datapath;monetdb start dataflow_analyzer;monetdb status;$JAVA_HOME/bin/java -jar $SIMULATION_DIR/dfa/DI-1.0.jar -daemon start
  # echo "cd $cpath;killall monetdb;killall monetdbd;killall mserver5;monetdbd set port=54321 $datapath; monetdbd start $datapath;monetdbd get all $datapath;monetdb start dataflow_analyzer;monetdb status;java -jar /home/users/vitorss/simulation/sedimentation/dfa/DI-1.0.jar -daemon start"
  # ssh $host "cd $cpath;monetdbd stop $datapath;killall monetdb;killall monetdbd;killall mserver5;monetdbd set port=54321 $datapath; monetdbd start $datapath;monetdbd get all $datapath;monetdb start dataflow_analyzer;monetdb status;java -jar /home/users/vitorss/simulation/sedimentation/dfa/DI-1.0.jar -daemon start  > db-out.txt 2> db-err.txt" &
done