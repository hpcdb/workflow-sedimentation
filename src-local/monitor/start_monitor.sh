#!/bin/sh
#set -x
host=$1
directory=`pwd`
echo "cd $directory;python $directory/../monitor/Main.py -h localhost -p 22000"
ssh $host "killall python;cd $directory;rm performance.log;source $directory/../monitor/env/bin/activate;python $directory/../monitor/Main.py -h localhost -p 22000;deactivate" &