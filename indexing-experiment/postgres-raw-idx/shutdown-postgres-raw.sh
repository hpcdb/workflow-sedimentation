#!/bin/bash
pkill postgres
pkill postgres
sleep 10
export POSTGRES_RAW_DIR=/home/users/vitorss/local/all/postgres-raw
cd $POSTGRES_RAW_DIR
rm -rf data
rm -rf db-output.log
rm -rf db-error.log
cd -