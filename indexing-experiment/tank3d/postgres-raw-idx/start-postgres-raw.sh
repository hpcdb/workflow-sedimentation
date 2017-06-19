#!/bin/bash
# LoboC
export POSTGRES_RAW_DIR=/home/users/vitorss/local/all/postgres-raw
cd $POSTGRES_RAW_DIR
rm -rf $POSTGRES_RAW_DIR/data
bin/initdb $POSTGRES_RAW_DIR/data
bin/postgres -D $POSTGRES_RAW_DIR/data >> db-output.log  2>> db-error.log &
sleep 5
# bin/createdb di_csv_idx_postgres_raw
cd -