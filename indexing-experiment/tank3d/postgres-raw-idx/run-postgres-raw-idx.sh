export RDI_DIR='../bin/RDI-1.0-jar-with-dependencies.jar'
export RDI="java -jar $RDI_DIR"

export BIN='/home/users/vitorss/local/all/postgres-raw/bin/'
export PG_DATA='/home/users/vitorss/local/all/postgres-raw/data'
export COLUMNS='[u:numeric:50,v:numeric:50,w:numeric:50,p:numeric:50,s:numeric:50,d:numeric:50,vtkValidPointMask:numeric,arc_length:numeric:50,Points0:numeric:50,Points1:numeric:50,Points2:numeric:50]'

export RDI_CMD="POSTGRES_RAW:INDEX di_csv_idx_postgres_raw $PWD \"*.csv\" $COLUMNS -bin=\"$BIN\" -delimiter=\",\" -pgData=\"$PG_DATA\" -option=[header:true]"

# Start Timestamp
STARTTIME=`date +%s.%N`

$RDI $RDI_CMD

# End timestamp
ENDTIME=`date +%s.%N`

# Convert nanoseconds to milliseconds
# crudely by taking first 3 decimal places
ELAPSED_TIME=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print $1"."substr($2,1,3)}'`
echo "Elapsed time for PostgresRaw Indexing: $ELAPSED_TIME"
