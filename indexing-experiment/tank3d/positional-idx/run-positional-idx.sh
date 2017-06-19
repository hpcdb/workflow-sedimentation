export RDI_DIR='../bin/RDI-1.0-jar-with-dependencies.jar'
export RDI="java -jar $RDI_DIR"
export COLUMNS='[u:numeric:50,v:numeric:50,w:numeric:50,p:numeric:50,s:numeric:50,d:numeric:50,vtkValidPointMask:numeric,arc_length:numeric:50,Points0:numeric:50,Points1:numeric:50,Points2:numeric:50]'

# Start Timestamp
STARTTIME=`date +%s.%N`

for file in *.csv
do
  export fileName="${file%.*}"
  echo "${fileName}_csv_idx"
  export RDI_CMD="CSV:INDEX ${fileName}_csv_idx $PWD $file $COLUMNS -delimiter=\",\""
  $RDI $RDI_CMD
done

# End timestamp
ENDTIME=`date +%s.%N`

# Convert nanoseconds to milliseconds
# crudely by taking first 3 decimal places
ELAPSED_TIME=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print $1"."substr($2,1,3)}'`
echo "Elapsed time for Positional Indexing: $ELAPSED_TIME"