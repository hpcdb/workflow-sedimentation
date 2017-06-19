export RDI_DIR='../bin/RDI-1.0-jar-with-dependencies.jar'
export RDI="java -jar $RDI_DIR"

export FASTBIT_BIN='/home/users/vitorss/local/all/fastbit/bin/'
export COLUMNS='[u:numeric:5,v:numeric:5,w:numeric:5,p:numeric:5,s:numeric:5,d:numeric:5,vtkValidPointMask:numeric,arc_length:numeric:5,Points0:numeric:5,Points1:numeric:5,Points2:numeric:5]'

export RDI_CMD="$PWD \"*.csv\" $COLUMNS -bin=\"$FASTBIT_BIN\" -delimiter=\",\""

export INDEX_TYPE="FASTBIT:INDEX"
export INDEX_NAME_BASE="fastbit_idx"

# Start Timestamp
STARTTIME=`date +%s.%N`

echo ${INDEX_NAME_BASE}_e.equality_c.uncompressdensebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.equality_c.uncompressdensebitmaps $RDI_CMD -option=[e:equality,c:uncompressdensebitmaps]
echo ${INDEX_NAME_BASE}_e.equality_c.uncompresslargebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.equality_c.uncompresslargebitmaps $RDI_CMD -option=[e:equality,c:uncompresslargebitmaps]
echo ${INDEX_NAME_BASE}_e.equality_c.uncompressall 
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.equality_c.uncompressall $RDI_CMD -option=[e:equality,c:uncompressall]
echo ${INDEX_NAME_BASE}_e.range_c.uncompressdensebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.range_c.uncompressdensebitmaps $RDI_CMD -option=[e:range,c:uncompressdensebitmaps]
echo ${INDEX_NAME_BASE}_e.range_c.uncompresslargebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.range_c.uncompresslargebitmaps $RDI_CMD -option=[e:range,c:uncompresslargebitmaps]
echo ${INDEX_NAME_BASE}_e.range_c.uncompressall 
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.range_c.uncompressall $RDI_CMD -option=[e:range,c:uncompressall]
echo ${INDEX_NAME_BASE}_e.interval_c.uncompressdensebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.interval_c.uncompressdensebitmaps $RDI_CMD -option=[e:interval,c:uncompressdensebitmaps]
echo ${INDEX_NAME_BASE}_e.interval_c.uncompresslargebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.interval_c.uncompresslargebitmaps $RDI_CMD -option=[e:interval,c:uncompresslargebitmaps]
echo ${INDEX_NAME_BASE}_e.interval_c.uncompressall
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.interval_c.uncompressall $RDI_CMD -option=[e:interval,c:uncompressall]
echo ${INDEX_NAME_BASE}_e.ncomp=ddd_c.uncompressdensebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.ncomp=ddd_c.uncompressdensebitmaps $RDI_CMD -option=[e:ncomp=ddd,c:uncompressdensebitmaps]
echo ${INDEX_NAME_BASE}_e.ncomp=ddd_c.uncompresslargebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.ncomp=ddd_c.uncompresslargebitmaps $RDI_CMD -option=[e:ncomp=ddd,c:uncompresslargebitmaps]
echo ${INDEX_NAME_BASE}_e.ncomp=ddd_c.uncompressall
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.ncomp=ddd_c.uncompressall $RDI_CMD -option=[e:ncomp=ddd,c:uncompressall]
echo ${INDEX_NAME_BASE}_e.equality-equality_c.uncompressdensebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.equality-equality_c.uncompressdensebitmaps $RDI_CMD -option=[e:equality-equality,c:uncompressdensebitmaps]
echo ${INDEX_NAME_BASE}_e.equality-equality_c.uncompresslargebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.equality-equality_c.uncompresslargebitmaps $RDI_CMD -option=[e:equality-equality,c:uncompresslargebitmaps]
echo ${INDEX_NAME_BASE}_e.equality-equality_c.uncompressall
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.equality-equality_c.uncompressall $RDI_CMD -option=[e:equality-equality,c:uncompressall]
echo ${INDEX_NAME_BASE}_e.interval-equality_c.uncompressdensebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.interval-equality_c.uncompressdensebitmaps $RDI_CMD -option=[e:interval-equality,c:uncompressdensebitmaps]
echo ${INDEX_NAME_BASE}_e.interval-equality_c.uncompresslargebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.interval-equality_c.uncompresslargebitmaps $RDI_CMD -option=[e:interval-equality,c:uncompresslargebitmaps]
echo ${INDEX_NAME_BASE}_e.interval-equality_c.uncompressall
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.interval-equality_c.uncompressall $RDI_CMD -option=[e:interval-equality,c:uncompressall]
echo ${INDEX_NAME_BASE}_e.range-equality_c.uncompressdensebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.range-equality_c.uncompressdensebitmaps $RDI_CMD -option=[e:range-equality,c:uncompressdensebitmaps]
echo ${INDEX_NAME_BASE}_e.range-equality_c.uncompresslargebitmaps
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.range-equality_c.uncompresslargebitmaps $RDI_CMD -option=[e:range-equality,c:uncompresslargebitmaps]
echo ${INDEX_NAME_BASE}_e.range-equality_c.uncompressall
$RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.range-equality_c.uncompressall $RDI_CMD -option=[e:range-equality,c:uncompressall]

# End timestamp
ENDTIME=`date +%s.%N`

# Convert nanoseconds to milliseconds
# crudely by taking first 3 decimal places
ELAPSED_TIME=`echo "$ENDTIME - $STARTTIME" | bc | awk -F"." '{print $1"."substr($2,1,3)}'`
echo "Elapsed time for FastBit Indexing: $ELAPSED_TIME"