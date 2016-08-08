#!/bin/bash
clear
monetdb-start-all
# delete files from the previous execution
./delete.sh
# Start daemon process to Data Ingestor
java -jar ../dfa/DI-1.0.jar -daemon start

