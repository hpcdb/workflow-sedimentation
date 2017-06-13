#!/bin/bash
cd fastbit-idx
# rm -rf *.csv
./delete.sh
cd ..
cd optimized-fastbit-idx
# rm -rf *.csv
./delete.sh
cd ..
cd positional-idx
# rm -rf *.csv
./delete.sh
cd ..
cd postgres-raw-idx
# rm -rf *.csv
./delete.sh
cd ..

rm indexing-sedimentation.o*
rm -rf output.log
rm -rf error.log