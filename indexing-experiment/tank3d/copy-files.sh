#!/bin/bash
cd fastbit-idx
./copy-files.sh
cd ..
cd optimized-fastbit-idx
./copy-files.sh
cd ..
cd positional-idx
./copy-files.sh
cd ..
cd postgres-raw-idx
./copy-files.sh
cd ..