#!/bin/bash
for EXP_DIR in {1..6}
do
	echo "Executing experiment $EXP_DIR"
	# organize input files
	cp ../monitoring/lock_necker3D_pc11-exp$EXP_DIR.in lock_necker3D_pc11.in
	# kill all previous processes
	killall sediment-opt
	# solver execution
	./delete.sh
	./execute-solver.sh
	# experiment directory
	rm -rf ../exp-$EXP_DIR
	mkdir ../exp-$EXP_DIR
	cp -rf * ../exp-$EXP_DIR
	clear
done
