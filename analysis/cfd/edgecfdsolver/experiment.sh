#!/bin/bash
# sequential
%=WFDIR%/../bin/edgecfdsolver < input.dat
# openmp
# %=WFDIR%/../bin/edgecfdsolver-openmp < input.dat

python %=WFDIR%/../bin/generate_erelation_solver.py %=MID% %=MESH%

# data selection from xmf files
meshName=`basename %=MESH%`
numberOfTimeSteps=`find . -maxdepth 1 -name "$meshName*.h5" | wc -l`
numberOfTimeSteps=`expr $numberOfTimeSteps - 1`

for timeStep in $(seq 0 1 `echo $numberOfTimeSteps`)
do
   /usr/bin/python %=WFDIR%/../bin/edit_Configuration_File-XDMF.py `basename %=MESH%` 1 `echo $timeStep` [0.5,0.5,0.0] [0.5,0.5,1.0] paraview.csv FALSE
   { time $PARAVIEW/bin/pvpython %=WFDIR%/../bin/XDMF-Filter.py config.properties ; } 2>> filter.log
done

# /usr/bin/python %=WFDIR%/../bin/edit_Configuration_File-XDMF.py `basename %=MESH%` 1 $numberOfTimeSteps [0.5,0.5,0.0] [0.5,0.5,1.0] paraview.csv TRUE
# { time $PARAVIEW/bin/pvpython %=WFDIR%/../bin/XDMF-Filter.py config.properties ; } 2>> filter.log

PATH=`pwd`
{ time /usr/bin/python %=WFDIR%/../bin/Stp-Extractor.py %=MESH%_1.stp stp.csv ; } 2>> filter.log

{ time /usr/bin/python %=WFDIR%/../bin/SolverJoinData.py stp.csv paraview.csv solver.csv ; } 2>> filter.log