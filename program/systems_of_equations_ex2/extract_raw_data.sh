export DFANALYZER_DIR="/home/vitor/Documents/dev/workflow-sedimentation/program/systems_of_equations_ex2/dfanalyzer"
java -jar $DFANALYZER_DIR/bin/RDE PROGRAM:EXTRACT extractor . "$PARAVIEW_DIR/bin/pvpython ./script/exodus_data_extraction.py 1" [u:NUMERIC,v:NUMERIC,w:NUMERIC,p:NUMERIC,x:NUMERIC,y:NUMERIC,z:NUMERIC]
