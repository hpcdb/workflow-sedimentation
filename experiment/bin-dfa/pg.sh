rm -rf di/*
rm -rf pg/*
echo "Dataflow - libMesh Sedimentation"
# Default mode
java -jar PG-1.0.jar -dataflow -tag sedimentation
java -jar PG-1.0.jar -transformation -dataflow sedimentation -tag equationSystems
java -jar PG-1.0.jar -program -dataflow sedimentation -transformation edgecfdpre -name EdgeCFD-Pre -filepath /Users/vitor/Documents/Repository/Thesis/EdgeCFD-trunk/workflow/bin/edgecfdpre

java -jar PG-1.0.jar -set -dataflow sedimentation -transformation edgecfdpre -tag iedgecfdpre -type input
java -jar PG-1.0.jar -set -dataflow sedimentation -transformation edgecfdpre -tag oedgecfdpre -type output

java -jar PG-1.0.jar -extractor -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -tag dl_pre -algorithm EXTRACTION:CSV
java -jar PG-1.0.jar -extractor -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -tag dl_in -algorithm EXTRACTION:CSV
java -jar PG-1.0.jar -extractor -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -tag dl_mat -algorithm EXTRACTION:CSV

java -jar PG-1.0.jar -extractor_combination -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -outer dl_pre -inner dl_in -keys {inn} -key_types {file}
java -jar PG-1.0.jar -extractor_combination -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -outer dl_pre -inner dl_mat -keys {mat} -key_types {file}

java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set iedgecfdpre -name mid -type numeric
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set iedgecfdpre -name mesh -type text
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set iedgecfdpre -name processes -type numeric

java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name mid -type numeric
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name mesh -type text
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name processes -type numeric

java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name edg -type file -extractor dl_pre
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name inn -type file -extractor dl_pre
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name mat -type file -extractor dl_pre
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name msh -type file -extractor dl_pre

java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name iprint -type numeric -extractor dl_in
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name dt -type numeric -extractor dl_in
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name etaf -type numeric -extractor dl_in

java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name visc -type numeric -extractor dl_mat
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name dens -type numeric -extractor dl_mat
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name kxx -type numeric -extractor dl_mat
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name kyy -type numeric -extractor dl_mat
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name kzz -type numeric -extractor dl_mat

echo "	- edgecfdsolver"
java -jar PG-1.0.jar -transformation -dataflow sedimentation -tag edgecfdsolver
java -jar PG-1.0.jar -program -dataflow sedimentation -transformation edgecfdsolver -name EdgeCFD-Solver -filepath /Users/vitor/Documents/Repository/Thesis/EdgeCFD-trunk/workflow/bin/edgecfdsolver

java -jar PG-1.0.jar -set -dataflow sedimentation -transformation edgecfdsolver -tag oedgecfdpre -type input -dependency edgecfdpre
java -jar PG-1.0.jar -set -dataflow sedimentation -transformation edgecfdsolver -tag oedgecfdsolver -type output

java -jar PG-1.0.jar -extractor -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -tag dl_paraview -algorithm EXTRACTION:CSV

java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name timestep -type numeric -extractor dl_paraview
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name pressure -type numeric -extractor dl_paraview
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name velocity_0 -type numeric -extractor dl_paraview
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name velocity_1 -type numeric -extractor dl_paraview
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name velocity_2 -type numeric -extractor dl_paraview
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name vtkvalidpointmask -type numeric -extractor dl_paraview
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name arc_length -type numeric -extractor dl_paraview
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name points_0 -type numeric -extractor dl_paraview
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name points_1 -type numeric -extractor dl_paraview
java -jar PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name points_2 -type numeric -extractor dl_paraview


echo "Dataflow ingestion"
java -jar PG-1.0.jar -ingest -dataflow edgecfd


