rm -rf prov/di/*
rm -rf prov/pg/*
echo "Dataflow - libMesh Sedimentation"
# Default mode
java -jar ../dfa/PG-1.0.jar -dataflow -tag sedimentation
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag meshRefinement
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation meshRefinement -name libmesh-sedimentation-opt::MeshRefinement -filepath /Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/libmesh-sedimentation

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshRefinement -tag imeshrefinement -type input
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshRefinement -tag omeshrefinement -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set imeshrefinement -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set imeshrefinement -name dim -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set imeshrefinement -name ncellx -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set imeshrefinement -name ncelly -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set imeshrefinement -name ncellz -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set imeshrefinement -name xmin -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set imeshrefinement -name ymin -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set imeshrefinement -name zmin -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set imeshrefinement -name xmax -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set imeshrefinement -name ymax -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set imeshrefinement -name zmax -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set imeshrefinement -name ref_interval -type numeric

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name r_fraction -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name c_fraction -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name max_h_level -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name hlevels -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name first_step_refinement -type text






# java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -tag dl_pre -algorithm EXTRACTION:CSV
# java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -tag dl_in -algorithm EXTRACTION:CSV
# java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -tag dl_mat -algorithm EXTRACTION:CSV

# java -jar ../dfa/PG-1.0.jar -extractor_combination -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -outer dl_pre -inner dl_in -keys {inn} -key_types {file}
# java -jar ../dfa/PG-1.0.jar -extractor_combination -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -outer dl_pre -inner dl_mat -keys {mat} -key_types {file}

# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set iedgecfdpre -name mid -type numeric
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set iedgecfdpre -name mesh -type text
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set iedgecfdpre -name processes -type numeric

# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name mid -type numeric
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name mesh -type text
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name processes -type numeric

# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name edg -type file -extractor dl_pre
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name inn -type file -extractor dl_pre
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name mat -type file -extractor dl_pre
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name msh -type file -extractor dl_pre

# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name iprint -type numeric -extractor dl_in
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name dt -type numeric -extractor dl_in
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name etaf -type numeric -extractor dl_in

# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name visc -type numeric -extractor dl_mat
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name dens -type numeric -extractor dl_mat
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name kxx -type numeric -extractor dl_mat
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name kyy -type numeric -extractor dl_mat
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdpre -set oedgecfdpre -name kzz -type numeric -extractor dl_mat

# echo "	- edgecfdsolver"
# java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag edgecfdsolver
# java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation edgecfdsolver -name EdgeCFD-Solver -filepath /Users/vitor/Documents/Repository/Thesis/EdgeCFD-trunk/workflow/bin/edgecfdsolver

# java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation edgecfdsolver -tag oedgecfdpre -type input -dependency edgecfdpre
# java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation edgecfdsolver -tag oedgecfdsolver -type output

# java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -tag dl_paraview -algorithm EXTRACTION:CSV

# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name timestep -type numeric -extractor dl_paraview
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name pressure -type numeric -extractor dl_paraview
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name velocity_0 -type numeric -extractor dl_paraview
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name velocity_1 -type numeric -extractor dl_paraview
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name velocity_2 -type numeric -extractor dl_paraview
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name vtkvalidpointmask -type numeric -extractor dl_paraview
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name arc_length -type numeric -extractor dl_paraview
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name points_0 -type numeric -extractor dl_paraview
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name points_1 -type numeric -extractor dl_paraview
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation edgecfdsolver -set oedgecfdsolver -name points_2 -type numeric -extractor dl_paraview


echo "Dataflow ingestion"
java -jar ../dfa/PG-1.0.jar -ingest -dataflow sedimentation


