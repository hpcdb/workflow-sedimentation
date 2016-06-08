rm -rf prov/di/*
rm -rf prov/pg/*

PGDIR=/Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/libmesh-sedimentation

echo "Dataflow - libMesh Sedimentation"
# Default mode
# Dataflow
java -jar ../dfa/PG-1.0.jar -dataflow -tag sedimentation

echo "Mesh Generation"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag meshGeneration
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation meshGeneration -name libmesh-sedimentation-opt::MeshGeneration -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshGeneration -tag imeshgeneration -type input
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshGeneration -tag omeshgeneration -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set imeshgeneration -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set imeshgeneration -name dim -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set imeshgeneration -name ncellx -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set imeshgeneration -name ncelly -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set imeshgeneration -name ncellz -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set imeshgeneration -name xmin -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set imeshgeneration -name ymin -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set imeshgeneration -name zmin -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set imeshgeneration -name xmax -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set imeshgeneration -name ymax -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set imeshgeneration -name zmax -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set imeshgeneration -name ref_interval -type numeric

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set omeshgeneration -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set omeshgeneration -name r_fraction -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set omeshgeneration -name c_fraction -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set omeshgeneration -name max_h_level -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshGeneration -set omeshgeneration -name hlevels -type numeric

echo "Create Equation Systems"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag createEquationSystems
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation createEquationSystems -name libmesh-sedimentation-opt::CreateEquationSystems -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation createEquationSystems -tag omeshgeneration -type input -dependency meshGeneration
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation createEquationSystems -tag ocreateequationsystems -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name Reynolds -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name Gr -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name Sc -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name Us -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name Diffusivity -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name xlock -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name alfa -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name theta -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name ex -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name ey -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name ez -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name c_factor -type numeric

echo "Get Maximum Iterations"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag getMaximumIterations
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation getMaximumIterations -name libmesh-sedimentation-opt::getMaximumIterations -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation getMaximumIterations -tag ocreateequationsystems -type input -dependency createEquationSystems
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation getMaximumIterations -tag ogetmaximumiterations -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterations -set ogetmaximumiterations -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterations -set ogetmaximumiterations -name dt -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterations -set ogetmaximumiterations -name tmax -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterations -set ogetmaximumiterations -name n_time_steps -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterations -set ogetmaximumiterations -name n_nonlinear_steps -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterations -set ogetmaximumiterations -name nonlinear_tolerance -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterations -set ogetmaximumiterations -name max_linear_iters -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterations -set ogetmaximumiterations -name max_r_steps -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterations -set ogetmaximumiterations -name write_interval -type numeric

echo "Solver Simulation to the Fluid"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag solverSimulationFluid
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation solverSimulationFluid -name libmesh-sedimentation-opt::SolverSimulationFluid -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation solverSimulationFluid -tag ogetmaximumiterations -type input -dependency getMaximumIterations
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation solverSimulationFluid -tag osolversimulationfluid -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFluid -set osolversimulationfluid -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFluid -set osolversimulationfluid -name t_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFluid -set osolversimulationfluid -name time -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFluid -set osolversimulationfluid -name r -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFluid -set osolversimulationfluid -name l -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFluid -set osolversimulationfluid -name n_linear_iterations -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFluid -set osolversimulationfluid -name final_linear_residual -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFluid -set osolversimulationfluid -name norm_delta -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFluid -set osolversimulationfluid -name norm_delta_u -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFluid -set osolversimulationfluid -name converged -type text

echo "Solver Simulation to the Sediments"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag solverSimulationSediments
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation solverSimulationSediments -name libmesh-sedimentation-opt::SolverSimulationSediments -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation solverSimulationSediments -tag osolversimulationfluid -type input -dependency solverSimulationFluid
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation solverSimulationSediments -tag osolversimulationsediments -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationSediments -set osolversimulationsediments -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationSediments -set osolversimulationsediments -name t_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationSediments -set osolversimulationsediments -name time -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationSediments -set osolversimulationsediments -name r -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationSediments -set osolversimulationsediments -name l -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationSediments -set osolversimulationsediments -name n_linear_iterations -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationSediments -set osolversimulationsediments -name final_linear_residual -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationSediments -set osolversimulationsediments -name norm_delta -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationSediments -set osolversimulationsediments -name norm_delta_u -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationSediments -set osolversimulationsediments -name converged -type text

echo "Mesh Refinement"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag meshRefinement
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation meshRefinement -name libmesh-sedimentation-opt::MeshRefinement -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshRefinement -tag osolversimulationsediments -type input -dependency solverSimulationSediments
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshRefinement -tag omeshrefinement -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name first_step_refinement -type text
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name before_n_active_elem -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name after_n_active_elem -type numeric

echo "Dataflow ingestion"
java -jar ../dfa/PG-1.0.jar -ingest -dataflow sedimentation

cp prov/pg/sedimentation/dataflow.json .