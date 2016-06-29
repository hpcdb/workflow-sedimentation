rm -rf prov/di/*
rm -rf prov/pg/*

# MacOS
PGDIR=/Users/vitor/Documents/Repository/Thesis/WorkflowSedimentation/sedimentation
# docker
# PGDIR=/experiment/libmesh-sedimentation
# Virtual Box - mint
# PGDIR=/media/sf_shared/libmesh-sedimentation

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
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation createEquationSystems -set ocreateequationsystems -name fopc -type numeric
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
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterations -set ogetmaximumiterations -name xdmf -type text

# 3D Analysis
# echo "Init Data Extraction"
# java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag initDataExtraction
# java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation initDataExtraction -name libmesh-sedimentation-opt::InitDataExtraction -filepath $PGDIR

# java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation initDataExtraction -tag ogetmaximumiterations -type input -dependency getMaximumIterations
# java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation initDataExtraction -tag oinitdataextraction -type output

# java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation initDataExtraction -set oinitdataextraction -tag irde -algorithm EXTRACTION:PROGRAM

# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation initDataExtraction -set oinitdataextraction -name simulationID -type numeric
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation initDataExtraction -set oinitdataextraction -name time_step -type numeric
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation initDataExtraction -set oinitdataextraction -name xdmf -type file
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation initDataExtraction -set oinitdataextraction -name u -type numeric -extractor irde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation initDataExtraction -set oinitdataextraction -name v -type numeric -extractor irde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation initDataExtraction -set oinitdataextraction -name p -type numeric -extractor irde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation initDataExtraction -set oinitdataextraction -name s -type numeric -extractor irde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation initDataExtraction -set oinitdataextraction -name d -type numeric -extractor irde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation initDataExtraction -set oinitdataextraction -name points0 -type numeric -extractor irde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation initDataExtraction -set oinitdataextraction -name points1 -type numeric -extractor irde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation initDataExtraction -set oinitdataextraction -name points2 -type numeric -extractor irde

echo "Initial Horizontal Line Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag iLine0Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation iLine0Extraction -name libmesh-sedimentation-opt::iLine0Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine0Extraction -tag ogetmaximumiterations -type input -dependency getMaximumIterations
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine0Extraction -tag oline0iextraction -type output

java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -tag iline0 -algorithm EXTRACTION:PROGRAM

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name u -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name v -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name w -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name p -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name s -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name d -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name vtkvalidpointmask -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name arc_length -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name points0 -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name points1 -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name points2 -type numeric -extractor iline0

echo "Initial Vertical Line 1 Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag iLine1Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation iLine1Extraction -name libmesh-sedimentation-opt::iLine1Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine1Extraction -tag ogetmaximumiterations -type input -dependency getMaximumIterations
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine1Extraction -tag oline1iextraction -type output

java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -tag iline1 -algorithm EXTRACTION:PROGRAM

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name u -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name v -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name w -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name p -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name s -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name d -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name vtkvalidpointmask -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name arc_length -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name points0 -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name points1 -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name points2 -type numeric -extractor iline1

echo "Initial Vertical Line 2 Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag iLine2Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation iLine2Extraction -name libmesh-sedimentation-opt::iLine2Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine2Extraction -tag ogetmaximumiterations -type input -dependency getMaximumIterations
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine2Extraction -tag oline2iextraction -type output

java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -tag iline2 -algorithm EXTRACTION:PROGRAM

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name u -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name v -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name w -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name p -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name s -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name d -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name vtkvalidpointmask -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name arc_length -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name points0 -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name points1 -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name points2 -type numeric -extractor iline2

echo "Initial Vertical Line 3 Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag iLine3Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation iLine3Extraction -name libmesh-sedimentation-opt::iLine3Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine3Extraction -tag ogetmaximumiterations -type input -dependency getMaximumIterations
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine3Extraction -tag oline3iextraction -type output

java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -tag iline3 -algorithm EXTRACTION:PROGRAM

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name u -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name v -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name w -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name p -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name s -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name d -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name vtkvalidpointmask -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name arc_length -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name points0 -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name points1 -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name points2 -type numeric -extractor iline3

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
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name t_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name before_n_active_elem -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name after_n_active_elem -type numeric

echo "Mesh Writer"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag meshWriter
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation meshWriter -name libmesh-sedimentation-opt::MeshWriter -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshWriter -tag osolversimulationsediments -type input -dependency solverSimulationSediments
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshWriter -tag omeshwriter -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshWriter -set omeshwriter -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshWriter -set omeshwriter -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshWriter -set omeshwriter -name xdmf -type file

# 2D Analysis
# echo "Data Extraction"
# java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag dataExtraction
# java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation dataExtraction -name libmesh-sedimentation-opt::DataExtraction -filepath $PGDIR

# java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation dataExtraction -tag omeshwriter -type input -dependency meshWriter
# java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation dataExtraction -tag odataextraction -type output

# java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation dataExtraction -set odataextraction -tag rde -algorithm EXTRACTION:PROGRAM

# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation dataExtraction -set odataextraction -name simulationID -type numeric
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation dataExtraction -set odataextraction -name time_step -type numeric
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation dataExtraction -set odataextraction -name xdmf -type file
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation dataExtraction -set odataextraction -name u -type numeric -extractor rde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation dataExtraction -set odataextraction -name v -type numeric -extractor rde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation dataExtraction -set odataextraction -name p -type numeric -extractor rde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation dataExtraction -set odataextraction -name s -type numeric -extractor rde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation dataExtraction -set odataextraction -name d -type numeric -extractor rde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation dataExtraction -set odataextraction -name points0 -type numeric -extractor rde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation dataExtraction -set odataextraction -name points1 -type numeric -extractor rde
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation dataExtraction -set odataextraction -name points2 -type numeric -extractor rde

# 3D Analysis
echo "Horizontal Line Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag line0Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation line0Extraction -name libmesh-sedimentation-opt::Line0Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line0Extraction -tag omeshwriter -type input -dependency meshWriter
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line0Extraction -tag oline0extraction -type output

java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line0Extraction -set oline0extraction -tag line0 -algorithm EXTRACTION:PROGRAM

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name u -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name v -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name w -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name p -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name s -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name d -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name vtkvalidpointmask -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name arc_length -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name points0 -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name points1 -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name points2 -type numeric -extractor line0

echo "Vertical Line 1 Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag line1Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation line1Extraction -name libmesh-sedimentation-opt::Line1Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line1Extraction -tag omeshwriter -type input -dependency meshWriter
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line1Extraction -tag oline1extraction -type output

java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line1Extraction -set oline1extraction -tag line1 -algorithm EXTRACTION:PROGRAM

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name u -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name v -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name w -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name p -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name s -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name d -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name vtkvalidpointmask -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name arc_length -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name points0 -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name points1 -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name points2 -type numeric -extractor line1

echo "Vertical Line 2 Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag line2Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation line2Extraction -name libmesh-sedimentation-opt::Line2Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line2Extraction -tag omeshwriter -type input -dependency meshWriter
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line2Extraction -tag oline2extraction -type output

java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line2Extraction -set oline2extraction -tag line2 -algorithm EXTRACTION:PROGRAM

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name u -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name v -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name w -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name p -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name s -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name d -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name vtkvalidpointmask -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name arc_length -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name points0 -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name points1 -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name points2 -type numeric -extractor line2

echo "Vertical Line 3 Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag line3Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation line3Extraction -name libmesh-sedimentation-opt::Line3Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line3Extraction -tag omeshwriter -type input -dependency meshWriter
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line3Extraction -tag oline3extraction -type output

java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line3Extraction -set oline3extraction -tag line3 -algorithm EXTRACTION:PROGRAM

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name u -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name v -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name w -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name p -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name s -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name d -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name vtkvalidpointmask -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name arc_length -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name points0 -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name points1 -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name points2 -type numeric -extractor line3

# echo "Compute Statistics"
# java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag computeStatistics
# java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation computeStatistics -name libmesh-sedimentation-opt::ComputeStatistics -filepath $PGDIR

# java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation computeStatistics -tag omeshwriter -type input -dependency meshWriter
# java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation computeStatistics -tag ostatistics -type output

# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeStatistics -set ostatistics -name simulationID -type numeric
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeStatistics -set ostatistics -name time_step -type numeric
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeStatistics -set ostatistics -name xdmf -type file
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeStatistics -set ostatistics -name sediments_amount -type numeric

echo "Mesh Aggregator"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag meshAggregator
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation meshAggregator -name libmesh-sedimentation-opt::MeshAggregator -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshAggregator -tag omeshwriter -type input -dependency meshWriter
# java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshAggregator -tag ostatistics -type input -dependency computeStatistics
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshAggregator -tag omeshaggregator -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshAggregator -set omeshaggregator -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshAggregator -set omeshaggregator -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshAggregator -set omeshaggregator -name n_processors -type numeric

echo "Dataflow ingestion"
java -jar ../dfa/PG-1.0.jar -ingest -dataflow sedimentation

cp prov/pg/sedimentation/dataflow.json .



