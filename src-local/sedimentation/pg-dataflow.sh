rm -rf prov/di/*
rm -rf prov/pg/*

# environments: xps-nacad,xps-home,inspiron-laptop, macos
environment="macos"
PGDIR=""

if [ "$environment" == "xps-nacad" ]; then
	PGDIR="/home/vitor/Documents/dev/workflow-sedimentation/src-local/sedimentation"
elif [ "$environment" == "xps-home" ]; then
	PGDIR="/home/vitor/Documents/dev/workflow-sedimentation/src-local/sedimentation"	
elif [ "$environment" == "inspiron-laptop" ]; then
	PGDIR="/media/vitor/data-linux/dev/workflow-sedimentation/src-local/sedimentation"	
elif [ "$environment" == "macos" ]; then
	PGDIR="/Users/vitor/Documents/repository/workflow-sedimentation/src-local/sedimentation"
fi

# Ubuntu - Inspiron
#PGDIR=/media/vitor/data-linux/dev/workflow-sedimentation/src-local/sedimentation
# docker
# PGDIR=/shared/experiment/libmesh-sedimentation
# Virtual Box - mint
#PGDIR=/home/vitor/dev/workflow-sedimentation/libmesh-sedimentation
# Ubuntu
# PGDIR=/home/vitor/Documents/dev/workflow-sedimentation/libmesh-sedimentation
# Stampede
#PGDIR=/work/03664/silva/experiments/sedimentation
# LoboC
# PGDIR=/home/users/vitorss/simulation/sedimentation/libmesh-sedimentation
# Mint
# PGDIR=/home/vitor/Documents/dev/workflow-sedimentation/src-local/libmesh-sedimentation

# Sedimentation Solver
dimension="2"
access="extraction"
cartridge="csv"

echo "Dataflow - libMesh Sedimentation"
# Default mode
# Dataflow
java -jar ../dfa/PG-1.0.jar -dataflow -tag sedimentation

echo "Input Mesh"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag inputMesh
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation inputMesh -name libmesh-sedimentation-opt::InputMesh -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation inputMesh -tag iinputmesh -type input
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation inputMesh -tag oinputmesh -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation inputMesh -set iinputmesh -name simulationID -type numeric

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation inputMesh -set oinputmesh -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation inputMesh -set oinputmesh -name dim -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation inputMesh -set oinputmesh -name mesh_file -type text
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation inputMesh -set oinputmesh -name restart_control -type text

echo "AMR Config"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag amrConfig
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation amrConfig -name libmesh-sedimentation-opt::AMRConfig -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation amrConfig -tag oinputmesh -type input -dependency inputMesh
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation amrConfig -tag oamrconfig -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation amrConfig -set oamrconfig -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation amrConfig -set oamrconfig -name r_fraction -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation amrConfig -set oamrconfig -name c_fraction -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation amrConfig -set oamrconfig -name max_h_level -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation amrConfig -set oamrconfig -name hlevels -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation amrConfig -set oamrconfig -name first_step_refinement -type text
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation amrConfig -set oamrconfig -name amrc_flow_transp -type text
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation amrConfig -set oamrconfig -name ref_interval -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation amrConfig -set oamrconfig -name max_r_steps -type numeric

echo "Create Equation Systems"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag createEquationSystems
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation createEquationSystems -name libmesh-sedimentation-opt::CreateEquationSystems -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation createEquationSystems -tag oinputmesh -type input -dependency inputMesh
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

echo "Time Step Control Config"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag timeStepControlConfig
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation timeStepControlConfig -name libmesh-sedimentation-opt::TimeStepControlConfig -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation timeStepControlConfig -tag oinputmesh -type input -dependency inputMesh
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation timeStepControlConfig -tag otimestepcontrolconfig -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name model_name -type text
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name dt_min -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name dt_max -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name tol_u -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name tol_s -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name kp -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name ki -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name kd -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name nsa_max -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name nsa_target_flow -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name nsa_target_transport -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name nsa_limit_flow -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name nsa_limit_transport -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name mult_factor_max -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name mult_factor_min -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name pc11_theta -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name alpha -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name k_exp -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name s_min -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name s_max -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name reduct_factor -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation timeStepControlConfig -set otimestepcontrolconfig -name complete_flow_norm -type text

echo "IO Config"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag ioConfig
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation ioConfig -name libmesh-sedimentation-opt::IOConfig -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation ioConfig -tag oinputmesh -type input -dependency inputMesh
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation ioConfig -tag oioconfig -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation ioConfig -set oioconfig -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation ioConfig -set oioconfig -name dpath -type text
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation ioConfig -set oioconfig -name rname -type text
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation ioConfig -set oioconfig -name write_interval -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation ioConfig -set oioconfig -name catalyst_interval -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation ioConfig -set oioconfig -name write_restart -type text

echo "Get Maximum Iterations to the Flow"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag getMaximumIterationsToFlow
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation getMaximumIterationsToFlow -name libmesh-sedimentation-opt::getMaximumIterationsToFlow -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation getMaximumIterationsToFlow -tag ocreateequationsystems -type input -dependency createEquationSystems
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation getMaximumIterationsToFlow -tag ogetmaximumiterationstoflow -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToFlow -set ogetmaximumiterationstoflow -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToFlow -set ogetmaximumiterationstoflow -name dt -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToFlow -set ogetmaximumiterationstoflow -name tmax -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToFlow -set ogetmaximumiterationstoflow -name n_time_steps -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToFlow -set ogetmaximumiterationstoflow -name n_nonlinear_steps -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToFlow -set ogetmaximumiterationstoflow -name nonlinear_tolerance -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToFlow -set ogetmaximumiterationstoflow -name max_linear_iters -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToFlow -set ogetmaximumiterationstoflow -name xdmf -type text

echo "Get Maximum Iterations to the Transport"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag getMaximumIterationsToTransport
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation getMaximumIterationsToTransport -name libmesh-sedimentation-opt::getMaximumIterationsToTransport -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation getMaximumIterationsToTransport -tag ogetmaximumiterationstoflow -type input -dependency getMaximumIterationsToFlow
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation getMaximumIterationsToTransport -tag ogetmaximumiterationstotransport -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToTransport -set ogetmaximumiterationstotransport -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToTransport -set ogetmaximumiterationstotransport -name dt -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToTransport -set ogetmaximumiterationstotransport -name tmax -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToTransport -set ogetmaximumiterationstotransport -name n_time_steps -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToTransport -set ogetmaximumiterationstotransport -name n_nonlinear_steps -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToTransport -set ogetmaximumiterationstotransport -name nonlinear_tolerance -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToTransport -set ogetmaximumiterationstotransport -name max_linear_iters -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation getMaximumIterationsToTransport -set ogetmaximumiterationstotransport -name xdmf -type text

# Analysis
echo "Initial Horizontal Line Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag iLine0Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation iLine0Extraction -name libmesh-sedimentation-opt::iLine0Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine0Extraction -tag ogetmaximumiterationstotransport -type input -dependency getMaximumIterationsToTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine0Extraction -tag oline0iextraction -type output

echo "Initial Vertical Line 1 Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag iLine1Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation iLine1Extraction -name libmesh-sedimentation-opt::iLine1Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine1Extraction -tag ogetmaximumiterationstotransport -type input -dependency getMaximumIterationsToTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine1Extraction -tag oline1iextraction -type output

echo "Initial Vertical Line 2 Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag iLine2Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation iLine2Extraction -name libmesh-sedimentation-opt::iLine2Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine2Extraction -tag ogetmaximumiterationstotransport -type input -dependency getMaximumIterationsToTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine2Extraction -tag oline2iextraction -type output

echo "Initial Vertical Line 3 Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag iLine3Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation iLine3Extraction -name libmesh-sedimentation-opt::iLine3Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine3Extraction -tag ogetmaximumiterationstotransport -type input -dependency getMaximumIterationsToTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iLine3Extraction -tag oline3iextraction -type output

if [ "$access" == "extraction" ]; then
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -tag iline0 -algorithm EXTRACTION:CSV
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -tag iline1 -algorithm EXTRACTION:CSV
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -tag iline2 -algorithm EXTRACTION:CSV
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -tag iline3 -algorithm EXTRACTION:CSV
elif [ "$access" == "indexing" ] && [ "$cartridge" == "csv" ]; then
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -tag iline0 -algorithm INDEXING:CSV
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -tag iline1 -algorithm INDEXING:CSV
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -tag iline2 -algorithm INDEXING:CSV
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -tag iline3 -algorithm INDEXING:CSV
elif [ "$access" == "indexing" ] && [ "$cartridge" == "fastbit" ]; then
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -tag iline0 -algorithm INDEXING:FASTBIT
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -tag iline1 -algorithm INDEXING:FASTBIT
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -tag iline2 -algorithm INDEXING:FASTBIT
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -tag iline3 -algorithm INDEXING:FASTBIT
elif [ "$access" == "indexing" ] && [ "$cartridge" == "optimized_fastbit" ]; then
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -tag iline0 -algorithm INDEXING:OPTIMIZED_FASTBIT
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -tag iline1 -algorithm INDEXING:OPTIMIZED_FASTBIT
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -tag iline2 -algorithm INDEXING:OPTIMIZED_FASTBIT
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -tag iline3 -algorithm INDEXING:OPTIMIZED_FASTBIT
fi

# line 0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name u -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name v -type numeric -extractor iline0
if [ "$dimension" == "3" ]; then
	java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name w -type numeric -extractor iline0
fi
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name p -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name s -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name d -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name r -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name vtkvalidpointmask -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name arc_length -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name points0 -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name points1 -type numeric -extractor iline0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine0Extraction -set oline0iextraction -name points2 -type numeric -extractor iline0

# line 1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name u -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name v -type numeric -extractor iline1
if [ "$dimension" == "3" ]; then
	java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name w -type numeric -extractor iline1
fi
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name p -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name s -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name d -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name r -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name vtkvalidpointmask -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name arc_length -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name points0 -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name points1 -type numeric -extractor iline1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine1Extraction -set oline1iextraction -name points2 -type numeric -extractor iline1

# line 2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name u -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name v -type numeric -extractor iline2
if [ "$dimension" == "3" ]; then
	java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name w -type numeric -extractor iline2
fi
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name p -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name s -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name d -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name r -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name vtkvalidpointmask -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name arc_length -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name points0 -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name points1 -type numeric -extractor iline2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine2Extraction -set oline2iextraction -name points2 -type numeric -extractor iline2

# line 3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name u -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name v -type numeric -extractor iline3
if [ "$dimension" == "3" ]; then
	java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name w -type numeric -extractor iline3
fi
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name p -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name s -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name d -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name r -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name vtkvalidpointmask -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name arc_length -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name points0 -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name points1 -type numeric -extractor iline3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iLine3Extraction -set oline3iextraction -name points2 -type numeric -extractor iline3

echo "Initial Visualization"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag iVisualization
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation iVisualization -name libmesh-sedimentation-opt::iVisualization -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iVisualization -tag ogetmaximumiterationstotransport -type input -dependency getMaximumIterationsToTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation iVisualization -tag oivisualization -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iVisualization -set oivisualization -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iVisualization -set oivisualization -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation iVisualization -set oivisualization -name png -type file

echo "Solver Simulation to the Flow"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag solverSimulationFlow
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation solverSimulationFlow -name libmesh-sedimentation-opt::SolverSimulationFlow -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation solverSimulationFlow -tag ogetmaximumiterationstotransport -type input -dependency getMaximumIterationsToTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation solverSimulationFlow -tag osolversimulationflow -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFlow -set osolversimulationflow -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFlow -set osolversimulationflow -name t_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFlow -set osolversimulationflow -name dt -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFlow -set osolversimulationflow -name time -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFlow -set osolversimulationflow -name r -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFlow -set osolversimulationflow -name flow_l -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFlow -set osolversimulationflow -name flow_n_linear_iterations -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFlow -set osolversimulationflow -name flow_final_linear_residual -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFlow -set osolversimulationflow -name flow_norm_delta -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFlow -set osolversimulationflow -name flow_norm_delta_u -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationFlow -set osolversimulationflow -name flow_converged -type text

echo "Solver Simulation to the Transport"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag solverSimulationTransport
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation solverSimulationTransport -name libmesh-sedimentation-opt::SolverSimulationTransport -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation solverSimulationTransport -tag osolversimulationflow -type input -dependency solverSimulationFlow
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation solverSimulationTransport -tag osolversimulationtransport -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationTransport -set osolversimulationtransport -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationTransport -set osolversimulationtransport -name t_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationTransport -set osolversimulationtransport -name dt -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationTransport -set osolversimulationtransport -name time -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationTransport -set osolversimulationtransport -name r -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationTransport -set osolversimulationtransport -name transport_l -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationTransport -set osolversimulationtransport -name transport_n_linear_iterations -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationTransport -set osolversimulationtransport -name transport_final_linear_residual -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationTransport -set osolversimulationtransport -name transport_norm_delta -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationTransport -set osolversimulationtransport -name transport_norm_delta_u -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation solverSimulationTransport -set osolversimulationtransport -name transport_converged -type text

echo "Compute Solution Change"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag computeSolutionChange
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation computeSolutionChange -name libmesh-sedimentation-opt::EvaluateTimeStepControl -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation computeSolutionChange -tag osolversimulationtransport -type input -dependency solverSimulationTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation computeSolutionChange -tag ocomputeSolutionChange -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeSolutionChange -set ocomputeSolutionChange -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeSolutionChange -set ocomputeSolutionChange -name t_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeSolutionChange -set ocomputeSolutionChange -name time -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeSolutionChange -set ocomputeSolutionChange -name dt -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeSolutionChange -set ocomputeSolutionChange -name n_flow_linear_iterations_total -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeSolutionChange -set ocomputeSolutionChange -name n_flow_nonlinear_iterations_total -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeSolutionChange -set ocomputeSolutionChange -name n_transport_linear_iterations_total -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeSolutionChange -set ocomputeSolutionChange -name n_transport_nonlinear_iterations_total -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeSolutionChange -set ocomputeSolutionChange -name solution_converged -type text
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeSolutionChange -set ocomputeSolutionChange -name error -type numeric

echo "Compute Time Step"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag computeTimeStep
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation computeTimeStep -name libmesh-sedimentation-opt::EvaluateTimeStepControl -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation computeTimeStep -tag ocomputeSolutionChange -type input -dependency computeSolutionChange
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation computeTimeStep -tag ocomputeTimeStep -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeTimeStep -set ocomputeTimeStep -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeTimeStep -set ocomputeTimeStep -name t_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeTimeStep -set ocomputeTimeStep -name time -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeTimeStep -set ocomputeTimeStep -name dt -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation computeTimeStep -set ocomputeTimeStep -name ts_converged -type text

echo "Mesh Refinement"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag meshRefinement
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation meshRefinement -name libmesh-sedimentation-opt::MeshRefinement -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshRefinement -tag osolversimulationtransport -type input -dependency solverSimulationTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshRefinement -tag omeshrefinement -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name first_step_refinement -type text
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name t_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name before_n_active_elem -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshRefinement -set omeshrefinement -name after_n_active_elem -type numeric

echo "Mesh Writer"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag meshWriter
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation meshWriter -name libmesh-sedimentation-opt::MeshWriter -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshWriter -tag osolversimulationtransport -type input -dependency solverSimulationTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshWriter -tag omeshwriter -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshWriter -set omeshwriter -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshWriter -set omeshwriter -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshWriter -set omeshwriter -name xdmf -type file

# Analysis
echo "Horizontal Line Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag line0Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation line0Extraction -name libmesh-sedimentation-opt::Line0Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line0Extraction -tag osolversimulationtransport -type input -dependency solverSimulationTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line0Extraction -tag oline0extraction -type output

echo "Vertical Line 1 Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag line1Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation line1Extraction -name libmesh-sedimentation-opt::Line1Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line1Extraction -tag osolversimulationtransport -type input -dependency solverSimulationTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line1Extraction -tag oline1extraction -type output

echo "Vertical Line 2 Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag line2Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation line2Extraction -name libmesh-sedimentation-opt::Line2Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line2Extraction -tag osolversimulationtransport -type input -dependency solverSimulationTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line2Extraction -tag oline2extraction -type output

echo "Vertical Line 3 Extraction"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag line3Extraction
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation line3Extraction -name libmesh-sedimentation-opt::Line3Extraction -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line3Extraction -tag osolversimulationtransport -type input -dependency solverSimulationTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation line3Extraction -tag oline3extraction -type output

if [ "$access" == "extraction" ]; then
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line0Extraction -set oline0extraction -tag line0 -algorithm EXTRACTION:CSV
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line1Extraction -set oline1extraction -tag line1 -algorithm EXTRACTION:CSV
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line2Extraction -set oline2extraction -tag line2 -algorithm EXTRACTION:CSV
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line3Extraction -set oline3extraction -tag line3 -algorithm EXTRACTION:CSV
elif [ "$access" == "indexing" ] && [ "$cartridge" == "csv" ]; then
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line0Extraction -set oline0extraction -tag line0 -algorithm INDEXING:CSV
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line1Extraction -set oline1extraction -tag line1 -algorithm INDEXING:CSV
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line2Extraction -set oline2extraction -tag line2 -algorithm INDEXING:CSV
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line3Extraction -set oline3extraction -tag line3 -algorithm INDEXING:CSV
elif [ "$access" == "indexing" ] && [ "$cartridge" == "fastbit" ]; then
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line0Extraction -set oline0extraction -tag line0 -algorithm INDEXING:FASTBIT
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line1Extraction -set oline1extraction -tag line1 -algorithm INDEXING:FASTBIT
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line2Extraction -set oline2extraction -tag line2 -algorithm INDEXING:FASTBIT
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line3Extraction -set oline3extraction -tag line3 -algorithm INDEXING:FASTBIT
elif [ "$access" == "indexing" ] && [ "$cartridge" == "optimized_fastbit" ]; then
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line0Extraction -set oline0extraction -tag line0 -algorithm INDEXING:OPTIMIZED_FASTBIT
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line1Extraction -set oline1extraction -tag line1 -algorithm INDEXING:OPTIMIZED_FASTBIT
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line2Extraction -set oline2extraction -tag line2 -algorithm INDEXING:OPTIMIZED_FASTBIT
	java -jar ../dfa/PG-1.0.jar -extractor -dataflow sedimentation -transformation line3Extraction -set oline3extraction -tag line3 -algorithm INDEXING:OPTIMIZED_FASTBIT
fi

# line 0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name u -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name v -type numeric -extractor line0
if [ "$dimension" == "3" ]; then
	java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name w -type numeric -extractor line0
fi
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name p -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name s -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name d -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name r -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name vtkvalidpointmask -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name arc_length -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name points0 -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name points1 -type numeric -extractor line0
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line0Extraction -set oline0extraction -name points2 -type numeric -extractor line0

# line 1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name u -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name v -type numeric -extractor line1
if [ "$dimension" == "3" ]; then
	java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name w -type numeric -extractor line1
fi
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name p -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name s -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name d -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name r -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name vtkvalidpointmask -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name arc_length -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name points0 -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name points1 -type numeric -extractor line1
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line1Extraction -set oline1extraction -name points2 -type numeric -extractor line1

# line 2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name u -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name v -type numeric -extractor line2
if [ "$dimension" == "3" ]; then
	java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name w -type numeric -extractor line2
fi
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name p -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name s -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name d -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name r -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name vtkvalidpointmask -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name arc_length -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name points0 -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name points1 -type numeric -extractor line2
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line2Extraction -set oline2extraction -name points2 -type numeric -extractor line2

# line 3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name u -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name v -type numeric -extractor line3
if [ "$dimension" == "3" ]; then
	java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name w -type numeric -extractor line3
fi
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name p -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name s -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name d -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name r -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name vtkvalidpointmask -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name arc_length -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name points0 -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name points1 -type numeric -extractor line3
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation line3Extraction -set oline3extraction -name points2 -type numeric -extractor line3

echo "Visualization"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag visualization
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation visualization -name libmesh-sedimentation-opt::Visualization -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation visualization -tag osolversimulationtransport -type input -dependency solverSimulationTransport
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation visualization -tag ovisualization -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation visualization -set ovisualization -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation visualization -set ovisualization -name time_step -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation visualization -set ovisualization -name png -type file

echo "Mesh Aggregator"
java -jar ../dfa/PG-1.0.jar -transformation -dataflow sedimentation -tag meshAggregator
java -jar ../dfa/PG-1.0.jar -program -dataflow sedimentation -transformation meshAggregator -name libmesh-sedimentation-opt::MeshAggregator -filepath $PGDIR

java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshAggregator -tag omeshwriter -type input -dependency meshWriter
# java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshAggregator -tag ostatistics -type input -dependency computeStatistics
java -jar ../dfa/PG-1.0.jar -set -dataflow sedimentation -transformation meshAggregator -tag omeshaggregator -type output

java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshAggregator -set omeshaggregator -name simulationID -type numeric
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshAggregator -set omeshaggregator -name xdmf -type file
java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshAggregator -set omeshaggregator -name n_processors -type numeric
# java -jar ../dfa/PG-1.0.jar -attribute -dataflow sedimentation -transformation meshAggregator -set omeshaggregator -name video -type file

echo "Dataflow ingestion"
java -jar ../dfa/PG-1.0.jar -ingest -dataflow sedimentation

cp prov/pg/sedimentation/dataflow.json .


