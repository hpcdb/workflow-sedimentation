# Introduction

This repository contains the source code of **libMesh-sedimentation**. It is a simulation solver built on top of a widely used parallel fine element framework, [libMesh](http://libmesh.github.io), which supports parallel simulation of multiscale, multiphysics applications. libMesh interfaces with several libraries for Computational Science and Engineering applications (e.g., PeTSc, Metis, Parmetis, LAPACK). Also, scientific visualization tools like [ParaView](https://www.paraview.org/), are typically used in these applications to gain insight from the computations. 

# User documentation for docker

## Clones
- docker pull vitorss/debian:sedimentation-vitor
- git clone git@github.com:hpcdb/workflow-sedimentation.git
- LOCAL_DIR=/path/to/repository/in/local/machine

## Run docker:
```
docker run -v $LOCAL_DIR:/shared -i -t vitorss/debian:sedimentation-vitor /bin/bash
```
Now, in docker:

## Adjust the Docker variables

### /shared/src-local/sedimentation/DfA.properties
```
# Docker
pg_dir=/shared/workflow-sedimentation/src-local/sedimentation/prov/pg
di_dir=/shared/workflow-sedimentation/src-local/sedimentation/prov/di
```

### /shared/src-local/sedimentation/execute-solver.sh
```
# Docker
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/opt/paraview/5.4/lib/paraview-5.4
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/paraview/5.4/lib/paraview-5.4
```

Change the end to:
```
time mpirun -np 2 ../../libmesh-sedimentation/sediment-opt -i lock_necker3D_pc11.in -m lock_necker3D.msh -e extraction.py -v visualization.py -o output -d /shared/workflow-sedimentation/src-local/sedimentation/output | tee -a "output-solver.log"
```

### /shared/workflow-sedimentation/src-local/sedimentation/pg-dataflow.sh
```
# Docker
PGDIR=/shared/workflow-sedimentation/src-local/libmesh-sedimentation
```

### /shared/workflow-sedimentation/src-local/sedimentation/provenance.in
```
# Docker
directory = /shared/workflow-sedimentation/src-local/sedimentation
pgFilePath = /shared/workflow-sedimentation/src-local/dfa/PG-1.0.jar
rdeFilePath = /shared/workflow-sedimentation/src-local/dfa/RDE-1.0.jar
rdiFilePath = /shared/workflow-sedimentation/src-local/dfa/RDI-1.0.jar
access=EXTRACTION
cartridge=csv
bin=/programs/fastbit/bin
extraArguments=[b:precision=2]
```

## Make the solver:

Uncomment `include config/config.docker`
```
export LIBMESH_DIR=/programs/libmesh/build
cd /shared/libmesh-sedimentation
rm -rf .depend
rm -rf .libs
make clean
make
```

## Cheat sheet of Execution commands:

```
monetdb-start-all
cd /shared/workflow-sedimentation/src-local/sedimentation
./pg-dataflow.sh --> generates the conceptual dataflow
./delete.sh --> deletes data from a previous execution
./start-data-ingestor.sh --> cleans the database and initializes MonetDB
./execute-solver.sh --> Executes the sedimentation solver
```
