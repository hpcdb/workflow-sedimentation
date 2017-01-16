#!/bin/bash
#mkdir sources
cd sources

echo "***************************************************"
echo " Installing zlib..."
echo "***************************************************"
#wget http://zlib.net/zlib-1.2.8.tar.gz
tar xzvf zlib-1.2.8.tar.gz
cd zlib-1.2.8/
./configure --prefix=${HOME}/local/intel/zlib
make CC=icc
make install
cd ..

echo "***************************************************"
echo " Installing cmake..."
echo "***************************************************"
#wget https://cmake.org/files/v3.6/cmake-3.6.0-rc4.tar.gz
tar xzvf cmake-3.6.0-rc4.tar.gz
cd cmake-3.6.0-rc4
./configure --prefix=${HOME}/local/all/cmake
gmake
gmake install
cd ..

echo "***************************************************"
echo " Installing szip..."
echo "***************************************************"
#wget http://www.hdfgroup.org/ftp/lib-external/szip/2.1/src/szip-2.1.tar.gz
tar xzvf szip-2.1.tar.gz
cd szip-2.1/
./configure --prefix=${HOME}/local/intel/szip CC=icc CXX=icpc FC=ifort
make
make install
cd ..

echo "***************************************************"
echo " Installing HDF5..."
echo "***************************************************"
#wget https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.18.tar
tar xjvf hdf5-1.8.18.tar.bz2
cd hdf5-1.8.18/
./configure --prefix=${HOME}/local/intel/hdf5 --with-zlib=${HOME}/local/intel/zlib --with-szlib=${HOME}/local/intel/szip CC=icc CXX=icpc FC=ifort --enable-fortran
make -j 4
make install
cd ..

echo "***************************************************"
echo " Installing PETSc..."
echo "***************************************************"
#wget  http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.2.tar.gz
tar xzvf petsc-lite-3.7.2.tar.gz
cd petsc-3.7.2
export PETSC_DIR=$PWD
export PETSC_ARCH=linux-intel
./configure --prefix=${HOME}/local/intel/petsc_lite --with-debugging=0 --with-blas-lapack-dir=/sw/apps/intel16/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64 --with-cc=mpiicc --with-cxx=mpiicpc --with-fc=mpiifort 
make PETSC_DIR=/home/user4/sources/petsc-3.7.2 PETSC_ARCH=arch-linux2-c-opt all
make PETSC_DIR=/home/user4/sources/petsc-3.7.2 PETSC_ARCH=arch-linux2-c-opt install
cd ..

export PETSC_DIR=${HOME}/local/intel/petsc_lite

echo "***************************************************"
echo " Installing libMesh..."
echo "***************************************************"
#wget https://github.com/libMesh/libmesh/releases/download/v1.0.0/libmesh-1.0.0.tar.gz
tar xzvf libmesh-1.0.0.tar.gz
cd libmesh-1.0.0
export I_MPI_F77=ifort
./configure --prefix=${HOME}/local/intel/libmesh --enable-parmesh HDF5_DIR=${HOME}/local/intel/hdf5 LDFLAGS=-L${HOME}/local/intel/zlib/lib CC=mpiicc CXX=mpiicpc FC=mpiifort
make -j 12
make install
cd ..

echo "***************************************************"
echo "Installing Mesa OS..."
echo "***************************************************"
#wget https://mesa.freedesktop.org/archive/current/mesa-11.0.7.tar.gz --no-check-certificate
tar -zxvf mesa-11.0.7.tar.gz
cd mesa-11.0.7
./configure CXXFLAGS="-O2 -g -DDEFAULT_SOFTWARE_DEPTH_BITS=31" CFLAGS="-O2 -g -DDEFAULT_SOFTWARE_DEPTH_BITS=31" --disable-xvmc --disable-glx --disable-dri --with-dri-drivers="" --with-gallium-drivers="swrast" --enable-texture-float --enable-gles2 --disable-egl --with-egl-platforms="" --enable-gallium-osmesa --enable-gallium-llvm=yes --prefix=${HOME}/local/all/mesaos
make
make install
cd ..

#echo "***************************************************"
#echo "Installing Python..."
#echo "***************************************************"
#wget https://www.python.org/ftp/python/3.5.2/Python-3.5.2.tgz --no-check-certificate
#wget https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz
#tar zxvf Python-3.5.2.tgz
#cd Python-3.5.2
#./configure --prefix=${HOME}/local/all/python
#make
#make install
#cd ..

echo "***************************************************"
echo "Installing ParaView..."
echo "***************************************************"
#wget http://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.1&type=source&os=all&downloadFile=ParaView-v5.1.2.tar.gz
cp -rf ~/ParaView-v5.1.2.tar.gz ~/sources/ParaView-v5.1.2.tar.gz
tar -zxvf ParaView-v5.1.2.tar.gz
cd ParaView-v5.1.2
mkdir build
cd build
~/local/all/cmake/bin/cmake ${HOME}/sources/ParaView-v5.1.2 -DBUILD_SHARED_LIBS=1 -DPARAVIEW_ENABLE_PYTHON=1 -DPARAVIEW_USE_MPI=1 -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_ENABLE_CATALYST=1 -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=1 -DPARAVIEW_BUILD_QT_GUI=0 -DVTK_OPENGL_HAS_OSMESA=1 -DVTK_USE_X=0 -DOSMESA_INCLUDE_DIR=${HOME}/local/all/mesaos/include -DOSMESA_LIBRARY=${HOME}/local/all/mesaos/lib/libOSMesa.so -DVTK_USE_OFFSCREEN=1 -DCMAKE_INSTALL_PREFIX=${HOME}/local/all/paraview -DOPENGL_INCLUDE_DIR= -DOPENGL_gl_LIBRARY= -DOPENGL_glu_LIBRARY= 
#/home/vitor/program/cmake/bin/cmake ~/Downloads/ParaView-v5.2.0 -DBUILD_SHARED_LIBS=1 -DPARAVIEW_ENABLE_PYTHON=1 -DPARAVIEW_USE_MPI=1 -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_ENABLE_CATALYST=1 -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=1 -DPARAVIEW_BUILD_QT_GUI=0 -DVTK_OPENGL_HAS_OSMESA=1 -DVTK_USE_X=0 -DOSMESA_INCLUDE_DIR=${HOME}/program/mesaos/include -DOSMESA_LIBRARY=$HOME/program/mesaos/lib/libOSMesa.so -DVTK_USE_OFFSCREEN=1 -DCMAKE_INSTALL_PREFIX=/home/vitor/program/paraview -DOPENGL_INCLUDE_DIR= -DOPENGL_gl_LIBRARY= -DOPENGL_glu_LIBRARY=
make -j 4
make install
cd ..