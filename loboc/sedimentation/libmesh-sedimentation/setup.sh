#!/bin/bash
cd ${HOME}/Documents/Program_Installations/src

echo "Installing szip..."
tar xzvf ${HOME}/Documents/Program_Installations/src/szip-2.1.tar.gz
cd szip-2.1/
./configure --prefix=${HOME}/Documents/Program_Installations/szip
make
make install
cd ..

echo "Installing zlib..."
tar xzvf ${HOME}/Documents/Program_Installations/src/zlib-1.2.8.tar.gz
cd zlib-1.2.8/
./configure --prefix=${HOME}/Documents/Program_Installations/zlib
make
make install
cd ..

echo "Installing HDF5..."
tar xjvf ${HOME}/Documents/Program_Installations/src/hdf5-1.8.17.tar.bz2
cd hdf5-1.8.17/
./configure --prefix=${HOME}/Documents/Program_Installations/hdf5 --with-zlib=${HOME}/Documents/Program_Installations/zlib --with-szlib=${HOME}/Documents/Program_Installations/szip
make
make install
cd ..

echo "Installing cmake..."
#wget https://cmake.org/files/v3.6/cmake-3.6.0-rc4.tar.gz
tar xzvf ${HOME}/Documents/Program_Installations/src/cmake-3.6.0-rc4.tar.gz
cd cmake-3.6.0-rc4
./configure --prefix=${HOME}/Documents/Program_Installations/cmake
gmake
gmake install
cd ..

echo "Installing PETSc..."
#wget  http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.2.tar.gz
tar xzvf ${HOME}/Documents/Program_Installations/src/petsc-lite-3.7.2.tar.gz
cd petsc-3.7.2
export PETSC_DIR=$PWD
./configure --prefix=${HOME}/Documents/Program_Installations/petsc-3.7.2 --with-debugging=0 --download-f2cblaslapack=1 --with-fc=0
make PETSC_DIR=/Users/vitor/Documents/Program_Installations/src/petsc-3.7.2 PETSC_ARCH=arch-darwin-c-opt all
make PETSC_DIR=/Users/vitor/Documents/Program_Installations/src/petsc-3.7.2 PETSC_ARCH=arch-darwin-c-opt install
cd ..

export PETSC_DIR=/Users/vitor/Documents/Program_Installations/3.7.2

echo "Installing libMesh..."
#wget https://github.com/libMesh/libmesh/releases/download/v1.0.0/libmesh-1.0.0.tar.gz
git clone git://github.com/libMesh/libmesh.git 
tar xzvf ${HOME}/Documents/Program_Installations/src/libmesh-1.0.0.tar.gz
cd libmesh-1.0.0
./configure --prefix=${HOME}/Documents/Program_Installations/libmesh --enable-parmesh HDF5_DIR=${HOME}/Documents/Program_Installations/hdf5 LDFLAGS=-L${HOME}/Documents/Program_Installations/zlib/lib --disable-fortran
make -j 4
make install
cd ..

echo "Installing Mesa OS..."
#wget https://mesa.freedesktop.org/archive/current/mesa-11.0.7.tar.gz --no-check-certificate
tar -zxvf ${HOME}/Documents/Program_Installations/src/mesa-11.0.7.tar.gz
cd mesa-11.0.7
./configure CXXFLAGS="-O2 -g -DDEFAULT_SOFTWARE_DEPTH_BITS=31" CFLAGS="-O2 -g -DDEFAULT_SOFTWARE_DEPTH_BITS=31" --disable-xvmc --disable-glx --disable-dri --with-dri-drivers="" --with-gallium-drivers="swrast" --enable-texture-float --enable-gles2 --disable-egl --with-egl-platforms="" --enable-gallium-osmesa --enable-gallium-llvm=yes --prefix=${HOME}/Documents/Program_Installations/mesaos

./configure CXXFLAGS="-O2 -g -DDEFAULT_SOFTWARE_DEPTH_BITS=31"     CFLAGS="-O2 -g -DDEFAULT_SOFTWARE_DEPTH_BITS=31"     --disable-xvmc     --disable-glx     --disable-dri     --with-dri-drivers=""     --with-gallium-drivers="swrast"     --enable-texture-float --enable-gles2     --disable-egl     --with-egl-platforms=""     --enable-gallium-osmesa     --enable-gallium-llvm=yes   --prefix=${HOME}/Documents/Program_Installations/mesaos

make
make install
cd ..

echo "Installing ParaView..."
#wget http://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.1&type=source&os=all&downloadFile=ParaView-v5.1.2.tar.gz
tar -zxvf ${HOME}/Documents/Program_Installations/src/ParaView-v5.1.2.tar.gz
cd ParaView-v5.1.2
mkdir build
cd build
cmake ${HOME}/Documents/Program_Installations/src/ParaView-v5.1.2 -DMACOS_RPATH=${HOME}/Documents/Program_Installations/paraview -DBUILD_SHARED_LIBS=1 -DPARAVIEW_ENABLE_PYTHON=1 -DPARAVIEW_USE_MPI=1 -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_ENABLE_CATALYST=1 -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=1 -DPARAVIEW_BUILD_QT_GUI=0 -DVTK_OPENGL_HAS_OSMESA=1 -DVTK_USE_X=0 -DOSMESA_INCLUDE_DIR=/usr/local/Cellar/mesalib-glw/7.2/include -DOSMESA_LIBRARY=/usr/local/Cellar/mesalib-glw/7.2/lib/libOSMesa.so -DVTK_USE_OFFSCREEN=1 -DCMAKE_INSTALL_PREFIX=${HOME}/Documents/Program_Installations/paraview -DOPENGL_INCLUDE_DIR= -DOPENGL_gl_LIBRARY= -DOPENGL_glu_LIBRARY= 

make -j 4
make install
cd ..


