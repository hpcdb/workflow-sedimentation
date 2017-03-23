# MacOS

# PETSc - Installation
./configure --prefix=/Users/vitor/Documents/Program_Installations/petsc-3.6.4/build --download-fblaslapack
make PETSC_DIR=/Users/vitor/Documents/Program_Installations/petsc-3.6.4 PETSC_ARCH=arch-darwin-c-debug all
make PETSC_DIR=/Users/vitor/Documents/Program_Installations/petsc-3.6.4 PETSC_ARCH=arch-darwin-c-debug install
export PETSC_DIR=/Users/vitor/Documents/Program_Installations/petsc-3.6.4

# libMesh - Installation
./configure --enable-parmesh --prefix=/Users/vitor/Documents/Program_Installations/libmesh/build --with-hdf5=/Users/vitor/Documents/Program_Installations/hdf5-1.8.9/hdf5 PETSC_DIR=/Users/vitor/Documents/Program_Installations/petsc-3.6.4
make -j 2
make install

# OSMesa Gallium llvmpipe state-tracker
 ./configure     CXXFLAGS="-O2 -g -DDEFAULT_SOFTWARE_DEPTH_BITS=31"     CFLAGS="-O2 -g -DDEFAULT_SOFTWARE_DEPTH_BITS=31"     --disable-xvmc     --disable-glx     --disable-dri     --with-dri-drivers=""     --with-gallium-drivers="swrast"     --enable-texture-float --enable-gles2     --disable-egl     --with-egl-platforms=""     --enable-gallium-osmesa     --enable-gallium-llvm=yes   --prefix=/opt/libs/mesa
make 
make install

# Paraview
BUILD_SHARED_LIBS: ON
CMAKE_BUILD_TYPE: Release
PARAVIEW_BUILD_QT_GUI: OFF
PARAVIEW_ENABLE_CATALYST: ON
PARAVIEW_ENABLE_PYTHON: ON
PARAVIEW_INSTALL_DEVELOPMENT_FILES:0N
PARAVIEW_USE_MPI: ON

VTK_OPENGL_HAS_OSMESA: ON
VTK_USE_X: OFF
VTK_USE_OFFSCREEN=ON

OSMESA_INCLUDE_DIR and OPENGL_INCLUDE_DIR: Make sure these are not set to the OpenGL directory location of the header files.
OSMESA_LIBRARY: Set this to the Mesa library.
OPENGL_gl_LIBRARY and OPENGL_glu_LIBRARY: These should not be set.


========================================================================
# Docker
# PETSc
./configure --prefix=/programs/petsc-3.7.1/build --download-fblaslapack
make PETSC_DIR=/programs/petsc-3.7.1 all
make PETSC_DIR=/programs/petsc-3.7.1 install
export PETSC_DIR=/programs/petsc-3.7.1

# libMesh
./configure --enable-parmesh --prefix=/programs/libmesh/build --with-hdf5=/programs/hdf5-1.8.17/build PETSC_DIR=/programs/petsc-3.7.1
make -j 2
make install

# MesaOS
 ./configure     CXXFLAGS="-O2 -g -DDEFAULT_SOFTWARE_DEPTH_BITS=31"     CFLAGS="-O2 -g -DDEFAULT_SOFTWARE_DEPTH_BITS=31"     --disable-xvmc     --disable-glx     --disable-dri     --with-dri-drivers=""     --with-gallium-drivers="swrast"     --enable-texture-float --enable-gles2     --disable-egl     --with-egl-platforms=""     --enable-gallium-osmesa     --enable-gallium-llvm=yes   --prefix=/opt/libs/mesa
make
make install

# ParaView
apt-get install libphonon-dev libphonon4 qt4-dev-tools libqt4-core libqt4-gui qt4-qmake libxt-dev g++ gcc cmake-curses-gui libqt4-opengl-dev mesa-common-dev openmpi-common openmpi-bin libopenmpi-dev python-dev openmpi-common openmpi-bin libopenmpi-dev
cmake /programs/ParaView-5.3.0-RC2-Qt5-OpenGL2-MPI-Linux-64bit -DBUILD_SHARED_LIBS=1 -DPARAVIEW_ENABLE_PYTHON=1 -DPARAVIEW_USE_MPI=1 -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_ENABLE_CATALYST=1 -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=1 -DPARAVIEW_BUILD_QT_GUI=0 -DVTK_OPENGL_HAS_OSMESA=1 -DVTK_USE_X=0 -DOSMESA_INCLUDE_DIR=/programs/mesa-build/include -DOSMESA_LIBRARY=/programs/mesa-build/lib/libOSMesa.so -DVTK_USE_OFFSCREEN=1 -DCMAKE_INSTALL_PREFIX=/programs/paraview-5.3.0 -DMPI_C_LIBRARIES=/usr/lib/libmpi.so -DMPI_C_INCLUDE_PATH=/usr/include/mpi -DOPENGL_INCLUDE_DIR= -DOPENGL_gl_LIBRARY= -DOPENGL_glu_LIBRARY= -DMPI_C_LIBRARIES=/usr/lib/libmpi.so -DMPI_C_INCLUDE_PATH=/usr/include/mpi
make
make install

==========================================================================
# Stampede

# PETSc - Installation
./configure --prefix=/work/03664/silva/programs/petsc-3.6.4/build --download-fblaslapack
make PETSC_DIR=/work/03664/silva/programs/petsc-3.6.4 PETSC_ARCH=arch-darwin-c-debug all
make PETSC_DIR=/work/03664/silva/programs/petsc-3.6.4 PETSC_ARCH=arch-darwin-c-debug install

export PETSC_DIR=/work/03664/silva/programs/petsc-3.6.4


# libMesh - Installation
./configure --enable-parmesh --prefix=/work/03664/silva/programs/libmesh/build --with-hdf5=/work/03664/silva/programs/hdf5-1.8.9/hdf5 PETSC_DIR=/work/03664/silva/programs/petsc-3.7.2/build
make -j 2
make install

========================================================================
# Virtual Box - Mint

# SZIP
./configure --prefix=/home/vitor/programs/szip
make
make install

# ZLIB
./configure --prefix=/home/vitor/programs/zlib
make
make install

# HDF5
./configure --prefix=/home/vitor/programs/hdf5 --enable-fortran --enable-cxx --with-szlib=/home/vitor/programs/szip
make
make install

# PETSc
./configure --prefix=/home/vitor/programs/petsc --download-fblaslapack 
make PETSC_DIR=/home/vitor/programs/petsc all
make PETSC_DIR=/home/vitor/programs/petsc install
export PETSC_DIR=/home/vitor/programs/petsc

# libMesh
./configure --enable-parmesh --prefix=/home/vitor/programs/libmesh-install --with-hdf5=/home/vitor/programs/hdf5 PETSC_DIR=/home/vitor/programs/petsc
make -j 2
make install

# MesaOS
 ./configure     CXXFLAGS="-O2 -g -DDEFAULT_SOFTWARE_DEPTH_BITS=31"     CFLAGS="-O2 -g -DDEFAULT_SOFTWARE_DEPTH_BITS=31"     --disable-xvmc     --disable-glx     --disable-dri     --with-dri-drivers=""     --with-gallium-drivers="swrast"     --enable-texture-float --enable-gles2     --disable-egl     --with-egl-platforms=""     --enable-gallium-osmesa     --enable-gallium-llvm=yes   --prefix=/opt/libs/mesa
make
make install

# ParaView
apt-get install libphonon-dev libphonon4 qt4-dev-tools libqt4-core libqt4-gui qt4-qmake libxt-dev g++ gcc cmake-curses-gui libqt4-opengl-dev mesa-common-dev openmpi-common openmpi-bin libopenmpi-dev python-dev openmpi-common openmpi-bin libopenmpi-dev
cmake /programs/ParaView-v5.3.0-RC2 -DBUILD_SHARED_LIBS=1 -DPARAVIEW_ENABLE_PYTHON=1 -DPARAVIEW_USE_MPI=1 -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_ENABLE_CATALYST=1 -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=1 -DPARAVIEW_BUILD_QT_GUI=0 -DVTK_OPENGL_HAS_OSMESA=1 -DVTK_USE_X=0 -DOSMESA_INCLUDE_DIR=/programs/mesa-build/include -DOSMESA_LIBRARY=/programs/mesa-build/lib/libOSMesa.so -DVTK_USE_OFFSCREEN=1 -DCMAKE_INSTALL_PREFIX=/programs/paraview-5.3.0 -DMPI_C_LIBRARIES=/usr/lib/libmpi.so -DMPI_C_INCLUDE_PATH=/usr/include/mpi -DOPENGL_INCLUDE_DIR= -DOPENGL_gl_LIBRARY= -DOPENGL_glu_LIBRARY= -DMPI_C_LIBRARIES=/usr/lib/libmpi.so -DMPI_C_INCLUDE_PATH=/usr/include/mpi
make
make install

--------------------------------------------

SETUP PYTHON ENVIRONMENT - LoboC
export PYTHONPATH=/home/user2/sources/ParaView-v5.1.0/build/lib:/home/user2/sources/ParaView-v5.1.0/build/lib/site-packages
export LD_LIBRARY_PATH=/home/user2/sources/ParaView-v5.1.0/build/lib:$LD_LIBRARY_PATH
export PATH=$PATH:/home/user2/sources/ParaView-v5.1.0/build/bin