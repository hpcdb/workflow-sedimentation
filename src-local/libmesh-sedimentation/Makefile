######################################################################
#
# Template libMesh application Makefile

# MacOS
#LIBMESH_DIR ?=/Users/vitor/Documents/Program_Installations/libmesh
# docker
# LIBMESH_DIR ?=/programs/libmesh/build
# Virtual Box
# LIBMESH_DIR ?=/home/vitor/src/libmesh-install
# Stampede
#LIBMESH_DIR ?=/work/03664/silva/programs/libmesh/build
# LoboC
#LIBMESH_DIR ?=/home/user4/local/intel/libmesh
#LIBMESH_DIR ?=/home/user2/local/gnu/libmesh-git
# LIBMESH_DIR ?=/home/user4/local/gnu/libmesh
# Mint
LIBMESH_DIR ?=/home/vitor/Documents/program/libmesh

# METHOD = dbg
# include the library options determined by configure
include $(LIBMESH_DIR)/Make.common

target     := ./libmesh-sedimentation-$(METHOD)


###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
srcfiles 	:= $(wildcard *.C dfanalyzer/*.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################

# MacOS
# PARAVIEW_DIR=/Users/vitor/Documents/Program_Installations/ParaView-v5.2.0/build/CMakeFiles/__macos_install
# docker
# PARAVIEW_DIR=/programs/paraview-install
# Stampede
#PARAVIEW_DIR=/work/03664/silva/programs/paraview
# LoboC
#PARAVIEW_DIR=/home/user4/local/all/paraview
#PARAVIEW_DIR=/home/user2/local/gnu/paraview/5.1
# PARAVIEW_DIR=/home/user4/local/all/paraview
# Mint
PARAVIEW_DIR=/home/vitor/Documents/program/paraview

# PARAVIEW_LIB_DIR=$(PARAVIEW_DIR)/lib/paraview-5.0
# PARAVIEW_LIBS =-L$(PARAVIEW_LIB_DIR) -lvtkPVPythonCatalyst-pv5.0 -lvtkPVCatalyst-pv5.0 -lvtkPVServerManagerApplication-pv5.0 -lvtkPVServerManagerCore-pv5.0 -lvtkPVServerImplementationCore-pv5.0 -lvtkPVClientServerCoreCore-pv5.0 -lvtkFiltersParallel-pv5.0 -lvtkFiltersModeling-pv5.0 -lvtkRenderingCore-pv5.0 -lvtkFiltersExtraction-pv5.0 -lvtkFiltersStatistics-pv5.0 -lvtkImagingFourier-pv5.0 -lvtkImagingCore-pv5.0 -lvtkalglib-pv5.0 -lvtkFiltersGeometry-pv5.0 -lvtkFiltersSources-pv5.0 -lvtkFiltersGeneral-pv5.0 -lvtkCommonComputationalGeometry-pv5.0 -lvtkCommonColor-pv5.0 -lvtkFiltersProgrammable-pv5.0 -lvtkPVVTKExtensionsCore-pv5.0 -lvtkPVCommon-pv5.0 -lvtkClientServer-pv5.0 -lvtkIOXMLParser-pv5.0 -lvtkFiltersCore-pv5.0 -lvtkParallelMPI-pv5.0 -lvtkParallelCore-pv5.0 -lvtkIOLegacy-pv5.0 -lvtkCommonExecutionModel-pv5.0 -lvtkCommonDataModel-pv5.0 -lvtkCommonSystem-pv5.0 -lvtkCommonTransforms-pv5.0 -lvtkCommonMisc-pv5.0 -lvtkCommonMath-pv5.0 -lprotobuf -lpthread  -lvtkPythonInterpreter-pv5.0 -lvtkWrappingPython27Core-pv5.0 -lvtkCommonCore-pv5.0 -lpython2.6 -lvtksys-pv5.0 -Wl,-rpath,$(PARAVIEW_LIB_DIR) 
# PARAVIEW_INCLUDE=-I$(PARAVIEW_DIR)/include/paraview-5.0
# libmesh_CXXFLAGS += $(PARAVIEW_INCLUDE) -I../../xdmf
# libmesh_LDFLAGS+=$(PARAVIEW_LIBS)

# PARAVIEW_LIB_DIR=$(PARAVIEW_DIR)/lib/paraview-5.0
# PARAVIEW_LIBS =-L$(PARAVIEW_LIB_DIR) -lvtkPVPythonCatalyst-pv5.0 -lvtkPVCatalyst-pv5.0 -lvtkPVServerManagerApplication-pv5.0 -lvtkPVServerManagerCore-pv5.0 -lvtkPVServerImplementationCore-pv5.0 -lvtkPVClientServerCoreCore-pv5.0 -lvtkFiltersParallel-pv5.0 -lvtkFiltersModeling-pv5.0 -lvtkRenderingCore-pv5.0 -lvtkFiltersExtraction-pv5.0 -lvtkFiltersStatistics-pv5.0 -lvtkImagingFourier-pv5.0 -lvtkImagingCore-pv5.0 -lvtkalglib-pv5.0 -lvtkFiltersGeometry-pv5.0 -lvtkFiltersSources-pv5.0 -lvtkFiltersGeneral-pv5.0 -lvtkCommonComputationalGeometry-pv5.0 -lvtkCommonColor-pv5.0 -lvtkFiltersProgrammable-pv5.0 -lvtkPVVTKExtensionsCore-pv5.0 -lvtkPVCommon-pv5.0 -lvtkClientServer-pv5.0 -lvtkIOXMLParser-pv5.0 -lvtkFiltersCore-pv5.0 -lvtkParallelMPI-pv5.0 -lvtkParallelCore-pv5.0 -lvtkIOLegacy-pv5.0 -lvtkCommonExecutionModel-pv5.0 -lvtkCommonDataModel-pv5.0 -lvtkCommonSystem-pv5.0 -lvtkCommonTransforms-pv5.0 -lvtkCommonMisc-pv5.0 -lvtkCommonMath-pv5.0 -lprotobuf -lpthread  -lvtkPythonInterpreter-pv5.0 -lvtkWrappingPython27Core-pv5.0 -lvtkCommonCore-pv5.0 -lpython2.7 -lvtksys-pv5.0 -Wl,-rpath,$(PARAVIEW_LIB_DIR) 
# PARAVIEW_INCLUDE=-I$(PARAVIEW_DIR)/include/paraview-5.0
# libmesh_CXXFLAGS += $(PARAVIEW_INCLUDE) -I../../xdmf
# libmesh_LDFLAGS+=$(PARAVIEW_LIBS)

# Virtual Box
# PARAVIEW_DIR=/home/vitor/src/paraview

# PARAVIEW_LIB_DIR=$(PARAVIEW_DIR)/lib/paraview-5.1
# PARAVIEW_LIBS =-L$(PARAVIEW_LIB_DIR) -lvtkPVPythonCatalyst-pv5.1 -lvtkPVCatalyst-pv5.1 -lvtkPVServerManagerApplication-pv5.1 -lvtkPVServerManagerCore-pv5.1 -lvtkPVServerImplementationCore-pv5.1 -lvtkPVClientServerCoreCore-pv5.1 -lvtkFiltersParallel-pv5.1 -lvtkFiltersModeling-pv5.1 -lvtkRenderingCore-pv5.1 -lvtkFiltersExtraction-pv5.1 -lvtkFiltersStatistics-pv5.1 -lvtkImagingFourier-pv5.1 -lvtkImagingCore-pv5.1 -lvtkalglib-pv5.1 -lvtkFiltersGeometry-pv5.1 -lvtkFiltersSources-pv5.1 -lvtkFiltersGeneral-pv5.1 -lvtkCommonComputationalGeometry-pv5.1 -lvtkCommonColor-pv5.1 -lvtkFiltersProgrammable-pv5.1 -lvtkPVVTKExtensionsCore-pv5.1 -lvtkPVCommon-pv5.1 -lvtkClientServer-pv5.1 -lvtkIOXMLParser-pv5.1 -lvtkFiltersCore-pv5.1 -lvtkParallelMPI-pv5.1 -lvtkParallelCore-pv5.1 -lvtkIOLegacy-pv5.1 -lvtkCommonExecutionModel-pv5.1 -lvtkCommonDataModel-pv5.1 -lvtkCommonSystem-pv5.1 -lvtkCommonTransforms-pv5.1 -lvtkCommonMisc-pv5.1 -lvtkCommonMath-pv5.1 -lprotobuf -lpthread  -lvtkPythonInterpreter-pv5.1 -lvtkWrappingPython27Core-pv5.1 -lvtkCommonCore-pv5.1 -lpython2.6 -lvtksys-pv5.1 -Wl,-rpath,$(PARAVIEW_LIB_DIR) 
# PARAVIEW_INCLUDE=-I$(PARAVIEW_DIR)/include/paraview-5.1
# libmesh_CXXFLAGS += $(PARAVIEW_INCLUDE) -I../../xdmf
# libmesh_LDFLAGS+=$(PARAVIEW_LIBS)

PARAVIEW_LIB_DIR=$(PARAVIEW_DIR)/lib/paraview-5.1
PARAVIEW_LIBS =-L$(PARAVIEW_LIB_DIR) -lvtkPVPythonCatalyst-pv5.1 -lvtkPVCatalyst-pv5.1 -lvtkPVServerManagerApplication-pv5.1 -lvtkPVServerManagerCore-pv5.1 -lvtkPVServerImplementationCore-pv5.1 -lvtkPVClientServerCoreCore-pv5.1 -lvtkFiltersParallel-pv5.1 -lvtkFiltersModeling-pv5.1 -lvtkRenderingCore-pv5.1 -lvtkFiltersExtraction-pv5.1 -lvtkFiltersStatistics-pv5.1 -lvtkImagingFourier-pv5.1 -lvtkImagingCore-pv5.1 -lvtkalglib-pv5.1 -lvtkFiltersGeometry-pv5.1 -lvtkFiltersSources-pv5.1 -lvtkFiltersGeneral-pv5.1 -lvtkCommonComputationalGeometry-pv5.1 -lvtkCommonColor-pv5.1 -lvtkFiltersProgrammable-pv5.1 -lvtkPVVTKExtensionsCore-pv5.1 -lvtkPVCommon-pv5.1 -lvtkClientServer-pv5.1 -lvtkIOXMLParser-pv5.1 -lvtkFiltersCore-pv5.1 -lvtkParallelMPI-pv5.1 -lvtkParallelCore-pv5.1 -lvtkIOLegacy-pv5.1 -lvtkCommonExecutionModel-pv5.1 -lvtkCommonDataModel-pv5.1 -lvtkCommonSystem-pv5.1 -lvtkCommonTransforms-pv5.1 -lvtkCommonMisc-pv5.1 -lvtkCommonMath-pv5.1 -lprotobuf -lpthread  -lvtkPythonInterpreter-pv5.1 -lvtkWrappingPython27Core-pv5.1 -lvtkCommonCore-pv5.1 -lpython2.7 -lvtksys-pv5.1 -Wl,-rpath,$(PARAVIEW_LIB_DIR) 
PARAVIEW_INCLUDE=-I$(PARAVIEW_DIR)/include/paraview-5.1
libmesh_CXXFLAGS += $(PARAVIEW_INCLUDE) -I../../xdmf
libmesh_LDFLAGS+=$(PARAVIEW_LIBS)

#PARAVIEW_LIB_DIR=$(PARAVIEW_DIR)/lib/paraview-5.2
#PARAVIEW_LIBS =-L$(PARAVIEW_LIB_DIR) -lvtkPVPythonCatalyst-pv5.2 -lvtkPVCatalyst-pv5.2 -lvtkPVServerManagerApplication-pv5.2 -lvtkPVServerManagerCore-pv5.2 -lvtkPVServerImplementationCore-pv5.2 -lvtkPVClientServerCoreCore-pv5.2 -lvtkFiltersParallel-pv5.2 -lvtkFiltersModeling-pv5.2 -lvtkRenderingCore-pv5.2 -lvtkFiltersExtraction-pv5.2 -lvtkFiltersStatistics-pv5.2 -lvtkImagingFourier-pv5.2 -lvtkImagingCore-pv5.2 -lvtkalglib-pv5.2 -lvtkFiltersGeometry-pv5.2 -lvtkFiltersSources-pv5.2 -lvtkFiltersGeneral-pv5.2 -lvtkCommonComputationalGeometry-pv5.2 -lvtkCommonColor-pv5.2 -lvtkFiltersProgrammable-pv5.2 -lvtkPVVTKExtensionsCore-pv5.2 -lvtkPVCommon-pv5.2 -lvtkClientServer-pv5.2 -lvtkIOXMLParser-pv5.2 -lvtkFiltersCore-pv5.2 -lvtkParallelMPI-pv5.2 -lvtkParallelCore-pv5.2 -lvtkIOLegacy-pv5.2 -lvtkCommonExecutionModel-pv5.2 -lvtkCommonDataModel-pv5.2 -lvtkCommonSystem-pv5.2 -lvtkCommonTransforms-pv5.2 -lvtkCommonMisc-pv5.2 -lvtkCommonMath-pv5.2 -lprotobuf -lpthread  -lvtkPythonInterpreter-pv5.2 -lvtkWrappingPython27Core-pv5.2 -lvtkCommonCore-pv5.2 -lpython2.7 -lvtksys-pv5.2 -Wl,-rpath,$(PARAVIEW_LIB_DIR) 
#PARAVIEW_INCLUDE=-I$(PARAVIEW_DIR)/include/paraview-5.2
#libmesh_CXXFLAGS += $(PARAVIEW_INCLUDE) -I../../xdmf
#libmesh_LDFLAGS+=$(PARAVIEW_LIBS)

.PHONY: dust clean distclean

###############################################################################
# Target:
#

all:: $(notdir $(target))

# Production rules:  how to make the target - depends on library configuration
$(notdir $(target)): $(objects)
	@echo "Linking "$@"..."
	@$(libmesh_LIBTOOL) --tag=CXX $(LIBTOOLFLAGS) --mode=link \
	  $(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS)


# Useful rules.
dust:
	@echo "Deleting old output and runtime files"
	@rm -f out*.m job_output.txt output.txt* *.gmv.* *.plt.* out*.xdr* out*.xda* PI* complete

clean: dust
	@rm -f $(objects) *.$(obj-suffix) *.lo

clobber: clean 
	@rm -f $(target)

distclean: clean
	@rm -rf *.o .libs .depend

echo:
	@echo srcfiles = $(srcfiles)
	@echo objects = $(objects)
	@echo target = $(target)

run: complete

complete: $(wildcard *.in)
#	@$(MAKE) dust
	@$(MAKE) -C $(dir $(target)) $(notdir $(target))
	@echo "***************************************************************"
	@echo "* Running App " $(notdir $(target))
	@echo "***************************************************************"
	@echo " "
	${LIBMESH_RUN} $(target) ${LIBMESH_OPTIONS} 2>&1 | tee output.txt
	@bzip2 -f output.txt
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running App " $(notdir $(target))
	@echo "***************************************************************"

gmv:
	@$(MAKE) -C $(LIBMESH_DIR)/roy/meshplot/ meshplot-$(METHOD)
	@for file in out.mesh.*; do ${LIBMESH_RUN} $(LIBMESH_DIR)/roy/meshplot/meshplot-$(METHOD) $$file out.soln.$${file##out.mesh.} out.gmv.$${file:9:4}; done

# include the dependency list
-include .depend

#
# Dependencies
#
.depend: $(srcfiles) $(LIBMESH_DIR)/include/libmesh/*.h
	@$(perl) $(LIBMESH_DIR)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(LIBMESH_DIR)/include $(wildcard $(LIBMESH_DIR)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend

###############################################################################
