######################################################################
#
# Template libMesh application Makefile
LIBMESH_DIR ?=/Users/vitor/Documents/Program_Installations/libmesh-git/build
# METHOD = dbg
# include the library options determined by configure
include $(LIBMESH_DIR)/Make.common

target     := ./libmesh-sedimentation-$(METHOD)


###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################

PARAVIEW_DIR=/Users/vitor/Documents/Program_Installations/ParaView-v5.0.1-source/build/CMakeFiles/__macos_install
PARAVIEW_LIB_DIR=$(PARAVIEW_DIR)/lib/paraview-5.0
PARAVIEW_LIBS =-L$(PARAVIEW_LIB_DIR) -lvtkPVPythonCatalyst-pv5.0 -lvtkPVCatalyst-pv5.0 -lvtkPVServerManagerApplication-pv5.0 -lvtkPVServerManagerCore-pv5.0 -lvtkPVServerImplementationCore-pv5.0 -lvtkPVClientServerCoreCore-pv5.0 -lvtkFiltersParallel-pv5.0 -lvtkFiltersModeling-pv5.0 -lvtkRenderingCore-pv5.0 -lvtkFiltersExtraction-pv5.0 -lvtkFiltersStatistics-pv5.0 -lvtkImagingFourier-pv5.0 -lvtkImagingCore-pv5.0 -lvtkalglib-pv5.0 -lvtkFiltersGeometry-pv5.0 -lvtkFiltersSources-pv5.0 -lvtkFiltersGeneral-pv5.0 -lvtkCommonComputationalGeometry-pv5.0 -lvtkCommonColor-pv5.0 -lvtkFiltersProgrammable-pv5.0 -lvtkPVVTKExtensionsCore-pv5.0 -lvtkPVCommon-pv5.0 -lvtkClientServer-pv5.0 -lvtkIOXMLParser-pv5.0 -lvtkFiltersCore-pv5.0 -lvtkParallelMPI-pv5.0 -lvtkParallelCore-pv5.0 -lvtkIOLegacy-pv5.0 -lvtkCommonExecutionModel-pv5.0 -lvtkCommonDataModel-pv5.0 -lvtkCommonSystem-pv5.0 -lvtkCommonTransforms-pv5.0 -lvtkCommonMisc-pv5.0 -lvtkCommonMath-pv5.0 -lprotobuf -lpthread  -lvtkPythonInterpreter-pv5.0 -lvtkWrappingPython27Core-pv5.0 -lvtkCommonCore-pv5.0 -lpython2.7 -lvtksys-pv5.0 -Wl,-rpath,$(PARAVIEW_LIB_DIR) -Wl,-rpath-link,$(PARAVIEW_LIB_DIR)
PARAVIEW_INCLUDE=-I$(PARAVIEW_DIR)/include/paraview-5.0

libmesh_CXXFLAGS += $(PARAVIEW_INCLUDE) -I../../xdmf

# libmesh_LDFLAGS  += $(shell ${PARAVIEW_CONGIG} --libs vtkPVPythonCatalyst --libs vtkWebPython --libs vtkPVCatalyst) 

# libmesh_LDFLAGS += /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkPVPythonCatalyst-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkPVCatalyst-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkPVServerManagerApplication-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkPVServerManagerCore-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkPVServerImplementationCore-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkPVClientServerCoreCore-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkFiltersParallel-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkFiltersModeling-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkRenderingCore-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkFiltersExtraction-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkFiltersStatistics-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkImagingFourier-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkImagingCore-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkalglib-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkFiltersGeometry-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkFiltersSources-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkFiltersGeneral-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkCommonComputationalGeometry-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkCommonColor-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkFiltersProgrammable-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkPVVTKExtensionsCore-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkPVCommon-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkClientServer-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkIOXMLParser-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkFiltersCore-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkParallelMPI-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkParallelCore-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkIOLegacy-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkIOCore-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkCommonExecutionModel-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkCommonDataModel-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkCommonSystem-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkCommonTransforms-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkCommonMisc-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkCommonMath-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libprotobuf.so -lpthread /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkPythonInterpreter-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkWrappingPython27Core-pv5.0.so.1 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtkCommonCore-pv5.0.so.1 -lpython2.7 /home/camata/local/src/ParaView-v5.0.1-source/Build/lib/libvtksys-pv5.0.so.1 -ldl -lmpi_cxx -Wl,-rpath,/home/camata/local/src/ParaView-v5.0.1-source/Build/lib -Wl,-rpath-link,/home/camata/local/src/ParaView-v5.0.1-source/Build/lib 
libmesh_LDFLAGS+=$(PARAVIEW_LIBS)

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
