######################################################################
#
# Template libMesh application Makefile
CATALYST=1
PROVENANCE=1
MACOS?=1

#include config/linda.dell

# include the library options determined by configure
include $(LIBMESH_DIR)/Make.common

TARGET_SUFFIX := $(METHOD)

###############################################################################
# File management.  This is where the source, header, and object files are
# defined

# source files
srcfiles 	:= $(wildcard *.C) $(wildcard src/common/*.C) $(wildcard src/flow/*.C) $(wildcard src/io/*.C) $(wildcard src/mesh_mov/*.C) $(wildcard src/transport/*.C) $(wildcard src/time_stepping/*.C)

libmesh_CXXFLAGS += -I./include/common -I./include/flow  -I./include/io -I./include/mesh_mov -I./include/transport -I./include/time_stepping -I./contrib 

srcfiles +=  $(wildcard src/provenance/*.C) $(wildcard contrib/dfanalyzer/*.C) $(wildcard contrib/rapidjson/*.C)
libmesh_CXXFLAGS +=  -I./contrib/dfanalyzer/*.h -I./include/provenance

ifeq (${PROVENANCE}, 1) 
libmesh_CXXFLAGS += -DPROVENANCE
endif

libmesh_LDFLAGS += -lcurl 

ifeq (${CATALYST}, 1) 
PARAVIEW_INCLUDE = $(shell /Users/vitor/Documents/program/paraview-5.4.0/bin/paraview-config --include vtkPVCatalyst vtkPVPythonCatalyst vtkParaViewWeb)
PARAVIEW_LIBS    = $(shell /Users/vitor/Documents/program/paraview-5.4.0/bin/paraview-config --libs vtkUtilitiesPythonInitializer vtkPVCatalyst vtkPVPythonCatalyst vtkParaViewWeb) 
libmesh_CXXFLAGS += $(PARAVIEW_INCLUDE) 
libmesh_CXXFLAGS +=-I./include/catalyst -DUSE_CATALYST
libmesh_LDFLAGS +=$(PARAVIEW_LIBS)
srcfiles +=  $(wildcard src/catalyst/*.C)
endif

ifeq (${BEDMOV}, 1) 
libmesh_CXXFLAGS += -DMESH_MOVIMENT
endif

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################

target     := ./sediment-$(TARGET_SUFFIX)
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
