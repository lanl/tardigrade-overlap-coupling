# Overlap coupling, tools for performing partial and total overlap
# coupling between a micromorphic macro-scale and micro micro scale.
#
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   July 10, 2020
#
# This is the common configuration file for all of the included makefiles

# Check for icc compiler or default to g++
ICC_EXIST:=$(shell which icc)
ifdef ICC_EXIST
    CXX=icc
    CFLAGS=-std=c++11 -Wall -Wextra -ansi -pedantic -I. -O3 -fmax-errors=5 -lyaml-cpp
else
    CXX=g++
    CFLAGS=-std=gnu++11 -Wall -Wextra -ansi -pedantic -I. -O3 -fmax-errors=5 -lyaml-cpp
endif

# Location of the Eigen library
EIGEN=-I$(abspath $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))/../eigen)

# Set the root directory
ROOTDIR:=$(abspath $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))/../)

# Location of the quickhull source code
QHULL=/projects/nathanm/micromorphic/micromorphic_library/external_repositories/quickhull/

#Python includes ( TODO: Make this more automatic )
PYTHON_INCLUDE:=/home/nathan/anaconda3/include
PYTHON_LIB:=/home/nathan/anaconda3/lib

# Add the location of the error_tools to the include and library
ERRORSOURCE = $(ROOTDIR)/error_tools/src/cpp/error_tools.cpp
ERRORHEADER = $(ROOTDIR)/error_tools/src/cpp/error_tools.h
INC=-I$(ROOTDIR)/error_tools/src/cpp
LIB=-L$(ROOTDIR)/error_tools/src/cpp

# Add the location of the vector_tools to the include and library
VECTORSOURCE = $(ROOTDIR)/vector_tools/src/cpp/vector_tools.cpp
VECTORHEADER = $(ROOTDIR)/vector_tools/src/cpp/vector_tools.h
INC+=-I$(ROOTDIR)/vector_tools/src/cpp
LIB+=-L$(ROOTDIR)/vector_tools/src/cpp

# Add the location of the input file processor definitions to the include and library
INPUTFILEPROCESSORSOURCE = $(ROOTDIR)/overlap_coupling/src/cpp/inputFileProcessor.cpp
INPUTFILEPROCESSORHEADER = $(ROOTDIR)/overlap_coupling/src/cpp/inputFileProcessor.h

# Add the location of the DOF projection definitions to the include and library
DOFPROJECTIONSOURCE = $(ROOTDIR)/overlap_coupling/src/cpp/DOFProjection.cpp
DOFPROJECTIONHEADER = $(ROOTDIR)/overlap_coupling/src/cpp/DOFProjection.h

# Add the filenames of the data-file interface to the include and library
DATAFILEINTERFACESOURCE = $(ROOTDIR)/overlap_coupling/src/cpp/dataFileInterface.cpp
DATAFILEINTERFACEHEADER = $(ROOTDIR)/overlap_coupling/src/cpp/dataFileInterface.h

# Add the filenames of the data-file interface to the include and library
ELEMENTSOURCE = $(ROOTDIR)/overlap_coupling/src/cpp/element.cpp
ELEMENTHEADER = $(ROOTDIR)/overlap_coupling/src/cpp/element.h

# Add the location of solver_tools to the include and library
SOLVERSOURCE = $(ROOTDIR)/solver_tools/src/cpp/solver_tools.cpp
SOLVERHEADER = $(ROOTDIR)/solver_tools/src/cpp/solver_tools.h
INC+=-I$(ROOTDIR)/solver_tools/src/cpp
LIB+=-L$(ROOTDIR)/solver_tools/src/cpp

# Add the location of the source and header files for volumeReconstruction
VOLUMERECONSTRUCTIONSOURCE = $(ROOTDIR)/overlap_coupling/src/cpp/volumeReconstruction.cpp
VOLUMERECONSTRUCTIONHEADER = $(ROOTDIR)/overlap_coupling/src/cpp/volumeReconstruction.h

#Add the location of the micromorphic tools files
MICROMORPHICTOOLSSOURCE = $(ROOTDIR)/micromorphic_tools/src/cpp/micromorphic_tools.cpp
MICROMORPHICTOOLSHEADER = $(ROOTDIR)/micromorphic_tools/src/cpp/micromorphic_tools.h
INC+=-I$(ROOTDIR)/micromorphic_tools/src/cpp
LIB+=-I$(ROOTDIR)/micromorphic_tools/src/cpp

#Add the location of the constitutive tools files
CONSTITUTIVETOOLSSOURCE = $(ROOTDIR)/constitutive_tools/src/cpp/constitutive_tools.cpp
CONSTITUTIVETOOLSHEADER = $(ROOTDIR)/constitutive_tools/src/cpp/constitutive_tools.h
INC+=-I$(ROOTDIR)/constitutive_tools/src/cpp
LIB+=-I$(ROOTDIR)/constitutive_tools/src/cpp

#Add the location of the balance equation file
BALANCEEQUATIONSSOURCE = $(ROOTDIR)/micromorphic_element/src/cpp/balance_equations.cpp
BALANCEEQUATIONSHEADER = $(ROOTDIR)/micromorphic_element/src/cpp/balance_equations.h
INC+=-I$(ROOTDIR)/micromorphic_element/src/cpp
LIB+=-I$(ROOTDIR)/micromorphic_element/src/cpp

# Add the location of the YAML headers and library
INC+= -I$(ROOTDIR)/yaml-cpp/include
LIB+= -L$(ROOTDIR)/yaml-cpp/build/

# Add the location of the XDMF headers, library, and associated files
INC += -I$(ROOTDIR)/xdmf -I$(ROOTDIR)/xdmf/build
INC += -I$(ROOTDIR)/xdmf/core -I$(ROOTDIR)/xdmf/build/core
INC += -I$(PYTHON_INCLUDE)/libxml2

LIB += -L$(ROOTDIR)/xdmf/build/lib -lXdmf -lXdmfCore -lXdmfUtils
LIB += -L$(PYTHON_LIB)/libxml2 -lxml2
LIB += -L$(PYTHON_LIB) -lhdf5 -ltiff

# Add the location of the Boost headers
INC+= -I$(BOOST_ROOT)

# Add the location of the voro++ libraries
INC+=-I$(ROOTDIR)/voro++/voro++
LIB+=-L$(ROOTDIR)/voro++/voro++

# The python command
PYTHON=/apps/anaconda3/bin/python
