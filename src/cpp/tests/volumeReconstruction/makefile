# volumeReconstruction, a collection of routines to reconstruct volume information
# from pointsets
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   July 3, 2020
#
# This is the makefile for test_volumeReconstruction.cpp. It generates the 
# tests for the reconstruction of volume information from pointsets

include ../../../../config.mk

# List of executables
EXECUTABLE=test_volumeReconstruction
TESTOBJECT=test_volumeReconstruction.o
TESTSOURCE=test_volumeReconstruction.cpp

OBJECT=volumeReconstruction.o
SOURCE=../../volumeReconstruction.cpp
HEADER=../../volumeReconstruction.h

ERROROBJECT=error_tools.o

OBJECTS = $(OBJECT) $(TESTOBJECT) $(ERROROBJECT)

#Set the mode due to the method variable
DBG =
ifeq ($(METHOD),dbg)
  DBG = -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -g
endif
ifeq ($(METHOD),oprof)
  DBG = -pg
endif

INC+=-I../..
# Makefile rules
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) $^ -o $@ $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++ -lXdmf -lXdmfCore -lXdmfUtils -ltiff -lxml2 -lhdf5

$(TESTOBJECT): $(TESTSOURCE) $(HEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++

$(OBJECT): $(SOURCE) $(HEADER) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11 -lXdmf -lXdmfCore -lXdmfUtils

$(ERROROBJECT): $(ERRORSOURCE) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)

.PHONY: all clean
