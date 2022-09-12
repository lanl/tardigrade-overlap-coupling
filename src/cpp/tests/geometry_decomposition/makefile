# geometry_decomposition, tools for decomposition geometry ( largely depreciated )
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   July 13, 2018
#
# This is the makefile for test_geometry_decomposition.cpp. It generates the 
# tests for the geometry_decomposition
#
# NOTE: This isn't actually used but we will maintain these tests for now.

include ../../../../config.mk

# List of executables
EXECUTABLES=test_geometry_decomposition
OBJECTS=geometry_decomposition.o test_geometry_decomposition.o
INC+=-I../..
LOCAL_CFLAGS=-ggdb

#Set the mode due to the method variable
DBG =
ifeq ($(METHOD),dbg)
  DBG = -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -g
endif
ifeq ($(METHOD),oprof)
  DBG = -pg
endif

# Makefile rules
all: $(EXECUTABLES)

test_geometry_decomposition: $(OBJECTS)
	$(CXX) $(LOCAL_CFLAGS) $(INC) $(LIB) $(DBG) $^ -o $@ $(CFLAGS) -lstdc++ -std=gnu++11

test_geometry_decomposition.o: test_geometry_decomposition.cpp ../../geometry_decomposition.h
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11

geometry_decomposition.o: ../../geometry_decomposition.cpp ../../geometry_decomposition.h $(VECTORTOOLSHEADER) $(VECTORTOOLSSOURCE)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11

clean:
	rm -f $(EXECUTABLES)
	rm $(OBJECTS)

.PHONY: all clean
