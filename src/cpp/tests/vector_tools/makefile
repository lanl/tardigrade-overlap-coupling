# vector_tools, the vector tools library ( test file is depreciated )
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   July 13, 2018
#
# This is the makefile for test_vector_tools.cpp. It generates the 
# test files for vector_tools.cpp

include ../../../../config.mk

# List of executables
EXECUTABLES=test_vector_tools
OBJECTS=test_vector_tools.o
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

test_vector_tools: $(OBJECTS)
	$(CXX) $(LOCAL_CFLAGS) $(INC) $(LIB) $(DBG) $^ -o $@ $(CFLAGS) -lstdc++ -std=gnu++11

test_vector_tools.o: test_vector_tools.cpp $(VECTORHEADER) $(VECTORSOURCE)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11

clean:
	rm -f $(EXECUTABLES)
	rm -f $(OBJECTS)

.PHONY: all clean
