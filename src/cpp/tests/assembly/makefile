# assembly, routines to assemble FEA systems
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   July 13, 2018
#
# This is the makefile for test_assembly.cpp. It generates the 
# test functions for the assembly routines
#
# NOTE: This file is not used in the current overlap coupling framework
#       but we will continue to track it for now.

include ../../../../config.mk

# List of executables
EXECUTABLES=test_assembly

ASSEMBLYOBJECT=assembly.o
ELEMENTOBJECT=element.o
ERROROBJECT=error_tools.o
OBJECTS=$(ELEMENTOBJECT) assembly.o test_assembly.o $(ERROROBJECT)
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

test_assembly: $(OBJECTS)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) $^ -o $@ $(CFLAGS) -lstdc++ -std=gnu++11

test_assembly.o: test_assembly.cpp ../../assembly.h $(ELEMENTHEADER) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11

$(ASSEMBLYOBJECT): ../../assembly.cpp ../../assembly.h
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11

$(ELEMENTOBJECT): $(ELEMENTSOURCE) $(ELEMENTHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11

$(ERROROBJECT): $(ERRORSOURCE) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

clean:
	rm -f $(EXECUTABLES)
	rm $(OBJECTS)

.PHONY: all clean
