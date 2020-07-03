# overlap_coupling, Routines for processing DNS and macroscale
# data for overlap coupling.
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   July 13, 2018
#
# This is the makefile for test_overlap_coupling.cpp. It generates test
# functions for overlap_coupling.cpp

include ../../../../config.mk

# List of executables
EXECUTABLES=test_overlap_coupling
OBJECTS=overlap_coupling.o test_overlap_coupling.o quickhull.o element.o assembly.o
INC+=-I../..
LOCAL_CFLAGS=-ggdb -O3 -fpic -Wall

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

test_overlap_coupling: $(OBJECTS)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) $^ -o $@ $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++

test_overlap_coupling.o: test_overlap_coupling.cpp ../../overlap_coupling.h ../../occonfiguration.h
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++

quickhull.o: $(QHULL)/QuickHull.cpp $(QHULL)/*.hpp
	$(CXX) $(CFLAGS) $(LOCAL_CFLAGS) $(INC) $(LIB) $(DBG) -c -o $@ $< -std=gnu++11

overlap_coupling.o: ../../overlap_coupling.cpp ../../overlap_coupling.h ../../geometry_decomposition.h ../../occonfiguration.h
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++

geometry_decomposition.o: geometry_decomposition.cpp geometry_decomposition.h
	$(CXX) $(CFLAGS) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< -lstdc++ -std=gnu++11

element.o: ../../element.cpp ../../element.h
	$(CXX) $(CFLAGS) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< -lstdc++ -std=gnu++11

assembly.o: ../../assembly.cpp ../../assembly.h
	$(CXX) $(CFLAGS) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< -lstdc++ -std=gnu++11

clean:
	rm -f $(EXECUTABLES)
	rm $(OBJECTS)

.PHONY: all clean
