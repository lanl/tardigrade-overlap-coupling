# Micromorphic filter, a variationally based filter for DNS
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   July 13, 2018
#
# This is the makefile for processDNS.cpp. It generates the 
# micromorphic filter

include ../../../../config.mk

# List of executables
EXECUTABLES=test_overlap_coupling
OBJECTS=overlap_coupling.o test_overlap_coupling.o
INC+=-I../..
LOCAL_CFLAGS=-ggdb

# Makefile rules
all: $(EXECUTABLES)

test_overlap_coupling: $(OBJECTS)
	$(CXX) $(CFLAGS) $(LOCAL_CFLAGS) $(VPPLIB) $(VPPINC) $(EIGEN) $(INC) $^ -o $@ -lvoro++

test_overlap_coupling.o: test_overlap_coupling.cpp ../../overlap_coupling.h
	$(CXX) $(CFLAGS) $(LOCAL_CFLAGS) $(VPPLIB) $(VPPINC) $(EIGEN) $(INC) -c -o $@ $< -lvoro++

overlap_coupling.o: ../../overlap_coupling.cpp ../../overlap_coupling.h
	$(CXX) $(CFLAGS) $(LOCAL_CFLAGS) $(VPPLIB) $(VPPINC) $(EIGEN) $(INC) -c -o $@ $< -lvoro++

clean:
	rm -f $(EXECUTABLES)
	rm $(OBJECTS)

.PHONY: all clean
