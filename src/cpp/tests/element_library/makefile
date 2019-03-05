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
EXECUTABLES=test_element_library
OBJECTS=element_library.o test_element_library.o

INC=-I../..
LOCAL_CFLAGS=-ggdb#-fprofile-arcs -ftest-coverage

# Makefile rules
all: $(EXECUTABLES)

test_element_library: $(OBJECTS)
	$(CXX) $(CFLAGS) $(LOCAL_CFLAGS) $(VPPLIB) $(VPPINC) $(EIGEN) $(INC) $^ -o $@ -lvoro++

test_element_library.o: test_element_library.cpp ../../element_library.h
	$(CXX) $(CFLAGS) $(LOCAL_CFLAGS) $(VPPLIB) $(VPPINC) $(EIGEN) $(INC) -o $@ -c $< -lvoro++

element_library.o: ../../element_library.cpp ../../element_library.h
	$(CXX) $(CFLAGS) $(LOCAL_CFLAGS) $(VPPLIB) $(VPPINC) $(EIGEN) $(INC) -o $@ -c $< -lvoro++

clean:
	rm -f $(EXECUTABLES)
	rm $(OBJECTS)

.PHONY: all clean
