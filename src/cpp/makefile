# Micromorphic filter, a variationally based filter for DNS
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   July 13, 2018
#
# This is the makefile for processDNS.cpp. It generates the 
# micromorphic filter

include ../../config.mk

# List of executables
OBJECTS=overlap_coupling.o
LOCAL_CFLAGS=-ggdb

# Makefile rules
all: $(OBJECTS)

overlap_coupling.o: overlap_coupling.cpp overlap_coupling.h
	$(CXX) $(CFLAGS) $(LOCAL_CFLAGS) $(VPPLIB) $(VPPINC) $(EIGEN) $(INC) -c -o $@ $< -lvoro++

clean:
	rm $(OBJECTS)

.PHONY: all clean
