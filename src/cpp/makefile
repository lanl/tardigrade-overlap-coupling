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
OBJECTS=overlap_coupling.o quickhull.o
#LOCAL_CFLAGS=-I. -O3 -flto -fpic -Wall
LOCAL_CFLAGS=-I. -O3 -fpic -Wall

#Set the mode due to the method variable
DBG = 
ifeq ($(METHOD),dbg)
  DBG = -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -g
endif
ifeq ($(METHOD),oprof)
  DBG = -pg
endif

#Define the library name
SONAME=liboverlap.so
SOVERSION=liboverlap.so.1

# Makefile rules
all: $(SONAME)

$(SONAME): $(SOVERSION)
	ln -s -f $(SOVERSION) $@

$(SOVERSION): $(OBJECTS)
	$(CC) $(STD) -shared -Wl,-soname,$(SOVERSION) -o $@ $^ $(CFLAGS) $(LOCAL_CFLAGS) $(ERRORFLG) $(DBG) $(INC) $(LIB) -std=gnu++11 -lvoro++

quickhull.o: $(QHULL)/QuickHull.cpp $(QHULL)/*.hpp
	$(CXX) $(CFLAGS) $(LOCAL_CFLAGS) $(INC) $(LIB) $(DBG) -c -o $@ $< -std=gnu++11

overlap_coupling.o: overlap_coupling.cpp overlap_coupling.h
	$(CXX) $(CFLAGS) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< -std=gnu++11 -lvoro++

clean:
	rm $(OBJECTS) $(SOVERSION) $(SONAME)

.PHONY: all clean
