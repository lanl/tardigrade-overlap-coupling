# DOFProjection, a collection of routines to project the degrees of freedom
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   June 10, 2020
#
# This is the makefile for test_DOFProjection.cpp. It generates the 
# tests for the projection of the degrees of freedom

include ../../../../config.mk

# List of executables
EXECUTABLE=test_inputFileProcessor
TESTOBJECT=test_inputFileProcessor.o
TESTSOURCE=test_inputFileProcessor.cpp

OBJECT=inputFileProcessor.o
SOURCE=../../inputFileProcessor.cpp
HEADER=../../inputFileProcessor.h

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
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) $^ -o $@ $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++

$(TESTOBJECT): $(TESTSOURCE) $(HEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++

$(OBJECT): $(SOURCE) $(HEADER) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(ERROROBJECT): $(ERRORSOURCE) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)

.PHONY: all clean
