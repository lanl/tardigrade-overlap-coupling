# dataFileInterface, a collection of routines to interface with datafiles
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   June 10, 2020
#
# This is the makefile for test_dataFileInterface.cpp. It generates the 
# tests for the interface with the datafile

include ../../../../config.mk

# List of executables
EXECUTABLE=test_dataFileInterface
TESTOBJECT=test_dataFileInterface.o
TESTSOURCE=test_dataFileInterface.cpp

OBJECT=dataFileInterface.o
SOURCE=../../dataFileInterface.cpp
HEADER=../../dataFileInterface.h

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
