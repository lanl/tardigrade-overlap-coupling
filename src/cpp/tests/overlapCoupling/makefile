# overlapCoupling, a collection of routines to performt the overlap coupling operations
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   June 10, 2020
#
# This is the makefile for test_overlapCoupling.cpp. It generates the 
# tests for the overlap coupling.

include ../../../../config.mk

# List of executables
EXECUTABLE=test_overlapCoupling
TESTOBJECT=test_overlapCoupling.o
TESTSOURCE=test_overlapCoupling.cpp

OBJECT=overlapCoupling.o
SOURCE=../../overlapCoupling.cpp
HEADER=../../overlapCoupling.h

INPUTFILEPROCESSOROBJECT=inputFileProcessor.o

DOFPROJECTIONOBJECT=DOFProjection.o

DATAFILEINTERFACEOBJECT=dataFileInterface.o

ERROROBJECT=error_tools.o

ELEMENTOBJECT=element.o

SUPPORTHEADERS  = $(ERRORHEADER) $(VECTORHEADER) $(ELEMENTHEADER) $(DOFPROJECTIONHEADER) $(DATAFILEINTERFACEHEADER)
SUPPORTHEADERS +=$(INPUTFILEPROCESSORHEADER)

SUPPORTOBJECTS = $(ERROROBJECT) $(DATAFILEINTERFACEOBJECT) $(DOFPROJECTIONOBJECT) $(INPUTFILEPROCESSOROBJECT) $(ELEMENTOBJECT)
OBJECTS = $(OBJECT) $(TESTOBJECT) $(SUPPORTOBJECTS)

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

$(INPUTFILEPROCESSOROBJECT): $(INPUTFILEPROCESSORSOURCE) $(INPUTFILEPROCESSORHEADER) $(ERRORHEADER) $(VECTORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++

$(DOFPROJECTIONOBJECT): $(DOFPROJECTIONSOURCE) $(DOFPROJECTIONHEADER) $(ERRORHEADER) $(VECTORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++

$(DATAFILEINTERFACEOBJECT): $(DATAFILEINTERFACESOURCE) $(DATAFILEINTERFACEHEADER) $(ERRORHEADER) $(VECTORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++

$(OBJECT): $(SOURCE) $(HEADER) $(SUPPORTHEADERS)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(ERROROBJECT): $(ERRORSOURCE) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(ELEMENTOBJECT): $(ELEMENTSOURCE) $(ELEMENTHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)

.PHONY: all clean
