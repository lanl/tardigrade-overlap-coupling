# inputFileProcessor, a collection of routines for processing input files
# for overlap coupling
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   June 10, 2020
#
# This is the makefile for test_inputFileProcessor.cpp. It generates the 
# tests for the processing the input file

include ../../../../config.mk

# List of executables
EXECUTABLE=test_inputFileProcessor
TESTOBJECT=test_inputFileProcessor.o
TESTSOURCE=test_inputFileProcessor.cpp

OBJECT=inputFileProcessor.o
SOURCE=../../inputFileProcessor.cpp
HEADER=../../inputFileProcessor.h

DATAFILEINTERFACEOBJECT=dataFileInterface.o

GENERATEXDMFDATAOBJECT=generateXMDFData.o

ERROROBJECT=error_tools.o

OBJECTS = $(OBJECT) $(TESTOBJECT) $(ERROROBJECT) $(DATAFILEINTERFACEOBJECT) $(GENERATEXDMFDATAOBJECT)

#Set the mode due to the method variable
DBG =
YAMLFLG = -lyaml-cpp
ifeq ($(METHOD),dbg)
  DBG = -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -g
  YAMLFLG = -lyaml-cppd
endif
ifeq ($(METHOD),oprof)
  DBG = -pg
endif

INC+=-I../..

# Makefile rules
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) $^ -o $@ $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++ -lXdmf -lXdmfCore -lXdmfUtils -ltiff -lxml2 -lhdf5 $(YAMLFLG)

$(TESTOBJECT): $(TESTSOURCE) $(HEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++

$(GENERATEXDMFDATAOBJECT): $(GENERATEXDMFDATASOURCE) $(GENERATEXDMFDATAHEADER) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++

$(DATAFILEINTERFACEOBJECT): $(DATAFILEINTERFACESOURCE) $(DATAFILEINTERFACEHEADER) $(ERRORHEADER) $(VECTORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=gnu++11 -lvoro++

$(OBJECT): $(SOURCE) $(HEADER) $(ERRORHEADER) $(VECTORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

$(ERROROBJECT): $(ERRORSOURCE) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)

.PHONY: all clean
