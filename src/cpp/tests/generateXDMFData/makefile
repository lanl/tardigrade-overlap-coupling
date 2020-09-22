# generateXMDFData, a collection of routines to generate XDMF datafiles
# 
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   September 21, 2020
#
# This is the makefile for test_generateXDMFData.cpp. It generates the 
# tests for generating XDMF datafiles

include ../../../../config.mk

# List of executables
EXECUTABLE=test_generateXDMFData
TESTOBJECT=test_generateXDMFData.o
TESTSOURCE=test_generateXDMFData.cpp

OBJECT=generateXDMFData.o
SOURCE=$(GENERATEXDMFDATASOURCE)
HEADER=$(GENERATEXDMFDATAHEADER)

DATAFILEINTERFACEOBJECT=dataFileInterface.o

ERROROBJECT=error_tools.o

OBJECTS = $(OBJECT) $(TESTOBJECT) $(ERROROBJECT) $(DATAFILEINTERFACEOBJECT)

#Set the mode due to the method variable
DBG =
YAMLFLG = -lyaml-cpp
ifeq ($(METHOD),dbg)
  DBG = -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
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

$(OBJECT): $(SOURCE) $(HEADER) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11 -lXdmf -lXdmfCore -lXdmfUtils

$(DATAFILEINTERFACEOBJECT): $(DATAFILEINTERFACESOURCE) $(DATAFILEINTERFACEHEADER) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11 -lXdmf -lXdmfCore -lXdmfUtils

$(ERROROBJECT): $(ERRORSOURCE) $(ERRORHEADER)
	$(CXX) $(LOCAL_CFLAGS) $(EIGEN) $(INC) $(LIB) $(DBG) -c -o $@ $< $(CFLAGS) -lstdc++ -std=c++11

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJECTS)

.PHONY: all clean
