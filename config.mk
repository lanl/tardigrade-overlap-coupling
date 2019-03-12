# Micromorphic filter, a variationally based filter for DNS
#
# Author: Nathan A. Miller (LANL / CU Boulder)
# Email:  nathanm@lanl.gov
# Date:   July 13, 2018
#
# This is the common configuration file for all of the included makefiles

# C++ compiler
CXX=/opt/moose/gcc-7.3.0/bin/g++

# Flags for the C++ compiler
CFLAGS=-std=gnu++11 -Wall -ansi -pedantic -O3 -I.

# Location of the Voro++ library
VPPLIB=-L/projects/nathanm/micromorphic/voro++/voro++-0.4.6/src
VPPINC=-I/projects/nathanm/micromorphic/voro++/voro++-0.4.6/src

# Location of the Eigen library
#EIGEN=-I/projects/nathanm/usr/local/include/eigen-git-mirror

# The python command
PYTHON=/apps/anaconda3/bin/python

# Additional includes
INC=-I/projects/nathanm/micromorphic/micromorphic_library/external_repositories/3d-quickhull
