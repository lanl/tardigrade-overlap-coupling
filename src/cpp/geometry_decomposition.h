/******************************************************************************
*                          geometry_decomposition.h                           *
===============================================================================
* A collection of functions and related utilities intended to help decompose  *
* simple geometric objects into easy to integrate subvolumes.                 *
*******************************************************************************
*/

#ifndef GEOMETRY_DECOMPOSITION_H
#define GEOMETRY_DECOMPOSITION_H

#include<stdio.h>
#include<iostream>
#include<exception>
#include<fstream>
#include<vector>
#include<map>
#include<math.h>
#include<assert.h>
#include<string.h>
#include<vector_tools.h>

namespace gDecomp{

    //Type definitions
    typedef double floatType;
    typedef std::vector< floatType > vectorType;
    typedef std::vector< vectorType > matrixType;

    //Geometry decomposition tools
    std::vector< matrixType > getTets(const vectorType &p, const matrixType &nodes);
    floatType getTetVolume(const matrixType &tet);
    int getUnitToTetMap(const matrixType &nodes, matrixType &A, vectorType &d);

}

#endif
