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
#define USE_EIGEN
#include<vector_tools.h>

namespace gDecomp{

    //Type definitions
    typedef double floatType;
    typedef std::vector< floatType > vectorType;
    typedef std::vector< vectorType > matrixType;
    typedef std::pair< vectorType, vectorType > faceType;

    //Geometry decomposition tools
    std::vector< matrixType > getTets(const vectorType &p, const matrixType &nodes);

    floatType getTetVolume(const matrixType &tet);

    int getUnitToTetMap(const matrixType &nodes, matrixType &A, vectorType &d);

    int getTetQuadrature(const unsigned int order, matrixType &points, vectorType &weights);

    int findPointsOnFace(const vectorType &n, const vectorType &q, const matrixType &points,
                         std::vector< unsigned int > &surfacePoints,
                         double tolr=1e-9, double tola=1e-9);

    int orderPlanarPoints(const matrixType &points, std::vector< unsigned int > &orderedIndices);

    int getFacePoints(const std::vector< faceType > &faces, const matrixType &points,
        std::vector< std::vector< unsigned int > > &indexFaces);

    int volumeToTets(const std::vector< faceType > &faces, const matrixType &points,
        std::vector< matrixType > &tets);

    int findMidpoints(const vectorType &p, const matrixType &points, matrixType &midpoints,
        double tolr=1e-9, double tola=1e-9);

    int findPointOfIntersection(const std::vector< faceType > &planes, vectorType &point, bool &solveFlag);

    int findAllPointsOfIntersection(const std::vector< faceType > &planes, matrixType &intersectionPoints,
        floatType tolr=1e-9, floatType tola=1e-9);

    bool isDuplicate(const vectorType &point, const matrixType &points);

    int determineInteriorPoints(const vectorType &pInside, const matrixType &points, const std::vector< faceType > &faces,
        std::vector< unsigned int > &interiorPoints, floatType tolr=1e-9, floatType tola=1e-9);

    int midpointsToFaces(const vectorType &p, const matrixType &midpoints, std::vector< faceType > &faces);

    int getVolumeSubdomainAsTets(const unsigned int index, const matrixType &domainPoints, const std::vector< faceType > &faces,
        std::vector< matrixType > &subdomainTets);

    int mapLocalTetPointsToGlobal(const matrixType &tet, const matrixType &localPoints, matrixType &globalPoints);
}

#endif
