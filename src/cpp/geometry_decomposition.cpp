/******************************************************************************
*                          geometry_decomposition.cpp                         *
===============================================================================
* A collection of functions and related utilities intended to help decompose  *
* simple geometric objects into easy to integrate subvolumes.                 *
*******************************************************************************
*/

#include "geometry_decomposition.h"

namespace gDecomp{

    std::vector< matrixType > getTets(const vectorType &p, const matrixType &nodes){
        /*!
         * Get the tetrahedra as lists of points which all use point p
         * as one of the vertices. Nodes should be an ordered list of 
         * points from a flat surface
         * 
         * :param const vectorType &p: Point p
         * :param const matrixType &nodes: A list of nodes on a surface (nnodes, ndim)
         */

        if (nodes.size() < 3){
            std::cerr << "Error: number of nodes must be at least three\n";
            assert(1==0);
        }

        //Compute the centroid of the points on the face
        vectorType faceCentroid;
        vectorTools::computeMean(nodes, faceCentroid);

        std::vector< matrixType > tets;
        tets.reserve(nodes.size());

        unsigned int i=0;
        for (auto it=nodes.begin(); it!=nodes.end()-1; it++, i++){
            tets.push_back({p, faceCentroid, nodes[i], nodes[i+1]});
        }

        tets.push_back({p, faceCentroid, nodes[nodes.size()-1], nodes[0]});

        return tets;
    }

    floatType getTetVolume(const matrixType &tet){
        /*!
         * Compute and return the volume of a tetrahedron
         * 
         * :param const matrixType &tet: The nodes of the tetrahedron
         */

        if (tet.size() != 4){
            std::cerr << "Error: A tetrahedron must be defined by 4 nodes not " << tet.size() << "\n";
            assert(1==0);
        }

        //Define the tetrahedron's sides
        vectorType s1 = tet[1] - tet[0];
        vectorType s2 = tet[2] - tet[0];
        vectorType s3 = tet[3] - tet[0];

        //Compute the volume
        return std::abs(vectorTools::dot(vectorTools::cross(s1, s2), s3))/6;
    }

    int getUnitToTetMap(const matrixType &nodes, matrixType &A, vectorType &d){
        /*!
         * Get the map from the unit tetrahedron to an arbitrary
         * tetrahedron with the given nodes.
         * 
         * :param const matrixType &nodes: The nodes of the tetrahedron 
         *     (nnoddes, ndim)
         * :param matrixType &A: The mapping matrix
         * :param vectorType &d: The shift vector
         */

        if (nodes.size() != 4){
            std::cerr << "Error: A tetrahedron must be defined by 4 nodes not " << nodes.size() << "\n";
            assert(1==0);
        }

        d = nodes[0];

        if (d.size() != 3){
            std::cerr << "Error: A tetrahedron must be defined in 3D not " << d.size() << "\n";
            assert(1==0);
        }

        A.resize(nodes[0].size(), {0, 0, 0});

        for (unsigned int i=1; i<4; i++){
            for (unsigned int j=0; j<3; j++){
                A[j][i-1] = nodes[i][j] - d[j];
            }
        }
        return 0;
    }
}
