/*!
===============================================================================
|                              overlap_coupling.h                             |
===============================================================================
| The header file for the overlap coupling library. This library provides the |
| classes, functions, and methods required to compute the required terms for  |
| the micro/meso-scale to macro-scale coupling following the micromorphic     |
| continuum mechanics framework.                                              |
===============================================================================
*/

#ifndef OVERLAP_COUPLING_H
#define OVERLAP_COUPLING_H

#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<math.h>
#include<assert.h>
#include<string.h>

#include "quickhull.h"

typedef std::vector< std::vector< double > > vecOfvec;

typedef std::map< std::vector< double >, std::vector< double > > planeMap;

namespace overlap{

    //!===
    //! | Classes
    //!===

    class OverlapCoupling{
        /*!
        The primary class for performing the overlap coupling operation. This 
        object will return the integration weights for volume integrals, area  
        weighted normal vectors for surface integrals, and also constructs 
        the shape function matrix which can be used to construct the projectors.
        */

        public:
            OverlapCoupling();
            OverlapCoupling(const vecOfvec &local_coordinates);

            //! > Interface to 3D-quickhull
            qh_vertex_t map_vector_to_quickhull(const std::vector< double > &vector) const;
            std::vector< double > map_quickhull_to_vector(const qh_vertex_t &vertex) const;
            void map_vectors_to_quickhull(const vecOfvec &vectors, std::vector< qh_vertex_t > &vertices) const;
            void map_quickhull_to_vectors(const std::vector< qh_vertex_t > &vertices, vecOfvec &vectors) const;
            void extract_mesh_info(const qh_mesh_t &mesh, vecOfvec &normals, vecOfvec &points) const;
            void compute_node_bounds(const vecOfvec &coordinates, planeMap &planes, const double tolr=1e-6, const double tola=1e-6) const;

            //! > Interface to defined quantities
            void get_element_planes(planeMap &planes) const;

        protected:
            vecOfvec local_coordinates;
            planeMap element_planes;
            planeMap dns_planes;

            void compute_element_bounds();
            void compute_dns_bounds(const vecOfvec &dns_coordinates);
            planeMap compute_unique_planes(const vecOfvec &normals, const vecOfvec &points, const double tolr=1e-6, const double tola=1e-6) const;


    };

    class ParsedData{
        /*!
        Class which stores data read from an input file. Used for testing.
        */

        public:
            vecOfvec global_nodes;
            vecOfvec local_nodes;
            std::vector< unsigned int > node_numbers;
            std::vector< double > volumes;
            std::vector< double > densities;
            vecOfvec coordinates;

            //! > Constructors
            ParsedData(){}

            ParsedData(vecOfvec _global_nodes, vecOfvec _local_nodes,
                       std::vector< unsigned int > _node_numbers, std::vector< double > _volumes,
                       std::vector< double > _densities, vecOfvec _coordinates){
                global_nodes = _global_nodes;
                local_nodes = _local_nodes;
                node_numbers = _node_numbers;
                volumes = _volumes;
                densities = _densities;
                coordinates = _coordinates;
            }
    };
    
   //!===
   //! | Functions
   //!===

    ParsedData read_data_from_file(std::string filename);

    std::vector< std::string > split(std::string s, char delimiter);

    template<typename Out>
    void split(const std::string &s, char delimiter, Out result);

    double dot(const std::vector< double > &a, const std::vector< double > &b);

    std::vector< double > cross(const std::vector< double > &a, const std::vector< double > &b);

    bool fuzzy_equals(const double a, const double b, const double tolr=1e-6, const double tola=1e-6);

    bool compare_vector_directions(const std::vector< double > &v1, const std::vector< double > &v2, const double tolr=1e-6, const double tola=1e-6);
}

#endif
