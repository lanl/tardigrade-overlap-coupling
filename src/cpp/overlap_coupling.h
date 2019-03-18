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

//#include "voro++.hh"

typedef std::vector< std::vector< double > > vecOfvec;
typedef std::map< std::vector< double >, std::vector< double > > planeMap;

#if CONVEXLIB == QUICKHULL
    #include "3d-quickhull/quickhull.h"
    typedef qh_vertex_t vertex_t;
    typedef qh_mesh_t mesh_t;
#elif CONVEXLIB == CONVHULL_3D
    #include "conhull_3d/convhull_3d.h"
    typedef ch_vertex vertex_t;
    typedef std::pair< std::vector< int >, std::vector< vertex_t > > mesh_t;
#elif CONVEXLIB == AKUUKKA
    #include "quickhull/QuickHull.hpp"
    typedef double FloatType;
    typedef quickhull::Vector3<FloatType> vertex_t;
    typedef quickhull::ConvexHull<FloatType> mesh_t;
#else
    #error CONVEXLIB must be defined. If defined, check that the value is supported.
#endif

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
            vertex_t map_vector_to_quickhull(const std::vector< double > &vector) const;
            std::vector< double > map_quickhull_to_vector(const vertex_t &vertex) const;
            void map_vectors_to_quickhull(const vecOfvec &vectors, std::vector< vertex_t > &vertices) const;
            void map_quickhull_to_vectors(const std::vector< vertex_t > &vertices, vecOfvec &vectors) const;
#if CONVEXLIB != AKUUKKA
            void extract_mesh_info(const mesh_t &mesh, vecOfvec &normals, vecOfvec &points) const;
#else
            void extract_mesh_info(mesh_t &mesh, vecOfvec &normals, vecOfvec &points) const;
#endif
            void compute_node_bounds(const vecOfvec &coordinates, planeMap &planes, 
                std::vector< double > &xbnds, std::vector< double > &ybnds, std::vector< double > &zbnds,
                const double tolr=1e-6, const double tola=1e-6);

            void compute_dns_bounds(const vecOfvec &dns_coordinates);

            //! > Interface to defined quantities
            const planeMap* get_element_planes() const;
            const vecOfvec* get_element_bounds() const;
            
            const planeMap* get_dns_planes() const;
            const vecOfvec* get_dns_bounds() const;

        protected:
            vecOfvec local_coordinates;
            planeMap element_planes;
            planeMap dns_planes;
            vecOfvec element_bounds;
            vecOfvec dns_bounds;

//            voro::container element_container;
//            voro::container dns_container;

            void compute_element_bounds();
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

    std::vector< double > normal_from_vertices(const vertex_t &p1, const vertex_t &p2, const vertex_t &p3);

    void print_vertex(const vertex_t &vertex);
    void print_vector(const std::vector< FloatType > &vector);
}

#endif
