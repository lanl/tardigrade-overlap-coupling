/*!
===============================================================================
|                            overlap_coupling.cpp                             |
===============================================================================
| The source file for the overlap coupling library. This library provides the |
| classes, functions, and methods required to compute the required terms for  |
| the micro/meso-scale to macro-scale coupling following the micromorphic     |
| continuum mechanics framework.                                              |
===============================================================================
*/

#include<sstream>
#include<iterator>

#define QUICKHULL_IMPLEMENTATION
#include "overlap_coupling.h"

namespace overlap{
    //!===
    //! | Classes
    //!===

    OverlapCoupling::OverlapCoupling(){}

    OverlapCoupling::OverlapCoupling(const vecOfvec &_local_coordinates){
        /*!
        The constructor for the filter.
        */

        local_coordinates = _local_coordinates;
        compute_element_bounds();
    }

    //! > Interface to 3D-quickhull

    qh_vertex_t OverlapCoupling::map_vector_to_quickhull(const std::vector< double > &vector) const{
        /*!
        Map a vector to a vertex which can be read using 3d-quickhull.

        :param std::vector< double > vector: The incoming vector to be mapped.
        */
        qh_vertex_t vertex;
        vertex.x = vector[0];
        vertex.y = vector[1];
        vertex.z = vector[2];
        return vertex;
    }

    std::vector< double > OverlapCoupling::map_quickhull_to_vector(const qh_vertex_t &vertex) const{
        /*!
        Map a quickhull vertex to std::vector

        :param qh_vertex_t vertex: The quickhull vertex type.
        */

        std::vector< double > vector(3,0);
        vector[0] = vertex.x;
        vector[1] = vertex.y;
        vector[2] = vertex.z;
        return vector;
    }

    void OverlapCoupling::map_vectors_to_quickhull(const vecOfvec &vectors, std::vector< qh_vertex_t > &vertices) const{
        /*!
        Map a collection of vectors to a vertex which can be read by 3D-quickhull
        
        :param vecOfvec vectors: a vector of std::vector< double > representing the vector values
        :param std::vector< qh_vertex_t > vertices: The output in a format readable by 3D-quickhull
        */

        vertices.reserve(vectors.size());
        for (unsigned int i=0; i<vectors.size(); i++){
            vertices.push_back(map_vector_to_quickhull(vectors[i]));
        }
        return;
    }

    void OverlapCoupling::map_quickhull_to_vectors(const std::vector< qh_vertex_t > &vertices, vecOfvec &vectors) const{
        /*!
        Map a collection of 3D-quickhull vertices to std::vectors

        :param std::vector< qh_vertex_t > vertices: The incoming quickhull vertices
        :param vecOfvec vectors: The outgoing vector of std::vectors
        */

        vectors.reserve(vertices.size());
        for (unsigned int i=0; i<vertices.size(); i++){
            vectors.push_back(map_quickhull_to_vector(vertices[i]));
        }
    }

    void OverlapCoupling::extract_mesh_info(const qh_mesh_t &mesh, vecOfvec &normals, vecOfvec &points) const{
        /*!
        Extract the required information from the quickhull mesh

        :param qh_mesh_t mesh: The output from qh_quickhull3d
        :param vecOfvec normals: The normal vectors of the hull's facets
        :param vecOfvec points: Points on the hull's facets
        */

        //!Extract the normal information
        std::vector< qh_vertex_t > _normals(mesh.normals, mesh.normals+mesh.nnormals);

        //!Get a single point of each facet corresponding to the normal.
        std::vector< qh_vertex_t > _points;
        _points.reserve(_normals.size());

        unsigned int index = 0;
        for (unsigned int i=0; i<mesh.nindices/3; i++){
            _points.push_back(mesh.vertices[mesh.indices[index]]);
            index += 3;
        }

        //!Map the quickhull representations to vectors
        map_quickhull_to_vectors(_normals, normals);
        map_quickhull_to_vectors(_points, points);
    }

    void OverlapCoupling::compute_element_bounds(){
        /*!
        Compute the bounds of the element by constructing its convex hull.
        */

        compute_node_bounds(local_coordinates, element_planes);

        return;
    }

    void OverlapCoupling::compute_node_bounds(const vecOfvec &coordinates, planeMap &planes, const double tolr, const double tola) const{
        /*!
        Compute the bounding planes for the provided coordinates.

        :param vecOfvec coordinates: A vector of std::vectors of the nodal coordinates.
        :param planeMap planes: The resulting planes which bound the coordinates.
        */
        std::cout << "derp\n";
        //!Map the coordinates to 3D-quickhull vertices
        std::vector< qh_vertex_t > vertices;
        map_vectors_to_quickhull(coordinates, vertices);
        
        //!Construct the mesh
        qh_mesh_t mesh = qh_quickhull3d(&vertices[0], vertices.size());
        qh_mesh_export(&mesh, "tmp2.obj");

        //!Extract the relevant information
        vecOfvec normals, points;
        extract_mesh_info(mesh, normals, points);

        //!Form the planes
        planes = compute_unique_planes(normals, points);

        return;
    }

    planeMap OverlapCoupling::compute_unique_planes(const vecOfvec &normals, const vecOfvec &points, const double tolr, const double tola) const{
        /*!
        Compute which normal vectors are unique. For a convex hull a unique normal indicates a unique plane.

        :param vecOfvec normals: A collection of unit normal vectors defining surfaces
        */

        planeMap _planes;
        planeMap::iterator it;

        for (unsigned int i=0; i<normals.size(); i++){

            for (it = _planes.begin(); it != _planes.end(); it++){
                if (compare_vector_directions(it->first, normals[i])){
                    break;
                }
            }
            if (it == _planes.end()){
                _planes.insert(it, std::pair< std::vector< double >, std::vector< double > >(normals[i], points[i]));
            }
        }
        return _planes;
    }

    void OverlapCoupling::get_element_planes(planeMap &planes) const{
        /*!
        Extract the element planes
        */

        planes = element_planes;
    }

    void compute_distances(const vecOfvec &normals, const vecOfvec &points, std::vector< double > &distances){
        /*!
        Compute the distances of planes from the origin

        :param vecOfvec normals: A collection of unit normal vectors defining surfaces
        :param vecOfvec points: A collection of points on the surfaces.
        :param std::vector< double > distances: The resulting distances.

        d = n_i p_i

        Note: Negative distances imply that the origin is outside of the bounding surface
        */

        //!Make sure the normals and points have the same size
        if (normals.size() != points.size()){
            std::cout << "Error: normals and points must have the same size\n";
            assert(1==0);
        }

        distances.reserve(normals.size());
        for (unsigned int i=0; i<normals.size(); i++){
            distances.push_back(dot(normals[i], points[i]));
        }
    }

    bool compare_vector_directions(const std::vector< double > &v1, const std::vector< double > &v2, const double tolr, const double tola){
        /*!
        Compare vectors to determine if they are in the same direction.

        :param std::vector< double > v1: The first vector
        :param std::vector< double > v2: The second vector

        The dot product is computed and checked if it is nearly equal to 1.
        */

        double factor = sqrt(dot(v1, v1)*dot(v2, v2));
        double result = dot(v1, v2)/factor;
        return fuzzy_equals(result, 1, tolr, tola);
    }

    //!===
    //! | Functions
    //!===

    ParsedData read_data_from_file(std::string filename){
        /*!
        Read in formatted data from a file. Used primarily for testing purposes.
        */

        std::ifstream file;
        std::string line;
        std::vector<std::string> split_line;
        std::vector< double > tmp;
        std::vector< std::vector< double > > global_nodes;
        std::vector< std::vector< double > > local_nodes;
        std::vector< unsigned int > node_numbers;
        std::vector< double > volumes;
        std::vector< double > densities;
        vecOfvec coordinates;

        //Open the file
        file.open(filename.c_str());

        //Skip past the header
        for (unsigned int i=0; i<3; i++){std::getline(file, line);}

        //Read in the nodes of the element
        for (unsigned int i=0; i<8; i++){
            std::getline(file, line);
            split_line = split(line, ' ');
            global_nodes.push_back({});
            local_nodes.push_back({});
            for (unsigned int j=0; j<3; j++){global_nodes[i].push_back(::atof(split_line[j].c_str()));}
            for (unsigned int j=3; j<6; j++){local_nodes[i].push_back(::atof(split_line[j].c_str()));}
        }

        //Split the strings
        while (std::getline(file, line)){
            split_line = split(line, ' ');
            node_numbers.push_back(::atoi(split_line[0].c_str()));
            volumes.push_back(::atof(split_line[1].c_str()));
            densities.push_back(::atof(split_line[2].c_str()));
            coordinates.push_back({});
            for (unsigned int i=3; i<6; i++){coordinates.back().push_back(::atof(split_line[i].c_str()));}
            
        }
        return ParsedData(global_nodes, local_nodes, node_numbers, volumes, densities, coordinates);
    }

    //!String parsing techniques from https://stackoverflow.com/questions/236129/how-do-i-iterate-over-the-words-of-a-string
    
    template<typename Out>
    void split(const std::string &s, char delimiter, Out result){
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delimiter)){
            *(result++) = item;
        }
    }

    std::vector< std::string > split(std::string s, char delimiter){
        /*!
        Parse the line to extract substrings.
        */

        std::vector< std::string > output;
        split(s, delimiter, std::back_inserter(output));
        return output;
    }

    double dot(const std::vector< double > &a, const std::vector< double > &b){
        /*!
        Compute the dot product between two vectors

        :param std::vector< double > a: The first vector
        :param std::vector< double > b: The second vector
        */

        if (a.size() != b.size()){
            std::cout << "Error: vectors must have the same size.\n";
            assert(1==0);
        }

        double result = 0;
        for (unsigned int i=0; i<a.size(); i++){
            result += a[i] * b[i];
        }
        return result;
    }

    std::vector< double > cross(const std::vector< double > &a, const std::vector< double > &b){
        /*!
        Compute the cross product between two (3d) vectors.

        :param std::vector< double > a: The first vector
        :param std::vector< double > b: The second vector
        */

        if ((a.size() != b.size()) || (a.size() != 3)){
            std::cout << "Error: check vector dimensions\n";
            assert(1==0);
        }

        std::vector< double > result(3,0);
        result[0] = a[1]*b[2] - a[2]*b[1];
        result[1] = -(a[0]*b[2] - a[2]*b[0]);
        result[2] = a[0]*b[1] - a[1]*b[0];
        return result;
    }

    bool fuzzy_equals(const double a, const double b, const double tolr, const double tola){
        /*!
        Compare two doubles to determine if they are equal.
        */

        double tol = fmin(tolr*fabs(a) + tola, tolr*fabs(b) + tola);
        return fabs(a-b)<tol;
    }


}
