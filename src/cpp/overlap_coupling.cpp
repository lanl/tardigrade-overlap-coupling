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

#include "occonfiguration.h"

#if CONVEXLIB == QUICKHULL
    #define QUICKHULL_IMPLEMENTATION
#elif CONVEXLIB == CONVHULL_3D
    #define CONVHULL_3D_ENABLE
#endif

#include "overlap_coupling.h"

namespace overlap{
    //!===
    //! | Classes
    //!===

    OverlapCoupling::OverlapCoupling(){}

    OverlapCoupling::OverlapCoupling(const vecOfvec &_local_coordinates, const vecOfvec &_gauss_points){
        /*!
        The constructor for the filter.
        */

        local_coordinates = _local_coordinates;
        gauss_points = _gauss_points;
        compute_element_bounds();
    }

    //! > Interface to 3D-quickhull

    vertex_t OverlapCoupling::map_vector_to_quickhull(const std::vector< double > &vector) const{
        /*!
        Map a vector to a vertex which can be read using 3d-quickhull.

        :param std::vector< double > vector: The incoming vector to be mapped.
        */
        vertex_t vertex;
        vertex.x = vector[0];
        vertex.y = vector[1];
        vertex.z = vector[2];
        return vertex;
    }

    std::vector< double > OverlapCoupling::map_quickhull_to_vector(const vertex_t &vertex) const{
        /*!
        Map a quickhull vertex to std::vector

        :param vertex_t vertex: The quickhull vertex type.
        */

        std::vector< double > vector(3,0);
        vector[0] = vertex.x;
        vector[1] = vertex.y;
        vector[2] = vertex.z;
        return vector;
    }

    void OverlapCoupling::map_vectors_to_quickhull(const vecOfvec &vectors, std::vector< vertex_t > &vertices) const{
        /*!
        Map a collection of vectors to a vertex which can be read by 3D-quickhull
        
        :param vecOfvec vectors: a vector of std::vector< double > representing the vector values
        :param std::vector< vertex_t > vertices: The output in a format readable by 3D-quickhull
        */

        vertices.reserve(vectors.size());
        for (unsigned int i=0; i<vectors.size(); i++){
            vertices.push_back(map_vector_to_quickhull(vectors[i]));
        }
        return;
    }

    void OverlapCoupling::map_quickhull_to_vectors(const std::vector< vertex_t > &vertices, vecOfvec &vectors) const{
        /*!
        Map a collection of 3D-quickhull vertices to std::vectors

        :param std::vector< vertex_t > vertices: The incoming quickhull vertices
        :param vecOfvec vectors: The outgoing vector of std::vectors
        */

        vectors.reserve(vertices.size());
        for (unsigned int i=0; i<vertices.size(); i++){
            vectors.push_back(map_quickhull_to_vector(vertices[i]));
        }
    }
#if CONVEXLIB != AKUUKKA
    void OverlapCoupling::extract_mesh_info(const mesh_t &mesh, vecOfvec &normals, vecOfvec &points) const{
#else
    void OverlapCoupling::extract_mesh_info(mesh_t &mesh, vecOfvec &normals, vecOfvec &points) const{
#endif
        /*!
        Extract the required information from the quickhull mesh

        :param mesh_t mesh: The output from the convex hull generating algorithm
        :param vecOfvec normals: The normal vectors of the hull's facets
        :param vecOfvec points: Points on the hull's facets
        */

        #if CONVEXLIB == QUICKHULL
            //!Extract the normal information
            std::vector< vertex_t > _normals(mesh.normals, mesh.normals+mesh.nnormals);

            //!Get a single point of each facet corresponding to the normal.
            std::vector< vertex_t > _points;
            _points.reserve(_normals.size());

            unsigned int index = 0;
            for (unsigned int i=0; i<mesh.nindices/3; i++){
                _points.push_back(mesh.vertices[mesh.indices[index]]);
                index += 3;
            }

            //!Map the convexhull representations to vectors
            map_quickhull_to_vectors(_normals, normals);
            map_quickhull_to_vectors(_points, points);

        #elif CONVEXLIB == CONVHULL_3D
            //!Extract the normal information and assign points
            normals.reserve(mesh.first.size()/3);
            points.reserve(mesh.first.size()/3);
            int i1, i2, i3, index;
            index = 0;
            for (unsigned int i=0; i<mesh.first.size()/3; i++){
                //Get the indices of the vertices which form a face
                i1 = mesh.first[index+0];
                i2 = mesh.first[index+1];
                i3 = mesh.first[index+2];

                //Construct the normal and get a point on the face.
                normals.push_back(normal_from_vertices(mesh.second[i1], mesh.second[i2], mesh.second[i3]));
                points.push_back(map_quickhull_to_vector(mesh.second[i1]));
                index += 3;
            }
        #elif CONVEXLIB == AKUUKKA
            //!Extract the normal information and assign points
            auto indexBuffer = mesh.getIndexBuffer();
            auto vertexBuffer = mesh.getVertexBuffer();

            normals.reserve(indexBuffer.size()/3);
            points.reserve(indexBuffer.size()/3);

            int i1, i2, i3;
            for (unsigned int i=0; i<indexBuffer.size()/3; i++){
                //Get the indices of the vertices which form a face
                i1 = indexBuffer[i*3+0];
                i2 = indexBuffer[i*3+1];
                i3 = indexBuffer[i*3+2];

                //Construct the normal and get a point on the face.
                normals.push_back(normal_from_vertices(vertexBuffer[i1], vertexBuffer[i2], vertexBuffer[i3]));
//                std::cout << "n: ";print_vector(normals[i]);
                points.push_back(map_quickhull_to_vector(vertexBuffer[i1]));
            }
            
        #endif

    }

    void OverlapCoupling::compute_element_bounds(){
        /*!
        Compute the bounds of the element by constructing its convex hull.
        */

        //!Compute the bounding planes
        element_bounds.resize(3);
        compute_node_bounds(local_coordinates, element_planes, element_bounds[0], element_bounds[1], element_bounds[2]);

        return;
    }

    void OverlapCoupling::compute_dns_bounds(const vecOfvec &DNScoordinates){
        /*!
        Compute the bounds of the DNS.

        :param vecOfvec DNScoordinates: A vector of std::vectors of the DNS point coordinates
        */

        //!Compute the bounding planes
        dns_bounds.resize(3);
        compute_node_bounds(DNScoordinates, dns_planes, dns_bounds[0], dns_bounds[1], dns_bounds[2]);
    }

    void OverlapCoupling::compute_node_bounds(const vecOfvec &coordinates, planeMap &planes,
        std::vector< double > &xbnds, std::vector< double > &ybnds, std::vector< double > &zbnds, const double tolr, const double tola){
        /*!
        Compute the bounding planes for the provided coordinates.

        :param vecOfvec coordinates: A vector of std::vectors of the nodal coordinates.
        :param planeMap planes: The resulting planes which bound the coordinates.
        :param std::vector< double > xbnds: The bounds of the nodes in the x direction
        :param std::vector< double > ybnds: The bounds of the nodes in the y direction
        :param std::vector< double > zbnds: The bounds of the nodes in the z direction
        :param double tolr: The relative tolerance
        :param double tola: The absolute tolerance
        */
        //!Map the coordinates to 3D-quickhull vertices
        std::vector< vertex_t > vertices;
        map_vectors_to_quickhull(coordinates, vertices);
        
        //!Construct the mesh
        #if CONVEXLIB == QUICKHULL
            mesh_t mesh = qh_quickhull3d(&vertices[0], vertices.size());
        #elif CONVEXLIB == CONVHULL_3D
            int *faceIndices = NULL;
            int nFaces;
            mesh_t mesh;
            convhull_3d_build(&vertices[0], vertices.size(), &faceIndices, &nFaces);
            std::cout << "nFaces: " << nFaces << "\n";
            mesh.first.reserve(3*nFaces);
            for (int *it = faceIndices; it != faceIndices + 3*nFaces; it++){
                mesh.first.push_back(*it);
            }
            free(faceIndices);
//            mesh.first.assign(faceIndices, faceIndices + 3*nFaces);
            mesh.second = vertices;
            std::cout << "mesh.first.size(): " << mesh.first.size() << "\n";
            std::cout << "mesh.second.size(): " << mesh.second.size() << "\n";
        #elif CONVEXLIB == AKUUKKA
            quickhull::QuickHull<FloatType> qh;
            mesh_t mesh = qh.getConvexHull(&vertices[0].x, vertices.size(), false, false);
        #endif

        //!Extract the relevant information
        vecOfvec normals, points;
        extract_mesh_info(mesh, normals, points);

        //!Form the planes
        planes = compute_unique_planes(normals, points);

        //!Find the bounding box
        xbnds.resize(2); ybnds.resize(2); zbnds.resize(2);
        xbnds[0] = xbnds[1] = planes.begin()->second[0];
        ybnds[0] = ybnds[1] = planes.begin()->second[1];
        zbnds[0] = zbnds[1] = planes.begin()->second[2];
        planeMap::iterator it;
        for (it = planes.begin(); it != planes.end(); it++){
            xbnds[0] = fmin(xbnds[0], it->second[0]);
            xbnds[1] = fmax(xbnds[1], it->second[0]);
            ybnds[0] = fmin(ybnds[0], it->second[1]);
            ybnds[1] = fmax(ybnds[1], it->second[1]);
            zbnds[0] = fmin(zbnds[0], it->second[2]);
            zbnds[1] = fmax(zbnds[1], it->second[2]);
        }
//        std::cout << "xbnds: " << xbnds[0] << " " << xbnds[1] << "\n";
//        std::cout << "ybnds: " << ybnds[0] << " " << ybnds[1] << "\n";
//        std::cout << "zbnds: " << zbnds[0] << " " << zbnds[1] << "\n";

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

    const planeMap* OverlapCoupling::get_element_planes() const{
        /*!
        Extract the element planes

        :param planeMap* planes: A pointer to the element planes.
        */

        return &element_planes;
    }

    const vecOfvec* OverlapCoupling::get_element_bounds() const{
        /*!
        Extract the element bounds

        :param vecOfvec* bounds: A pointer to the element bounds.
        */

        return &element_bounds;
    }

    const planeMap* OverlapCoupling::get_dns_planes() const{
        /*!
        Extract the dns planes

        :param planeMap* planes: A pointer to the dns planes.
        */

        return &dns_planes;
    }

    const vecOfvec* OverlapCoupling::get_dns_bounds() const{
        /*!
        Extract the dns bounds

        :param vecOfvec* bounds: A pointer to the dns planes.
        */

        return &dns_bounds;
    }

//    void OverlapCoupling::construct_element_container(){
//        /*!
//        Construct the element container. Filling it with the gauss points.
//        */
//
//        element_container = 
//    }

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
        vecOfvec global_nodes;
        vecOfvec local_nodes;
        vecOfvec local_gpts;
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

        //Read in the gauss points of the element
        for (unsigned int i=0; i<8; i++){
            std::getline(file, line);
            split_line = split(line, ' ');
            local_gpts.push_back({});
            for (unsigned int j=0; j<3; j++){local_gpts[i].push_back(::atof(split_line[j].c_str()));}
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

        return ParsedData(global_nodes, local_nodes, local_gpts, node_numbers, volumes, densities, coordinates);
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

    std::vector< double > normal_from_vertices(const vertex_t &p1, const vertex_t &p2, const vertex_t &p3){
        /*!
        Compute the normal vector from three vertices.

        :param vertex_t p1: The center vertex of the triangle
        :param vertex_t p2: The first vertex CCW from the center vertex
        :param vertex_t p3: the second vertex CCW from the center vertex
        */

//        std::cout << "p1: ", print_vertex(p1);
//        std::cout << "p2: ", print_vertex(p2);
//        std::cout << "p3: ", print_vertex(p3);

        std::vector< double > v1(3), v2(3), normal;
        v1[0] = p2.x - p1.x;
        v1[1] = p2.y - p1.y;
        v1[2] = p2.z - p1.z;

        v2[0] = p3.x - p1.x;
        v2[1] = p3.y - p1.y;
        v2[2] = p3.z - p1.z;

        //Compute the normal
        normal = cross(v1, v2);

        //Normalize the normal vector
        double mag = sqrt(dot(normal, normal));
        for (unsigned int i=0; i<normal.size(); i++){normal[i]/=mag;}
        return normal;
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

    void print_vertex(const vertex_t &vertex){
        /*!
        Print the value of a vertex to the terminal (debugging tool)

        :param vertex_t vertex: The vertex to print
        */

        printf("%1.6f %1.6f %1.6f", vertex.x, vertex.y, vertex.z);
    }

    void print_vector(const std::vector< FloatType > &vector){
        /*!
        Print the value of a vector to the terminal (debugging tool)

        :param std::vector< FloatType > vector: The vector to print
        */

        printf("%1.6f %1.6f %1.6f", vector[0], vector[1], vector[2]);
    }

    void print_planeMap(const planeMap &planes){
        /*!
        Print the value of a planeMap to the terminal (debugging tool)

        :param planeMap planes: The planemap to print

        */
        planeMap::const_iterator it;
        int padlen = 30;
        std::string str1 = "normals";
        std::string str2 = "points";
        int prel1 = std::ceil(0.5*(padlen - str1.length()));
        int postl1 = std::floor(0.5*(padlen - str1.length()));
        int prel2 = std::ceil(0.5*(padlen - str2.length()));
        int postl2 = std::floor(0.5*(padlen - str2.length()));

        printf("%*s%s%*s|%*s%s%*s\n", prel1, "", "normals", postl1, "", prel2, "", "points", postl2, "");
        for (it = planes.begin(); it != planes.end(); it++){
            printf("%+1.6f %+1.6f %+1.6f | %+1.6f %+1.6f %+1.6f\n",
                   it->first[0], it->first[1], it->first[2],
                   it->second[0], it->second[1], it->second[2]);
        }
    }

    
}
