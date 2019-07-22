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

//Indicate that the library is being compiled
#define OVERLAP_LIBCOMPILE

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
        initialize(_local_coordinates, _gauss_points);
//        local_coordinates = _local_coordinates;
//        gauss_points = _gauss_points;
//        compute_element_bounds();
//        construct_gauss_domains();
    }

    void OverlapCoupling::initialize(const vecOfvec &_local_coordinates, const vecOfvec &_gauss_points){
        /*!
        Initialize the overlap coupling object.

        :param const vecOfvec &_local_coordinates: The local coordinates of the domain's nodes
        :param const vecOfvec &_gauss_points: The local coordinates of the gauss points for the domain.
        */
        local_coordinates = _local_coordinates;
        gauss_points = _gauss_points;
        compute_element_bounds();
        construct_gauss_domains();
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
/*        std::cout << "dns_bounds:\n";
        for (unsigned int i=0; i<dns_bounds.size(); i++){
            std::cout << "    ";
            for (unsigned int j=0; j<dns_bounds[i].size(); j++){
                std::cout << dns_bounds[i][j] << " ";
            }
            std::cout << "\n";
        }
        double min_x, max_x = DNScoordinates[0][0];
        double min_y, max_y = DNScoordinates[0][1];
        double min_z, max_z = DNScoordinates[0][2];

        for (unsigned int i=1; i<DNScoordinates.size(); i++){
            min_x = fmin(min_x, DNScoordinates[i][0]);
            min_y = fmin(min_y, DNScoordinates[i][1]);
            min_z = fmin(min_z, DNScoordinates[i][2]);
            max_x = fmax(max_x, DNScoordinates[i][0]);
            max_y = fmax(max_y, DNScoordinates[i][1]);
            max_z = fmax(max_z, DNScoordinates[i][2]);
        }

        std::cout << "bounds x: " << min_x << " " << max_y << "\n";
        std::cout << "bounds y: " << min_y << " " << max_y << "\n";
        std::cout << "bounds z: " << min_z << " " << max_z << "\n";
*/
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

    void OverlapCoupling::construct_gauss_domains(){
        /*!
        Construct the Gauss domains by using a Voronoi cell representation of the volumes associated with a Gauss point.
        */

        //Map the planes to voro::wall_plane objects
        std::vector< voro::wall_plane > vplanes;
//        std::cout << "element_planes:\n";
//        print_planeMap(element_planes);
        map_planes_to_voro(element_planes, vplanes);

        //Construct the container
        std::vector< unsigned int > gpt_nums(gauss_points.size());
        for (unsigned int i=0; i<gpt_nums.size(); i++){gpt_nums[i] = i;}
//        std::cout << "gauss_points:\n";
//        print_matrix(gauss_points);
/*        std::cout << "  element_bounds:\n";
        for (unsigned int i=0; i<element_bounds.size(); i++){
            std::cout << "    ";
            for (unsigned int j=0; j<element_bounds[i].size(); j++){
                std::cout << element_bounds[i][j] << " ";
            }
            std::cout << "\n";
        }
*/
        voro::container *container = construct_container(gpt_nums, gauss_points, element_bounds, vplanes);

        //Iterate over the gauss points
        voro::voronoicell_neighbor c;
        voro::c_loop_all cl(*container);

        //Loop over the contained points
        std::vector< double > cell_normals;
        std::vector< int > face_vertices;
        std::vector< int > planes;
        std::vector< double > vertices;
        std::vector< double > centroid(3);
        std::vector< double > areas;
        double x, y, z;
        int ifv = 0, index=0;
        vecOfvec normals;
        vecOfvec points;

        gauss_domains.resize(gauss_points.size());

        if (cl.start()) do if (container->compute_cell(c, cl)){
            //Extract the required values from the point
            cl.pos(x, y, z);
            c.normals( cell_normals );
            c.face_vertices( face_vertices );
            c.face_areas( areas );
            c.vertices(x, y, z, vertices);

            //Prepare to populate the required values
            ifv = 0;
            normals.resize(cell_normals.size()/3);
            points.resize(cell_normals.size()/3);
            planes.resize(cell_normals.size()/3);
            //Loop over the domain's faces
            for (unsigned int i=0; i<cell_normals.size()/3; i++){
                //Extract the normal
                normals[i] = {cell_normals[3*i+0],
                              cell_normals[3*i+1],
                              cell_normals[3*i+2]};

                //Find the centroid of each of the face
                find_face_centroid(face_vertices, vertices, ifv, points[i]);

                //Set the face number
                planes[i] = i;

                ifv += face_vertices[ifv]+1;
            }

            //Compute the centroid of the domain
            c.centroid(centroid[0], centroid[1], centroid[2]);
            centroid[0] += x;
            centroid[1] += y;
            centroid[2] += z;

            //Add the point
            gauss_domains[index] = MicroPoint(c.volume(), centroid, planes, 
                                              areas, normals, points);
            index++;

        } while (cl.inc());

        //Free the memory associated with the container
        delete(container);
    }

    void OverlapCoupling::compute_weights(const std::vector< unsigned int > &numbers, const vecOfvec &positions,
                                          std::vector< integrateMap > &points){
        /*!
        Compute the weights of the DNS points for their integration of the gauss domains along with other quantities which 
        will be required.

        :param std::vector< unsigned int > numbers: The number (index) of the DNS point.
        :param vecOfvec positions: The positions of the DNS points in local coordinates.
        :param std::vector< std::vector< MicroPoint > > points: The outgoing weights and other required quantities.
        */

        const MicroPoint* mp;
        voro::container* container;

        //Compute the bounds of the DNS
        compute_dns_bounds(positions);

        //Iterate through the gauss domains
        points.resize(gauss_domains.size());
        for (unsigned int gd=0; gd<gauss_domains.size(); gd++){

            //Construct the Voro++ planes from the gauss domain
            mp = &gauss_domains[gd];
            std::vector< voro::wall_plane > planes;
            map_domain_to_voro(*mp, planes);
            map_planes_to_voro(dns_planes, planes, planes.size());

            //Construct the container
            container = construct_container(numbers, positions, element_bounds, planes);
            
            //Evaluate the point information
            evaluate_container_information(container, points[gd], boundary_node_volumes);

            delete(container);
        }

        for (unsigned int gd=0; gd<gauss_domains.size(); gd++){
//            std::cout << "gd: " << gd << "\n";

            for (auto mpit = points[gd].begin(); mpit != points[gd].end(); mpit++){
                if (boundary_node_volumes.find(mpit->first) != boundary_node_volumes.end()){
                    mpit->second.weight =mpit->second.volume/boundary_node_volumes[mpit->first];
                }
//                std::cout << mpit->first << ": " << mpit->second.weight << "\n";
            }

        }
    }

    const std::vector< MicroPoint >* OverlapCoupling::get_gauss_domains() const{
        /*
        Get a pointer to the gauss domains
        */
        return &gauss_domains;
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

    void MicroPoint::print() const{
        /*!
        Print the contents of the MicroPoint to the terminal (debugging tool)
        */

        std::cout << "MicroPoint:\n";

        std::cout << "  volume: " << volume << "\n";

        std::cout << "  coordinates: "; print_vector(coordinates);

        std::cout << "  planes:";
        for (unsigned int i=0; i<planes.size(); i++){
            printf("%10d", planes[i]);
        }
        std::cout << "\n";

//        std::cout << "  areas: ";
//        for (unsigned int i=0; i<areas.size(); i++){
//            printf("  %1.6f", areas[i]);
//        }
//        std::cout << "\n";
//
//        std::cout << "  normals:\n";
//        for (unsigned int i=0; i<normals.size(); i++){
//            std::cout << "          ";
//            print_vector(normals[i]);
//        }

        std::cout << "  das:\n";
        for (unsigned int i=0; i<das.size(); i++){
            std::cout << "          ";
            print_vector(das[i]);
        }

        std::cout << "  face centroids:\n";
        for (unsigned int i=0; i<face_centroids.size(); i++){
            std::cout << "          ";
            print_vector(face_centroids[i]);
        }
    }

    std::vector< double > MicroPoint::normal(unsigned int i) const{
        /*
        Return the normal vector of the external face

        :param unsigned int i: The index of the face to compute the normal for
        */

        std::vector< double > n = das[i];
        double area = sqrt(dot(n, n));
        for (unsigned int j=0; j<n.size(); j++){n[j]/=area;}
        return n;
    }

    double MicroPoint::area(unsigned int i) const{
        /*
        Compute the area of the external face

        :param unsigned int i: The index of the face to compute the area for
        */

        std::vector< double > n = das[i];
        return sqrt(dot(n, n));
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

        printf("%+1.6f %+1.6f %+1.6f", vertex.x, vertex.y, vertex.z);
    }

    void print_vector(const std::vector< FloatType > &vector){
        /*!
        Print the value of a vector to the terminal (debugging tool)

        :param std::vector< FloatType > vector: The vector to print
        */

        for (unsigned int i=0; i<vector.size(); i++){
            printf("%+1.6f ",vector[i]);
        }
        std::cout << "\n";
    }

    void print_matrix(const std::vector< std::vector< FloatType > > &matrix){
        /*!
        Print the value of a matrix to the terminal (debugging tool)
        */
        for (unsigned int i=0; i<matrix.size(); i++){
            print_vector(matrix[i]);
        }
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

    void add_planes_to_container(std::vector< voro::wall_plane > &planes, voro::container &container){
        /*!
        Add the planes as defined to the voro::container object

        :param std::vector< voro::wall_plane > planes: The planes to be added to the container
        :param voro::container container: The voro++ container object.
        */

        for (unsigned int i = 0; i<planes.size(); i++){
            container.add_wall(planes[i]);
        }
    }

    voro::container* construct_container(const std::vector< unsigned int > &point_numbers, const vecOfvec &point_coords,
                                         const vecOfvec &bounds, std::vector< voro::wall_plane > &planes, double expand){
        /*!
        Returns the pointer to a new voro::container formed by the walls in planes and containing the points in point_coords.

        NOTE: delete this memory after using!

        :param std::vector< unsigned int > point_numbers: The point id numbers
        :param vecOfvec point_coords: The coordinates of the points
        :param vecOfvec bounds: The bounds of the domains
        :param std::vector< voro::wall_plane > planes: The definitions of the bounding planes.
        :param double expand: The amount to expand the bounds by. Ensures that points on the surface aren't id'd as being outside
        */

        //Set up the pre-container class. This will try to optimize the fitting process
        int nx, ny, nz;
        voro::pre_container pcontainer(bounds[0][0]-expand, bounds[0][1]+expand,
                                       bounds[1][0]-expand, bounds[1][1]+expand,
                                       bounds[2][0]-expand, bounds[2][1]+expand,
                                       false, false, false);

        //Add the points to the pre-container
        if (point_numbers.size() != point_coords.size()){
            std::cout << "Error: The point indices and coordinates must have the same length\n";
            assert(1==0);
        }
        for (unsigned int i=0; i<point_numbers.size(); i++){
            pcontainer.put(point_numbers[i], point_coords[i][0], point_coords[i][1], point_coords[i][2]);
        }
        pcontainer.guess_optimal(nx, ny, nz);

        //Set up the container using the guess as to the best options from the pre-container
        voro::container* container = new voro::container(bounds[0][0]-expand, bounds[0][1]+expand,
                                                         bounds[1][0]-expand, bounds[1][1]+expand,
                                                         bounds[2][0]-expand, bounds[2][1]+expand,
                                                         nx, ny, nz, false, false, false, 8);
        pcontainer.setup(*container);

        //Add the additional bounding planes to the container
        std::vector< voro::wall_plane > ps;
        add_planes_to_container(planes, *container);

        return container;
    }

    void evaluate_container_information(voro::container* container, integrateMap &points, std::map< unsigned int, FloatType> &boundary_node_volumes){
        /*!
        Compute required container information (volumes, surface areas, etc.) and return them.

        The planes which describe the container's bounds are desired to have negative id numbers starting 
        at -1 and progressing onwards. It will be assumed that any negative number for a neighboring point 
        is actually a surface on the outside of the domain. The planes will be id'd by adding 1 and taking 
        the negative.

        :param voro::container* container: A pointer to the Voro++ container class to be investigated
        :param std::vector< integrateMap > points: A vector of point information containers.
        :param std::map< unsigned int, FloatType > boundary_node_volumes: The total volume of the cells on gauss domain boundaries
        */

        voro::voronoicell_neighbor c;
        voro::c_loop_all cl(*container);
        
        std::vector< int > neighbors;
        std::vector< double > face_areas;
        std::vector< double > cell_normals;
        std::vector< int > face_orders;
        std::vector< int > face_vertices;
        double x, y, z;
        std::vector< double > vertices;
        std::vector< double > centroid(3);

        std::vector< int > planes; 
        std::vector< double > areas;
        vecOfvec normals;
        vecOfvec face_centroids;
        unsigned int index = 0;
        unsigned int index_order = 0;
        bool is_boundary = false;

        //Define the iterators
        std::vector< int >::iterator viit;
        std::vector< double >::iterator vdit;

//        double vtot = 0;

        //Loop over the contained points
        if (cl.start()) do if (container->compute_cell(c, cl)){
            cl.pos(x, y, z);
            c.neighbors( neighbors );
            c.face_areas( face_areas );
            c.normals( cell_normals );
            c.face_vertices( face_vertices );
            c.vertices(x, y, z, vertices);

//            std::cout << "face_vertices: ";
//            for (unsigned int i=0; i<face_vertices.size(); i++){printf("%4d",face_vertices[i]);}
//            std::cout << "\n";

            planes.resize(0);
            areas.resize(0);
            normals.resize(0);
            face_centroids.resize(0);

            viit = neighbors.begin();
            vdit = face_areas.begin();

            index = 0;
            index_order = 0;
            is_boundary = false;

            while (viit != neighbors.end()){
//                std::cout << "index_order: " << index_order << "\n";
                if (*viit < 0){
                    planes.push_back(-(*viit+1));
                    areas.push_back(*vdit);
                    normals.push_back({cell_normals[index+0],
                                       cell_normals[index+1],
                                       cell_normals[index+2]});
                    find_face_centroid(face_vertices, vertices, index_order, centroid);
                    face_centroids.push_back({centroid[0], centroid[1], centroid[2]});
                    is_boundary = true;

                }
                viit++;
                vdit++;
                index += 3;
                index_order += face_vertices[index_order]+1;
            }
            if (is_boundary){
                auto bit = boundary_node_volumes.find(cl.pid());
                if (bit != boundary_node_volumes.end()){
                    boundary_node_volumes[cl.pid()] += c.volume();
                }
                else{
                    boundary_node_volumes.insert( std::pair< unsigned int, FloatType >(cl.pid(), c.volume()));
                }
            }

            c.centroid(centroid[0], centroid[1], centroid[2]);
            centroid[0] += x;
            centroid[1] += y;
            centroid[2] += z;
//            std::cout << "cl.pid(): " << cl.pid() << "\n";
            points.insert(std::pair< unsigned int, MicroPoint>(cl.pid(), MicroPoint(c.volume(), centroid, planes, areas, normals, face_centroids)));
//            points[cl.pid()].print();
//            vtot += c.volume();

        } while (cl.inc());
//        std::cout << "vtot: " << vtot << "\n";
    }

    void find_face_centroid(const std::vector< int > &face_vertices, const std::vector< double > &vertices, const int &index, std::vector< double > &centroid){
        /*!
        Find the centroid of the given face

        :param std::vector< int > face_vertices: The indices of the vertices for the faces of the cell
        :param std::vector< double > vertices: The coordinates of the vertices for the given cell
        :param int index: The index corresponding to the desired face
        :param std::vector< double > centroid: The centroid of the face
        */

        int k, l, n = face_vertices[index];
        centroid.resize(3);
        centroid[0] = centroid[1] = centroid[2] = 0;
        for (k=0; k<n; k++){
            l = 3*face_vertices[index+k+1];
//            printf("%+1.6f %+1.6f %+1.6f\n", vertices[l+0], vertices[l+1], vertices[l+2]);
            centroid[0] += vertices[l+0]/n;
            centroid[1] += vertices[l+1]/n;
            centroid[2] += vertices[l+2]/n;
        }

/*        vecOfvec edges;
        edges.reserve(n);
        int l0, l1;
        for (k=1; k<n; k++){
            l0 = 3*face_vertices[index+1];
            l1 = 3*face_vertices[index+k+1];
            edges.push_back({vertices[l1+0]-vertices[l0+0],
                             vertices[l1+1]-vertices[l0+1],
                             vertices[l1+2]-vertices[l0+2]});
        }

        std::cout << "distances:\n";
        std::vector< double > normal = cross(edges[1], edges[0]);
        for (unsigned int i=0; i<edges.size(); i++){
            l1 = 3*face_vertices[index+k+1];
            std::cout << dot(normal, edges[i]) << "\n";
        }

        std::cout << "centroid distance: ";
        l0 = 3*face_vertices[index+1];
        std::cout << dot(normal, {centroid[0]-vertices[l0+0], centroid[1]-vertices[l0+1], centroid[2]-vertices[l0+2]});
        std::cout << "\n";
*/
    }
    
    void map_planes_to_voro(const planeMap &planes, std::vector< voro::wall_plane > &vplanes, int j){
        /*!
        Map planes to voro::wall_plane objects.

        :param planeMap planes: The planeMap object which stores the plane information (normal, point on plane)
        :param std::vector< voro::wall_plane > vplanes: The vector of Voro++ plane objects.
        :param int j: The initial value of the planes to use
        */

        //Map the planes to voro::wall_plane objects
        vplanes.reserve(planes.size());
        planeMap::const_iterator it;
        double distance;

//        std::cout << "map planes to voro:\n";

        for (it=planes.begin(); it!=planes.end(); it++){
            distance = overlap::dot(it->first, it->second);
//            std::cout << "  normal: " << it->first[0] << " " << it->first[1] << " " << it->first[2] << "\n";
//            std::cout << "  face centroid: " << it->second[0] << " " << it->second[1] << " " << it->second[2] << "\n";
//            std::cout << "  distance: " << distance << "\n";
            vplanes.push_back(voro::wall_plane(it->first[0], it->first[1], it->first[2], distance, -(j+1)));
            j++;
        }
    }

    void map_domain_to_voro(const MicroPoint &domain, std::vector< voro::wall_plane > &vplanes){
        /*!
        Map a domain (here represented by a MicroPoint) to a std::vector of Voro++ wall_plane objects
    
        :param MicroPoint domain: The MicroPoint representation of a domain.
        :param std::vector< voro::wall_plane > vplanes: A std::vector of voro::wall_plane objects which can be added to a voro::container
        */

        unsigned int j=1, n = domain.das.size();
        double distance;
        std::vector< double > normal(3);
        vplanes.reserve(n);

//        std::cout << "map_domain_to_voro:\n";

        for (unsigned int i=0; i<n; i++){
            normal = domain.normal(i);
            distance = overlap::dot(normal, domain.face_centroids[i]);
//            std::cout << "  normal: " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
//            std::cout << "  face centroid: " << domain.face_centroids[i][0] << " " << domain.face_centroids[i][1] << " " << domain.face_centroids[i][2] << "\n";
//            std::cout << "  distance: " << distance << "\n";
            vplanes.push_back(voro::wall_plane(normal[0], normal[1], normal[2], distance, -j));
            j++;
        }
    }

    void apply_nansons_relation(const std::vector< double > &N, const double &JdA, const vecOfvec &Finv, std::vector< double > &nda){
        /*!
        Apply Nanson's relation to the differential areas.

        Assumes 3D

        :param std::vector< double > N: The normal in the reference configuration
        :param double JdA: The Jacobian of transformation times the differential area in the reference configuration
        :param vecOfvec Finv: The inverse of the deformation gradient which transforms between the reference and current configurations
        :param std::vector< double > nda: The normal vector in the current configuration weighted by the area
        */

        if ((N.size() != 3) || (Finv.size() != 3)){
            std::cout << "Error: This implementation only works for 3D\n";
            assert(1==0);
        }

        for (unsigned int i=0; i<Finv.size(); i++){
            if (Finv[i].size() != 3){
                std::cout << "Error: This implementation only works for 3D\n";
                assert(1==0);
            }
        }

        nda.resize(3);

        for (unsigned int i=0; i<3; i++){
            nda[i] = 0;
            for (unsigned int j=0; j<3; j++){
                nda[i] += JdA*N[j]*Finv[j][i];
            }
        }
    }

    void perform_volume_integration( const std::map< unsigned int, double > &values, const std::vector< integrateMap > &weights, std::vector< double > &result){
        /*
        Perform the volume integration of a scalar value. Returns the integrated value at the gauss points.

        :param map::< unsigned int, double > values: The map object which defines the values to be integrated at each of the micro-nodes
        :param std::vector< integrateMap > weights: The weights of each of the nodes in true space for each gauss point.
        :param std::vector< double > result: The integrated result.
        */

        //Set up an iterator for the value map
        std::map<unsigned int, double >::const_iterator itv;
        integrateMap::const_iterator itiM;

        //Initialize the result vector
        result = std::vector< double >(weights.size(), 0);
        
        //Loop over the gauss points
        for (unsigned int gp=0; gp<weights.size(); gp++){

            for (itiM=weights[gp].begin(); itiM!=weights[gp].end(); itiM++){
                
                //Find the value of the function at the node
                itv = values.find(itiM->first);
                if (itv == values.end()){
                    std::cout << "Error: node " << itiM->first << " not found in values\n";
                    assert(1==0);
                }

                //Add the term to the integral.
                result[gp] += itv->second*itiM->second.volume;
            }
        }
    }

    void perform_position_weighted_volume_integration( const std::map< unsigned int, double > &values, const std::vector< integrateMap > &weights, std::vector< std::vector< double > > &result){
        /*
        Perform the volume integration of a scalar value multiplied by the position of the point.
        Returns the integrated value at the gauss points.

        :param map::< unsigned int, double > values: The map object which defines the values to be integrated at each of the micro-nodes
        :param std::vector< integrateMap > weights: The weights of each of the nodes in true space for each gauss point.
        :param std::vector< std::vector< double > > result: The integrated result.
        */

        //Set up an iterator for the value map
        std::map<unsigned int, double >::const_iterator itv;
        integrateMap::const_iterator itiM;

        //Initialize the result vector
        result = std::vector< std::vector< double > >(weights.size());
        //Loop over the gauss points
        for (unsigned int gp=0; gp<weights.size(); gp++){
            result[gp] = std::vector<double>(weights[gp].begin()->second.coordinates.size(), 0);

            for (itiM=weights[gp].begin(); itiM!=weights[gp].end(); itiM++){
                
                //Find the value of the function at the node
                itv = values.find(itiM->first);
                if (itv == values.end()){
                    std::cout << "Error: node " << itiM->first << " not found in values\n";
                    assert(1==0);
                }

                //Add the term to the integral.
                for (unsigned int i=0; i<result[gp].size(); i++){
                    result[gp][i] += itv->second*itiM->second.coordinates[i]*itiM->second.volume;
                }
            }
        }
    }

    void perform_volume_integration( const std::map< unsigned int, std::vector< double > > &values, const std::vector< integrateMap > &weights, std::vector< std::vector< double > > &result){
        /*
        Perform the volume integration of a vector value. Returns the integrated value at the gauss points.

        :param map::< unsigned int, std::vector< double > > values: The map object which defines the values to be integrated at each of the micro-nodes
        :param std::vector< integrateMap > weights: The weights of each of the nodes in true space for each gauss point.
        :param std::vector< std::vector< double > > result: The integrated result.
        */

        //Set up an iterator for the value map
        std::map<unsigned int, std::vector< double >>::const_iterator itv;
        integrateMap::const_iterator itiM;

        //Initialize the result vector
        result = std::vector< std::vector< double > >(weights.size());
        
        //Loop over the gauss points
        for (unsigned int gp=0; gp<weights.size(); gp++){

//            result[gp].resize(values.begin()->second.size());
            result[gp] = std::vector<double>(values.begin()->second.size(), 0);

            for (itiM=weights[gp].begin(); itiM!=weights[gp].end(); itiM++){
                
                //Find the value of the function at the node
                itv = values.find(itiM->first);
                if (itv == values.end()){
                    std::cout << "Error: node " << itiM->first << " not found in values\n";
                    assert(1==0);
                }

                //Add the terms to the integral.
                if (itv->second.size() != result[gp].size()){
                    std::cout << "Error: result and value must have the same size\n";
                }
                for (unsigned int i=0; i<result[gp].size(); i++){
                    result[gp][i] += itv->second[i]*itiM->second.volume;
                }
            }
        }
    }

    void perform_surface_integration( const std::map< unsigned int, double > &values, const std::vector< integrateMap > &weights, std::vector< std::map< unsigned int, double > > &result){
        /*!
        Perform the surface integration of a scalar value. Returns the integrated value at each of the surface of the gauss points

        :param std::map< unsigned int, double > &values: The values at each of the micro-points
        :param std::vector< integrateMap > weights: The weights of each of the nodes in true space for each gauss point
        :param std::vector< std::map< unsigned int, double > > result: The result of the integration over each of the faces
        */

        //Set up an iterator for the value map
        std::map< unsigned int, double >::const_iterator itv;
        integrateMap::const_iterator itiM;
        std::map< unsigned int, double >::iterator itr;

        //Initialize the result vector
        result.resize(weights.size());

        //Loop over the gauss points
        for (unsigned int gp=0; gp<weights.size(); gp++){
        
            //Loop over the micro-node weights    
            for (itiM=weights[gp].begin(); itiM!=weights[gp].end(); itiM++){

                //Find the value of the function at the node
                itv = values.find(itiM->first);
                if (itv == values.end()){
                    std::cout << "Error: node " << itiM->first << " not found in values\n";
                    assert(1==0);
                }

                //Loop over the planes
                for (unsigned int j=0; j<itiM->second.planes.size(); j++){
                    itr = result[gp].find(itiM->second.planes[j]);

                    //Insert the plane if new
                    if (itr == result[gp].end()){
                        result[gp].insert( std::pair< unsigned int, double >(itiM->second.planes[j], itv->second*itiM->second.area(j)));
                    }
                    //Add to the plane if it exists already
                    else{
                        itr->second += itv->second*itiM->second.area(j);
                    }
                }
            }
        }
    }

    void perform_surface_integration( const std::map< unsigned int, std::vector< double > > &values, const std::vector< integrateMap > &weights, std::vector< std::map< unsigned int, std::vector< double > > > &result){
        /*!
        Perform the surface integration of a scalar value. Returns the integrated value at each of the surface of the gauss points

        :param std::map< unsigned int, double > &values: The values at each of the micro-points
        :param std::vector< integrateMap > weights: The weights of each of the nodes in true space for each gauss point
        :param std::vector< std::map< unsigned int, double > > result: The result of the integration over each of the faces
        */

        //Set up an iterator for the value map
        std::map< unsigned int, std::vector< double > >::const_iterator itv;
        integrateMap::const_iterator itiM;
        std::map< unsigned int, std::vector< double > >::iterator itr;

        //Initialize the result vector
        result.resize(weights.size());

        //Loop over the gauss points
        for (unsigned int gp=0; gp<weights.size(); gp++){
        
            //Loop over the micro-node weights    
            for (itiM=weights[gp].begin(); itiM!=weights[gp].end(); itiM++){

                //Find the value of the function at the node
                itv = values.find(itiM->first);
                if (itv == values.end()){
                    std::cout << "Error: node " << itiM->first << " not found in values\n";
                    assert(1==0);
                }

                //Loop over the planes
                for (unsigned int j=0; j<itiM->second.planes.size(); j++){
                    itr = result[gp].find(itiM->second.planes[j]);

                    //Compute the value of the vector
                    std::vector< double > vec = itv->second;
                    for (unsigned int i=0; i<vec.size(); i++){vec[i] *= itiM->second.area(j);}

                    //Insert the plane if new
                    if (itr == result[gp].end()){
                        result[gp].insert( std::pair< unsigned int, std::vector< double > >(itiM->second.planes[j], vec));
                    }
                    //Add to the plane if it exists already
                    else{
                        for (unsigned int i=0; i<vec.size(); i++){ itr->second[i] += vec[i];}
                    }
                }
            }
        }
    }

    void construct_triplet_list(const std::map< unsigned int, unsigned int >* macro_node_to_col_map,
                                const std::map< unsigned int, unsigned int >* dns_node_to_row_map,
                                const std::vector< unsigned int > &macro_node_ids,
                                const std::vector< FloatType > &cg,
                                const vecOfvec &psis,
                                const integrateMap &dns_weights,
                                const std::map< unsigned int, unsigned int>* micro_node_elcount,
                                bool share_ghost_free_boundary_nodes,
                                bool macro_elem_is_ghost,
                                unsigned int num_micro_free,
                                std::vector< T > &tripletList,
                                unsigned int num_macro_dof,
                                unsigned int num_micro_dof){
        /*!
        Add the contributions of the nodes contained within a quadrature domain to the shape-function matrix
        triplet list.
        
        :param std::map< unsigned int, unsigned int >* macro_node_to_col_map: A pointer to the map from the macro node id numbers to the location in the shape-function matrix. The value will be scaled by 12 since there are assumed to be 12 DOF for the macro nodes (i.e. 3D isothermal behavior)
        :param std::map< unsigned int, unsigned int >* micro_node_to_row_map: A pointer to the map from the micro node id numbers to the location in the shape-function matrix. The value will be scaled by 3 since there are assumed to be 3 DOF for the micro nodes (i.e. 3D isothermal behavior)
        :param std::vector< unsigned int > macro_node_ids: The id numbers of the macro-scale nodes
        :param std::vector< FloatType > cg: The center of gravity of the macro node.
        :param vecOfvec psis: The shape function values for each of the nodes at the cg.
        :param integrateMap dns_weights: The weights and locations of the micro-nodes in the macro element. All positions, volumes, and das should be in true space (i.e. not in local/master coordinates)
        :param std::map< unsigned int, unsigned int>* micro_node_elcount: The number of elements a given node is inside. Assumed 1 if not in map.
        :param bool share_ghost_free_boundary_nodes: Boolean indicating if micro-nodes on a boundary between a free and ghost macro element should be shared
        :param bool macro_elem_is_ghost: Boolean indicating if the macro-element is a ghost element or not
        :param std::vector< T > tripletList: The list of triplets used to construct the sparse matrix representation of the shape-function matrix for this quadrature domain
        :param unsigned int num_macro_dof: The number of macro-scale DOF per node. Note that the code is only set up to deal with 12 but the values are in place for future expansion.
        :param unsigned int num_micro_dof: The number of micro-scale DOF per node. Note that the code is only set up to deal with 3 but the values are in place for future expansion.
        */

        //Initialize the variables and iterators
        
        std::map< unsigned int, unsigned int >::const_iterator it;
        integrateMap::const_iterator iMit;
        std::vector< FloatType > xi(3,0);
        unsigned int row0, col0;
        double psi_n;
        double weight;

        //Reserve the memory required for the tripletList
        tripletList.reserve(tripletList.size() + 8*num_macro_dof*dns_weights.size());

//        std::cout << "cg: ";
//        for (unsigned int i=0; i<3; i++){
//            std::cout << cg[i] << " ";
//        }
//        std::cout << "\n";

        //Loop over the macro nodes
        for (unsigned int n=0; n<psis.size(); n++){
            //Set the initial index of the column
            it = macro_node_to_col_map->find(macro_node_ids[n]);
            if (it != macro_node_to_col_map->end()){
                col0 = num_macro_dof*it->second;
                psi_n = psis[n][0];
            }
            else{
                std::cout << "Error: Macro node not found in macro_node_to_col map\n";
                assert(1==0);
            }

            //Loop over the micro-nodes and add the contributions to the shape-function matrix
            for (iMit=dns_weights.begin(); iMit!=dns_weights.end(); iMit++){

                //Set the initial index of the row
                it = dns_node_to_row_map->find(iMit->first);
                if (it != dns_node_to_row_map->end()){
                    row0 = num_micro_dof*it->second;
                }
                else{
                    std::cout << "Error: Micro node not found in micro_node_to_row map\n";
                    assert(1==0);
                }

                //Check if macro elements id'd as free can contain free micro elements and if the current node fails that test
                if ((!share_ghost_free_boundary_nodes) && (!macro_elem_is_ghost) && (it->second < num_micro_free)){
                    continue;
                }

                //Set the weight
                weight = iMit->second.weight;
                auto mnelit = micro_node_elcount->find(iMit->first);
                if (mnelit != micro_node_elcount->end()){
                    weight /= mnelit->second;
                }

                //Compute the xi vector
                for (unsigned int j=0; j<3; j++){
                    xi[j] = iMit->second.coordinates[j] - cg[j];
                }

/*                if (it->second == 450){
                    std::cout << "row0: " << row0 << "\n";
                    std::cout << "col0: " << col0 << "\n";
                    std::cout << "share_ghost_free_boundary_nodes: " << share_ghost_free_boundary_nodes << "\n";
                    std::cout << "macro_elem_is_ghost: " << macro_elem_is_ghost << "\n";
                    std::cout << "num_micro_free: " << num_micro_free << "\n";
                    std::cout << "iMit->second.coordinates: ";
                    for (unsigned int i=0; i<iMit->second.coordinates.size(); i++){
                        std::cout << iMit->second.coordinates[i] << " ";
                    }
	            std::cout << "\n";
                    std::cout << "xi: ";
                    for (unsigned int i=0; i<xi.size(); i++){
                        std::cout << xi[i] << " ";
                    }
	            std::cout << "\n";
                    std::cout << "psi_n: " << psi_n << "\n";
                    std::cout << "dns_weight: " << iMit->second.weight << "\n";
                    std::cout << "weight: " << weight << "\n";
                    if (mnelit != micro_node_elcount->end()){
                        std::cout << "number of neighboring elements: " << mnelit->second << "\n";
                    }
                    
                }
*/
                //Add the values to the matrix
                tripletList.push_back(T(row0 + 0, col0 +  0, weight*psi_n));
                tripletList.push_back(T(row0 + 1, col0 +  1, weight*psi_n));
                tripletList.push_back(T(row0 + 2, col0 +  2, weight*psi_n));
                tripletList.push_back(T(row0 + 0, col0 +  3, weight*psi_n*xi[0]));
                tripletList.push_back(T(row0 + 1, col0 +  4, weight*psi_n*xi[1]));
                tripletList.push_back(T(row0 + 2, col0 +  5, weight*psi_n*xi[2]));
                tripletList.push_back(T(row0 + 1, col0 +  6, weight*psi_n*xi[2]));
                tripletList.push_back(T(row0 + 0, col0 +  7, weight*psi_n*xi[2]));
                tripletList.push_back(T(row0 + 0, col0 +  8, weight*psi_n*xi[1]));
                tripletList.push_back(T(row0 + 2, col0 +  9, weight*psi_n*xi[1]));
                tripletList.push_back(T(row0 + 2, col0 + 10, weight*psi_n*xi[0]));
                tripletList.push_back(T(row0 + 1, col0 + 11, weight*psi_n*xi[0]));
            }
        }
        return;
    }

    SpMat form_sparsematrix(const std::vector< T > &tripletList, unsigned int nrows, unsigned int ncols, const bool &ignore_dup){
        /*!
        Form a sparse matrix from a triplet list

        :param std::vector< T > tripletlist: A list of triplets to use to construct the matrix
        :param bool ignore_dup: Flag to indicate whether duplicates should be ignored
        :param SpMat result: The resulting sparse matrix
        */

        //Initialize the matrix
        SpMat result(nrows, ncols);

        //Construct the matrix following the logic on ignoring duplicates
        if (ignore_dup){
            result.setFromTriplets(tripletList.begin(), tripletList.end(), [] (const double&, const double &b) { return b; });
        }
        else{
            result.setFromTriplets(tripletList.begin(), tripletList.end());
        }
        result.makeCompressed();
        return result;
    }

    SpMat extract_block(const SpMat &A, unsigned int start_row, unsigned int start_col, unsigned int nrows, unsigned int ncols){
        /*!
        Extract a block from a larger sparse matrix
        */

        SpMat block = A.block(start_row, start_col, nrows, ncols);
        block.makeCompressed();
        return block;
    }

//    QRsolver form_solver(SpMat &A){
//        /*!
//        Form the QR solver.
//        */
//        QRsolver Asolver;
//        A.makeCompressed();
//        Asolver.compute(A);
//        if (Asolver.info()!=Eigen::Success){
//            std::cout << "Error: Least squares solution to update degrees of freedom failed.\n";
//            assert(1==0);
//        }
//        return Asolver;
//    }

    void solve_for_projector(const SpMat &A, const SpMat &B, SpMat &X){
        /*!
        Solve for the projector by solving the sparse matrix linear algebra problem using QR decomposition.
    
        AX = B

        :param SpMat A: The left-hand-side matrix
        :param SpMat B: The right-hand-side matrix
        :param SpMat X: The solution matrix.
        */

        //Define the solver
        QRsolver solver;

        //Perform the decomposition
        solver.compute(A);
        if( solver.info() != Eigen::Success){
            std::cout << "Error: Least squares solution to solving for the projector failed\n";
            assert(1==0);
        }
        X = solver.solve(B);
    }

    Projector::Projector(){
        /*!
        Constructor
        */
        tripletList.resize(0);
    }

    Projector::Projector(unsigned int _num_macro_dof, unsigned int _num_micro_dof,
                         unsigned int _num_macro_ghost, unsigned int _num_macro_free, 
                         unsigned int _num_micro_ghost, unsigned int _num_micro_free){
        /*!
        Constructor
        */

        num_macro_dof = _num_macro_dof;
        num_micro_dof = _num_micro_dof;
        num_macro_ghost = _num_macro_ghost;
        num_macro_free = _num_macro_free;
        num_micro_ghost = _num_micro_ghost;
        num_micro_free = _num_micro_free;
        tripletList.resize(0);
    }

    Projector::Projector(const Projector &p){
        /*!
        The copy constructor of Projector.
        */

        num_macro_dof = p.num_macro_dof;
        num_micro_dof = p.num_micro_dof;
        num_macro_ghost = p.num_macro_ghost;
        num_macro_free = p.num_macro_free;
        num_micro_ghost = p.num_micro_ghost;
        num_micro_free = p.num_micro_free;
        tripletList = p.tripletList;
    }

    void Projector::initialize(unsigned int _num_macro_dof, unsigned int _num_micro_dof,
                               unsigned int _num_macro_ghost, unsigned int _num_macro_free,
                               unsigned int _num_micro_ghost, unsigned int _num_micro_free){
        /*!
        Initialize the projector object

        :param unsigned int _num_macro_dof: The number of macro-scale degrees of freedom
        :param unsigned int _num_micro_dof: The number of micro-scale degrees of freedom
        :param unsigned int _num_macro_ghost: The number of macro-scale ghost nodes in the overlap domain
        :param unsigned int _num_macro_free: The number of macro-scale free nodes in the overlap domain
        :param unsigned int _num_micro_ghost: The number of micro-scale ghost nodes in the overlap domain
        :param unsigned int _num_micro_free: The number of micro-scale free nodes in the overlap domain

        :return: Initialize the projector object
 
        :rtype: void
        */

        num_macro_dof = _num_macro_dof;
        num_micro_dof = _num_micro_dof;
        num_macro_ghost = _num_macro_ghost;
        num_macro_free = _num_macro_free;
        num_micro_ghost = _num_micro_ghost;
        num_micro_free = _num_micro_free;
        tripletList.resize(0);
    }

    void Projector::add_shapefunction_terms(const std::map< unsigned int, unsigned int >* macro_node_to_col_map,
                                            const std::map< unsigned int, unsigned int >* micro_node_to_row_map,
                                            const std::vector< unsigned int > &macro_node_ids,
                                            const std::vector< FloatType > &cg,
                                            const vecOfvec &psis,
                                            const integrateMap &dns_weights,
                                            const std::map< unsigned int, unsigned int>* micro_node_elcount,
                                            bool share_ghost_free_boundary_nodes,
                                            bool macro_elem_is_ghost,
                                            unsigned int num_micro_free){//,
//                                            unsigned int num_macro_dof,
//                                            unsigned int num_micro_dof){
        /*!
        Add the contributions of the nodes contained within a quadrature domain to the shape-function matrix
        triplet list.

        Note: This principally exists to isolate the version of Eigen used here.
        
        :param std::map< unsigned int, unsigned int >* macro_node_to_col_map: A pointer to the map from the macro node id numbers to the location in the shape-function matrix. The value will be scaled by 12 since there are assumed to be 12 DOF for the macro nodes (i.e. 3D isothermal behavior)
        :param std::map< unsigned int, unsigned int >* micro_node_to_row_map: A pointer to the map from the micro node id numbers to the location in the shape-function matrix. The value will be scaled by 3 since there are assumed to be 3 DOF for the micro nodes (i.e. 3D isothermal behavior)
        :param std::vector< unsigned int > macro_node_ids: The id numbers of the macro-scale nodes
        :param std::vector< FloatType > cg: The center of gravity of the macro node.
        :param vecOfvec psis: The shape function values for each of the nodes at the cg.
        :param integrateMap dns_weights: The weights and locations of the micro-nodes in the macro element. All positions, volumes, and das should be in true space (i.e. not in local/master coordinates)
        :param std::map< unsigned int, unsigned int>* micro_node_elcount: The number of elements a given node is inside. Assumed 1 if not in map.
        :param bool share_ghost_free_boundary_nodes: Boolean indicating if micro-nodes on a boundary between a free and ghost macro element should be shared
        :param bool macro_elem_is_ghost: Boolean indicating if the macro-element is a ghost element or not
        :param unsigned int num_macro_dof: The number of macro-scale DOF per node. Note that the code is only set up to deal with 12 but the values are in place for future expansion.
        :param unsigned int num_micro_dof: The number of micro-scale DOF per node. Note that the code is only set up to deal with 3 but the values are in place for future expansion.
        */

        construct_triplet_list(macro_node_to_col_map, micro_node_to_row_map, macro_node_ids,
                               cg, psis, dns_weights,
                               micro_node_elcount, share_ghost_free_boundary_nodes,
                               macro_elem_is_ghost, num_micro_free,
                               tripletList, num_macro_dof, num_micro_dof);
    }

    void Projector::form_shapefunction_matrix(unsigned int nrows, unsigned int ncols){
        /*!
        Form the shapefunction matrix
        :param unsigned int nrows: The number of rows in the matrix
        :param unsigned int ncols: The number of columns in the matrix
        */

        shapefunction = SpMat(nrows, ncols);
        shapefunction.setFromTriplets(tripletList.begin(), tripletList.end());
        shapefunction.makeCompressed();
    }

    int Projector::form_BDhQsolver(){
        /*!
        Form the solver for the BDhQ projector.
        */

        std::cout << "  Performing NQDh QR decomposition\n";
        BDhQsolver.compute(shapefunction.block( 0, num_macro_dof*num_macro_free,
                                                num_micro_dof*num_micro_free, num_macro_dof*num_macro_ghost));
        if (BDhQsolver.info() != Eigen::Success){
            return 1;
        }
        return 0;
    }

    int Projector::form_NQDh_PR_transpose_solver(){
        /*!
        Form the solver for problems of the form:

        BDhQ.transpose x = b

        where

        NQDh P = Q R

        NQDh = Q R P.transpose

        NQDh BDhQ = I

        Q R P.transpose BDhQ = I

        BDhQ = P R.inverse Q.transpose
        BDhQ.transpose = Q R.inverse_transpose P.transpose

        So for problems of the form
        x = BDhQ.transpose b

        x = Q R.inverse_transpose P.transpose b
        P R.transpose Q.transpose x = b
        
        So to solve for x we first must solve the linear equation
        P R.transpose xp = b

        And then:
        x = Q xp
        */

        std::cout << "  Performing NQDh transpose QR decomposition\n";

        unsigned int rows = BDhQsolver.matrixR().rows();
        unsigned int cols = BDhQsolver.matrixR().cols();

        if (rows < cols){
            //Return error code because there are more macro scale dof than micro scale
            return 1;
        }

        //Extract the transpose of R
//        SpMat matrixR_transpose = BDhQsolver.matrixR().block(0, 0,
//                                                             num_macro_dof*num_macro_ghost, num_macro_dof*num_macro_ghost).transpose();
        SpMat matrixR_transpose = BDhQsolver.matrixR().block(0, 0, cols, cols);

        //Perform the permutation
        SpMat PR_transpose = BDhQsolver.colsPermutation()*matrixR_transpose;
        NQDh_PR_transpose_solver.compute(PR_transpose);
        if (NQDh_PR_transpose_solver.info() != Eigen::Success){
            //Return error code because the decomposition failed
            return 2;
        }
        return 0;
    }

    void Projector::solve_BDhQ(const std::vector< FloatType > &Qvec, std::vector< FloatType > &Dhvec) const{
        /*!
        Solve an equation of the type:

        Dh = BDhQ Q

        :param const std::vector< FloatType > Q: The incoming Q vector.
        :param std::vector< FloatType > Dh: The outgoing Dh vector.
        */

        Eigen::Map<const EigVec> Q(Qvec.data(), Qvec.size(), 1);

        if (Dhvec.size() != num_macro_dof*num_macro_ghost){
            Dhvec = std::vector< double >(num_macro_dof*num_macro_ghost, 0);
        }

        Eigen::Map< EigVec > Dh(Dhvec.data(), Dhvec.size(), 1);
        Dh = BDhQsolver.solve(Q);
//        Dhvec.resize(Dh.size());
//        Eigen::VectorXd::Map(&Dhvec[0], Dh.size()) = Dh;
    }

    void Projector::solve_BDhQtranspose(const std::vector< FloatType > &bvec, std::vector< FloatType > &xvec) const{
        /*!
        Solve an equation of the type

        x = BDhQ.transpose b

        :param const std::vector< FloatType > bvec: The right-hand-side vector
        :param std::vector< FloatType > xvec: The solution vector
        */

        Eigen::Map<const EigVec> b(bvec.data(), bvec.size(), 1);
        xvec.resize(num_micro_dof*num_micro_free);
        Eigen::Map<EigVec> x(xvec.data(), xvec.size(), 1);
//        if (xvec.size() != num_micro_dof*num_micro_free){
//            xvec = std::vector< double >(num_micro_dof*num_micro_free, 0);
//        }
//        Eigen::Map< EigVec > x(xvec.data(), xvec.size(), 1);
//
//        EigVec xp = NQDh_PR_transpose_solver.solve(b);
//        x = BDhQsolver.matrixQ()*xp;

        EigVec xp = NQDh_PR_transpose_solver.solve(b);
        Eigen::MatrixXd matrixQ = BDhQsolver.matrixQ();
        x = matrixQ.block(0, 0, num_micro_dof*num_micro_free, num_macro_dof*num_macro_ghost)*xp;
    }

    void Projector::solve_BQhDtranspose(const std::vector< FloatType > &bvec, std::vector< FloatType > &xvec) const{
        /*!
        Solve an equation of the type

        x = BQhD b

        where
        BQhD = NQhd + NQhDh BDhD

        and we assume BDhD = 0

        :param const std::vector< FloatType > bvec: The right-hand-side vector
        :param std::vector< FloatType > xvec: The solution vector.

        */
        Eigen::Map<const Eigen::MatrixXd> b(bvec.data(), 1, bvec.size());
        xvec.resize(num_macro_dof*num_macro_free);
        Eigen::Map<Eigen::MatrixXd> x(xvec.data(), 1, xvec.size());

//        assert(-1==2);

        x = b*shapefunction.block(num_micro_dof*num_micro_free, 0,
                                  num_micro_dof*num_micro_ghost, num_macro_dof*num_macro_free);
    }

    void Projector::solve_BQhQtranspose(const std::vector< FloatType > &bvec, std::vector< FloatType > &xvec) const{
        /*!
        Solve an equation of the type

        x = BQhQtranspose b

        where
        BQhQ = NQhD BDhQ -> BDhQtranspose = BDhQtranspose NQhDtranspose

        -> x = BDhQtranspose NQhDtranspose b
        */

        Eigen::Map< const Eigen::MatrixXd > b(bvec.data(), 1, bvec.size());

        std::vector< FloatType > bstarvec(num_macro_dof*num_macro_ghost, 0);
        Eigen::Map< Eigen::MatrixXd > bstar(bstarvec.data(), 1, bstarvec.size());

        //Transform using the right-hand-size using the shapefunction matrix
        bstar = b*shapefunction.block(num_micro_dof*num_micro_free, num_macro_dof*num_macro_free,
                                      num_micro_dof*num_micro_ghost, num_macro_dof*num_macro_ghost);

        //Solve for the transformed matrix
        solve_BDhQtranspose(bstarvec, xvec);
    }

    int Projector::_run_tests(bool solve_for_projectors){
        /*!
        Run some informal verification tests.
        */

        EigVec Dhtmp = EigVec::Zero(num_macro_dof*num_macro_ghost, 1);
        for (unsigned int i=0; i<num_macro_ghost; i++){
            Dhtmp(num_macro_dof*i + 0) =  0.32;
            Dhtmp(num_macro_dof*i + 1) =  1.00;
            Dhtmp(num_macro_dof*i + 2) = -3.42;
        }

        EigVec Dtmp = EigVec::Zero(num_macro_dof*num_macro_free, 1);
        for (unsigned int i=0; i<num_macro_free; i++){
            Dtmp(num_macro_dof*i + 0) =  0.32;
            Dtmp(num_macro_dof*i + 1) =  1.00;
            Dtmp(num_macro_dof*i + 2) = -3.42;
        }

        //Test if the macro-scale values are interpolated correctly
        SpMat NQDh = shapefunction.block( 0, num_macro_dof*num_macro_free,
                                          num_micro_dof*num_micro_free, num_macro_dof*num_macro_ghost);
        EigVec Qtmp = NQDh*Dhtmp;

        bool xtest, ytest, ztest;

        for (unsigned int i=0; i<num_micro_free; i++){

            xtest = fuzzy_equals(Qtmp(num_micro_dof*i + 0),  0.32);
            ytest = fuzzy_equals(Qtmp(num_micro_dof*i + 1),  1.00);
            ztest = fuzzy_equals(Qtmp(num_micro_dof*i + 2), -3.42);
            if (!(xtest && ytest && ztest)){
                std::cout << "i: " << i << "\n";
                std::cout << "num_micro_free: " << num_micro_free << "\n";
                std::cout << "Qtmp(" << num_micro_dof*i + 0 << "): " << Qtmp(num_micro_dof*i + 0) << "\n";
                std::cout << "Qtmp(" << num_micro_dof*i + 1 << "): " << Qtmp(num_micro_dof*i + 1) << "\n";
                std::cout << "Qtmp(" << num_micro_dof*i + 2 << "): " << Qtmp(num_micro_dof*i + 2) << "\n";
                std::cout << "Test 1 failed: Micro-dof not expected value\n";
                return 1;
                assert(1==0);
            }
        }

        //Test if the solver solves correctly
        EigVec Dhans = BDhQsolver.solve(Qtmp);

        if (!fuzzy_equals((Dhans - Dhtmp).norm(), 0)){
            std::cout << "Test 2 failed\n";
            return 2;
            assert(2==0);
        }

        //Test of the wrapper for the solver solves correctly
        std::vector< double > Dhvec;
        std::vector< double > Qtmpvec(Qtmp.rows(), 0);
        for (unsigned int i=0; i<Qtmp.rows(); i++){
            Qtmpvec[i] = Qtmp[i];
        }
        solve_BDhQ(Qtmpvec, Dhvec);
        for (unsigned int i=0; i<Dhvec.size(); i++){
            if (!fuzzy_equals(Dhvec[i], Dhans[i])){
                std::cout << "Test 3 failed\n";
                return 3;
            }
        }

        //Make sure that there are no non-zero terms in NQD
        SpMat NQD   = shapefunction.block( 0, 0, num_micro_dof*num_micro_free, num_macro_dof*num_macro_free);
        if (!fuzzy_equals(NQD.norm(), 0)){
            std::cout << "Test 4 failed\n";
            return 4;
            assert(3==0);
        }

        //Make sure that the interpolation of NQhD and NQhDh is carried out correctly
        SpMat NQhD  = shapefunction.block( num_micro_dof*num_micro_free, 0, num_micro_dof*num_micro_ghost,  num_macro_dof*num_macro_free);
        SpMat NQhDh = shapefunction.block( num_micro_dof*num_micro_free, num_macro_dof*num_macro_free, num_micro_dof*num_micro_ghost, num_macro_dof*num_macro_ghost);
        EigVec Qhtmp = NQhD*Dtmp + NQhDh*Dhtmp;

        std::vector< double > sum_shape_fxn(12, 0);
        unsigned int dofnum = 450;
        EigVec NQhDrow0 = NQhD.row(num_micro_dof*dofnum);
        std::cout << "NQhD.row(" << num_micro_dof*dofnum << "):\n";
        for (unsigned int i=0; i<num_macro_free; i++){
            for (unsigned int j=0; j<12; j++){
                std::cout << NQhDrow0[12*i + j] << " ";
                sum_shape_fxn[j] += NQhDrow0[12*i + j];
            }
            std::cout << "\n";
        }
        std::cout << "sum of values from NQhD: ";
        for (unsigned int j=0; j<12; j++){
            std::cout << sum_shape_fxn[j] << " ";
        }
        std::cout << "\n";

        std::cout << "NQhDh.row(" << num_micro_dof*dofnum << "):\n";
        EigVec NQhDhrow0 = NQhDh.row(num_micro_dof*dofnum);
        std::vector< double > sum_shape_fxn2(12, 0);
        for (unsigned int i=0; i<num_macro_ghost; i++){
            for (unsigned int j=0; j<12; j++){
                std::cout << NQhDhrow0[12*i + j] << " ";
                sum_shape_fxn2[j] += NQhDhrow0[12*i + j];
            }
           std::cout << "\n";
        }
        std::cout << "\n";

        std::cout << "sum of values from NQhDh: ";
        for (unsigned int j=0; j<12; j++){
            std::cout << sum_shape_fxn2[j] << " ";
        }
        std::cout << "\n";

        std::cout << "sum of values from NQhD and NQhDh: ";
        for (unsigned int j=0; j<12; j++){
            std::cout << sum_shape_fxn[j] + sum_shape_fxn2[j] << " ";
        }
        std::cout << "\n";

        std::cout << "D:\n";
        for (unsigned int i=0; i<num_macro_free; i++){
            for (unsigned int j=0; j<12; j++){
                std::cout << Dtmp[12*i + j] << " ";
            }
            std::cout << "\n";
        }

        for (unsigned int i=0; i<num_micro_ghost; i++){
            xtest = fuzzy_equals(Qhtmp(num_micro_dof*i + 0),  0.32);
            ytest = fuzzy_equals(Qhtmp(num_micro_dof*i + 1),  1.00);
            ztest = fuzzy_equals(Qhtmp(num_micro_dof*i + 2), -3.42);
            if (!(xtest && ytest && ztest)){
                std::cout << "i: " << i << "\n";
                std::cout << "num_micro_ghost: " << num_micro_ghost << "\n";
                std::cout << "Qtmp(" << num_micro_dof*i + 0 << "): " << Qhtmp(num_micro_dof*i + 0) << "\n";
                std::cout << "Qtmp(" << num_micro_dof*i + 1 << "): " << Qhtmp(num_micro_dof*i + 1) << "\n";
                std::cout << "Qtmp(" << num_micro_dof*i + 2 << "): " << Qhtmp(num_micro_dof*i + 2) << "\n";
                std::cout << "Test 5 failed: Micro-dof not expected value\n";
                return 5;
            }
        }

        if (solve_for_projectors){
            SpMat NQDh_transpose = NQDh.transpose();

            std::cout << "Null rows/columns of NQDh: size: " << NQDh.rows() << ", " << NQDh.cols() << "\n";
            for (int k=0; k<NQDh.rows(); ++k){
                if (NQDh.row(k).norm() < 1e-8){
                    std::cout << "row: " << k << "\n";
                }
            }

            for (int k=0; k<NQDh.cols(); ++k){
                if (NQDh.col(k).norm() < 1e-8){
                    std::cout << "col: " << k << "\n";
                }
            }

            std::cout << "Null rows/columns of NQDh_transpose: size: " << NQDh_transpose.rows() << ", " << NQDh_transpose.cols() << "\n";
            for (int k=0; k<NQDh_transpose.rows(); ++k){
                if (NQDh_transpose.row(k).norm() < 1e-8){
                    std::cout << "row: " << k << "\n";
                }
            }

            for (int k=0; k<NQDh_transpose.cols(); ++k){
                if (NQDh_transpose.col(k).norm() < 1e-8){
                    std::cout << "col: " << k << "\n";
                }
            }

            Dhans = NQDh_transpose*Qtmp;
            std::cout << "Dhans:\n";
            for (unsigned int i=0; i<num_macro_ghost; i++){
                for (unsigned int j=0; j<12; j++){
                    std::cout << Dhans[12*i + j] << " ";
                }
                std::cout << "\n";
            }

            //Make sure that the BDhQtranspose solver is working as expected
            std::vector< double > Qvec;
            solve_BDhQtranspose(Dhvec, Qvec);
            std::cout << "Qvec size: " << Qvec.size() << "\n";
//            for (unsigned int i=0; i < Qvec.size(); i++){
//                std::cout << Qvec[i] << "\n";
//            }

            if (Qvec.size() != num_micro_dof*num_micro_free){
                std::cout << "Test 6 failed: BDhQtranspose solver returned a vector of improper size.";
                return 6;
            }

            std::vector< double > Qhvec(num_micro_dof*num_micro_ghost, 0);
            for (unsigned int i=0; i<num_micro_ghost; i++){
                Qhvec[i+0] =  1.2;
                Qhvec[i+1] =  2.3;
                Qhvec[i+2] = -3.4;
            }


            std::vector< double > Dvec;
            solve_BQhDtranspose(Qhvec, Dvec);
            if (Dvec.size() != num_macro_dof*num_macro_free){
                std::cout << "Test 7 failed: BQhDtranspose solver returned a vector of improper size.";
            }

            solve_BQhQtranspose(Qhvec, Qvec);
            if (Qvec.size() != num_micro_dof*num_micro_free){
                std::cout << "Test 8 failed: BQhQtranspose solver returned a vector of improper size.";
            }

            assert(-23==-24);

//            EigVec QprojDh_result = NQDh_PR_transpose_solver.solve(Dhans);
//            std::cout << "Q->Dh:\n" << QprojDh_result << "\n";
//            std::cout << "Qtmp_result shape: " << QprojDh_result.rows() << ", " << QprojDh_result.cols() << "\n";
//            std::cout << "NQDH shape: " << NQDh.rows() << ", " << NQDh.cols() << "\n";
//            std::cout << "R shape: " << BDhQsolver.matrixR().rows() << ", " << BDhQsolver.matrixR().cols() << "\n";
//            std::cout << "Q shape: " << BDhQsolver.matrixQ().rows() << ", " << BDhQsolver.matrixQ().cols() << "\n";
//            Eigen::MatrixXd matrixQ = BDhQsolver.matrixQ();
//            Eigen::MatrixXd testing_things = matrixQ.block(0, 0,
//                                                 matrixQ.rows(),
//                                                 BDhQsolver.matrixR().cols())*QprojDh_result;
//            std::cout << "testing_things:\n" << testing_things << "\n"; 

//            Eigen::MatrixXd testing_things = BDhQsolver.matrixQ()*NQDh_PR_transpose_solver.solve(Dhans);
            //std::cout << "Qtmp_result shape: " << Qtmp_result.rows() << ", " << Qtmp_result.cols() << "\n";
//            std::cout << "It worked, but why can't I save it?\n";
//            assert(-1==1);

//            if ((Qtmp - Qtmp_result).norm() > 1e-8){
//                std::cout << "Qtmp_result:\n" << Qtmp_result << "\n";
//                mooseWarning("Test 5 failed: BDhQ_transpose solver returning unexpected values");
//            }
        }

        return 0;
        std::cout << "All tests passed\n";
        assert(-1==1);

    }

    void Projector::project_dof(const std::vector< double > &Dvec, const std::vector< double > &Qvec,
                                std::vector< double > &Dhvec, std::vector< double > &Qhvec) const{
        /*!
        Perform the projection operation projecting from free to ghost degrees of freedom

        :param std::vector< double > Dvec: The macro-scale (micromorphic) free degrees of freedom
        :param std::vector< double > Qvec: The micro-scale (DNS) free degrees of freedom
        :param std::vector< double > Dhvec: The macro-scale (micromorphic) ghost degrees of freedom
        :param std::vector< double > Qhvec: The micro-scale (DNS) free degrees of freedom
        */

        if (Dhvec.size() != num_macro_dof*num_macro_ghost){
            Dhvec = std::vector< double >(num_macro_dof*num_macro_ghost, 0);
        }
        if (Qhvec.size() != num_micro_dof*num_micro_ghost){
	    Qhvec = std::vector< double >(num_micro_dof*num_micro_ghost, 0);
        }

        //Form the Eigen maps
        Eigen::Map<const EigVec> D(Dvec.data(), Dvec.size(), 1);
        Eigen::Map<const EigVec> Q(Qvec.data(), Qvec.size(), 1);
        Eigen::Map<EigVec> Dh(Dhvec.data(), num_macro_dof*num_macro_ghost, 1);
        Eigen::Map<EigVec> Qh(Qhvec.data(), num_micro_dof*num_micro_ghost, 1);

        //Solve for Dh
        Dh = BDhQsolver.solve(Q);

        //Solve for Qh (NQhD * D + NQhDh * Dh
        Qh = shapefunction.block( num_micro_dof*num_micro_free, 0, num_micro_dof*num_micro_ghost, num_macro_dof*num_macro_free) * D
       +shapefunction.block( num_micro_dof*num_micro_free, num_macro_dof*num_macro_free, num_micro_dof*num_micro_ghost, num_macro_dof*num_macro_ghost) * Dh;
        return;
        
    }

//    MicromorphicFilter::MicromorphicFilter(elib::Element* _element, bool _shared_dof_material){
//        /*!
//        Initialize the micromorphic filter
//
//        :param elib::Element &element: The underlying element over which the filter is computed.
//        */
//
//        element = _element;
//        shared_dof_material = _shared_dof_material;
//
//        elib::vecOfvec gauss_points(element->qrule.size());
//        for (unsigned int gp=0; gp<element->qrule.size(); gp++){
//            gauss_points[gp] = element->qrule[gp].first;
//        }
//
//        material_overlap.initialize(_element->local_node_coordinates, gauss_points); 
//    }

    MicromorphicFilter::MicromorphicFilter(const std::string &element_type, const elib::vecOfvec &nodes,
		                           const elib::quadrature_rule &qrule, bool _shared_dof_material){
        /*!
	 * Initialize the micromorphic filter
	 *
	 * :param const std::string &element_type: The name of the underlying element type
	 * :param const elib::vecOfvec &nodes: The nodes which make up the element type (i.e. their coordinates)
	 * :param const elib::quadrature_rule &qrule: The quadrature rule of the element.
	 * :param bool _shared_dof_material: A flag which indicates if the material points should be used to 
	 *                                   determine the macro-scale nodal coordinate information.
	 */

	 element = elib::build_element_from_string(element_type, nodes, qrule);
	 shared_dof_material = _shared_dof_material;

	 elib::vecOfvec gauss_points(element->qrule.size());
	 for (unsigned int gp=0; gp<element->qrule.size(); gp++){
             gauss_points[gp] = element->qrule[gp].first;
	 }
	 
	 material_overlap.initialize(element->local_node_coordinates, gauss_points);
    }

    bool MicromorphicFilter::add_micro_dof_point(const unsigned int &id, const elib::vec &coordinates){
        /*!
        Check if a dof point is inside the filter and add it if it is.

        :param unsigned int &id: The micro point id number
        :param elib::vec &coordinates: The global coordinates of the micro dof point.
        */

        if (element->bounding_box_contains_point(coordinates)){

            elib::vec xi;
            element->compute_local_coordinates(coordinates, xi);
            if (element->local_point_inside(xi)){
                dof_id_numbers.push_back(id);
                micro_dof_local_coordinates.push_back(xi);
                return true;
            }
        }
        return false;
    }

    bool MicromorphicFilter::add_micro_material_point(const unsigned int &id, const elib::vec &coordinates){
        /*!
        * Check if a material point is inside the filter and add it if it is.
        * 
        * :param unsigned int &id: The micro point id number
        * :param elib::vec &coordinates: The global coordinates of the micro material point.
        */

        if (element->bounding_box_contains_point(coordinates)){

            elib::vec xi;
            element->compute_local_coordinates(coordinates, xi);
            if (element->local_point_inside(xi)){
                material_id_numbers.push_back(id);
                micro_material_local_coordinates.push_back(xi);
                return true;
            }
        }
        return false;
    }

    int MicromorphicFilter::construct_integrators(){
        /*!
        Construct the filter's integrators.
        */

        construct_material_point_integrator();
        return 0;
    }

    int MicromorphicFilter::construct_material_point_integrator(){
        /*!
        Construct the integrator for the material points.

        */

        material_overlap.compute_weights(material_id_numbers, micro_material_local_coordinates, material_weights);
        return 0;
    }

    int MicromorphicFilter::compute_mass_properties(const std::map< unsigned int, double> &micro_density){
        /*!
        Compute the mass properties of the homogenized DNS.

        :param std::map< unsigned int , elib::vec> &density:
        */

        compute_volume();
        compute_density(micro_density);
        compute_centers_of_mass(micro_density);
        return 0;
    }

    int MicromorphicFilter::compute_volume(){
        /*!
        Compute the volume of each of the filter's gauss domains
        */

        integrateMap::const_iterator itiM;

        volume = std::vector< double >(material_weights.size(), 0);
        for (unsigned int gp=0; gp<material_weights.size(); gp++){

            for (itiM=material_weights[gp].begin(); itiM!=material_weights[gp].end(); itiM++){

                //Add the term to the integral.
                volume[gp] += itiM->second.volume;
            }
        }
        return 0;
    }

    int MicromorphicFilter::compute_density(const std::map< unsigned int, double> &micro_density){
        /*!
        Compute the density

        :param std::map< unsigned int , elib::vec> &density:
        */

        perform_volume_integration(micro_density, material_weights, density);

        for (unsigned int gp=0; gp<density.size(); gp++){
            density[gp] /= volume[gp];
        }
        return 0;
    }

    int MicromorphicFilter::compute_centers_of_mass(const std::map< unsigned int, double> &micro_density){
        /*!
        Compute the centers of mass of the material.

        :param std::map< unsigned int , elib::vec> &density:
        */

        perform_position_weighted_volume_integration(micro_density, material_weights, center_of_mass);
        
        for (unsigned int gp=0; gp<density.size(); gp++){
            for (unsigned int i=0; i<center_of_mass[gp].size(); i++){
                center_of_mass[gp][i] /= volume[gp]*density[gp];
            }
        }
        return 0;
    }

    int MicromorphicFilter::print(){
        /*!
	 * Print out the information contained in the filter
	 */

	std::cout << "DOF id numbers:\n";
        elib::print(dof_id_numbers);
	std::cout << "Material Point id numbers:\n";
	elib::print(material_id_numbers);
	std::cout << "DOF local coordinates:\n";
	elib::print(micro_dof_local_coordinates);
	std::cout << "Material Point local coordinates:\n";
	elib::print(micro_material_local_coordinates);

	return 0;
    }
}
