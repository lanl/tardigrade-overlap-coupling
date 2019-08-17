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

    void OverlapCoupling::map_vectors_to_quickhull(const std::map< unsigned int, std::vector< FloatType > > &vectors, std::vector< vertex_t > &vertices) const{
        /*!
         * Map a collection of vectors to a vertex which can be read by 3D-quickhull
         * 
         * :param std::map< unsigned int, std::vector< FloatType > > &vectors: a map of std::vector< double > representing the vector values
         * :param std::vector< vertex_t > vertices: The output in a format readable by 3D-quickhull
         */

        vertices.reserve(vectors.size());
        for (auto vector=vectors.begin(); vector!=vectors.end(); vector++){
            vertices.push_back(map_vector_to_quickhull(vector->second));
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
    void OverlapCoupling::extract_mesh_info(mesh_t &mesh, vecOfvec &normals, vecOfvec &points) const {//, vecOfvec &facepoints) const{
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
                std::vector< FloatType > normal;
                int nfv_result = normal_from_vertices(mesh.second[i1], mesh.second[i2], mesh.secont[i3], normal);
                if (!(nfv_result > 0)){
                    normals.push_back(normal);//normal_from_vertices(mesh.second[i1], mesh.second[i2], mesh.second[i3]));
                    points.push_back(map_quickhull_to_vector(mesh.second[i1]));
                }
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
                std::vector< FloatType > normal;
                int nfv_result = normal_from_vertices(vertexBuffer[i1], vertexBuffer[i2], vertexBuffer[i3], normal);

                if (!(nfv_result > 0)){
//                    normals.push_back(normal_from_vertices(vertexBuffer[i1], vertexBuffer[i2], vertexBuffer[i3]));
//                std::cout << "n: ";print_vector(normals[i]);
                    normals.push_back(normal);
                    points.push_back(map_quickhull_to_vector(vertexBuffer[i1]));
                    //facepoints.push_back(map_quickhull_to_vector(vertexBuffer[i1]));
                    //facepoints.push_back(map_quickhull_to_vector(vertexBuffer[i2]));
                    //facepoints.push_back(map_quickhull_to_vector(vertexBuffer[i3]));
                }
            }
            
        #endif

    }

    void OverlapCoupling::compute_element_bounds(){
        /*!
         * Compute the bounds of the element by constructing its convex hull.
         */

        //!Compute the bounding planes
        element_bounds.resize(3);
        std::map< unsigned int, std::vector< FloatType > > local_coordinate_map;
        for (unsigned int i=0; i<local_coordinates.size(); i++){local_coordinate_map.emplace(i, local_coordinates[i]);}

        compute_node_bounds(local_coordinate_map, element_planes, element_bounds[0], element_bounds[1], element_bounds[2]);

        return;
    }

    void OverlapCoupling::compute_dns_bounds(const std::map< unsigned int, std::vector< FloatType > > &DNScoordinates, bool use_dns_bounds){
        /*!
         * Compute the bounds of the DNS.
         * 
         * :param std::map< unsigned int, std::vector< FloatType > > &DNScoordinates: The map of std::vectors of the DNS point coordinates
         */

        //!Compute the bounding planes
        if (use_dns_bounds){
            dns_bounds.resize(3);
            compute_node_bounds(DNScoordinates, dns_planes, dns_bounds[0], dns_bounds[1], dns_bounds[2]);
        }
        else{
            dns_planes = element_planes;
            dns_bounds = element_bounds;
        }
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

    void OverlapCoupling::compute_node_bounds(const std::map< unsigned int, std::vector< FloatType > > &coordinates, planeMap &planes,
        std::vector< double > &xbnds, std::vector< double > &ybnds, std::vector< double > &zbnds, const double tolr, const double tola){
        /*!
        Compute the bounding planes for the provided coordinates.

        :param std::map< unsigned int, std::vector< FloatType > > &coordinates: A map of std::vectors of the nodal coordinates.
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
//        std::cout << "coordinates:\n[";
//        for (unsigned int i=0; i<coordinates.size(); i++){
//            std::cout << "[" << coordinates[i][0] << ", " << coordinates[i][1] << ", " << coordinates[i][2] << "],\n";
//        }
//        std::cout << "]\n";
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
        vecOfvec normals, points;//, face_points;
        extract_mesh_info(mesh, normals, points);//, face_points);

/*        std::cout << "normals:\n[";
        for (unsigned int i=0; i<normals.size(); i++){
            std::cout << "[" << normals[i][0] << ", " << normals[i][1] << ", " << normals[i][2] << "],\n";
        }
        std::cout << "]\n";
        std::cout << "points:\n[";
        for (unsigned int i=0; i<points.size(); i++){
            std::cout << "[" << points[i][0] << ", " << points[i][1] << ", " << points[i][2] << "],\n";
        }
        std::cout << "]\n";
        std::cout << "face_points:\n[";
        for (unsigned int i=0; i<face_points.size(); i++){
            std::cout << "[" << face_points[i][0] << ", " << face_points[i][1] << ", " << face_points[i][2] << "],\n";
        }
        std::cout << "]\n";
*/
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
        map_planes_to_voro(element_planes, vplanes);

        //Add the planes to the external surface id vector
        external_face_ids.resize(gauss_points.size());
        for (unsigned int gp=0; gp<gauss_points.size(); gp++){
            external_face_ids[gp].resize(element_planes.size());
            for (unsigned int i=0; i<element_planes.size(); i++){
                external_face_ids[gp][i] = i;
            }
        }

        //Construct the container
//        std::vector< unsigned int > gpt_nums(gauss_points.size());
        std::map< unsigned int, std::vector< FloatType > > gpt_map;
        for (unsigned int i=0; i<gauss_points.size(); i++){gpt_map.emplace(i, gauss_points[i]);}
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
        voro::container *container = construct_container(gpt_map, element_bounds, vplanes);

        //Iterate over the gauss points
        voro::voronoicell_neighbor c;
        voro::c_loop_all cl(*container);

        //Loop over the contained points
        std::vector< int > neighbors;
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
        domain_vertices.resize(gauss_points.size());
        vertex_planes.resize(gauss_points.size());

        unsigned int plane_num = 0;

        if (cl.start()) do if (container->compute_cell(c, cl)){
            //Extract the required values from the point
//            std::cout << "cl.pid(): " << cl.pid() << "\n";
            index = cl.pid();
//            std::cout << "id: " << index << "\n";
            cl.pos(x, y, z);
//            std::cout << "pos: " << x << ", " << y << ", " << z << "\n";
            c.normals( cell_normals );
            c.neighbors( neighbors );
            c.face_vertices( face_vertices );
            c.face_areas( areas );
            c.vertices(x, y, z, vertices);

            //Set the number of vertices for the current gauss domain
            vertex_planes[index].resize(vertices.size()/3);
            for (unsigned int i=0; i<vertex_planes[index].size(); i++){vertex_planes[index][i].reserve(3);}

            //Prepare to populate the required values
            ifv = 0;
            normals.resize(cell_normals.size()/3);
            points.resize(cell_normals.size()/3);
            planes.resize(cell_normals.size()/3);
            //Loop over the domain's faces
            for (unsigned int i=0; i<cell_normals.size()/3; i++){
                if (neighbors[i] < 0){
                    plane_num = - (neighbors[i] + 1);
                }
                else{
                    plane_num = neighbors[i] + element_planes.size();
                }
//                std::cout << "plane_num: " << plane_num << "\n";
                //Extract the normal
                normals[i] = {cell_normals[3*i+0],
                              cell_normals[3*i+1],
                              cell_normals[3*i+2]};

//                std::cout << "normal: "; elib::print(normals[i]);

                //Find the centroid of each of the faces
                find_face_centroid(face_vertices, vertices, ifv, points[i]);

                //Set the face number
                planes[i] = plane_num;

                //Add the plane to the vertices which are included
                for (int j=0; j<face_vertices[ifv]; j++){vertex_planes[index][face_vertices[ifv+1+j]].push_back(plane_num);}

                plane_num++;

                ifv += face_vertices[ifv]+1;
            }

            //Compute the centroid of the domain
            c.centroid(centroid[0], centroid[1], centroid[2]);
            centroid[0] += x;
            centroid[1] += y;
            centroid[2] += z;

            //Add the point
            gauss_domains[index] = MicroPoint(c.volume(), centroid, {x, y, z}, planes, 
                                              areas, normals, points);

            //Add the coordinates of the domain vertices
            domain_vertices[index].resize(vertices.size()/3);
            for (unsigned int i=0; i<vertices.size()/3; i++){
                domain_vertices[index][i] = {vertices[3*i+0],
                                             vertices[3*i+1],
                                             vertices[3*i+2]};
            }
//            std::cout << "gauss domain " << index << "\n"; gauss_domains[index].print();
//            assert(1==0);
//            index++;

        } while (cl.inc());
//        assert(1==0);

        //Free the memory associated with the container
        delete(container);

    }

    void OverlapCoupling::compute_weights(const std::map< unsigned int, std::vector< FloatType > > &positions,
                                          std::vector< integrateMap > &points, bool use_dns_bounds){
        /*!
         * Compute the weights of the DNS points for their integration of the gauss domains along with other quantities which 
         * will be required.
         * 
         * NOTE: All quantities will be in the coordinate system provided to the filter. This means that if you use,
         *      in the finite element context, the local coordinates of the element, all of the properties will be 
         *      in that reference frame. Mapping these quantities to the global coordinate frame will be required 
         *      to perform integration.
         * 
         * :param std::map< unsigned int, std::vector< FloatType > > positions: Map from the DNS point number to it's local position.
         * :param std::vector< std::vector< MicroPoint > > points: The outgoing weights and other required quantities.
         * :param bool use_dns_bounds: Flag on whether the DNS bounds should be used. When used in the full overlap 
         *     context this may not be desired.
        */

        const MicroPoint* mp;
        voro::container* container;

        std::map< int, std::pair< std::vector< FloatType >, std::vector< FloatType > > > bounding_faces;

        //Compute the bounds of the DNS
        compute_dns_bounds(positions, use_dns_bounds);

        //Iterate through the gauss domains
        points.resize(gauss_domains.size());
        for (unsigned int gd=0; gd<gauss_domains.size(); gd++){

            //Construct the Voro++ planes from the gauss domain
            mp = &gauss_domains[gd];
            std::vector< voro::wall_plane > planes;
//            std::cout << "gd: " << gd << "\n";
//            (*mp).print();
            map_domain_to_voro(*mp, planes);
            bounding_faces.clear();
            for (unsigned int i=0; i<(*mp).planes.size(); i++){
                bounding_faces.emplace(-((*mp).planes[i]+1),
                                       std::pair< std::vector< FloatType >, std::vector< FloatType > >((*mp).normal(i), (*mp).face_centroids[i]));
            }

            if (use_dns_bounds){
                int bni=0;
                for (auto dpit=dns_planes.begin(); dpit!=dns_planes.end(); dpit++){
                    bounding_faces.emplace(-(planes.size() + bni), std::pair< std::vector< FloatType > , std::vector< FloatType > >(dpit->first, dpit->second));
                    external_face_ids[gd].push_back(planes.size() + bni);
                    bni++;
                }
                
                map_planes_to_voro(dns_planes, planes, planes.size());
            }

//            for (auto tmpit=bounding_faces.begin(); tmpit!=bounding_faces.end(); tmpit++){
//                std::cout << tmpit->first << ": normal: "; elib::print(tmpit->second.first);
//                std::cout << "   : point: "; elib::print(tmpit->second.second);
//            }

            //Construct the container
            container = construct_container(positions, element_bounds, planes);
            
            //Evaluate the point information
            evaluate_container_information(positions, bounding_faces, container, points[gd], boundary_node_volumes);
//            assert(1==0);

            if (points[gd].size() == 0 ){
                gauss_domains[gd].print();
                print_planeMap(dns_planes);
                assert(-23==-22);
            }

            delete(container);
        }

//        assert(1==0);

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

    const std::vector< vecOfvec >* OverlapCoupling::get_domain_vertices() const{
        /*
        Get a pointer to the coordinates vertices which comprise the current gauss domain
        */
        return &domain_vertices;
    }

    const std::vector< std::vector< std::vector< unsigned int > > >* OverlapCoupling::get_vertex_planes() const{
        /*
        Get a pointer to the indices of the planes which go through the domain vertices
        */
        return &vertex_planes;
    }

    const std::vector< std::vector< unsigned int > >* OverlapCoupling::get_external_face_ids() const{
        /*!
         * Get a pointer to the external face ids vector
         */
        return &external_face_ids;
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

    void OverlapCoupling::print_element() const{
        /*!
         * Print a selection of properties of the element
         */

        std::cout << "OverlapCoupling Object\n";
        std::cout << "element_bounds:\n";
        elib::print(element_bounds);
        std::cout << "element planes:\n";
        print_planeMap(element_planes);
        std::cout << "Gauss domains\n";
        for (unsigned int i=0; i<gauss_domains.size(); i++){
            gauss_domains[i].print();
        }
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

        std::cout << "  particle coordinates: "; print_vector(particle_coordinates);

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

    bool fuzzy_equals(const std::vector< double > &a, const std::vector< double > &b, const double tolr, const double tola){
        /*!
        Compare two vectors of doubles to determine if they are equal.
        */

        if (a.size() != b.size()){
            return false;
        }
        else{
            for (unsigned int i=0; i<a.size(); i++){
                bool tmp = fuzzy_equals(a[i], b[i], tolr, tola);
                if (!tmp){
                    return false;
                }
            }
        }
        return true;
    }

    int normal_from_vertices(const vertex_t &p1, const vertex_t &p2, const vertex_t &p3, std::vector< double > &normal, double tolr, double tola){
        /*!
         * Compute the normal vector from three vertices. If the resulting area of the two vectors is smaller than the tolerance
         * then the return value will be greater than zero.
         * 
         * :param vertex_t p1: The center vertex of the triangle
         * :param vertex_t p2: The first vertex CCW from the center vertex
         * :param vertex_t p3: the second vertex CCW from the center vertex
         * :param std::vector< double > &normal: The normal vector
         * :param double tol: The tolerance.
         */

//        std::cout << "p1: ", print_vertex(p1);
//        std::cout << "p2: ", print_vertex(p2);
//        std::cout << "p3: ", print_vertex(p3);

        std::vector< double > v1(3), v2(3);
        v1[0] = p2.x - p1.x;
        v1[1] = p2.y - p1.y;
        v1[2] = p2.z - p1.z;

        v2[0] = p3.x - p1.x;
        v2[1] = p3.y - p1.y;
        v2[2] = p3.z - p1.z;

        //Normalize v1 and v2
        double v1_mag = sqrt(dot(v1, v1));
        double v2_mag = sqrt(dot(v2, v2));
        for (unsigned int i=0; i<v1.size(); i++){v1[i]/=v1_mag;}
        for (unsigned int i=0; i<v2.size(); i++){v2[i]/=v2_mag;}

        //Compute the normal
        normal = cross(v1, v2);

        //Normalize the normal vector
        double narea = dot(normal, normal);
        double mag = sqrt(narea);
        for (unsigned int i=0; i<normal.size(); i++){normal[i]/=mag;}

        //Check for numerical issues
        double tol = 1*tolr + tola;
        if (narea < tol ){
            //The area of the rhombus between the two vectors is very small.
            //They are likely either nearly in the same direction or in 
            //opposite directions and the result is suspect.
            return 1;
        }
        
        return 0;
    }

    bool compare_vector_directions(const std::vector< double > &v1, const std::vector< double > &v2, const double tolr, const double tola, bool opposite_is_unique){
        /*!
         * Compare vectors to determine if they are in the same direction.
         * 
         * :param std::vector< double > v1: The first vector
         * :param std::vector< double > v2: The second vector
         * :param const double tolr: The relative tolerance
         * :param const double tola: The absolute tolerance
         * :param bool opposite_is_unique: If true, vectors in opposite directions are different.
         * 
         * The dot product is computed and checked if it is nearly equal to 1.
         */

        double factor = sqrt(dot(v1, v1)*dot(v2, v2));
        double result = dot(v1, v2)/factor;
        if (!opposite_is_unique){
//            std::cout << "result: " << result << "\n";
            result = std::abs(result);
        }
        return fuzzy_equals(result, 1, tolr, tola);
    }

    void print_vertex(const vertex_t &vertex){
        /*!
         * Print the value of a vertex to the terminal (debugging tool)
         * 
         * :param vertex_t vertex: The vertex to print
         */

        printf("%+1.6f %+1.6f %+1.6f", vertex.x, vertex.y, vertex.z);
    }

    void print_vector(const std::vector< FloatType > &vector){
        /*!
         * Print the value of a vector to the terminal (debugging tool)
         * 
         * :param std::vector< FloatType > vector: The vector to print
         */

        for (unsigned int i=0; i<vector.size(); i++){
            printf("%+1.6f ",vector[i]);
        }
        std::cout << "\n";
    }

    void print_matrix(const std::vector< std::vector< FloatType > > &matrix){
        /*!
         * Print the value of a matrix to the terminal (debugging tool)
         * 
         * const vecOfvec &matrix: The matrix to be printed
         */
        for (unsigned int i=0; i<matrix.size(); i++){
            print_vector(matrix[i]);
        }
    }

    void print_planeMap(const planeMap &planes){
        /*!
         * Print the value of a planeMap to the terminal (debugging tool)
         * 
         * :param planeMap planes: The planemap to print
         * 
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

    void print_coordinateMap(const std::map< unsigned int, std::vector< FloatType > > &coordinates){
        /*!
         * Print the coordinate map (id, coordinates) to the terminal
         * 
         * :param const std::map< unsigned int, std::vector< FloatType > > &coordinates: The coordinate map
         */
        
        std::map< unsigned int, std::vector< FloatType > >::const_iterator it;
        int padlen = 30;
        std::string str1 = "id";
        std::string str2 = "coordinates";
        int prel1 = std::ceil(0.5*(padlen - str1.length()));
        int postl1 = std::floor(0.5*(padlen - str1.length()));
        int prel2 = std::ceil(0.5*(padlen - str2.length()));
        int postl2 = std::floor(0.5*(padlen - str2.length()));

        printf("%*s%s%*s|%*s%s%*s\n", prel1, "", "id", postl1, "", prel2, "", "coordinates", postl2, "");
        for (it = coordinates.begin(); it != coordinates.end(); it++){
            printf("%+6d | %+1.6f %+1.6f %+1.6f\n",
                   it->first, it->second[0], it->second[1], it->second[2]);
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

    voro::container* construct_container(const std::map< unsigned int, std::vector< FloatType > >  &point_coords,
                                         const vecOfvec &bounds, std::vector< voro::wall_plane > &planes, double expand){
        /*!
        Returns the pointer to a new voro::container formed by the walls in planes and containing the points in point_coords.

        NOTE: delete this memory after using!

        :param std::map< unsigned int, std::vector< FloatType > > &point_coords: The map to the coordinates of the points
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

//        //Add the points to the pre-container
//        if (point_numbers.size() != point_coords.size()){
//            std::cout << "Error: The point indices and coordinates must have the same length\n";
//            assert(1==0);
//        }
        for (auto point = point_coords.begin(); point != point_coords.end(); point++){//unsigned int i=0; i<point_numbers.size(); i++){
            pcontainer.put(point->first, point->second[0], point->second[1], point->second[2]);
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

    void evaluate_container_information(const std::map< unsigned int, std::vector< FloatType > > &positions,
                                        const std::map< int, std::pair< std::vector< FloatType >, std::vector< FloatType > > > &bounding_faces,
                                        voro::container* container, integrateMap &points, 
                                        std::map< unsigned int, FloatType> &boundary_node_volumes){
        /*!
        Compute required container information (volumes, surface areas, etc.) and return them.

        The planes which describe the container's bounds are desired to have negative id numbers starting 
        at -1 and progressing onwards. It will be assumed that any negative number for a neighboring point 
        is actually a surface on the outside of the domain. The planes will be id'd by adding 1 and taking 
        the negative.

        :param const std::map< unsigned int, std::vector< FloatType > > &positions: A map from the dns node number to its local position
        :param const std::map< unsigned int, std::pair< std::vector< FloatType >, std::vector< FloatType > > > &bounding_faces: The bounding surfaces
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
//                std::cout << *viit << ", ";

                //Accept points cut by planes which will have a negative index
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
                else{

                    //Check if the current face is coincident with one of the bounding faces
                    std::vector< FloatType > tmp_normal({cell_normals[index+0],
                                                         cell_normals[index+1],
                                                         cell_normals[index+2]});

                    find_face_centroid(face_vertices, vertices, index_order, centroid);

                    //Iterate through the bounding faces
                    int temp_id = 0;
                    bool dircomp, planecmp;
                    for (auto face_it=bounding_faces.begin(); face_it!=bounding_faces.end(); face_it++){
                        dircomp = compare_vector_directions(tmp_normal, face_it->second.first);
                        planecmp = point_on_surface(centroid, face_it->second.first, face_it->second.second);

                        if ((dircomp) && (planecmp)){
                            temp_id = face_it->first;
                            break;
                        }
                    }
                    if (temp_id < 0){
                        planes.push_back(-(temp_id+1));
                        areas.push_back(*vdit);
                        normals.push_back(tmp_normal);
                        face_centroids.push_back(centroid);
                        is_boundary = true;

                    }
/*
                
                    //Check if a neighboring point is not included in the domain. This means the point is also on a boundary
                    auto nit = positions.find(*viit);
                    if (nit != positions.end()){
                        if (!container->point_inside(nit->second[0], nit->second[1], nit->second[2])){

                            //Get the normal and centroid information
                            std::vector< FloatType > tmp_normal({cell_normals[index+0],
                                                                 cell_normals[index+1],
                                                                 cell_normals[index+2]});

                            find_face_centroid(face_vertices, vertices, index_order, centroid);

                            //Determine the face's id number
                            int temp_id=0;
                            bool dircomp, planecmp;
//                            std::cout << "neighbor location: "; elib::print(nit->second);
//                            std::cout << "face centroid: "; elib::print(centroid);
//                            std::cout << "face normal:   "; elib::print(tmp_normal);
                            for (auto face_it=bounding_faces.begin(); face_it!=bounding_faces.end(); face_it++){
//                                std::cout << " plane id: " << face_it->first << "\n";
//                                std::cout << " plane normal: "; elib::print(face_it->second.first);
//                                std::cout << " plane point:  "; elib::print(face_it->second.second); std::cout << "\n";

                                dircomp = compare_vector_directions(tmp_normal, face_it->second.first);
                                planecmp = point_on_surface(centroid, face_it->second.first, face_it->second.second);

                                if ((dircomp) && (planecmp)){
                                    temp_id = face_it->first;
                                    break;
                                }
                            }
                            if (temp_id == 0){
                                std::cout << "Error: bounding volume has no associated normal\n";
                                //assert(1==0);
                            }
                            else{
                                planes.push_back(-(temp_id+1));
                                areas.push_back(*vdit);
                                normals.push_back(tmp_normal);
                                find_face_centroid(face_vertices, vertices, index_order, centroid);
                                face_centroids.push_back(centroid);
                                is_boundary = true;
                            }
                        }
                    }

*/
                }

                viit++;
                vdit++;
                index += 3;
                index_order += face_vertices[index_order]+1;
            }
//            std::cout << "\n";
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
            points.insert(std::pair< unsigned int, MicroPoint>(cl.pid(), MicroPoint(c.volume(), centroid, {x, y, z}, planes, areas, normals, face_centroids)));
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

        unsigned int n = domain.das.size();
//        unsigned int j=1, n = domain.das.size();
        double distance;
        std::vector< double > normal(3);
        vplanes.reserve(n);

//        std::cout << "\nmap_domain_to_voro:\n";

        for (unsigned int i=0; i<n; i++){
            normal = domain.normal(i);
            distance = overlap::dot(normal, domain.face_centroids[i]);
//            std::cout << "  true plane id: " << domain.planes[i] << "\n";
//            std::cout << "  voro plane id: " << -(domain.planes[i]+1) << "\n";
//            std::cout << "  normal: " << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
//            std::cout << "  face centroid: " << domain.face_centroids[i][0] << " " << domain.face_centroids[i][1] << " " << domain.face_centroids[i][2] << "\n";
//            std::cout << "  distance: " << distance << "\n";
            vplanes.push_back(voro::wall_plane(normal[0], normal[1], normal[2], distance, -(domain.planes[i]+1)));//j));
//            j++;
        }
//        assert(1==0);
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
        * Perform the volume integration of a scalar value multiplied by the position of the point.
        * Returns the integrated value at the gauss points.
        *
        * :param map::< unsigned int, double > values: The map object which defines the values to be integrated at each of the micro-nodes
        * :param std::vector< integrateMap > weights: The weights of each of the nodes in true space for each gauss point.
        * :param std::vector< std::vector< double > > result: The integrated result.
        */

        //Set up an iterator for the value map
        std::map<unsigned int, double >::const_iterator itv;
//        std::map<unsigned int, elib::vec>::const_iterator itp;
        integrateMap::const_iterator itiM;

        //Initialize the result vector
        result = std::vector< std::vector< double > >(weights.size());
        //Loop over the gauss points
        for (unsigned int gp=0; gp<weights.size(); gp++){
            auto tmpit = weights[gp].begin();
            if (weights[gp].size()==0){
                std::cout << "Warning: gauss point with no micro-scale points detected.\n";
                continue;
            }
            result[gp] = std::vector<double>(tmpit->second.coordinates.size(), 0);
//            result[gp] = std::vector<double>(weights[gp].begin()->second.coordinates.size(), 0);

            for (itiM=weights[gp].begin(); itiM!=weights[gp].end(); itiM++){

                //Find the value of the function at the node
                itv = values.find(itiM->first);
                if (itv == values.end()){
                    std::cout << "Error: node " << itiM->first << " not found in values\n";
                    assert(1==0);
                }
//                itp = positions.find(itiM->first);
//                if (itp == positions.end()){
//                    std::cout << "Error node " << itiM->first << " not found in positions\n";
//                    assert(1==0);
//                }

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

    void compute_surface_area_normal(const std::vector< integrateMap > &weights,
                                     scalar_surface_map &surface_area, vector_surface_map &surface_normal){
        /*!
         * Compute the surface areas and normals of each of the faces in the gauss domains
         * 
         * :param const std::vector< integrateMap > &weights: The weights of each of the nodes in true-space for each gauss point
         * :param scalar_surface_map &surface_area: The computed surface areas
         * :param vector_surface_map &surface_normal: The computed surface normals
         */

        integrateMap::const_iterator itiM;
        std::map< unsigned int, double >::iterator itr_sa;
        std::map< unsigned int, std::vector< double > >::iterator itr_sn;

        //Initialize the surface_area vector
        surface_area.resize(weights.size());
        surface_normal.resize(weights.size());

        //Loop over the gauss domains
        for (unsigned int gp=0; gp<weights.size(); gp++){
            
            //Loop over the micro-node weights
            for (itiM=weights[gp].begin(); itiM!=weights[gp].end(); itiM++){

                //Loop over the planes
                for (unsigned int j=0; j<itiM->second.planes.size(); j++){

                    itr_sa = surface_area[gp].find(itiM->second.planes[j]);
                    itr_sn = surface_normal[gp].find(itiM->second.planes[j]);

                    //Insert the plane if new
                    if (itr_sa == surface_area[gp].end()){
                        surface_area[gp].emplace(itiM->second.planes[j], itiM->second.area(j));
                        surface_normal[gp].emplace(itiM->second.planes[j], itiM->second.das[j]);
                    }
                    //Add to the plane if it exists already
                    else{
                        itr_sa->second += itiM->second.area(j);
                        for (unsigned int i=0; i<itr_sn->second.size(); i++){itr_sn->second[i] += itiM->second.das[j][i];}
                    }
                }
            }

//            std::cout << "gp: " << gp << "\n";
            for (auto ittmp = surface_normal[gp].begin(); ittmp!=surface_normal[gp].end(); ittmp++){
                auto area = surface_area[gp].find(ittmp->first);
                if (area == surface_area[gp].end()){
                    std::cout << "Error: surface area for normal plane not found\n";
                    assert(1==0);
                }
                else{
//                    std::cout << " face: " << ittmp->first << "area: " << area->second << " , da: "; elib::print(ittmp->second);
                    for (unsigned int i=0; i<ittmp->second.size(); i++){ittmp->second[i] /= area->second;}
                }
            }
        }


        return;
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
         * Perform the surface integration of a scalar value. Returns the integrated value at each of the surface of the gauss points
         * 
         * :param std::map< unsigned int, double > &values: The values at each of the micro-points
         * :param std::vector< integrateMap > weights: The weights of each of the nodes in true space for each gauss point
         * :param std::vector< std::map< unsigned int, std::vector< double > > > result: The result of the integration over each of the faces
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

    void perform_symmetric_tensor_surface_traction_integration(const std::map< unsigned int, std::vector< double > > &tensor, 
                                                               const std::vector< integrateMap > &weights,
                                                               std::vector< std::map< unsigned int, std::vector< double > > > &result){
        /*!
         * Integrate the fluxes of the symmetric tensor over the surfaces (assumes 3D) return will be a traction weighted by the area
         * 
         * :param const std::map< unsigned int, std::vector< double > > &tensor: The symmetric tensors in voigt notation i.e.
         *     tensor = t11, t22, t33, t23, t13, t12
         * :param std::vector< integrateMap > weights: The weights of each of the nodes in true space for each gauss point
         * :param std::vector< std::map< unsigned int, double > > result: The result of the integration over each of the faces
         */
        
        //Set up an iterator for the value map
        std::map< unsigned int, std::vector< double > >::const_iterator itv; //The tensor value iterator
        integrateMap::const_iterator itiM; //The integration map iterator
        std::map< unsigned int, std::vector< double > >::iterator itr; //The result iterator

        //Initialize the result vector
        result.resize(weights.size());

        //Loop over the gauss points
        for (unsigned int gp=0; gp<weights.size(); gp++){
        
           //Loop over the micro-node weights
           for (itiM=weights[gp].begin(); itiM!=weights[gp].end(); itiM++){

               //Find the value of the function at the current node
               itv = tensor.find(itiM->first);
               if (itv == tensor.end()){
                   std::cout << "Error: node " << itiM->first << " not found in values\n";
                   assert(1==0);
               }

               //Loop over the planes
               for (unsigned int j=0; j<itiM->second.planes.size(); j++){
                   itr = result[gp].find(itiM->second.planes[j]);

                   //Compute the traction - Assumes 3D
                   std::vector< double > traction = std::vector< FloatType >(3, 0);
                   traction[0] = itiM->second.das[j][0]*itv->second[0]
                                +itiM->second.das[j][1]*itv->second[5]
                                +itiM->second.das[j][2]*itv->second[4];

                   traction[1] = itiM->second.das[j][0]*itv->second[5]
                                +itiM->second.das[j][1]*itv->second[1]
                                +itiM->second.das[j][2]*itv->second[3];

                   traction[2] = itiM->second.das[j][0]*itv->second[4]
                                +itiM->second.das[j][1]*itv->second[3]
                                +itiM->second.das[j][2]*itv->second[2];

                   //Insert the plane if new
                   if (itr == result[gp].end()){
                       result[gp].emplace(itiM->second.planes[j], traction);
                   }
                   //Add to the plane if it exists already
                   else{
                       for (unsigned int i=0; i<traction.size(); i++){ itr->second[i] += traction[i];}
                   }
               }
           }
        }
    }

    void perform_symmetric_tensor_surface_couple_traction_integration(
             const std::map< unsigned int, std::vector< double > > &tensor,
             const std::vector< integrateMap > &weights,
             const vecOfvec &centers_of_mass,
             std::vector< std::map< unsigned int, std::vector< double > > > &result){
        /*!
         * Integrate the couple fluxes of the symmetric tensor over the surfaces (assumes 3D) return will be a couple traction
         * weighted by the area.
         * 
         * :param const std::map< unsigned int, std::vector< double > > &tensor: The symmetric tensors in voigt notation i.e.
         *     tensor = t11, t22, t33, t23, t13, t12
         * :param const std::vector< integrateMap > weights: The weights of each of the nodes in true space for each gauss point
         * :param const vecOfvec &centers_of_mass: The centers of mass of the gauss points.
         * :param std::vector< std::map< unsigned int, double > > result: The result of the integration over each of the faces
         */
        
        //Set up an iterator for the value map
        std::map< unsigned int, std::vector< double > >::const_iterator itv; //The tensor value iterator
        integrateMap::const_iterator itiM; //The integration map iterator
        std::map< unsigned int, std::vector< double > >::iterator itr; //The result iterator

        //Initialize the result vector
        result.resize(weights.size());

        //Initialize the center of mass and xi vectors
        std::vector< double > center_of_mass;
        std::vector< double > xi;

        //Loop over the gauss points
        for (unsigned int gp=0; gp<weights.size(); gp++){
            center_of_mass = centers_of_mass[gp];
            xi.resize(center_of_mass.size());

            //Loop over the micro-node weights
            for (itiM=weights[gp].begin(); itiM!=weights[gp].end(); itiM++){

                //Find the value of the tensor at the current node
                itv = tensor.find(itiM->first);
                if (itv == tensor.end()){
                    std::cout << "Error: node " << itiM->first << " not found in values\n";
                    assert(1==0);
                }


                for (unsigned int j=0; j<itiM->second.planes.size(); j++){
                    itr = result[gp].find(itiM->second.planes[j]);

                    //Compute the xi vector
                    for (unsigned int i=0; i<xi.size(); i++){xi[i] = itiM->second.face_centroids[j][i] - center_of_mass[i];}

                    //Compute the couple traction - Assumes 3D
                    std::vector< double > couple_traction(9, 0);
                    couple_traction[0] = itiM->second.das[j][0]*itv->second[0]*xi[0]
                                       + itiM->second.das[j][1]*itv->second[5]*xi[0]
                                       + itiM->second.das[j][2]*itv->second[4]*xi[0];

                    couple_traction[1] = itiM->second.das[j][0]*itv->second[5]*xi[1]
                                       + itiM->second.das[j][1]*itv->second[1]*xi[1]
                                       + itiM->second.das[j][2]*itv->second[3]*xi[1];

                    couple_traction[2] = itiM->second.das[j][0]*itv->second[4]*xi[2]
                                       + itiM->second.das[j][1]*itv->second[3]*xi[2]
                                       + itiM->second.das[j][2]*itv->second[2]*xi[2];

                    couple_traction[3] = itiM->second.das[j][0]*itv->second[5]*xi[2]
                                       + itiM->second.das[j][1]*itv->second[1]*xi[2]
                                       + itiM->second.das[j][2]*itv->second[3]*xi[2];

                    couple_traction[4] = itiM->second.das[j][0]*itv->second[0]*xi[2]
                                       + itiM->second.das[j][1]*itv->second[5]*xi[2]
                                       + itiM->second.das[j][2]*itv->second[4]*xi[2];

                    couple_traction[5] = itiM->second.das[j][0]*itv->second[0]*xi[1]
                                       + itiM->second.das[j][1]*itv->second[5]*xi[1]
                                       + itiM->second.das[j][2]*itv->second[4]*xi[1];

                    couple_traction[6] = itiM->second.das[j][0]*itv->second[4]*xi[1]
                                       + itiM->second.das[j][1]*itv->second[3]*xi[1]
                                       + itiM->second.das[j][2]*itv->second[2]*xi[1];

                    couple_traction[7] = itiM->second.das[j][0]*itv->second[4]*xi[0]
                                       + itiM->second.das[j][1]*itv->second[3]*xi[0]
                                       + itiM->second.das[j][2]*itv->second[2]*xi[0];

                    couple_traction[8] = itiM->second.das[j][0]*itv->second[5]*xi[0]
                                       + itiM->second.das[j][1]*itv->second[1]*xi[0]
                                       + itiM->second.das[j][2]*itv->second[3]*xi[0];

                   //Insert the plane if new
                   if (itr == result[gp].end()){
                       result[gp].emplace(itiM->second.planes[j], couple_traction);
                   }
                   //Add to the plane if it exists already
                   else{
                       for (unsigned int i=0; i<couple_traction.size(); i++){ itr->second[i] += couple_traction[i];}
                   }
                }
            }
        }
        return; 
    }

    void construct_couple_least_squares(const std::vector< std::map< unsigned int, std::vector< double > > > &surface_normals,
                                        const std::vector< std::map< unsigned int, std::vector< double > > > &surface_couples,
                                        Eigen::MatrixXd &A, Eigen::MatrixXd &b){
        /*!
         * Construct the normal matrix which can project the couple stresses at the Gauss points to the couple traction vector b 
         * on each of the surfaces. Assumes a 3D couple stress.
         * 
         * :param const std::vector< std::map< unsigned int, std::vector< double > > > &surface_normals: The vector of maps from the 
         *     gauss domain's face number to the normal of that face.
         * :param const std::vector< std::map< unsigned int, std::vector< double > > > &surface_couples: The vector of maps from the 
         *     gauss domain's face number to the surface couple traction on that face.
         * :param Eigen::MatrixXd &A: The normal matrix
         * :param Eigen::MatrixXd &b: The couple traction vector
         */

        //Determine the size of the A matrix
        unsigned int nrows, ncols;
        unsigned int nstress = 27;
        unsigned int dim = 9;

        //Find the faces that are unique
        std::vector< std::map< unsigned int, std::vector< double > > > unique_normals(surface_normals.size());
        for (unsigned int gp=0; gp<surface_normals.size(); gp++){
            id_unique_vectors(surface_normals[gp], unique_normals[gp], false);

//            std::cout << "gp " << gp << "\n";
//            for (auto un=unique_normals[gp].begin(); un!=unique_normals[gp].end(); un++){
//                std::cout << un->first << ": "; print_vector(un->second);
//            }
        }

        nrows = 0;
        ncols = nstress*surface_normals.size(); //27 components of the couple stress for each gauss point
        for (unsigned int gp=0; gp<surface_normals.size(); gp++){
//            nrows += dim*surface_normals[gp].size();
            nrows += dim*unique_normals[gp].size();

        }

        if (surface_couples.size() != surface_normals.size()){
            std::cerr << "Error: surface_couples should have the same size as surface_normals\n";
            std::cerr << "       surface_normals.size(): " << surface_normals.size() << "\n";
            std::cerr << "       surface_couples.size(): " << surface_couples.size() << "\n";
            assert(1==0);
        }

        //Resize A and b
        A = Eigen::MatrixXd::Zero(nrows, ncols);
        b = Eigen::MatrixXd::Zero(nrows, 1);

        //Iterate over the gauss points
        unsigned int row0, col0;
        row0 = col0 = 0;

//        for (unsigned int gp=0; gp<surface_normals.size(); gp++){
        for (unsigned int gp=0; gp<unique_normals.size(); gp++){

//            for (auto face = surface_normals[gp].begin(); face != surface_normals[gp].end(); face++){
            for (auto face = unique_normals[gp].begin(); face != unique_normals[gp].end(); face++){
                //Set the values in the A matrix

                A(row0 + 0, col0 +  0) = face->second[0];
                A(row0 + 0, col0 +  9) = face->second[1];
                A(row0 + 0, col0 + 18) = face->second[2];
                A(row0 + 1, col0 +  1) = face->second[0];
                A(row0 + 1, col0 + 10) = face->second[1];
                A(row0 + 1, col0 + 19) = face->second[2];
                A(row0 + 2, col0 +  2) = face->second[0];
                A(row0 + 2, col0 + 11) = face->second[1];
                A(row0 + 2, col0 + 20) = face->second[2];
                A(row0 + 3, col0 +  3) = face->second[0];
                A(row0 + 3, col0 + 12) = face->second[1];
                A(row0 + 3, col0 + 21) = face->second[2];
                A(row0 + 4, col0 +  4) = face->second[0];
                A(row0 + 4, col0 + 13) = face->second[1];
                A(row0 + 4, col0 + 22) = face->second[2];
                A(row0 + 5, col0 +  5) = face->second[0];
                A(row0 + 5, col0 + 14) = face->second[1];
                A(row0 + 5, col0 + 23) = face->second[2];
                A(row0 + 6, col0 +  6) = face->second[0];
                A(row0 + 6, col0 + 15) = face->second[1];
                A(row0 + 6, col0 + 24) = face->second[2];
                A(row0 + 7, col0 +  7) = face->second[0];
                A(row0 + 7, col0 + 16) = face->second[1];
                A(row0 + 7, col0 + 25) = face->second[2];
                A(row0 + 8, col0 +  8) = face->second[0];
                A(row0 + 8, col0 + 17) = face->second[1];
                A(row0 + 8, col0 + 26) = face->second[2];

                auto couple = surface_couples[gp].find(face->first);
                if (couple == surface_couples[gp].end()){
                    std::cerr << "Error: surface couple for face " << face->first << " not found\n";
                    assert(1==0);
                }

                for (unsigned int i=0; i<dim; i++){
                    b(row0 + i, 0) = couple->second[i];
                }

                row0 += dim; //Increment row0
            }
            col0 += nstress; //Increment col0
        }
        return;
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
        tripletList.reserve(tripletList.size() + psis.size()*num_macro_dof*dns_weights.size());

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
//            std::cout << "dns_weights.size(): " << dns_weights.size() << "\n";
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
//                    std::cout << "dns node, weight, elcount: " << iMit->first << ", " << iMit->second.weight << ", " << mnelit->second << "\n";
                    
                }

                //Compute the xi vector
                for (unsigned int j=0; j<3; j++){
//                    xi[j] = iMit->second.coordinates[j] - cg[j];
                      xi[j] = iMit->second.particle_coordinates[j] - cg[j];
                }

/*                bool diff_centroid_particle = !fuzzy_equals(iMit->second.coordinates, iMit->second.particle_coordinates);
                if (diff_centroid_particle){//(it->second == 450){
                    std::cout << "id: " << it->second << "\n";
                    std::cout << "row0: " << row0 << "\n";
                    std::cout << "col0: " << col0 << "\n";
                    std::cout << "share_ghost_free_boundary_nodes: " << share_ghost_free_boundary_nodes << "\n";
                    std::cout << "macro_elem_is_ghost: " << macro_elem_is_ghost << "\n";
                    std::cout << "num_micro_free: " << num_micro_free << "\n";
                    std::cout << "iMit->second.coordinates: "; elib::print(iMit->second.coordinates);
                    std::cout << "iMit->second.particle_coordinates: "; elib::print(iMit->second.particle_coordinates);
                    
                    std::cout << "cg: "; elib::print(cg);
                    std::cout << "xi: "; elib::print(xi);
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

    MicromorphicFilter::MicromorphicFilter(const unsigned int _id, const std::string &element_type,
                                           const std::vector< unsigned int > &global_node_ids, const elib::vecOfvec &nodes,
                                           const elib::quadrature_rule &qrule, const unsigned int num_macro_dof,
                                           bool _shared_dof_material, bool _use_dns_bounds){
        /*!
	 * Initialize the micromorphic filter
	 *
         * :param unsigned int id: The id number of the filter
	 * :param const std::string &element_type: The name of the underlying element type
         * :param const std::vector< unsigned int > global_node_ids: The global node id numbers
	 * :param const elib::vecOfvec &nodes: The nodes which make up the element type (i.e. their coordinates)
	 * :param const elib::quadrature_rule &qrule: The quadrature rule of the element.
         * :param const unsigned int num_macro_dof: The number of macro-scale degrees of freedom.
	 * :param bool _shared_dof_material: A flag which indicates if the material points should be used to 
	 *                                   determine the macro-scale nodal coordinate information.
         * :param bool _use_dns_bounds: A flag which indicates whether the DNS bounds should be used for the overlap.
	 */

        filter_id = _id;
        element = elib::build_element_from_string(element_type, global_node_ids, nodes, qrule);
        shared_dof_material = _shared_dof_material;
        use_dns_bounds = _use_dns_bounds;

        elib::vecOfvec gauss_points(element->qrule.size());
        for (unsigned int gp=0; gp<element->qrule.size(); gp++){
            gauss_points[gp] = element->qrule[gp].first;
        }

        material_overlap.initialize(element->local_node_coordinates, gauss_points);

        if (!shared_dof_material){
            dof_overlap.initialize(element->local_node_coordinates, gauss_points);
        }

        dof_values.resize(nodes.size());
        for (unsigned int i=0; i<nodes.size(); i++){
            dof_values[i] = std::vector< FloatType >(num_macro_dof, 0);
        }

        filter_dim = nodes[0].size();

    }

    bool MicromorphicFilter::add_micro_dof_point(const unsigned int &id, const elib::vec &coordinates, FloatType tol){
        /*!
        Check if a dof point is inside the filter and add it if it is.

        :param unsigned int &id: The micro point id number
        :param elib::vec &coordinates: The global coordinates of the micro dof point.
        */

        if (element->bounding_box_contains_point(coordinates)){

            elib::vec xi;
            int clc_result = element->compute_local_coordinates(coordinates, xi);
            if ((element->local_point_inside(xi, tol)) and (clc_result == 0)){
//                dof_id_numbers.push_back(id);
//                micro_dof_local_coordinates.push_back(xi);
                micro_dof_local_coordinates.emplace(id, xi);
                return true;
            }
        }
        return false;
    }

    bool MicromorphicFilter::add_micro_material_point(const unsigned int &id, const elib::vec &coordinates, FloatType tol){
        /*!
        * Check if a material point is inside the filter and add it if it is.
        * 
        * :param unsigned int &id: The micro point id number
        * :param elib::vec &coordinates: The global coordinates of the micro material point.
        */

        if (element->bounding_box_contains_point(coordinates)){

            elib::vec xi;
            element->compute_local_coordinates(coordinates, xi);
            if (element->local_point_inside(xi, tol)){
                micro_material_local_coordinates.emplace(id, xi);
//                material_id_numbers.push_back(id);
//                micro_material_local_coordinates.push_back(xi);
                return true;
            }
        }
        return false;
    }

    int MicromorphicFilter::construct_integrators(bool update_shapefunction_matrix){
        /*!
        * Construct the filter's integrators.
        * 
        * :param bool update_shapefunction_matrix: Flag indicating if the shape-function matrix
        *                                          should be updated.
        */
        
//        std::cout << "material_id_numbers.size(): " << material_id_numbers.size() << "\n";

        construct_material_point_integrator();
        if ((!shared_dof_material) && (update_shapefunction_matrix)){
            construct_dof_point_integrator();
        }
        return 0;
    }

    int MicromorphicFilter::construct_material_point_integrator(){
        /*!
        * Construct the integrator for the material points.
        */

//        material_overlap.compute_weights(material_id_numbers, micro_material_local_coordinates, material_weights, use_dns_bounds);
        material_overlap.compute_weights(micro_material_local_coordinates, material_weights, use_dns_bounds);

        //Compute the shapefunctions at the gauss domain face centroids
        compute_face_centroid_shapefunctions();

        //Transform the volumes and normals
        elib::vecOfvec jacobian;
        elib::vecOfvec invjacobian;
        double je;
        for (unsigned int gp=0; gp<material_weights.size(); gp++){
/*            std::cout << "points in gauss point " << gp << ": " << material_weights[gp].size() << "\n";
            if (material_weights[gp].size()<10){

                const std::vector< MicroPoint >* gd = material_overlap.get_gauss_domains();
                (*gd)[gp].print();

                for (unsigned int i=0; i<micro_material_local_coordinates.size(); i++){
                    elib::print(micro_material_local_coordinates[i]);
                }
                assert(1==0);
            }
*/
            for (auto it=material_weights[gp].begin(); it!=material_weights[gp].end(); it++){
//                std::cout << "material point: " << it->first << "\n";
//                it->second.print();

                element->get_jacobian(it->second.coordinates, element->local_node_coordinates, jacobian);
                elib::determinant_3x3(jacobian, je);
                elib::invert(jacobian, invjacobian);

                //Transform the volume
                it->second.volume *= je;

                //Transform the normals
                for (unsigned int j=0; j<it->second.das.size(); j++){
                    overlap::apply_nansons_relation(it->second.normal(j), je*it->second.area(j), invjacobian, it->second.das[j]);
                }

                //Transform the coordinates to true-space
                elib::vec lc = it->second.coordinates;
                element->interpolate(element->nodes, lc, it->second.coordinates);
                lc = it->second.particle_coordinates;
                element->interpolate(element->nodes, lc, it->second.particle_coordinates);

                //Transform the face centroids to true-space
                for (unsigned int n=0; n<it->second.face_centroids.size(); n++){
                    lc = it->second.face_centroids[n];
                    element->interpolate(element->nodes, lc, it->second.face_centroids[n]);
                }

//                std::cout << "true space\n";
//                it->second.print();
            }
//            assert(2==3);
        }
        return 0;
    }

    int MicromorphicFilter::construct_dof_point_integrator(){
        /*!
        * Construct the integrator for the degree of freedom points.
        */
//        dof_overlap.compute_weights(dof_id_numbers, micro_dof_local_coordinates, dof_weights, use_dns_bounds);
        dof_overlap.compute_weights(micro_dof_local_coordinates, dof_weights, use_dns_bounds);

        //Transform the volumes and normals
        elib::vecOfvec jacobian;
        elib::vecOfvec invjacobian;
        double je;
        for (unsigned int gp=0; gp<dof_weights.size(); gp++){
            for (auto it=dof_weights[gp].begin(); it!=dof_weights[gp].end(); it++){
                element->get_jacobian(it->second.coordinates, element->local_node_coordinates, jacobian);
                elib::determinant_3x3(jacobian, je);
                elib::invert(jacobian, invjacobian);

                //Transform the volume
                it->second.volume *= je;

                //Transform the normals
                for (unsigned int j=0; j<it->second.das.size(); j++){
                    overlap::apply_nansons_relation(it->second.normal(j), je*it->second.area(j), invjacobian, it->second.das[j]);
                }

                //Transform the coordinates to true-space
                elib::vec lc = it->second.coordinates;
                element->interpolate(element->nodes, lc, it->second.coordinates);
                lc = it->second.particle_coordinates;
                element->interpolate(element->nodes, lc, it->second.particle_coordinates);

                //Transform the face centroids to true-space
                for (unsigned int n=0; n<it->second.face_centroids.size(); n++){
                    lc = it->second.face_centroids[n];
                    element->interpolate(element->nodes, lc, it->second.face_centroids[n]);
                }
            }
        }
        return 0;
    }

    int MicromorphicFilter::compute_mass_properties(const std::map< unsigned int, double> &micro_density){
        /*!
         * Compute the mass properties of the homogenized DNS.
         *
         * :param std::map< unsigned int , double> &density: The densities of the micro points
         */

        compute_volume();
        compute_surface_area_normal();
        compute_density(micro_density);
        compute_centers_of_mass(micro_density);
        compute_com_shapefunction_gradients();
        return 0;
    }

    int MicromorphicFilter::compute_stress_properties(const std::map< unsigned int, std::vector< double > > &micro_stress){
        /*!
         * Compute the stress properties of the homogenized DNS.
         * 
         * :param std::map< unsigned int, std::vector< double > > &micro_stress: The stresses of the micro-points.
         */

         //Compute the symmetric microstress
         compute_symmetric_microstress(micro_stress);

         //Compute the Cauchy stress
         compute_traction(micro_stress);
         compute_vertices_cauchy_stress();
//         construct_cauchy_least_squares();
         construct_linear_momentum_surface_external_force();
         construct_linear_momentum_least_squares_matrix();
         //TODO: Add construction of body force term
         //TODO: Add construction of kinetic force term
         construct_linear_momentum_b_vector();
         construct_weight_constraints();
         compute_cauchy_stress();

         //Compute the couple stress
         compute_couple_traction(micro_stress);
         compute_vertices_couple_stress();
         construct_couple_least_squares();
         construct_first_moment_surface_external_couple();
         construct_first_moment_symm_cauchy_couple();
         construct_first_moment_constraint_matrix();
         //TODO: Add construction of body couple term
         //TODO: Add construction of kinetic couple term
         construct_first_moment_d_vector();
         compute_couple_stress();
//         std::cout << "couple_stress:\n";
//         for (unsigned int gp=0; gp<couple_stress.size(); gp++){
//             std::cout << "gp " << gp << ": "; print_vector(couple_stress[gp]);
//         }
//         assert(1==0);

         return 0;
    }

    int MicromorphicFilter::construct_linear_momentum_b_vector(){
        /*!
         * Construct the d vector for the constrained least squares calculation 
         * of the cauchy stress
         */

        overlap::construct_linear_momentum_b_vector(linear_momentum_A.rows(),
                                                    surface_external_force, body_external_force, kinetic_force,
                                                    linear_momentum_b);
        return 0;
    }

    int MicromorphicFilter::construct_first_moment_d_vector(){
        /*!
         * Construct the d vector for the constrained first moment of momentum calculation
         * of the couple stress
         */

        overlap::construct_first_moment_d_vector(first_moment_C.rows(),
                                                 surface_external_couple, symm_cauchy_couple,
                                                 body_external_couple, kinetic_couple,
                                                 first_moment_d);
        return 0;
    }

    int MicromorphicFilter::construct_weight_constraints(){
        /*!
         * Construct the constraint matrix and vector for the weights
         * 
         * TODO: Add a constraint that the weights must be greater than zero
         */

        const std::vector< vecOfvec >* domain_vertices = material_overlap.get_domain_vertices();

        unsigned int ncol = 0;
        for (unsigned int gp=0; gp<(*domain_vertices).size(); gp++){
            ncol += (*domain_vertices)[gp].size();
        }

        weight_constraints = Eigen::MatrixXd::Zero((*domain_vertices).size(), ncol);
        linear_momentum_d = Eigen::MatrixXd::Zero((*domain_vertices).size(), 1);

        unsigned int col0 = 0;

        for (unsigned int gp=0; gp<(*domain_vertices).size(); gp++){
            for (unsigned int i=0; i<(*domain_vertices)[gp].size(); i++){
                weight_constraints(gp, col0 + i) = 1.;
            }
            linear_momentum_d(gp) = 1;
            col0 += (*domain_vertices)[gp].size();
        }
        return 0;
    }

    int MicromorphicFilter::compute_cauchy_stress(){
        /*!
         * Compute the Cauchy stress using a constrained least squares technique
         */

        Eigen::MatrixXd x;
//        solve_constrained_least_squares(linear_momentum_A, linear_momentum_b, linear_momentum_C, linear_momentum_d, x);
        Eigen::MatrixXd w;
        solve_constrained_least_squares(linear_momentum_A, linear_momentum_b, weight_constraints, linear_momentum_d, w);

        //Convert the weights to the Cauchy stresses
        cauchy_stress.resize(symmetric_microstress.size());
        unsigned int nstress;
        double weight;
        unsigned int col0 = 0;
        const std::vector< vecOfvec >* domain_vertices = material_overlap.get_domain_vertices();

        for (unsigned int gp=0; gp<cauchy_stress.size(); gp++){
            nstress = vertex_cauchy[gp][0].size();
            cauchy_stress[gp] = std::vector< FloatType >(nstress, 0);
            for (unsigned int v=0; v<(*domain_vertices)[gp].size(); v++){
                weight = w(col0 + v);
                for (unsigned int i=0; i<nstress; i++){
                    cauchy_stress[gp][i] += weight*vertex_cauchy[gp][v][i];
                }
            }
            col0 += (*domain_vertices)[gp].size();
        }

//        //Extract the Cauchy stresses
//        cauchy_stress.resize(symmetric_microstress.size());
//        unsigned int nstress = x.rows()/cauchy_stress.size();
//        for (unsigned int gp=0; gp<cauchy_stress.size(); gp++){
//            cauchy_stress[gp].resize(nstress);
//            for (unsigned int i=0; i<nstress; i++){
//                cauchy_stress[gp][i] = x(nstress*gp + i);
//            }
//        }
        return 0;
    }

    int MicromorphicFilter::compute_couple_stress(){
        /*!
         * Compute the couple stress using a constrained least squares technique
         */

        Eigen::MatrixXd x;

        std::ofstream file("couple_lsq_test.txt");
        if (file.is_open()){
            file << "A:\n" << first_moment_A << "\n";
            file << "b:\n" << first_moment_b << "\n";
            file << "C:\n" << first_moment_C << "\n";
            file << "d:\n" << first_moment_d << "\n";
            file.flush();
        }
        assert(1==0);

        solve_constrained_least_squares(first_moment_A, first_moment_b, first_moment_C, first_moment_d, x);

        //Extract the couple stresses assuming 3d
        couple_stress.resize(cauchy_stress.size());
        unsigned int nstress = x.rows()/couple_stress.size();
        for (unsigned int gp=0; gp<couple_stress.size(); gp++){
            couple_stress[gp].resize(nstress);
            for (unsigned int i=0; i<nstress; i++){
                couple_stress[gp][i] = x(nstress*gp + i);
            }
        }
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

    int MicromorphicFilter::compute_surface_area_normal(){
        /*!
         * Compute the surface areas and average normals of each of the filter's gauss domains
         */

        overlap::compute_surface_area_normal(material_weights, surface_area, surface_normal);
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
         * Compute the centers of mass of the material.
         * 
         * :param std::map< unsigned int , elib::vec> &density:
         */

//        std::cout << "entering position weighted volume integration\n";
        perform_position_weighted_volume_integration(micro_density, material_weights, center_of_mass);

//        std::cout << "material_weights[0]:\n";
//        for (auto p = material_weights[0].begin(); p!=material_weights[0].end(); p++){
//            std::cout << p->first << " ";
//            p->second.print();
//
//        }
        
//        std::cout << "computing local coordinates\n";
        local_center_of_mass.resize(center_of_mass.size());
        com_shapefunction_values.resize(center_of_mass.size());
        for (unsigned int gp=0; gp<density.size(); gp++){
            for (unsigned int i=0; i<center_of_mass[gp].size(); i++){
                center_of_mass[gp][i] /= volume[gp]*density[gp];
            }
            element->compute_local_coordinates(center_of_mass[gp], local_center_of_mass[gp]);
            element->get_shape_functions(local_center_of_mass[gp], com_shapefunction_values[gp]);
        }


        return 0;
    }

    int MicromorphicFilter::compute_symmetric_microstress(const std::map< unsigned int, std::vector< double > > &micro_cauchy){
        /*!
         * Compute the symmetric microstress (i.e. the volume average of the micro-scale's Cauchy stress)
         * 
         * :param const std::map< unsigned int, std::vector< double > > &micro_cauchy: The micro-scale cauchy stresses. It is 
         *     assumed they are symmetric and stored in the form (s11, s22, s33, s23, s13, s12)
         */
        
        perform_volume_integration(micro_cauchy, material_weights, symmetric_microstress);

        //Divide by the volume
        for (unsigned int gp=0; gp<symmetric_microstress.size(); gp++){
            for (unsigned int i=0; i<symmetric_microstress[gp].size(); i++){symmetric_microstress[gp][i] /= volume[gp];}
        }

        return 0;
    }

    int MicromorphicFilter::compute_traction(const std::map< unsigned int, std::vector< double > > &micro_cauchy){
        /*!
         * Compute the traction of the micro-cauchy stress through the each of the surfaces of the gauss points.
         * 
         * :param const std::map< unsigned int, std::vector< double > > &micro_cauchy: The micro-scale cauchy stresses. It is 
         *     assumed they are symmetric and stored in the form (s11, s22, s33, s23, s13, s12)
         */

        perform_symmetric_tensor_surface_traction_integration(micro_cauchy, material_weights, traction);

        //Divide by the surface area
        for (unsigned int gp=0; gp<traction.size(); gp++){
            for (auto face=traction[gp].begin(); face!=traction[gp].end(); face++){
                auto area = surface_area[gp].find(face->first);
                if (area == surface_area[gp].end()){
                    std::cerr << "Error: face " << face->first << "not found in surface areas\n";
                    assert(1==0);
                }
                for (unsigned int i=0; i<face->second.size(); i++){
                    face->second[i] /= area->second;
                }
            }
        }
        return 0;
    }

    int MicromorphicFilter::compute_vertices_cauchy_stress(){
        /*!
         * Compute the Cauchy stress at the Gauss domain vertices
         */

        const std::vector< std::vector< std::vector< unsigned int > > > *vertex_planes = material_overlap.get_vertex_planes();
        const std::vector< MicroPoint > *gauss_domains = material_overlap.get_gauss_domains();
        vecOfvec vertex_normals, vertex_tractions;

        vertex_cauchy.clear();
        vertex_cauchy.resize((*vertex_planes).size());

        //Iterate through the gauss points
        for (unsigned int gp=0; gp<(*vertex_planes).size(); gp++){

           overlap::compute_vertices_cauchy_stress((*vertex_planes)[gp], (*gauss_domains)[gp], traction[gp], vertex_cauchy[gp]);

        }
        return 0;
    }

    int MicromorphicFilter::compute_vertices_couple_stress(){
        /*!
         * Compute the higher order stress at the Gauss domain vertices
         */
        const std::vector< std::vector< std::vector< unsigned int > > > *vertex_planes = material_overlap.get_vertex_planes();
        const std::vector< MicroPoint > *gauss_domains = material_overlap.get_gauss_domains();
        vecOfvec vertex_normals;

        vertex_hostress.clear();
        vertex_hostress.resize((*vertex_planes).size());

        //Iterate through the gauss points
        for (unsigned int gp=0; gp<(*vertex_planes).size(); gp++){
            overlap::compute_vertices_couple_stress((*vertex_planes)[gp], (*gauss_domains)[gp], couple_traction[gp], vertex_hostress[gp]);
        }
        return 0;
    }

    int MicromorphicFilter::compute_couple_traction(const std::map< unsigned int, std::vector< double > > &micro_cauchy){
        /*!
         * Compute the couple traction of the micro-cauchy stress through the each of the surfaces of the gauss points.
         * 
         * :param const std::map< unsigned int, std::vector< double > > &micro_cauchy: The micro-scale cauchy stresses. It is 
         *     assumed they are symmetric and stored in the form (s11, s22, s33, s23, s13, s12)
         */

        perform_symmetric_tensor_surface_couple_traction_integration(micro_cauchy, material_weights, center_of_mass, couple_traction);

        //Divide by the surface area
        for (unsigned int gp=0; gp<couple_traction.size(); gp++){
//            std::cout << "gp: " << gp << "\n";
            for (auto face=couple_traction[gp].begin(); face!=couple_traction[gp].end(); face++){
                auto area = surface_area[gp].find(face->first);
                if (area == surface_area[gp].end()){
                    std::cerr << "Error: face " << face->first << "not found in surface areas\n";
                    assert(1==0);
                }
                for (unsigned int i=0; i<face->second.size(); i++){
                    face->second[i] /= area->second;
                }
//                std::cout << " " << face->first << ": "; print_vector(face->second);
            }
        }

        return 0;
    }

//    int MicromorphicFilter::construct_cauchy_least_squares(){
//        /*!
//         * Construct the normal matrix that when dotted with the vector of the Cauchy stress at the gauss domain 
//         * CGs will produce the surface traction vector b;
//         */
//
//        overlap::construct_cauchy_least_squares(surface_normal, traction, linear_momentum_A, linear_momentum_b);
//        return 0;
//    }

    int MicromorphicFilter::construct_couple_least_squares(){
        /*!
         * Construct the normal matrix that when dotted with the vector of the couple stress at the gauss domain 
         * CGs will produce the surface couple vector b;
         */

        overlap::construct_couple_least_squares(surface_normal, couple_traction, first_moment_A, first_moment_b);
        return 0;
    }
    int MicromorphicFilter::construct_linear_momentum_surface_external_force(){
        /*!
         * Construct the external force applied on the surfaces of the gauss domains.
         */
         const std::vector< std::vector< unsigned int > > *external_face_ids = material_overlap.get_external_face_ids();
         overlap::construct_linear_momentum_surface_external_force(face_shapefunctions, traction,
                                                                   surface_area, *external_face_ids,
                                                                   surface_external_force);
         return 0;
    }

    int MicromorphicFilter::construct_first_moment_surface_external_couple(){
        /*!
         * Construct the external couple applied on the surfaces of the gauss domains
         */
         const std::vector< std::vector< unsigned int > > *external_face_ids = material_overlap.get_external_face_ids();
         overlap::construct_first_moment_surface_external_couple(face_shapefunctions, couple_traction,
                                                                 surface_area, *external_face_ids,
                                                                 surface_external_couple);
         return 0;
    }

    int MicromorphicFilter::construct_first_moment_symm_cauchy_couple(){
        /*!
         * Construct the couple resulting from the difference between the Cauchy stress and symmetric microstress
         */
         overlap::construct_first_moment_symm_cauchy_couple(com_shapefunction_values,
                                                            symmetric_microstress, cauchy_stress,
                                                            volume,
                                                            symm_cauchy_couple);
         return 0;
    }

    int MicromorphicFilter::construct_linear_momentum_least_squares_matrix(){
        /*!
         * Construct the least_squares matrix for the balance of linear momentum.
         */

        overlap::construct_linear_momentum_least_squares_matrix(com_shapefunction_gradients, volume, vertex_cauchy, linear_momentum_A);
        return 0;
    }

    int MicromorphicFilter::construct_first_moment_constraint_matrix(){
        /*!
         * Construct the constraint matrix for the balance of the first moment of momentum.
         */

        overlap::construct_first_moment_constraint_matrix(com_shapefunction_gradients, volume, first_moment_C);
        return 0;
    }

    int MicromorphicFilter::compute_face_centroid_shapefunctions(){
        /*!
         * Compute the shapefunction values at the centroids of the gauss-domain faces
         */

        //Get a reference to the gauss domains
        const std::vector< MicroPoint >* gauss_domains = material_overlap.get_gauss_domains();

        //Initialize the face_shapefunctions vector
        face_shapefunctions.resize(gauss_domains->size());

        //Loop over the gauss domains
        unsigned int index = 0;
        for (auto domain = gauss_domains->begin(); domain!=gauss_domains->end(); domain++){
            for (unsigned int fc=0; fc<(*domain).face_centroids.size(); fc++){
                std::vector< FloatType > fc_shapefunctions;
                element->get_shape_functions((*domain).face_centroids[fc], fc_shapefunctions);
                face_shapefunctions[index].emplace((*domain).planes[fc], fc_shapefunctions);
            }
            index++;
        }
        return 0;
    }

    int MicromorphicFilter::compute_com_shapefunction_gradients(){
        /*!
         * Compute the global gradients of the shapefunctions at the centers of mass
         */
        
        com_shapefunction_gradients.resize(local_center_of_mass.size());

        for (unsigned int lcom=0; lcom<local_center_of_mass.size(); lcom++){
            element->get_global_shapefunction_gradients(local_center_of_mass[lcom], com_shapefunction_gradients[lcom]);
        }
        return 0;
    }

    int MicromorphicFilter::add_shapefunction_matrix_contribution(const std::map< unsigned int, unsigned int > &macro_node_to_col,
                                                                  const std::map< unsigned int, unsigned int > &micro_node_to_row,
                                                                  const std::vector< unsigned int > &macro_node_ids,
                                                                  const std::map< unsigned int, unsigned int > &micro_node_elcount,
                                                                  const unsigned int num_macro_dof, const unsigned int num_micro_dof,
                                                                  const unsigned int num_micro_free,
                                                                  std::vector< T > &tripletList){
        /*!
         * Add the contribution of this filter to the shape-function matrix
         * :param const std::map< unsigned int, unsigned int > &macro_node_to_col: The map from the macro nodes to the column of the 
         *                                                                         matrix (divided by num_macro_dof) which corresponds
         * :param const std::map< unsigned int, unsigned int > &micro_node_to_row: The map from the macro nodes to the row of the matrix
         *                                                                         (divided by num_micro_dof) which corresponds.
         * :param const std::vector< unsigned int > &macro_node_ids: The id numbers of the macro nodes.
         * :param const std::map< unsigned int, unsigned int > &micro_node_elcount: The number of elements the micro node is contained within
         * :param const unsigned int num_macro_dof: The number of degrees of freedom in the macro-scale
         * :param const unsigned int num_micro_dof: The number of degrees of freedom in the micro-scale
         * :param const unsigned int num_micro_free: The number of free micro degrees of freedom
         * :param std::vector< T > &tripletList: The resulting terms of the sparse matrix.
         */

        elib::vecOfvec phis;
        elib::vecOfvec cg_phis;

        //Compute the shape-function values at the centers of mass
        get_cg_phis(cg_phis);

//TEMP
//        std::cout << "cg_phis:\n";
//        elib::print(cg_phis);
//
//        std::cout << "local_com:\n";
//        elib::print(local_center_of_mass);
//
//        std::cout << "gpts:\n";
//        elib::print(element->qrule);
//
//        elib::print(macro_node_ids);
//        std::cout << "macro_node_elcount:\n";
//        for (auto it=micro_node_elcount.begin(); it!=micro_node_elcount.end(); it++){
//            std::cout << it->first << ": " << it->second << "\n";
//        }
//        std::cout << "num_macro_dof: " << num_macro_dof << "\n";
//        std::cout << "num_micro_dof: " << num_micro_dof << "\n";
//        std::cout << "num_micro_free: " << num_micro_free << "\n";
//ENDTEMP

        //Iterate over the filter's gauss points
        if (shared_dof_material){
            for (unsigned int gp=0; gp<material_weights.size(); gp++){
                //Re-cast the cg_phis at the gauss point to the shape expected by construct triplet list
                phis.resize(cg_phis[gp].size());
                for (unsigned int i=0; i<phis.size(); i++){phis[i].resize(1); phis[i][0] = cg_phis[gp][i];}

                //Add the terms to the triplet list
                overlap::construct_triplet_list(&macro_node_to_col, &micro_node_to_row, macro_node_ids,
                                                center_of_mass[gp], phis, material_weights[gp],
                                                &micro_node_elcount, false, true, num_micro_free,
                                                        tripletList, num_macro_dof, num_micro_dof);
            }
        }
        else{

            for (unsigned int gp=0; gp<dof_weights.size(); gp++){
                //Re-cast the cg_phis at the gauss point to the shape expected by construct triplet list
                phis.resize(cg_phis[gp].size());
                for (unsigned int i=0; i<phis.size(); i++){phis[i].resize(1); phis[i][0] = cg_phis[gp][i];}

                //Add the terms to the triplet list
                overlap::construct_triplet_list(&macro_node_to_col, &micro_node_to_row, macro_node_ids,
                                                center_of_mass[gp], phis, dof_weights[gp],
                                                &micro_node_elcount, false, true, num_micro_free,
                                                tripletList, num_macro_dof, num_micro_dof);
            }
        }
        return 0;
    }

    int MicromorphicFilter::print(const bool show_microscale_info){
        /*!
	 * Print out the information contained in the filter
         *
         * :param const bool show_microscale_info: Print out the micro-scale points 
         *     contained in the filter. This can result in huge amounds of output.
         *     defaults to false.
	 */

        if (show_microscale_info){
            std::cout << "DOF information: (id, coordinates)\n";
            print_coordinateMap(micro_dof_local_coordinates);
//            elib::print(dof_id_numbers);
            std::cout << "Material Point information (id, coordinates):\n";
            print_coordinateMap(micro_material_local_coordinates);
//            std::cout << "DOF local coordinates:\n";
//            elib::print(micro_dof_local_coordinates);
//            std::cout << "Material Point local coordinates:\n";
//            elib::print(micro_material_local_coordinates);
        }
        std::cout << "Element planes:\n";
        const planeMap *ep = material_overlap.get_element_planes();
        const planeMap *dnsp = material_overlap.get_dns_planes();
        print_planeMap(*ep);
        std::cout << "DNS planes:\n";
        print_planeMap(*dnsp);

        std::cout << "**************************\n";
        std::cout << "*** ELEMENT PROPERTIES ***\n";
        std::cout << "**************************\n";
        elib::print(*element);

        //Print the degree of freedom properties
        std::cout << "*********************************\n";
        std::cout << "*** DEGREES OF FREEDOM VALUES ***\n";
        std::cout << "*********************************\n";
        elib::print(dof_values);
        //Print the mass properties
        print_mass_properties();
	return 0;
    }

    int MicromorphicFilter::print_mass_properties(){
        /*!
         * Print out the mass properties of the filter
         */

        std::cout << "***********************\n";
        std::cout << "*** MASS PROPERTIES ***\n";
        std::cout << "***********************\n";
        for (unsigned int gp=0; gp<material_weights.size(); gp++){
            std::cout << " Gauss Point " << gp << "\n";
            std::cout << "  volume:  " << volume[gp] << "\n";
            std::cout << "  density: " << density[gp] << "\n";
            std::cout << "  C. of mass: "; elib::print(center_of_mass[gp]);
            std::cout << "  local C. of mass: "; elib::print(local_center_of_mass[gp]);
        }
        return 0;
    }

    int MicromorphicFilter::get_cg_phis(elib::vecOfvec &cg_phis){
        /*!
         * Compute the shape function values at the centers of mass.
         */
         cg_phis.resize(local_center_of_mass.size());
         for (unsigned int gp=0; gp<local_center_of_mass.size(); gp++){
             element->get_shape_functions(local_center_of_mass[gp], cg_phis[gp]);
         }
         return 0;
    }

    const unsigned int MicromorphicFilter::id(){
        /*!
         * Return the filter's id number
         */

         return filter_id;
    }

    const unsigned int MicromorphicFilter::dim(){
        /*!
         * Return the filter's spatial dimension
         */
        return filter_dim;
    }

    const std::vector< unsigned int >* MicromorphicFilter::get_element_global_node_ids(){
        /*!
         * Get the global node ids of the underlying Element object.
         */
        return element->get_global_node_ids();
    }

    unsigned int MicromorphicFilter::get_dns_point_gauss_domain(const unsigned int dns_id){
        /*!
         * Get the gauss point associated with the provided dns point
         * returns -1 if the point is not in the filter.
         *
         * :param const unsigned int dns_id: The id of the dns point.
         */

        //Search in the material weights integrator
        for (unsigned int i=0; i<material_weights.size(); i++){
            auto it = material_weights[i].find(dns_id);
            if (it != material_weights[i].end()){
                return i;
            }
        }

        //Search in the dof weights integrator
        for (unsigned int i=0; i<dof_weights.size(); i++){
            auto it = dof_weights[i].find(dns_id);
            if (it != dof_weights[i].end()){
                return i;
            }
        }

        return -1;
    }

    const std::vector< FloatType >* MicromorphicFilter::get_center_of_mass(const unsigned int &gp_id){
        /*!
         * get the center of mass for the gauss point indicated. Returns 
         * NULL if the gp_id is out of range.
         *
         * :param const int &gp_id: The local gauss point number
         */

        if (gp_id > center_of_mass.size()){
            std::cerr << "Error: Gauss point " << gp_id << " is out of range.\n";
            return NULL;
        }
        return &center_of_mass[gp_id];
    }

    const std::string MicromorphicFilter::element_type(){
        /*!
         * Get the name of the element type
         * 
         * :param std::string &eltype: The name of the element type
         */
        
        return element->name;
    }

    int MicromorphicFilter::update_element_node_position(const unsigned int n){
        /*!
         * Update the position of node n using the dof values of the filter.
         *
         * Note: The first filter_dim values are assumed to be the nodal displacements.
         * 
         * :param const unsigned int n: The local node number
         */
        std::vector< FloatType > displacement = std::vector< FloatType >(&dof_values[n][0], &dof_values[n][0]+filter_dim);
        return update_element_node_position(n, displacement);
    }

    int MicromorphicFilter::update_element_node_position(const unsigned int n, const elib::vec &displacement){
        /*!
         * Update the nodal position of node n in the element
         * 
         * :param const unsigned int n: The local node number
         * :param const elib::vec &displacement: The displacement of the node from the reference state
         */

       return element->update_node_position(n, displacement);

//       if (element->reference_nodes[n].size() != displacement.size()){
//           std::cerr << "Error: local node " << n << " has a dimension of " << element->reference_nodes[n].size() << ".\n";
//           std::cerr << "       the nodal displacement has a dimension of " << displacement.size() << ".\n";
//           return 1;
//       }
//       for (unsigned int i=0; i<element->reference_nodes[n].size(); i++){
//           element->nodes[n][i] = element->reference_nodes[n][i] + displacement[i];
//       }
//       return 0;

    }

    int MicromorphicFilter::update_element_node_positions(const elib::vecOfvec &displacements){
        /*!
         * Update the nodal positions of the elements to reflect a movement of the filter.
         * 
         * :param elib::vecOfvec &displacements: The displacements at the nodes from the reference coordinates
         */

         return element->update_node_positions(displacements);

//         if (element->nodes.size() != displacements.size()){
//             std::cerr << "Error: " << displacements.size() << " nodal displacements provided to an element which has " << element->nodes.size() << "\n";
//             return 1;
//         }
//         for (unsigned int n=0; n<element->nodes.size(); n++){
//             int uenp_result = update_element_node_positions(n, displacements[n]);
//             if (uenp_result > 0){
//                 return uenp_result;
//             }
//         }
//         return 0;
    }

    int MicromorphicFilter::update_dof_values(const unsigned int n,
                                              const std::vector< FloatType > &new_dof_values){
        /*!
         * Update the degree of freedom values for node n
         * 
         * :param const unsigned int n: The node to update the degree of freedom values for
         * :param const std::vector< FloatType > &new_dof_values: The new values of the degree of freedom
         */
        dof_values[n] = new_dof_values;
        return 0;
    }

    int MicromorphicFilter::update_dof_values(const vecOfvec &new_dof_values){
        /*!
         * Update the degree of freedom values for the filter.
         * 
         * :param std::vector< double > &new_dof_values: The new values of the degrees of freedom
         */

        if (new_dof_values.size() != element->nodes.size()){
            std::cerr << "Error: new degrees of freedom must be defined for all nodes to use this function.\n";
            std::cerr << "       Individual degrees of freedom can be updated using\n";
            std::cerr << "        MicromorphicFilter::update_dof_values(local_node_number, new_dof_at_node)\n";
            return 1;
        }

        for (unsigned int n=0; n<new_dof_values.size(); n++){
            update_dof_values(n, new_dof_values[n]);
        }
        return 0;
    }

    int MicromorphicFilter::write_to_file(std::ofstream &file){
        /*!
         * Write the filter data to the provided file.
         * 
         * :param std::ofstream &file: The output file
         * :param unsigned int &id: The filter's id number
        */
        
        file << "*MICROMORPHIC FILTER, " << filter_id << "\n";
        file << " *ELEMENT\n";
        file << "  *NODES\n";
        for (unsigned int n=0; n<element->nodes.size(); n++){
            file << "   ";
            for (unsigned int i=0; i<element->nodes[n].size(); i++){
                file << element->nodes[n][i] << ", ";
            }
            file << "\n";
        }
        //Write out the nodal values
        file << " *DOF VALUES\n";
        for (unsigned int n=0; n<dof_values.size(); n++){
            file << "  ";
            for (unsigned int i=0; i<dof_values[n].size(); i++){
                file << dof_values[n][i] << ", ";
            }
            file << "\n";
        }
        //Write out the Gauss point values
        file << " *GAUSS POINT INFORMATION\n";
        for (unsigned int gp=0; gp<material_weights.size(); gp++){
            file << "  *VOLUME, " << volume[gp] << "\n";
            file << "  *SURFACE AREAS (plane, area)\n";
            for (auto it=surface_area[gp].begin(); it!=surface_area[gp].end(); it++){
                file << "   " << it->first << ", " << it->second << "\n";
            }
            file << "  *SURFACE NORMALS (plane, N1, N2, ...)\n";
            for (auto it=surface_normal[gp].begin(); it!=surface_normal[gp].end(); it++){
                file << "   " << it->first;
                for (unsigned int i=0; i<it->second.size(); i++){
                    file << ", " << it->second[i];
                }
                file << "\n";
            }
            file << "  *DENSITY, " << density[gp] << "\n";
            file << "  *LOCAL MASS CENTER, ";
            for (unsigned int i=0; i<local_center_of_mass[gp].size(); i++){
                file << local_center_of_mass[gp][i] << ", ";
            }
            file << "\n";
            file << "  *GLOBAL MASS CENTER, ";
            for (unsigned int i=0; i<center_of_mass[gp].size(); i++){
                file << center_of_mass[gp][i] << ", ";
            }
            file << "\n";
            file << "  *SYMMETRIC MICROSTRESS\n";
            file << "   ";
            for (unsigned int i=0; i<symmetric_microstress[gp].size(); i++){
                if (i==0){
                    file << symmetric_microstress[gp][i];
                }
                else{
                    file << ", " << symmetric_microstress[gp][i];
                }
            }
            file << "\n";
            file << "  *CAUCHY STRESS\n";
            file << "   ";
            for (unsigned int i=0; i<cauchy_stress[gp].size(); i++){
                if (i==0){
                    file << cauchy_stress[gp][i];
                }
                else{
                    file << ", " << cauchy_stress[gp][i];
                }
            }
            file << "\n";
            
        }
        file.flush();
        return 0;
    }

    int MicromorphicFilter::clear_microscale(){
        /*!
         * Clear any stored micro-scale information. This 
         * should be done prior to building a new integrator 
         * after an old one has been built. 
         */

        //dof_id_numbers.clear();
        //material_id_numbers.clear();
        micro_dof_local_coordinates.clear();
        micro_material_local_coordinates.clear();
        dof_weights.clear();
        material_weights.clear();

        //Clear average quantities

        //Mass and geometry
        volume.clear();
        surface_area.clear();
        surface_normal.clear();
        density.clear();
        local_center_of_mass.clear();
        center_of_mass.clear();

        //Stresses
        symmetric_microstress.clear();
        cauchy_stress.clear();
        couple_stress.clear();

        //Tractions
        traction.clear();
        couple_traction.clear();

        //Force vectors
        surface_external_force.clear();
        body_external_force.clear();
        kinetic_force.clear();

        //Force vectors
        surface_external_couple.clear();
        symm_cauchy_couple.clear();
        body_external_couple.clear();
        kinetic_couple.clear();
        return 0;
    }

    std::vector< double > subtract(const std::vector< double > &a, const std::vector< double > &b){
        /*!
         * Subtract two vectors c = a - b
         * 
         * :param std::vector< double > &a: The leading vector
         * :param std::vector< double > &b: The tailing vector
         */
        
        if (a.size() != b.size()){
            std::cerr << "Error: vectors of different sizes cannot be subtracted.\n";
            assert(1==0);
        }

        std::vector< double > c(a.size());
        for (unsigned int i=0; i<a.size(); i++){
            c[i] = a[i] - b[i];
        }
        return c;
    }

    bool point_on_surface(const std::vector< double > &p, const std::vector< double > &n, const std::vector< double > &a){
        /*!
         * Determine whether the point p is on the surface defined by the normal n and the point on the surface a.
         * 
         * :param const std::vector< double > &p: The query point
         * :param const std::vector< double > &n: The surface normal
         * :param const std::vector< double > &a: A point on the surface
         */
        
        std::vector< double > d = subtract(p, a);
        double distance = dot(d, n);
        return fuzzy_equals(distance, 0);
    }

//    void construct_cauchy_least_squares(const std::vector< std::map< unsigned int, std::vector< double > > > &surface_normals,
//                                        const std::vector< std::map< unsigned int, std::vector< double > > > &surface_tractions,
//                                        Eigen::MatrixXd &A, Eigen::MatrixXd &b){
//        /*!
//         * Construct the normal matrix which can project the Cauchy stresses at the Gauss points to the traction vector b on each of the 
//         * surfaces. Assumes a 3D Cauchy stress.
//         * 
//         * :param const std::vector< std::map< unsigned int, std::vector< double > > > &surface_normals: The vector of maps from the 
//         *     gauss domain's face number to the normal of that face.
//         * :param const std::vector< std::map< unsigned int, std::vector< double > > > &surface_tractions: The vector of maps from the 
//         *     gauss domain's face number to the surface traction on that face.
//         * :param Eigen::MatrixXd &A: The normal matrix
//         * :param Eigen::MatrixXd &b: The traction vector
//         */
//        
//        //Determine the size of the A matrix
//        unsigned int nrows, ncols;
//        unsigned int nstress = 9;
//        unsigned int dim = 3;
//        nrows = 0;
//        ncols = nstress*surface_normals.size(); //9 components of the non-symmetric cauchy stress for each gauss point
//        for (unsigned int gp=0; gp<surface_normals.size(); gp++){
//            nrows += dim*surface_normals[gp].size();
//
//        }
//
//        if (surface_tractions.size() != surface_normals.size()){
//            std::cerr << "Error: surface_tractions should have the same size as surface_normals\n";
//            std::cerr << "       surface_normals.size(): " << surface_normals.size() << "\n";
//            std::cerr << "       surface_tractions.size(): " << surface_tractions.size() << "\n";
//            assert(1==0);
//        }
//
//        //Resize A
//        A = Eigen::MatrixXd::Zero(nrows, ncols);
//        b = Eigen::MatrixXd::Zero(nrows, 1);
//
//        //Iterate over the gauss points
//        unsigned int row0, col0;
//        row0 = col0 = 0;
//
//        for (unsigned int gp=0; gp<surface_normals.size(); gp++){
//
//            for (auto face = surface_normals[gp].begin(); face != surface_normals[gp].end(); face++){
//                //Set the values in the A matrix
//                A(row0 + 0, col0 + 0) = face->second[0];
//                A(row0 + 0, col0 + 7) = face->second[2];
//                A(row0 + 0, col0 + 8) = face->second[1];
//                A(row0 + 1, col0 + 1) = face->second[1];
//                A(row0 + 1, col0 + 5) = face->second[0];
//                A(row0 + 1, col0 + 6) = face->second[2];
//                A(row0 + 2, col0 + 2) = face->second[2];
//                A(row0 + 2, col0 + 3) = face->second[1];
//                A(row0 + 2, col0 + 4) = face->second[0];
//
//                auto traction = surface_tractions[gp].find(face->first);
//                if (traction == surface_tractions[gp].end()){
//                    std::cerr << "Error: surface traction for face " << face->first << " not found\n";
//                    assert(1==0);
//                }
//
//                b(row0 + 0, 0) = traction->second[0];
//                b(row0 + 1, 0) = traction->second[1];
//                b(row0 + 2, 0) = traction->second[2];
//
//                row0 += dim; //Increment row0
//            }
//            col0 += nstress; //Increment col0
//        }
//
//    }

    void construct_linear_momentum_surface_external_force(
             const std::vector< std::map< unsigned int, std::vector< double > > > &face_shapefunctions,
             const std::vector< std::map< unsigned int, std::vector< double > > > &face_tractions,
             const std::vector< std::map< unsigned int, double > > &face_areas,
             const std::vector< std::vector< unsigned int > > &external_face_ids,
             std::vector< double > &surface_external_force){
        /*!
         * Construct the surface force acting on the nodes of the finite element
         * 
         * :param const std::vector< std::map< unsigned int, double > > &face_shapefunctions: The shapefunctions at the 
         *     gauss domain centroids.
         * :param const std::vector< std::map< unsigned int, std::vector< double > > > &face_tractions: The tractions on the 
         *     gauss domain faces.
         * :param const std::vector< std::map< unsigned int, double > > &face_areas: The areas of each of the gauss domain faces in 
         *     global coordinates.
         * :param const std::vector< std::vector< unsigned int > > &external_face_ids: The id numbers of faces which are the boundaries 
         *     of the filter.
         * :param std::vector< double > &surface_external_force: The external surface force vector
         */

        unsigned int dim = 3;
        surface_external_force.clear();

        //Make sure the incoming vectors have a non-zero size
        unsigned int ngpts = face_shapefunctions.size();
        if (ngpts == 0){
            std::cerr << "Error: no gauss points defined in face_shapefunctions\n";
            assert(1==0);
        }
        if (ngpts != face_tractions.size()){
            std::cerr << "Error: face_tractions doesn't have as many gauss points as face_shapefunctions\n";
            std::cerr << "       face_tractions.size(): " << face_tractions.size() << "\n";
            std::cerr << "       face_shapefunctions.size(): " << ngpts << "\n";
            assert(1==0);
        }
        if (ngpts != face_areas.size()){
            std::cerr << "Error: face_areas doesn't have as many gauss points as face_shapefunctions\n";
            std::cerr << "       face_areas.size(): " << face_areas.size() << "\n";
            std::cerr << "       face_shapefunctions.size(): " << ngpts << "\n";
            assert(1==0);
        }

        //Initialize the surface external force vector.
        surface_external_force = std::vector< double >(dim*face_shapefunctions[0].begin()->second.size(), 0);

        //Loop over the gauss points
        for (unsigned int gp=0; gp<face_shapefunctions.size(); gp++){
            //Loop over the faces
            for (auto face=face_shapefunctions[gp].begin(); face!=face_shapefunctions[gp].end(); face++){
                //check if the face is an external face
                auto extface = std::find(external_face_ids[gp].begin(), external_face_ids[gp].end(), face->first);
                if (extface == external_face_ids[gp].end()){
                    continue;
                }

                auto traction = face_tractions[gp].find(face->first);
                auto area = face_areas[gp].find(face->first);

                if (traction == face_tractions[gp].end()){
                    std::cerr << "Error: Face " << face->first << " not found in tractions.\n";
                    assert(1==0);
                }
                if (area == face_areas[gp].end()){
                    std::cerr << "Error: Face " << face->first << " not found in areas.\n";
                    assert(1==0);
                }

                //Iterate through the element's nodes
                for (unsigned int n=0; n<face->second.size(); n++){
                    //Loop over the traction's indices
                    for (unsigned int i=0; i<traction->second.size(); i++){
                        //N_node(x)*traction(x)*area
                        surface_external_force[dim*n + i] += face->second[n]*traction->second[i]*area->second;
                    }
                }
            }
        }
        return;
    }

    void construct_first_moment_surface_external_couple(
             const std::vector< std::map< unsigned int, std::vector< double > > > &face_shapefunctions,
             const std::vector< std::map< unsigned int, std::vector< double > > > &face_couples,
             const std::vector< std::map< unsigned int, double > > &face_areas,
             const std::vector< std::vector< unsigned int > > &external_face_ids,
             std::vector< double > &surface_external_couple){
        /*!
         * Construct the surface force acting on the nodes of the finite element
         * 
         * :param const std::vector< std::map< unsigned int, double > > &face_shapefunctions: The shapefunctions at the 
         *     gauss domain centroids.
         * :param const std::vector< std::map< unsigned int, std::vector< double > > > &face_couples: The couple tractions on the 
         *     gauss domain faces.
         * :param const std::vector< std::map< unsigned int, double > > &face_areas: The areas of each of the gauss domain faces in 
         *     global coordinates.
         * :param const std::vector< std::vector< unsigned int > > &external_face_ids: The id numbers of faces which are the boundaries 
         *     of the filter.
         * :param std::vector< double > &surface_external_force: The external surface couple vector
         */

        unsigned int ncouple = 9;
        surface_external_couple.clear();

        //Make sure the incoming vectors have a non-zero size
        unsigned int ngpts = face_shapefunctions.size();
        if (ngpts == 0){
            std::cerr << "Error: no gauss points defined in face_shapefunctions\n";
            assert(1==0);
        }
        if (ngpts != face_couples.size()){
            std::cerr << "Error: face_couples doesn't have as many gauss points as face_shapefunctions\n";
            std::cerr << "       face_couples.size(): " << face_couples.size() << "\n";
            std::cerr << "       face_shapefunctions.size(): " << ngpts << "\n";
            assert(1==0);
        }
        if (ngpts != face_areas.size()){
            std::cerr << "Error: face_areas doesn't have as many gauss points as face_shapefunctions\n";
            std::cerr << "       face_areas.size(): " << face_areas.size() << "\n";
            std::cerr << "       face_shapefunctions.size(): " << ngpts << "\n";
            assert(1==0);
        }

        //Initialize the surface external force vector.
        surface_external_couple = std::vector< double >(ncouple*face_shapefunctions[0].begin()->second.size(), 0);

        //Loop over the gauss points
        for (unsigned int gp=0; gp<face_shapefunctions.size(); gp++){
            //Loop over the faces
            for (auto face=face_shapefunctions[gp].begin(); face!=face_shapefunctions[gp].end(); face++){
                //check if the face is an external face
                auto extface = std::find(external_face_ids[gp].begin(), external_face_ids[gp].end(), face->first);
                if (extface == external_face_ids[gp].end()){
                    continue;
                }

                auto couple = face_couples[gp].find(face->first);
                auto area = face_areas[gp].find(face->first);

                if (couple == face_couples[gp].end()){
                    std::cerr << "Error: Face " << face->first << " not found in couples.\n";
                    assert(1==0);
                }
                if (area == face_areas[gp].end()){
                    std::cerr << "Error: Face " << face->first << " not found in areas.\n";
                    assert(1==0);
                }

                //Iterate through the element's nodes
                for (unsigned int n=0; n<face->second.size(); n++){

                    //Assign the couple's indices
                    for (unsigned int i=0; i<ncouple; i++){
                        surface_external_couple[ncouple*n + i] += face->second[n]*couple->second[i]*area->second;
                    }
                }
            }
        }
        return;
    }

    void construct_first_moment_symm_cauchy_couple(const vecOfvec &com_shapefunctions,
                                                   const vecOfvec &symmetric_microstress, const vecOfvec &cauchy_stress,
                                                   const std::vector< double > &volume,
                                                   std::vector< double > &symm_cauchy_couple){
        /*!
         * Compute the couple resulting from the difference between the Cauchy and symmetric microstress
         * 
         * :param const vecOfvec &com_shapefunctions: The shapefunction values at the centers of mass
         * :param const vecOfvec &symmetric_microstress: The values of the symmetric micro-stress at the gauss points
         * :param const vecOfvec &cauchy_stress: The values of the cauchy stress at the gauss points
         * :param const std::vector< double > &volume: The volume associated with the gauss point.
         * :param std::vector< double > &symm_cauchy_couple: The couple resulting from the difference between 
         *      the Cauchy and symmetric microstress.
         */

        unsigned int ncouple = 9;
        symm_cauchy_couple.clear();

        //Make sure the incoming vectors have a non-zero size
        unsigned int ngpts = com_shapefunctions.size();
        if (ngpts == 0){
            std::cerr << "Error: no gauss points defined in com_shapefunctions\n";
            assert(1==0);
        }
        if (ngpts != symmetric_microstress.size()){
            std::cerr << "Error: symmetric_microstress doesn't have as many gauss points as com_shapefunctions\n";
            std::cerr << "       symmetric_microstress.size(): " << symmetric_microstress.size() << "\n";
            std::cerr << "       com_shapefunctions.size(): " << ngpts << "\n";
            assert(1==0);
        }
        if (ngpts != cauchy_stress.size()){
            std::cerr << "Error: cauchy_stress doesn't have as many gauss points as com_shapefunctions\n";
            std::cerr << "       cauchy_stress.size(): " << cauchy_stress.size() << "\n";
            std::cerr << "       com_shapefunctions.size(): " << ngpts << "\n";
            assert(1==0);
        }
        if (ngpts != volume.size()){
            std::cerr << "Error: volume doesn't have as many gauss points as com_shapefunctions\n";
            std::cerr << "       volume.size(): " << volume.size() << "\n";
            std::cerr << "       com_shapefunctions.size(): " << ngpts << "\n";
            assert(1==0);
        }

        //Initialize symm_cauchy_couple
        symm_cauchy_couple = std::vector< double >(com_shapefunctions[0].size()*ncouple, 0);

        //Iterate through the gauss points
        for (unsigned int gp=0; gp<ngpts; gp++){

            //Iterate through the element's nodes
            for (unsigned int n=0; n<com_shapefunctions[gp].size(); n++){
                //Assign the couple's indices. Require transpose of Voigt so doing this explicitly
                symm_cauchy_couple[ncouple*n + 0] += com_shapefunctions[gp][n]*(cauchy_stress[gp][0] - symmetric_microstress[gp][0]);
                symm_cauchy_couple[ncouple*n + 1] += com_shapefunctions[gp][n]*(cauchy_stress[gp][1] - symmetric_microstress[gp][1]);
                symm_cauchy_couple[ncouple*n + 2] += com_shapefunctions[gp][n]*(cauchy_stress[gp][2] - symmetric_microstress[gp][2]);
                symm_cauchy_couple[ncouple*n + 3] += com_shapefunctions[gp][n]*(cauchy_stress[gp][6] - symmetric_microstress[gp][6]);
                symm_cauchy_couple[ncouple*n + 4] += com_shapefunctions[gp][n]*(cauchy_stress[gp][7] - symmetric_microstress[gp][7]);
                symm_cauchy_couple[ncouple*n + 5] += com_shapefunctions[gp][n]*(cauchy_stress[gp][8] - symmetric_microstress[gp][8]);
                symm_cauchy_couple[ncouple*n + 6] += com_shapefunctions[gp][n]*(cauchy_stress[gp][3] - symmetric_microstress[gp][3]);
                symm_cauchy_couple[ncouple*n + 7] += com_shapefunctions[gp][n]*(cauchy_stress[gp][4] - symmetric_microstress[gp][4]);
                symm_cauchy_couple[ncouple*n + 8] += com_shapefunctions[gp][n]*(cauchy_stress[gp][5] - symmetric_microstress[gp][5]);
             }

        }

    }

    void construct_linear_momentum_least_squares_matrix(const std::vector< vecOfvec > &cg_shapefunction_gradients,
                                                        const std::vector< FloatType > &volume,
                                                        const std::vector< vecOfvec > &vertex_cauchy_stress,
                                                        Eigen::MatrixXd &C){
        /*!
         * Construct the linear momentum least squares matrix i.e. the divergence of the Cauchy stress
         * 
         * :param const std::vector< vecOfvec > &cg_shapefunction_gradients: The global gradients of the shapefunctions 
         *     at the gauss points.
         * :param const std::vector< FloatType > &volume: The volume of the gauss domains
         * :param const std::vector< vecOfvec > &vertex_cauchy_stress: The cauchy stress at the vertices of the gauss domains
         * :param Eigen::MatrixXd &C: The constraint matrix.
         */
        
        //Assume the problem is 3D
        unsigned int dim=3;
        unsigned int nstress=9;

        //Make sure at least one gauss domain is defined
        unsigned int ngpts = cg_shapefunction_gradients.size();
        if (ngpts == 0){
            std::cerr << "Error: At least one Gauss point must be defined.\n";
            assert(1==0);
        }

        if (ngpts != volume.size()){
            std::cerr << "Error: The number of shapefunction gradients and volumes is not equal.\n";
            std::cerr << "       cg_shapefunction_gradients.size(): " << ngpts << "\n";
            std::cerr << "       volume.size(): " << volume.size() << "\n";
            assert(1==0);
        }

        //Initialize the C matrix
        C = Eigen::MatrixXd::Zero(dim*cg_shapefunction_gradients[0].size(), nstress*cg_shapefunction_gradients.size());

        //Initialize the C2 matrix
        unsigned int a2cols = 0;
        for (unsigned int gp=0; gp<vertex_cauchy_stress.size(); gp++){
            a2cols += vertex_cauchy_stress[gp].size();
        }
        Eigen::MatrixXd C2 = Eigen::MatrixXd::Zero(C.cols(), a2cols);


        //Loop over the gauss domains
        for (unsigned int gp=0; gp<cg_shapefunction_gradients.size(); gp++){
            //Loop over the nodes
            for (unsigned int n=0; n<cg_shapefunction_gradients[gp].size(); n++){
                C(dim*n+0, nstress*gp + 0) += cg_shapefunction_gradients[gp][n][0]*volume[gp];
                C(dim*n+0, nstress*gp + 7) += cg_shapefunction_gradients[gp][n][2]*volume[gp];
                C(dim*n+0, nstress*gp + 8) += cg_shapefunction_gradients[gp][n][1]*volume[gp];
                C(dim*n+1, nstress*gp + 1) += cg_shapefunction_gradients[gp][n][1]*volume[gp];
                C(dim*n+1, nstress*gp + 5) += cg_shapefunction_gradients[gp][n][0]*volume[gp];
                C(dim*n+1, nstress*gp + 6) += cg_shapefunction_gradients[gp][n][2]*volume[gp];
                C(dim*n+2, nstress*gp + 2) += cg_shapefunction_gradients[gp][n][2]*volume[gp];
                C(dim*n+2, nstress*gp + 3) += cg_shapefunction_gradients[gp][n][1]*volume[gp];
                C(dim*n+2, nstress*gp + 4) += cg_shapefunction_gradients[gp][n][0]*volume[gp];
            }

        }

        //Construct C2
        unsigned int row0, col0;
        row0 = col0 = 0;
        for (unsigned int gp=0; gp<vertex_cauchy_stress.size(); gp++){

            for (unsigned int v=0; v<vertex_cauchy_stress[gp].size(); v++){

                for (unsigned int i=0; i<vertex_cauchy_stress[gp][v].size(); i++){
                    C2(row0+i, col0+v) = vertex_cauchy_stress[gp][v][i];
                }

            }

            row0 += nstress;
            col0 += vertex_cauchy_stress[gp].size(); //Increment col0

        }

        //Incorporate the 
        C *= C2;

        return;
    }

    void construct_first_moment_constraint_matrix(const std::vector< vecOfvec > &cg_shapefunction_gradients,
                                                  const std::vector< FloatType > &volume,
                                                  Eigen::MatrixXd &C){
        /*!
         * Construct the first moment of momentum constraint matrix i.e. the divergence of the couple stress
         * 
         * :param const std::vector< vecOfvec > &cg_shapefunction_gradients: The global gradients of the shapefunctions 
         *     at the gauss points.
         * :param const std::vector< FloatType > &volume: The volume of the gauss domains
         * :param Eigen::MatrixXd &C: The constraint matrix.
         */
        
        //Assume the problem is 3D
        unsigned int ncouple=9;
        unsigned int nstress=27;

        //Make sure at least one gauss domain is defined
        unsigned int ngpts = cg_shapefunction_gradients.size();
        if (ngpts == 0){
            std::cerr << "Error: At least one Gauss point must be defined.\n";
            assert(1==0);
        }

        if (ngpts != volume.size()){
            std::cerr << "Error: The number of shapefunction gradients and volumes is not equal.\n";
            std::cerr << "       cg_shapefunction_gradients.size(): " << ngpts << "\n";
            std::cerr << "       volume.size(): " << volume.size() << "\n";
            assert(1==0);
        }

        //Initialize the C matrix
        C = Eigen::MatrixXd::Zero(ncouple*cg_shapefunction_gradients[0].size(), nstress*cg_shapefunction_gradients.size());

        //Loop over the gauss domains
        for (unsigned int gp=0; gp<ngpts; gp++){
            //Loop over the nodes
            for (unsigned int n=0; n<cg_shapefunction_gradients[gp].size(); n++){

                C(ncouple*n + 0, nstress*gp +  0) = cg_shapefunction_gradients[gp][n][0]*volume[gp];
                C(ncouple*n + 0, nstress*gp +  9) = cg_shapefunction_gradients[gp][n][1]*volume[gp];
                C(ncouple*n + 0, nstress*gp + 18) = cg_shapefunction_gradients[gp][n][2]*volume[gp];
                C(ncouple*n + 1, nstress*gp +  1) = cg_shapefunction_gradients[gp][n][0]*volume[gp];
                C(ncouple*n + 1, nstress*gp + 10) = cg_shapefunction_gradients[gp][n][1]*volume[gp];
                C(ncouple*n + 1, nstress*gp + 19) = cg_shapefunction_gradients[gp][n][2]*volume[gp];
                C(ncouple*n + 2, nstress*gp +  2) = cg_shapefunction_gradients[gp][n][0]*volume[gp];
                C(ncouple*n + 2, nstress*gp + 11) = cg_shapefunction_gradients[gp][n][1]*volume[gp];
                C(ncouple*n + 2, nstress*gp + 20) = cg_shapefunction_gradients[gp][n][2]*volume[gp];
                C(ncouple*n + 3, nstress*gp +  3) = cg_shapefunction_gradients[gp][n][0]*volume[gp];
                C(ncouple*n + 3, nstress*gp + 12) = cg_shapefunction_gradients[gp][n][1]*volume[gp];
                C(ncouple*n + 3, nstress*gp + 21) = cg_shapefunction_gradients[gp][n][2]*volume[gp];
                C(ncouple*n + 4, nstress*gp +  4) = cg_shapefunction_gradients[gp][n][0]*volume[gp];
                C(ncouple*n + 4, nstress*gp + 13) = cg_shapefunction_gradients[gp][n][1]*volume[gp];
                C(ncouple*n + 4, nstress*gp + 22) = cg_shapefunction_gradients[gp][n][2]*volume[gp];
                C(ncouple*n + 5, nstress*gp +  5) = cg_shapefunction_gradients[gp][n][0]*volume[gp];
                C(ncouple*n + 5, nstress*gp + 14) = cg_shapefunction_gradients[gp][n][1]*volume[gp];
                C(ncouple*n + 5, nstress*gp + 23) = cg_shapefunction_gradients[gp][n][2]*volume[gp];
                C(ncouple*n + 6, nstress*gp +  6) = cg_shapefunction_gradients[gp][n][0]*volume[gp];
                C(ncouple*n + 6, nstress*gp + 15) = cg_shapefunction_gradients[gp][n][1]*volume[gp];
                C(ncouple*n + 6, nstress*gp + 24) = cg_shapefunction_gradients[gp][n][2]*volume[gp];
                C(ncouple*n + 7, nstress*gp +  7) = cg_shapefunction_gradients[gp][n][0]*volume[gp];
                C(ncouple*n + 7, nstress*gp + 16) = cg_shapefunction_gradients[gp][n][1]*volume[gp];
                C(ncouple*n + 7, nstress*gp + 25) = cg_shapefunction_gradients[gp][n][2]*volume[gp];
                C(ncouple*n + 8, nstress*gp +  8) = cg_shapefunction_gradients[gp][n][0]*volume[gp];
                C(ncouple*n + 8, nstress*gp + 17) = cg_shapefunction_gradients[gp][n][1]*volume[gp];
                C(ncouple*n + 8, nstress*gp + 26) = cg_shapefunction_gradients[gp][n][2]*volume[gp];
            }
        }

        return;
    }

    void solve_constrained_least_squares(const Eigen::MatrixXd &A, const Eigen::MatrixXd &b,
                                         const Eigen::MatrixXd &C, const Eigen::MatrixXd &d,
                                         Eigen::MatrixXd &x){
        /*!
         * Solve a least squares problem to minimize ||Ax - b|| while subject to the equality 
         * constraint Cx = d
         * 
         * :param const Eigen::MatrixXd &A: The least-squares matrix (nequations, nvariables)
         * :param const Eigen::MatrixXd &b: The solution vector for the least-squares matrix (nvariables, 1)
         * :param const Eigen::MatrixXd &C: The constraint matrix (nconstraints, nvariables)
         * :param const Eigen::MatrixXd &d: The solution vector for the constraint matrix (nvariables, 1)
         * :param Eigen::MatrixXd &x: The optimal solution for the free variables.
         */

        unsigned int nvariables = A.cols();
        unsigned int nconstraints = C.rows();

        if (nconstraints > nvariables){
            std::cerr << "Error: more constraints than variables. Least squares should be\n";
            std::cerr << "       performed on the constraint matrix.\n";
            assert(1==0);
        }

        //Construct the M matrix
        Eigen::MatrixXd M = Eigen::MatrixXd::Zero(nvariables+nconstraints, nvariables+nconstraints);
        M.block(0, 0, nvariables, nvariables) = 2*A.transpose()*A;
        M.block(0, nvariables, nvariables, nconstraints) = C.transpose();
        M.block(nvariables, 0, nconstraints, nvariables) = C;

        //Form the RHS vector
        Eigen::MatrixXd RHS = Eigen::MatrixXd::Zero(nvariables + nconstraints, 1);
        RHS.block(0, 0, nvariables, 1) = 2*A.transpose()*b;
        RHS.block(nvariables, 0, nconstraints, 1) = d;

        //Solve the system
        Eigen::MatrixXd solution = M.colPivHouseholderQr().solve(RHS);
        double relative_error = (M*solution - RHS).norm()/std::max((M*solution).norm(), RHS.norm());

        if (relative_error > 1e-6){
              double numerator = (C*solution.block(0, 0, nvariables, 1) - d).norm();
              double T1norm = (C*solution.block(0, 0, nvariables, 1)).norm();
              double T2norm = (d.norm());
              double constraint_error = numerator/std::max(T1norm, T2norm);
//            double constraint_error = (C*solution.block(nvariables, 0, nconstraints, 1) - d).norm()/std::max((C*solution.block(nvariables, 0, nconstraints, 1)).norm(), d.norm());
            std::cerr << "Warning: Relative error is larger than threshold.\n";
            std::cerr << "         Relative error  : " << relative_error << "\n";
            std::cerr << "         Constraint error: " << constraint_error << "\n";
        }

        //Return the variables
        x = solution.block(0, 0, nvariables, 1);
        return;
    }

    void construct_linear_momentum_b_vector(const unsigned int nconstraints,
                                            const std::vector< FloatType > &surface_external_force,
                                            const std::vector< FloatType > &body_external_force,
                                            const std::vector< FloatType > &kinetic_force,
                                            Eigen::MatrixXd &d){
        /*!
         * Construct the b vector for contrained least squares given the equation
         * \int_{\Omega} N_{,j} sigma_{ji} dv = \int_{\partial \Omega} N n_{j} \sigma_{ji} da - \int_{\Omega} \rho\left(b_i - a_i\right) dv\\
         * :param const unsigned int nconstraints: The number of constraint equations.
         * :param const std::vector< FloatType > &surface_external_force: The external tractions acting on the body.
         * :param const std::vector< FloatType > &body_external_force: The body force external forces acting on the body.
         * :param const std::vector< FloatType > &kinetic_force: The kinetic force acting on the body.
         * :param Eigen::MatrixXd &d: The b vector to be constructed.
         */

        d = Eigen::MatrixXd::Zero(nconstraints, 1);

        if (surface_external_force.size() != nconstraints){
            std::cout << "Error: the external tractions on the body and the constraint matrix must be defined.\n";
            assert(1==0);
        }

        //Add the surface external force to d
        d += Eigen::Map< const EigVec >(surface_external_force.data(), nconstraints, 1);

        //Add the body external force to d
        if (body_external_force.size() == nconstraints){
            d += Eigen::Map< const EigVec >(body_external_force.data(), nconstraints, 1);
        }
        else if (body_external_force.size() > 0){
            std::cout << "Error: The external body force vector is not the same size as the surface traction vector.\n";
            assert(1==0);
        }

        //Subtract the kinetic force from d
        if (kinetic_force.size() == nconstraints){
            d -= Eigen::Map< const EigVec>(kinetic_force.data(), nconstraints, 1);
        }
        else if (kinetic_force.size() > 0){
            std::cout << "Error: The kinetic force vector is not the same size as the surface traction vector.\n";
            assert(1==0);
        }

        return;
    }

    void construct_first_moment_d_vector(const unsigned int nconstraints,
                                         const std::vector< FloatType > &surface_external_couple,
                                         const std::vector< FloatType > &symm_cauchy_couple,
                                         const std::vector< FloatType > &body_external_couple,
                                         const std::vector< FloatType > &kinetic_couple,
                                         Eigen::MatrixXd &d){
        /*!
         * Construct the d vector for contrained least squares given the equation
         * \int_{\Omega} N_{,k} m_{kij} dv = \int_{\partial \Omega} N n_{k} m_{kij} da + \int_{Omega} N (\sigma_{ji} - s_{ji}\right) dv + \int_{\Omega} \rho\left(l_{ij} - \omega_{ij}\right) dv\\
         * :param const unsigned int nconstraints: The number of constraint equations.
         * :param const std::vector< FloatType > &surface_external_couple: The external couples acting on the body.
         * :param const std::vector< FloatType > &symm_cauchy_couple: The couple resulting from the difference between the Cauchy
         *     and symmetric microstress.
         * :param const std::vector< FloatType > &body_external_couple: The body couple external forces acting on the body.
         * :param const std::vector< FloatType > &kinetic_couple: The kinetic couple acting on the body.
         * :param Eigen::MatrixXd &d: The d vector to be constructed.
         */

        d = Eigen::MatrixXd::Zero(nconstraints, 1);

        if (surface_external_couple.size() != nconstraints){
            std::cout << "Error: the external couples on the body and the constraint matrix must be defined.\n";
            assert(1==0);
        }
        if (symm_cauchy_couple.size() != nconstraints){
            std::cout << "Error: the cauchy - symmetric microstress couples on the body and the constraint matrix must be defined.\n";
            assert(1==0);
        }

        //Add the surface external force to d
        d += Eigen::Map< const EigVec >(surface_external_couple.data(), nconstraints, 1);

        //Add the cauchy-symmetric microstress couple
        d += Eigen::Map< const EigVec >(symm_cauchy_couple.data(), nconstraints, 1);

        //Add the body external couple to d
        if (body_external_couple.size() == nconstraints){
            assert(1==2);
            d += Eigen::Map< const EigVec >(body_external_couple.data(), nconstraints, 1);
        }
        else if (body_external_couple.size() > 0){
            std::cout << "Error: The external body couple vector is not the same size as the surface external couple vector.\n";
            assert(1==0);
        }

        //Subtract the kinetic couple from d
        if (kinetic_couple.size() == nconstraints){
            assert(2==3);
            d -= Eigen::Map< const EigVec>(kinetic_couple.data(), nconstraints, 1);
        }
        else if (kinetic_couple.size() > 0){
            std::cout << "Error: The kinetic couple vector is not the same size as the surface external couple vector.\n";
            assert(1==0);
        }

        return;
    }

    void id_unique_vectors(const std::map< unsigned int, std::vector< double > > &vectors,
                           std::map< unsigned int, std::vector< double > > &unique,
                           double tolr, double tola, bool opposite_is_unique){
        /*!
         * Identify a subset of the incoming vectors that are in unique directions. Note that 
         * if opposite_is_unique is false (the default) vectors in opposite directions are id'd 
         * as equivalent since they are linear combinations of each-other.
         * 
         * :param const std::map< unsigned int, std::vector< double > > &vectors: The set of vectors to compare.
         * :param std::map< unsigned int, std::vector< double > > &unique: A set of unique vectors
         * :param double tolr: The relative tolerance
         * :param double tola: The absolute tolerance
         * :param bool opposite_is_unique: If false, vectors facing in exactly opposite directions are not unique.
         */
        
        unique.clear();
        if (vectors.size() == 0){
            return;
        }

        bool is_unique;
        for (auto vector=vectors.begin(); vector!=vectors.end(); vector++){
            is_unique = true;

            for (auto u = unique.begin(); u!=unique.end(); u++){
//                std::cout << "vector: "; print_vector(vector->second);
//                std::cout << "u: "; print_vector(u->second);
//                std::cout << "opposite_is_unique: " << opposite_is_unique << "\n";
                if (compare_vector_directions(vector->second, u->second, tolr, tola, opposite_is_unique)){
                    is_unique=false;
                    break;
                }
            }
            if (is_unique){
                unique.emplace(vector->first, vector->second);
            }
        }
        return;
    }

    void compute_vertex_cauchy_stress(const vecOfvec &normals, const vecOfvec &tractions, std::vector< double > &cauchy_stress){
        /*!
         * Compute the Cauchy stress from the traction associated with the provided normals
         * 
         * :param const vecOfvec &normals: The normal vectors
         * :param const vecOfvec &traction: The surface tractions
         * :param std::vector< double > &cauchy_stress: The re-constructed Cauchy stress 
         */

        //Construct the A matrix and b vector
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3*tractions.size(), 9);
        Eigen::MatrixXd b(3*tractions.size(), 1);

        //Set up a map for the cauchy stress
        cauchy_stress.resize(9);
        Eigen::Map< Eigen::MatrixXd > x(cauchy_stress.data(), 9, 1);

        unsigned int row0=0;
        for (auto n=normals.begin(); n!=normals.end(); n++){
            A(row0+0, 0) = (*n)[0];
            A(row0+0, 7) = (*n)[2];
            A(row0+0, 8) = (*n)[1];
            A(row0+1, 1) = (*n)[1];
            A(row0+1, 5) = (*n)[0];
            A(row0+1, 6) = (*n)[2];
            A(row0+2, 2) = (*n)[2];
            A(row0+2, 3) = (*n)[1];
            A(row0+2, 4) = (*n)[0];

            b(row0+0, 0) = tractions[row0/3][0];
            b(row0+1, 0) = tractions[row0/3][1];
            b(row0+2, 0) = tractions[row0/3][2];

            row0 += 3;
        }

        //Perform the least squares solution
        x = A.colPivHouseholderQr().solve(b);

        return;
    }

    void compute_vertex_couple_stress(const vecOfvec &normals, const vecOfvec &couples, std::vector< double > &couple_stress){
        /*!
         * Compute the higher order stress from the traction associated with the provided normals
         * 
         * :param const vecOfvec &normals: The normal vectors
         * :param const vecOfvec &couples: The surface couple tractions
         * :param std::vector< double > &couple_stress: The re-constructed higher order stress 
         */

        //Construct the A matrix and b vector
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(9*couples.size(), 27);
        Eigen::MatrixXd b(9*couples.size(), 1);

        //Set up a map for the couple stress
        couple_stress.resize(27);
        Eigen::Map< Eigen::MatrixXd > x(couple_stress.data(), couple_stress.size(), 1);

        unsigned int row0=0;
        for (auto n=normals.begin(); n!=normals.end(); n++){
            for (unsigned int i=0; i<9; i++){

                A(row0+i, 0+i) = (*n)[0];
                A(row0+i, 9+i) = (*n)[1];
                A(row0+i,18+i) = (*n)[2];

                b(row0+i, 0) = couples[row0/9][i];
            }
            row0 += 9;
        }

        //Perform the least squares solution
        x = A.colPivHouseholderQr().solve(b);

        return;
    }


    void compute_vertices_cauchy_stress(const std::vector< std::vector< unsigned int > > &vertex_planes,
                                        const MicroPoint &domain,
                                        const std::map< unsigned int, std::vector< FloatType > > &tractions,
                                        vecOfvec &vertex_cauchy){
        /*!
         * Compute the cauchy stress at a series of vertices and their associated planes
         * 
         * :param std::vector< std::vector< unsigned int > > &vertex_planes: The id numbers of the planes associated with the vertex
         * :param std::map< unsigned int, std::vector< FloatType > > &planes: The normal - point definitions of the planes
         * :param std::map< unsigned int, std::vector< FloatType > > &tractions: The tractions associated with each plane
         * :param vecOvec &vertex_cauchy: The cauchy stress at the vertices
         */

        vecOfvec vertex_normals; //The normals associated with the vertices
        vecOfvec vertex_tractions; //The tractions associated with the vertices
        
        //Resize the Cauchy stress vector
        vertex_cauchy.resize(vertex_planes.size());

        //Iterate through the vertices
        for (unsigned int v=0; v<vertex_planes.size(); v++){

            //Assemble the normals and tractions
            vertex_normals.clear(); 
            vertex_normals.resize(vertex_planes[v].size());

            vertex_tractions.clear();
            vertex_tractions.resize(vertex_planes[v].size());

            //Assemble the tractions and normals
            for (unsigned int f=0; f<vertex_planes[v].size(); f++){

                auto face = std::find(domain.planes.begin(), domain.planes.end(), vertex_planes[v][f]);
                if (face == domain.planes.end()){
                    std::cerr << "Error: vertex plane not found in element_planes\n";
                    assert(1==0);
                }

                auto f_traction = tractions.find(vertex_planes[v][f]);
                if (f_traction == tractions.end()){
                    std::cerr << "Error: traction not found in traction\n";
                    assert(1==0);
                }

                //Save the normal
                vertex_normals[f] = domain.normal( std::distance( domain.planes.begin(), face) );

                //Save the traction
                vertex_tractions[f] = f_traction->second;
            }

            compute_vertex_cauchy_stress(vertex_normals, vertex_tractions, vertex_cauchy[v]);
        }
        return;
    }

    void compute_vertices_couple_stress(const std::vector< std::vector< unsigned int > > &vertex_planes,
                                        const MicroPoint &domain,
                                        const std::map< unsigned int, std::vector< FloatType > > &couple_tractions,
                                        vecOfvec &vertex_hostress){
        /*!
         * Compute the higher order stress at a series of vertices and their associated planes
         * 
         * :param std::vector< std::vector< unsigned int > > &vertex_planes: The id numbers of the planes associated with the vertex
         * :param std::map< unsigned int, std::vector< FloatType > > &planes: The normal - point definitions of the planes
         * :param std::map< unsigned int, std::vector< FloatType > > &couple_tractions: The couple_tractions associated with each plane
         * :param vecOvec &vertex_hostress: The higher order stress at the vertices
         */

        vecOfvec vertex_normals; //The normals associated with the vertices
        vecOfvec vertex_couples; //The couple tractions associated with the vertices

        //Resize the couple stress vector
        vertex_hostress.resize(vertex_planes.size());

        //Iterate through the vertices
        for (unsigned int v=0; v<vertex_planes.size(); v++){

            //Assemble the normals and couple tractions
            vertex_normals.clear(); 
            vertex_normals.resize(vertex_planes[v].size());

            vertex_couples.clear();
            vertex_couples.resize(vertex_planes[v].size());

            //Assemble the tractions and normals
            for (unsigned int f=0; f<vertex_planes[v].size(); f++){

                auto face = std::find(domain.planes.begin(), domain.planes.end(), vertex_planes[v][f]);
                if (face == domain.planes.end()){
                    std::cerr << "Error: vertex plane not found in element_planes\n";
                    assert(1==0);
                }

                auto f_traction = couple_tractions.find(vertex_planes[v][f]);
                if (f_traction == couple_tractions.end()){
                    std::cerr << "Error: couple traction not found in couple_tractions\n";
                    assert(1==0);
                }

                //Save the normal
                vertex_normals[f] = domain.normal( std::distance( domain.planes.begin(), face) );

                //Save the couple traction
                vertex_couples[f] = f_traction->second;
            }

            compute_vertex_couple_stress(vertex_normals, vertex_couples, vertex_hostress[v]);
        }
    }
}
