//!The test file for overlap_coupling.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>

#include "occonfiguration.h"
#include "overlap_coupling.h"

bool fuzzy_equals(double a, double b, double tolr=1e-6, double tola=1e-6){
    /*!
    Compare two doubles to determine if they are equal.
    */

    double tol = fmin(tolr*fabs(a) + tola, tolr*fabs(b) + tola);
    return fabs(a-b)<tol;
}

bool fuzzy_equals(std::vector< double > a, std::vector< double > b, double tolr=1e-6, double tola=1e-6){
    /*!
    Compare two vectors to determine if they are equal
    */

    if (a.size() != b.size()){
        std::cout << "Error: vectors must have the same size.\n";
        assert(1==0);
    }

    for (unsigned int i=0; i<a.size(); i++){
        if (!fuzzy_equals(a[i], b[i], tolr, tola)){
            return false;
        }
    }
    return true;
}

//!From voro++ documentation
// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}
//!End from voro++ documentation

template<typename T>
void print_vector(std::vector< T > vector){
    /*!
    Print the vector to the terminal
    */
    for (unsigned int i=0; i<vector.size(); i++){
        std::cout << vector[i] << " ";
    }
    std::cout << "\n";
    return;
}

template<typename T>
void print_matrix(std::vector< std::vector< T > > matrix){
    /*!
    Print the matrix to the terminal.
    */
    for (unsigned int i=0; i<matrix.size(); i++){
        print_vector(matrix[i]);
    }
}

void test_map_vector_to_quickhull(std::ofstream &results){
    /*!
    Test mapping a std::vector to a 3D-quickhull vertex.
    */

    std::vector< double > a(3,0);
    a[1] = 1;
    a[2] = 2;

    overlap::OverlapCoupling oc;
    vertex_t result = oc.map_vector_to_quickhull(a);

    if (!((result.x == a[0]) && (result.y == a[1]) && (result.z == a[2]))){
        results << "test_map_vector_to_quickhull & False\n";
    }
    else{
        results << "test_map_vector_to_quickhull & True\n";
    }
}

void test_map_quickhull_to_vector(std::ofstream &results){
    /*!
    Test mapping a 3D-quickhull vertex to a std::vector
    */

    vertex_t v;
    v.x = 1.2;
    v.y = 3.7;
    v.z = -1.2;

    overlap::OverlapCoupling oc;
    std::vector< double > result = oc.map_quickhull_to_vector(v);

    if (!(fuzzy_equals(result[0],v.x) && (fuzzy_equals(result[1],v.y)) && (fuzzy_equals(result[2], v.z)))){
        results << "test_map_quickhull_to_vector & False\n";
        return;
    }

    results << "test_map_quickhull_to_vector & True\n";
    return;
}

void test_map_vectors_to_quickhull(std::ofstream &results){
    /*!
    Test mapping a collection of std::vectors to a 3D-quickhull vertex.
    */

    std::vector< double > a(3,0), b(3,-1);
    a[1] = 1;
    a[2] = 2;
    b[1] = .32;
    b[2] = 7.8;

    overlap::vecOfvec in;
    in.push_back(a);
    in.push_back(b);

    overlap::OverlapCoupling oc;
    std::vector< vertex_t > result;
    oc.map_vectors_to_quickhull(in, result);

    for (unsigned int i=0; i<in.size(); i++){
        if (!fuzzy_equals(in[i][0], result[i].x) || !fuzzy_equals(in[i][1], result[i].y) || !fuzzy_equals(in[i][2], result[i].z)){
//            print_vector(in[i]);
//            std::cout << result[i].x << " " << result[i].y << " " << result[i].z << "\n";
            results << "test_map_vectors_to_quickhull & False\n";
            return;
        }
    }

    results << "test_map_vectors_to_quickhull & True\n";
    return;
}

void test_map_quickhull_to_vectors(std::ofstream &results){
    /*!
    Test mapping a collection of 3D-quickhull vertices to a vector of std::vectors
    */

    vertex_t v1, v2;
    v1.x = 1.;
    v1.y = 2.;
    v1.z = 3.;

    v2.x = .27;
    v2.y = 1.23;
    v2.z = -2.1;

    std::vector< vertex_t > in;
    in.push_back(v1);
    in.push_back(v2);

    overlap::OverlapCoupling oc;
    overlap::vecOfvec result;
    oc.map_quickhull_to_vectors(in, result);

    for (unsigned int i=0; i<in.size(); i++){
        if (!(fuzzy_equals(result[i][0], in[i].x) && fuzzy_equals(result[i][1], in[i].y) && fuzzy_equals(result[i][2], in[i].z))){
            results << "test_map_quickhull_to_vectors & False\n";
            return;
        }
    }

    results << "test_map_quickhull_to_vectors & True\n";
    return;
}

void test_dot(std::ofstream &results){
    /*!
    Test the computation of the dot product of two vectors.
    */

    std::vector< double > a(3,0), b(3,-1);
    a[1] = 1;
    a[2] = 2;
    b[1] = 0.32;
    b[2] = 7.8;

    double result = overlap::dot(a, b);
    
    double answer = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];

    if (!fuzzy_equals(result, answer)){
        results << "test_dot & False\n";
        return;
    }
    results << "test_dot & True\n";
    return;
    
}

void test_cross(std::ofstream &results){
    /*!
    Test the computation of the cross product of two vectors.
    */

    std::vector< double > a(3,0), b(3,-1);
    a[1] = 1;
    a[2] = 2;
    b[1] = 0.32;
    b[2] = 7.8;

    std::vector< double > result = overlap::cross(a, b);

    if ((!fuzzy_equals(overlap::dot(result, a), 0)) || (!fuzzy_equals(overlap::dot(result, b), 0))){
        results << "test_cross & False\n";
        return;
    }
    results << "test_cross & True\n";
    return;
}

void test_fuzzy_equals(std::ofstream &results){
    /*!
    Test the comparison of two values using the fuzzy (tolerant) comparison.
    */

    double a = 1;
    double b = 1;

    if (!overlap::fuzzy_equals(a, b)){
        results << "test_fuzzy_equals (test 1) & False\n";
        return;
    }

    a += 1e-3;
    if (overlap::fuzzy_equals(a, b)){
        results << "test_fuzzy_equals (test 2) & False\n";
        return;
    }

    a = -1;
    b = -1;

    if (!overlap::fuzzy_equals(a, b)){
        results << "test_fuzzy_equals (test 3) & False\n";
        return;
    }

    results << "test_fuzzy_equals & True\n";
    return;
}

void test_compare_vector_directions(std::ofstream &results){
    /*!
    Test the comparison of two vector directions for equality.
    */
//    std::cout << "Testing the comparision of vector directions...\n";

    std::vector< double > a(3, 1), b(3, 3);
    
    if (!overlap::compare_vector_directions(a, b)){
        results << "test_compare_vector_directions (test 1) & False\n";
        return;
    }

    a[0] += 1;
    if (overlap::compare_vector_directions(a, b)){
        results << "test_compare_vector_directions (test 2) & False\n";
        return;
    }

    a[0] = 1;
    a[1] = 2;
    a[2] = 3;
    b[0] = -1;
    b[1] = -2;
    b[2] = -3;

    if (overlap::compare_vector_directions(a, b)){
        results << "test_compare_vector_directions (test 3) & False\n";
        assert(1==0);
        return;
    }

    results << "test_compare_vector_directions & True\n";
    return;
}

void test_compute_element_bounds(std::ofstream &results){
    /*!
    Test the computation of the element bounds (also tests compute_unique_planes)
    */
//    std::cout << "Testing the computation of the element bounds...\n";
    overlap::ParsedData data = overlap::read_data_from_file("overlap.txt");
    overlap::OverlapCoupling oc(data.local_nodes, data.local_gpts);
    const overlap::planeMap *element_planes = oc.get_element_planes();
    const overlap::vecOfvec *element_bounds = oc.get_element_bounds();

    if (element_planes->size()!=6){
        results << "test_compute_element_bounds (test 1)& False\n";
        return;
    }

    //!Assumes the underlying element is a hex
    overlap::planeMap::const_iterator it;
    for (it=element_planes->begin(); it!=element_planes->end(); it++){
        for (unsigned int i=0; i<it->first.size(); i++){
//            print_vector(it->first);
//            print_vector(it->second);
            if (fuzzy_equals(fabs(it->first[i]), 1)){
                if (!fuzzy_equals(it->first[i], it->second[i])){
                    results << "test_compute_element_bounds (test 2) & False\n";
                    return;
                }
            }
        }
    }

    //!Check the bounds (assuming it is a hex)
    std::vector< double > answer(2, 1);
    answer[0] = -1;
    for (unsigned int i=0; i<3; i++){
        if (!fuzzy_equals(answer, (*element_bounds)[i])){
            results << "test_compute_element_bounds(test 3) & False\n";
            return;
        }
    }

    results << "test_compute_element_bounds & True\n";
    return;
}

void test_compute_node_bounds(std::ofstream &results){
    /*!
    Test the computation of the node bounds.
    */
//    std::cout << "Testing the computation of the node bounds...\n";

    overlap::ParsedData data = overlap::read_data_from_file("overlap.txt");
    overlap::OverlapCoupling oc = overlap::OverlapCoupling(data.local_nodes, data.local_gpts);
    overlap::planeMap dns_planes;
    std::vector< vertex_t > vertices;
    oc.map_vectors_to_quickhull(data.coordinates, vertices);

    overlap::vecOfvec bounds;
    bounds.resize(3);
    oc.compute_node_bounds(data.coordinates, dns_planes, bounds[0], bounds[1], bounds[2], 1e-9, 1e-9);

    if (dns_planes.size() != 6){
        results << "test_compute_node_bounds (test 1) & False\n";
        return;
    }

    overlap::vecOfvec answer(3);
    for (unsigned int i=0; i<answer.size(); i++){
        answer[i] = std::vector< double >(2,1);
    }
    answer[0][0] = 0;
    answer[1][0] = answer[2][0] = -1;

    for (unsigned int i=0; i<3; i++){
        if (!(fuzzy_equals(answer[i], bounds[i]))){
            results << "test_compute_node_bounds (test 2) & False\n";
            return;
        }
    }

    results << "test_compute_node_bounds & True\n";
    return;
}

void test_extract_mesh_info(std::ofstream &results){
    /*!
    Test the extraction of the mesh information
    */
//    std::cout << "Testing the extraction of the mesh information...\n";

    overlap::ParsedData data = overlap::read_data_from_file("overlap.txt");
    overlap::OverlapCoupling oc;
    std::vector< vertex_t > vertices;
    oc.map_vectors_to_quickhull(data.local_nodes, vertices);
    #if CONVEXLIB == QUICKHULL
        mesh_t mesh = qh_quickhull3d(&vertices[0], vertices.size());
    #elif CONVEXLIB == CONVHULL_3D
        mesh_t mesh;
        int *faceIndices = NULL;
        int nFaces;
        convhull_3d_build(&vertices[0], vertices.size(), &faceIndices, &nFaces);
        mesh.first.assign(faceIndices, faceIndices + nFaces);
        mesh.second = vertices;
    #elif CONVEXLIB == AKUUKKA
        quickhull::QuickHull<FloatType> qh;
        mesh_t mesh = qh.getConvexHull(vertices, false, false);
    #endif

    overlap::vecOfvec normals;
    overlap::vecOfvec points;

    oc.extract_mesh_info(mesh, normals, points);

    #if CONVEXLIB == QUICKHULL

        if (mesh.nnormals != normals.size()){
            results << "test_extract_mesh_info (test 1) & False\n";
            return;
        }
        if (mesh.nnormals != points.size()){
            results << "test_extract_mesh_info (test 2) & False\n";
            return;
        }

    #elif CONVEXLIB == CONVHULL_3D
        if (mesh.first.size()/3 != normals.size()){
            results << "test_extract_mesh_info (test 1) & False\n";
            return;
        }
        if (mesh.first.size()/3 != points.size()){
            results << "test_extract_mesh_info (test 2) & False\n";
            return;
        }
    #elif CONVEXLIB == AKUUKKA
        auto indexBuffer = mesh.getIndexBuffer();

        if (indexBuffer.size()/3 != normals.size()){
            results << "test_extract_mesh_info (test 1) & False\n";
            return;
        }
        if (indexBuffer.size()/3 != points.size()){
            results << "test_extract_mesh_info (test 2) & False\n";
            return;
        }
    #endif

    unsigned int index = 0;
    vertex_t temp_vertex;
    for (unsigned int i=0; i<normals.size(); i++){
        #if CONVEXLIB == QUICKHULL
            if (!(fuzzy_equals(normals[i][0], mesh.normals[i].x) && fuzzy_equals(normals[i][1], mesh.normals[i].y) && fuzzy_equals(normals[i][2], mesh.normals[i].z))){
                results << "test_extract_mesh_info (test 3) & False\n";
                return;
            }
            temp_vertex = mesh.vertices[mesh.indices[index]];

            if (!(fuzzy_equals(points[i][0], temp_vertex.x) && fuzzy_equals(points[i][1], temp_vertex.y) && fuzzy_equals(points[i][2], temp_vertex.z))){
                results << "test_extract_mesh_info (test 4) & False\n";
                return;
            }
//        #elif CONVEXLIB == CONVHULL_3D
//            if (!(fuzzy_equals(normals[i][0], mesh.normals[i].x) && fuzzy_equals(normals[i][1], mesh.normals[i].y) && fuzzy_equals(normals[i][2], mesh.normals[i].z))){
//                results << "test_extract_mesh_info (test 3) & False\n";
//                return;
//            }
//            temp_vertex = mesh.vertices[mesh.indices[index]];
//
//            if (!(fuzzy_equals(points[i][0], temp_vertex.x) && fuzzy_equals(points[i][1], temp_vertex.y) && fuzzy_equals(points[i][2], temp_vertex.z))){
//                results << "test_extract_mesh_info (test 4) & False\n";
//                return;
//            }

        #endif


        index += 3;

    }

    results << "test_extract_mesh_info & True\n";
    return;
}

void test_normal_from_vertices(std::ofstream &results){
    /*!
    Test the computation of a normal from a set of three vertices that 
    define a plane.
    */
//    std::cout << "Testing normal from vertices...\n";
    vertex_t v1, v2, v3;
    v1.x = 1;
    v1.y = 0;
    v1.z = 0;

    v2.x = 0;
    v2.y = 1;
    v2.z = 0;

    v3.x = 0;
    v3.y = 0;
    v3.z = 1;

    std::vector< double > result = overlap::normal_from_vertices(v1, v2, v3);

    std::vector< double > answer(3, 1./sqrt(3));

    if (!fuzzy_equals(answer, result)){
        results << "test_normal_from_vertices & False\n";
        return;
    }

    results << "test_normal_from_vertices & True\n";
    return;
}

void test_compute_dns_bounds(std::ofstream &results){
    /*!
    Test the computation of the bounding planes for the DNS point positions.
    */
//    std::cout << "Testing the computation of the dns bounds...\n";

    overlap::ParsedData data = overlap::read_data_from_file("overlap.txt");
    overlap::OverlapCoupling oc = overlap::OverlapCoupling(data.local_nodes, data.local_gpts);
    oc.compute_dns_bounds(data.coordinates);
    const overlap::planeMap *dns_planes = oc.get_dns_planes();
    const overlap::vecOfvec *dns_bounds = oc.get_dns_bounds();

    //!Compare bounds to expected values
    overlap::vecOfvec answer(3);
    for (unsigned int i=0; i<answer.size(); i++){
        answer[i] = std::vector< double >(2,1);
    }

    answer[0][0] = 0;
    answer[1][0] = answer[2][0] = -1;

    for (unsigned int i=0; i<answer.size(); i++){
        if (!fuzzy_equals(answer[i], dns_bounds->operator[](i))){
            results << "test_compute_dns_bounds (test 1) & False\n";
            return;
        }
    }
    //!Assumes the underlying dns has a hexahedral domain
    overlap::planeMap::const_iterator it;
    for (it=dns_planes->begin(); it!=dns_planes->end(); it++){
        for (unsigned int i=0; i<it->first.size(); i++){
            if (fuzzy_equals(fabs(it->first[i]), 1)){
                if (fuzzy_equals(it->first[i], 1) && !fuzzy_equals(it->second[i], (*dns_bounds)[i][1])){
                    results << "test_compute_dns_bounds (test 2) & False\n";
                    return;
                }
                else if (fuzzy_equals(it->first[i], -1) && !fuzzy_equals(it->second[i], (*dns_bounds)[i][0])){
                    results << "test_compute_dns_bounds (test 2) & False\n";
                    return;
                }
            }
        }
    }

    results << "test_compute_dns_bounds & True\n";
    return;
}

void test_construct_container(std::ofstream &results){
    /*!
    Tests the construction of a voro++ container class object. Also tests add_planes_to_container,
    evaluate_container_information, find_face_centroid, and map_planes_to_voro.
    */

    //Set the number of particles
    const int particles=64;

    const double x_min=-2, x_max=2;
    const double y_min=-2, y_max=2;
    const double z_min=-2, z_max=2;

    //Set the bounds
    overlap::vecOfvec bounds;
    bounds.resize(3);
    for (unsigned int i=0; i<bounds.size(); i++){
        bounds[i].resize(2);
        bounds[i][0] = -2;
        bounds[i][1] =  2;
    }

    //Build the bounding planes (tetrahedron)
    overlap::planeMap planes;
    planes.insert(std::pair< std::vector< double >, std::vector< double > >({ 1, 1, 1}, { 1, 0, 0}));
    planes.insert(std::pair< std::vector< double >, std::vector< double > >({-1,-1, 1}, {-1, 0, 0}));
    planes.insert(std::pair< std::vector< double >, std::vector< double > >({ 1,-1,-1}, { 0, 0,-1}));
    planes.insert(std::pair< std::vector< double >, std::vector< double > >({-1, 1,-1}, { 0, 1, 0}));

    //Construct voro++ planes from the definitions
    std::vector< voro::wall_plane > vplanes;
    overlap::map_planes_to_voro(planes, vplanes);
//    vplanes.reserve(planes.size());
//    double distance;
//    overlap::planeMap::iterator it;
//    int j=1;
//    for (it=planes.begin(); it!=planes.end(); it++){
//        distance = sqrt(overlap::dot(it->first, it->second));
//        vplanes.push_back(voro::wall_plane(it->first[0], it->first[1], it->first[2], distance, -j));
//        j++;
//    }

    //Define the point coordinates
    std::vector< unsigned int > point_numbers;
    overlap::vecOfvec point_coords;

    point_numbers.reserve(particles);
    point_coords.reserve(particles);

    double x, y, z;
    for (unsigned int i=0; i<particles; i++){
        point_numbers.push_back(i);
        x = x_min + rnd()*(x_max - x_min);
        y = y_min + rnd()*(y_max - y_min);
        z = z_min + rnd()*(z_max - z_min);
        point_coords.push_back({x, y, z});
    }

    //Construct the container
    voro::container* container = overlap::construct_container(point_numbers, point_coords, bounds, vplanes);
    overlap::integrateMap points;
    overlap::integrateMap::iterator itiM;
    overlap::evaluate_container_information( container, points);
    double result_d = 0;
    for (itiM=points.begin(); itiM!=points.end(); itiM++){
        result_d += itiM->second.volume;
    }
    double answer_d = 8./3;

    //Check that the volume is what was expected    
    if (!fuzzy_equals(result_d, answer_d)){
        results << "test_construct_container (test 1) & False\n";
        delete(container);
        return;
    }

    //Check that the surface areas are what was expected
    std::vector< double > sub_surface_areas(4, 0);
    for (itiM=points.begin(); itiM!=points.end(); itiM++){
        for (unsigned int j=0; j<itiM->second.areas.size(); j++){
            sub_surface_areas[itiM->second.planes[j]] += itiM->second.areas[j];
        }
    }

    answer_d = sqrt(12);
    std::vector< double >::iterator vdit = sub_surface_areas.begin();
    while (vdit != sub_surface_areas.end()){
        if (!fuzzy_equals(*vdit, answer_d)){
            results << "test_construct_container (test 2) & False\n";
            delete(container);
            return;
        }
        vdit++;
    }

    //Check that the normals for each plane are consistent with expectations and that the centroids are on the plane
    overlap::planeMap::iterator it;
    std::vector< double > normal(3);
    for (itiM=points.begin(); itiM!=points.end(); itiM++){
        for (unsigned int i=0; i<itiM->second.planes.size(); i++){
            it = planes.begin();
            for (int j=0; j<itiM->second.planes[i]; j++){it++;}
            normal = it->first;
            for (unsigned int j=0; j<normal.size(); j++){normal[j] /= sqrt(overlap::dot(it->first, it->first));}
            if (!fuzzy_equals(normal, itiM->second.normals[i])){
                results << "test_construct_container (test 3) & False\n";
                delete(container);
                return;
            }
            if (!fuzzy_equals(overlap::dot(normal, it->second), overlap::dot(normal, itiM->second.face_centroids[i]))){
                results << "test_construct_container (test 4) & False\n";
                delete(container);
                return;
            }
        }
    }

    results << "test_construct_container & True\n";
    delete(container);
    return;
}

void test_construct_gauss_domains(std::ofstream &results){
    /*!
    Test the construction of the gauss domains.
    */

    overlap::ParsedData data = overlap::read_data_from_file("overlap.txt");
    overlap::OverlapCoupling oc = overlap::OverlapCoupling(data.local_nodes, data.local_gpts);
    const std::vector< overlap::MicroPoint > *gauss_domains = oc.get_gauss_domains();
    std::vector< double > pos(3,0);

    //Macro-element assumed to be a fully integrated linear hex
    for (unsigned int i=0; i<gauss_domains->size(); i++){
    
        //See if the volume is 1 as expected    
        if (!fuzzy_equals((*gauss_domains)[i].volume, 1.)){
            results << "test_construct_gauss_domains (test 1) & False\n";
            return;
        }

        //Make sure that the centroids are located where they are expected to be
        for (unsigned int j=0; j<(*gauss_domains)[i].coordinates.size(); j++){
            if (!fuzzy_equals(fabs((*gauss_domains)[i].coordinates[j]), 0.5)){
                results << "test_construct_gauss_domains (test 2) & False\n";
                return;
            }

            if (!fuzzy_equals((*gauss_domains)[i].coordinates[j]/fabs((*gauss_domains)[i].coordinates[j]),
                              data.local_gpts[i][j]/fabs(data.local_gpts[i][j]))){
                results << "test_construct_gauss_domains (test 2) & False\n";
                return;
            }
        }

        //Make sure the surface areas are all 1
        for (unsigned int j=0; j<(*gauss_domains)[i].areas.size(); j++){
            if (!fuzzy_equals((*gauss_domains)[i].areas[j], 1.)){
                results << "test_construct_gauss_domains (test 3) & False\n";
                return;
            }
        }

        //Make sure the centroid is contained within the surfaces
        for (unsigned int j=0; j<(*gauss_domains)[i].normals.size(); j++){
            pos[0] = (*gauss_domains)[i].coordinates[0] - (*gauss_domains)[i].face_centroids[j][0];
            pos[1] = (*gauss_domains)[i].coordinates[1] - (*gauss_domains)[i].face_centroids[j][1];
            pos[2] = (*gauss_domains)[i].coordinates[2] - (*gauss_domains)[i].face_centroids[j][2];

            if (overlap::dot((*gauss_domains)[i].normals[j], pos)>0){
                results << "test_constructed_gauss_domains (test 4) & False\n";
                return;
            }
        }

    
    }

    results << "test_construct_gauss_domains & True\n";
    return;
}

void test_compute_weights(std::ofstream &results){
    /*!
    Test to make sure that the computation of weights and other required quantities is performed correctly.

    Also tests OverlapCoupling::map_domain_to_voro.
    */

    overlap::ParsedData data = overlap::read_data_from_file("overlap.txt");
    overlap::OverlapCoupling oc = overlap::OverlapCoupling(data.local_nodes, data.local_gpts);
    std::vector< overlap::integrateMap > points;
    oc.compute_weights(data.node_numbers, data.coordinates, points);
}

int main(){
    /*!
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or False 
    if the test passes or fails respectively.
    */

    //Open the results file
    std::ofstream results;
    results.open("results.tex");

    overlap::ParsedData data = overlap::read_data_from_file("overlap.txt");

    //Tests for the interface to the hull building routines
    test_map_vector_to_quickhull(results);
    test_map_vectors_to_quickhull(results);
    test_map_quickhull_to_vector(results);
    test_map_quickhull_to_vectors(results);

    //Test for the computations of the bounds
    test_extract_mesh_info(results);
    test_compute_element_bounds(results);
    test_compute_node_bounds(results);
    test_compute_dns_bounds(results);
    test_construct_gauss_domains(results);

    //Tests for the interface to Voro++
    test_construct_container(results);

    //Test misc. functions
    test_dot(results);
    test_cross(results);
    test_fuzzy_equals(results);
    test_compare_vector_directions(results);
    test_normal_from_vertices(results);

    //Close the results file
    results.close();
}
