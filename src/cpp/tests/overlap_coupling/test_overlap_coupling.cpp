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

    vecOfvec in;
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
    vecOfvec result;
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
    const planeMap *element_planes = oc.get_element_planes();
    const vecOfvec *element_bounds = oc.get_element_bounds();

    if (element_planes->size()!=6){
        results << "test_compute_element_bounds (test 1)& False\n";
        return;
    }

    //!Assumes the underlying element is a hex
    planeMap::const_iterator it;
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
    planeMap dns_planes;
    std::vector< vertex_t > vertices;
    oc.map_vectors_to_quickhull(data.coordinates, vertices);

    vecOfvec bounds;
    bounds.resize(3);
    oc.compute_node_bounds(data.coordinates, dns_planes, bounds[0], bounds[1], bounds[2], 1e-9, 1e-9);

    if (dns_planes.size() != 6){
        results << "test_compute_node_bounds (test 1) & False\n";
        return;
    }

    vecOfvec answer(3);
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

    vecOfvec normals;
    vecOfvec points;

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
        auto vertexBuffer = mesh.getVertexBuffer();

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
    const planeMap *dns_planes = oc.get_dns_planes();
    const vecOfvec *dns_bounds = oc.get_dns_bounds();

    //!Compare bounds to expected values
    vecOfvec answer(3);
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
    planeMap::const_iterator it;
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

    //Tests for the interface to 3D-quickhull
    test_map_vector_to_quickhull(results);
    test_map_vectors_to_quickhull(results);
    test_map_quickhull_to_vector(results);
    test_map_quickhull_to_vectors(results);

    //Test for the computations of the bounds
    test_extract_mesh_info(results);
    test_compute_element_bounds(results);
    test_compute_node_bounds(results);
    test_compute_dns_bounds(results);

    //Test misc. functions
    test_dot(results);
    test_cross(results);
    test_fuzzy_equals(results);
    test_compare_vector_directions(results);
    test_normal_from_vertices(results);

    //Close the results file
    results.close();
}
