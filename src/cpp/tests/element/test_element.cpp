/*!
Tests for element.h and element.cpp
*/

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>

#include<element.h>

bool fuzzy_equals(double a, double b, double tolr=1e-6, double tola=1e-6){
    /*!
    Compare two doubles to determine if they are equal.
    */

    double tol = fmin(tolr*fabs(a) + tola, tolr*fabs(b) + tola);
    return fabs(a-b)<tol;
}

bool fuzzy_equals(elib::vec a, elib::vec b, double tolr=1e-6, double tola=1e-6){
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

bool fuzzy_equals(elib::vecOfvec A, elib::vecOfvec B, double tolr=1e-6, double tola=1e-6){
    /*!
    Compare two matrices to determine if they are equal
    */

    if (A.size() != B.size()){
        std::cout << "Error: matrices must have the same size.\n";
        assert(1==0);
    }

    for (unsigned int i=0; i<A.size(); i++){
        if (!fuzzy_equals(A[i], B[i], tolr, tola)){
            return false;
        }
    }
    return true;
}

void print(elib::vec a){
    /*!
    Print the vector to the terminal
    */

    for (unsigned int i=0; i<a.size(); i++){
        std::cout << a[i] << " ";
    }
    std::cout << "\n";
}

void print(elib::vecOfvec A){
    /*!
    Print the matrix to the terminal
    */

    for (unsigned int i=0; i<A.size(); i++){
        print(A[i]);
    }
}

void get_scalar_field_definition(elib::vec &a){
    /*!
    Return the values required to define a linear scalar field

    :param elib::vec &a: The parameters for the scalar field
    */

    a = {0.1, 0, 0};//-0.2, 0.3};
}

void get_vector_field_definition(elib::vecOfvec &A, elib::vec &b){
    /*!
    Return the values required to define a linear vector field

    :param elib::vecOfvec &A: A reference to the linear mapping
    :param elib::vec &b: A reference to the additive portion
    */

    A = {{0.69646919, 0.28613933, 0.22685145},
         {0.55131477, 0.71946897, 0.42310646},
         {0.9807642,  0.68482974, 0.4809319 },
         {0.39211752, 0.34317802, 0.72904971}};

    b = {0.43857224, 0.0596779, 0.39804426, 0.73799541};

    return;
}

double scalar_field(const elib::vec &x){
    /*!
    Compute the value of a scalar field at location x

    :param elib::vec x: The global coordinates
    */

    elib::vec a;
    get_scalar_field_definition(a);

    double value = 0;
    for (unsigned int i=0; i<x.size(); i++){
        value += a[i]*x[i];
    }
    return value;
}

elib::vec vector_field(const elib::vec &x){
    /*!
    Compute the value of a vector field at location x

    :param elib::vec x: The global coordinates
    */

    elib::vecOfvec A;
    elib::vec b;
    get_vector_field_definition(A, b);

    elib::vec out(A.size(), 0);

    for (unsigned int i=0; i<A.size(); i++){
        out[i] = b[i];
        for (unsigned int j=0; j<x.size(); j++){
            out[i] += A[i][j]*x[j];
        }
    }

    return out;
}

void define_hex8_fully_integrated_quadrature(elib::quadrature_rule &qrule){
    /*!
    Define the quadrature rule for a fully integrated hexehedral element.

    :param elib::quadrature_rule &qrule: A reference to the quadrature rule.
    */

    elib::vecOfvec quadrature_points = {{-0.57735027, -0.57735027, -0.57735027},
                                        { 0.57735027, -0.57735027, -0.57735027},
                                        { 0.57735027,  0.57735027, -0.57735027},
                                        {-0.57735027,  0.57735027, -0.57735027},
                                        {-0.57735027, -0.57735027,  0.57735027},
                                        { 0.57735027, -0.57735027,  0.57735027},
                                        { 0.57735027,  0.57735027,  0.57735027},
                                        {-0.57735027,  0.57735027,  0.57735027}};

    elib::vec quadrature_weights = {1, 1, 1, 1, 1, 1, 1, 1};

    qrule.resize(8);
    for (unsigned int i=0; i<8; i++){
        qrule[i] = {quadrature_points[i], quadrature_weights[i]};
    }
    return;
}

int test_Hex8_get_shape_functions(std::ofstream &results){
    /*!
    Test the computation of the shape functions for a Hex8 Element
    */

    // Define the element's nodes
    elib::vecOfvec nodes = {{0, 0, 0},
                            {1, 0, 0},
                            {1, 1, 0},
                            {0, 1, 0},
                            {0, 0, 1},
                            {1, 0, 1},
                            {1, 1, 1},
                            {0, 1, 1}};

    // Define the element's quadrature rule
    elib::quadrature_rule qrule;
    define_hex8_fully_integrated_quadrature(qrule);

    // Construct the element
    elib::Hex8 element(nodes, qrule);

    // Initialize the shape function vector
    elib::vec shape_functions;

    // Make sure the shape functions return the expected value at the center
    element.get_shape_functions({0, 0, 0}, shape_functions);
    elib::vec answer = {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125};
    
    if (!fuzzy_equals(answer, shape_functions)){
        results << "test_Hex8_get_shape_functions (test 1) & False\n";
        return 1;
    }

    // Check the shape function values at the nodes
    for (unsigned int n=0; n<element.local_node_coordinates.size(); n++){
        element.get_shape_functions(element.local_node_coordinates[n], shape_functions);
        for (unsigned int m=0; m<element.local_node_coordinates.size(); m++){
            if (m==n){
                if (!fuzzy_equals(shape_functions[m], 1)){
                    results << "test_Hex8_get_shape_functions (test 2a) & False\n";
                    return 1;
                }
            }
            else{
                if (!fuzzy_equals(shape_functions[m], 0)){
                    results << "test_Hex8_get_shape_functions (test 2b) & False\n";
                    return 1;
                }
            }
        }
    }

    results << "test_Hex8_get_shape_functions & True\n";
    return 0;
}

int test_Hex8_get_local_grad_shape_functions(std::ofstream &results){
    /*!
    Test the computation of the local gradients of the shape functions for a Hex8 Element
    */

    // Define the element's nodes
    elib::vecOfvec nodes = {{0, 0, 0},
                            {1, 0, 0},
                            {1, 1, 0},
                            {0, 1, 0},
                            {0, 0, 1},
                            {1, 0, 1},
                            {1, 1, 1},
                            {0, 1, 1}};

    // Define the element's quadrature rule
    elib::quadrature_rule qrule;
    define_hex8_fully_integrated_quadrature(qrule);

    // Construct the element
    elib::Hex8 element(nodes, qrule);

    //Compute a numeric gradient of the shape functions
    double eps = 1e-6;
    elib::vec sf0, sfpx, sfpy, sfpz;
    element.get_shape_functions({    0.1,     -0.2,     0.3},  sf0);
    element.get_shape_functions({0.1+eps,     -0.2,     0.3}, sfpx);
    element.get_shape_functions({    0.1, -0.2+eps,     0.3}, sfpy);
    element.get_shape_functions({    0.1,     -0.2, 0.3+eps}, sfpz);

    elib::vecOfvec answer;
    answer.resize(8);
    for (unsigned int i=0; i<8; i++){
        answer[i].resize(3);
        answer[i][0] = (sfpx[i] - sf0[i])/eps;
        answer[i][1] = (sfpy[i] - sf0[i])/eps;
        answer[i][2] = (sfpz[i] - sf0[i])/eps;
    }

    elib::vecOfvec local_grad_shape_functions;
    element.get_local_grad_shape_functions({0.1, -0.2, 0.3}, local_grad_shape_functions);

    if (!fuzzy_equals(answer, local_grad_shape_functions)){
        results << "test_Hex8_get_local_grad_shape_functions & False\n";
        return 1;
    }

    results << "test_Hex8_get_local_grad_shape_functions & True\n";
    return 0;
}

int test_interpolate(elib::Element &element, std::ofstream &results){
    /*!
    Test whether interpolation is performed correctly on the element.

    :param elib::Element element: The element to be tested
    :param std::ofstream &results: The output file to write the results to
    */

    //Compute the global coordinates of the nodes via interpolation
    elib::vec value;
    for (unsigned int n=0; n<element.local_node_coordinates.size(); n++){
        element.interpolate(element.nodes, element.local_node_coordinates[n], value);
        if (!fuzzy_equals(value, element.nodes[n])){
            results << element.name.c_str() << "_test_interpolate (test 1) & False\n";
            return 1;
        }
    }

    //Interpolate a constant scalar
    elib::vec scalar(element.nodes.size(), 1);
    double scalar_result;
    element.interpolate(scalar, {-.2, .8, .5}, scalar_result);
    if (!fuzzy_equals(scalar_result, 1)){
        results << element.name.c_str() << "_test_interpolate (test 2) & False\n";
        return 1;
    }

    //Interpolate a variable linear scalar field
    double scalar_answer;
    elib::vec scalar_nodal_values(element.nodes.size());
    for (unsigned int n=0; n<element.nodes.size(); n++){
        scalar_nodal_values[n] = scalar_field(element.nodes[n]);
    }

    elib::vec xi = {-0.2, 0.4, 0.8};
    elib::vec x;
    element.interpolate(element.nodes, xi, x);

    element.interpolate(scalar_nodal_values, xi, scalar_result);
    scalar_answer = scalar_field(x);
    if (!fuzzy_equals(scalar_result, scalar_answer)){
        results << element.name.c_str() << "_test_interpolate (test 3) & False\n";
        return 1;
    }

    results << element.name.c_str() << "_test_interpolate & True\n";
    return 0;
}

int test_element_functionality(elib::Element &element, std::ofstream &results){
    /*!
    Test the provided element's functionality

    :param elib::Element element: The element to be tested
    :param std::ofstream &results: The output file to write the results to

    */

    test_interpolate(element, results);

    return 0;
}


int test_Hex8_functionality(std::ofstream &results){
    /*!
    Test the Hex8 element's functionality
    */

    // Define the element's nodes
    elib::vecOfvec nodes = {{0, 0, 0},
                            {1, 0, 0},
                            {1, 1, 0},
                            {0, 1, 0},
                            {0, 0, 1},
                            {1, 0, 1},
                            {1, 1, 1},
                            {0, 1, 1}};

    // Define the element's quadrature rule
    elib::quadrature_rule qrule;
    define_hex8_fully_integrated_quadrature(qrule);

    // Construct the element
    elib::Hex8 element(nodes, qrule);

    return test_element_functionality(element, results);
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

    //Hex8 tests
    test_Hex8_get_shape_functions(results);
    test_Hex8_get_local_grad_shape_functions(results);
    test_Hex8_functionality(results);

    //Close the results file
    results.close();

    return 0;
}

