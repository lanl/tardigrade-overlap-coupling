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

    a = {0.1, -0.2, 0.3};
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

void get_linear_transformation_definition(elib::vecOfvec &A, elib::vec &b){
    /*!
    Get the definition of a linear transformation applied to the provided vector
    */

    A = {{ 0.26921601, -0.28725274,  0.01841124},
         { 0.19559688,  0.01621845, -1.43394978},
         { 0.33276929,  0.22285938,  0.82795953}};

    b = {1.23409356, 0.50251371, 0.41645453};
    return;
}

void linear_transform(const elib::vec &v, elib::vec &w){
    /*!
    Apply a linear transformation to a vector of a given dimension.

    :param const elib::vec &v: The un-transformed vector
    :param const elib::vec &w: The transformed vector
    */

    elib::vecOfvec A;
    elib::vec b;
    get_linear_transformation_definition(A, b);

    w.resize(v.size());
    for (unsigned int i=0; i<w.size(); i++){
        w[i] = b[i];
        for (unsigned int j=0; j<v.size(); j++){
            w[i] += A[i][j]*v[j];
        }
    }
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

    TODO: Generalize to non-3D in local coordinates elements.

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

    //Interpolate a variable vector field
    elib::vecOfvec vector_nodal_values(element.nodes.size());
    for (unsigned int n=0; n<element.nodes.size(); n++){
        vector_nodal_values[n] = vector_field(element.nodes[n]);
    }

    elib::vec vector_result, vector_answer;
    element.interpolate(vector_nodal_values, xi, vector_result);
    vector_answer = vector_field(x);

    if (!fuzzy_equals(vector_result, vector_answer)){
        results << element.name.c_str() << "_test_interpolate (test 4) & False\n";
        return 1;
    }

    results << element.name.c_str() << "_test_interpolate & True\n";
    return 0;
}

int test_get_local_gradient(elib::Element &element, std::ofstream &results){
    /*!
    Test the computation of the gradient with respect to the local coordinates

    TODO: Generalize to non-3D in local coordinates elements

    :param elib::Element element: The element to be tested
    :param std::ofstream &results: The output file to write the results to
    */

    double eps = 1e-6;
    elib::vec scalar_answer(3, 0), scalar_result(3, 0);
    double sg0, sgpx, sgpy, sgpz;

    //Form the scalar field at the nodes
    elib::vec scalar_nodal_values(element.nodes.size());
    for (unsigned int n=0; n<element.nodes.size(); n++){
        scalar_nodal_values[n] = scalar_field(element.nodes[n]);
    }

    //Interpolate the field
    element.interpolate(scalar_nodal_values, {    -0.2,     0.4,     0.64}, sg0);
    element.interpolate(scalar_nodal_values, {-0.2+eps,     0.4,     0.64}, sgpx);
    element.interpolate(scalar_nodal_values, {    -0.2, 0.4+eps,     0.64}, sgpy);
    element.interpolate(scalar_nodal_values, {    -0.2,     0.4, 0.64+eps}, sgpz);

    //Compute the numeric gradient
    scalar_answer[0] = (sgpx - sg0)/eps;
    scalar_answer[1] = (sgpy - sg0)/eps;
    scalar_answer[2] = (sgpz - sg0)/eps;

    //Compute the element result
    element.get_local_gradient(scalar_nodal_values, {0.1, -0.2, 0.3}, scalar_result);

    if (!fuzzy_equals(scalar_answer, scalar_result)){
        results << element.name.c_str() << "_test_get_local_gradient (test 1) & False\n";
	return 1;
    } 

    elib::vecOfvec vector_answer, vector_result;
    elib::vecOfvec vector_nodal_values(element.nodes.size());
    elib::vec vg0, vgpx, vgpy, vgpz;

    for (unsigned int n=0; n<element.nodes.size(); n++){
        vector_nodal_values[n] = vector_field(element.nodes[n]);
    }

    //Interpolate the field
    elib::vec local_coordinates = {-0.2, 0.4, 0.64};
    elib::vec perturbed_coordinates(local_coordinates.size());
    elib::vecOfvec perturbation_matrix(element.local_node_coordinates[0].size()+1);

    //Interpolate at the nominal position
    element.interpolate(vector_nodal_values, local_coordinates, perturbation_matrix[0]);

    //Perturb the local coordinates and interpolate
    for (unsigned int i=0; i<perturbation_matrix.size()-1; i++){
        for (unsigned int j=0; j<perturbation_matrix.size()-1; j++){
            if (i == j){
                perturbed_coordinates[j] = local_coordinates[j] + eps;
    	    }
	    else{
                perturbed_coordinates[j] = local_coordinates[j];
	    }
	}
	element.interpolate(vector_nodal_values, perturbed_coordinates, perturbation_matrix[i+1]);
    }

    //Approximate the derivative using finite differences
    vector_answer.resize(vector_nodal_values[0].size());
    for (unsigned int i=0; i<vector_answer.size(); i++){
	vector_answer[i].resize(perturbation_matrix.size()-1);
        for (unsigned int j=1; j<perturbation_matrix.size(); j++){
            vector_answer[i][j-1] = (perturbation_matrix[j][i] - perturbation_matrix[0][i])/eps;
	}
    }
    element.get_local_gradient(vector_nodal_values, local_coordinates, vector_result);

    if (!fuzzy_equals(vector_answer, vector_result)){
        results << element.name.c_str() << "_test_get_local_gradient (test 2) & False\n";
	return 1;
    }

    results << element.name.c_str() << "_test_get_local_gradient & True\n";
    return 0;
}

int test_get_global_gradient(elib::Element &element, std::ofstream &results){
    /*!
    Test the computation of the global gradient.

    TODO: Generalize to non-3D elements.

    :param elib::Element element: The element to be tested
    :param std::ofstream &results: The output file to write the results to
    */

    //Compute a set of reference coordinates using a linear transformation
    elib::vecOfvec reference_coordinates(element.nodes.size());
    for (unsigned int n=0; n<element.nodes.size(); n++){
        reference_coordinates[n].resize(element.nodes.size());
        linear_transform(element.nodes[n], reference_coordinates[n]);
    }

    //Evaluate the scalar field
    elib::vec scalar_nodal_current_values(element.nodes.size());
    elib::vec scalar_nodal_reference_values(element.nodes.size());

    for (unsigned int n=0; n<element.nodes.size(); n++){
        scalar_nodal_current_values[n] = scalar_field(element.nodes[n]);
        scalar_nodal_reference_values[n] = scalar_field(reference_coordinates[n]);
    }

    //Compute the two gradients
    elib::vec grad_scalar_current, grad_scalar_reference;
    element.get_global_gradient(scalar_nodal_current_values, {0.1, 0.2, 0.3}, grad_scalar_current);
    element.get_global_gradient(scalar_nodal_reference_values, {0.1, 0.2, 0.3}, reference_coordinates, grad_scalar_reference);

    //Get the expected gradient
    elib::vec scalar_answer;
    get_scalar_field_definition(scalar_answer);

    if (!fuzzy_equals(grad_scalar_current, grad_scalar_reference) || !fuzzy_equals(grad_scalar_current, scalar_answer)){
        results << element.name.c_str() << "_test_get_global_gradient (test 1) & False\n";
        return 1;
    }

    //Evaluate the vector field
    elib::vecOfvec vector_nodal_current_values(element.nodes.size());
    elib::vecOfvec vector_nodal_reference_values(element.nodes.size());

    for (unsigned int n=0; n<element.nodes.size(); n++){
        vector_nodal_current_values[n] = vector_field(element.nodes[n]);
        vector_nodal_reference_values[n] = vector_field(reference_coordinates[n]);
    }

    //Compute the two gradients
    elib::vecOfvec grad_vector_current, grad_vector_reference;
    element.get_global_gradient(vector_nodal_current_values, {0.1, 0.2, 0.3}, grad_vector_current);
    element.get_global_gradient(vector_nodal_reference_values, {0.1, 0.2, 0.3}, reference_coordinates, grad_vector_reference);

    //Get the expected gradient
    elib::vecOfvec vector_answer;
    elib::vec b;
    get_vector_field_definition(vector_answer, b);

    if (!fuzzy_equals(grad_vector_current, grad_vector_reference) || !fuzzy_equals(grad_vector_current, vector_answer)){
        results << element.name.c_str() << "_test_get_global_gradient (test 2) & False\n";
        return 1;
    }


    results << element.name.c_str() << "_test_get_global_gradient & True\n";
    return 0;
}

int test_compute_local_coordinates(elib::Element &element, std::ofstream &results){
    /*!
    Test for the computation of an element's local coordinates given it's global coordinates.

    TODO: Generalize to non-3D elements.

    :param elib::Element element: The element to be tested
    :param std::ofstream &results: The output file to write the results to
    */

    elib::vec x = {0.25, 0.75, .14};
    elib::vec xi, result;

    element.compute_local_coordinates(x, xi);

    element.interpolate(element.nodes, xi, result);

    if (!fuzzy_equals(result, x)){
        results << element.name.c_str() << "_test_compute_local_coordinates & False\n";
        return 1;
    }

    results << element.name.c_str() << "_test_compute_local_coordinates & True\n";
    return 0;
}

int test_element_functionality(elib::Element &element, std::ofstream &results){
    /*!
    Test the provided element's functionality

    :param elib::Element element: The element to be tested
    :param std::ofstream &results: The output file to write the results to

    */

    test_interpolate(element, results);
    test_get_local_gradient(element, results);
    test_get_global_gradient(element, results);
    test_compute_local_coordinates(element, results);

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

int test_invert(std::ofstream &results){
    /*!
    Test the computation of the matrix inverse
    */

    elib::vecOfvec A = {{2, 3,  5},
                        {3, 6,  7},
                        {5, 7, 10}};

    elib::vecOfvec Ainv;
    elib::invert(A, Ainv);

    elib::vecOfvec answer = {{1, 0, 0},
                             {0, 1, 0},
                             {0, 0, 1}};
    elib::vecOfvec result;

    result.resize(A.size());
    for (unsigned int i=0; i<A.size(); i++){
        result[i].resize(A[i].size());
        for (unsigned int j=0; j<A[i].size(); j++){
            result[i][j] = 0;
            for (unsigned int k=0; k<A[i].size(); k++){
                result[i][j] += A[i][k]*Ainv[k][j];
            }
        }
    }

    if (!fuzzy_equals(result, answer)){
        results << "test_invert & False\n";
        return 1;
    }

    results << "test_invert & True\n";
    return 0;
}

int test_solve(std::ofstream &results){
    /*!
    Test the computation of the matrix solve
    */

    elib::vecOfvec A = {{2, 3,  5},
                        {3, 6,  7},
                        {5, 7, 10}};

    elib::vec answer = {1, 2, 3};

    elib::vec b(A.size(), 0);

    for (unsigned int i=0; i<A.size(); i++){
        for (unsigned int j=0; j<A[i].size(); j++){
            b[i] += A[i][j]*answer[j];
        }
    }

    elib::vec result;
    elib::solve(A, b, result);

    if (!fuzzy_equals(answer, result)){
        results << "test_solve & False\n";
        return 1;
    }

    results << "test_solve & True\n";
    return 0;
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

    //Eigen tool tests
    test_invert(results);
    test_solve(results);

    //Close the results file
    results.close();

    return 0;
}

