/*!
Tests for element.h and element.cpp
*/

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>
#include<memory>

#include<element.h>

#define BOOST_TEST_MODULE test_element
#include <boost/test/included/unit_test.hpp>

typedef elib::uitype uitype;
typedef elib::uivec uivec;
typedef elib::vecOfuivec uimat;
typedef elib::errorNode errorNode;
typedef elib::errorOut errorOut;

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

    for (uitype i=0; i<a.size(); i++){
        BOOST_CHECK( fuzzy_equals( a[i], b[i], tolr, tola ) );
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

    for (uitype i=0; i<A.size(); i++){
        BOOST_CHECK( fuzzy_equals( A[i], B[i], tolr, tola ) );
    }
    return true;
}

void print(elib::vec a){
    /*!
    Print the vector to the terminal
    */

    for (uitype i=0; i<a.size(); i++){
        std::cout << a[i] << " ";
    }
    std::cout << "\n";
}

void print(elib::vecOfvec A){
    /*!
    Print the matrix to the terminal
    */

    for (uitype i=0; i<A.size(); i++){
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
    for (uitype i=0; i<x.size(); i++){
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

    for (uitype i=0; i<A.size(); i++){
        out[i] = b[i];
        for (uitype j=0; j<x.size(); j++){
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
    for (uitype i=0; i<w.size(); i++){
        w[i] = b[i];
        for (uitype j=0; j<v.size(); j++){
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
    for (uitype i=0; i<8; i++){
        qrule[i] = {quadrature_points[i], quadrature_weights[i]};
    }
    return;
}

BOOST_AUTO_TEST_CASE( testHex8_get_shape_functions ){
    /*!
    Test the computation of the shape functions for a Hex8 Element
    */

    // Define the element's nodes
    std::vector< uitype > node_ids  = {1, 2, 3, 4, 5, 6, 7, 8};

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
    elib::Hex8 element(node_ids, nodes, qrule);

    // Initialize the shape function vector
    elib::vec shape_functions;

    // Make sure the shape functions return the expected value at the center
    element.get_shape_functions({0, 0, 0}, shape_functions);
    elib::vec answer = {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125};
    
    BOOST_CHECK( fuzzy_equals( answer, shape_functions ) );

    // Check the shape function values at the nodes
    for (uitype n=0; n<element.local_node_coordinates.size(); n++){
        element.get_shape_functions(element.local_node_coordinates[n], shape_functions);
        for (uitype m=0; m<element.local_node_coordinates.size(); m++){
            if (m==n){
                BOOST_CHECK( fuzzy_equals( shape_functions[m], 1 ) );
            }
            else{
                BOOST_CHECK( fuzzy_equals( shape_functions[m], 0 ) );
            }
        }
    }

}

BOOST_AUTO_TEST_CASE( testHex8_get_local_grad_shape_functions ){
    /*!
    Test the computation of the local gradients of the shape functions for a Hex8 Element
    */

    // Define the element's nodes
    std::vector< uitype > node_ids  = {1, 2, 3, 4, 5, 6, 7, 8};

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
    elib::Hex8 element(node_ids, nodes, qrule);

    //Compute a numeric gradient of the shape functions
    double eps = 1e-6;
    elib::vec sf0, sfpx, sfpy, sfpz;
    element.get_shape_functions({    0.1,     -0.2,     0.3},  sf0);
    element.get_shape_functions({0.1+eps,     -0.2,     0.3}, sfpx);
    element.get_shape_functions({    0.1, -0.2+eps,     0.3}, sfpy);
    element.get_shape_functions({    0.1,     -0.2, 0.3+eps}, sfpz);

    elib::vecOfvec answer;
    answer.resize(8);
    for (uitype i=0; i<8; i++){
        answer[i].resize(3);
        answer[i][0] = (sfpx[i] - sf0[i])/eps;
        answer[i][1] = (sfpy[i] - sf0[i])/eps;
        answer[i][2] = (sfpz[i] - sf0[i])/eps;
    }

    elib::vecOfvec local_grad_shape_functions;
    element.get_local_grad_shape_functions({0.1, -0.2, 0.3}, local_grad_shape_functions);

    BOOST_CHECK( fuzzy_equals( answer, local_grad_shape_functions ) );

    //Test for a distorted element
    nodes = {{3.13443, -0.61357,  1.90472},
             {4.24588,  1.41151,  3.82988},
             {3.97724,  1.34621,  4.43285},
             {2.86579, -0.678866, 2.50769},
             {3.95241, -0.996794, 1.71353},
             {5.06385,  1.02829,  3.63869},
             {4.79521,  0.96299,  4.24166},
             {3.68377, -1.06209,  2.3165}};

    element = elib::Hex8(node_ids, nodes, qrule);
    element.get_shape_functions({    0.1,     -0.2,     0.3},  sf0);
    element.get_shape_functions({0.1+eps,     -0.2,     0.3}, sfpx);
    element.get_shape_functions({    0.1, -0.2+eps,     0.3}, sfpy);
    element.get_shape_functions({    0.1,     -0.2, 0.3+eps}, sfpz);

    for (uitype i=0; i<8; i++){
        answer[i][0] = (sfpx[i] - sf0[i])/eps;
        answer[i][1] = (sfpy[i] - sf0[i])/eps;
        answer[i][2] = (sfpz[i] - sf0[i])/eps;
    }

    element.get_local_grad_shape_functions({0.1, -0.2, 0.3}, local_grad_shape_functions);

    BOOST_CHECK( fuzzy_equals( answer, local_grad_shape_functions ) );

}

BOOST_AUTO_TEST_CASE( testHex8_local_point_inside ){
    /*!
    Test the determination of a point in local coordinates is inside the Hex8 element.
    */

    // Define the element's nodes
    std::vector< uitype > node_ids  = {1, 2, 3, 4, 5, 6, 7, 8};

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
    elib::Hex8 element(node_ids, nodes, qrule);

    // Set the xi value
    elib::vec xi(3, 0);
    BOOST_CHECK( element.local_point_inside( xi ) );

    for (uitype i=0; i<xi.size(); i++){
        xi[i] = 2;
        BOOST_CHECK( !element.local_point_inside( xi ) );
        xi[i] = -2;
        BOOST_CHECK( !element.local_point_inside( xi ) );
        xi[i] = 0;
    }

}

int test_interpolate( elib::Element &element ){
    /*!
    Test whether interpolation is performed correctly on the element.

    TODO: Generalize to non-3D in local coordinates elements.

    :param elib::Element element: The element to be tested
    */

    //Compute the global coordinates of the nodes via interpolation
    elib::vec value;
    for (uitype n=0; n<element.local_node_coordinates.size(); n++){
        element.interpolate(element.nodes, element.local_node_coordinates[n], value);
        BOOST_CHECK( fuzzy_equals( value, element.nodes[n] ) );
    }

    //Interpolate a constant scalar
    elib::vec scalar(element.nodes.size(), 1);
    double scalar_result;
    element.interpolate(scalar, {-.2, .8, .5}, scalar_result);
    BOOST_CHECK( fuzzy_equals( scalar_result, 1 ) );

    //Interpolate a variable linear scalar field
    double scalar_answer;
    elib::vec scalar_nodal_values(element.nodes.size());
    for (uitype n=0; n<element.nodes.size(); n++){
        scalar_nodal_values[n] = scalar_field(element.nodes[n]);
    }

    elib::vec xi = {-0.2, 0.4, 0.8};
    elib::vec x;
    element.interpolate(element.nodes, xi, x);

    element.interpolate(scalar_nodal_values, xi, scalar_result);
    scalar_answer = scalar_field(x);
    BOOST_CHECK( fuzzy_equals( scalar_result, scalar_answer ) );

    //Interpolate a variable vector field
    elib::vecOfvec vector_nodal_values(element.nodes.size());
    for (uitype n=0; n<element.nodes.size(); n++){
        vector_nodal_values[n] = vector_field(element.nodes[n]);
    }

    elib::vec vector_result, vector_answer;
    element.interpolate(vector_nodal_values, xi, vector_result);
    vector_answer = vector_field(x);

    BOOST_CHECK( fuzzy_equals( vector_result, vector_answer ) );

    return 0;

}

int test_get_global_shapefunction_gradients( elib::Element &element, elib::vec &local_test_point ){
    /*!
     * Test the computation of the gradient of the shape functions w.r.t. the global coordinates.
     * 
     * :param elib::Element element: The element to be tested
     * :param elib::vec local_test_point: A point in local coordinates to use as the test point.
     */
    
    double eps = 1e-6;

    elib::vec global_test_point;
    element.interpolate(element.nodes, local_test_point, global_test_point);
    
    elib::vec delta = global_test_point;
    elib::vec xtmp, xi, xip, xim, N0, Ntmpp, Ntmpm;
    elib::vecOfvec dNdx_num(element.nodes.size());
    for (uitype n=0; n<dNdx_num.size(); n++){
        dNdx_num[n] = elib::vec(element.nodes[n].size(), 0);
    }

    //Compute the initial value of xi
    element.compute_local_coordinates(global_test_point, xi);

    //Compute the initial values of the shape functions
    element.get_shape_functions(xi, N0);

    for (uitype i=0; i<element.nodes[0].size(); i++){

        //Perturb the global coordinates positively
        delta[i] *= 1+eps;

        //Compute the local coordinates of the positively perturbed global coordinates
        element.compute_local_coordinates(delta, xip);
        
        //Compute the new shape-function values
        element.get_shape_functions(xip, Ntmpp);

        //Remove the positive perturbation
        delta[i] /= 1+eps;

        //Perturb the global coordinates negatively
        delta[i] *= 1-eps;

        //Compute the local coordinates of the negatively perturbed global coordinates
        element.compute_local_coordinates(delta, xim);

        //Compute the new shape-function values
        element.get_shape_functions(xim, Ntmpm);

        //Set the values of the estimated gradient
        for (uitype n=0; n<Ntmpp.size(); n++){dNdx_num[n][i] = (Ntmpp[n] - Ntmpm[n])/(2*global_test_point[i]*eps);}

        //Remove the negative perturbation
        delta[i] /= 1-eps;
    }

    elib::vecOfvec dNdx;
    element.get_global_shapefunction_gradients(local_test_point, dNdx);

    BOOST_CHECK( fuzzy_equals( dNdx_num, dNdx ) );

    return 0;

}

int test_get_local_gradient( elib::Element &element ){
    /*!
    Test the computation of the gradient with respect to the local coordinates

    TODO: Generalize to non-3D in local coordinates elements

    :param elib::Element element: The element to be tested
    */

    double eps = 1e-6;
    elib::vec scalar_answer(3, 0), scalar_result(3, 0);
    double sg0, sgpx, sgpy, sgpz, sgmx, sgmy, sgmz;

    //Form the scalar field at the nodes
    elib::vec scalar_nodal_values(element.nodes.size());
    for (uitype n=0; n<element.nodes.size(); n++){
        scalar_nodal_values[n] = scalar_field(element.nodes[n]);
    }

    //Interpolate the field
    element.interpolate(scalar_nodal_values, {    -0.2,     0.4,     0.64}, sg0);
    element.interpolate(scalar_nodal_values, {-0.2+eps,     0.4,     0.64}, sgpx);
    element.interpolate(scalar_nodal_values, {    -0.2, 0.4+eps,     0.64}, sgpy);
    element.interpolate(scalar_nodal_values, {    -0.2,     0.4, 0.64+eps}, sgpz);
    element.interpolate(scalar_nodal_values, {-0.2-eps,     0.4,     0.64}, sgmx);
    element.interpolate(scalar_nodal_values, {    -0.2, 0.4-eps,     0.64}, sgmy);
    element.interpolate(scalar_nodal_values, {    -0.2,     0.4, 0.64-eps}, sgmz);

    //Compute the numeric gradient
    scalar_answer[0] = (sgpx - sgmx)/(2*eps);
    scalar_answer[1] = (sgpy - sgmy)/(2*eps);
    scalar_answer[2] = (sgpz - sgmz)/(2*eps);

    //Compute the element result
    element.get_local_gradient(scalar_nodal_values, {-0.2, 0.4, 0.64}, scalar_result);

    BOOST_CHECK( fuzzy_equals( scalar_answer, scalar_result ) ); 

    elib::vecOfvec vector_answer, vector_result;
    elib::vecOfvec vector_nodal_values(element.nodes.size());
    elib::vec vg0, vgpx, vgpy, vgpz;

    for (uitype n=0; n<element.nodes.size(); n++){
        vector_nodal_values[n] = vector_field(element.nodes[n]);
    }

    //Interpolate the field
    elib::vec local_coordinates = {-0.2, 0.4, 0.64};
    elib::vec perturbed_coordinates(local_coordinates.size());
    elib::vecOfvec perturbation_matrix(element.local_node_coordinates[0].size()+1);

    //Interpolate at the nominal position
    element.interpolate(vector_nodal_values, local_coordinates, perturbation_matrix[0]);

    //Perturb the local coordinates and interpolate
    for (uitype i=0; i<perturbation_matrix.size()-1; i++){
        for (uitype j=0; j<perturbation_matrix.size()-1; j++){
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
    for (uitype i=0; i<vector_answer.size(); i++){
	vector_answer[i].resize(perturbation_matrix.size()-1);
        for (uitype j=1; j<perturbation_matrix.size(); j++){
            vector_answer[i][j-1] = (perturbation_matrix[j][i] - perturbation_matrix[0][i])/eps;
	}
    }
    element.get_local_gradient(vector_nodal_values, local_coordinates, vector_result);

    BOOST_CHECK( fuzzy_equals(vector_answer, vector_result ) );

    return 0;
}

int test_get_global_gradient( elib::Element &element ){
    /*!
    Test the computation of the global gradient.

    TODO: Generalize to non-3D elements.

    :param elib::Element element: The element to be tested
    */

    //Compute a set of reference coordinates using a linear transformation
    elib::vecOfvec reference_coordinates(element.nodes.size());
    for (uitype n=0; n<element.nodes.size(); n++){
        reference_coordinates[n].resize(element.nodes.size());
        linear_transform(element.nodes[n], reference_coordinates[n]);
    }

    //Evaluate the scalar field
    elib::vec scalar_nodal_current_values(element.nodes.size());
    elib::vec scalar_nodal_reference_values(element.nodes.size());

    for (uitype n=0; n<element.nodes.size(); n++){
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

    BOOST_CHECK( fuzzy_equals( grad_scalar_current, grad_scalar_reference ) && fuzzy_equals( grad_scalar_current, scalar_answer ) );

    //Evaluate the vector field
    elib::vecOfvec vector_nodal_current_values(element.nodes.size());
    elib::vecOfvec vector_nodal_reference_values(element.nodes.size());

    for (uitype n=0; n<element.nodes.size(); n++){
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

    BOOST_CHECK( fuzzy_equals( grad_vector_current, grad_vector_reference ) && fuzzy_equals( grad_vector_current, vector_answer ) );

    return 0;
}

int test_compute_local_coordinates(elib::Element &element, elib::vec xtest, bool isoutside ){
    /*!
    Test for the computation of an element's local coordinates given it's global coordinates.

    TODO: Generalize to non-3D elements.

    :param elib::Element element: The element to be tested
    */

    elib::vec xi, result;

    std::unique_ptr< errorNode > clc_result;
    clc_result.reset( element.compute_local_coordinates(xtest, xi) );

    if ( clc_result ){
        BOOST_CHECK( isoutside );
        return 0;
    }

    element.interpolate( element.nodes, xi, result );

    BOOST_CHECK( fuzzy_equals( result, xtest ) );

    return 0;
}

int test_get_jacobian( elib::Element &element ){
    /*!
    Test the computation of the element's jacobian of transformation.

    TODO: Generalize to non-3D elements.

    :param elib::Element element: The element to be tested
    */

    //Compute a set of reference coordinates using a linear transformation
    elib::vecOfvec reference_coordinates(element.nodes.size());
    for (uitype n=0; n<element.nodes.size(); n++){
        reference_coordinates[n].resize(element.nodes.size());
        linear_transform(element.nodes[n], reference_coordinates[n]);
    }

    elib::vecOfvec result;
    element.get_jacobian({0.2, -0.3, 0.4}, reference_coordinates, result);

    elib::vecOfvec A, answer;
    elib::vec b;
    get_linear_transformation_definition(A, b);
    elib::invert(A, answer);

    BOOST_CHECK( fuzzy_equals( answer, result ) );

    return 0;
}

int test_bounding_box_contains_point( elib::Element &element ){
    /*!
    Test if the bounding box point detection works correctly.

    :param elib::Element element: The element to be tested
    */

    double delta = 0.1;
    elib::vec x;

    x = element.bounding_box[0];
    BOOST_CHECK( element.bounding_box_contains_point( x ) );

    //Check if points smaller than the bounds will be detected as being outside of the element
    for (uitype i=0; i<element.bounding_box[0].size(); i++){
        x = element.bounding_box[0];
        x[i] -= delta;
        BOOST_CHECK( !element.bounding_box_contains_point( x ) );
    } 

    x = element.bounding_box[1];
    BOOST_CHECK( element.bounding_box_contains_point( x ) );

    //Check if points larger than the bounds will be detected as being outside of the element
    for (uitype i=0; i<element.bounding_box[1].size(); i++){
        x = element.bounding_box[1];
        x[i] += delta;
        BOOST_CHECK( !element.bounding_box_contains_point( x ) );
    }

    return 0;
}

int test_contains_point( elib::Element &element ){
    /*!
    Test if the element can identify if a global point is contained inside.

    :param elib::Element element: The element to be tested
    */

    for (uitype n=0; n<element.nodes.size(); n++){
        BOOST_CHECK( element.contains_point( element.nodes[ n ] ) );
    }

    return 0;
}

int test_build_element_from_string( elib::Element &element ){
    /*!
    Test if the build element from string function is working correctly.

    :param elib::Element element: The element to be tested
    */

    std::unique_ptr<elib::Element> new_element = elib::build_element_from_string(element.name, element.global_node_ids, element.nodes, element.qrule);

    std::string result = typeid(*new_element).name();
    std::string answer = typeid(element).name();

    BOOST_CHECK( std::strcmp( result.c_str( ), answer.c_str( ) ) == 0);

    return 0;
}

int test_element_functionality( elib::Element &element, elib::vec &global_test_point, bool isoutside ){
    /*!
    Test the provided element's functionality

    :param elib::Element element: The element to be tested
    :param elib::vec &global_test_point: A point in global space to use for testing
    :param bool isoutside: Flag indicating if the test point is truly outside the element or not.

    */

    test_interpolate(element);
    if (!isoutside){
        test_get_local_gradient(element);
    }
    test_get_global_gradient( element );
    test_compute_local_coordinates( element, global_test_point, isoutside );
    test_get_jacobian( element );
    test_bounding_box_contains_point( element );
    test_contains_point( element );
    test_build_element_from_string( element );
    if ( !isoutside ){
        elib::vec local_test_point;
        errorOut clc_result = element.compute_local_coordinates( global_test_point, local_test_point );
        if ( clc_result ){
            test_get_global_shapefunction_gradients( element, local_test_point );
        }
    }

    return 0;
}

BOOST_AUTO_TEST_CASE( testHex8_functionality ){
    /*!
    Test the Hex8 element's functionality
    */

    // First test a non-distorted element

    // Define the element's nodes
    elib::vecOfvec nodes = {{0, 0, 0},
                            {1, 0, 0},
                            {1, 1, 0},
                            {0, 1, 0},
                            {0, 0, 1},
                            {1, 0, 1},
                            {1, 1, 1},
                            {0, 1, 1}};

    // Define the element's node ids
    std::vector< uitype > node_ids = {1, 2, 3, 4, 5, 6, 7, 8};

    // Define the element's quadrature rule
    elib::quadrature_rule qrule;
    define_hex8_fully_integrated_quadrature(qrule);

    // Construct the element
    elib::Hex8 element(node_ids, nodes, qrule);

    // Define a test point
    elib::vec xtest = {0.25, 0.75, .14};
    bool isoutside = false;

    int tef_return = test_element_functionality( element, xtest, isoutside );

    // Test a distorted element
    nodes = {{0.516905, 0.391528, 0.293894 },
             {0.86161, 0.442245, 0.178099 },
             {1.10153, 0.877418, 0.274955 },
             {0.846862, 0.78123, 0.445236 },
             {0.315421, 0.42434, 0.676207 },
             {0.720471, 0.459122, 0.606603 },
             {0.869162, 0.915384, 0.665252 },
             {0.52575, 0.848709, 0.771187}};

    element = elib::Hex8(node_ids, nodes, qrule);

    xtest = {0.672, 0.636, 0.368};
    isoutside = true;

    tef_return = test_element_functionality( element, xtest, isoutside );

    // Test another distorted element
    nodes = {{3.13443, -0.61357,  1.90472},
             {4.24588,  1.41151,  3.82988},
             {3.97724,  1.34621,  4.43285},
             {2.86579, -0.678866, 2.50769},
             {3.95241, -0.996794, 1.71353},
             {5.06385,  1.02829,  3.63869},
             {4.79521,  0.96299,  4.24166},
             {3.68377, -1.06209,  2.3165}};

    element = elib::Hex8(node_ids, nodes, qrule);

    xtest = {4.38002, 0.56885, 3.65742};
    isoutside = false;

    tef_return = test_element_functionality( element, xtest, isoutside );

}

BOOST_AUTO_TEST_CASE( testHex8_point_on_surface ){
    /*!
     * Test if a point is on the surface of the Hex8 element
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

    // Define the element's node ids
    std::vector< uitype > node_ids = {1, 2, 3, 4, 5, 6, 7, 8};

    // Define the element's quadrature rule
    elib::quadrature_rule qrule;
    define_hex8_fully_integrated_quadrature(qrule);

    // Construct the element
    elib::Hex8 element(node_ids, nodes, qrule);

    //Define the point
    elib::vec point = { 0, 0, 0 };

    std::vector< elib::uitype > answer1 = { 0, 2, 4 };

    std::vector< elib::uitype > result;

    BOOST_CHECK( element.point_on_surface( point, result, 1e-9 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer1 ) );

    point = { -0.1, 0, 0 };

    BOOST_CHECK( !element.point_on_surface( point, result, 1e-9 ) );

    BOOST_CHECK( element.point_on_surface( point, result, 3e-1 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer1 ) );

    point = { -0.1, 0.5, 0.5 };

    BOOST_CHECK( element.point_on_surface( point, result, 3e-1 ) );

    std::vector< elib::uitype > answer2 = { 0 };

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer2 ) );

}

BOOST_AUTO_TEST_CASE( testInvert ){
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
    for (uitype i=0; i<A.size(); i++){
        result[i].resize(A[i].size());
        for (uitype j=0; j<A[i].size(); j++){
            result[i][j] = 0;
            for (uitype k=0; k<A[i].size(); k++){
                result[i][j] += A[i][k]*Ainv[k][j];
            }
        }
    }

    BOOST_CHECK( fuzzy_equals( result, answer ) );

}

BOOST_AUTO_TEST_CASE( testSolve ){
    /*!
    Test the computation of the matrix solve
    */

    elib::vecOfvec A = {{2, 3,  5},
                        {3, 6,  7},
                        {5, 7, 10}};

    elib::vec answer = {1, 2, 3};

    elib::vec b(A.size(), 0);

    for (uitype i=0; i<A.size(); i++){
        for (uitype j=0; j<A[i].size(); j++){
            b[i] += A[i][j]*answer[j];
        }
    }

    elib::vec result;
    elib::solve(A, b, result);

    BOOST_CHECK( fuzzy_equals( answer, result ) );

}

BOOST_AUTO_TEST_CASE( testHex8_transform_local_vector ){
    /*!
     * Test the transformation of a local vector to the global reference frame
     *
     */

    elib::vecOfvec referenceNodes = { { 0, 0, 0 },
                                      { 1, 0, 0 },
                                      { 1, 1, 0 },
                                      { 0, 1, 0 },
                                      { 0, 0, 1 },
                                      { 1, 0, 1 },
                                      { 1, 1, 1 },
                                      { 0, 1, 1 } };

    elib::vecOfvec displacements = { { 0.        ,  0.        ,  0.        },
                                     {-0.81824397,  0.33884637, -1.0510223 },
                                     {-0.71099902,  0.2144174 , -0.81869947},
                                     { 0.10724495, -0.12442897,  0.23232284},
                                     { 0.9034876 ,  0.03718896, -0.43120554},
                                     { 0.08524363,  0.37603533, -1.48222784},
                                     { 0.19248858,  0.25160636, -1.24990501},
                                     { 1.01073255, -0.08724002, -0.1988827 } };

    // Define the element's node ids
    std::vector< uitype > node_ids = {1, 2, 3, 4, 5, 6, 7, 8};

    // Define the element's quadrature rule
    elib::quadrature_rule qrule;
    define_hex8_fully_integrated_quadrature(qrule);

    // Construct the element
    elib::Hex8 element(node_ids, referenceNodes, qrule);

    element.update_node_positions( displacements );

    // local vector 1
    elib::vec local_vector_1 = { 2, 0, 0 };

    elib::vec answer1 = { 1, 0, 0 };

    elib::vec result1;
    element.transform_local_vector( element.local_node_coordinates[ 0 ], local_vector_1, result1, false );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer1, result1 ) );

    elib::vec local_vector_2 = { 2, 2, 2 };
    elib::vec answer2 = { 1, 1, 1 };
    elib::vec result2;

    element.transform_local_vector( element.local_node_coordinates[ 1 ], local_vector_2, result2, false );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer2, result2 ) );

    elib::vec answer3 = referenceNodes[ 1 ] + displacements[ 1 ];
    elib::vec result3;
    element.transform_local_vector( element.local_node_coordinates[ 0 ], local_vector_1, result3, true );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer1, result1 ) );

    elib::vec answer4 = referenceNodes[ 6 ] + displacements[ 6 ];
    elib::vec result4;

    element.transform_local_vector( element.local_node_coordinates[ 0 ], local_vector_2, result4, true );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer4, result4 ) );
}
