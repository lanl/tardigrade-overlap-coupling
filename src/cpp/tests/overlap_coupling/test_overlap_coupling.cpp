/*!
Unit tests for overlap_coupling.cpp
*/

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>

#include "overlap_coupling.h"

template< typename T1, typename T2 >
bool double_compare(const T1 &a, const T2 &b, const double tolr = 1e-9, const double tola = 1e-9){
    /*!
    Compare two floats to see if the are equal
    */

    double tol = fmin(tolr*a + tola, tolr*b + tola);
    return tol>(a - b);
}

void test_vector_operators(std::ofstream &results){
    /*!
    Test the overloaded operators in Vector
    */

    std::vector< double > a(3), a2(3), b(3);
    a[0] = 1;
    a[1] = 2;
    a[2] = 5;
    a2[0] = a[0];
    a2[1] = a[1];
    a2[2] = a[2];
    b[0] = 4;
    b[1] = -1;
    b[2] = 3.5;

    overlap::Vector v1(a);
    overlap::Vector v2(a2);
    overlap::Vector v3(b);

    //Test element access
    double result = v1(2);
    bool answer = double_compare(result, 5);

    if (!answer){
        results << "test_vector_operators (component access) & False\n";
        return;
    }

    //Test vector comparison
    answer = (v1 == v2);
    if (!answer){
        results << "test_vector_operators (equality 1) & False\n";
        return;
    }

    answer = !(v1 == v3);
    if (!answer){
        results << "test_vector_operators (equality 2) & False\n";
        return;
    }

    //Test addition
    overlap::Vector r = v1 + v2;
    std::vector< double > _va(3);
    _va[0] = a[0] + a2[0];
    _va[1] = a[1] + a2[1];
    _va[2] = a[2] + a2[2];
    overlap::Vector va(_va);

    answer = (va == r);
    if (!answer){
        results << "test_vector_operators (addition 1) & False\n";
        return;
    }

    r = v1 + v3;
    _va[0] = a[0] + b[0];
    _va[1] = a[1] + b[1];
    _va[2] = a[2] + b[2];

    va = overlap::Vector(_va);
    
    answer = (va == r);
    if (!answer){
        results << "test_vector_operators (addition 2) & False\n";
        return;
    }

    _va[0] = a[0] + .42;
    _va[1] = a[1] + .42;
    _va[2] = a[2] + .42;

    va = overlap::Vector(_va);
    r = .42 + v1;

    answer = (va == r);

    if (!answer){
        results << "test_vector_operators (addition 3) & False\n";
        return;
    }

    _va[0] = 1.2+a[0];
    _va[1] = 1.2+a[1];
    _va[2] = 1.2+a[2];
    va = overlap::Vector(_va);
    r = a;
    r += 1.2;

    answer = (va == r);
    if (!answer){
        results << "test_vector_operators (addition 4) & False\n";
        return;
    }

    _va[0] = a[0] + b[0];
    _va[1] = a[1] + b[1];
    _va[2] = a[2] + b[2];
    va = overlap::Vector(_va);
    r = a;
    r += b;

    answer = (va == r);
    if (!answer){
        results << "test_vector_operators (addition 5) & False\n";
        return;
    }

    //Test negation
    _va[0] = -a[0];
    _va[1] = -a[1];
    _va[2] = -a[2];
    va = overlap::Vector(_va);

    r = -v1;
    answer = (va == r);

    if (!answer){
        results << "test_vector_operators (negation) & False\n";
        return;
    }

    //Test subtraction
    _va[0] = a[0] - b[0];
    _va[1] = a[1] - b[1];
    _va[2] = a[2] - b[2];
    va = overlap::Vector(_va);

    r = v1 - v3;

    answer = (va == r);
    if (!answer){
        results << "test_vector_operators (subtraction) & False\n";
        return;
    }

    //Test multiplication
    _va[0] = 3.4*a[0];
    _va[1] = 3.4*a[1];
    _va[2] = 3.4*a[2];

    va = overlap::Vector(_va);
    r = 3.4*v1;
    
    answer = (va == r);
    if (!answer){
        results << "test_vector_operators (multiplication 1) & False\n";
        return;
    }

    _va[0] = a[0]*b[0];
    _va[1] = a[1]*b[1];
    _va[2] = a[2]*b[2];

    va = overlap::Vector(_va);
    r = v1*v3;
    
    answer = (va == r);
    if (!answer){
        results << "test_vector_operators (multiplication 2) & False\n";
        return;
    }

    _va[0] = 1.2*a[0];
    _va[1] = 1.2*a[1];
    _va[2] = 1.2*a[2];
    va = overlap::Vector(_va);
    r = a;
    r *= 1.2;

    answer = (va == r);
    if (!answer){
        results << "test_vector_operators (multiplication 3) & False\n";
        return;
    }

    _va[0] = a[0]/7.5;
    _va[1] = a[1]/7.5;
    _va[2] = a[2]/7.5;
    va = overlap::Vector(_va);
    r = v1/7.5;

    answer = (va == r);
    if (!answer){
        results << "test_vector_operators (division) & False\n";
        return;
    }

    //All tests passed!
    results << "test_vector_operators & True\n";
}

void test_vector_methods(std::ofstream &results){
    /*!
    Test the methods implemented in vector
    */

    //!Test the dyadic product
    std::vector< double > a(3);
    std::vector< double > b(5);
    a[0] = 1.3;
    a[1] = 3.4;
    a[2] = -2.1;

    b[0] = .12;
    b[1] = -3.42;
    b[2] = 1.42;
    b[3] = 9.21;
    b[4] = 1.2;

    overlap::Vector v1 = overlap::Vector(a);
    overlap::Vector v2 = overlap::Vector(b);

    std::vector< overlap::Vector > dpr = v1.dyadic_product(v2);

    std::vector< std::vector< double > > dpvec_answer(a.size());
    std::vector< overlap::Vector > dpanswer(a.size());
    for (unsigned int i=0; i<a.size(); i++){
        dpvec_answer[i].resize(b.size());
        for (unsigned int j=0; j<b.size(); j++){
            dpvec_answer[i][j] = a[i]*b[j];
        }
        dpanswer[i] = overlap::Vector(dpvec_answer[i]);

        if (!(dpanswer[i] == dpr[i])){
            results << "test_vector_methods (dyadic product) & False\n";
            return;
        }
    }

    //!Test the sum
    double rdouble = v1.sum();
    
    if (!overlap::fuzzy_compare(rdouble, 1.3+3.4-2.1)){
        results << "test_vector_methods (sum) & False\n";
        return;
    }

    //!Test the product
    rdouble = v1.product();

    if (!overlap::fuzzy_compare(rdouble, -1.3*3.4*2.1)){
        results << "test_vector_methods (product) & False\n";
        return;
    }

    results << "test_vector_methods & True\n";
    return;
}

void test_Hex8(std::ofstream &results){
    /*!
    Test the implementation of the linearly integrated isoparametric hexahedron element.
    */

    //Define the global coordinates of the nodes
    std::vector< overlap::Vector > global_coordinates(8);

    std::vector< double > vec(3,0);
    global_coordinates[0] = overlap::Vector(vec);

    vec[0] = 1;
    global_coordinates[1] = overlap::Vector(vec);

    vec[1] = 2;
    global_coordinates[2] = overlap::Vector(vec);

    vec[0] = 0;
    global_coordinates[3] = overlap::Vector(vec);

    vec[1] = 0;
    vec[2] = 3;
    global_coordinates[4] = overlap::Vector(vec);

    vec[0] = 1;
    global_coordinates[5] = overlap::Vector(vec);

    vec[1] = 2;
    global_coordinates[6] = overlap::Vector(vec);

    vec[0] = 0;
    global_coordinates[7] = overlap::Vector(vec);

    overlap::Hex8 element(global_coordinates);

    //!Check the shape functions
    vec = std::vector< double >(3, 0);
    overlap::Vector position(vec);
    double dresult;
    for (unsigned int m=0; m<global_coordinates.size(); m++){
        position = element.get_local_coordinates(m);
        for (unsigned int n=0; n<global_coordinates.size(); n++){
            dresult = element.shape_function(n, position);
            if (m==n){
                if (!overlap::fuzzy_compare(dresult, 1)){
                    results << "test_Hex8 (shape_functions) & False\n";
                    return;
                }
            }
            else{
                if (!overlap::fuzzy_compare(dresult, 0)){
                    results << "test_Hex8 (shape_functions) & False\n";
                    return;
                }
            }
        }
    }
    

    //!Check the interpolate function
    overlap::Vector Vec_result;
    position = std::vector< double >(3, 0);
    element.interpolate(global_coordinates, position, Vec_result);

    vec[0] = 0.5;
    vec[1] = 1.0;
    vec[2] = 1.5;
    overlap::Vector Vec_answer(vec);

    if (!(Vec_result == Vec_answer)){
        results << "test_Hex8 (interpolate) & False\n";
        Vec_result.print();
        Vec_answer.print();
        return;
    }

    results << "test_Hex8 & True\n";
}

int main(){
    /*!
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or 
    False if the test passes or fails respectively.
    */

    //Open the results file
    std::ofstream results;
    results.open ("results.tex");

    test_vector_operators(results);
    test_vector_methods(results);

    //!Test the supported element types
    test_Hex8(results);

    //Close the results file
    results.close();
}
