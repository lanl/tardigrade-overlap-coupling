#include<vector>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector_tools.h>

typedef double floatType;
typedef std::vector< floatType > vectorType;
typedef std::vector< vectorType > matrixType;

bool fuzzy_equals(double a, double b, double tolr=1e-6, double tola=1e-6){
    /*!
    Compare two doubles to determine if they are equal.
    */

    floatType tol = fmin(tolr*fabs(a) + tola, tolr*fabs(b) + tola);
    return fabs(a-b)<tol;
}

bool fuzzy_equals(vectorType a, vectorType b, double tolr=1e-6, double tola=1e-6){
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

bool fuzzy_equals(matrixType A, matrixType B, double tolr=1e-6, double tola=1e-6){
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

void print(vectorType a){
    /*!
    Print the vector to the terminal
    */

    for (unsigned int i=0; i<a.size(); i++){
        std::cout << a[i] << " ";
    }
    std::cout << "\n";
}

void print(matrixType A){
    /*!
    Print the matrix to the terminal
    */

    for (unsigned int i=0; i<A.size(); i++){
        print(A[i]);
    }
}

int test_addition_operators(std::ofstream &results){
    /*!
     * Test the addition operators
     * 
     * :param std::ofstream &results: The output file
     */

    vectorType a = { 1, 2, 3};
    vectorType b = {-2, 7, 2};
    vectorType c;

    a += b;

    if (!fuzzy_equals(a, {-1, 9, 5})){
        results << "test_addition_operators (test 1) & False\n";
        return 1;
    }

    c = a + b;

    if (!fuzzy_equals(c, {-3, 16, 7})){
        results << "test_addition_operators (test 2) & False\n";
        return 1;
    }

    //All tests passed    
    results << "test_addition_operators & True\n";
    return 0;
}

int test_subtraction_operators(std::ofstream &results){
    /*!
     * Test the subtraction operators
     * 
     * :param std::ofstream &results: The output file
     */

    vectorType a = { 1, 2, 3};
    vectorType b = {-2, 7, 2};
    vectorType c;

    if (!fuzzy_equals(-a, {-1, -2, -3})){
        results << "test_subtraction_operators (test 1) & False\n";
        return 1;
    }

    a -= b;

    if (!fuzzy_equals(a, {3, -5, 1})){
        results << "test_subtraction_operators (test 2) & False\n";
        return 1;
    }

    c = a - b;

    if (!fuzzy_equals(c, {5,-12, -1})){
        results << "test_subtraction_operators (test 3) & False\n";
        return 1;
    }

    //All tests passed    
    results << "test_subtraction_operators & True\n";
    return 0;
}

int test_multiplication_operators(std::ofstream &results){
    /*!
     * Test the multiplication operators
     * 
     * :param std::ofstream &results: The output file
     */

    vectorType a = { 1, 2, 3};
    vectorType b, c;


    a *= 2;

    if (!fuzzy_equals(a, {2, 4, 6})){
        results << "test_multiplication_operators (test 1) & False\n";
        return 1;
    }

    b = 3*a;
    c = a*3;

    if (!fuzzy_equals(b, c) && !fuzzy_equals(b, {6, 12, 18})){
        results << "test_multiplication_operators (test 2) & False\n";
        return 1;
    }

    //All tests passed    
    results << "test_multiplication_operators & True\n";
    return 0;
}

int test_division_operators(std::ofstream &results){
    /*!
     * Test the division operators
     * 
     * :param std::ofstream &results: The output file
     */

    vectorType a = { 1, 2, 3};
    vectorType b;


    a /= 2;

    if (!fuzzy_equals(a, {0.5, 1, 1.5})){
        results << "test_division_operators (test 1) & False\n";
        return 1;
    }

    b = a/2;

    if (!fuzzy_equals(b, {0.25, 0.5, 0.75})){
        results << "test_division_operators (test 2) & False\n";
        return 1;
    }

    //All tests passed    
    results << "test_division_operators & True\n";
    return 0;
}

int test_computeMean(std::ofstream &results){
    /*!
     * Test the computation of the mean of a vector of vectors
     * 
     * :param std::ofstream &results: The output file
     */

    matrixType A = {{ 1,  2, 3.0, 4},
                    {-4, 13, 0.4, 5},
                    { 2,  6, 1.0, 7}};

    vectorType answer = {-1./3, 7, 8.8/6, 5 + 1/3.};
    vectorType result;
    vectorTools::computeMean(A, result);

    if (!fuzzy_equals(result, answer)){
        results << "test_computeMean (test 1) & False\n";
        return 1;
    }

    results << "test_computeMean & True\n";
    return 0;
}

int test_cross(std::ofstream &results){
    /*!
     * Test the computation of the cross product of two vectors
     * 
     * :param std::ofstream &results: The output file
     */

    vectorType a = { 1, 2};
    vectorType b = {-1, 7};
    vectorType c;

    vectorTools::cross(a, b, c);

    if (!fuzzy_equals(c, {0, 0, 9})){
        results << "test_cross (test 1) & False\n";
        return 1;
    }

    a = {1, 2, 3};
    b = {-1, 7, -3};

    vectorTools::cross(a, b, c);

    if (!fuzzy_equals(c, {-27, 0, 9})){
        results << "test_cross (test 2) & False\n";
        return 1;
    }

    results << "test_cross & True\n";
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

    //Test the operator overloading
    test_addition_operators(results);
    test_subtraction_operators(results);
    test_multiplication_operators(results);
    test_division_operators(results);

    //Test the utility functions
    test_computeMean(results);
    test_cross(results);

    //Close the results file
    results.close();

    return 0;

}
