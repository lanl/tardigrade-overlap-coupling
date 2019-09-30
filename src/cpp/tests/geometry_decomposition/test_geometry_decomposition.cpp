#include<vector>
#include<iostream>
#include<fstream>
#include<math.h>
#include<geometry_decomposition.h>

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

int test_getTets(std::ofstream &results){
    /*!
     * Test the creation of a collection of tetrahedra that describe 
     * the volume associated with ordered points on a plane that 
     * represent the boundary of a convex polyhedra and some 
     * center point.
     * 
     * :param std::ofstream &results: The output file
     */

    matrixType nodes = {{-1, -1, 1},
                        { 1, -1, 1},
                        { 1,  1, 1},
                        {-1,  1, 1}};

    vectorType centroid = {0, 0, 0};

    std::vector< matrixType > tets = gDecomp::getTets(centroid, nodes);

    unsigned int i=0, j=1;
    for (auto tet = tets.begin(); tet!=tets.end(); tet++, i++, j++){
        if (!fuzzy_equals((*tet)[0], centroid)){
            results << "test_getTets (test 1) & False\n";
            return 0;
        }
        if (!fuzzy_equals((*tet)[1], {0, 0, 1})){
            results << "test_getTets (test 2) & False\n";
            return 0;
        }
        if (!fuzzy_equals((*tet)[2], nodes[i])){
            results << "test_getTets (test 3) & False\n";
        }
        if (j >= nodes.size()){
            j = 0;
        }
        if (!fuzzy_equals((*tet)[3], nodes[j])){
            results << "test_getTets (test 4) & False\n";
        }
    }

    results << "test_getTets & True\n";
    return 0;
}

int test_getTetVolume(std::ofstream &results){
    /*!
     * Compute the volume of a tetrahedron
     * 
     * :param std::ofstream &results: The output file
     */

    matrixType tet = {{ 0, 0, 0},
                      { 1, 0, 0},
                      { 0, 1, 0},
                      { 0, 0, 1}};

    double volume = gDecomp::getTetVolume(tet);
    
    if (!fuzzy_equals(volume, 1./6)){
        results << "test_getTetVolume (test 1) & False\n";
        return 1;
    }

    results << "test_getTetVolume & True\n";
    return 0;
}

int test_getUnitToTetMap(std::ofstream &results){
    /*!
     * Test the computation of the map between the unit tetrahedra and an 
     * arbitrary tetrahedron.
     * 
     * :param const matrixType &nodes: The nodes of the arbitrary tetrahedron.
     * :param matrixType &A: The mapping matrix A
     * :param vectorType &d: The shift vector d
     */

    matrixType localNodes = {{ 0, 0, 0},
                             { 1, 0, 0},
                             { 0, 1, 0},
                             { 0, 0, 1}};

    matrixType nodes = {{  1,  4, 2},
                        {  6,  4, 1},
                        { 10,  3, 5},
                        {  3, -1, 4}};

    matrixType A;
    vectorType d;
    gDecomp::getUnitToTetMap(nodes, A, d);

    unsigned int i=0;
    for (auto ln=localNodes.begin(); ln!=localNodes.end(); ln++, i++){
        if (!fuzzy_equals(nodes[i], vectorTools::dot(A, *ln) + d)){
            results << "test_getUnitToTetMap (test 1) & False\n";
            return 1;
        }
    }

    results << "test_getUnitToTetMap & True\n";
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

    //Test the utility functions
    test_getTets(results);
    test_getTetVolume(results);
    test_getUnitToTetMap(results);

    //Close the results file
    results.close();

    return 0;

}
