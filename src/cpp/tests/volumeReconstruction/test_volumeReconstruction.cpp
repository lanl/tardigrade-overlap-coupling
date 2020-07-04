//!The test file for volumeReconstruction.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>
#define USE_EIGEN
#include<vector_tools.h>

#include<volumeReconstruction.h>

typedef volumeReconstruction::errorNode errorNode; //!Redefinition for the error node
typedef volumeReconstruction::errorOut errorOut; //!Redefinition for a pointer to the error node
typedef volumeReconstruction::floatType floatType; //!Define the float values type.
typedef volumeReconstruction::floatVector floatVector; //! Define a vector of floats
typedef volumeReconstruction::floatMatrix floatMatrix; //!Define a matrix of floats
typedef volumeReconstruction::uIntVector uIntVector; //!Define a vector of unsigned ints

int test_dualContouring_constructor( std::ofstream &results ){
    /*!
     * Tests for the constructors of dualContouring
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError( ) ){
        
        dc.getError( )->print( );

        results << "test_dualContouring_constructor (test 1) & False\n";
        return 1;
    
    }

    std::shared_ptr< volumeReconstruction::volumeReconstructionBase > yR
        = volumeReconstruction::volumeReconstructionBase( yf ).create( );

    if ( yR->getError( ) ){

        yR->getError( )->print( );

        results << "test_dualContouring_constructor (test 2) & False\n";
        return 1;

    }

   results << "test_dualContouring_constructor & True\n";
   return 0;
}

int test_dualContouring_loadPoints( std::ofstream &results ){
    /*!
     * Test for loading the points into the object
     */

    floatVector points = { 1, 2, 3, 4, 5, 6 };

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_loadPoints & False\n";
        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){

        error->print( );
        results << "test_dualContouring_loadPoints (test 1) & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( *dc.getPoints( ), points ) ){

        results << "test_dualContouring_loadPoints (test 2) & False\n";
        return 1;

    }

    floatVector points2 = { 1, 2, 3, 4, 5 };

    error = dc.loadPoints( &points2 );

    if ( !error ){

        error->print( );
        results << "test_dualContouring_loadPoints (test 3) & False\n";
        return 1;

    }

    results << "test_dualContouring_loadPoints & True\n";
    return 0;
}

int test_dualContouring_loadFunction( std::ofstream &results ){
    /*!
     * Test for loading the function into the object
     */

    floatVector points = { 1, 2, 3, 4, 5, 6 };

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_loadFunction & False\n";
        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){

        error->print( );
        results << "test_dualContouring_loadFunction & False\n";
        return 1;

    }

    floatVector function = { -1, -2, 1, 7, 8, 10 };

    error = dc.loadFunction( &function );

    if ( error ){

        error->print( );
        results << "test_dualContouring_loadFunction & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( *dc.getFunction( ), function ) ){

        results << "test_dualContouring_loadFunction (test 1) & False\n";
        return 1;

    }

    floatVector function2 = { 2 };

    if ( !dc.loadFunction( &function2 ) ){

        results << "test_dualContouring_loadFunction (test 2) & False\n";
        return 1;

    }

    results << "test_dualContouring_loadFunction & True\n";
    return 0;
}

int test_KDNode_constructor( std::ofstream &results ){
    /*!
     * Test the KD tree creation
     *
     * :param std::ofstream &results: The output file
     */

    floatVector points = { 2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2 };
    uIntVector  ownedIndices = { 0, 2, 4, 6, 8, 10 };
    unsigned int dim = 2;

    volumeReconstruction::KDNode tree( &points, ownedIndices, 0, dim );

    results << "test_KDNode_constructor & True\n";
    return 0;
}

int test_KDNode_getIndex( std::ofstream &results ){
    /*!
     * Test getting the index of the KDNode
     *
     * :param std::ofstream &results: The output file
     */

    floatVector points = { 2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2 };
    uIntVector  ownedIndices = { 0, 2, 4, 6, 8, 10 };
    unsigned int dim = 2;

    volumeReconstruction::KDNode tree( &points, ownedIndices, 0, dim );

    if ( !( tree.getIndex( ), 10 ) ){
        
        results << "test_KDNode_getIndex (test 1) & True\n";
        return 1;

    }

    results << "test_KDNode_getIndex & True\n";
    return 0;
}

int test_KDNode_getMinimumValueDimension( std::ofstream &results ){
    /*!
     * Test getting the upper bound of the KD tree
     *
     * :param std::ofstream &results: The output file
     */

    floatVector points = { 2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2 };
    uIntVector  ownedIndices = { 0, 2, 4, 6, 8, 10 };
    unsigned int dim = 2;

    volumeReconstruction::KDNode tree( &points, ownedIndices, 0, dim );

    floatType answer1 = 2;
    floatType answer2 = 1;

    floatType result = tree.getMinimumValueDimension( 0 );

    if ( !vectorTools::fuzzyEquals( result, answer1 ) ){
        results << "test_KDNode_getMinimumValueDimension (test 1) & False\n";
        return 1;
    }

    result = tree.getMinimumValueDimension( 1 );

    if ( !vectorTools::fuzzyEquals( result, answer2 ) ){
        results << "test_KDNode_getMinimumValueDimesnion (test 2) & False\n";
        return 1;
    }

    floatVector points2 = { 1, 10, 50, 50, 10, 30, 35, 90, 55, 1, 60, 80, 25, 40, 70, 70, 51, 75 };
    uIntVector ownedIndices2 = { 0, 2, 4, 6, 8, 10, 12, 14, 16 };

    volumeReconstruction::KDNode tree2( &points2, ownedIndices2, 0, dim );

    floatType answer3 = 1;
    floatType answer4 = 1;

    result = tree2.getMinimumValueDimension( 0 );

    if ( !vectorTools::fuzzyEquals( result, answer3 ) ){
        results << "test_KDNode_getMinimumValueDimension (test 3) & False\n";
        return 1;
    }

    result = tree2.getMinimumValueDimension( 1 );

    if ( !vectorTools::fuzzyEquals( result, answer4 ) ){
        results << "test_KDNode_getMinimumValueDimension (test 4) & False\n";
        return 1;
    }

    results << "test_KDNode_getMinimumValueDimension & True\n";
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

    test_dualContouring_constructor( results );
    test_dualContouring_loadPoints( results );
    test_dualContouring_loadFunction( results );
    test_KDNode_constructor( results );
    test_KDNode_getIndex( results );
    test_KDNode_getMinimumValueDimension( results );

    //Close the results file
    results.close();
}
