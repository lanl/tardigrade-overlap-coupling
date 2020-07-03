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

    YAML::Node yf = YAML::LoadFile( "reconstructionConfiguration.yaml" );
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

    YAML::Node yf = YAML::LoadFile( "reconstructionConfiguration.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_loadPoints (test 1) & False\n";
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

    //Close the results file
    results.close();
}
