//!The test file for dataFileInterface.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>
#define USE_EIGEN
#include<vector_tools.h>

#include<generateXDMFData.h>

typedef dataFileInterface::errorNode errorNode; //!Redefinition for the error node
typedef dataFileInterface::errorOut errorOut; //!Redefinition for a pointer to the error node
typedef dataFileInterface::floatType floatType; //!Define the float values type.
typedef dataFileInterface::floatVector floatVector; //! Define a vector of floats
typedef dataFileInterface::floatMatrix floatMatrix; //!Define a matrix of floats
typedef dataFileInterface::uIntType uIntType; //!Define the unsigned int type
typedef dataFileInterface::uIntVector uIntVector; //!Define a vector of unsigned ints
typedef dataFileInterface::stringVector stringVector; //!Define a vector of strings

int test_fileGenerator_constructor( std::ofstream &results ){
    /*!
     * Test the generateXDMFData constructor
     *
     * :param std::ofstream &results: The output file
     */

    //Make sure the default constructor runs
    fileGenerator::fileGenerator fG;

    if ( fG.getError( ) ){

        fG.getError( )->print( );
        results << "test_fileGenerator_constructor (test 1) & False\n";
        return 1;

    }

    fG = fileGenerator::fileGenerator( "bad_file" );

    if ( !fG.getError( ) ){
        results << "test_fileGenerator_constructor (test 2) & False\n";
        return 1;
    }

    fG = fileGenerator::fileGenerator( "testYAML.yaml" );

    if ( fG.getError( ) ){
        fG.getError( )->print( );
        results << "test_fileGenerator_constructor (test 3) & False\n";
        return 1;
    }

    results << "test_fileGenerator_constructor & True\n";
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

    test_fileGenerator_constructor( results );

    //Close the results file
    results.close();
}
