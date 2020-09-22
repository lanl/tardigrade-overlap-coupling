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

bool compareFiles( const std::string &fn1, const std::string &fn2 ){
    /*!
     * Compare two files to determine if they are exactly the same
     *
     * :param const std::string &fn1: The first filename
     * :param const std::string &fn2: The second filename
     */

    std::ifstream f1, f2;
    f1.open( fn1 );
    f2.open( fn2 );

    //Check if the files can be opened
    if ( ( !f1.good( ) ) || ( !f2.good( ) ) ){

        return false;

    }

    //Check if the files are the same size
    if ( f1.tellg( ) != f2.tellg( ) ){

        return false;

    }

    //Rewind the files
    f1.seekg( 0 );
    f2.seekg( 0 );

    //Compare the lines
    std::istreambuf_iterator<char> begin1(f1);
    std::istreambuf_iterator<char> begin2(f2);

    return std::equal(begin1,std::istreambuf_iterator<char>(),begin2); //Second argument is end-of-range iterator

}

int test_fileGenerator_constructor( std::ofstream &results ){
    /*!
     * Test the generateXDMFData constructor
     *
     * :param std::ofstream &results: The output file
     */

    remove( "xdmf_out.xdmf" );

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

    std::ifstream output_file;
    output_file.open( "xdmf_out.xdmf" );

    if ( !output_file.good( ) ){
        results << "test_fileGenerator_constructor (test 4) & False\n";
        return 1;
    }

    remove( "xdmf_out.xdmf" );

    results << "test_fileGenerator_constructor & True\n";
    return 0;
}

int test_fileGenerator_build( std::ofstream &results ){
    /*!
     * Test the function which builds the XDMF file
     *
     * :param std::ofstream &results: The output file
     */

    remove( "xdmf_out.xdmf" );
    remove( "xdmf_out.h5" );

    fileGenerator::fileGenerator fG( "testYAML.yaml" );

    if ( fG.getError( ) ){
        fG.getError( )->print( );
        results << "test_fileGenerator_build & False\n";
        return 1;
    }

    if ( fG.build( ) ){

        fG.getError( )->print( );
        results << "test_fileGenerator_build & False\n";
        return 1;

    }

    if ( *fG.getCurrentIncrement( ) != 1 ){

        results << "test_fileGenerator_build (test 1) & False\n";
        return 1;

    }

    if ( !compareFiles( "xdmf_out.xdmf", "xdmf_answer.xdmf" ) ){

        results << "test_fileGenerator_build (test 2) & False\n";
        return 1;

    }

    remove( "xdmf_out.xdmf" );
    remove( "xdmf_out.h5" );

    //Tests for errors
    fG = fileGenerator::fileGenerator( "badYAML.yaml" );

    if ( fG.getError( ) ){

        fG.getError( )->print( );
        results << "test_fileGenerator_build & False\n";
        return 1;

    }

    if ( !fG.build( ) ){

        results << "test_fileGenerator_build (test 3) & False\n";
        return 1;

    }

    remove( "xdmf_out.xdmf" );
    remove( "xdmf_out.h5" );

    results << "test_fileGenerator_build & True\n";
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
    test_fileGenerator_build( results );

    //Close the results file
    results.close();
}
