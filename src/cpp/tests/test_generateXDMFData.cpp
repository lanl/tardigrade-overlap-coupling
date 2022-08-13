//!The test file for dataFileInterface.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>
#define USE_EIGEN
#include<vector_tools.h>

#include<generateXDMFData.h>

#define BOOST_TEST_MODULE test_dataFileInterface
#include <boost/test/included/unit_test.hpp>

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

BOOST_AUTO_TEST_CASE( testfileGenerator_constructor ){
    /*!
     * Test the generateXDMFData constructor
     *
     */

    remove( "xdmf_out.xdmf" );

    //Make sure the default constructor runs
    fileGenerator::fileGenerator fG;

    BOOST_CHECK( !fG.getError( ) );

    fG = fileGenerator::fileGenerator( "bad_file" );

    BOOST_CHECK( fG.getError( ) );

    fG = fileGenerator::fileGenerator( "generateXDMFData_testYAML.yaml" );

    BOOST_CHECK( !fG.getError( ) );

    std::ifstream output_file;
    output_file.open( "xdmf_out.xdmf" );

    BOOST_CHECK( output_file.good( ) );

    remove( "xdmf_out.xdmf" );
    remove( "xdmf_out.h5" );

}

BOOST_AUTO_TEST_CASE( testfileGenerator_build ){
    /*!
     * Test the function which builds the XDMF file
     *
     * :param std::ofstream &results: The output file
     */

    remove( "xdmf_out.xdmf" );
    remove( "xdmf_out.h5" );

    fileGenerator::fileGenerator fG( "generateXDMFData_testYAML.yaml" );

    BOOST_CHECK( !fG.getError( ) );

    BOOST_CHECK( !fG.build( ) );

    BOOST_CHECK( *fG.getCurrentIncrement( ) == 1 );

    BOOST_CHECK( compareFiles( "xdmf_out.xdmf", "generateXDMFData_xdmf_answer.xdmf" ) );

    remove( "xdmf_out.xdmf" );
    remove( "xdmf_out.h5" );

    //Tests for errors
    fG = fileGenerator::fileGenerator( "generateXDMFData_badYAML.yaml" );

    BOOST_CHECK( !fG.getError( ) );

    BOOST_CHECK( fG.build( ) );

    remove( "xdmf_out.xdmf" );
    remove( "xdmf_out.h5" );

}
