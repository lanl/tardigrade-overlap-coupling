//!The test file for dataFileInterface.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>
#define USE_EIGEN
#include<vector_tools.h>

#include<dataFileInterface.h>

typedef dataFileInterface::errorNode errorNode; //!Redefinition for the error node
typedef dataFileInterface::errorOut errorOut; //!Redefinition for a pointer to the error node
typedef dataFileInterface::floatType floatType; //!Define the float values type.
typedef dataFileInterface::floatVector floatVector; //! Define a vector of floats
typedef dataFileInterface::floatMatrix floatMatrix; //!Define a matrix of floats
typedef dataFileInterface::uIntVector uIntVector; //!Define a vector of unsigned ints

int test_XDMFDataFile_constructor( std::ofstream &results ){
    /*!
     * Test the interface with the XDMF file format
     * constructor
     *
     * :param std::ofstream &results: The output file
     */

    std::unique_ptr<dataFileInterface::dataFileBase> df;
    df = dataFileInterface::dataFileBase().Create( "XDMF" );

    if ( !df->_error ){
        df->_error->print();
        results << "test_XDMFDataFile_constructor & False\n";
        return 1;
    }

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    df = dataFileInterface::dataFileBase( yf["filetest1"] ).Create( "XDMF" );

    if ( df->_error ){
        df->_error->print();
        results << "test_XDMFDataFile_constructor & False\n";
        return 1;
    }

    if ( df->_filename.compare( "macroscale_xdmf.xdmf" ) != 0 ){
        results << "test_XDMFDataFile_constructor (test 3) & False\n";
        return 1;
    }

    if ( df->_mode.compare( "read" ) != 0 ){
        results << "test_XDMFDataFile_constructor (test 4) & False\n";
        return 1;
    }

    df = dataFileInterface::dataFileBase( yf[ "filetest2" ] ).Create( "XDMF" );

    if ( !df->_error ){
        results << "test_XDMFDataFile_constructor & False\n";
        return 1;
    }

    results << "test_XDMFDataFile_constructor & True\n";
    return 0;
}

int test_XDMFDataFile_readMesh( std::ofstream &results ){
    /*!
     * Test the interface with the XDFM file format
     * in a read capacity.
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdf( yf[ "filetest1" ] );

    floatVector nodePositionsAnswer = { 1, 0, 1,
                                        1, 0, 0,
                                        0, 0, 0,
                                        0, 0, 1,
                                        1, 1, 1,
                                        1, 1, 0,
                                        0, 1, 0,
                                        0, 1, 1,
                                        0, 1, 2,
                                        1, 1, 2,
                                        0, 0, 2,
                                        1, 0, 2,
                                        0, 0, 3,
                                        0, 1, 3,
                                        1, 1, 3,
                                        1, 0, 3 };

    floatVector nodePositionsResult;

    errorOut error = xdf.readMesh( 1, nodePositionsResult );

    if ( error ){
        error->print( );
        results << "test_XDMFDataFile_readMesh & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( nodePositionsAnswer, nodePositionsResult ) ){
        results << "test_XDMFDataFile_readMesh (test 1) & False\n";
        return 1;
    }

    results << "test_XDMFDataFile_readMesh & True\n";
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

    test_XDMFDataFile_constructor( results );
    test_XDMFDataFile_readMesh( results );

    //Close the results file
    results.close();
}
