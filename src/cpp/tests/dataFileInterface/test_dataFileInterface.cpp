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

    std::shared_ptr<dataFileInterface::dataFileBase> df;
    df = dataFileInterface::dataFileBase().create( "XDMF" );

    if ( !df->_error ){
        df->_error->print();
        results << "test_XDMFDataFile_constructor & False\n";
        return 1;
    }

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    df = dataFileInterface::dataFileBase( yf["filetest1"] ).create( "XDMF" );

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

    if ( df->_mode.compare( "read" ) != 0 ){
        results << "test_XDMFDataFile_constructor (test 4) & False\n";
        return 1;
    }

    df = dataFileInterface::dataFileBase( yf["filetest1"] ).create( );

    if ( !df ){
        results << "test_XDMFDataFile_constructor (NULL) & False\n";
        return 1;
    }

    if ( df->_error ){
        df->_error->print();
        results << "test_XDMFDataFile_constructor & False\n";
        return 1;
    }

    if ( df->_filename.compare( "macroscale_xdmf.xdmf" ) != 0 ){
        results << "test_XDMFDataFile_constructor (test 5) & False\n";
        return 1;
    }

    if ( df->_mode.compare( "read" ) != 0 ){
        results << "test_XDMFDataFile_constructor (test 6) & False\n";
        return 1;
    }

    df = dataFileInterface::dataFileBase( yf[ "filetest2" ] ).create( "XDMF" );

    if ( !df->_error ){
        results << "test_XDMFDataFile_constructor & False\n";
        return 1;
    }

    results << "test_XDMFDataFile_constructor & True\n";
    return 0;
}

int test_XDMFDataFile_readMesh( std::ofstream &results ){
    /*!
     * Test the interface with the mesh for the XDMF file format.
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

int test_XDMFDataFile_getNumIncrements( std::ofstream &results ){
    /*!
     * Test the interface with the XDMF file to get the number of
     * temporal increments.
     *
     * :param std::ofstream &results: The output file.
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdf( yf[ "filetest1" ] );

    unsigned int numIncrementsAnswer = 11;
    unsigned int numIncrementsResult;

    errorOut error = xdf.getNumIncrements( numIncrementsResult );
    
    if ( error ){
        error->print( );
        results << "test_XDMFDataFile_getNumIncrements (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( numIncrementsResult, numIncrementsAnswer ) ){
        results << "test_XDMFDataFile_getNumIncrements (test 1) & False\n";
        return 1;
    }

    results << "test_XDMFDataFile_getNumIncrements & True\n";
    return 0;
}

int test_XDMFDataFile_getDomainNodes( std::ofstream &results ){
    /*!
     * Get the nodes from a domain.
     *
     * :param std::ofstream &results The output file
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdf( yf[ "filetest1" ] );

    uIntVector domainNodesAnswer = { 2, 3, 6, 7, 8, 10, 12, 13 };

    uIntVector domainNodesResult;
    std::string domainName = "left";
    errorOut error = xdf.getDomainNodes( 0, domainName, domainNodesResult );

    if ( error ){
        error->print( );
        results << "test_XDMFDataFile_getDomainNodes & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( domainNodesResult, domainNodesAnswer ) ){
        results << "test_XDMFDataFile_getDomainNodes (test 1) & False\n";
        return 1;
    }

    domainName = "free";
    error = xdf.getDomainNodes( 0, domainName, domainNodesResult );

    if ( !error ){
        results << "test_XDMFDataFile_getDomainNodes (test 2) & False\n";
        return 1;
    }

    results << "test_XDMFDataFile_getDomainNodes & True\n";
    return 0;
}

int test_XDMFDataFile_getNumNodes( std::ofstream &results ){
    /*!
     * Test the function to extract the number of nodes in the domain
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdf( yf[ "filetest1" ] );

    unsigned int answer = 16;
    unsigned int result;

    errorOut error = xdf.getNumNodes( 0, result );

    if ( error ){
        error->print( );
        results << "test_XDMFDataFile_getNumNodes & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_XDMFDataFile_getNumNodes (test 1) & False\n";
        return 1;
    }

    results << "test_XDMFDataFile_getNumNodes & True\n";
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
    test_XDMFDataFile_getNumIncrements( results );
    test_XDMFDataFile_readMesh( results );
    test_XDMFDataFile_getDomainNodes( results );
    test_XDMFDataFile_getNumNodes( results );

    //Close the results file
    results.close();
}
