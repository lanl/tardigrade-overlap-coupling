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
typedef dataFileInterface::uIntType uIntType; //!Define the unsigned int type
typedef dataFileInterface::uIntVector uIntVector; //!Define a vector of unsigned ints
typedef dataFileInterface::stringVector stringVector; //!Define a vector of strings

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

    if ( df->_filename.compare( "../testFiles/macroscale_xdmf.xdmf" ) != 0 ){
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

    if ( df->_filename.compare( "../testFiles/macroscale_xdmf.xdmf" ) != 0 ){
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

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    df = dataFileInterface::dataFileBase( yf[ "filetest3" ] ).create( "XDMF" );

    if ( df->_error ){
        results << "test_XDMFDataFile_constructor & False\n";
        return 1;
    }

    std::ifstream infile( "test_output.xdmf" );
    if ( !infile.good( ) ){
        results << "test_XDMFDataFile_constructor (test 7) & False\n";
        return 1;
    }

    infile = std::ifstream( "test_output.h5" );
    if ( !infile.good( ) ){
        results << "test_XDMFDataFile_constructor (test 8) & False\n";
        return 1;
    }

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

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
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

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

    errorOut error = xdmf.readMesh( 1, nodePositionsResult );

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
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    unsigned int numIncrementsAnswer = 2;
    unsigned int numIncrementsResult;

    errorOut error = xdmf.getNumIncrements( numIncrementsResult );
    
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

int test_XDMFDataFile_getSubDomainNodes( std::ofstream &results ){
    /*!
     * Get the nodes from a domain.
     *
     * :param std::ofstream &results The output file
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    uIntVector domainNodesAnswer = { 2, 3, 6, 7, 8, 10, 12, 13 };

    uIntVector domainNodesResult;
    std::string domainName = "left";
    std::unique_ptr< errorNode > error;
    error.reset( xdmf.getSubDomainNodes( 0, domainName, domainNodesResult ) );

    if ( error ){
        error->print( );
        results << "test_XDMFDataFile_getSubDomainNodes & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( domainNodesResult, domainNodesAnswer ) ){
        results << "test_XDMFDataFile_getSubDomainNodes (test 1) & False\n";
        return 1;
    }

    domainName = "free";
    error.reset( xdmf.getSubDomainNodes( 0, domainName, domainNodesResult ) );

    if ( !error ){
        results << "test_XDMFDataFile_getSubDomainNodes (test 2) & False\n";
        return 1;
    }

    results << "test_XDMFDataFile_getSubDomainNodes & True\n";
    return 0;
}

int test_XDMFDataFile_getNumNodes( std::ofstream &results ){
    /*!
     * Test the function to extract the number of nodes in the domain
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    unsigned int answer = 16;
    unsigned int result;

    std::unique_ptr< errorNode > error;
    error.reset( xdmf.getNumNodes( 0, result ) );

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

int test_XDMFDataFile_getSetNames( std::ofstream &results ){
    /*!
     * Test the function to extract the names of the sets
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    std::vector< std::string > answer = { "free_nodes", "ghost_nodes",
                                          "left", "right", "bottom", "top", "back", "front", "all",
                                          "non_overlapped_nodes", "non_overlapped_elements",
                                          "free_elements", "ghost_elements" };
    std::vector< std::string > result;

    errorOut error = xdmf.getSetNames( 1, result );

    if ( error ){
        error->print( );
        results << "test_XDMFDataFile_getSetNames & False\n";
        return 1;
    }

    if ( answer.size() != result.size() ){
        results << "test_XDMFDataFile_getSetNames (test 1) & False\n";
        return 1;
    }

    for ( unsigned int i = 0; i < result.size( ); i++ ){

        if ( answer[ i ].compare( result[ i ] ) != 0 ){

            results << "test_XDMFDataFile_getSetNames (test 2) & False\n";
            return 1;

        }

    }

    results << "test_XDMFDataFile_getSetNames & True\n";
    return 0;
}

int test_XDMFDataFile_getSolutionData( std::ofstream &results ){
    /*!
     * Test the function to extract the solution data
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    floatVector answer = { -0.001, -0.001, -0.001,
                           -0.001, -0.001, -0.001,
                           -0.001, -0.001, -0.001,
                           -0.001, -0.001, -0.001,
                           -0.001, -0.001, -0.001,
                           -0.001 };

    floatVector result;
    errorOut error = xdmf.getSolutionData( 1, "disp_z", "Node", result );

    if ( error ){
        error->print( );
        results << "test_XDMFDataFile_getSolutionData & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_XDMFDataFile_getSolutionData (test 1) & False\n";
        return 1;
    }

    results << "test_XDMFDataFile_getSolutionData & True\n";
    return 0;

}

int test_XDMFDataFile_getMeshData( std::ofstream &results ){
    /*!
     * Test the function to get the mesh data
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "testConfig_polyhedron.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    floatVector nodePositionAnswer = { 1, 0, 1, 1, 0, 0, 0, 0, 0,
                                       0, 0, 1, 1, 1, 1, 1, 1, 0,
                                       0, 1, 0, 0, 1, 1, 0, 1, 2,
                                       1, 1, 2, 0, 0, 2, 1, 0, 2,
                                       0, 0, 3, 0, 1, 3, 1, 1, 3,
                                       1, 0, 3 };

    uIntVector connectivityAnswer = { 9,  0,  1,  2,  3,  4, 5, 6,  7,
                                      9,  8,  7,  4,  9, 10, 3, 0, 11,
                                      9, 12, 13, 14, 15, 10, 8, 9, 11 };

    uIntVector connectivityCellIndicesAnswer = { 0, 9, 18 };

    unsigned int cellCountAnswer = 3;

    floatVector nodePositionResult;
    uIntVector connectivityResult, connectivityCellIndicesResult;
    unsigned int cellCountResult;

    errorOut error = xdmf.getMeshData( 1, nodePositionResult, connectivityResult, connectivityCellIndicesResult, cellCountResult );

    if ( error ){
        error->print( );
        results << "test_XDMFDataFile_getMeshData & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( nodePositionAnswer, nodePositionResult ) ){
        results << "test_XDMFDataFile_getMeshData (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( connectivityAnswer, connectivityResult ) ){
        results << "test_XDMFDataFile_getMeshData (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( cellCountAnswer, cellCountResult ) ){
        results << "test_XDMFDataFile_getMeshData (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( connectivityCellIndicesAnswer, connectivityCellIndicesResult ) ){
        results << "test_XDMFDataFile_getMeshData (test 4) & False\n";
        return 1;
    }

    results << "test_XDMFDataFile_getMeshData & True\n";
    return 0;
}

int test_XDMFDataFile_getMeshData2( std::ofstream &results ){
    /*!
     * Second test the function to get the mesh data
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    floatVector nodePositionAnswer = { 1, 0, 1, 1, 0, 0, 0, 0, 0,
                                       0, 0, 1, 1, 1, 1, 1, 1, 0,
                                       0, 1, 0, 0, 1, 1, 0, 1, 2,
                                       1, 1, 2, 0, 0, 2, 1, 0, 2,
                                       0, 0, 3, 0, 1, 3, 1, 1, 3,
                                       1, 0, 3 };

    uIntVector connectivityAnswer = { 9,  0,  1,  2,  3,  4, 5, 6,  7,
                                      9,  8,  7,  4,  9, 10, 3, 0, 11,
                                      9, 12, 13, 14, 15, 10, 8, 9, 11 };

    uIntVector connectivityCellIndicesAnswer = { 0, 9, 18 };

    unsigned int cellCountAnswer = 3;

    floatVector nodePositionResult;
    uIntVector connectivityResult, connectivityCellIndicesResult;
    unsigned int cellCountResult;

    errorOut error = xdmf.getMeshData( 1, nodePositionResult, connectivityResult, connectivityCellIndicesResult, cellCountResult );

    if ( error ){
        error->print( );
        results << "test_XDMFDataFile_getMeshData2 & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( nodePositionAnswer, nodePositionResult ) ){
        results << "test_XDMFDataFile_getMeshData2 (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( connectivityAnswer, connectivityResult ) ){
        results << "test_XDMFDataFile_getMeshData2 (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( cellCountAnswer, cellCountResult ) ){
        results << "test_XDMFDataFile_getMeshData2 (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( connectivityCellIndicesAnswer, connectivityCellIndicesResult ) ){
        results << "test_XDMFDataFile_getMeshData2 (test 4) & False\n";
        return 1;
    }

    results << "test_XDMFDataFile_getMeshData2 & True\n";
    return 0;
}

int test_XDMFDataFile_getNumSubDomainNodes( std::ofstream &results ){
    /*!
     * Test the determination of the number of nodes are in a given domain
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    unsigned int numSubDomainNodesAnswer = 8;

    unsigned int numSubDomainNodesResult;
    std::string domainName = "left";
    std::unique_ptr< errorNode > error;
    error.reset( xdmf.getNumSubDomainNodes( 0, domainName, numSubDomainNodesResult ) );

    if ( error ){
        error->print( );
        results << "test_XDMFDataFile_getNumSubDomainNodes & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( numSubDomainNodesResult, numSubDomainNodesAnswer ) ){
        results << "test_XDMFDataFile_getNumSubDomainNodes (test 1) & False\n";
        return 1;
    }

    domainName = "free";
    error.reset( xdmf.getNumSubDomainNodes( 0, domainName, numSubDomainNodesResult ) );

    if ( !error ){
        results << "test_XDMFDataFile_getNumSubDomainNodes (test 2) & False\n";
        return 1;
    }

    results << "test_XDMFDataFile_getNumSubDomainNodes & True\n";
    return 0;
}

int test_XDMFDataFile_getSolutionVectorDataFromComponents( std::ofstream &results ){
    /*!
     * Test the extraction of vector solution data from a file where the 
     * components are separated.
     *
     * :param std::ofstream &results: The output file.
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    floatVector answer = { 0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001,
                           0., 0., -0.001 };

    floatVector result;
    stringVector componentNames = { "disp_x", "disp_y", "disp_z" };
    errorOut error = xdmf.getSolutionVectorDataFromComponents( 1, componentNames, "Node", result );

    if ( error ){
        error->print( );
        results << "test_XDMFDataFile_getSolutionVectorDataFromComponents & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_XDMFDataFile_getSolutionVectorDataFromComponents (test 1) & False\n";
        return 1;
    }


    results << "test_XDMFDataFile_getSolutionVectorDataFromComponents & True\n";
    return 0;
}

int test_XDMFDataFile_getIncrementTime( std::ofstream &results ){
    /*!
     * Test the extraction of the timestamp for a given increment
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    floatType answer1 = 0.;
    floatType result;
    errorOut error = xdmf.getIncrementTime( 0, result );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_getIncrementTime & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( result, answer1 ) ){

        results << "test_XDMFDataFile_getIncrementTime (test 1) & False\n";
        return 1;

    }

    floatType answer2 = 1;
    error = xdmf.getIncrementTime( 1, result );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_getIncrementTime & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( result, answer2 ) ){

        results << "test_XDMFDataFile_getIncrementTime (test 1) & False\n";
        return 1;

    }

    results << "test_XDMFDataFile_getIncrementTime & True\n";
    return 0;
}

int test_XDMFDataFile_writeIncrementMeshData( std::ofstream &results ){
    /*!
     * Write the increment mesh data
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest3" ] );

    floatType timeAnswer = 0.0;

    uIntType reference_increment = 0;

    uIntVector nodeIdsAnswer = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

    floatVector nodePositionsAnswer = { 1, 0, 1, 1, 0, 0, 0, 0, 0,
                                        0, 0, 1, 1, 1, 1, 1, 1, 0,
                                        0, 1, 0, 0, 1, 1, 0, 1, 2,
                                        1, 1, 2, 0, 0, 2, 1, 0, 2,
                                        0, 0, 3, 0, 1, 3, 1, 1, 3,
                                        1, 0, 3 };

    uIntVector elementIdsAnswer = { 1, 2, 3 };

    uIntVector connectivityAnswer = { 16, 6,
                                      4, 0, 3, 2, 1,
                                      4, 0, 1, 5, 4,
                                      4, 1, 2, 6, 5,
                                      4, 2, 3, 7, 6,
                                      4, 3, 0, 4, 7,
                                      4, 4, 5, 6, 7,
                                      16, 6,
                                      4, 8, 9, 4, 7,
                                      4, 8, 7, 3, 10,
                                      4, 7, 4, 0, 3,
                                      4, 4, 9, 11, 0,
                                      4, 9, 8, 10, 11,
                                      4, 10, 3, 0, 11,
                                      16, 6,
                                      4, 12, 15, 14, 13,
                                      4, 12, 13, 8, 10,
                                      4, 13, 14, 9, 8,
                                      4, 14, 15, 11, 9,
                                      4, 15, 12, 10, 11,
                                      4, 10, 8, 9, 11 };

    uIntVector cellIndicesAnswer = { 0, 32, 64 };
    uIntType cellCountsAnswer = 3;

    uIntType increment;
    uIntType collectionNumber = 0;

    errorOut error = xdmf.initializeIncrement( timeAnswer, reference_increment, collectionNumber, increment );

    if ( error ){

        error->print( );
        results << "test_writeIncrementMeshData & False\n";
        return 1;

    }

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    error = xdmf.writeIncrementMeshData( increment, collectionNumber, nodeIdsAnswer, { { } }, { { } }, nodePositionsAnswer,
                                        elementIdsAnswer, { { } }, { { } }, connectivityAnswer );

    if ( error ){

        error->print( );
        results << "test_writeIncrementMeshData & False\n";
        return 1;

    }

    //Read in the mesh data to determine if things have been stored correctly
    YAML::Node af = YAML::Load( "mode: read\nfilename: test_output.xdmf\ncell_id_variable_name: ELEMID\n" );
    dataFileInterface::XDMFDataFile xdmf_result( af );

    if ( xdmf_result._error ){

        xdmf_result._error->print( );
        results << "test_writeIncrementMeshData & False\n";
        return 1;

    }

    //Check if the timestep was stored correctly
    floatType scalarResult;
    error = xdmf_result.getIncrementTime( increment, scalarResult );

    if ( error ){
        error->print( );
        results << "test_writeIncrementMeshData & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( scalarResult, 0.0 ) ){
        results << "test_writeIncrementMeshData (test 1) & False\n";
        return 1;
    }

    //Check if the mesh information is stored correctly
    uIntVector nodeIdsResult;
    error = xdmf_result.getNodeIds( increment, "NODEID", nodeIdsResult );

    if ( error ){

        error->print( );
        results << "test_writeIncrementMeshData & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( nodeIdsResult, nodeIdsAnswer ) ){

        results << "test_writeIncrementMeshData (test 2) & False\n";
        return 1;

    }

    floatVector nodePositionsResult;
    uIntVector connectivityResult;
    uIntVector cellIndicesResult;
    uIntType cellCountsResult;

    error = xdmf_result.getMeshData( increment, nodePositionsResult, connectivityResult, cellIndicesResult, cellCountsResult );

    if ( error ){

        error->print( );
        results << "test_writeIncrementMeshData & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( nodePositionsAnswer, nodePositionsResult ) ){

        results << "test_writeIncrementMeshData (test 3) & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( connectivityAnswer, connectivityResult ) ){

        results << "test_writeIncrementMeshData (test 4) & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( cellIndicesAnswer, cellIndicesResult ) ){

        vectorTools::print( cellIndicesAnswer );
        vectorTools::print( cellIndicesResult );
        results << "test_writeIncrementMeshData (test 5) & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( cellCountsAnswer, cellCountsResult ) ){

        results << "test_writeIncrementMeshData (test 6) & False\n";
        return 1;

    }


    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    results << "test_writeIncrementMeshData & True\n";
    return 0;

}

int test_XDMFDataFile_getNodeIds( std::ofstream &results ){
    /*!
     * Extract the node ids from the domain
     */


    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    if ( xdmf._error ){

        xdmf._error->print( );
        results << "test_XDMFDataFile_getNodeIds & False\n";
        return 1;
    }

    uIntVector nodeIdAnswer = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

    uIntVector nodeIdResult;
    errorOut error = xdmf.getNodeIds( 0, "NODEID", nodeIdResult );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_getNodeIds & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( nodeIdResult, nodeIdAnswer ) ){

        vectorTools::print( nodeIdResult );
        vectorTools::print( nodeIdAnswer );
        results << "test_XDMFDataFile_getNodeIds (test 1) & False\n";
        return 1;

    }

    results << "test_XDMFDataFile_getNodeIds & True\n";
    return 0;

}

int test_XDMFDataFile_getCellIds( std::ofstream &results ){
    /*!
     * Extract the cell ids from the domain
     */


    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    if ( xdmf._error ){

        xdmf._error->print( );
        results << "test_XDMFDataFile_getCellIds & False\n";
        return 1;
    }

    uIntVector cellIdAnswer = { 0, 1, 2 };

    uIntVector cellIdResult;
    errorOut error = xdmf.getCellIds( 0, "ELEMID", cellIdResult );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_getCellIds & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( cellIdResult, cellIdAnswer ) ){

        vectorTools::print( cellIdResult );
        vectorTools::print( cellIdAnswer );
        results << "test_XDMFDataFile_getCellIds (test 1) & False\n";
        return 1;

    }

    results << "test_XDMFDataFile_getCellIds & True\n";
    return 0;

}

int test_XDMFDataFile_initializeIncrement( std::ofstream &results ){
    /*!
     * Initialize an increment in an output XDMF data file
     *
     * :param std::ofstream &results: The output file
     */

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest3" ] );

    if ( xdmf._error ){

        xdmf._error->print( );
        results << "test_XDMFDataFile_initializeIncrement & False\n";
        return 1;
    }

    uIntType incrementAnswer1 = 0;
    uIntType incrementResult;

    errorOut error = xdmf.initializeIncrement( 0.0, 0, 0, incrementResult );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_initializeIncrement & False\n";
        return 1;

    }

    if ( incrementResult != incrementAnswer1 ){

        results << "test_XDMFDataFile_initializeIncrement (test 1) & False\n";
        return 1;

    }

    uIntType incrementAnswer2 = 1;

    error = xdmf.initializeIncrement( 0.1, 0, 0, incrementResult );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_initializeIncrement & False\n";
        return 1;

    }

    if ( incrementResult != incrementAnswer2 ){

        results << "test_XDMFDataFile_initializeIncrement (test 2) & False\n";
        return 1;

    }

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    results << "test_XDMFDataFile_initializeIncrement & True\n";
    return 0;
}

int test_XDMFDataFile_addRootCollection( std::ofstream &results ){
    /*!
     * Test adding a root collection to the datafile
     *
     * :param std::ofstream &results: The output file
     */

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest3" ] );

    if ( xdmf._error ){

        xdmf._error->print( );
        results << "test_XDMFDataFile_addRootCollection & False\n";
        return 1;
    }

    uIntType collectionNumberResult;
    errorOut error = xdmf.addRootCollection( "TEST", "Test collection info", collectionNumberResult );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_addRootCollection & False\n";
        return 1;

    }

    if ( collectionNumberResult != 1 ){

        results << "test_XDMFDataFile_addRootCollection (test 1) & False\n";
        return 1;

    }

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    results << "test_XDMFDataFile_addRootCollection & True\n";
    return 0;
}

int test_XDMFDataFile_writeScalarSolutionData( std::ofstream &results ){
    /*!
     * Test writing a scalar solution data to the output file
     *
     * :param std::ofstream &results: The output file
     */

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest3" ] );

    if ( xdmf._error ){

        xdmf._error->print( );
        results << "test_XDMFDataFile_writeScalarSolutionData & False\n";
        return 1;
    }

    uIntType increment;
    uIntType collectionNumber = 0;
    errorOut error = xdmf.initializeIncrement( 0.0, collectionNumber, collectionNumber, increment );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeScalarSolutionData & False\n";
        return 1;

    }

    uIntVector nodeIds = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

    floatVector nodePositions = { 1, 0, 1, 1, 0, 0, 0, 0, 0,
                                  0, 0, 1, 1, 1, 1, 1, 1, 0,
                                  0, 1, 0, 0, 1, 1, 0, 1, 2,
                                  1, 1, 2, 0, 0, 2, 1, 0, 2,
                                  0, 0, 3, 0, 1, 3, 1, 1, 3,
                                  1, 0, 3 };

    uIntVector elementIds = { 1, 2, 3 };

    uIntVector connectivity = { 16, 6,
                                4, 0, 3, 2, 1,
                                4, 0, 1, 5, 4,
                                4, 1, 2, 6, 5,
                                4, 2, 3, 7, 6,
                                4, 3, 0, 4, 7,
                                4, 4, 5, 6, 7,
                                16, 6,
                                4, 8, 9, 4, 7,
                                4, 8, 7, 3, 10,
                                4, 7, 4, 0, 3,
                                4, 4, 9, 11, 0,
                                4, 9, 8, 10, 11,
                                4, 10, 3, 0, 11,
                                16, 6,
                                4, 12, 15, 14, 13,
                                4, 12, 13, 8, 10,
                                4, 13, 14, 9, 8,
                                4, 14, 15, 11, 9,
                                4, 15, 12, 10, 11,
                                4, 10, 8, 9, 11 };

    error = xdmf.writeIncrementMeshData( increment, collectionNumber, nodeIds, { { } }, { { } }, nodePositions,
                                         elementIds, { { } }, { { } }, connectivity );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeScalarSolutionData & False\n";
        return 1;

    }

    floatVector nodeDataAnswer = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5 };

    error = xdmf.writeScalarSolutionData( increment, 0, "TEST_DATA", "NODE", nodeDataAnswer );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeScalarSolutionData & False\n";
        return 1;

    }

    floatVector bigNodeDataAnswer( 1000, 1 );

    error = xdmf.writeScalarSolutionData( increment, 0, "BIG_TEST_DATA", "NODE", bigNodeDataAnswer );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeScalarSolutionData & False\n";
        return 1;

    }

    floatVector elementDataAnswer = { -1, -2, -3 };

    error = xdmf.writeScalarSolutionData( increment, 0, "TEST_DATA_", "CeLl", elementDataAnswer );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeScalarSolutionData & False\n";
        return 1;

    }

    //Read in the mesh data to determine if things have been stored correctly
    YAML::Node af = YAML::Load( "mode: read\nfilename: test_output.xdmf\ncell_id_variable_name: ELEMID\n" );
    dataFileInterface::XDMFDataFile xdmf_result( af );

    if ( xdmf_result._error ){

        xdmf_result._error->print( );
        results << "test_writeSolutionData & False\n";
        return 1;

    }

    floatVector nodeDataResult;
    error = xdmf_result.getSolutionData( increment, "TEST_DATA", "Node", nodeDataResult );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeScalarSolutionData & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( nodeDataResult, nodeDataAnswer ) ){

        results << "test_XDMFDataFile_writeScalarSolutionData (test 1) & False\n";
        return 1;

    }

    floatVector elementDataResult;
    error = xdmf_result.getSolutionData( increment, "TEST_DATA_", "Cell", elementDataResult );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeScalarSolutionData & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( elementDataResult, elementDataAnswer ) ){

        results << "test_XDMFDataFile_writeScalarSolutionData (test 2) & False\n";
        return 1;

    }

    floatVector bigNodeDataResult;
    error = xdmf_result.getSolutionData( increment, "BIG_TEST_DATA", "Node", bigNodeDataResult );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeScalarSolutionData (test 3) & False\n";
        return 1;

    }

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    results << "test_XDMFDataFile_writeScalarSolutionData & True\n";
    return 0;
}

int test_XDMFDataFile_writeSolutionData( std::ofstream &results ){
    /*!
     * Test writing the solution data to the XDMF output file.
     *
     * :param std::ofstream &results: The output file
     */

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest3" ] );

    if ( xdmf._error ){

        xdmf._error->print( );
        results << "test_XDMFDataFile_writeSolutionData & False\n";
        return 1;
    }

    uIntType increment;
    uIntType collectionNumber = 0;
    errorOut error = xdmf.initializeIncrement( 0.0, 0, collectionNumber, increment );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeScalarSolutionData & False\n";
        return 1;

    }

    uIntVector nodeIds = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

    floatVector nodePositions = { 1, 0, 1, 1, 0, 0, 0, 0, 0,
                                  0, 0, 1, 1, 1, 1, 1, 1, 0,
                                  0, 1, 0, 0, 1, 1, 0, 1, 2,
                                  1, 1, 2, 0, 0, 2, 1, 0, 2,
                                  0, 0, 3, 0, 1, 3, 1, 1, 3,
                                  1, 0, 3 };

    uIntVector elementIds = { 1, 2, 3 };

    uIntVector connectivity = { 16, 6,
                                4, 0, 3, 2, 1,
                                4, 0, 1, 5, 4,
                                4, 1, 2, 6, 5,
                                4, 2, 3, 7, 6,
                                4, 3, 0, 4, 7,
                                4, 4, 5, 6, 7,
                                16, 6,
                                4, 8, 9, 4, 7,
                                4, 8, 7, 3, 10,
                                4, 7, 4, 0, 3,
                                4, 4, 9, 11, 0,
                                4, 9, 8, 10, 11,
                                4, 10, 3, 0, 11,
                                16, 6,
                                4, 12, 15, 14, 13,
                                4, 12, 13, 8, 10,
                                4, 13, 14, 9, 8,
                                4, 14, 15, 11, 9,
                                4, 15, 12, 10, 11,
                                4, 10, 8, 9, 11 };

    error = xdmf.writeIncrementMeshData( increment, collectionNumber, nodeIds, { { } }, { { } }, nodePositions,
                                         elementIds, { { } }, { { } }, connectivity );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeSolutionData & False\n";
        return 1;

    }

    floatVector nodeDataAnswer =
        {
            0.0, 0.1, 0.2,
            0.3, 0.4, 0.5,
            0.6, 0.7, 0.8,
            0.9, 1.0, 1.1,
            1.2, 1.3, 1.4,
            1.5, 1.6, 1.7,
            1.8, 1.9, 2.0,
            2.1, 2.2, 2.3,
            2.4, 2.5, 2.6,
            2.7, 2.8, 2.9,
            3.0, 3.1, 3.2,
            3.3, 3.4, 3.5,
            3.6, 3.7, 3.8,
            3.9, 4.0, 4.1,
            4.2, 4.3, 4.4
        };

    error = xdmf.writeSolutionData( increment, 0, { "TEST_DATA_1", "TEST_DATA_2", "TEST_DATA_3" }, "NODE", nodeDataAnswer );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeSolutionData & False\n";
        return 1;

    }

    floatVector elementDataAnswer = { -1, -2, -3,
                                      -4, -5, -6 };

    error = xdmf.writeSolutionData( increment, 0, { "TEST_DATA_1_", "TEST_DATA_2_" }, "CeLl", elementDataAnswer );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeSolutionData & False\n";
        return 1;

    }

    //Read in the mesh data to determine if things have been stored correctly
    YAML::Node af = YAML::Load( "mode: read\nfilename: test_output.xdmf\ncell_id_variable_name: ELEMID\n" );
    dataFileInterface::XDMFDataFile xdmf_result( af );

    if ( xdmf_result._error ){

        xdmf_result._error->print( );
        results << "test_writeSolutionData & False\n";
        return 1;

    }

    for ( uIntType i = 0; i < 3; i++ ){

        floatVector nodeDataResult;
        std::string name = "TEST_DATA_" + std::to_string( i + 1 );
        error = xdmf_result.getSolutionData( increment, name, "Node", nodeDataResult );
    
        if ( error ){
    
            error->print( );
            results << "test_XDMFDataFile_writeSolutionData & False\n";
            return 1;
    
        }
    
        uIntType indx = 0;
        for ( unsigned int j = i; j < nodeDataAnswer.size( ); j += 3, indx++ ){

            if ( !vectorTools::fuzzyEquals( nodeDataResult[ indx ], nodeDataAnswer[ j ] ) ){
    
                results << "test_XDMFDataFile_writeSolutionData (test 1) & False\n";
                return 1;
    
            }

        }

    }

    for ( uIntType i = 0; i < 2; i++ ){

        floatVector elementDataResult;
        std::string name = "TEST_DATA_" + std::to_string( i + 1 ) + "_";
        error = xdmf_result.getSolutionData( increment, name, "Cell", elementDataResult );
    
        if ( error ){
    
            error->print( );
            results << "test_XDMFDataFile_writeSolutionData & False\n";
            return 1;
    
        }

        uIntType indx = 0;
        for ( unsigned int j = i; j < elementDataAnswer.size( ); j += 2, indx++ ){
    
            if ( !vectorTools::fuzzyEquals( elementDataResult[ indx ], elementDataAnswer[ j ] ) ){
      
                results << "test_XDMFDataFile_writeSolutionData (test 2) & False\n";
                return 1;
        
            }

        }

    }

    dataFileInterface::XDMFDataFile xdmf2( yf[ "filetest3" ] );

    error = xdmf2.initializeIncrement( 1.0, 0, collectionNumber, increment );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeSolutionData & False\n";
        return 1;

    }

    if ( increment != 1 ){

        results << "test_XDMFDataFile_writeSolutionData (test 3) & False\n";
        return 1;

    }

    error = xdmf2.writeIncrementMeshData( increment, collectionNumber, { }, { { } }, { }, { }, { }, { { } }, { }, { } );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeSolutionData & False\n";
        return 1;

    }

    floatVector nodeDataAnswer2 = nodeDataAnswer + 1.;

    error = xdmf2.writeSolutionData( increment, 0, { "TEST_DATA_1", "TEST_DATA_2", "TEST_DATA_3" }, "NODE", nodeDataAnswer2 );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeSolutionData & False\n";
        return 1;

    }

    floatVector elementDataAnswer2 = elementDataAnswer - 2.;

    error = xdmf2.writeSolutionData( increment, 0, { "TEST_DATA_1_", "TEST_DATA_2_" }, "CeLl", elementDataAnswer2 );

    if ( error ){

        error->print( );
        results << "test_XDMFDataFile_writeSolutionData & False\n";
        return 1;

    }

    dataFileInterface::XDMFDataFile xdmf_result2( af );
    
    if ( xdmf_result2._error ){

        xdmf_result2._error->print( );
        results << "test_writeSolutionData & False\n";
        return 1;

    }

    for ( uIntType i = 0; i < 3; i++ ){

        floatVector nodeDataResult;
        std::string name = "TEST_DATA_" + std::to_string( i + 1 );
        error = xdmf_result2.getSolutionData( increment, name, "Node", nodeDataResult );
    
        if ( error ){
    
            error->print( );
            results << "test_XDMFDataFile_writeSolutionData & False\n";
            return 1;
    
        }
    
        uIntType indx = 0;
        for ( unsigned int j = i; j < nodeDataAnswer.size( ); j += 3, indx++ ){

            if ( !vectorTools::fuzzyEquals( nodeDataResult[ indx ], nodeDataAnswer2[ j ] ) ){
    
                results << "test_XDMFDataFile_writeSolutionData (test 4) & False\n";
                return 1;
    
            }

        }

    }

    for ( uIntType i = 0; i < 2; i++ ){

        floatVector elementDataResult;
        std::string name = "TEST_DATA_" + std::to_string( i + 1 ) + "_";
        error = xdmf_result2.getSolutionData( increment, name, "Cell", elementDataResult );
    
        if ( error ){
    
            error->print( );
            results << "test_XDMFDataFile_writeSolutionData & False\n";
            return 1;
    
        }

        uIntType indx = 0;
        for ( unsigned int j = i; j < elementDataAnswer.size( ); j += 2, indx++ ){
    
            if ( !vectorTools::fuzzyEquals( elementDataResult[ indx ], elementDataAnswer2[ j ] ) ){
      
                results << "test_XDMFDataFile_writeSolutionData (test 5) & False\n";
                return 1;
        
            }

        }

    }

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    results << "test_XDMFDataFile_writeSolutionData & True\n";
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
    test_XDMFDataFile_getNumSubDomainNodes( results );
    test_XDMFDataFile_getNodeIds( results );
    test_XDMFDataFile_getCellIds( results );
    test_XDMFDataFile_getSubDomainNodes( results );
    test_XDMFDataFile_getNumNodes( results );
    test_XDMFDataFile_getSetNames( results );
    test_XDMFDataFile_getSolutionData( results );
    test_XDMFDataFile_getSolutionVectorDataFromComponents( results );
    test_XDMFDataFile_getMeshData( results );
    test_XDMFDataFile_getMeshData2( results );
    test_XDMFDataFile_getIncrementTime( results );

    test_XDMFDataFile_initializeIncrement( results );
    test_XDMFDataFile_addRootCollection( results );
    test_XDMFDataFile_writeIncrementMeshData( results );
    test_XDMFDataFile_writeScalarSolutionData( results );
    test_XDMFDataFile_writeSolutionData( results );

    //Close the results file
    results.close();
}
