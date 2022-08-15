//!The test file for dataFileInterface.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>
#define USE_EIGEN
#include<vector_tools.h>

#include<dataFileInterface.h>

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

BOOST_AUTO_TEST_CASE( testXDMFDataFile_constructor ){
    /*!
     * Test the interface with the XDMF file format
     * constructor
     *
     */

    std::shared_ptr<dataFileInterface::dataFileBase> df;
    df = dataFileInterface::dataFileBase().create( "XDMF" );

    BOOST_CHECK( df->_error );

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    df = dataFileInterface::dataFileBase( yf["filetest1"] ).create( "XDMF" );

    BOOST_CHECK( !df->_error );

    BOOST_CHECK( df->_filename.compare( "testFiles/macroscale_xdmf.xdmf" ) == 0 );

    BOOST_CHECK( df->_mode.compare( "read" ) == 0 );

    df = dataFileInterface::dataFileBase( yf["filetest1"] ).create( );

    BOOST_CHECK( df );

    BOOST_CHECK( !df->_error );

    BOOST_CHECK( df->_filename.compare( "testFiles/macroscale_xdmf.xdmf" ) == 0 );

    BOOST_CHECK( df->_mode.compare( "read" ) == 0 );

    df = dataFileInterface::dataFileBase( yf[ "filetest2" ] ).create( "XDMF" );

    BOOST_CHECK( df->_error );

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    df = dataFileInterface::dataFileBase( yf[ "filetest3" ] ).create( "XDMF" );

    BOOST_CHECK( !df->_error );

    std::ifstream infile( "test_output.xdmf" );
    BOOST_CHECK( infile.good( ) );

    infile = std::ifstream( "test_output.h5" );
    BOOST_CHECK( infile.good( ) );

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_readMesh ){
    /*!
     * Test the interface with the mesh for the XDMF file format.
     *
     */

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
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

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( nodePositionsAnswer, nodePositionsResult ) );

}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_getNumIncrements ){
    /*!
     * Test the interface with the XDMF file to get the number of
     * temporal increments.
     *
     */

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    unsigned int numIncrementsAnswer = 2;
    unsigned int numIncrementsResult;

    errorOut error = xdmf.getNumIncrements( numIncrementsResult );
    
    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( numIncrementsResult, numIncrementsAnswer ) );

}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_getSubDomainNodes ){
    /*!
     * Get the nodes from a domain.
     *
     */

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    uIntVector domainNodesAnswer = { 2, 3, 6, 7, 8, 10, 12, 13 };

    uIntVector domainNodesResult;
    std::string domainName = "left";
    std::unique_ptr< errorNode > error;
    error.reset( xdmf.getSubDomainNodes( 0, domainName, domainNodesResult ) );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( domainNodesResult, domainNodesAnswer ) );

    domainName = "free";
    error.reset( xdmf.getSubDomainNodes( 0, domainName, domainNodesResult ) );

    BOOST_CHECK( error );

}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_getNumNodes ){
    /*!
     * Test the function to extract the number of nodes in the domain
     *
     */

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    unsigned int answer = 16;
    unsigned int result;

    std::unique_ptr< errorNode > error;
    error.reset( xdmf.getNumNodes( 0, result ) );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer, result ) );

}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_getSetNames ){
    /*!
     * Test the function to extract the names of the sets
     *
     */

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    std::vector< std::string > answer = { "free_nodes", "ghost_nodes",
                                          "left", "right", "bottom", "top", "back", "front", "all",
                                          "non_overlapped_nodes", "non_overlapped_elements",
                                          "free_elements", "ghost_elements" };
    std::vector< std::string > result;

    errorOut error = xdmf.getSetNames( 1, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( answer.size() == result.size() );

    for ( unsigned int i = 0; i < result.size( ); i++ ){

        BOOST_CHECK( answer[ i ].compare( result[ i ] ) == 0 );

    }

}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_getSolutionData ){
    /*!
     * Test the function to extract the solution data
     *
     */

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    floatVector answer = { -0.001, -0.001, -0.001,
                           -0.001, -0.001, -0.001,
                           -0.001, -0.001, -0.001,
                           -0.001, -0.001, -0.001,
                           -0.001, -0.001, -0.001,
                           -0.001 };

    floatVector result;
    errorOut error = xdmf.getSolutionData( 1, "disp_z", "Node", result );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer, result ) );


}

// NOTE: Commenting this out because we don't actually use this representation
//       of the macro-scale anymore. If we did make a mesh that didn't have
//       an explicit type in XDMF or if we hadn't put it in yet, then this
//       would be the function to test with
//BOOST_AUTO_TEST_CASE( testXDMFDataFile_getMeshData_polyhedron ){
//    /*!
//     * Test the function to get the mesh data
//     *
//     */
//
//    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig_polyhedron.yaml" );
//    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );
//
//    floatVector nodePositionAnswer = { 1, 0, 1, 1, 0, 0, 0, 0, 0,
//                                       0, 0, 1, 1, 1, 1, 1, 1, 0,
//                                       0, 1, 0, 0, 1, 1, 0, 1, 2,
//                                       1, 1, 2, 0, 0, 2, 1, 0, 2,
//                                       0, 0, 3, 0, 1, 3, 1, 1, 3,
//                                       1, 0, 3 };
//
//    uIntVector connectivityAnswer = { 9,  0,  1,  2,  3,  4, 5, 6,  7,
//                                      9,  8,  7,  4,  9, 10, 3, 0, 11,
//                                      9, 12, 13, 14, 15, 10, 8, 9, 11 };
//
//    uIntVector connectivityCellIndicesAnswer = { 0, 9, 18 };
//
//    unsigned int cellCountAnswer = 3;
//
//    floatVector nodePositionResult;
//    uIntVector connectivityResult, connectivityCellIndicesResult;
//    unsigned int cellCountResult;
//
//    errorOut error = xdmf.getMeshData( 1, nodePositionResult, connectivityResult, connectivityCellIndicesResult, cellCountResult );
//
//    BOOST_CHECK( !error );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( nodePositionAnswer, nodePositionResult ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( connectivityAnswer, connectivityResult ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( cellCountAnswer, cellCountResult ) );
//
//    BOOST_CHECK( vectorTools::fuzzyEquals( connectivityCellIndicesAnswer, connectivityCellIndicesResult ) );
//
//}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_getMeshData ){
    /*!
     * Test the function to get the mesh data
     *
     */

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
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

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( nodePositionAnswer, nodePositionResult ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( connectivityAnswer, connectivityResult ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( cellCountAnswer, cellCountResult ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( connectivityCellIndicesAnswer, connectivityCellIndicesResult ) );

}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_getNumSubDomainNodes ){
    /*!
     * Test the determination of the number of nodes are in a given domain
     *
     */

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    unsigned int numSubDomainNodesAnswer = 8;

    unsigned int numSubDomainNodesResult;
    std::string domainName = "left";
    std::unique_ptr< errorNode > error;
    error.reset( xdmf.getNumSubDomainNodes( 0, domainName, numSubDomainNodesResult ) );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( numSubDomainNodesResult, numSubDomainNodesAnswer ) );

    domainName = "free";
    error.reset( xdmf.getNumSubDomainNodes( 0, domainName, numSubDomainNodesResult ) );

    BOOST_CHECK( error );

}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_getSolutionVectorDataFromComponents ){
    /*!
     * Test the extraction of vector solution data from a file where the 
     * components are separated.
     *
     */

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
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

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( answer, result ) );

}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_getIncrementTime ){
    /*!
     * Test the extraction of the timestamp for a given increment
     *
     */

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    floatType answer1 = 0.;
    floatType result;
    errorOut error = xdmf.getIncrementTime( 0, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer1 ) );

    floatType answer2 = 1;
    error = xdmf.getIncrementTime( 1, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer2 ) );

}

int writeIncrementMeshData( const floatType &timeAnswer, const uIntType &reference_increment,
                            const uIntType &collectionNumber, uIntType &increment,
                            const uIntVector &nodeIdsAnswer, const floatVector &nodePositionsAnswer,
                            const uIntVector &elementIdsAnswer, const uIntVector &connectivityAnswer ){
    /*!
     * Write the increment's mesh data for use in the tests
     * 
     */

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest3" ] );

    errorOut error = xdmf.initializeIncrement( timeAnswer, reference_increment, collectionNumber, increment );

    BOOST_CHECK( !error );

    error = xdmf.writeIncrementMeshData( increment, collectionNumber, nodeIdsAnswer, { { } }, { { } }, nodePositionsAnswer,
                                        elementIdsAnswer, { { } }, { { } }, connectivityAnswer );

    BOOST_CHECK( !error );

    return 0;
}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_writeIncrementMeshData ){
    /*!
     * Write the increment mesh data
     *
     */

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

    int result_int = writeIncrementMeshData( timeAnswer, reference_increment, collectionNumber, increment,
                                             nodeIdsAnswer, nodePositionsAnswer, elementIdsAnswer, connectivityAnswer );

    BOOST_CHECK(result_int == 0);

    //Read in the mesh data to determine if things have been stored correctly
    YAML::Node af = YAML::Load( "mode: read\nfilename: test_output.xdmf\ncell_id_variable_name: ELEMID\n" );
    dataFileInterface::XDMFDataFile xdmf_result( af );

    BOOST_CHECK( !xdmf_result._error );

    //Check if the timestep was stored correctly
    floatType scalarResult;
    errorOut error = xdmf_result.getIncrementTime( increment, scalarResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( scalarResult, 0.0 ) );

    //Check if the mesh information is stored correctly
    uIntVector nodeIdsResult;
    error = xdmf_result.getNodeIds( increment, "NODEID", nodeIdsResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( nodeIdsResult, nodeIdsAnswer ) );

    floatVector nodePositionsResult;
    uIntVector connectivityResult;
    uIntVector cellIndicesResult;
    uIntType cellCountsResult;

    error = xdmf_result.getMeshData( increment, nodePositionsResult, connectivityResult, cellIndicesResult, cellCountsResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( nodePositionsAnswer, nodePositionsResult ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( connectivityAnswer, connectivityResult ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( cellIndicesAnswer, cellIndicesResult ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( cellCountsAnswer, cellCountsResult ) );


    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );


}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_getNodeIds ){
    /*!
     * Extract the node ids from the domain
     */


    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    BOOST_CHECK( !xdmf._error );

    uIntVector nodeIdAnswer = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

    uIntVector nodeIdResult;
    errorOut error = xdmf.getNodeIds( 0, "NODEID", nodeIdResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( nodeIdResult, nodeIdAnswer ) );

}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_getCellIds ){
    /*!
     * Extract the cell ids from the domain
     */


    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest1" ] );

    BOOST_CHECK( !xdmf._error );

    uIntVector cellIdAnswer = { 0, 1, 2 };

    uIntVector cellIdResult;
    errorOut error = xdmf.getCellIds( 0, "ELEMID", cellIdResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( cellIdResult, cellIdAnswer ) );

}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_initializeIncrement ){
    /*!
     * Initialize an increment in an output XDMF data file
     *
     */

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest3" ] );

    BOOST_CHECK( !xdmf._error );

    uIntType incrementAnswer1 = 0;
    uIntType incrementResult;

    errorOut error = xdmf.initializeIncrement( 0.0, 0, 0, incrementResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( incrementResult == incrementAnswer1 );

    uIntType incrementAnswer2 = 1;

    error = xdmf.initializeIncrement( 0.1, 0, 0, incrementResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( incrementResult == incrementAnswer2 );

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_addRootCollection ){
    /*!
     * Test adding a root collection to the datafile
     *
     */

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest3" ] );

    BOOST_CHECK( !xdmf._error );

    uIntType collectionNumberResult;
    errorOut error = xdmf.addRootCollection( "TEST", "Test collection info", collectionNumberResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( collectionNumberResult == 1 );

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

}

int writeScalarSolutionData( const uIntType &collectionNumber, uIntType &increment, const uIntVector &nodeIds,
                             const floatVector &nodePositions, const uIntVector &elementIds, const uIntVector &connectivity,
                             const floatVector &nodeDataAnswer, const floatVector &bigNodeDataAnswer, const floatVector &elementDataAnswer ){
    /*!
     * Write the scalar solution data to the datafile for the test
     * 
     * :param const uIntType &collectionNumber: The collection number
     * :param uIntType &increment: The increment
     * :param const uIntVector &nodeIds: The node id numbers
     * :param const floatVector &nodePositions: The nodal positions
     * :param const uIntVector &elementIds: The element id numbers
     * :param const uIntVector &connectivity: The connectivity vector
     * :param const floatVector &nodeDataAnswer: The nodal data vector
     * :param const floatVector &bigNodeDataAnswer: The large nodal data vector
     * :param const floatVector &elementDataAnswer: The element data answer vector
     */

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest3" ] );

    BOOST_CHECK( !xdmf._error );

    errorOut error = xdmf.initializeIncrement( 0.0, collectionNumber, collectionNumber, increment );

    BOOST_CHECK( !error );

    error = xdmf.writeIncrementMeshData( increment, collectionNumber, nodeIds, { { } }, { { } }, nodePositions,
                                         elementIds, { { } }, { { } }, connectivity );

    BOOST_CHECK( !error );

    error = xdmf.writeScalarSolutionData( increment, 0, "TEST_DATA", "NODE", nodeDataAnswer );

    BOOST_CHECK( !error );

    error = xdmf.writeScalarSolutionData( increment, 0, "BIG_TEST_DATA", "NODE", bigNodeDataAnswer );

    BOOST_CHECK( !error );

    error = xdmf.writeScalarSolutionData( increment, 0, "TEST_DATA_", "CeLl", elementDataAnswer );

    BOOST_CHECK( !error );

    return 0;
}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_writeScalarSolutionData ){
    /*!
     * Test writing a scalar solution data to the output file
     *
     */

    uIntType increment;
    uIntType collectionNumber = 0;

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

    floatVector nodeDataAnswer = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5 };

    floatVector bigNodeDataAnswer( 1000, 1 );

    floatVector elementDataAnswer = { -1, -2, -3 };

    int return_int = writeScalarSolutionData( collectionNumber, increment, nodeIds, nodePositions, elementIds,
                                              connectivity, nodeDataAnswer, bigNodeDataAnswer, elementDataAnswer );

    BOOST_CHECK( return_int == 0 );

    //Read in the mesh data to determine if things have been stored correctly
    YAML::Node af = YAML::Load( "mode: read\nfilename: test_output.xdmf\ncell_id_variable_name: ELEMID\n" );
    dataFileInterface::XDMFDataFile xdmf_result( af );

    BOOST_CHECK( !xdmf_result._error );

    floatVector nodeDataResult;
    errorOut error = xdmf_result.getSolutionData( increment, "TEST_DATA", "Node", nodeDataResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( nodeDataResult, nodeDataAnswer ) );

    floatVector elementDataResult;
    error = xdmf_result.getSolutionData( increment, "TEST_DATA_", "Cell", elementDataResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( elementDataResult, elementDataAnswer ) );

    floatVector bigNodeDataResult;
    error = xdmf_result.getSolutionData( increment, "BIG_TEST_DATA", "Node", bigNodeDataResult );

    BOOST_CHECK( !error );

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

}

int writeSolutionData( uIntType &increment, const uIntType &collectionNumber, const uIntVector &nodeIds, const floatVector &nodePositions,
                       const uIntVector &elementIds, const uIntVector &connectivity, const floatVector &nodeDataAnswer,
                       const floatVector &elementDataAnswer ){
    /*!
     * Write the solution data to an output file
     * 
     * :param uIntType &increment: The increment number
     * :param const uIntType &collectionNumber: The collection number
     * :param const uIntVctor &nodeIds: The node ids of the solution
     * :param const floatVector &nodePositions: The nodal positions
     * :param const uIntVector &elementIds: The element ID numbers
     * :param const uIntVector &connectivity: The connectivity vector
     * :param const floatVector &nodeDataAnswer: The node data to be written to the file
     * :param const floatVector &elementDataAnswer: The element data to be written to the file
     */

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    YAML::Node yf = YAML::LoadFile( "dataFileInterface_testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdmf( yf[ "filetest3" ] );

    errorOut error = xdmf.initializeIncrement( 0.0, 0, collectionNumber, increment );

    BOOST_CHECK( !error );

    error = xdmf.writeIncrementMeshData( increment, collectionNumber, nodeIds, { { } }, { { } }, nodePositions,
                                         elementIds, { { } }, { { } }, connectivity );

    BOOST_CHECK( !error );

    error = xdmf.writeSolutionData( increment, 0, { "TEST_DATA_1", "TEST_DATA_2", "TEST_DATA_3" }, "NODE", nodeDataAnswer );

    BOOST_CHECK( !error );

    error = xdmf.writeSolutionData( increment, 0, { "TEST_DATA_1_", "TEST_DATA_2_" }, "CeLl", elementDataAnswer );

    BOOST_CHECK( !error );

    return 0;

}

BOOST_AUTO_TEST_CASE( testXDMFDataFile_writeSolutionData ){
    /*!
     * Test writing the solution data to the XDMF output file.
     *
     */

    uIntType increment;
    uIntType collectionNumber = 0;

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

    floatVector elementDataAnswer = { -1, -2, -3,
                                      -4, -5, -6 };

    int return_value = writeSolutionData( increment, collectionNumber, nodeIds, nodePositions, elementIds,
                                          connectivity, nodeDataAnswer, elementDataAnswer );

    BOOST_CHECK( return_value == 0 );

    //Read in the mesh data to determine if things have been stored correctly
    YAML::Node af = YAML::Load( "mode: read\nfilename: test_output.xdmf\ncell_id_variable_name: ELEMID\n" );
    dataFileInterface::XDMFDataFile xdmf_result( af );

    BOOST_CHECK( !xdmf_result._error );

    errorOut error;

    for ( uIntType i = 0; i < 3; i++ ){

        floatVector nodeDataResult;
        std::string name = "TEST_DATA_" + std::to_string( i + 1 );
        error = xdmf_result.getSolutionData( increment, name, "Node", nodeDataResult );
    
        BOOST_CHECK( !error );
    
        uIntType indx = 0;
        for ( unsigned int j = i; j < nodeDataAnswer.size( ); j += 3, indx++ ){

            BOOST_CHECK( vectorTools::fuzzyEquals( nodeDataResult[ indx ], nodeDataAnswer[ j ] ) );

        }

    }

    for ( uIntType i = 0; i < 2; i++ ){

        floatVector elementDataResult;
        std::string name = "TEST_DATA_" + std::to_string( i + 1 ) + "_";
        error = xdmf_result.getSolutionData( increment, name, "Cell", elementDataResult );
    
        BOOST_CHECK( !error );

        uIntType indx = 0;
        for ( unsigned int j = i; j < elementDataAnswer.size( ); j += 2, indx++ ){
    
            BOOST_CHECK( vectorTools::fuzzyEquals( elementDataResult[ indx ], elementDataAnswer[ j ] ) );

        }

    }

    floatVector nodeDataAnswer2 = nodeDataAnswer + 1.;

    floatVector elementDataAnswer2 = elementDataAnswer - 2.;

    return_value = writeSolutionData( increment, collectionNumber, { }, { }, { }, { }, nodeDataAnswer2,
                                      elementDataAnswer2 );

    BOOST_CHECK( return_value <= 0 );

    dataFileInterface::XDMFDataFile xdmf_result2( af );
    
    BOOST_CHECK( !xdmf_result2._error );

    for ( uIntType i = 0; i < 3; i++ ){

        floatVector nodeDataResult;
        std::string name = "TEST_DATA_" + std::to_string( i + 1 );
        error = xdmf_result2.getSolutionData( increment, name, "Node", nodeDataResult );
    
        BOOST_CHECK( !error );
    
        uIntType indx = 0;
        for ( unsigned int j = i; j < nodeDataAnswer.size( ); j += 3, indx++ ){

            BOOST_CHECK( vectorTools::fuzzyEquals( nodeDataResult[ indx ], nodeDataAnswer2[ j ] ) );

        }

    }

    for ( uIntType i = 0; i < 2; i++ ){

        floatVector elementDataResult;
        std::string name = "TEST_DATA_" + std::to_string( i + 1 ) + "_";
        error = xdmf_result2.getSolutionData( increment, name, "Cell", elementDataResult );
    
        BOOST_CHECK( !error );

        uIntType indx = 0;
        for ( unsigned int j = i; j < elementDataAnswer.size( ); j += 2, indx++ ){
    
            BOOST_CHECK( vectorTools::fuzzyEquals( elementDataResult[ indx ], elementDataAnswer2[ j ] ) );

        }

    }

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

}
