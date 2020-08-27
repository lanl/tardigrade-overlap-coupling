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

    unsigned int numIncrementsAnswer = 2;
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

int test_XDMFDataFile_getSubDomainNodes( std::ofstream &results ){
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
    errorOut error = xdf.getSubDomainNodes( 0, domainName, domainNodesResult );

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
    error = xdf.getSubDomainNodes( 0, domainName, domainNodesResult );

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

int test_XDMFDataFile_getSetNames( std::ofstream &results ){
    /*!
     * Test the function to extract the names of the sets
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "testConfig.yaml" );
    dataFileInterface::XDMFDataFile xdf( yf[ "filetest1" ] );

    std::vector< std::string > answer = { "free_nodes", "ghost_nodes",
                                          "left", "right", "bottom", "top", "back", "front",
                                          "non_overlapped_elements", "free_elements", "ghost_elements" };
    std::vector< std::string > result;

    errorOut error = xdf.getSetNames( 1, result );

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
    dataFileInterface::XDMFDataFile xdf( yf[ "filetest1" ] );

    floatVector answer = { -0.001, -0.001, -0.001,
                           -0.001, -0.001, -0.001,
                           -0.001, -0.001, -0.001,
                           -0.001, -0.001, -0.001,
                           -0.001, -0.001, -0.001,
                           -0.001 };

    floatVector result;
    errorOut error = xdf.getSolutionData( 1, "disp_z", "Node", result );

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
    dataFileInterface::XDMFDataFile xdf( yf[ "filetest1" ] );

    floatVector nodePositionAnswer = { 1, 0, 1, 1, 0, 0, 0, 0, 0,
                                       0, 0, 1, 1, 1, 1, 1, 1, 0,
                                       0, 1, 0, 0, 1, 1, 0, 1, 2,
                                       1, 1, 2, 0, 0, 2, 1, 0, 2,
                                       0, 0, 3, 0, 1, 3, 1, 1, 3,
                                       1, 0, 3 };

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

    uIntVector connectivityCellIndicesAnswer = { 0, 32, 64 };

    unsigned int cellCountAnswer = 3;

    floatVector nodePositionResult;
    uIntVector connectivityResult, connectivityCellIndicesResult;
    unsigned int cellCountResult;

    errorOut error = xdf.getMeshData( 1, nodePositionResult, connectivityResult, connectivityCellIndicesResult, cellCountResult );

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
    dataFileInterface::XDMFDataFile xdf( yf[ "filetest1" ] );

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

    errorOut error = xdf.getMeshData( 1, nodePositionResult, connectivityResult, connectivityCellIndicesResult, cellCountResult );

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
    dataFileInterface::XDMFDataFile xdf( yf[ "filetest1" ] );

    unsigned int numSubDomainNodesAnswer = 8;

    unsigned int numSubDomainNodesResult;
    std::string domainName = "left";
    errorOut error = xdf.getNumSubDomainNodes( 0, domainName, numSubDomainNodesResult );

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
    error = xdf.getNumSubDomainNodes( 0, domainName, numSubDomainNodesResult );

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
    dataFileInterface::XDMFDataFile xdf( yf[ "filetest1" ] );

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
    errorOut error = xdf.getSolutionVectorDataFromComponents( 1, componentNames, "Node", result );

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
    dataFileInterface::XDMFDataFile xdf( yf[ "filetest1" ] );

    floatType answer1 = 0.;
    floatType result;
    errorOut error = xdf.getIncrementTime( 0, result );

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
    error = xdf.getIncrementTime( 1, result );

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
    dataFileInterface::XDMFDataFile xdf( yf[ "filetest3" ] );

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

    errorOut error = xdf.initializeIncrement( timeAnswer, reference_increment, collectionNumber, increment );

    if ( error ){

        error->print( );
        results << "test_writeIncrementMeshData & False\n";
        return 1;

    }

    std::remove( "test_output.xdmf" );
    std::remove( "test_output.h5" );

    error = xdf.writeIncrementMeshData( increment, collectionNumber, nodeIdsAnswer, { { } }, { { } }, nodePositionsAnswer,
                                        elementIdsAnswer, { { } }, { { } }, connectivityAnswer );

    if ( error ){

        error->print( );
        results << "test_writeIncrementMeshData & False\n";
        return 1;

    }

    //Read in the mesh data to determine if things have been stored correctly
    YAML::Node af = YAML::Load( "mode: read\nfilename: test_output.xdmf\ncell_id_variable_name: ELEMID\n" );
    dataFileInterface::XDMFDataFile xdf_result( af );

    if ( xdf_result._error ){

        xdf_result._error->print( );
        results << "test_writeIncrementMeshData & False\n";
        return 1;

    }

    //Check if the timestep was stored correctly
    floatType scalarResult;
    error = xdf_result.getIncrementTime( increment, scalarResult );

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
    error = xdf_result.getNodeIds( increment, "NODEID", nodeIdsResult );

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

    error = xdf_result.getMeshData( increment, nodePositionsResult, connectivityResult, cellIndicesResult, cellCountsResult );

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
    dataFileInterface::XDMFDataFile xdf( yf[ "filetest1" ] );

    if ( xdf._error ){

        xdf._error->print( );
        results << "test_XDMFDataFile_getNodeIds & False\n";
        return 1;
    }

    uIntVector nodeIdAnswer = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

    uIntVector nodeIdResult;
    errorOut error = xdf.getNodeIds( 0, "NODEID", nodeIdResult );

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
    test_XDMFDataFile_getSubDomainNodes( results );
    test_XDMFDataFile_getNumNodes( results );
    test_XDMFDataFile_getSetNames( results );
    test_XDMFDataFile_getSolutionData( results );
    test_XDMFDataFile_getSolutionVectorDataFromComponents( results );
    test_XDMFDataFile_getMeshData( results );
    test_XDMFDataFile_getMeshData2( results );
    test_XDMFDataFile_getIncrementTime( results );

    test_XDMFDataFile_writeIncrementMeshData( results );

    //Close the results file
    results.close();
}
