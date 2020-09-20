//!The test file for dataFileInterface.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>
#include <boost/algorithm/string.hpp>
#define USE_EIGEN
#include<vector_tools.h>

#include<overlapCoupling.h>

typedef overlapCoupling::errorNode errorNode; //!Redefinition for the error node
typedef overlapCoupling::errorOut errorOut; //!Redefinition for a pointer to the error node
typedef overlapCoupling::floatType floatType; //!Define the float values type.
typedef overlapCoupling::floatVector floatVector; //! Define a vector of floats
typedef overlapCoupling::floatMatrix floatMatrix; //!Define a matrix of floats
typedef overlapCoupling::uIntVector uIntVector; //!Define a vector of unsigned ints
typedef overlapCoupling::SparseMatrix SparseMatrix; //!Define a sparse matrix
typedef DOFProjection::T T; //!Define the triplet

errorOut readMatrixFromFile( const std::string filename, floatVector &data, Eigen::MatrixXd &matrix ){
    /*!
     * Read a matrix of doubles from a file
     */

    std::ifstream file( filename );
    std::string line;
    std::vector< std::string > splitLine;

    unsigned int rows = 0;
    unsigned int cols = 0;

    if ( file.is_open( ) ){

        while ( std::getline( file, line ) ){


            boost::algorithm::split( splitLine, line, boost::is_any_of( "," ) );
    
            if ( ( splitLine.size( ) != rows ) && ( rows != 0 ) ){
    
                return new errorNode( "readMatrixFromFile", "The matrix is not a consistent matrix" );
    
            }
            else{

                rows = splitLine.size( );

            }
    
            for ( auto v = splitLine.begin( ); v != splitLine.end( ); v++ ){
    
                data.push_back( std::stod( *v ) );
    
            }
    
            cols++;
    
        }

    }
    else{

        return new errorNode( "readMatrixFromFile", "Can't open " + filename );

    }

    if ( cols == 0 ){

        return new errorNode( "readMatrixFromFile", "There are no columns in the matrix" );

    }

    if ( rows == 0 ){

        return new errorNode( "readMatrixFromFile", "there are no rows in the matrix" );

    }

    matrix = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>( data.data(), rows, cols );

    return NULL;

}

int test_overlapCoupling_constructor( std::ostream &results ){
    /*!
     * Test the constructor to make sure that the code generates properly
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "../testFiles/testConfig.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_constructor & False\n";
        return 1;
    }

    results << "test_overlapCoupling_constructor & True\n";
    return 0;
}

int test_overlapCoupling_initializeCoupling( std::ofstream &results ){
    /*!
     * Test the initialization of the coupling
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string filename = "../testFiles/testConfig.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_initializeCoupling & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_initializeCoupling & False\n";
        return 1;
    }

    if ( !std::ifstream( "reference_information.xdmf" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling (test 1) & False\n";
        return 1;

    }

    if ( !std::ifstream( "reference_information.h5" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling (test 2) & False\n";
        return 1;

    }

    if ( !std::ifstream( "homogenized_response.xdmf" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling (test 3) & False\n";
        return 1;

    }

    if ( !std::ifstream( "homogenized_response.h5" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling (test 4) & False\n";
        return 1;

    }

    std::string xdmf_filename = "reference_information.xdmf";
    shared_ptr< XdmfReader > reader = XdmfReader::New( );
    shared_ptr< XdmfDomain > _readDomain = shared_dynamic_cast< XdmfDomain >( reader->read( xdmf_filename ) );
    shared_ptr< XdmfUnstructuredGrid > _readGrid = _readDomain->getUnstructuredGrid( 0 );

    SparseMatrix N;
    overlapCoupling::readSparseMatrixFromXDMF( _readGrid, "N", N );

    Eigen::MatrixXd macroD( N.cols( ), 1 );

    for ( unsigned int i = 0; i < N.cols( ); i+=12 ){

        macroD( i +  0, 0 ) = 1;
        macroD( i +  1, 0 ) = 2;
        macroD( i +  2, 0 ) = 3;
        macroD( i +  3, 0 ) = 0;
        macroD( i +  4, 0 ) = 0;
        macroD( i +  5, 0 ) = 0;
        macroD( i +  6, 0 ) = 0;
        macroD( i +  7, 0 ) = 0;
        macroD( i +  8, 0 ) = 0;
        macroD( i +  9, 0 ) = 0;
        macroD( i + 10, 0 ) = 0;
        macroD( i + 11, 0 ) = 0;

    }

    Eigen::MatrixXd R = N * macroD;

    for ( unsigned int i = 0; i < R.rows( ); i+=3 ){

        if ( !vectorTools::fuzzyEquals( R( i + 0, 0 ), 1. ) ){

            results << "test_overlapCoupling_initializeCoupling (test 5) & False\n";
            return 1;

        }
        if ( !vectorTools::fuzzyEquals( R( i + 1, 0 ), 2. ) ){

            results << "test_overlapCoupling_initializeCoupling (test 6) & False\n";
            return 1;

        }
        if ( !vectorTools::fuzzyEquals( R( i + 2, 0 ), 3. ) ){

            results << "test_overlapCoupling_initializeCoupling (test 7) & False\n";
            return 1;

        }

    }

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_initializeCoupling & True\n";
    return 0;
}

int test_overlapCoupling_getReferenceFreeMicroDomainMasses( std::ofstream &results ){
    /*!
     * Test the extraction of the reference free micro-domain mass
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string filename = "../testFiles/testConfig.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_getReferenceFreeMicroDomainMasses & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_getReferenceFreeMicroDomainMasses & False\n";
        return 1;
    }

    floatVector referenceFreeMicroDomainMassesAnswer = { 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125 };
    referenceFreeMicroDomainMassesAnswer *= 2;

    const floatVector *referenceFreeMicroDomainMassesResult = oc.getReferenceFreeMicroDomainMasses( );

    if ( !vectorTools::fuzzyEquals( referenceFreeMicroDomainMassesAnswer, *referenceFreeMicroDomainMassesResult ) ){
        vectorTools::print( referenceFreeMicroDomainMassesAnswer );
        vectorTools::print( *referenceFreeMicroDomainMassesResult );
        std::cout << referenceFreeMicroDomainMassesResult->size( ) << "\n";
        results << "test_overlapCoupling_getReferenceFreeMicroDomainMasses (test 1) & False\n";
        return 1;
    }

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_getReferenceFreeMicroDomainMasses & True\n";
    return 0;
}

int test_overlapCoupling_getReferenceGhostMicroDomainMasses( std::ofstream &results ){
    /*!
     * Test the extraction of the reference ghost micro-domain mass
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string filename = "../testFiles/testConfig.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_getReferenceGhostMicroDomainMasses & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_getReferenceGhostMicroDomainMasses & False\n";
        return 1;
    }

    floatVector referenceGhostMicroDomainMassesAnswer = { 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125 };
    referenceGhostMicroDomainMassesAnswer *= 2;

    const floatVector *referenceGhostMicroDomainMassesResult = oc.getReferenceGhostMicroDomainMasses( );

    if ( !vectorTools::fuzzyEquals( referenceGhostMicroDomainMassesAnswer, *referenceGhostMicroDomainMassesResult ) ){
        results << "test_overlapCoupling_getReferenceGhostMicroDomainMasses (test 1) & False\n";
        return 1;
    }

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_getReferenceGhostMicroDomainMasses & True\n";
    return 0;
}
int test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass( std::ofstream &results ){
    /*!
     * Test the extraction of the reference free micro-domain centers of mass
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string filename = "../testFiles/testConfig.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass & False\n";
        return 1;
    }

    floatVector referenceFreeMicroDomainCentersOfMassAnswer = { 0.75, 0.25, 2.75,
                                                                0.75, 0.25, 2.25,
                                                                0.25, 0.25, 2.75,
                                                                0.25, 0.25, 2.25,
                                                                0.75, 0.75, 2.75,
                                                                0.75, 0.75, 2.25,
                                                                0.25, 0.75, 2.75,
                                                                0.25, 0.75, 2.25 };

    const floatVector *referenceFreeMicroDomainCentersOfMassResult = oc.getReferenceFreeMicroDomainCentersOfMass( );

    if ( !vectorTools::fuzzyEquals( referenceFreeMicroDomainCentersOfMassAnswer, *referenceFreeMicroDomainCentersOfMassResult ) ){
        vectorTools::print( *referenceFreeMicroDomainCentersOfMassResult );
        vectorTools::print(  referenceFreeMicroDomainCentersOfMassAnswer );
        results << "test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass (test 1) & False\n";
        return 1;
    }

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass & True\n";
    return 0;
}

int test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass( std::ofstream &results ){
    /*!
     * Test the extraction of the reference ghost micro-domain centers of mass
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string filename = "../testFiles/testConfig.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass & False\n";
        return 1;
    }

    floatVector referenceGhostMicroDomainCentersOfMassAnswer = { 0.75, 0.25, 1.75,
                                                                 0.75, 0.25, 1.25,
                                                                 0.25, 0.25, 1.75,
                                                                 0.25, 0.25, 1.25,
                                                                 0.75, 0.75, 1.75,
                                                                 0.75, 0.75, 1.25,
                                                                 0.25, 0.75, 1.75,
                                                                 0.25, 0.75, 1.25 };

    const floatVector *referenceGhostMicroDomainCentersOfMassResult = oc.getReferenceGhostMicroDomainCentersOfMass( );

    if ( !vectorTools::fuzzyEquals( referenceGhostMicroDomainCentersOfMassAnswer, *referenceGhostMicroDomainCentersOfMassResult ) ){
        vectorTools::print( *referenceGhostMicroDomainCentersOfMassResult );
        vectorTools::print(  referenceGhostMicroDomainCentersOfMassAnswer );
        results << "test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass (test 1) & False\n";
        return 1;
    }

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass & True\n";
    return 0;
}

int test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions( std::ofstream &results ){
    /*!
     * Test of the retrieval of the reference free micro domain centers of mass shape functions
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string filename = "../testFiles/testConfig.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_ & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions & False\n";
        return 1;
    }

    const floatVector referenceFreeMicroDomainCenterOfMassShapeFunctionsAnswer =
        {
            0.140625, 0.046875, 0.140625, 0.421875, 0.046875, 0.015625, 0.046875, 0.140625,
            0.046875, 0.015625, 0.046875, 0.140625, 0.140625, 0.046875, 0.140625, 0.421875,
            0.421875, 0.140625, 0.046875, 0.140625, 0.140625, 0.046875, 0.015625, 0.046875,
            0.140625, 0.046875, 0.015625, 0.046875, 0.421875, 0.140625, 0.046875, 0.140625,
            0.046875, 0.140625, 0.421875, 0.140625, 0.015625, 0.046875, 0.140625, 0.046875,
            0.015625, 0.046875, 0.140625, 0.046875, 0.046875, 0.140625, 0.421875, 0.140625,
            0.140625, 0.421875, 0.140625, 0.046875, 0.046875, 0.140625, 0.046875, 0.015625,
            0.046875, 0.140625, 0.046875, 0.015625, 0.140625, 0.421875, 0.140625, 0.046875
        };

    const floatVector *referenceFreeMicroDomainCenterOfMassShapeFunctionsResult =
        oc.getReferenceFreeMicroDomainCenterOfMassShapeFunctions( );

    if ( !vectorTools::fuzzyEquals( referenceFreeMicroDomainCenterOfMassShapeFunctionsAnswer,
                                   *referenceFreeMicroDomainCenterOfMassShapeFunctionsResult ) ){
        vectorTools::print( *referenceFreeMicroDomainCenterOfMassShapeFunctionsResult );
        results << "test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions (test 1) & False\n";
        return 1;
    } 

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions & True\n";
    return 0;

}

int test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions( std::ofstream &results ){
    /*!
     * Test of the retrieval of the reference ghost micro domain centers of mass shape functions
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string filename = "../testFiles/testConfig.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_ & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions & False\n";
        return 1;
    }

    const floatVector referenceGhostMicroDomainCenterOfMassShapeFunctionsAnswer =
        {
            0.046875, 0.015625, 0.046875, 0.140625, 0.140625, 0.046875, 0.140625, 0.421875,
            0.015625, 0.046875, 0.140625, 0.046875, 0.046875, 0.140625, 0.421875, 0.140625,
            0.140625, 0.046875, 0.015625, 0.046875, 0.421875, 0.140625, 0.046875, 0.140625,
            0.046875, 0.140625, 0.046875, 0.015625, 0.140625, 0.421875, 0.140625, 0.046875,
            0.140625, 0.046875, 0.140625, 0.421875, 0.046875, 0.015625, 0.046875, 0.140625,
            0.046875, 0.140625, 0.421875, 0.140625, 0.015625, 0.046875, 0.140625, 0.046875,
            0.421875, 0.140625, 0.046875, 0.140625, 0.140625, 0.046875, 0.015625, 0.046875,
            0.140625, 0.421875, 0.140625, 0.046875, 0.046875, 0.140625, 0.046875, 0.015625
        };

    const floatVector *referenceGhostMicroDomainCenterOfMassShapeFunctionsResult =
        oc.getReferenceGhostMicroDomainCenterOfMassShapeFunctions( );

    if ( !vectorTools::fuzzyEquals( referenceGhostMicroDomainCenterOfMassShapeFunctionsAnswer,
                                   *referenceGhostMicroDomainCenterOfMassShapeFunctionsResult ) ){
        vectorTools::print( *referenceGhostMicroDomainCenterOfMassShapeFunctionsResult );
        results << "test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions (test 1) & False\n";
        return 1;
    } 

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions & True\n";
    return 0;

}

int test_overlapCoupling_processIncrement( std::ofstream &results ){
    /*!
     * Test the initialization of the coupling
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    remove( "macroscale_dof.xdmf" );
    remove( "macroscale_dof.h5" );

    remove( "microscale_dof.xdmf" );
    remove( "microscale_dof.h5" );

    std::string filename = "../testFiles/testConfig.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_processIncrement & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_processIncrement & False\n";
        return 1;
    }

    error = oc.processIncrement( 1, 1 );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_processIncrement & False\n";
        return 1;
    }

    const floatVector *projectedGhostMacroDisplacement = oc.getProjectedGhostMacroDisplacement( );
    const floatVector *projectedGhostMicroDisplacement = oc.getProjectedGhostMicroDisplacement( );

    if ( projectedGhostMacroDisplacement->size( ) == 0 ){
        results << "test_overlapCoupling_processIncrement (test 1) & False\n";
        return 1;
    }

    if ( projectedGhostMicroDisplacement->size( ) == 0 ){
        results << "test_overlapCoupling_processIncrement (test 2) & False\n";
        return 1;
    }

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    remove( "macroscale_dof.xdmf" );
    remove( "macroscale_dof.h5" );

    remove( "microscale_dof.xdmf" );
    remove( "microscale_dof.h5" );

    results << "test_overlapCoupling_processIncrement & True\n";
    return 0;
}

int test_overlapCoupling_processLastIncrements( std::ofstream &results ){
    /*!
     * Test processing the last increments
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    remove( "macroscale_dof.xdmf" );
    remove( "macroscale_dof.h5" );

    remove( "microscale_dof.xdmf" );
    remove( "microscale_dof.h5" );

    std::string filename = "../testFiles/testConfig.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_processLastIncrements & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_processLastIncrements & False\n";
        return 1;
    }

    error = oc.processLastIncrements( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_processLastIncrements & False\n";
        return 1;
    }

    const floatVector *projectedGhostMacroDisplacement = oc.getProjectedGhostMacroDisplacement( );
    const floatVector *projectedGhostMicroDisplacement = oc.getProjectedGhostMicroDisplacement( );

    if ( projectedGhostMacroDisplacement->size( ) == 0 ){
        results << "test_overlapCoupling_processLastIncrements (test 1) & False\n";
        return 1;
    }

    if ( projectedGhostMicroDisplacement->size( ) == 0 ){
        results << "test_overlapCoupling_processLastIncrements (test 2) & False\n";
        return 1;
    }

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    remove( "macroscale_dof.xdmf" );
    remove( "macroscale_dof.h5" );

    remove( "microscale_dof.xdmf" );
    remove( "microscale_dof.h5" );

    results << "test_overlapCoupling_processLastIncrements & True\n";
    return 0;
}

int test_MADOutlierDetection( std::ofstream &results ){
    /*!
     * Test the computation of outliers using the maximum absolution deviation 
     * detection metric.
     * 
     * :param std::ofstream &results: The output file.
     */

    floatVector x = { 0.70154526, 0.00265005, 0.29766985, 0.0570927 , 0.12136678 };
    uIntVector outliers;

    errorOut error = overlapCoupling::MADOutlierDetection( x, outliers, 5 );

    if ( error ){

        error->print( );
        results << "test_MADOutlierDetection & False\n";
        return 1;

    }

    if ( outliers.size() > 0 ){
        results << "test_MADOutlierDetection (test 1) & False\n";
        return 1;
    }

    overlapCoupling::MADOutlierDetection( x, outliers, 4 );

    if ( !vectorTools::fuzzyEquals( outliers, { 0 } ) ){

        results << "test_MADOutlierDetection (test 1) & False\n";
        return 1;

    }

    results << "test_MADOutlierDetection & True\n";
    return 0;
}

int test_formMicromorphicElementMassMatrix( std::ofstream &results ){
    /*!
     * Test the formation of the mass matrix for a single micromorphic element.
     *
     * :param std::ofstream &results: The output file
     */

    std::unique_ptr< elib::Element > element;

    floatMatrix reference_nodes =
        {
            { 0, 0, 0 },
            { 1, 0, 0 },
            { 1, 1, 0 },
            { 0, 1, 0 },
            { 0, 0, 1 },
            { 1, 0, 1 },
            { 1, 1, 1 },
            { 0, 1, 1 }
        };

    floatVector degreeOfFreedomValues =
        {
            -0.04066559,  0.0390943 , -0.00232655, -0.03264781,  0.02111065,
             0.00079745, -0.00409562,  0.03073762,  0.01799572, -0.02048283,
             0.00914281, -0.04062647,  0.02985775, -0.02334215,  0.01175122,
            -0.02032996,  0.0489607 , -0.01085405, -0.00867754,  0.00475848,
             0.04859963,  0.00442464,  0.0211827 ,  0.04006121, -0.03019419,
            -0.04428046,  0.04353155, -0.03999713,  0.03122015,  0.01574864,
            -0.0025961 , -0.0232323 , -0.03535742, -0.00621033, -0.00485358,
             0.00917265, -0.03230599,  0.0228427 , -0.00198795,  0.04252169,
             0.04223892,  0.01351459,  0.01295129, -0.04424972,  0.0322117 ,
            -0.03558341, -0.01547058,  0.01382653, -0.0379305 ,  0.04802716,
            -0.03400389,  0.0348995 , -0.02256884, -0.04231619,  0.00654967,
             0.01205778, -0.01045064, -0.01303744, -0.00900963, -0.01305943,
            -0.04141735, -0.02224897,  0.03300644, -0.03487759, -0.04996436,
            -0.00305152, -0.00471929, -0.04222132, -0.03447746, -0.03298123,
            -0.02164259, -0.03886441,  0.02980088,  0.00453406,  0.03759968,
             0.00967625,  0.04231608, -0.001698  ,  0.01340712, -0.03988135,
            -0.04814852,  0.0196056 , -0.01305974, -0.0333702 , -0.03428086,
            -0.04671768,  0.02380142, -0.02290295,  0.02739323,  0.00481225,
             0.0424606 , -0.00402194, -0.00948801,  0.0395755 ,  0.02496865,
            -0.04143367
        };

    floatMatrix nodeDisplacement;
    nodeDisplacement.reserve( 8 );
    for ( unsigned int i = 0; i < 8; i++ ){

        nodeDisplacement.push_back( floatVector( degreeOfFreedomValues.begin( ) + 12 * i,
                                                 degreeOfFreedomValues.begin( ) + 12 * i + 3 ) );

    }

    auto qrule = elib::default_qrules.find( "Hex8" );
    element = elib::build_element_from_string( "Hex8", { 10, 7, 3, 9, 1, 8, 13, 4 }, reference_nodes, qrule->second );
    element->update_node_positions( nodeDisplacement );

    overlapCoupling::DOFMap nodeIDToIndex =
        {
            { 10, 0 },
            {  7, 2 },
            {  3, 4 },
            {  9, 3 },
            {  1, 6 },
            {  8, 1 },
            { 13, 7 },
            {  4, 5 }
        };

    floatVector density =
        { 693.53490713, 1765.4802207 ,   91.36052518,  666.64526727,
          51.16415254,  398.63874113,  702.24020488, 1190.92397094 };

    floatVector momentOfInertia =
        {
            0.44595488,  0.13676299, -0.2525482 ,  0.13676299, -0.19746052,
           -0.30931581, -0.2525482 , -0.30931581,  0.17223508, -0.31771288,
            0.3408996 ,  0.11782056,  0.3408996 ,  0.38466266,  0.1586499 ,
            0.11782056,  0.1586499 ,  0.19104204,  0.11479176,  0.17524618,
            0.04975924,  0.17524618, -0.17921051,  0.15571073,  0.04975924,
            0.15571073, -0.0283223 , -0.28271601,  0.18506962,  0.26725315,
            0.18506962,  0.17534291, -0.12381001,  0.26725315, -0.12381001,
            0.33227658,  0.32864468, -0.16755297, -0.25824399, -0.16755297,
            0.10339622, -0.40619467, -0.25824399, -0.40619467, -0.19037855,
            0.27080415,  0.35140531, -0.00926281,  0.35140531, -0.1955969 ,
           -0.0311474 , -0.00926281, -0.0311474 , -0.31288414,  0.41621005,
            0.10722768,  0.16218443,  0.10722768, -0.00765149,  0.06275192,
            0.16218443,  0.06275192, -0.31463761,  0.49814788, -0.14514796,
            0.24525217, -0.14514796, -0.46891292,  0.10017765,  0.24525217,
            0.10017765, -0.21952903
        };

    floatVector answerData;
    Eigen::MatrixXd answer;

    errorOut error = readMatrixFromFile( "mass_matrix_answer.csv", answerData, answer );

    if ( error ){

        error->print( );
        results << "test_formMicromorphicElementMassMatrix & False\n";
        return 1;

    }

    std::vector< DOFProjection::T > coefficients;

    error = overlapCoupling::formMicromorphicElementMassMatrix( element, degreeOfFreedomValues,
                                                                momentOfInertia, density, &nodeIDToIndex, coefficients );
    DOFProjection::SparseMatrix result( 8 * 12, 8 * 12 );
    result.setFromTriplets( coefficients.begin( ), coefficients.end( ) );

    if ( error ){

        error->print( );
        results << "test_formMicromorphicElementMassMatrix & False\n";
        return 1;

    }

    if ( !answer.isApprox( result.toDense( ), 1e-5 ) ){

        results << "test_formMicromorphicElementMassMatrix (test 1) & False\n";
        return 1;

    }

    results << "test_formMicromorphicElementMassMatrix & True\n";
    return 0;
}

int test_computeMicromorphicElementRequiredValues( std::ofstream &results ){
    /*!
     * Test the computation of the default required values from the micromorphic element
     *
     * :param std::ofstream &results: The output file
     */

    std::unique_ptr< elib::Element > element;

    floatMatrix reference_nodes =
        {
            { 0, 0, 0 },
            { 1, 0, 0 },
            { 1, 1, 0 },
            { 0, 1, 0 },
            { 0, 0, 1 },
            { 1, 0, 1 },
            { 1, 1, 1 },
            { 0, 1, 1 }
        };

    floatVector degreeOfFreedomValues =
        {
            -0.1       , -0.1       , -0.1       , -0.03264781,  0.02111065,
             0.00079745, -0.00409562,  0.03073762,  0.01799572, -0.02048283,
             0.00914281, -0.04062647,  0.1       , -0.1       , -0.1       ,
            -0.02032996,  0.0489607 , -0.01085405, -0.00867754,  0.00475848,
             0.04859963,  0.00442464,  0.0211827 ,  0.04006121,  0.1       ,
             0.1       , -0.1       , -0.03999713,  0.03122015,  0.01574864,
            -0.0025961 , -0.0232323 , -0.03535742, -0.00621033, -0.00485358,
             0.00917265, -0.1       ,  0.1       , -0.1       ,  0.04252169,
             0.04223892,  0.01351459,  0.01295129, -0.04424972,  0.0322117 ,
            -0.03558341, -0.01547058,  0.01382653, -0.1       , -0.1       ,
             0.1       ,  0.0348995 , -0.02256884, -0.04231619,  0.00654967,
             0.01205778, -0.01045064, -0.01303744, -0.00900963, -0.01305943,
             0.1       , -0.1       ,  0.1       , -0.03487759, -0.04996436,
            -0.00305152, -0.00471929, -0.04222132, -0.03447746, -0.03298123,
            -0.02164259, -0.03886441,  0.1       ,  0.1       ,  0.1       ,
             0.00967625,  0.04231608, -0.001698  ,  0.01340712, -0.03988135,
            -0.04814852,  0.0196056 , -0.01305974, -0.0333702 , -0.1       ,
             0.1       ,  0.1       , -0.02290295,  0.02739323,  0.00481225,
             0.0424606 , -0.00402194, -0.00948801,  0.0395755 ,  0.02496865,
            -0.04143367
        };

    floatMatrix nodeDisplacement;
    nodeDisplacement.reserve( 8 );
    for ( unsigned int i = 0; i < 8; i++ ){

        nodeDisplacement.push_back( floatVector( degreeOfFreedomValues.begin( ) + 12 * i,
                                                 degreeOfFreedomValues.begin( ) + 12 * i + 3 ) );

    }

    auto qrule = elib::default_qrules.find( "Hex8" );
    element = elib::build_element_from_string( "Hex8", { 10, 7, 3, 9, 1, 8, 13, 4 }, reference_nodes, qrule->second );
    element->update_node_positions( nodeDisplacement );

    floatMatrix reshapedDOFValues = vectorTools::inflate( degreeOfFreedomValues, 8, 12 );
    floatVector uQptResult, XiQptResult, shapeFunctionsResult, deformationGradientResult;
    floatMatrix gradShapeFunctionsResult;
    floatType JResult, JxwResult;

    floatType JAnswer = 1.728;
    floatType JxwReferenceAnswer = 0.125;
    floatType JxwCurrentAnswer = 0.216;
    floatType tmp = 0.0577350269;

    floatMatrix shapeFunctionAnswer
        {
            { 0.490563,   0.131446,   0.0352208,  0.131446,   0.131446,   0.0352208,  0.00943739, 0.0352208 }, 
            { 0.131446,   0.490563,   0.131446,   0.0352208,  0.0352208,  0.131446,   0.0352208,  0.00943739 },
            { 0.0352208,  0.131446,   0.490563,   0.131446,   0.00943739, 0.0352208,  0.131446,   0.0352208 },
            { 0.131446,   0.0352208,  0.131446,   0.490563,   0.0352208,  0.00943739, 0.0352208,  0.131446 },
            { 0.131446,   0.0352208,  0.00943739, 0.0352208,  0.490563,   0.131446,   0.0352208,  0.131446 },
            { 0.0352208,  0.131446,   0.0352208,  0.00943739, 0.131446,   0.490563,   0.131446,   0.0352208 },
            { 0.00943739, 0.0352208,  0.131446,   0.0352208,  0.0352208,  0.131446,   0.490563,   0.131446 },
            { 0.0352208,  0.00943739, 0.0352208,  0.131446,   0.131446,   0.0352208,  0.131446,   0.490563 },
        };

    floatMatrix gradShapeFunctionsReferenceAnswer
        {
            { -0.622008, -0.622008, -0.622008, 0.622008, -0.166667, -0.166667, 0.166667, 0.166667, -0.0446582, -0.166667, 0.622008, -0.166667, -0.166667, -0.166667, 0.622008, 0.166667, -0.0446582, 0.166667, 0.0446582, 0.0446582, 0.0446582, -0.0446582, 0.166667, 0.166667 },
            { -0.622008, -0.166667, -0.166667, 0.622008, -0.622008, -0.622008, 0.166667, 0.622008, -0.166667, -0.166667, 0.166667, -0.0446582, -0.166667, -0.0446582, 0.166667, 0.166667, -0.166667, 0.622008, 0.0446582, 0.166667, 0.166667, -0.0446582, 0.0446582, 0.0446582 },
            { -0.166667, -0.166667, -0.0446582, 0.166667, -0.622008, -0.166667, 0.622008, 0.622008, -0.622008, -0.622008, 0.166667, -0.166667, -0.0446582, -0.0446582, 0.0446582, 0.0446582, -0.166667, 0.166667, 0.166667, 0.166667, 0.622008, -0.166667, 0.0446582, 0.166667 },
            { -0.166667, -0.622008, -0.166667, 0.166667, -0.166667, -0.0446582, 0.622008, 0.166667, -0.166667, -0.622008, 0.622008, -0.622008, -0.0446582, -0.166667, 0.166667, 0.0446582, -0.0446582, 0.0446582, 0.166667, 0.0446582, 0.166667, -0.166667, 0.166667, 0.622008 },
            { -0.166667, -0.166667, -0.622008, 0.166667, -0.0446582, -0.166667, 0.0446582, 0.0446582, -0.0446582, -0.0446582, 0.166667, -0.166667, -0.622008, -0.622008, 0.622008, 0.622008, -0.166667, 0.166667, 0.166667, 0.166667, 0.0446582, -0.166667, 0.622008, 0.166667 },
            { -0.166667, -0.0446582, -0.166667, 0.166667, -0.166667, -0.622008, 0.0446582, 0.166667, -0.166667, -0.0446582, 0.0446582, -0.0446582, -0.622008, -0.166667, 0.166667, 0.622008, -0.622008, 0.622008, 0.166667, 0.622008, 0.166667, -0.166667, 0.166667, 0.0446582 },
            { -0.0446582, -0.0446582, -0.0446582, 0.0446582, -0.166667, -0.166667, 0.166667, 0.166667, -0.622008, -0.166667, 0.0446582, -0.166667, -0.166667, -0.166667, 0.0446582, 0.166667, -0.622008, 0.166667, 0.622008, 0.622008, 0.622008, -0.622008, 0.166667, 0.166667 },
            { -0.0446582, -0.166667, -0.166667, 0.0446582, -0.0446582, -0.0446582, 0.166667, 0.0446582, -0.166667, -0.166667, 0.166667, -0.622008, -0.166667, -0.622008, 0.166667, 0.166667, -0.166667, 0.0446582, 0.622008, 0.166667, 0.166667, -0.622008, 0.622008, 0.622008 }
        };

    floatMatrix gradShapeFunctionsCurrentAnswer
        {
            { -0.51834, -0.51834, -0.51834, 0.51834, -0.138889, -0.138889, 0.138889, 0.138889, -0.0372152, -0.138889, 0.51834, -0.138889, -0.138889, -0.138889, 0.51834, 0.138889, -0.0372152, 0.138889, 0.0372152, 0.0372152, 0.0372152, -0.0372152, 0.138889, 0.138889 },
            { -0.51834, -0.138889, -0.138889, 0.51834, -0.51834, -0.51834, 0.138889, 0.51834, -0.138889, -0.138889, 0.138889, -0.0372152, -0.138889, -0.0372152, 0.138889, 0.138889, -0.138889, 0.51834, 0.0372152, 0.138889, 0.138889, -0.0372152, 0.0372152, 0.0372152 },
            { -0.138889, -0.138889, -0.0372152, 0.138889, -0.51834, -0.138889, 0.51834, 0.51834, -0.51834, -0.51834, 0.138889, -0.138889, -0.0372152, -0.0372152, 0.0372152, 0.0372152, -0.138889, 0.138889, 0.138889, 0.138889, 0.51834, -0.138889, 0.0372152, 0.138889 },
            { -0.138889, -0.51834, -0.138889, 0.138889, -0.138889, -0.0372152, 0.51834, 0.138889, -0.138889, -0.51834, 0.51834, -0.51834, -0.0372152, -0.138889, 0.138889, 0.0372152, -0.0372152, 0.0372152, 0.138889, 0.0372152, 0.138889, -0.138889, 0.138889, 0.51834 },
            { -0.138889, -0.138889, -0.51834, 0.138889, -0.0372152, -0.138889, 0.0372152, 0.0372152, -0.0372152, -0.0372152, 0.138889, -0.138889, -0.51834, -0.51834, 0.51834, 0.51834, -0.138889, 0.138889, 0.138889, 0.138889, 0.0372152, -0.138889, 0.51834, 0.138889 },
            { -0.138889, -0.0372152, -0.138889, 0.138889, -0.138889, -0.51834, 0.0372152, 0.138889, -0.138889, -0.0372152, 0.0372152, -0.0372152, -0.51834, -0.138889, 0.138889, 0.51834, -0.51834, 0.51834, 0.138889, 0.51834, 0.138889, -0.138889, 0.138889, 0.0372152 },
            { -0.0372152, -0.0372152, -0.0372152, 0.0372152, -0.138889, -0.138889, 0.138889, 0.138889, -0.51834, -0.138889, 0.0372152, -0.138889, -0.138889, -0.138889, 0.0372152, 0.138889, -0.51834, 0.138889, 0.51834, 0.51834, 0.51834, -0.51834, 0.138889, 0.138889 },
            { -0.0372152, -0.138889, -0.138889, 0.0372152, -0.0372152, -0.0372152, 0.138889, 0.0372152, -0.138889, -0.138889, 0.138889, -0.51834, -0.138889, -0.51834, 0.138889, 0.138889, -0.138889, 0.0372152, 0.51834, 0.138889, 0.138889, -0.51834, 0.51834, 0.51834 },
        };

    floatVector deformationGradientAnswer = { 1.2, 0.0, 0.0,
                                              0.0, 1.2, 0.0,
                                              0.0, 0.0, 1.2 };

    floatMatrix uQptAnswer =
        {
            { -tmp, -tmp, -tmp },
            {  tmp, -tmp, -tmp },
            {  tmp,  tmp, -tmp },
            { -tmp,  tmp, -tmp },
            { -tmp, -tmp,  tmp },
            {  tmp, -tmp,  tmp },
            {  tmp,  tmp,  tmp },
            { -tmp,  tmp,  tmp },
        };

    floatMatrix XiQptAnswer =
        {
            {  9.88136118e-01,  2.00813257e-02, -4.22070732e-03,  7.77908101e-04,  1.00864936e+00,  1.48284518e-02, -1.56589841e-02,  3.87462985e-03,  9.82616852e-01 },
            {  9.78744967e-01,  2.67710907e-02, -4.57959127e-03,  -4.19704309e-03,  9.95194617e-01,  1.60081991e-02, -6.32184786e-03,  7.02386931e-03,  1.00887017e+00 },
            { 9.81712207e-01,  3.26011268e-02,  7.54296900e-03,   2.29737042e-03,  9.77737555e-01, -1.40649506e-02, -5.17738385e-03, -2.99266479e-03,  1.00281440e+00 },
            {  1.00882498e+00,  3.31486739e-02,  7.47592280e-03,  1.14078262e-02,  9.77539770e-01,  1.15950117e-02, -1.56866541e-02, -3.97893908e-03,  9.96610936e-01 },
            {  1.00597891e+00, -6.26624679e-03, -2.02399680e-02,  8.23385344e-03,  1.00086213e+00, -7.72359125e-03, -8.68680450e-03, -3.08540289e-03,  9.78507688e-01 },
            {  9.83113378e-01, -1.22727139e-02, -7.82937249e-03,  5.49537395e-04,  9.75961052e-01, -1.88694411e-02, -1.46165030e-02, -9.84911837e-03,  9.77660786e-01 },
            {  9.93597014e-01,  2.45120434e-02,  7.93844761e-05,  1.15392913e-02,  9.70627254e-01, -3.13985669e-02,  7.91831580e-03, -6.63726565e-03,  9.75335260e-01 },
            {  9.96234380e-01,  2.21313473e-02, -1.27548343e-03,  2.46714378e-02,  9.87375422e-01, -9.49015591e-03,  1.35403120e-02,  6.90293147e-03,  9.73290037e-01 }
        };

    errorOut error = NULL;
    for ( auto qpt = element->qrule.begin( ); qpt != element->qrule.end( ); qpt++ ){

        error = overlapCoupling::computeMicromorphicElementRequiredValues( element, qpt, 3, reshapedDOFValues, true,
                                                                           shapeFunctionsResult, gradShapeFunctionsResult,
                                                                           deformationGradientResult,
                                                                           JResult, JxwResult, uQptResult, XiQptResult );

        if ( error ){
        
            error->print( );
            results << "test_computeMicromorphicElementRequiredValues & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( shapeFunctionAnswer[ qpt - element->qrule.begin( ) ], shapeFunctionsResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 1) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( gradShapeFunctionsReferenceAnswer[ qpt - element->qrule.begin( ) ],
                                        vectorTools::appendVectors( gradShapeFunctionsResult ) ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 2) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( deformationGradientAnswer, deformationGradientResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 3) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( JResult, JAnswer ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 4) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( JxwResult, JxwReferenceAnswer ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 5) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( uQptAnswer[ qpt - element->qrule.begin( ) ], uQptResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 6) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( XiQptAnswer[ qpt - element->qrule.begin( ) ], XiQptResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 7) & False\n";
            return 1;

        }

        error = overlapCoupling::computeMicromorphicElementRequiredValues( element, qpt, 3, reshapedDOFValues, false,
                                                                           shapeFunctionsResult, gradShapeFunctionsResult,
                                                                           deformationGradientResult,
                                                                           JResult, JxwResult, uQptResult, XiQptResult );

        if ( error ){
        
            error->print( );
            results << "test_computeMicromorphicElementRequiredValues & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( shapeFunctionAnswer[ qpt - element->qrule.begin( ) ], shapeFunctionsResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 8) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( gradShapeFunctionsCurrentAnswer[ qpt - element->qrule.begin( ) ],
                                        vectorTools::appendVectors( gradShapeFunctionsResult ) ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 9) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( deformationGradientAnswer, deformationGradientResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 10) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( JResult, JAnswer ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 11) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( JxwResult, JxwCurrentAnswer ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 12) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( uQptAnswer[ qpt - element->qrule.begin( ) ], uQptResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 13) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( XiQptAnswer[ qpt - element->qrule.begin( ) ], XiQptResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 14) & False\n";
            return 1;

        }

    }

    results << "test_computeMicromorphicElementRequiredValues & True\n";
    return 0;

}

int test_computeMicromorphicElementInternalForceVector( std::ofstream &results ){
    /*!
     * Test the computation of the micromorphic internal force vector
     *
     * :param std::ofstream &results: The output file
     */

    std::unique_ptr< elib::Element > element;

    floatMatrix reference_nodes =
        {
            { 0, 0, 0 },
            { 1, 0, 0 },
            { 1, 1, 0 },
            { 0, 1, 0 },
            { 0, 0, 1 },
            { 1, 0, 1 },
            { 1, 1, 1 },
            { 0, 1, 1 }
        };

    floatVector degreeOfFreedomValues =
        {
            -0.04066559,  0.0390943 , -0.00232655, -0.03264781,  0.02111065,
             0.00079745, -0.00409562,  0.03073762,  0.01799572, -0.02048283,
             0.00914281, -0.04062647,  0.02985775, -0.02334215,  0.01175122,
            -0.02032996,  0.0489607 , -0.01085405, -0.00867754,  0.00475848,
             0.04859963,  0.00442464,  0.0211827 ,  0.04006121, -0.03019419,
            -0.04428046,  0.04353155, -0.03999713,  0.03122015,  0.01574864,
            -0.0025961 , -0.0232323 , -0.03535742, -0.00621033, -0.00485358,
             0.00917265, -0.03230599,  0.0228427 , -0.00198795,  0.04252169,
             0.04223892,  0.01351459,  0.01295129, -0.04424972,  0.0322117 ,
            -0.03558341, -0.01547058,  0.01382653, -0.0379305 ,  0.04802716,
            -0.03400389,  0.0348995 , -0.02256884, -0.04231619,  0.00654967,
             0.01205778, -0.01045064, -0.01303744, -0.00900963, -0.01305943,
            -0.04141735, -0.02224897,  0.03300644, -0.03487759, -0.04996436,
            -0.00305152, -0.00471929, -0.04222132, -0.03447746, -0.03298123,
            -0.02164259, -0.03886441,  0.02980088,  0.00453406,  0.03759968,
             0.00967625,  0.04231608, -0.001698  ,  0.01340712, -0.03988135,
            -0.04814852,  0.0196056 , -0.01305974, -0.0333702 , -0.03428086,
            -0.04671768,  0.02380142, -0.02290295,  0.02739323,  0.00481225,
             0.0424606 , -0.00402194, -0.00948801,  0.0395755 ,  0.02496865,
            -0.04143367
        };

    floatMatrix nodeDisplacement;
    nodeDisplacement.reserve( 8 );
    for ( unsigned int i = 0; i < 8; i++ ){

        nodeDisplacement.push_back( floatVector( degreeOfFreedomValues.begin( ) + 12 * i,
                                                 degreeOfFreedomValues.begin( ) + 12 * i + 3 ) );

    }

    auto qrule = elib::default_qrules.find( "Hex8" );
    element = elib::build_element_from_string( "Hex8", { 10, 7, 3, 9, 1, 8, 13, 4 }, reference_nodes, qrule->second );
    element->update_node_positions( nodeDisplacement );

    overlapCoupling::DOFMap nodeIDToIndex =
        {
            { 10, 0 },
            {  7, 2 },
            {  3, 4 },
            {  9, 3 },
            {  1, 6 },
            {  8, 1 },
            { 13, 7 },
            {  4, 5 }
        };

    floatVector cauchyStress =
        {
            -0.45969764,  1.14822033, -1.36295921, -0.5000321 , -1.42325377,
            -0.90204189, -0.44663969,  0.28638596,  1.33786465,  1.16783921,
            -0.98509841, -0.34726126, -0.59238156, -1.46301493, -0.53019081,
            -0.39287991,  0.51670525, -0.95313864,  0.81999322,  0.5800708 ,
             1.2254044 ,  1.15398799, -0.20563341,  0.12685255,  1.15753279,
             0.57312611,  0.34437528,  1.01836067, -0.59267923,  0.03311643,
            -0.53763152, -0.0644638 ,  0.89123469,  0.99872816, -0.08496691,
            -0.84814803,  0.39031588,  0.44741445,  0.72010183,  0.52602873,
            -0.33027138,  1.15610447, -0.41258865,  0.99682975,  0.98593404,
             0.40067375, -1.3895129 , -1.40599158,  0.31018865, -1.30229246,
            -0.15995105, -1.13487559, -0.57113566,  0.81241244,  1.15210152,
             1.13622788,  0.95633037, -0.58679069,  0.90543097, -1.19180477,
            -1.09231537, -0.93865242,  1.24863015, -0.2600532 ,  1.11112689,
             0.90341929, -0.02295448,  0.42265386, -0.1380203 , -0.04614231,
            -0.94391283,  1.49444364
        };

    floatVector symmetricMicroStress =
        {
            0.10117987,  1.37582714, -0.97206764, -1.40673876,  0.47691406,
            0.09568587,  0.12415184, -0.6585871 , -0.63620366, -1.00348131,
           -0.05910748, -1.22504474,  0.95867523,  0.02521843,  0.62481998,
           -0.33498708,  0.72102278, -0.39559034, -1.37320129, -1.32955232,
            0.5171092 ,  1.13703248, -1.35457667,  0.53214959, -0.02469652,
            0.70267216,  0.03722338, -1.31772695, -0.21335467,  0.81220832,
           -0.98121639, -0.887004  , -0.46146348,  1.48095631,  0.90664032,
            0.9653973 , -0.68733076, -0.46874014,  0.59222427,  0.26895499,
           -1.4239325 , -0.74967125, -0.64881934, -1.10759892, -1.33468558,
           -0.34723361, -0.74556723, -0.01161972, -0.93394231,  0.40266023,
            0.39595606,  0.10547423,  1.43015732,  1.07281307, -0.01329458,
            1.27192279,  0.03659976,  1.46190328,  0.55927909, -0.78215141,
           -0.51928948, -1.18045737,  0.03809132, -0.60385602, -0.00554656,
            0.19209442,  1.02908222,  0.82778923, -1.38683747, -0.25118759,
            0.28963075, -0.85042722
        };

    floatVector higherOrderStress =
        {
            0.39229651, -1.27188227, -0.75057189, -1.36700197, -1.15241392,
           -0.83755545, -0.21498426,  1.02131234, -1.24288282, -0.76007716,
           -0.92450045, -1.04121465,  1.14264551, -0.94492367,  1.02005295,
            0.13905738, -1.08039572, -1.27358044,  0.7775386 , -1.49857357,
            0.68073573, -0.9120657 , -1.196205  , -0.27346835, -0.86683224,
            0.11928698, -1.40885654, -0.25191655, -1.21035354, -0.13516862,
            0.42551714, -1.367081  ,  1.41603026, -0.29139481,  0.66615639,
            1.45493764, -1.00795599, -1.08847277, -0.40932263,  1.37012762,
            1.3090307 , -1.0598115 ,  1.42470902, -1.10917598, -0.15906364,
            1.17699552, -0.74454968,  0.37418966,  0.5419793 ,  0.52773656,
            0.53722857,  0.60228796,  0.15701482, -1.40756007, -0.84260176,
            0.46713455, -0.9109566 , -0.4587495 , -0.27443065,  1.0929497 ,
           -0.56260743, -0.37428778,  1.1547629 , -1.26992331,  0.02180244,
           -0.99848932, -1.10047995, -0.73077205, -1.48752174,  0.36737644,
           -0.58238242,  1.01795009, -0.66254986, -1.10319654,  0.66144406,
            0.0965511 , -0.49245866, -0.35159224, -0.35986417,  0.52077983,
            0.37951793,  1.02467304,  1.12864357, -0.07594849,  1.31868883,
            0.58521494, -0.03745289, -0.31158896, -0.17807641, -0.91852362,
            0.03084061,  0.92307018,  0.0202181 ,  1.29077688,  1.30206557,
            0.0921027 ,  1.27553382,  0.34262177, -1.48643965,  0.61201724,
            1.49347647,  0.64445638,  0.28889739, -1.43068598,  0.03570288,
           -0.08223876, -0.44111068, -1.09584718, -0.19206467,  0.46059502,
           -0.40809829,  1.08533245,  0.51658791,  0.7905599 , -0.96952095,
            1.23554638,  0.87045523, -0.84229178, -1.44913607, -0.18983496,
            0.10093825, -1.1185102 ,  0.44234559,  0.60141309, -0.89067448,
            0.31832896, -0.55737341,  1.12297842, -1.23493082, -1.41979439,
            0.18887879, -0.70628017,  0.14308575, -0.39197145, -0.5690881 ,
            0.57602306,  0.99915749,  0.86749036, -0.43421454,  1.43900802,
            0.29587861,  1.48878773, -0.08768882, -0.14582274,  1.1333723 ,
            0.24309208,  0.23585109,  0.81010256, -0.27329013, -0.23483857,
           -0.25537111,  1.31371579, -0.36190462, -0.67439408,  1.15628843,
           -1.09453436, -0.16136072, -0.76631781, -1.3900636 ,  0.51565904,
            0.51245367,  1.41764668, -0.58441173,  1.34112171,  1.31309759,
           -0.39597866, -0.14284022,  1.41928088, -0.10742838,  0.9797636 ,
            1.0535908 , -0.63374945, -1.13702384,  0.82982881, -0.50849777,
           -0.00453001,  1.19397083,  0.44670022,  1.28185961, -0.45946417,
            1.09819494,  0.75844966,  0.30450907, -1.46525729, -1.07253823,
            1.31673755,  0.56305741,  0.7642357 , -0.42590002,  0.30657098,
           -1.2917611 ,  0.94766539,  0.89394059, -0.30803095, -0.16673604,
           -1.22619946, -0.46748525, -1.229485  , -0.73181422,  1.213918  ,
           -0.17310261,  0.28313866, -1.43347574,  0.99141871, -0.79494399,
           -0.0173043 ,  0.80792229, -0.79151497,  1.02627251, -1.09927188,
            0.0545878 ,  0.99842225,  0.3788543 , -0.01641077, -0.60538728,
            0.91039346
        };


    floatVector answer =
        {
            0.02739017,  0.08038737, -0.08497879,  0.00601456,  0.08946546,
            0.27811672,  0.01057751,  0.38096787, -0.13814203, -0.04606262,
            0.00897415,  0.18528464, -0.03866123,  0.07169987,  0.13703867,
           -0.00299104,  0.18306985,  0.04642505, -0.03120374,  0.18365814,
            0.08300798,  0.0682409 ,  0.07193369, -0.0009522 ,  0.2231008 ,
            0.16225866, -0.13121907, -0.03476567,  0.01903948,  0.12197495,
           -0.03377065,  0.0310203 ,  0.33407419, -0.33467015,  0.26452521,
            0.06854623, -0.27872297, -0.1310955 , -0.22403472, -0.34850733,
           -0.29972842, -0.022039  ,  0.18845184,  0.06034523, -0.01170348,
            0.2754791 , -0.12785453, -0.138628  ,  0.22756336, -0.09145956,
           -0.11057031, -0.47555049, -0.00796385, -0.26483904,  0.15577536,
            0.10029588, -0.04724648, -0.0130081 , -0.15922499,  0.03647805,
           -0.0681886 , -0.25736625,  0.03305347, -0.23813408,  0.16353921,
           -0.21829181, -0.1012466 , -0.2838822 ,  0.09664438,  0.08055681,
           -0.25891319, -0.2240477 , -0.19875004,  0.1424041 ,  0.19372511,
            0.05390138, -0.03151748, -0.05120333, -0.36223144,  0.00697449,
           -0.431687  , -0.0595205 , -0.2103422 , -0.43627243,  0.10626851,
            0.0231713 ,  0.18698564, -0.15977423,  0.12031186,  0.28120184,
           -0.01767001, -0.20123309,  0.29285515, -0.07107898,  0.31030247,
           -0.16779977
        };

    Eigen::MatrixXd internalForceVectorAnswer = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>( answer.data(), 12 * 8, 1 );
    Eigen::MatrixXd internalForceVectorResult = Eigen::MatrixXd::Zero( 12 * 8, 1 );

    errorOut error = overlapCoupling::formMicromorphicElementInternalForceVector( element, degreeOfFreedomValues,
                                                                                  cauchyStress, symmetricMicroStress, higherOrderStress,
                                                                                  &nodeIDToIndex, internalForceVectorResult );

    if ( error ){

        error->print( );
        results << "test_computeMicromorphicElementInternalForceVector & False\n";
        return 1;

    }

    if ( !internalForceVectorAnswer.isApprox( internalForceVectorResult, 1e-5 ) ){

        results << "test_computeMicromorphicElementInternalForceVector (test 1) & False\n";
        return 1;

    }

    results << "test_computeMicromorphicElementInternalForceVector & True\n";
    return 0;
}

int test_writeSparseMatrixToXDMF( std::ofstream &results ){
    /*!
     * Test writing a sparse matrix to XDMF file
     *
     * :param std::ofstream &results: The output file
     */

    shared_ptr< XdmfDomain > domain = XdmfDomain::New( );
    shared_ptr< XdmfUnstructuredGrid > grid = XdmfUnstructuredGrid::New( );

    std::string filename = "test_output_file";

    std::string xdmf_filename = filename + ".xdmf";
    std::string h5_filename = filename + ".h5";
    shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( h5_filename, true );
    shared_ptr< XdmfWriter > writer = XdmfWriter::New( xdmf_filename, heavyWriter );

    SparseMatrix A1( 3, 4 );
    std::vector< T > triplets;
    triplets.reserve( 7 );
    triplets.push_back( T( 0, 0, 1.0 ) );
    triplets.push_back( T( 0, 3, 1.5 ) );
    triplets.push_back( T( 2, 1, 7.0 ) );
    triplets.push_back( T( 1, 2, 5.0 ) );
    triplets.push_back( T( 0, 1, 2.0 ) );
    triplets.push_back( T( 1, 0, 1.6 ) );
    triplets.push_back( T( 0, 2, 3.0 ) );

    A1.setFromTriplets( triplets.begin( ), triplets.end( ) );

    domain->insert( grid );

    std::string matrixName = "A_MATRIX";
    errorOut error = overlapCoupling::writeSparseMatrixToXDMF( A1, matrixName, writer, domain, grid );

    if ( error ){

        error->print( );
        results << "test_writeSparseMatrixToXDMF & False\n";
        return 1;

    }

    shared_ptr< XdmfReader > reader = XdmfReader::New( );
    shared_ptr< XdmfDomain > _readDomain = shared_dynamic_cast< XdmfDomain >( reader->read( "test_output_file.xdmf" ) );
    shared_ptr< XdmfUnstructuredGrid > _readGrid = _readDomain->getUnstructuredGrid( 0 );

    SparseMatrix A1_result;
    error = overlapCoupling::readSparseMatrixFromXDMF( _readGrid, matrixName, A1_result );
    
    if ( error ){

        error->print( );
        results << "test_writeSparseMatrixToXDMF & False\n";
        return 1;

    }

    if ( !A1.isApprox( A1_result ) ){

        results << "test_writeSparseMatrixToXDMF (test 1) & False\n";
        return 1;

    }

    remove( h5_filename.c_str( ) );
    remove( xdmf_filename.c_str( ) );

    results << "test_writeSparseMatrixToXDMF & True\n";
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

    test_overlapCoupling_constructor( results );
    test_overlapCoupling_initializeCoupling( results );
    test_overlapCoupling_processIncrement( results );
    test_overlapCoupling_processLastIncrements( results );
    test_overlapCoupling_getReferenceFreeMicroDomainMasses( results );
    test_overlapCoupling_getReferenceGhostMicroDomainMasses( results );
    test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass( results );
    test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass( results );
//    test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions( results );
//    test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions( results );
    test_MADOutlierDetection( results );
    test_formMicromorphicElementMassMatrix( results );
    test_computeMicromorphicElementRequiredValues( results );
    test_computeMicromorphicElementInternalForceVector( results );
    test_writeSparseMatrixToXDMF( results );

    //Close the results file
    results.close();
}
