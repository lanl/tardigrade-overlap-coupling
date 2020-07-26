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

    results << "test_overlapCoupling_initializeCoupling & True\n";
    return 0;
}

int test_overlapCoupling_getReferenceFreeMicroDomainMasses( std::ofstream &results ){
    /*!
     * Test the extraction of the reference free micro-domain mass
     *
     * :param std::ofstream &results: The output file.
     */

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
    referenceFreeMicroDomainMassesAnswer *= 2000;

    const floatVector *referenceFreeMicroDomainMassesResult = oc.getReferenceFreeMicroDomainMasses( );

    if ( !vectorTools::fuzzyEquals( referenceFreeMicroDomainMassesAnswer, *referenceFreeMicroDomainMassesResult ) ){
        vectorTools::print( referenceFreeMicroDomainMassesAnswer );
        vectorTools::print( *referenceFreeMicroDomainMassesResult );
        std::cout << referenceFreeMicroDomainMassesResult->size( ) << "\n";
        results << "test_overlapCoupling_getReferenceFreeMicroDomainMasses (test 1) & False\n";
        return 1;
    }

    results << "test_overlapCoupling_getReferenceFreeMicroDomainMasses & True\n";
    return 0;
}

int test_overlapCoupling_getReferenceGhostMicroDomainMasses( std::ofstream &results ){
    /*!
     * Test the extraction of the reference ghost micro-domain mass
     *
     * :param std::ofstream &results: The output file.
     */

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
    referenceGhostMicroDomainMassesAnswer *= 2000;

    const floatVector *referenceGhostMicroDomainMassesResult = oc.getReferenceGhostMicroDomainMasses( );

    if ( !vectorTools::fuzzyEquals( referenceGhostMicroDomainMassesAnswer, *referenceGhostMicroDomainMassesResult ) ){
        results << "test_overlapCoupling_getReferenceGhostMicroDomainMasses (test 1) & False\n";
        return 1;
    }

    results << "test_overlapCoupling_getReferenceGhostMicroDomainMasses & True\n";
    return 0;
}
int test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass( std::ofstream &results ){
    /*!
     * Test the extraction of the reference free micro-domain centers of mass
     *
     * :param std::ofstream &results: The output file.
     */

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

    results << "test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass & True\n";
    return 0;
}

int test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass( std::ofstream &results ){
    /*!
     * Test the extraction of the reference ghost micro-domain centers of mass
     *
     * :param std::ofstream &results: The output file.
     */

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

    results << "test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass & True\n";
    return 0;
}

int test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions( std::ofstream &results ){
    /*!
     * Test of the retrieval of the reference free micro domain centers of mass shape functions
     */

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

    results << "test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions & True\n";
    return 0;

}

int test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions( std::ofstream &results ){
    /*!
     * Test of the retrieval of the reference ghost micro domain centers of mass shape functions
     */

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

    results << "test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions & True\n";
    return 0;

}

int test_overlapCoupling_processIncrement( std::ofstream &results ){
    /*!
     * Test the initialization of the coupling
     *
     * :param std::ofstream &results: The output file.
     */

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

    error = oc.processIncrement( 1, 1 );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_initializeCoupling & False\n";
        return 1;
    }

    const floatVector *projectedGhostMacroDisplacement = oc.getProjectedGhostMacroDisplacement( );
    const floatVector *projectedGhostMicroDisplacement = oc.getProjectedGhostMicroDisplacement( );

    if ( projectedGhostMacroDisplacement->size( ) == 0 ){
        results << "test_overlapCoupling_initializeCoupling (test 1) & False\n";
        return 1;
    }

    if ( projectedGhostMicroDisplacement->size( ) == 0 ){
        results << "test_overlapCoupling_initializeCoupling (test 2) & False\n";
        return 1;
    }

    results << "test_overlapCoupling_initializeCoupling & True\n";
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
    test_overlapCoupling_getReferenceFreeMicroDomainMasses( results );
    test_overlapCoupling_getReferenceGhostMicroDomainMasses( results );
    test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass( results );
    test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass( results );
//    test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions( results );
//    test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions( results );
    test_MADOutlierDetection( results );
    test_formMicromorphicElementMassMatrix( results );

    //Close the results file
    results.close();
}
