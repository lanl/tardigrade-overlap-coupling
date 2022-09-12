//!The test file for volumeReconstruction.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>
#define USE_EIGEN
#include<vector_tools.h>

#include<volumeReconstruction.h>

#define BOOST_TEST_MODULE test_volumeReconstruction
#include <boost/test/included/unit_test.hpp>

typedef volumeReconstruction::errorNode errorNode; //!Redefinition for the error node
typedef volumeReconstruction::errorOut errorOut; //!Redefinition for a pointer to the error node
typedef volumeReconstruction::floatType floatType; //!Define the float values type.
typedef volumeReconstruction::floatVector floatVector; //! Define a vector of floats
typedef volumeReconstruction::floatMatrix floatMatrix; //!Define a matrix of floats
typedef volumeReconstruction::intMatrix intMatrix; //!Define a matrix of ints
typedef volumeReconstruction::uIntType uIntType; //!Define the unsigned int type
typedef volumeReconstruction::uIntVector uIntVector; //!Define a vector of unsigned ints

BOOST_AUTO_TEST_CASE( testDualContouring_constructor ){
    /*!
     * Tests for the constructors of dualContouring
     *
     */

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    std::shared_ptr< volumeReconstruction::volumeReconstructionBase > yR
        = volumeReconstruction::volumeReconstructionBase( yf ).create( );

    BOOST_CHECK( !yR->getError( ) );

}

BOOST_AUTO_TEST_CASE( testDualContouring_loadPoints ){
    /*!
     * Test for loading the points into the object
     *
     */

    floatVector points = { 1, 2, 3, 4, 5, 6 };

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    std::unique_ptr< errorNode > error;
    error.reset( dc.loadPoints( &points ) );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( *dc.getPoints( ), points ) );

    floatVector points2 = { 1, 2, 3, 4, 5 };


    error.reset( dc.loadPoints( &points2 ) );

    BOOST_CHECK( error );

}

BOOST_AUTO_TEST_CASE( testDualContouring_loadFunction ){
    /*!
     * Test for loading the function into the object
     *
     */

    floatVector points = { 1, 2, 3, 4, 5, 6 };

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    std::unique_ptr< errorNode > error;
    error.reset( dc.loadPoints( &points ) );

    BOOST_CHECK( !error );

    floatVector function = { -1, 10 };

    error.reset( dc.loadFunction( &function ) );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( *dc.getFunction( ), function ) );

    floatVector function2 = { 2 };

    error.reset( dc.loadFunction( &function2 ) );
    BOOST_CHECK( error );

}

BOOST_AUTO_TEST_CASE( testDualContouring_getFunctionValue ){
    /*!
     * Test for getting a particular value of the function
     *
     */

    floatVector points = { 1, 2, 3, 4, 5, 6 };

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    errorOut error = dc.loadPoints( &points );

    BOOST_CHECK( !error );

    floatVector function = { -1, 10 };

    error = dc.loadFunction( &function );

    BOOST_CHECK( !error );

    floatType result = 0;
    error = dc.getFunctionValue( 0, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( function [ 0 ], result ) );

    volumeReconstruction::dualContouring dc2( yf );

    BOOST_CHECK( !dc2.getError( ) );

    error = dc2.loadPoints( &points );

    BOOST_CHECK( !error );

    result = 0;
    error = dc2.getFunctionValue( 0, result );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( 1., result ) );

}

BOOST_AUTO_TEST_CASE( testKDNode_constructor ){
    /*!
     * Test the KD tree creation
     *
     */

    floatVector points = { 2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2 };
    uIntVector  ownedIndices = { 0, 2, 4, 6, 8, 10 };
    unsigned int dim = 2;

    volumeReconstruction::KDNode tree( &points, ownedIndices, 0, dim );

}

BOOST_AUTO_TEST_CASE( testKDNode_getIndex ){
    /*!
     * Test getting the index of the KDNode
     *
     */

    floatVector points = { 2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2 };
    uIntVector  ownedIndices = { 0, 2, 4, 6, 8, 10 };
    unsigned int dim = 2;

    volumeReconstruction::KDNode tree( &points, ownedIndices, 0, dim );

    BOOST_CHECK( ( tree.getIndex( ), 10 ) );

}

BOOST_AUTO_TEST_CASE( testKDNode_getMinimumValueDimension ){
    /*!
     * Test getting the minimum of a dimension of the KD tree
     *
     */

    floatVector points = { 2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2 };
    uIntVector  ownedIndices = { 0, 2, 4, 6, 8, 10 };
    unsigned int dim = 2;

    volumeReconstruction::KDNode tree( &points, ownedIndices, 0, dim );

    floatType answer1 = 2;
    floatType answer2 = 1;

    floatType result = tree.getMinimumValueDimension( 0 );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer1 ) );

    result = tree.getMinimumValueDimension( 1 );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer2 ) );

    floatVector points2 = { 1, 10, 50, 50, 10, 30, 35, 90, 55, 1, 60, 80, 25, 40, 70, 70, 51, 75 };
    uIntVector ownedIndices2 = { 0, 2, 4, 6, 8, 10, 12, 14, 16 };

    volumeReconstruction::KDNode tree2( &points2, ownedIndices2, 0, dim );

    floatType answer3 = 1;
    floatType answer4 = 1;

    result = tree2.getMinimumValueDimension( 0 );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer3 ) );

    result = tree2.getMinimumValueDimension( 1 );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer4 ) );

    floatVector points3 =
        {
            0.88993453,  0.76248591, -0.93419017,  0.16825499, -0.64930721,
           -0.02736989,  0.3818072 , -0.56789234,  0.19636903, -0.01416032,
           -0.0726343 , -0.17899304,  0.58038292,  0.49538273, -0.16969211,
           -0.01199516,  0.83431932, -0.5944871 ,  0.68423515, -0.78823811,
           -0.30881808,  0.57273371, -0.3135023 , -0.0364649 ,  0.79618035,
           -0.36185261, -0.1701662 ,  0.4627748 , -0.22219916, -0.87057911,
           -0.34563548, -0.32495209,  0.90188535, -0.16083641, -0.67036656,
            0.04278752, -0.72313702, -0.62661373,  0.04437484, -0.62416924,
            0.3790917 , -0.99593042, -0.4127237 , -0.33142047,  0.91994519,
           -0.20164988,  0.48012869, -0.21532596,  0.92378245,  0.03305427,
            0.64008042,  0.36990853,  0.20280656,  0.82268547,  0.70533715,
            0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782
        };

    uIntVector ownedIndices3 = { 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57 };

    volumeReconstruction::KDNode tree3( &points3, ownedIndices3, 0, 3 );

    floatType answer5 = -0.72313702;
    floatType answer6 = -0.78823811;
    floatType answer7 = -0.99593042;

    BOOST_CHECK( vectorTools::fuzzyEquals( tree3.getMinimumValueDimension( 0 ), answer5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( tree3.getMinimumValueDimension( 1 ), answer6 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( tree3.getMinimumValueDimension( 2 ), answer7 ) );

}

BOOST_AUTO_TEST_CASE( testKDNode_getMaximumValueDimension ){
    /*!
     * Test getting the maximum of a dimension of the KD tree
     *
     */

    floatVector points = { 2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2 };
    uIntVector  ownedIndices = { 0, 2, 4, 6, 8, 10 };
    unsigned int dim = 2;

    volumeReconstruction::KDNode tree( &points, ownedIndices, 0, dim );

    floatType answer1 = 9;
    floatType answer2 = 7;

    floatType result = tree.getMaximumValueDimension( 0 );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer1 ) );

    result = tree.getMaximumValueDimension( 1 );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer2 ) );

    floatVector points2 = { 1, 10, 50, 50, 10, 30, 35, 90, 55, 1, 60, 80, 25, 40, 70, 70, 51, 75 };
    uIntVector ownedIndices2 = { 0, 2, 4, 6, 8, 10, 12, 14, 16 };

    volumeReconstruction::KDNode tree2( &points2, ownedIndices2, 0, dim );

    floatType answer3 = 70;
    floatType answer4 = 90;

    result = tree2.getMaximumValueDimension( 0 );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer3 ) );

    result = tree2.getMaximumValueDimension( 1 );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer4 ) );

    floatVector points3 =
        {
            0.88993453,  0.76248591, -0.93419017,  0.16825499, -0.64930721,
           -0.02736989,  0.3818072 , -0.56789234,  0.19636903, -0.01416032,
           -0.0726343 , -0.17899304,  0.58038292,  0.49538273, -0.16969211,
           -0.01199516,  0.83431932, -0.5944871 ,  0.68423515, -0.78823811,
           -0.30881808,  0.57273371, -0.3135023 , -0.0364649 ,  0.79618035,
           -0.36185261, -0.1701662 ,  0.4627748 , -0.22219916, -0.87057911,
           -0.34563548, -0.32495209,  0.90188535, -0.16083641, -0.67036656,
            0.04278752, -0.72313702, -0.62661373,  0.04437484, -0.62416924,
            0.3790917 , -0.99593042, -0.4127237 , -0.33142047,  0.91994519,
           -0.20164988,  0.48012869, -0.21532596,  0.92378245,  0.03305427,
            0.64008042,  0.36990853,  0.20280656,  0.82268547,  0.70533715,
            0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782
        };

    uIntVector ownedIndices3 = { 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57 };

    volumeReconstruction::KDNode tree3( &points3, ownedIndices3, 0, 3 );

    floatType answer5 = 0.92378245;
    floatType answer6 = 0.83431932;
    floatType answer7 = 0.91994519;

    BOOST_CHECK( vectorTools::fuzzyEquals( tree3.getMaximumValueDimension( 0 ), answer5 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( tree3.getMaximumValueDimension( 1 ), answer6 ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( tree3.getMaximumValueDimension( 2 ), answer7 ) );
}

BOOST_AUTO_TEST_CASE( testKDNode_getPointsInRange ){
    /*!
     * Get all of the points in a range
     *
     */

    floatVector points = { 2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2 };
    uIntVector  ownedIndices = { 0, 2, 4, 6, 8, 10 };
    unsigned int dim = 2;

    volumeReconstruction::KDNode tree( &points, ownedIndices, 0, dim );

    uIntVector answer = { 10, 2 };

    floatVector upperBound = { 7.5, 5.0 };
    floatVector lowerBound = { 3.5, 1.0 };

    uIntVector result;

    tree.getPointsInRange( upperBound, lowerBound, result );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

    floatVector points2 = { 1, 10, 50, 50, 10, 30, 35, 90, 55, 1, 60, 80, 25, 40, 70, 70, 51, 75 };
    uIntVector ownedIndices2 = { 0, 2, 4, 6, 8, 10, 12, 14, 16 };

    volumeReconstruction::KDNode tree2( &points2, ownedIndices2, 0, dim );

    answer = { 2, 12, 16 };

    upperBound = { 52, 76 };
    lowerBound = { 22,  5 };

    result.clear( );

    tree2.getPointsInRange( upperBound, lowerBound, result );

    BOOST_CHECK( vectorTools::fuzzyEquals( result, answer ) );

    floatVector points3 =
        {
            0.88993453,  0.76248591, -0.93419017,  0.16825499, -0.64930721,
           -0.02736989,  0.3818072 , -0.56789234,  0.19636903, -0.01416032,
           -0.0726343 , -0.17899304,  0.58038292,  0.49538273, -0.16969211,
           -0.01199516,  0.83431932, -0.5944871 ,  0.68423515, -0.78823811,
           -0.30881808,  0.57273371, -0.3135023 , -0.0364649 ,  0.79618035,
           -0.36185261, -0.1701662 ,  0.4627748 , -0.22219916, -0.87057911,
           -0.34563548, -0.32495209,  0.90188535, -0.16083641, -0.67036656,
            0.04278752, -0.72313702, -0.62661373,  0.04437484, -0.62416924,
            0.3790917 , -0.99593042, -0.4127237 , -0.33142047,  0.91994519,
           -0.20164988,  0.48012869, -0.21532596,  0.92378245,  0.03305427,
            0.64008042,  0.36990853,  0.20280656,  0.82268547,  0.70533715,
            0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782
        };

    uIntVector ownedIndices3 = { 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57 };

    volumeReconstruction::KDNode tree3( &points3, ownedIndices3, 0, 3 );

    answer = { 21, 9, 3 };

    upperBound = {  0.75,  0.35,  0.0 };
    lowerBound = { -0.70, -0.70, -0.7 };

    result.clear( );

    tree3.getPointsInRange( upperBound, lowerBound, result );

    for ( auto a = answer.begin( ); a != answer.end( ); a++ ){

        BOOST_CHECK( std::any_of( result.begin( ), result.end( ),
                     [&]( uIntType r ){ return vectorTools::fuzzyEquals( r, *a ); } ) );

    }

}

BOOST_AUTO_TEST_CASE( testKDNode_getPointsWithinRadiusOfOrigin ){
    /*!
     * Get all of the points within a given radius of the origin
     *
     */

    floatVector points = { 2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2 };
    uIntVector  ownedIndices = { 0, 2, 4, 6, 8, 10 };
    unsigned int dim = 2;
    floatVector origin = { 4.5, 2.1 };
    floatType radius = 3.0;

    volumeReconstruction::KDNode tree( &points, ownedIndices, 0, dim );

    uIntVector answer = { 0, 2, 10 };

    uIntVector result;

    tree.getPointsWithinRadiusOfOrigin( origin, radius, result );

    for ( auto a = answer.begin( ); a != answer.end( ); a++ ){

        BOOST_CHECK( std::any_of( result.begin( ), result.end( ),
                           [&]( uIntType r ){ return vectorTools::fuzzyEquals( r, *a ); } ) );

    }

    floatVector points2 = { 1, 10, 50, 50, 10, 30, 35, 90, 55, 1, 60, 80, 25, 40, 70, 70, 51, 75 };
    uIntVector ownedIndices2 = { 0, 2, 4, 6, 8, 10, 12, 14, 16 };

    volumeReconstruction::KDNode tree2( &points2, ownedIndices2, 0, dim );

    origin = { 24, 65 };

    radius = 32;

    answer = { 2, 6, 12, 16 };

    result.clear( );

    tree2.getPointsWithinRadiusOfOrigin( origin, radius, result );

    for ( auto a = answer.begin( ); a != answer.end( ); a++ ){

        BOOST_CHECK( std::any_of( result.begin( ), result.end( ),
                                  [&]( uIntType r ){ return vectorTools::fuzzyEquals( r, *a ); } ) );

    }

    floatVector points3 =
        {
            0.88993453,  0.76248591, -0.93419017,  0.16825499, -0.64930721,
           -0.02736989,  0.3818072 , -0.56789234,  0.19636903, -0.01416032,
           -0.0726343 , -0.17899304,  0.58038292,  0.49538273, -0.16969211,
           -0.01199516,  0.83431932, -0.5944871 ,  0.68423515, -0.78823811,
           -0.30881808,  0.57273371, -0.3135023 , -0.0364649 ,  0.79618035,
           -0.36185261, -0.1701662 ,  0.4627748 , -0.22219916, -0.87057911,
           -0.34563548, -0.32495209,  0.90188535, -0.16083641, -0.67036656,
            0.04278752, -0.72313702, -0.62661373,  0.04437484, -0.62416924,
            0.3790917 , -0.99593042, -0.4127237 , -0.33142047,  0.91994519,
           -0.20164988,  0.48012869, -0.21532596,  0.92378245,  0.03305427,
            0.64008042,  0.36990853,  0.20280656,  0.82268547,  0.70533715,
            0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782,
           -0.75918987, -0.19216117,  0.82956824,  0.18528101, -0.52393903,
            0.70272336, -0.45594571,  0.18255914, -0.03938857, -0.78883385,
            0.85656866,  0.75744535, -0.77273419, -0.14567324,  0.1562869 ,
            0.67723809, -0.63066704, -0.10145066, -0.94238038,  0.58534946,
            0.29427379,  0.87440634,  0.81035162,  0.0477631 ,  0.82727563,
            0.74320746,  0.22778598,  0.58000421,  0.57222056,  0.09287177,
            0.47295155,  0.72360676, -0.41622328,  0.82540552,  0.32944547,
           -0.29634423,  0.19641306,  0.41306438,  0.36451089, -0.31113366,
           -0.93875619, -0.53893408, -0.36899228, -0.50751604,  0.42847307,
           -0.72351422, -0.1097771 , -0.51914444, -0.8399524 , -0.37058809,
           -0.53075988, -0.75968017, -0.13391175, -0.42047612,  0.93090777,
            0.87886273, -0.90105205,  0.37213616, -0.98616276,  0.92134844,
            0.98173517, -0.45791212, -0.41495029, -0.12509758, -0.25764766,
           -0.39942542, -0.07300332, -0.50951107,  0.47335924,  0.9226246 ,
           -0.8225396 ,  0.83189879,  0.97164758,  0.07520523,  0.96401949,
            0.09115673, -0.15473381, -0.12517667,  0.53923746,  0.17375421,
           -0.90864801,  0.13164519, -0.81746094, -0.04091139,  0.70761535,
            0.27949176, -0.72591527,  0.62191923,  0.58579313, -0.63652057,
           -0.82922138, -0.20341235, -0.36896014,  0.93666234,  0.99606496,
            0.02604581,  0.682005  , -0.83573812, -0.94639137, -0.93744807,
            0.41838408,  0.53512206,  0.73833179,  0.43511637,  0.07707317,
            0.94551536,  0.06136697, -0.88529597, -0.37223707, -0.13058219,
            0.19189113,  0.91949737,  0.41992734,  0.73392705, -0.81952681,
            0.97131812, -0.90971185, -0.83583963, -0.39223689, -0.2724691 ,
           -0.88024076, -0.2430722 ,  0.83700368, -0.35641362,  0.60074684,
            0.47984156,  0.7821192 ,  0.92905141, -0.66941659,  0.92052141,
            0.66587504,  0.64540513,  0.32427648,  0.20626864,  0.98728317,
           -0.46970216, -0.9715591 ,  0.52876951,  0.00966212,  0.4786293 ,
           -0.08921145,  0.93996844, -0.86556617, -0.43238781, -0.67206028,
           -0.05996167, -0.23409636, -0.79355771, -0.78246022, -0.52978658,
           -0.81887169,  0.32057883, -0.1984804 ,  0.20959359,  0.89333809,
           -0.26993128, -0.41085675, -0.2522169 ,  0.44736448,  0.62315547,
           -0.7794221 , -0.18013406,  0.26102072, -0.18430868,  0.99247162,
           -0.6383439 , -0.59488566,  0.36539787,  0.90600975, -0.81800276,
           -0.73782801,  0.79265261,  0.28872337, -0.16722855, -0.6239667 ,
           -0.45177053, -0.196289  ,  0.84749471, -0.41257737,  0.41538694,
           -0.57353971, -0.0891511 , -0.66645294,  0.6060595 ,  0.06013818,
            0.43417399,  0.85628154, -0.25946277, -0.17124759, -0.27269205,
           -0.65624534,  0.76907186, -0.15572072, -0.85328073,  0.64179789,
            0.7853888 , -0.73223644,  0.67797289, -0.66645228,  0.01685972,
           -0.27934479, -0.95963138, -0.78014927,  0.62397642, -0.41740632,
            0.82117605, -0.02459754,  0.15318692, -0.5935911 , -0.70041395,
           -0.79377562,  0.93858127, -0.02915152,  0.75699326,  0.39613993,
            0.82603044,  0.96104565,  0.2936643 ,  0.48070712,  0.83404655,
            0.14211884,  0.00691666,  0.80661045,  0.88463205, -0.88978482,
           -0.28513985, -0.30738596, -0.20413492,  0.21292329, -0.26989904,
           -0.85370154,  0.72208774, -0.55094191,  0.37445557,  0.54146363,
            0.57375785,  0.3679136 ,  0.65563202, -0.06865393, -0.59308587
        };

    uIntVector ownedIndices3( points3.size( ) / 3 );
    for ( unsigned int i = 0; i < points3.size( ) / 3; i++ ){

        ownedIndices3[ i ] = 3 * i;

    }

    volumeReconstruction::KDNode tree3( &points3, ownedIndices3, 0, 3 );

    answer = { 63, 117, 252 };

    origin = { 0.16131353, -0.68809742,  0.91193468 };
    radius = 0.45;

    result.clear( );

    tree3.getPointsWithinRadiusOfOrigin( origin, radius, result );

    for ( auto a = answer.begin( ); a != answer.end( ); a++ ){

        BOOST_CHECK( std::any_of( result.begin( ), result.end( ),
                                  [&]( uIntType r ){ return vectorTools::fuzzyEquals( r, *a ); } ) );

    }

    floatVector points4 =
        {
            -0.960002, -0.986542, -0.996311,
            -0.96004, -0.98658, -0.49883,
            -0.96004, -0.98658, 0.495371,
            -0.960002, -0.986542, 0.992852,
            -0.96004, -0.490606, -0.996349,
            -0.960117, -0.490606, -0.49883,
            -0.960117, -0.490606, 0.495371,
            -0.96004, -0.490606, 0.99289,
            -0.96004, 0.500508, -0.996349,
            -0.960117, 0.500508, -0.49883,
            -0.960117, 0.500508, 0.495371,
            -0.96004, 0.500508, 0.99289,
            -0.960002, 0.996444, -0.996311,
            -0.96004, 0.996482, -0.49883,
            -0.96004, 0.996482, 0.495371,
            -0.960002, 0.996444, 0.992852,
            -0.47429, -0.98658, -0.996349,
            -0.47429, -0.986658, -0.49883,
            -0.47429, -0.986658, 0.495371,
            -0.47429, -0.98658, 0.99289,
            -0.47429, -0.490606, -0.996428,
            -0.47429, -0.490606, 0.992969,
            -0.47429, 0.500508, -0.996428,
            -0.359348, 0.616812, 0.61188,
            -0.474518, 0.50028, 0.992869,
            -0.47429, 0.996482, -0.996349,
            -0.47429, 0.996561, -0.49883,
            -0.474517, 0.996461, 0.495143,
            -0.474541, 0.996424, 0.992832,
            0.496394, -0.98658, -0.996349,
            0.496394, -0.986658, -0.49883,
            0.496394, -0.986658, 0.495371,
            0.496394, -0.98658, 0.99289,
            0.496394, -0.490606, -0.996428,
            0.496394, -0.490606, 0.992969,
            0.496394, 0.500508, -0.996428,
            0.381452, 0.616812, 0.61188,
            0.496622, 0.50028, 0.992869,
            0.496394, 0.996482, -0.996349,
            0.496394, 0.996561, -0.49883,
            0.496621, 0.996461, 0.495143,
            0.496645, 0.996424, 0.992832,
            0.982106, -0.986542, -0.996311,
            0.982144, -0.98658, -0.49883,
            0.982144, -0.98658, 0.495371,
            0.982106, -0.986542, 0.992852,
            0.982144, -0.490606, -0.996349,
            0.982221, -0.490606, -0.49883,
            0.982221, -0.490606, 0.495371,
            0.982144, -0.490606, 0.99289,
            0.982144, 0.500508, -0.996349,
            0.982221, 0.500508, -0.49883,
            0.982221, 0.500508, 0.495371,
            0.982144, 0.500508, 0.99289,
            0.982106, 0.996444, -0.996311,
            0.982144, 0.996482, -0.49883,
            0.982144, 0.996482, 0.495371,
            0.982106, 0.996444, 0.992852
        };

    uIntVector ownedIndices4( points4.size( ) / 3 );

    for ( unsigned int i = 0; i < ( points4.size( ) / 3 ); i++ ){

        ownedIndices4[ i ] = 3 * i;

    }

    volumeReconstruction::KDNode tree4( &points4, ownedIndices4, 0, 3 );

    answer = { 0,  3, 12, 15, 48, 51, 60 };

    origin = { -0.960002, -0.986542, -0.996311 };
    radius = 1.0;

    result.clear( );

    tree4.getPointsWithinRadiusOfOrigin( origin, radius, result );

    for ( auto a = answer.begin( ); a != answer.end( ); a++ ){

        BOOST_CHECK( std::any_of( result.begin( ), result.end( ),
                                  [&]( uIntType r ){ return vectorTools::fuzzyEquals( r, *a ); } ) );

    }

}

BOOST_AUTO_TEST_CASE( testDualContouring_evaluate ){
    /*!
     * Test the dualContouring evaluate function. This prepares the 
     * volume reconstruction object to perform the other functions
     *
     */

    floatVector points =
        {
             0.88993453,  0.76248591, -0.93419017,  0.16825499, -0.64930721,
        -0.02736989,  0.3818072 , -0.56789234,  0.19636903, -0.01416032,
        -0.0726343 , -0.17899304,  0.58038292,  0.49538273, -0.16969211,
        -0.01199516,  0.83431932, -0.5944871 ,  0.68423515, -0.78823811,
        -0.30881808,  0.57273371, -0.3135023 , -0.0364649 ,  0.79618035,
        -0.36185261, -0.1701662 ,  0.4627748 , -0.22219916, -0.87057911,
        -0.34563548, -0.32495209,  0.90188535, -0.16083641, -0.67036656,
         0.04278752, -0.72313702, -0.62661373,  0.04437484, -0.62416924,
         0.3790917 , -0.99593042, -0.4127237 , -0.33142047,  0.91994519,
        -0.20164988,  0.48012869, -0.21532596,  0.92378245,  0.03305427,
         0.64008042,  0.36990853,  0.20280656,  0.82268547,  0.70533715,
         0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782,
        -0.75918987, -0.19216117,  0.82956824,  0.18528101, -0.52393903,
         0.70272336, -0.45594571,  0.18255914, -0.03938857, -0.78883385,
         0.85656866,  0.75744535, -0.77273419, -0.14567324,  0.1562869 ,
         0.67723809, -0.63066704, -0.10145066, -0.94238038,  0.58534946,
         0.29427379,  0.87440634,  0.81035162,  0.0477631 ,  0.82727563,
         0.74320746,  0.22778598,  0.58000421,  0.57222056,  0.09287177,
         0.47295155,  0.72360676, -0.41622328,  0.82540552,  0.32944547,
        -0.29634423,  0.19641306,  0.41306438,  0.36451089, -0.31113366,
        -0.93875619, -0.53893408, -0.36899228, -0.50751604,  0.42847307,
        -0.72351422, -0.1097771 , -0.51914444, -0.8399524 , -0.37058809,
        -0.53075988, -0.75968017, -0.13391175, -0.42047612,  0.93090777,
         0.87886273, -0.90105205,  0.37213616, -0.98616276,  0.92134844,
         0.98173517, -0.45791212, -0.41495029, -0.12509758, -0.25764766,
        -0.39942542, -0.07300332, -0.50951107,  0.47335924,  0.9226246 ,
        -0.8225396 ,  0.83189879,  0.97164758,  0.07520523,  0.96401949,
         0.09115673, -0.15473381, -0.12517667,  0.53923746,  0.17375421,
        -0.90864801,  0.13164519, -0.81746094, -0.04091139,  0.70761535,
         0.27949176, -0.72591527,  0.62191923,  0.58579313, -0.63652057,
        -0.82922138, -0.20341235, -0.36896014,  0.93666234,  0.99606496,
         0.02604581,  0.682005  , -0.83573812, -0.94639137, -0.93744807,
         0.41838408,  0.53512206,  0.73833179,  0.43511637,  0.07707317,
         0.94551536,  0.06136697, -0.88529597, -0.37223707, -0.13058219,
         0.19189113,  0.91949737,  0.41992734,  0.73392705, -0.81952681,
         0.97131812, -0.90971185, -0.83583963, -0.39223689, -0.2724691 ,
        -0.88024076, -0.2430722 ,  0.83700368, -0.35641362,  0.60074684,
         0.47984156,  0.7821192 ,  0.92905141, -0.66941659,  0.92052141,
         0.66587504,  0.64540513,  0.32427648,  0.20626864,  0.98728317,
        -0.46970216, -0.9715591 ,  0.52876951,  0.00966212,  0.4786293 ,
        -0.08921145,  0.93996844, -0.86556617, -0.43238781, -0.67206028,
        -0.05996167, -0.23409636, -0.79355771, -0.78246022, -0.52978658,
        -0.81887169,  0.32057883, -0.1984804 ,  0.20959359,  0.89333809,
        -0.26993128, -0.41085675, -0.2522169 ,  0.44736448,  0.62315547,
        -0.7794221 , -0.18013406,  0.26102072, -0.18430868,  0.99247162,
        -0.6383439 , -0.59488566,  0.36539787,  0.90600975, -0.81800276,
        -0.73782801,  0.79265261,  0.28872337, -0.16722855, -0.6239667 ,
        -0.45177053, -0.196289  ,  0.84749471, -0.41257737,  0.41538694,
        -0.57353971, -0.0891511 , -0.66645294,  0.6060595 ,  0.06013818,
         0.43417399,  0.85628154, -0.25946277, -0.17124759, -0.27269205,
        -0.65624534,  0.76907186, -0.15572072, -0.85328073,  0.64179789,
         0.7853888 , -0.73223644,  0.67797289, -0.66645228,  0.01685972,
        -0.27934479, -0.95963138, -0.78014927,  0.62397642, -0.41740632,
         0.82117605, -0.02459754,  0.15318692, -0.5935911 , -0.70041395,
        -0.79377562,  0.93858127, -0.02915152,  0.75699326,  0.39613993,
         0.82603044,  0.96104565,  0.2936643 ,  0.48070712,  0.83404655,
         0.14211884,  0.00691666,  0.80661045,  0.88463205, -0.88978482,
        -0.28513985, -0.30738596, -0.20413492,  0.21292329, -0.26989904,
        -0.85370154,  0.72208774, -0.55094191,  0.37445557,  0.54146363,
         0.57375785,  0.3679136 ,  0.65563202, -0.06865393, -0.59308587
        };

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    errorOut error = dc.loadPoints( &points );

    BOOST_CHECK( !error );

    error = dc.evaluate( );

    BOOST_CHECK( !error );

    const floatVector upperBoundsAnswer = { 0.98173517, 0.99606496, 0.99247162 };

    const floatVector lowerBoundsAnswer = { -0.95963138, -0.98616276, -0.99593042 };

    const floatType medianNeighborhoodDistanceAnswer = 0.398116;

    const floatVector* upperBoundsResult = dc.getUpperBounds( );
    const floatVector* lowerBoundsResult = dc.getLowerBounds( );
    const floatType*   medianNeighborhoodDistanceResult = dc.getMedianNeighborhoodDistance( );

    BOOST_CHECK( vectorTools::fuzzyEquals( *upperBoundsResult, upperBoundsAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *lowerBoundsResult, lowerBoundsAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *medianNeighborhoodDistanceResult, medianNeighborhoodDistanceAnswer ) );

}

BOOST_AUTO_TEST_CASE( testDualContouringInternalPointResidual ){
    /*!
     * Test the dual contouring internal point residual
     *
     */

    floatVector X =
        { 
            0.85131874, 0.6459241 , 0.40004273, 0.32050015, 0.28067341,
            0.03255095, 0.84674781, 0.74372308, 0.21725097, 0.15472211,
            0.50591758, 0.11292911, 0.53883793, 0.28100951, 0.27188129
        };

    floatMatrix floatArgs =
        {
            { 0.73032719,  1.35644613,  0.95486593 },
            { 0.23311812, -0.34546654, -0.33838301 },
            { 0.13504783,  0.5675009 ,  0.01237095 },
            { 0.902437  ,  0.40287295,  0.21996613 },
            { 0.2710092 ,  0.93281803,  0.96473859 },
            { 0.95151858,  0.77561579,  0.17862172 },
            { 0.71134891,  0.12861138,  0.42750324 },
            { 0.13878509,  0.47367194,  0.65040988 }
        };

    intMatrix intArgs =
        {
            { 3, 3 },
        };

    floatVector residualAnswer =
        {
             2.93458538,  2.60123634,  0.94198874,  0.09917692,  0.28399522,
             0.0073519 , -0.91251968, -0.41798652, -0.11813294, -0.22371189,
             0.63174447,  0.55376364, -0.09878124,  0.43826662,  0.69122776            
        };

    floatVector jacobianAnswerVec =
        {
            4.43066618,  0.895239  ,  0.56433305,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  1.        ,
            0.        ,  0.        , -1.        , -0.        , -0.        ,
            0.895239  ,  3.84248584,  0.50160452,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            1.        ,  0.        , -0.        , -1.        , -0.        ,
            0.56433305,  0.50160452,  3.63769775,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  1.        , -0.        , -0.        , -1.        ,
            0.        ,  0.        ,  0.        ,  0.30944423,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.64100031,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  1.01183516,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.56134682,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.22585822,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.06510191,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        , -1.07767587, -0.        , -0.        ,  0.        ,
            0.        ,  0.        , -1.69349563, -0.        , -0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        , -0.        , -0.56201902, -0.        ,  0.        ,
            0.        ,  0.        , -0.        , -1.48744617, -0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        , -0.        , -0.        , -0.54376258,  0.        ,
            0.        ,  0.        , -0.        , -0.        , -0.43450193,
           -1.        , -0.        , -0.        , -0.64100031, -0.        ,
           -0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
           -0.        , -1.        , -0.        , -0.        , -0.56134682,
           -0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
           -0.        , -0.        , -1.        , -0.        , -0.        ,
           -0.06510191,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        , -1.69349563, -0.        , -0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  1.        ,  0.        ,  0.        ,  0.        ,
            0.        , -0.        , -1.48744617, -0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  1.        ,  0.        ,  0.        ,
            0.        , -0.        , -0.        , -0.43450193,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.
        };

    floatVector residualResult;
    floatMatrix jacobian;
    floatMatrix floatOuts;
    intMatrix   intOuts;

    errorOut error = volumeReconstruction::dualContouringInternalPointResidual( X, floatArgs, intArgs, residualResult, jacobian,
                                                                                floatOuts, intOuts );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( residualResult, residualAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( vectorTools::appendVectors( jacobian ), jacobianAnswerVec ) );

}

BOOST_AUTO_TEST_CASE( testDualContouring_performVolumeIntegration ){
    /*!
     * Test volume integration over the reconstructed domain
     *
     */

    floatVector points =
        {
             0.88993453,  0.76248591, -0.93419017,  0.16825499, -0.64930721,
            -0.02736989,  0.3818072 , -0.56789234,  0.19636903, -0.01416032,
            -0.0726343 , -0.17899304,  0.58038292,  0.49538273, -0.16969211,
            -0.01199516,  0.83431932, -0.5944871 ,  0.68423515, -0.78823811,
            -0.30881808,  0.57273371, -0.3135023 , -0.0364649 ,  0.79618035,
            -0.36185261, -0.1701662 ,  0.4627748 , -0.22219916, -0.87057911,
            -0.34563548, -0.32495209,  0.90188535, -0.16083641, -0.67036656,
             0.04278752, -0.72313702, -0.62661373,  0.04437484, -0.62416924,
             0.3790917 , -0.99593042, -0.4127237 , -0.33142047,  0.91994519,
            -0.20164988,  0.48012869, -0.21532596,  0.92378245,  0.03305427,
             0.64008042,  0.36990853,  0.20280656,  0.82268547,  0.70533715,
             0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782,
            -0.75918987, -0.19216117,  0.82956824,  0.18528101, -0.52393903,
             0.70272336, -0.45594571,  0.18255914, -0.03938857, -0.78883385,
             0.85656866,  0.75744535, -0.77273419, -0.14567324,  0.1562869 ,
             0.67723809, -0.63066704, -0.10145066, -0.94238038,  0.58534946,
             0.29427379,  0.87440634,  0.81035162,  0.0477631 ,  0.82727563,
             0.74320746,  0.22778598,  0.58000421,  0.57222056,  0.09287177,
             0.47295155,  0.72360676, -0.41622328,  0.82540552,  0.32944547,
            -0.29634423,  0.19641306,  0.41306438,  0.36451089, -0.31113366,
            -0.93875619, -0.53893408, -0.36899228, -0.50751604,  0.42847307,
            -0.72351422, -0.1097771 , -0.51914444, -0.8399524 , -0.37058809,
            -0.53075988, -0.75968017, -0.13391175, -0.42047612,  0.93090777,
             0.87886273, -0.90105205,  0.37213616, -0.98616276,  0.92134844,
             0.98173517, -0.45791212, -0.41495029, -0.12509758, -0.25764766,
            -0.39942542, -0.07300332, -0.50951107,  0.47335924,  0.9226246 ,
            -0.8225396 ,  0.83189879,  0.97164758,  0.07520523,  0.96401949,
             0.09115673, -0.15473381, -0.12517667,  0.53923746,  0.17375421,
            -0.90864801,  0.13164519, -0.81746094, -0.04091139,  0.70761535,
             0.27949176, -0.72591527,  0.62191923,  0.58579313, -0.63652057,
            -0.82922138, -0.20341235, -0.36896014,  0.93666234,  0.99606496,
             0.02604581,  0.682005  , -0.83573812, -0.94639137, -0.93744807,
             0.41838408,  0.53512206,  0.73833179,  0.43511637,  0.07707317,
             0.94551536,  0.06136697, -0.88529597, -0.37223707, -0.13058219,
             0.19189113,  0.91949737,  0.41992734,  0.73392705, -0.81952681,
             0.97131812, -0.90971185, -0.83583963, -0.39223689, -0.2724691 ,
            -0.88024076, -0.2430722 ,  0.83700368, -0.35641362,  0.60074684,
             0.47984156,  0.7821192 ,  0.92905141, -0.66941659,  0.92052141,
             0.66587504,  0.64540513,  0.32427648,  0.20626864,  0.98728317,
            -0.46970216, -0.9715591 ,  0.52876951,  0.00966212,  0.4786293 ,
            -0.08921145,  0.93996844, -0.86556617, -0.43238781, -0.67206028,
            -0.05996167, -0.23409636, -0.79355771, -0.78246022, -0.52978658,
            -0.81887169,  0.32057883, -0.1984804 ,  0.20959359,  0.89333809,
            -0.26993128, -0.41085675, -0.2522169 ,  0.44736448,  0.62315547,
            -0.7794221 , -0.18013406,  0.26102072, -0.18430868,  0.99247162,
            -0.6383439 , -0.59488566,  0.36539787,  0.90600975, -0.81800276,
            -0.73782801,  0.79265261,  0.28872337, -0.16722855, -0.6239667 ,
            -0.45177053, -0.196289  ,  0.84749471, -0.41257737,  0.41538694,
            -0.57353971, -0.0891511 , -0.66645294,  0.6060595 ,  0.06013818,
             0.43417399,  0.85628154, -0.25946277, -0.17124759, -0.27269205,
            -0.65624534,  0.76907186, -0.15572072, -0.85328073,  0.64179789,
             0.7853888 , -0.73223644,  0.67797289, -0.66645228,  0.01685972,
            -0.27934479, -0.95963138, -0.78014927,  0.62397642, -0.41740632,
             0.82117605, -0.02459754,  0.15318692, -0.5935911 , -0.70041395,
            -0.79377562,  0.93858127, -0.02915152,  0.75699326,  0.39613993,
             0.82603044,  0.96104565,  0.2936643 ,  0.48070712,  0.83404655,
             0.14211884,  0.00691666,  0.80661045,  0.88463205, -0.88978482,
            -0.28513985, -0.30738596, -0.20413492,  0.21292329, -0.26989904,
            -0.85370154,  0.72208774, -0.55094191,  0.37445557,  0.54146363,
             0.57375785,  0.3679136 ,  0.65563202, -0.06865393, -0.59308587
        };

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    errorOut error = dc.loadPoints( &points );

    BOOST_CHECK( !error );

    floatVector functionValues( points.size( ) );
    for ( unsigned int i = 0; i < points.size( ); i+=3 ){

        functionValues[ i + 0 ] = 1;
        functionValues[ i + 1 ] = 2;
        functionValues[ i + 2 ] = 3;

    }

    floatVector integratedVolumeResult;
    floatVector integratedVolumeAnswer = { 6.52002 , 13.04003 , 19.560051 };

    error = dc.performVolumeIntegration( functionValues, 3, integratedVolumeResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( integratedVolumeResult, integratedVolumeAnswer ) );

}

BOOST_AUTO_TEST_CASE( testDualContouring_performRelativePositionVolumeIntegration ){
    /*!
     * Test volume integration over the reconstructed domain utilizing the
     * relative position to compute a dyadic product.
     *
     */

    floatVector points =
        {
             0.88993453,  0.76248591, -0.93419017,  0.16825499, -0.64930721,
            -0.02736989,  0.3818072 , -0.56789234,  0.19636903, -0.01416032,
            -0.0726343 , -0.17899304,  0.58038292,  0.49538273, -0.16969211,
            -0.01199516,  0.83431932, -0.5944871 ,  0.68423515, -0.78823811,
            -0.30881808,  0.57273371, -0.3135023 , -0.0364649 ,  0.79618035,
            -0.36185261, -0.1701662 ,  0.4627748 , -0.22219916, -0.87057911,
            -0.34563548, -0.32495209,  0.90188535, -0.16083641, -0.67036656,
             0.04278752, -0.72313702, -0.62661373,  0.04437484, -0.62416924,
             0.3790917 , -0.99593042, -0.4127237 , -0.33142047,  0.91994519,
            -0.20164988,  0.48012869, -0.21532596,  0.92378245,  0.03305427,
             0.64008042,  0.36990853,  0.20280656,  0.82268547,  0.70533715,
             0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782,
            -0.75918987, -0.19216117,  0.82956824,  0.18528101, -0.52393903,
             0.70272336, -0.45594571,  0.18255914, -0.03938857, -0.78883385,
             0.85656866,  0.75744535, -0.77273419, -0.14567324,  0.1562869 ,
             0.67723809, -0.63066704, -0.10145066, -0.94238038,  0.58534946,
             0.29427379,  0.87440634,  0.81035162,  0.0477631 ,  0.82727563,
             0.74320746,  0.22778598,  0.58000421,  0.57222056,  0.09287177,
             0.47295155,  0.72360676, -0.41622328,  0.82540552,  0.32944547,
            -0.29634423,  0.19641306,  0.41306438,  0.36451089, -0.31113366,
            -0.93875619, -0.53893408, -0.36899228, -0.50751604,  0.42847307,
            -0.72351422, -0.1097771 , -0.51914444, -0.8399524 , -0.37058809,
            -0.53075988, -0.75968017, -0.13391175, -0.42047612,  0.93090777,
             0.87886273, -0.90105205,  0.37213616, -0.98616276,  0.92134844,
             0.98173517, -0.45791212, -0.41495029, -0.12509758, -0.25764766,
            -0.39942542, -0.07300332, -0.50951107,  0.47335924,  0.9226246 ,
            -0.8225396 ,  0.83189879,  0.97164758,  0.07520523,  0.96401949,
             0.09115673, -0.15473381, -0.12517667,  0.53923746,  0.17375421,
            -0.90864801,  0.13164519, -0.81746094, -0.04091139,  0.70761535,
             0.27949176, -0.72591527,  0.62191923,  0.58579313, -0.63652057,
            -0.82922138, -0.20341235, -0.36896014,  0.93666234,  0.99606496,
             0.02604581,  0.682005  , -0.83573812, -0.94639137, -0.93744807,
             0.41838408,  0.53512206,  0.73833179,  0.43511637,  0.07707317,
             0.94551536,  0.06136697, -0.88529597, -0.37223707, -0.13058219,
             0.19189113,  0.91949737,  0.41992734,  0.73392705, -0.81952681,
             0.97131812, -0.90971185, -0.83583963, -0.39223689, -0.2724691 ,
            -0.88024076, -0.2430722 ,  0.83700368, -0.35641362,  0.60074684,
             0.47984156,  0.7821192 ,  0.92905141, -0.66941659,  0.92052141,
             0.66587504,  0.64540513,  0.32427648,  0.20626864,  0.98728317,
            -0.46970216, -0.9715591 ,  0.52876951,  0.00966212,  0.4786293 ,
            -0.08921145,  0.93996844, -0.86556617, -0.43238781, -0.67206028,
            -0.05996167, -0.23409636, -0.79355771, -0.78246022, -0.52978658,
            -0.81887169,  0.32057883, -0.1984804 ,  0.20959359,  0.89333809,
            -0.26993128, -0.41085675, -0.2522169 ,  0.44736448,  0.62315547,
            -0.7794221 , -0.18013406,  0.26102072, -0.18430868,  0.99247162,
            -0.6383439 , -0.59488566,  0.36539787,  0.90600975, -0.81800276,
            -0.73782801,  0.79265261,  0.28872337, -0.16722855, -0.6239667 ,
            -0.45177053, -0.196289  ,  0.84749471, -0.41257737,  0.41538694,
            -0.57353971, -0.0891511 , -0.66645294,  0.6060595 ,  0.06013818,
             0.43417399,  0.85628154, -0.25946277, -0.17124759, -0.27269205,
            -0.65624534,  0.76907186, -0.15572072, -0.85328073,  0.64179789,
             0.7853888 , -0.73223644,  0.67797289, -0.66645228,  0.01685972,
            -0.27934479, -0.95963138, -0.78014927,  0.62397642, -0.41740632,
             0.82117605, -0.02459754,  0.15318692, -0.5935911 , -0.70041395,
            -0.79377562,  0.93858127, -0.02915152,  0.75699326,  0.39613993,
             0.82603044,  0.96104565,  0.2936643 ,  0.48070712,  0.83404655,
             0.14211884,  0.00691666,  0.80661045,  0.88463205, -0.88978482,
            -0.28513985, -0.30738596, -0.20413492,  0.21292329, -0.26989904,
            -0.85370154,  0.72208774, -0.55094191,  0.37445557,  0.54146363,
             0.57375785,  0.3679136 ,  0.65563202, -0.06865393, -0.59308587
        };

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    errorOut error = dc.loadPoints( &points );

    BOOST_CHECK( !error );

    floatVector functionValues( points.size( ) );
    for ( unsigned int i = 0; i < points.size( ); i+=3 ){

        functionValues[ i + 0 ] = 1;
        functionValues[ i + 1 ] = 2;
        functionValues[ i + 2 ] = 3;

    }

    floatVector integratedVolumeResult;
    floatVector integratedVolumeAnswer = { 0.48573563, -0.32615985,  0.0715282 ,
                                           0.97147225, -0.6523191 ,  0.14305681,
                                           1.45720788, -0.97847915,  0.21458421 };

    floatVector origin = { 0., 0., 0. };

    error = dc.performRelativePositionVolumeIntegration( functionValues, 3, origin, integratedVolumeResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( integratedVolumeResult, integratedVolumeAnswer ) );

}

BOOST_AUTO_TEST_CASE( testDualContouring_performSurfaceIntegration ){
    /*!
     * Test surface integration over the reconstructed domain
     *
     */

    floatVector points =
        {
             0.88993453,  0.76248591, -0.93419017,  0.16825499, -0.64930721,
            -0.02736989,  0.3818072 , -0.56789234,  0.19636903, -0.01416032,
            -0.0726343 , -0.17899304,  0.58038292,  0.49538273, -0.16969211,
            -0.01199516,  0.83431932, -0.5944871 ,  0.68423515, -0.78823811,
            -0.30881808,  0.57273371, -0.3135023 , -0.0364649 ,  0.79618035,
            -0.36185261, -0.1701662 ,  0.4627748 , -0.22219916, -0.87057911,
            -0.34563548, -0.32495209,  0.90188535, -0.16083641, -0.67036656,
             0.04278752, -0.72313702, -0.62661373,  0.04437484, -0.62416924,
             0.3790917 , -0.99593042, -0.4127237 , -0.33142047,  0.91994519,
            -0.20164988,  0.48012869, -0.21532596,  0.92378245,  0.03305427,
             0.64008042,  0.36990853,  0.20280656,  0.82268547,  0.70533715,
             0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782,
            -0.75918987, -0.19216117,  0.82956824,  0.18528101, -0.52393903,
             0.70272336, -0.45594571,  0.18255914, -0.03938857, -0.78883385,
             0.85656866,  0.75744535, -0.77273419, -0.14567324,  0.1562869 ,
             0.67723809, -0.63066704, -0.10145066, -0.94238038,  0.58534946,
             0.29427379,  0.87440634,  0.81035162,  0.0477631 ,  0.82727563,
             0.74320746,  0.22778598,  0.58000421,  0.57222056,  0.09287177,
             0.47295155,  0.72360676, -0.41622328,  0.82540552,  0.32944547,
            -0.29634423,  0.19641306,  0.41306438,  0.36451089, -0.31113366,
            -0.93875619, -0.53893408, -0.36899228, -0.50751604,  0.42847307,
            -0.72351422, -0.1097771 , -0.51914444, -0.8399524 , -0.37058809,
            -0.53075988, -0.75968017, -0.13391175, -0.42047612,  0.93090777,
             0.87886273, -0.90105205,  0.37213616, -0.98616276,  0.92134844,
             0.98173517, -0.45791212, -0.41495029, -0.12509758, -0.25764766,
            -0.39942542, -0.07300332, -0.50951107,  0.47335924,  0.9226246 ,
            -0.8225396 ,  0.83189879,  0.97164758,  0.07520523,  0.96401949,
             0.09115673, -0.15473381, -0.12517667,  0.53923746,  0.17375421,
            -0.90864801,  0.13164519, -0.81746094, -0.04091139,  0.70761535,
             0.27949176, -0.72591527,  0.62191923,  0.58579313, -0.63652057,
            -0.82922138, -0.20341235, -0.36896014,  0.93666234,  0.99606496,
             0.02604581,  0.682005  , -0.83573812, -0.94639137, -0.93744807,
             0.41838408,  0.53512206,  0.73833179,  0.43511637,  0.07707317,
             0.94551536,  0.06136697, -0.88529597, -0.37223707, -0.13058219,
             0.19189113,  0.91949737,  0.41992734,  0.73392705, -0.81952681,
             0.97131812, -0.90971185, -0.83583963, -0.39223689, -0.2724691 ,
            -0.88024076, -0.2430722 ,  0.83700368, -0.35641362,  0.60074684,
             0.47984156,  0.7821192 ,  0.92905141, -0.66941659,  0.92052141,
             0.66587504,  0.64540513,  0.32427648,  0.20626864,  0.98728317,
            -0.46970216, -0.9715591 ,  0.52876951,  0.00966212,  0.4786293 ,
            -0.08921145,  0.93996844, -0.86556617, -0.43238781, -0.67206028,
            -0.05996167, -0.23409636, -0.79355771, -0.78246022, -0.52978658,
            -0.81887169,  0.32057883, -0.1984804 ,  0.20959359,  0.89333809,
            -0.26993128, -0.41085675, -0.2522169 ,  0.44736448,  0.62315547,
            -0.7794221 , -0.18013406,  0.26102072, -0.18430868,  0.99247162,
            -0.6383439 , -0.59488566,  0.36539787,  0.90600975, -0.81800276,
            -0.73782801,  0.79265261,  0.28872337, -0.16722855, -0.6239667 ,
            -0.45177053, -0.196289  ,  0.84749471, -0.41257737,  0.41538694,
            -0.57353971, -0.0891511 , -0.66645294,  0.6060595 ,  0.06013818,
             0.43417399,  0.85628154, -0.25946277, -0.17124759, -0.27269205,
            -0.65624534,  0.76907186, -0.15572072, -0.85328073,  0.64179789,
             0.7853888 , -0.73223644,  0.67797289, -0.66645228,  0.01685972,
            -0.27934479, -0.95963138, -0.78014927,  0.62397642, -0.41740632,
             0.82117605, -0.02459754,  0.15318692, -0.5935911 , -0.70041395,
            -0.79377562,  0.93858127, -0.02915152,  0.75699326,  0.39613993,
             0.82603044,  0.96104565,  0.2936643 ,  0.48070712,  0.83404655,
             0.14211884,  0.00691666,  0.80661045,  0.88463205, -0.88978482,
            -0.28513985, -0.30738596, -0.20413492,  0.21292329, -0.26989904,
            -0.85370154,  0.72208774, -0.55094191,  0.37445557,  0.54146363,
             0.57375785,  0.3679136 ,  0.65563202, -0.06865393, -0.59308587
        };

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    errorOut error = dc.loadPoints( &points );

    BOOST_CHECK( !error );

    floatVector functionValues( points.size( ) );
    for ( unsigned int i = 0; i < points.size( ); i+=3 ){

        functionValues[ i + 0 ] = 1;
        functionValues[ i + 1 ] = 2;
        functionValues[ i + 2 ] = 3;

    }

    floatVector integratedSurfaceResult;
    floatVector integratedSurfaceAnswer = { 22.3709189 , 44.74184179, 67.11275769 };

    error = dc.performSurfaceIntegration( functionValues, 3, integratedSurfaceResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( integratedSurfaceResult, integratedSurfaceAnswer ) );

}

BOOST_AUTO_TEST_CASE( testDualContouring_performPositionWeightedSurfaceIntegration ){
    /*!
     * Test position weighted surface integration over the reconstructed domain
     *
     */

    floatVector points =
        {
             0.88993453,  0.76248591, -0.93419017,  0.16825499, -0.64930721,
            -0.02736989,  0.3818072 , -0.56789234,  0.19636903, -0.01416032,
            -0.0726343 , -0.17899304,  0.58038292,  0.49538273, -0.16969211,
            -0.01199516,  0.83431932, -0.5944871 ,  0.68423515, -0.78823811,
            -0.30881808,  0.57273371, -0.3135023 , -0.0364649 ,  0.79618035,
            -0.36185261, -0.1701662 ,  0.4627748 , -0.22219916, -0.87057911,
            -0.34563548, -0.32495209,  0.90188535, -0.16083641, -0.67036656,
             0.04278752, -0.72313702, -0.62661373,  0.04437484, -0.62416924,
             0.3790917 , -0.99593042, -0.4127237 , -0.33142047,  0.91994519,
            -0.20164988,  0.48012869, -0.21532596,  0.92378245,  0.03305427,
             0.64008042,  0.36990853,  0.20280656,  0.82268547,  0.70533715,
             0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782,
            -0.75918987, -0.19216117,  0.82956824,  0.18528101, -0.52393903,
             0.70272336, -0.45594571,  0.18255914, -0.03938857, -0.78883385,
             0.85656866,  0.75744535, -0.77273419, -0.14567324,  0.1562869 ,
             0.67723809, -0.63066704, -0.10145066, -0.94238038,  0.58534946,
             0.29427379,  0.87440634,  0.81035162,  0.0477631 ,  0.82727563,
             0.74320746,  0.22778598,  0.58000421,  0.57222056,  0.09287177,
             0.47295155,  0.72360676, -0.41622328,  0.82540552,  0.32944547,
            -0.29634423,  0.19641306,  0.41306438,  0.36451089, -0.31113366,
            -0.93875619, -0.53893408, -0.36899228, -0.50751604,  0.42847307,
            -0.72351422, -0.1097771 , -0.51914444, -0.8399524 , -0.37058809,
            -0.53075988, -0.75968017, -0.13391175, -0.42047612,  0.93090777,
             0.87886273, -0.90105205,  0.37213616, -0.98616276,  0.92134844,
             0.98173517, -0.45791212, -0.41495029, -0.12509758, -0.25764766,
            -0.39942542, -0.07300332, -0.50951107,  0.47335924,  0.9226246 ,
            -0.8225396 ,  0.83189879,  0.97164758,  0.07520523,  0.96401949,
             0.09115673, -0.15473381, -0.12517667,  0.53923746,  0.17375421,
            -0.90864801,  0.13164519, -0.81746094, -0.04091139,  0.70761535,
             0.27949176, -0.72591527,  0.62191923,  0.58579313, -0.63652057,
            -0.82922138, -0.20341235, -0.36896014,  0.93666234,  0.99606496,
             0.02604581,  0.682005  , -0.83573812, -0.94639137, -0.93744807,
             0.41838408,  0.53512206,  0.73833179,  0.43511637,  0.07707317,
             0.94551536,  0.06136697, -0.88529597, -0.37223707, -0.13058219,
             0.19189113,  0.91949737,  0.41992734,  0.73392705, -0.81952681,
             0.97131812, -0.90971185, -0.83583963, -0.39223689, -0.2724691 ,
            -0.88024076, -0.2430722 ,  0.83700368, -0.35641362,  0.60074684,
             0.47984156,  0.7821192 ,  0.92905141, -0.66941659,  0.92052141,
             0.66587504,  0.64540513,  0.32427648,  0.20626864,  0.98728317,
            -0.46970216, -0.9715591 ,  0.52876951,  0.00966212,  0.4786293 ,
            -0.08921145,  0.93996844, -0.86556617, -0.43238781, -0.67206028,
            -0.05996167, -0.23409636, -0.79355771, -0.78246022, -0.52978658,
            -0.81887169,  0.32057883, -0.1984804 ,  0.20959359,  0.89333809,
            -0.26993128, -0.41085675, -0.2522169 ,  0.44736448,  0.62315547,
            -0.7794221 , -0.18013406,  0.26102072, -0.18430868,  0.99247162,
            -0.6383439 , -0.59488566,  0.36539787,  0.90600975, -0.81800276,
            -0.73782801,  0.79265261,  0.28872337, -0.16722855, -0.6239667 ,
            -0.45177053, -0.196289  ,  0.84749471, -0.41257737,  0.41538694,
            -0.57353971, -0.0891511 , -0.66645294,  0.6060595 ,  0.06013818,
             0.43417399,  0.85628154, -0.25946277, -0.17124759, -0.27269205,
            -0.65624534,  0.76907186, -0.15572072, -0.85328073,  0.64179789,
             0.7853888 , -0.73223644,  0.67797289, -0.66645228,  0.01685972,
            -0.27934479, -0.95963138, -0.78014927,  0.62397642, -0.41740632,
             0.82117605, -0.02459754,  0.15318692, -0.5935911 , -0.70041395,
            -0.79377562,  0.93858127, -0.02915152,  0.75699326,  0.39613993,
             0.82603044,  0.96104565,  0.2936643 ,  0.48070712,  0.83404655,
             0.14211884,  0.00691666,  0.80661045,  0.88463205, -0.88978482,
            -0.28513985, -0.30738596, -0.20413492,  0.21292329, -0.26989904,
            -0.85370154,  0.72208774, -0.55094191,  0.37445557,  0.54146363,
             0.57375785,  0.3679136 ,  0.65563202, -0.06865393, -0.59308587
        };

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    errorOut error = dc.loadPoints( &points );

    BOOST_CHECK( !error );

    floatVector functionValues( points.size( ) / 3 );
    for ( unsigned int i = 0; i < points.size( ) / 3; i++ ){

        functionValues[ i + 0 ] = 1;

    }

    floatVector integratedSurfaceResult;
    floatVector integratedSurfaceAnswer = { 0.903828 , -0.3204276, -0.3170567 };

    error = dc.performPositionWeightedSurfaceIntegration( functionValues, 1, integratedSurfaceResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( integratedSurfaceResult, integratedSurfaceAnswer ) );

}


BOOST_AUTO_TEST_CASE( testDualContouring_performSurfaceFluxIntegration ){
    /*!
     * Test surface flux integration over the reconstructed domain
     *
     */

    floatVector points =
        {
             0.88993453,  0.76248591, -0.93419017,  0.16825499, -0.64930721,
            -0.02736989,  0.3818072 , -0.56789234,  0.19636903, -0.01416032,
            -0.0726343 , -0.17899304,  0.58038292,  0.49538273, -0.16969211,
            -0.01199516,  0.83431932, -0.5944871 ,  0.68423515, -0.78823811,
            -0.30881808,  0.57273371, -0.3135023 , -0.0364649 ,  0.79618035,
            -0.36185261, -0.1701662 ,  0.4627748 , -0.22219916, -0.87057911,
            -0.34563548, -0.32495209,  0.90188535, -0.16083641, -0.67036656,
             0.04278752, -0.72313702, -0.62661373,  0.04437484, -0.62416924,
             0.3790917 , -0.99593042, -0.4127237 , -0.33142047,  0.91994519,
            -0.20164988,  0.48012869, -0.21532596,  0.92378245,  0.03305427,
             0.64008042,  0.36990853,  0.20280656,  0.82268547,  0.70533715,
             0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782,
            -0.75918987, -0.19216117,  0.82956824,  0.18528101, -0.52393903,
             0.70272336, -0.45594571,  0.18255914, -0.03938857, -0.78883385,
             0.85656866,  0.75744535, -0.77273419, -0.14567324,  0.1562869 ,
             0.67723809, -0.63066704, -0.10145066, -0.94238038,  0.58534946,
             0.29427379,  0.87440634,  0.81035162,  0.0477631 ,  0.82727563,
             0.74320746,  0.22778598,  0.58000421,  0.57222056,  0.09287177,
             0.47295155,  0.72360676, -0.41622328,  0.82540552,  0.32944547,
            -0.29634423,  0.19641306,  0.41306438,  0.36451089, -0.31113366,
            -0.93875619, -0.53893408, -0.36899228, -0.50751604,  0.42847307,
            -0.72351422, -0.1097771 , -0.51914444, -0.8399524 , -0.37058809,
            -0.53075988, -0.75968017, -0.13391175, -0.42047612,  0.93090777,
             0.87886273, -0.90105205,  0.37213616, -0.98616276,  0.92134844,
             0.98173517, -0.45791212, -0.41495029, -0.12509758, -0.25764766,
            -0.39942542, -0.07300332, -0.50951107,  0.47335924,  0.9226246 ,
            -0.8225396 ,  0.83189879,  0.97164758,  0.07520523,  0.96401949,
             0.09115673, -0.15473381, -0.12517667,  0.53923746,  0.17375421,
            -0.90864801,  0.13164519, -0.81746094, -0.04091139,  0.70761535,
             0.27949176, -0.72591527,  0.62191923,  0.58579313, -0.63652057,
            -0.82922138, -0.20341235, -0.36896014,  0.93666234,  0.99606496,
             0.02604581,  0.682005  , -0.83573812, -0.94639137, -0.93744807,
             0.41838408,  0.53512206,  0.73833179,  0.43511637,  0.07707317,
             0.94551536,  0.06136697, -0.88529597, -0.37223707, -0.13058219,
             0.19189113,  0.91949737,  0.41992734,  0.73392705, -0.81952681,
             0.97131812, -0.90971185, -0.83583963, -0.39223689, -0.2724691 ,
            -0.88024076, -0.2430722 ,  0.83700368, -0.35641362,  0.60074684,
             0.47984156,  0.7821192 ,  0.92905141, -0.66941659,  0.92052141,
             0.66587504,  0.64540513,  0.32427648,  0.20626864,  0.98728317,
            -0.46970216, -0.9715591 ,  0.52876951,  0.00966212,  0.4786293 ,
            -0.08921145,  0.93996844, -0.86556617, -0.43238781, -0.67206028,
            -0.05996167, -0.23409636, -0.79355771, -0.78246022, -0.52978658,
            -0.81887169,  0.32057883, -0.1984804 ,  0.20959359,  0.89333809,
            -0.26993128, -0.41085675, -0.2522169 ,  0.44736448,  0.62315547,
            -0.7794221 , -0.18013406,  0.26102072, -0.18430868,  0.99247162,
            -0.6383439 , -0.59488566,  0.36539787,  0.90600975, -0.81800276,
            -0.73782801,  0.79265261,  0.28872337, -0.16722855, -0.6239667 ,
            -0.45177053, -0.196289  ,  0.84749471, -0.41257737,  0.41538694,
            -0.57353971, -0.0891511 , -0.66645294,  0.6060595 ,  0.06013818,
             0.43417399,  0.85628154, -0.25946277, -0.17124759, -0.27269205,
            -0.65624534,  0.76907186, -0.15572072, -0.85328073,  0.64179789,
             0.7853888 , -0.73223644,  0.67797289, -0.66645228,  0.01685972,
            -0.27934479, -0.95963138, -0.78014927,  0.62397642, -0.41740632,
             0.82117605, -0.02459754,  0.15318692, -0.5935911 , -0.70041395,
            -0.79377562,  0.93858127, -0.02915152,  0.75699326,  0.39613993,
             0.82603044,  0.96104565,  0.2936643 ,  0.48070712,  0.83404655,
             0.14211884,  0.00691666,  0.80661045,  0.88463205, -0.88978482,
            -0.28513985, -0.30738596, -0.20413492,  0.21292329, -0.26989904,
            -0.85370154,  0.72208774, -0.55094191,  0.37445557,  0.54146363,
             0.57375785,  0.3679136 ,  0.65563202, -0.06865393, -0.59308587
        };

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    errorOut error = dc.loadPoints( &points );

    BOOST_CHECK( !error );

    floatVector functionValues( 2 * points.size( ) );
    for ( unsigned int i = 0; i < 2 * points.size( ); i+=6 ){

        functionValues[ i + 0 ] = 1;
        functionValues[ i + 1 ] = 2;
        functionValues[ i + 2 ] = 3;
        functionValues[ i + 3 ] = 4;
        functionValues[ i + 4 ] = 5;
        functionValues[ i + 5 ] = 6;

    }

    floatVector integratedSurfaceResult;
    floatVector integratedSurfaceAnswer = { 0., 0. };

    error = dc.performSurfaceFluxIntegration( functionValues, 6, integratedSurfaceResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( integratedSurfaceResult, integratedSurfaceAnswer ) );

}

BOOST_AUTO_TEST_CASE( testDualContouring_performRelativePositionSurfaceFluxIntegration ){
    /*!
     * Test surface flux integration over the reconstructed domain
     *
     */

    floatVector points =
        {
             0.88993453,  0.76248591, -0.93419017,  0.16825499, -0.64930721,
            -0.02736989,  0.3818072 , -0.56789234,  0.19636903, -0.01416032,
            -0.0726343 , -0.17899304,  0.58038292,  0.49538273, -0.16969211,
            -0.01199516,  0.83431932, -0.5944871 ,  0.68423515, -0.78823811,
            -0.30881808,  0.57273371, -0.3135023 , -0.0364649 ,  0.79618035,
            -0.36185261, -0.1701662 ,  0.4627748 , -0.22219916, -0.87057911,
            -0.34563548, -0.32495209,  0.90188535, -0.16083641, -0.67036656,
             0.04278752, -0.72313702, -0.62661373,  0.04437484, -0.62416924,
             0.3790917 , -0.99593042, -0.4127237 , -0.33142047,  0.91994519,
            -0.20164988,  0.48012869, -0.21532596,  0.92378245,  0.03305427,
             0.64008042,  0.36990853,  0.20280656,  0.82268547,  0.70533715,
             0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782,
            -0.75918987, -0.19216117,  0.82956824,  0.18528101, -0.52393903,
             0.70272336, -0.45594571,  0.18255914, -0.03938857, -0.78883385,
             0.85656866,  0.75744535, -0.77273419, -0.14567324,  0.1562869 ,
             0.67723809, -0.63066704, -0.10145066, -0.94238038,  0.58534946,
             0.29427379,  0.87440634,  0.81035162,  0.0477631 ,  0.82727563,
             0.74320746,  0.22778598,  0.58000421,  0.57222056,  0.09287177,
             0.47295155,  0.72360676, -0.41622328,  0.82540552,  0.32944547,
            -0.29634423,  0.19641306,  0.41306438,  0.36451089, -0.31113366,
            -0.93875619, -0.53893408, -0.36899228, -0.50751604,  0.42847307,
            -0.72351422, -0.1097771 , -0.51914444, -0.8399524 , -0.37058809,
            -0.53075988, -0.75968017, -0.13391175, -0.42047612,  0.93090777,
             0.87886273, -0.90105205,  0.37213616, -0.98616276,  0.92134844,
             0.98173517, -0.45791212, -0.41495029, -0.12509758, -0.25764766,
            -0.39942542, -0.07300332, -0.50951107,  0.47335924,  0.9226246 ,
            -0.8225396 ,  0.83189879,  0.97164758,  0.07520523,  0.96401949,
             0.09115673, -0.15473381, -0.12517667,  0.53923746,  0.17375421,
            -0.90864801,  0.13164519, -0.81746094, -0.04091139,  0.70761535,
             0.27949176, -0.72591527,  0.62191923,  0.58579313, -0.63652057,
            -0.82922138, -0.20341235, -0.36896014,  0.93666234,  0.99606496,
             0.02604581,  0.682005  , -0.83573812, -0.94639137, -0.93744807,
             0.41838408,  0.53512206,  0.73833179,  0.43511637,  0.07707317,
             0.94551536,  0.06136697, -0.88529597, -0.37223707, -0.13058219,
             0.19189113,  0.91949737,  0.41992734,  0.73392705, -0.81952681,
             0.97131812, -0.90971185, -0.83583963, -0.39223689, -0.2724691 ,
            -0.88024076, -0.2430722 ,  0.83700368, -0.35641362,  0.60074684,
             0.47984156,  0.7821192 ,  0.92905141, -0.66941659,  0.92052141,
             0.66587504,  0.64540513,  0.32427648,  0.20626864,  0.98728317,
            -0.46970216, -0.9715591 ,  0.52876951,  0.00966212,  0.4786293 ,
            -0.08921145,  0.93996844, -0.86556617, -0.43238781, -0.67206028,
            -0.05996167, -0.23409636, -0.79355771, -0.78246022, -0.52978658,
            -0.81887169,  0.32057883, -0.1984804 ,  0.20959359,  0.89333809,
            -0.26993128, -0.41085675, -0.2522169 ,  0.44736448,  0.62315547,
            -0.7794221 , -0.18013406,  0.26102072, -0.18430868,  0.99247162,
            -0.6383439 , -0.59488566,  0.36539787,  0.90600975, -0.81800276,
            -0.73782801,  0.79265261,  0.28872337, -0.16722855, -0.6239667 ,
            -0.45177053, -0.196289  ,  0.84749471, -0.41257737,  0.41538694,
            -0.57353971, -0.0891511 , -0.66645294,  0.6060595 ,  0.06013818,
             0.43417399,  0.85628154, -0.25946277, -0.17124759, -0.27269205,
            -0.65624534,  0.76907186, -0.15572072, -0.85328073,  0.64179789,
             0.7853888 , -0.73223644,  0.67797289, -0.66645228,  0.01685972,
            -0.27934479, -0.95963138, -0.78014927,  0.62397642, -0.41740632,
             0.82117605, -0.02459754,  0.15318692, -0.5935911 , -0.70041395,
            -0.79377562,  0.93858127, -0.02915152,  0.75699326,  0.39613993,
             0.82603044,  0.96104565,  0.2936643 ,  0.48070712,  0.83404655,
             0.14211884,  0.00691666,  0.80661045,  0.88463205, -0.88978482,
            -0.28513985, -0.30738596, -0.20413492,  0.21292329, -0.26989904,
            -0.85370154,  0.72208774, -0.55094191,  0.37445557,  0.54146363,
             0.57375785,  0.3679136 ,  0.65563202, -0.06865393, -0.59308587
        };

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    errorOut error = dc.loadPoints( &points );

    BOOST_CHECK( !error );

    floatVector functionValues( 2 * points.size( ) );
    for ( unsigned int i = 0; i < 2 * points.size( ); i+=6 ){

        functionValues[ i + 0 ] = 1;
        functionValues[ i + 1 ] = 2;
        functionValues[ i + 2 ] = 3;
        functionValues[ i + 3 ] = 4;
        functionValues[ i + 4 ] = 5;
        functionValues[ i + 5 ] = 6;

    }

    floatVector integratedSurfaceResult;
    floatVector integratedSurfaceAnswer = { 6.81813769, 20.45440861, 34.09068637, 13.63627627, 27.27255183, 40.90882794 };
    floatVector origin = { 0., 0., 0. };

    error = dc.performRelativePositionSurfaceFluxIntegration( functionValues, 6, origin, integratedSurfaceResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( integratedSurfaceResult, integratedSurfaceAnswer ) );

}

BOOST_AUTO_TEST_CASE( testDualContouring_exportConfiguration ){
    /*!
     * Test the export of the configuration file
     *
     */

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    YAML::Node config = dc.exportConfiguration( );

    BOOST_CHECK( config[ "type" ] );
    BOOST_CHECK( config[ "type" ].as< std::string >( ).compare( "dual_contouring" ) == 0 );

}

BOOST_AUTO_TEST_CASE( testDualContouring_getSurfaceSubdomains ){
    /*!
     * Test getting sub-domains on the surface of the body
     *
     */

    floatVector points =
        {
             0.88993453,  0.76248591, -0.93419017,  0.16825499, -0.64930721,
            -0.02736989,  0.3818072 , -0.56789234,  0.19636903, -0.01416032,
            -0.0726343 , -0.17899304,  0.58038292,  0.49538273, -0.16969211,
            -0.01199516,  0.83431932, -0.5944871 ,  0.68423515, -0.78823811,
            -0.30881808,  0.57273371, -0.3135023 , -0.0364649 ,  0.79618035,
            -0.36185261, -0.1701662 ,  0.4627748 , -0.22219916, -0.87057911,
            -0.34563548, -0.32495209,  0.90188535, -0.16083641, -0.67036656,
             0.04278752, -0.72313702, -0.62661373,  0.04437484, -0.62416924,
             0.3790917 , -0.99593042, -0.4127237 , -0.33142047,  0.91994519,
            -0.20164988,  0.48012869, -0.21532596,  0.92378245,  0.03305427,
             0.64008042,  0.36990853,  0.20280656,  0.82268547,  0.70533715,
             0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782,
            -0.75918987, -0.19216117,  0.82956824,  0.18528101, -0.52393903,
             0.70272336, -0.45594571,  0.18255914, -0.03938857, -0.78883385,
             0.85656866,  0.75744535, -0.77273419, -0.14567324,  0.1562869 ,
             0.67723809, -0.63066704, -0.10145066, -0.94238038,  0.58534946,
             0.29427379,  0.87440634,  0.81035162,  0.0477631 ,  0.82727563,
             0.74320746,  0.22778598,  0.58000421,  0.57222056,  0.09287177,
             0.47295155,  0.72360676, -0.41622328,  0.82540552,  0.32944547,
            -0.29634423,  0.19641306,  0.41306438,  0.36451089, -0.31113366,
            -0.93875619, -0.53893408, -0.36899228, -0.50751604,  0.42847307,
            -0.72351422, -0.1097771 , -0.51914444, -0.8399524 , -0.37058809,
            -0.53075988, -0.75968017, -0.13391175, -0.42047612,  0.93090777,
             0.87886273, -0.90105205,  0.37213616, -0.98616276,  0.92134844,
             0.98173517, -0.45791212, -0.41495029, -0.12509758, -0.25764766,
            -0.39942542, -0.07300332, -0.50951107,  0.47335924,  0.9226246 ,
            -0.8225396 ,  0.83189879,  0.97164758,  0.07520523,  0.96401949,
             0.09115673, -0.15473381, -0.12517667,  0.53923746,  0.17375421,
            -0.90864801,  0.13164519, -0.81746094, -0.04091139,  0.70761535,
             0.27949176, -0.72591527,  0.62191923,  0.58579313, -0.63652057,
            -0.82922138, -0.20341235, -0.36896014,  0.93666234,  0.99606496,
             0.02604581,  0.682005  , -0.83573812, -0.94639137, -0.93744807,
             0.41838408,  0.53512206,  0.73833179,  0.43511637,  0.07707317,
             0.94551536,  0.06136697, -0.88529597, -0.37223707, -0.13058219,
             0.19189113,  0.91949737,  0.41992734,  0.73392705, -0.81952681,
             0.97131812, -0.90971185, -0.83583963, -0.39223689, -0.2724691 ,
            -0.88024076, -0.2430722 ,  0.83700368, -0.35641362,  0.60074684,
             0.47984156,  0.7821192 ,  0.92905141, -0.66941659,  0.92052141,
             0.66587504,  0.64540513,  0.32427648,  0.20626864,  0.98728317,
            -0.46970216, -0.9715591 ,  0.52876951,  0.00966212,  0.4786293 ,
            -0.08921145,  0.93996844, -0.86556617, -0.43238781, -0.67206028,
            -0.05996167, -0.23409636, -0.79355771, -0.78246022, -0.52978658,
            -0.81887169,  0.32057883, -0.1984804 ,  0.20959359,  0.89333809,
            -0.26993128, -0.41085675, -0.2522169 ,  0.44736448,  0.62315547,
            -0.7794221 , -0.18013406,  0.26102072, -0.18430868,  0.99247162,
            -0.6383439 , -0.59488566,  0.36539787,  0.90600975, -0.81800276,
            -0.73782801,  0.79265261,  0.28872337, -0.16722855, -0.6239667 ,
            -0.45177053, -0.196289  ,  0.84749471, -0.41257737,  0.41538694,
            -0.57353971, -0.0891511 , -0.66645294,  0.6060595 ,  0.06013818,
             0.43417399,  0.85628154, -0.25946277, -0.17124759, -0.27269205,
            -0.65624534,  0.76907186, -0.15572072, -0.85328073,  0.64179789,
             0.7853888 , -0.73223644,  0.67797289, -0.66645228,  0.01685972,
            -0.27934479, -0.95963138, -0.78014927,  0.62397642, -0.41740632,
             0.82117605, -0.02459754,  0.15318692, -0.5935911 , -0.70041395,
            -0.79377562,  0.93858127, -0.02915152,  0.75699326,  0.39613993,
             0.82603044,  0.96104565,  0.2936643 ,  0.48070712,  0.83404655,
             0.14211884,  0.00691666,  0.80661045,  0.88463205, -0.88978482,
            -0.28513985, -0.30738596, -0.20413492,  0.21292329, -0.26989904,
            -0.85370154,  0.72208774, -0.55094191,  0.37445557,  0.54146363,
             0.57375785,  0.3679136 ,  0.65563202, -0.06865393, -0.59308587
        };

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    errorOut error = dc.loadPoints( &points );

    BOOST_CHECK( !error );

    uIntVector subdomainNodeCountAnswer = { 31, 50, 26, 31, 33, 37, 22, 42, 24, 32 };

    uIntVector subdomainNodesAnswer = { 
        244, 241, 240, 225, 224, 221, 220, 149, 148, 147, 146, 141, 140,
        121, 120, 115, 114, 113, 112,  45,  44,  13,  12,  11,  10,   9,
          8,   3,   2,   1,   0, 315, 314, 307, 306, 305, 304, 291, 290,
        287, 286, 275, 271, 270, 267, 266, 211, 210, 207, 206, 205, 204,
        191, 190, 187, 186, 185, 184, 169, 168,  99,  98,  87,  86,  79,
         78,  77,  76,  75,  74,  73,  72,  66,  65,  64,  63,  62,  61,
         60,  57,  56, 327, 326, 325, 324, 319, 318, 299, 298, 295, 294,
        279, 219, 218, 215, 214, 199, 198, 197, 171, 111, 110, 107, 106,
         81,  80,  67, 311, 310, 303, 302, 301, 300, 283, 282, 281, 280,
        261, 260, 203, 202, 201, 200, 183, 182, 181, 180, 161,  91,  90,
         89,  88,  85,  84,  83,  82,  71,  70, 263, 262, 247, 246, 243,
        242, 227, 226, 223, 222, 163, 162, 160, 145, 144, 143, 142, 119,
        118, 117, 116,  55,  54,  51,  50,  17,  16,  15,  14,   7,   6,
          5,   4, 272, 257, 256, 253, 252, 249, 248, 245, 237, 236, 233,
        232, 229, 228, 151, 150, 132, 131, 130, 127, 126, 123, 122,  47,
         46,  37,  36,  35,  34,  30,  29,  28,  22,  21,  20,  19,  18,
        312, 309, 308, 289, 288, 285, 284, 269, 268, 265, 264, 209, 208,
        189, 188, 167, 166, 165, 164,  92,  53,  52, 274, 259, 258, 255,
        254, 251, 250, 239, 238, 235, 234, 231, 230, 170, 159, 158, 157,
        139, 138, 137, 136, 129, 128, 125, 124,  59,  58,  49,  48,  43,
         42,  41,  40,  39,  33,  32,  31,  27,  26,  25,  24,  23, 278,
        277, 276, 273, 196, 194, 179, 178, 177, 176, 175, 174, 173, 172,
        156, 155, 154, 153, 152, 135, 134, 133,  69,  38, 323, 322, 321,
        320, 317, 316, 313, 297, 296, 293, 292, 217, 216, 213, 212, 195,
        193, 192, 109, 108, 105, 104, 103, 102, 101, 100,  97,  96,  95,
         94,  93,  68
    };

    uIntVector subdomainNodeCountResult;
    uIntVector subdomainNodesResult;

    floatType minDistance = 1.0;

    error = dc.getSurfaceSubdomains( minDistance, subdomainNodeCountResult, subdomainNodesResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( subdomainNodeCountResult, subdomainNodeCountAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( subdomainNodesResult, subdomainNodesAnswer ) );

}

BOOST_AUTO_TEST_CASE( testDualContouring_getBoundaryInformation ){
    /*!
     * Check that the boundary points and cells can be extracted
     *
     */

    floatVector points =
        {
             0.88993453,  0.76248591, -0.93419017,  0.16825499, -0.64930721,
            -0.02736989,  0.3818072 , -0.56789234,  0.19636903, -0.01416032,
            -0.0726343 , -0.17899304,  0.58038292,  0.49538273, -0.16969211,
            -0.01199516,  0.83431932, -0.5944871 ,  0.68423515, -0.78823811,
            -0.30881808,  0.57273371, -0.3135023 , -0.0364649 ,  0.79618035,
            -0.36185261, -0.1701662 ,  0.4627748 , -0.22219916, -0.87057911,
            -0.34563548, -0.32495209,  0.90188535, -0.16083641, -0.67036656,
             0.04278752, -0.72313702, -0.62661373,  0.04437484, -0.62416924,
             0.3790917 , -0.99593042, -0.4127237 , -0.33142047,  0.91994519,
            -0.20164988,  0.48012869, -0.21532596,  0.92378245,  0.03305427,
             0.64008042,  0.36990853,  0.20280656,  0.82268547,  0.70533715,
             0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782,
            -0.75918987, -0.19216117,  0.82956824,  0.18528101, -0.52393903,
             0.70272336, -0.45594571,  0.18255914, -0.03938857, -0.78883385,
             0.85656866,  0.75744535, -0.77273419, -0.14567324,  0.1562869 ,
             0.67723809, -0.63066704, -0.10145066, -0.94238038,  0.58534946,
             0.29427379,  0.87440634,  0.81035162,  0.0477631 ,  0.82727563,
             0.74320746,  0.22778598,  0.58000421,  0.57222056,  0.09287177,
             0.47295155,  0.72360676, -0.41622328,  0.82540552,  0.32944547,
            -0.29634423,  0.19641306,  0.41306438,  0.36451089, -0.31113366,
            -0.93875619, -0.53893408, -0.36899228, -0.50751604,  0.42847307,
            -0.72351422, -0.1097771 , -0.51914444, -0.8399524 , -0.37058809,
            -0.53075988, -0.75968017, -0.13391175, -0.42047612,  0.93090777,
             0.87886273, -0.90105205,  0.37213616, -0.98616276,  0.92134844,
             0.98173517, -0.45791212, -0.41495029, -0.12509758, -0.25764766,
            -0.39942542, -0.07300332, -0.50951107,  0.47335924,  0.9226246 ,
            -0.8225396 ,  0.83189879,  0.97164758,  0.07520523,  0.96401949,
             0.09115673, -0.15473381, -0.12517667,  0.53923746,  0.17375421,
            -0.90864801,  0.13164519, -0.81746094, -0.04091139,  0.70761535,
             0.27949176, -0.72591527,  0.62191923,  0.58579313, -0.63652057,
            -0.82922138, -0.20341235, -0.36896014,  0.93666234,  0.99606496,
             0.02604581,  0.682005  , -0.83573812, -0.94639137, -0.93744807,
             0.41838408,  0.53512206,  0.73833179,  0.43511637,  0.07707317,
             0.94551536,  0.06136697, -0.88529597, -0.37223707, -0.13058219,
             0.19189113,  0.91949737,  0.41992734,  0.73392705, -0.81952681,
             0.97131812, -0.90971185, -0.83583963, -0.39223689, -0.2724691 ,
            -0.88024076, -0.2430722 ,  0.83700368, -0.35641362,  0.60074684,
             0.47984156,  0.7821192 ,  0.92905141, -0.66941659,  0.92052141,
             0.66587504,  0.64540513,  0.32427648,  0.20626864,  0.98728317,
            -0.46970216, -0.9715591 ,  0.52876951,  0.00966212,  0.4786293 ,
            -0.08921145,  0.93996844, -0.86556617, -0.43238781, -0.67206028,
            -0.05996167, -0.23409636, -0.79355771, -0.78246022, -0.52978658,
            -0.81887169,  0.32057883, -0.1984804 ,  0.20959359,  0.89333809,
            -0.26993128, -0.41085675, -0.2522169 ,  0.44736448,  0.62315547,
            -0.7794221 , -0.18013406,  0.26102072, -0.18430868,  0.99247162,
            -0.6383439 , -0.59488566,  0.36539787,  0.90600975, -0.81800276,
            -0.73782801,  0.79265261,  0.28872337, -0.16722855, -0.6239667 ,
            -0.45177053, -0.196289  ,  0.84749471, -0.41257737,  0.41538694,
            -0.57353971, -0.0891511 , -0.66645294,  0.6060595 ,  0.06013818,
             0.43417399,  0.85628154, -0.25946277, -0.17124759, -0.27269205,
            -0.65624534,  0.76907186, -0.15572072, -0.85328073,  0.64179789,
             0.7853888 , -0.73223644,  0.67797289, -0.66645228,  0.01685972,
            -0.27934479, -0.95963138, -0.78014927,  0.62397642, -0.41740632,
             0.82117605, -0.02459754,  0.15318692, -0.5935911 , -0.70041395,
            -0.79377562,  0.93858127, -0.02915152,  0.75699326,  0.39613993,
             0.82603044,  0.96104565,  0.2936643 ,  0.48070712,  0.83404655,
             0.14211884,  0.00691666,  0.80661045,  0.88463205, -0.88978482,
            -0.28513985, -0.30738596, -0.20413492,  0.21292329, -0.26989904,
            -0.85370154,  0.72208774, -0.55094191,  0.37445557,  0.54146363,
             0.57375785,  0.3679136 ,  0.65563202, -0.06865393, -0.59308587
        };

    YAML::Node yf = YAML::LoadFile( "volumeReconstruction_dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    errorOut error = dc.loadPoints( &points );

    BOOST_CHECK( !error );

    error = dc.evaluate( );

    BOOST_CHECK( !error );

    const uIntVector boundaryCellsAnswer
        =
        {
             1,   2,   3,   4,   8,   9,  10,  11,  12,  15,  16,  17,  18,
            19,  22,  23,  24,  25,  26,  28,  29,  30,  31,  32,  35,  36,
            37,  38,  39,  50,  51,  52,  53,  57,  58,  59,  60,  61,  63,
            64,  68,  70,  71,  72,  74,  75,  77,  78,  79,  80,  81,  84,
            85,  86,  87,  88,  99, 100, 101, 102, 103, 105, 106, 109, 110,
           112, 113, 117, 119, 120, 123, 124, 126, 127, 128, 129, 130, 134,
           135, 136, 147, 148, 149, 150, 151, 152, 154, 155, 159, 161, 162,
           166, 168, 169, 172, 173, 175, 176, 178, 179, 180, 182, 183, 184,
           185, 196, 197, 198, 199, 200, 201, 203, 204, 207, 208, 210, 211,
           214, 215, 217, 222, 224, 227, 228, 229, 231, 232, 233, 234, 235,
           245, 246, 247, 248, 249, 250, 252, 253, 254, 255, 256, 257, 259,
           260, 261, 262, 263, 264, 266, 267, 268, 269, 270, 271, 273, 274,
           275, 276, 277, 278, 280, 281, 282, 283, 284
        };

    const floatVector boundaryPointsAnswer
        =
        {
            9.81978157e-01,  9.40322163e-01,  5.68747397e-01,  9.81978146e-01,
            8.99162053e-01,  4.90738908e-01,  9.81977762e-01,  9.55152728e-01,
            1.20890063e-01,  9.81977519e-01,  8.41360993e-01, -2.10611106e-02,
            9.81977847e-01,  8.82521472e-01, -3.92249300e-01,  9.81979448e-01,
            8.00508426e-01, -5.85222104e-01,  9.81954444e-01,  9.14267059e-01,
           -8.57547279e-01,  9.81954430e-01,  8.31602564e-01, -9.40475616e-01,
            9.81977963e-01,  5.76093790e-01,  9.48248922e-01,  9.81978513e-01,
            4.34743612e-01,  8.87174313e-01,  9.81978230e-01,  6.41467059e-01,
            6.78323586e-01,  9.81978838e-01,  4.51757589e-01,  4.31900227e-01,
            9.81978511e-01,  5.95660849e-01,  1.02364376e-01,  9.81978185e-01,
            3.63822353e-01, -1.12890329e-01,  9.81978794e-01,  5.27509451e-01,
           -4.29288619e-01,  9.81977865e-01,  4.26535716e-01, -6.58068398e-01,
            9.81978506e-01,  6.16870196e-01, -8.96930221e-01,  9.81977450e-01,
            4.51099458e-01, -9.52603465e-01,  9.81978436e-01,  9.93479725e-02,
            9.31645684e-01,  9.81977599e-01, -6.99940179e-02,  8.39115529e-01,
            9.81978285e-01,  7.56921516e-02,  5.48220799e-01,  9.81977576e-01,
           -4.02119510e-03,  3.48060460e-01,  9.81977805e-01,  1.40205746e-01,
            1.64851845e-02,  9.81977566e-01, -4.14774828e-02, -1.04864662e-01,
            9.81977017e-01,  8.11653891e-02, -4.63019954e-01,  9.81976773e-01,
           -4.42845457e-02, -6.00061335e-01,  9.81976565e-01,  1.35315043e-01,
           -8.50783023e-01,  9.81976857e-01, -3.69296023e-02, -8.94358393e-01,
            9.81976894e-01, -3.79089937e-01,  6.17207056e-01,  9.81977922e-01,
           -5.03462032e-01,  4.56499814e-01,  9.81978365e-01, -3.10076519e-01,
            6.60033724e-02,  9.81978811e-01, -6.12771414e-01, -1.04998181e-01,
            9.81977767e-01, -4.69622312e-01, -4.12907596e-01,  9.81977514e-01,
           -6.50755774e-01, -6.39308635e-01,  9.81948735e-01, -8.88373828e-01,
            9.70915427e-01,  9.81948829e-01, -9.46219614e-01,  8.88007317e-01,
            9.81978368e-01, -8.26962459e-01,  6.34827038e-01,  9.81978192e-01,
           -8.67100455e-01,  4.07933986e-01,  9.81978191e-01, -8.46294923e-01,
            5.86052083e-02,  9.81980143e-01, -9.65608873e-01, -1.02107089e-01,
            9.81980509e-01, -9.13218826e-01, -4.92111518e-01,  9.81978200e-01,
           -9.34020364e-01, -6.39555168e-01,  9.81975986e-01, -8.85265259e-01,
           -8.88572611e-01,  9.81927747e-01, -9.37609487e-01, -9.12839547e-01,
            7.43467743e-01,  9.86448965e-01,  4.92379525e-01,  7.34764365e-01,
            9.55633447e-01,  6.31473062e-01,  6.99658649e-01, -3.62403730e-01,
            9.32107705e-01,  7.91647578e-01, -4.78945027e-01,  9.22671789e-01,
            8.16342107e-01, -3.61353191e-01, -9.34044989e-01,  8.01197717e-01,
           -5.18683124e-01, -9.77796058e-01,  3.35080463e-01,  9.27501342e-01,
           -9.47255166e-01,  2.75540839e-01,  8.44068550e-01, -9.24366550e-01,
            1.96993146e-01,  5.74229003e-01,  8.51650515e-01,  1.83334581e-01,
            4.39393890e-01,  9.36439759e-01,  2.94249855e-01,  6.06159994e-01,
           -9.73290103e-01,  2.41016567e-01,  4.52371754e-01, -8.95788815e-01,
            2.24135853e-01,  9.17470369e-02, -9.18677347e-01,  1.30303882e-01,
           -2.91786392e-02, -8.77712569e-01,  2.49054013e-01, -8.40896969e-01,
           -9.48350247e-01,  2.81324864e-01, -9.22984290e-01, -8.54696968e-01,
           -2.51999445e-01,  9.12454780e-01, -5.87988525e-01, -3.23354438e-01,
            8.71741515e-01, -4.67117949e-01, -2.65092507e-01,  5.42648357e-01,
           -8.21963203e-01, -2.52692860e-01,  3.85507888e-01, -8.89057942e-01,
           -2.06149783e-01,  7.81553554e-02, -8.38725491e-01, -2.65505328e-01,
           -3.82197159e-02, -9.45846059e-01, -2.02661726e-01, -3.43521649e-01,
           -9.45846025e-01, -1.78107592e-01, -5.92286614e-01, -8.75497653e-01,
           -2.60653019e-01, -9.14693913e-01,  9.20580476e-01, -1.93284440e-01,
           -9.48936131e-01,  9.28315595e-01, -8.39943166e-01,  9.89927313e-01,
            6.32728986e-01, -7.64310550e-01,  9.76342175e-01,  5.09521193e-01,
           -7.64020000e-01,  8.68296570e-01, -8.29366476e-01, -7.34927372e-01,
            7.98786976e-01, -8.80410711e-01, -7.36732307e-01,  5.26067315e-01,
           -6.44522633e-01, -7.99101891e-01,  4.71763278e-01, -5.70369838e-01,
           -6.41053717e-01,  5.57636084e-01, -8.80410906e-01, -6.63993184e-01,
            3.66609701e-01, -8.62317987e-01, -6.68539695e-01,  8.31128371e-02,
           -9.78086465e-01, -7.14124855e-01, -7.62410326e-02, -9.62564719e-01,
           -7.67000763e-01, -9.36226947e-01, -3.73604833e-02, -7.14464960e-01,
           -9.66645167e-01,  9.92570303e-02, -9.59883309e-01,  9.65385049e-01,
            3.82743210e-01, -9.59874424e-01,  9.47734126e-01,  5.50321336e-01,
           -9.59877607e-01,  9.39511627e-01, -6.05808957e-03, -9.59883138e-01,
            9.08581920e-01,  1.46094603e-01, -9.59814712e-01,  9.48634068e-01,
           -9.54870370e-01, -9.59874946e-01,  8.63010407e-01, -9.08509171e-01,
           -9.59874278e-01,  6.44503709e-01,  4.27118873e-01, -9.59873990e-01,
            4.48783581e-01,  5.59314888e-01, -9.59874162e-01,  7.04578407e-01,
           -1.92869675e-02, -9.59874186e-01,  4.58151118e-01,  9.41591769e-02,
           -9.59874323e-01,  1.01733457e-01,  8.53733657e-01, -9.59874604e-01,
           -8.46515354e-02,  9.42866685e-01, -9.59873633e-01,  7.85480639e-03,
            3.86742537e-01, -9.59874012e-01, -6.84158582e-02,  6.02602180e-01,
           -9.59874545e-01,  7.15656902e-02, -7.77390989e-02, -9.59874279e-01,
           -1.83010494e-01,  8.26482464e-02, -9.59873946e-01,  7.98602789e-02,
           -6.38943574e-01, -9.59874598e-01, -4.13706110e-02, -4.39843291e-01,
           -9.59875043e-01, -4.07148290e-01,  8.42011487e-01, -9.59874659e-01,
           -5.62663049e-01,  8.91865312e-01, -9.59873924e-01, -4.41357985e-01,
            3.80707745e-01, -9.59874673e-01, -5.93705425e-01,  5.30984025e-01,
           -9.59874303e-01, -4.93368089e-01, -4.36234271e-02, -9.59873568e-01,
           -6.75444915e-01,  3.85764289e-02, -9.59874640e-01, -4.37739168e-01,
           -5.92313158e-01, -9.59874610e-01, -6.35160461e-01, -3.56480639e-01,
           -9.59874272e-01, -8.75112162e-01,  3.75455872e-01, -9.59874513e-01,
           -9.36424985e-01,  5.58423915e-01, -9.59874627e-01, -8.78550689e-01,
           -4.86605342e-01, -9.59873855e-01, -9.39090334e-01, -3.93448881e-01,
            9.17671617e-01,  9.96313781e-01,  4.26694890e-01,  8.61769548e-01,
            9.96313140e-01,  5.72468127e-01,  8.70757811e-01,  9.96313146e-01,
           -7.84698033e-02,  8.06451656e-01,  9.96314089e-01,  1.49326536e-01,
            8.89841551e-01,  9.96313021e-01, -5.93496729e-01,  7.78621271e-01,
            9.96312961e-01, -4.21668405e-01,  9.00986212e-01,  9.96280881e-01,
           -9.13220649e-01,  8.08874920e-01,  9.96313796e-01, -8.42142903e-01,
            9.08967768e-01,  7.93162344e-01,  9.00817732e-01,  8.33054687e-01,
            7.75702609e-01,  9.45289401e-01,  8.36357120e-01, -2.18248912e-01,
            9.00190147e-01,  7.66010128e-01, -2.23924418e-01,  8.39577987e-01,
            9.32180297e-01, -1.94448839e-01, -8.50607254e-01,  8.48902058e-01,
           -2.19071222e-01, -9.52427862e-01,  9.28346764e-01, -7.50009521e-01,
            9.61535926e-01,  8.61994801e-01, -6.84945067e-01,  9.83284474e-01,
            9.17037261e-01, -7.86994563e-01, -9.71858794e-01,  8.84476646e-01,
           -7.04129464e-01, -9.53475104e-01,  9.00840207e-01, -9.86358128e-01,
            9.09755659e-01,  8.67870479e-01, -9.86411429e-01,  8.05589989e-01,
            9.48978924e-01, -9.86410912e-01,  5.56948596e-01,  8.75973889e-01,
           -9.86410305e-01,  3.81555111e-01,  9.08975109e-01, -9.86413800e-01,
            7.18282033e-02,  8.40445181e-01, -9.86414661e-01, -6.07326871e-02,
            9.13447906e-01, -9.86411273e-01, -3.82111630e-01,  8.47047251e-01,
           -9.86407573e-01, -6.24998881e-01,  9.15529306e-01, -9.86360892e-01,
           -9.06709354e-01,  8.34548245e-01, -9.86364821e-01, -9.89995695e-01,
            5.75404080e-01,  9.96313525e-01, -1.69815808e-02,  4.81717939e-01,
            9.96313156e-01,  9.47614151e-02,  5.45601489e-01,  9.96313439e-01,
           -6.21390491e-01,  4.06690453e-01,  9.96313060e-01, -4.44500005e-01,
            6.08200126e-01,  9.96313473e-01, -9.25070933e-01,  4.56146974e-01,
            9.96313176e-01, -8.76146844e-01,  5.79184409e-01,  7.92520658e-01,
            8.60500468e-01,  4.20043637e-01,  6.95516684e-01,  9.07931066e-01,
            5.31395986e-01,  8.95473569e-01,  4.94084796e-01,  3.95827465e-01,
            8.05811670e-01,  6.20225995e-01,  6.25225572e-01, -9.86411594e-01,
            8.88498841e-01,  4.49008281e-01, -9.86411653e-01,  8.56062181e-01,
            6.18648154e-01, -9.86410881e-01,  5.97759845e-01,  3.78977913e-01,
           -9.86410049e-01,  4.75409191e-01,  5.67773721e-01, -9.86410335e-01,
            1.26643566e-01,  3.61704139e-01, -9.86410786e-01,  1.85590234e-02,
            5.67974943e-01, -9.86411413e-01, -4.11704481e-01,  4.00540799e-01,
           -9.86355650e-01, -5.99815240e-01,  5.98632046e-01, -9.86355738e-01,
           -8.96395981e-01,  4.71277104e-01, -9.86355502e-01, -9.02525952e-01,
            1.17177091e-01,  9.96312628e-01, -1.22604988e-01, -8.75728430e-02,
            9.96312416e-01, -1.29311176e-02,  1.45817198e-01,  9.96312671e-01,
           -6.17236495e-01, -5.00100235e-02,  9.96312607e-01, -4.07084474e-01,
            7.43770712e-02,  7.88769624e-01,  4.21020976e-01, -2.43275126e-02,
            7.48954329e-01,  5.74191987e-01,  2.82021640e-02,  3.36455667e-01,
            8.40059227e-01, -3.47579577e-02,  2.31206847e-01,  8.96339742e-01,
           -2.37025498e-02, -1.77047605e-01, -9.55213458e-01, -1.53011522e-01,
           -1.56055731e-01, -9.04881080e-01,  4.69189754e-02, -7.99222808e-01,
           -8.78001483e-01, -3.99389850e-02, -7.47986112e-01, -9.25830136e-01,
            1.22977383e-01, -9.86377867e-01,  9.60236848e-01, -5.01728581e-02,
           -9.86377876e-01,  8.95878749e-01,  3.61877999e-02, -9.86410831e-01,
            6.69353428e-01, -8.57901456e-02, -9.86410654e-01,  3.97707717e-01,
            5.09303854e-02, -9.86410324e-01,  9.56789360e-02, -4.12038296e-02,
           -9.86410818e-01, -7.56669651e-02,  1.02135780e-01, -9.86355517e-01,
           -3.40897637e-01, -2.29257486e-02, -9.86355691e-01, -5.60428863e-01,
           -4.51334141e-01,  9.96312914e-01, -1.44142213e-01, -5.98376038e-01,
            9.96313551e-01,  4.52385969e-02, -3.29465968e-01,  8.34895452e-01,
            4.12061156e-01, -5.46395875e-01,  8.76366431e-01,  6.36927500e-01,
           -4.35685582e-01,  7.43725342e-01, -6.10425288e-01, -5.85140475e-01,
            7.45889853e-01, -4.59724994e-01, -4.76106274e-01,  7.27232316e-01,
           -9.29084311e-01, -5.64337947e-01,  7.14929188e-01, -8.13315846e-01,
           -3.47824716e-01,  2.55229834e-01,  8.54474046e-01, -5.36419012e-01,
            2.32062323e-01,  9.50854666e-01, -4.00205904e-01, -1.97726219e-01,
           -9.30324633e-01, -5.39527166e-01, -1.86173679e-01, -9.80657097e-01,
           -4.06440490e-01, -8.16289920e-01,  9.84937645e-01, -5.89756512e-01,
           -8.11326852e-01,  9.50210849e-01, -4.22550369e-01, -9.86411167e-01,
            5.38824000e-01, -5.56349367e-01, -9.86410719e-01,  3.32837863e-01,
           -3.11723359e-01, -9.86410700e-01,  2.03756782e-02, -4.82916724e-01,
           -9.86410179e-01, -3.17985433e-02, -3.37458723e-01, -9.86410471e-01,
           -3.47372413e-01, -5.70261328e-01, -9.86410392e-01, -5.39025359e-01,
           -8.75845777e-01,  9.96318482e-01,  4.08361057e-01, -9.47240533e-01,
            9.96318118e-01,  4.80359825e-01, -8.09857588e-01,  9.96310459e-01,
           -2.29491473e-02, -8.88483022e-01,  9.96315120e-01,  1.58591749e-01,
           -8.70065367e-01,  9.96299558e-01, -9.45096978e-01, -9.13723461e-01,
            9.96299559e-01, -9.03825957e-01, -7.78585507e-01,  7.13150026e-01,
           -8.34011455e-01, -8.84691424e-01,  7.93500814e-01, -9.49780228e-01,
           -8.21779197e-01,  3.05396651e-01,  8.61721974e-01, -8.88467451e-01,
            2.27309567e-01,  9.03587525e-01, -8.63894298e-01,  2.57169044e-01,
           -7.17525822e-01, -9.16138114e-01,  3.13706514e-01, -4.84171042e-01,
           -8.10642126e-01, -7.19710338e-01,  9.57993541e-01, -8.96344089e-01,
           -7.68128121e-01,  8.57138311e-01, -8.63787256e-01, -8.74914894e-01,
            7.86903801e-03, -9.14719685e-01, -8.47359897e-01, -1.17252882e-01,
           -8.24329467e-01, -9.86410969e-01,  5.08861206e-01, -8.62184078e-01,
           -9.86410748e-01,  3.88935839e-01, -8.13408551e-01, -9.86408915e-01,
           -4.23940444e-01, -9.10195796e-01, -9.86410263e-01, -5.91035481e-01,
            9.26076087e-01,  9.40321523e-01,  7.14520634e-01,  8.53066164e-01,
            9.09506520e-01,  7.97931682e-01,  9.04380225e-01,  8.59787091e-01,
           -9.96148968e-01,  8.23412899e-01,  9.42451564e-01, -9.96148956e-01,
            9.06064881e-01,  5.58634055e-01,  9.92720591e-01,  8.80312714e-01,
            4.55871996e-01,  9.92721004e-01,  8.73458777e-01,  3.88303015e-01,
           -9.96179860e-01,  7.95885627e-01,  5.82258280e-01, -9.96179968e-01,
            9.56225717e-01,  1.37936091e-01,  9.92720706e-01,  8.10605469e-01,
           -2.01295471e-02,  9.92720005e-01,  8.98698618e-01, -6.15519855e-02,
           -9.96179001e-01,  7.90179653e-01,  4.78962172e-02, -9.96180027e-01,
            9.11630024e-01, -4.38744537e-01,  8.17829703e-01,  8.57999541e-01,
           -5.44009482e-01,  9.00923240e-01,  9.49416859e-01, -5.39944622e-01,
           -8.51654871e-01,  8.99620345e-01, -3.36730808e-01, -8.32224381e-01,
            9.15596772e-01, -8.23309373e-01,  9.92663976e-01,  8.34487923e-01,
           -8.81154914e-01,  9.92663656e-01,  9.00946687e-01, -9.37613415e-01,
           -9.96125887e-01,  8.36054022e-01, -8.39385306e-01, -9.96179530e-01,
            5.94463293e-01,  9.86448850e-01,  2.98798688e-01,  3.99499806e-01,
            9.26288511e-01,  2.45048391e-01,  6.10235096e-01,  8.73672646e-01,
           -9.96179533e-01,  4.87133615e-01,  9.27501639e-01, -9.96179255e-01,
            6.48161747e-01,  5.75936576e-01,  9.92720411e-01,  4.56545943e-01,
            4.15026225e-01,  9.92720565e-01,  5.81574167e-01,  4.27013900e-01,
           -9.96179419e-01,  4.76890961e-01,  6.35763793e-01, -9.96178996e-01,
            5.43076662e-01,  1.36896008e-01,  9.92720456e-01,  3.21398275e-01,
           -3.65483398e-02,  9.92720332e-01,  5.13084683e-01, -1.48377034e-01,
           -9.96179774e-01,  3.82052347e-01,  3.67853841e-02, -9.96179058e-01,
            5.22575122e-01, -3.57892012e-01,  9.92720294e-01,  4.18983370e-01,
           -5.86305114e-01,  9.92720455e-01,  5.65522952e-01, -5.68156823e-01,
           -9.96179577e-01,  3.53368931e-01, -4.49580904e-01, -9.96179395e-01,
            5.85305271e-01, -8.54662412e-01,  9.92720233e-01,  4.09012648e-01,
           -9.59865792e-01,  9.92721086e-01,  6.42300639e-01, -8.88136476e-01,
           -9.96179309e-01,  4.55048589e-01, -8.05993557e-01, -9.96179074e-01,
            8.20045190e-02,  9.36152604e-01,  1.50684988e-01, -1.07549696e-01,
            8.88295616e-01,  2.00516375e-01,  5.57068924e-02,  8.29022178e-01,
           -8.57271793e-01, -5.61722234e-02,  9.12454843e-01, -7.98140547e-01,
            7.11692259e-02,  6.02072860e-01,  7.64547094e-01, -6.22435607e-02,
            4.70841734e-01,  7.52955612e-01,  9.42409143e-02,  4.36186801e-01,
           -7.88668068e-01, -7.23597438e-02,  5.74928669e-01, -7.99074598e-01,
            1.05094551e-01,  6.94135178e-02,  9.92720446e-01,  7.68307825e-04,
           -1.39988940e-01,  9.92720608e-01,  9.94910017e-04, -8.18676529e-03,
           -8.27380191e-01, -5.19487716e-02,  9.65539570e-02, -7.61224221e-01,
            1.13480280e-01, -4.74363900e-01,  9.92720655e-01, -3.95266921e-02,
           -7.05279511e-01,  9.92720619e-01,  7.88544194e-02, -5.72810080e-01,
           -9.96178846e-01, -7.33527539e-02, -3.64513523e-01, -9.96178404e-01,
            8.96837751e-02, -9.22424106e-01,  9.92721273e-01, -2.01341995e-02,
           -9.48936121e-01,  9.92673694e-01,  7.09421407e-02, -9.22984725e-01,
           -7.77647739e-01, -8.89817497e-02, -8.81365812e-01, -8.00952204e-01,
           -4.24168954e-01,  9.48456065e-01,  2.26281079e-01, -5.47380643e-01,
            9.34871196e-01,  2.84654850e-01, -4.72809331e-01,  8.73906025e-01,
           -3.16417655e-01, -5.98028521e-01,  9.14619810e-01, -2.84019840e-01,
           -3.74670206e-01,  6.59194182e-01,  7.87466848e-01, -5.97291239e-01,
            5.48831101e-01,  8.32704768e-01, -4.01358401e-01,  3.98450547e-01,
           -9.96179210e-01, -5.52822044e-01,  5.69939213e-01, -9.96179371e-01,
           -3.27678497e-01, -1.36312158e-03,  9.92720524e-01, -6.08213619e-01,
           -8.46746800e-02,  9.92720710e-01, -4.04826590e-01, -2.66671757e-02,
           -9.96178522e-01, -4.94136586e-01,  1.02650555e-01, -9.96179223e-01,
           -4.32773786e-01, -4.54216110e-01,  9.92720757e-01, -6.22501031e-01,
           -5.72956721e-01,  9.92720482e-01, -4.18435863e-01, -6.67520016e-01,
           -8.14216438e-01, -5.77885152e-01, -4.48791190e-01, -8.69043765e-01,
           -4.43969041e-01, -9.09730845e-01,  8.85853680e-01, -6.27000028e-01,
           -9.44006530e-01,  8.08179762e-01, -4.15652244e-01, -9.44736699e-01,
           -7.75386274e-01, -6.09637940e-01, -9.05362108e-01, -7.91154242e-01,
           -9.11337922e-01,  9.89926949e-01,  7.04727755e-01, -9.23980474e-01,
            9.72280938e-01,  7.49880765e-01, -8.16465004e-01,  8.57816954e-01,
           -2.59288869e-01, -8.81252173e-01,  9.39506966e-01, -1.87598986e-01,
           -8.07678094e-01,  8.68296571e-01, -7.88095454e-01, -8.53769029e-01,
            7.82659619e-01, -7.92740398e-01, -8.41033288e-01,  8.79137764e-01,
           -9.96179484e-01, -9.16156618e-01,  9.48634067e-01, -9.96141391e-01,
           -8.57292077e-01,  7.61942477e-01,  8.50526223e-01, -8.93185739e-01,
            5.41675538e-01,  7.83162808e-01, -8.51345706e-01,  5.28300747e-01,
           -3.37015058e-01, -8.95081462e-01,  6.53818102e-01, -2.43129543e-01,
           -8.83642387e-01,  7.06636810e-02,  9.92720403e-01, -9.55049147e-01,
           -8.84496092e-02,  9.92720654e-01, -8.15683760e-01, -8.08944416e-02,
           -9.27415979e-01, -9.07630131e-01,  2.33228091e-02, -8.72298354e-01,
           -8.69347622e-01, -3.62528581e-01,  9.92720687e-01, -8.74172695e-01,
           -5.14245267e-01,  9.92720542e-01, -8.18765265e-01, -6.26092572e-01,
           -8.20398514e-01, -8.67927618e-01, -4.20725528e-01, -8.46531066e-01,
           -8.58489462e-01, -8.94020319e-01,  7.71681701e-01, -9.22019902e-01,
           -9.36425206e-01,  6.78349282e-01, -8.11251535e-01, -9.66645537e-01,
            2.50131742e-01, -9.08941426e-01, -9.05332872e-01,  2.09447845e-01,
           -8.17933191e-01, -9.08671950e-01, -1.62482403e-01, -8.63086610e-01,
           -9.39088985e-01, -2.26353844e-01, -8.61034269e-01, -9.47037150e-01,
           -7.21888364e-01, -9.10712304e-01, -8.86496439e-01, -6.96305309e-01
        };

    const uIntVector *boundaryCellsResult = dc.getBoundaryIDs( );
    const floatVector *boundaryPointsResult = dc.getBoundaryPoints( );

    BOOST_CHECK( boundaryCellsResult );

    BOOST_CHECK( boundaryPointsResult );

    BOOST_CHECK( vectorTools::fuzzyEquals( *boundaryCellsResult, boundaryCellsAnswer ) );
    
    BOOST_CHECK( vectorTools::fuzzyEquals( *boundaryPointsResult, boundaryPointsAnswer ) );

}

BOOST_AUTO_TEST_CASE( testDualContouring_planes ){
    /*!
     * Test of adding planes to the dual-contouring reconstruction
     * 
     */

    floatVector points = {  3.2       , -0.5       ,  2.        ,  3.2817976 , -0.53832241,
                            1.98088073,  3.36359519, -0.57664481,  1.96176147,  3.1731361 ,
                           -0.50652963,  2.06029672,  3.2549337 , -0.54485204,  2.04117746,
                            3.3367313 , -0.58317445,  2.02205819,  3.14627221, -0.51305926,
                            2.12059345,  3.2280698 , -0.55138167,  2.10147418,  3.3098674 ,
                           -0.58970408,  2.08235492,  3.27069372, -0.36755489,  2.13373739,
                            3.35249132, -0.40587729,  2.11461812,  3.43428891, -0.4441997 ,
                            2.09549886,  3.24382982, -0.37408452,  2.19403412,  3.32562742,
                           -0.41240692,  2.17491485,  3.40742502, -0.45072933,  2.15579558,
                            3.21696593, -0.38061415,  2.25433084,  3.29876352, -0.41893656,
                            2.23521157,  3.38056112, -0.45725896,  2.21609231,  3.34138744,
                           -0.23510977,  2.26747478,  3.42318504, -0.27343218,  2.24835551,
                            3.50498263, -0.31175459,  2.22923625,  3.31452354, -0.2416394 ,
                            2.32777151,  3.39632114, -0.27996181,  2.30865224,  3.47811874,
                           -0.31828422,  2.28953297,  3.28765965, -0.24816904,  2.38806823,
                            3.36945724, -0.28649144,  2.36894896,  3.45125484, -0.32481385,
                            2.3498297};

    // Bounding plane information
    floatMatrix planePoints = { { 3.24382982, -0.37408452,  2.19403412},
                                { 3.29876352, -0.41893656,  2.23521157},
                                { 3.39632114, -0.27996181,  2.30865224} };

    floatMatrix planeNormals = { { -0.72388679,  0.64185401, -0.25300465},
                                 { -0.14037132, -0.66539974,  0.73317057},
                                 {  0.46110124,  0.83651484,  0.29602124} };

    // Set the filename
    std::string filename = "volumeReconstruction_dualContouring.yaml";

    YAML::Node yf = YAML::LoadFile( filename );

    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    errorOut error = dc.loadPoints( &points );

    BOOST_CHECK( !error );

    error = dc.addBoundingPlanes( planePoints, planeNormals );

    BOOST_CHECK( !error );

    error = dc.evaluate( );

    BOOST_CHECK( !error );

    floatVector boundaryPointAnswer = { 3.50504245, -0.33433854,  2.29253706,  3.50504243, -0.36811677,
                                        2.25480982,  3.50504244, -0.33823237,  2.18681146,  3.5050424 ,
                                       -0.36208101,  2.17586859,  3.4576968 , -0.34798511,  2.10535246,
                                        3.45611512, -0.38859457,  2.06363645,  3.47073551, -0.45838326,
                                        2.19953628,  3.4658851 , -0.49471646,  2.16388414,  3.46322261,
                                       -0.46258832,  2.09331844,  3.45045409, -0.50542083,  2.05823178,
                                        3.43994932, -0.55856558,  2.09206125,  3.43699295, -0.56559238,
                                        2.05203932,  3.44666247, -0.54534534,  2.00106873,  3.43370102,
                                       -0.56951614,  1.98034456,  3.32062494, -0.32819962,  2.14224419,
                                        3.31348903, -0.36436467,  2.18650646,  3.3163666 , -0.34522613,
                                        2.06605195,  3.31786444, -0.35991807,  2.07776709,  3.31538884,
                                       -0.4410281 ,  1.99418113,  3.33269969, -0.48471752,  1.97350221,
                                        3.20823181, -0.34609435,  2.28041505,  3.21197541, -0.36652317,
                                        2.30581811,  3.19665801, -0.45012023,  2.17758136,  3.19581603,
                                       -0.47576551,  2.1988168 ,  3.19877663, -0.46955473,  2.06875221,
                                        3.17946035, -0.50573768,  2.0980532 ,  3.21174738, -0.46186382,
                                        1.99028867,  3.18716673, -0.50158929,  2.01045077,  3.16782435,
                                       -0.56147801,  2.05887397,  3.19730134, -0.56997996,  2.0868294 ,
                                        3.18687304, -0.55755542,  1.98187028,  3.18684694, -0.58126119,
                                        2.00099249,  3.49534101, -0.29792221,  2.24057559,  3.47288857,
                                       -0.29772517,  2.29161523,  3.4860701 , -0.30265645,  2.15640509,
                                        3.47636866, -0.29612452,  2.17244197,  3.49018246, -0.39170176,
                                        2.30339375,  3.47665903, -0.41950117,  2.27128209,  3.49151897,
                                       -0.41976482,  2.21175528,  3.48862703, -0.41093981,  2.17830262,
                                        3.4096442 , -0.2941614 ,  2.24978047,  3.34731011, -0.30473436,
                                        2.30630129,  3.41702394, -0.29408351,  2.15235564,  3.37275192,
                                       -0.29231743,  2.17865449,  3.39726249, -0.32860293,  2.07587679,
                                        3.36339912, -0.31527376,  2.10423019,  3.418331  , -0.41171925,
                                        2.30024704,  3.36671035, -0.42225831,  2.27583904,  3.41013849,
                                       -0.5365973 ,  2.19206473,  3.36611126, -0.54768571,  2.16643526,
                                        3.41929866, -0.52523292,  1.98243232,  3.38900754, -0.51044604,
                                        1.99422638,  3.3962912 , -0.58976318,  2.08852087,  3.37067026,
                                       -0.58976329,  2.05402447,  3.39841057, -0.58976328,  1.99681875,
                                        3.35994258, -0.58976336,  1.97818231,  3.26938541, -0.31065916,
                                        2.26986499,  3.23732198, -0.31321102,  2.31170621,  3.28360026,
                                       -0.39946617,  2.29895107,  3.2408026 , -0.40130693,  2.28041289,
                                        3.28015359, -0.40414464,  2.15965189,  3.24085632, -0.41171876,
                                        2.22129184,  3.28584651, -0.39830163,  2.04901789,  3.25527156,
                                       -0.40636316,  2.08664042,  3.27332653, -0.42783457,  1.98238707,
                                        3.25393223, -0.42176774,  2.01096757,  3.28706332, -0.52461749,
                                        2.18731381,  3.23502297, -0.52995676,  2.16120388,  3.3031541 ,
                                       -0.58976326,  2.08350565,  3.25130052, -0.58976327,  2.0454221 ,
                                        3.2803487 , -0.58976338,  1.99730448,  3.24537211, -0.58976325,
                                        1.98083031,  3.48259   , -0.3341415 ,  2.34357671,  3.46773005,
                                       -0.3573743 ,  2.34258594,  3.50215046, -0.353256  ,  2.14241593,
                                        3.48317813, -0.32796185,  2.13479989,  3.39569452, -0.33153093,
                                        2.34714528,  3.35981889, -0.35807648,  2.33691804,  3.39300222,
                                       -0.39179603,  2.02850449,  3.35022997, -0.3585553 ,  2.03769855,
                                        3.40714685, -0.4672019 ,  2.24531741,  3.36203172, -0.48543943,
                                        2.2349984 ,  3.41714752, -0.47741371,  2.01136754,  3.36680314,
                                       -0.43778267,  2.01132224,  3.39369103, -0.57280818,  2.13341443,
                                        3.35570264, -0.57983499,  2.13196185,  3.39523303, -0.56951622,
                                        1.96170813,  3.37516094, -0.5494038 ,  1.96170808,  3.27937105,
                                       -0.33383177,  2.33791527,  3.25477307, -0.3646824 ,  2.32435629,
                                        3.27419176, -0.37193878,  2.24814641,  3.24029525, -0.34354249,
                                        2.23857382,  3.27879756, -0.44233661,  2.24698186,  3.24785638,
                                       -0.47042624,  2.22492673,  3.27619039, -0.49604568,  1.96170812,
                                        3.23114169, -0.46793065,  1.96170817,  3.27608252, -0.56005179,
                                        2.13887276,  3.24915492, -0.56997995,  2.12491295,  3.30489391,
                                       -0.56965082,  1.96170809,  3.25650239, -0.54594498,  1.96170813 };

    BOOST_CHECK( vectorTools::fuzzyEquals( *dc.getBoundaryPoints( ), boundaryPointAnswer ) );

}

BOOST_AUTO_TEST_CASE( testDualContouring_localBoundary ){
    /*!
     * Test of adding the local boundary to the dual-contouring reconstruction
     * 
     */

    floatVector points = {  3.2       , -0.5       ,  2.        ,  3.2817976 , -0.53832241,
                            1.98088073,  3.36359519, -0.57664481,  1.96176147,  3.1731361 ,
                           -0.50652963,  2.06029672,  3.2549337 , -0.54485204,  2.04117746,
                            3.3367313 , -0.58317445,  2.02205819,  3.14627221, -0.51305926,
                            2.12059345,  3.2280698 , -0.55138167,  2.10147418,  3.3098674 ,
                           -0.58970408,  2.08235492,  3.27069372, -0.36755489,  2.13373739,
                            3.35249132, -0.40587729,  2.11461812,  3.43428891, -0.4441997 ,
                            2.09549886,  3.24382982, -0.37408452,  2.19403412,  3.32562742,
                           -0.41240692,  2.17491485,  3.40742502, -0.45072933,  2.15579558,
                            3.21696593, -0.38061415,  2.25433084,  3.29876352, -0.41893656,
                            2.23521157,  3.38056112, -0.45725896,  2.21609231,  3.34138744,
                           -0.23510977,  2.26747478,  3.42318504, -0.27343218,  2.24835551,
                            3.50498263, -0.31175459,  2.22923625,  3.31452354, -0.2416394 ,
                            2.32777151,  3.39632114, -0.27996181,  2.30865224,  3.47811874,
                           -0.31828422,  2.28953297,  3.28765965, -0.24816904,  2.38806823,
                            3.36945724, -0.28649144,  2.36894896,  3.45125484, -0.32481385,
                            2.3498297};

    uIntVector global_node_ids = { 0, 1, 2, 3, 4, 5, 6, 7 };

    floatMatrix reference_nodes = { { 0, 0, 0 },
                                    { 1, 0, 0 },
                                    { 1, 1, 0 },
                                    { 0, 1, 0 },
                                    { 0, 0, 1 },
                                    { 1, 0, 1 },
                                    { 1, 1, 1 },
                                    { 0, 1, 1 } };

    floatMatrix current_nodes = { { 3.2       , -0.5       ,  2.        },
                                  { 3.9069372 ,  0.82445115,  3.3373739 },
                                  { 3.63829824,  0.75915482,  3.94034115},
                                  { 2.93136104, -0.56529632,  2.60296725},
                                  { 4.01797596, -0.88322407,  1.80880734},
                                  { 4.72491316,  0.44122707,  3.14618124},
                                  { 4.45627421,  0.37593075,  3.74914849},
                                  { 3.74933701, -0.94852039,  2.41177459} };

    std::string element_name = "Hex8";

    auto qrule = elib::default_qrules.find( element_name );

    std::unique_ptr< elib::Element > element = elib::build_element_from_string( element_name, global_node_ids, current_nodes, qrule->second );

    element->reference_nodes = reference_nodes;

    // Set the filename
    std::string filename = "volumeReconstruction_dualContouring.yaml";

    YAML::Node yf = YAML::LoadFile( filename );

    volumeReconstruction::dualContouring dc( yf );

    BOOST_CHECK( !dc.getError( ) );

    errorOut error = dc.loadPoints( &points );

    BOOST_CHECK( !error );

    error = dc.reconstructInLocalDomain( element );

    BOOST_CHECK( !error );

    error = dc.evaluate( );

    BOOST_CHECK( !error );

    floatVector boundaryPointAnswer =
        {
            3.51011836, -0.26532355,  2.36828031,  3.50658028, -0.25459436,
            2.34975918,  3.44605973, -0.22643716,  2.39243639,  3.43565531,
           -0.21337013,  2.37428845,  3.37630876, -0.19581583,  2.40912407,
            3.36759107, -0.18634552,  2.39252426,  3.5442135 , -0.26921112,
            2.30130302,  3.52281627, -0.25534824,  2.29808593,  3.47477863,
           -0.22766445,  2.32839416,  3.44727712, -0.21851907,  2.31713469,
            3.39943611, -0.19475445,  2.34642918,  3.38252767, -0.18832789,
            2.33464372,  3.43232728, -0.31796119,  2.42356113,  3.42584907,
           -0.30726524,  2.40379427,  3.36324563, -0.29044004,  2.44136913,
            3.35435013, -0.2719122 ,  2.42229264,  3.30013207, -0.25741441,
            2.44798999,  3.29692184, -0.2515022 ,  2.41585128,  3.55359801,
           -0.35444114,  2.3249919 ,  3.53505203, -0.33729285,  2.32061244,
            3.37429684, -0.37286478,  2.36276326,  3.40687148, -0.37454238,
            2.34016904,  3.29566998, -0.33179759,  2.37834318,  3.3304426 ,
           -0.34393915,  2.35984492,  3.24279017, -0.31330412,  2.38829517,
            3.265695  , -0.31236286,  2.36984565,  3.48727937, -0.39877834,
            2.25062218,  3.53055736, -0.4202009 ,  2.21687121,  3.36939079,
           -0.43637498,  2.30986538,  3.35480044, -0.4232575 ,  2.29219314,
            3.29590721, -0.40789019,  2.31738173,  3.28135103, -0.38725307,
            2.28855759,  3.23446366, -0.37480333,  2.30908538,  3.22904867,
           -0.37308012,  2.27875426,  3.48057687, -0.47268178,  2.21409646,
            3.46485214, -0.45797892,  2.19671511,  3.31045091, -0.49279839,
            2.23734037,  3.35186959, -0.49788518,  2.21835779,  3.2294351 ,
           -0.45421989,  2.25383853,  3.26480166, -0.46763301,  2.23570205,
            3.18300788, -0.42625822,  2.2701923 ,  3.20431637, -0.42921425,
            2.2499627 ,  3.42559194, -0.52204109,  2.13360354,  3.46091519,
           -0.5401211 ,  2.09692467,  3.3118167 , -0.55878799,  2.17911221,
            3.29049966, -0.54464034,  2.16898817,  3.22699497, -0.52769841,
            2.18703753,  3.20255574, -0.50757317,  2.18032159,  3.15730868,
           -0.49766108,  2.19237228,  3.15342727, -0.49333995,  2.17653681,
            3.41488876, -0.59592195,  2.09253983,  3.39584641, -0.58208147,
            2.06989987,  3.26563629, -0.58550226,  2.13923786,  3.29798313,
           -0.58996618,  2.10146683,  3.19327883, -0.5482546 ,  2.14668969,
            3.22862591, -0.55637832,  2.11458736,  3.13700251, -0.52387081,
            2.16544418,  3.1627192 , -0.52560403,  2.13028419,  3.36971967,
           -0.60893788,  2.04333687,  3.40681209, -0.62309811,  2.02557418,
            3.29274378, -0.57316637,  2.06215401,  3.32120861, -0.58079141,
            2.03936271,  3.23118489, -0.54017823,  2.0648221 ,  3.26295593,
           -0.55250236,  2.05015978,  3.17219979, -0.51473819,  2.08481121,
            3.20531841, -0.52118867,  2.05145143,  3.31839381, -0.56818233,
            2.00811564,  3.35231785, -0.57710075,  1.98047563,  3.25713608,
           -0.53834954,  2.01923121,  3.28410914, -0.54640447,  1.99997837,
            3.20750051, -0.51267474,  2.02399392,  3.22754715, -0.51751148,
            2.00643601,  3.46153562, -0.27006201,  2.39398874,  3.47334253,
           -0.29472871,  2.38064639,  3.39320886, -0.23117505,  2.41502891,
            3.40514666, -0.2640707 ,  2.39806013,  3.33270642, -0.21081214,
            2.41936581,  3.3392021 , -0.23310433,  2.40712271,  3.54155873,
           -0.28006344,  2.26494976,  3.51236297, -0.27574808,  2.26327654,
            3.466119  , -0.24188751,  2.28515325,  3.42870841, -0.23601115,
            2.28332403,  3.39295808, -0.20621615,  2.30352002,  3.36788882,
           -0.21246058,  2.29305365,  3.38237177, -0.33521838,  2.42572647,
            3.39939251, -0.36504865,  2.39339429,  3.30985555, -0.29681056,
            2.44900322,  3.32250348, -0.33089295,  2.41359371,  3.25928924,
           -0.27778814,  2.45026547,  3.26011722, -0.29865645,  2.43220222,
            3.50254955, -0.35186833,  2.33453133,  3.52665454, -0.39060937,
            2.29290281,  3.56094914, -0.3776502 ,  2.26643258,  3.5200506 ,
           -0.38412418,  2.23688942,  3.50604875, -0.32476268,  2.21758375,
            3.46108532, -0.34348982,  2.19197898,  3.42998894, -0.28859887,
            2.23584255,  3.38690166, -0.30940215,  2.20871264,  3.36055186,
           -0.26305886,  2.24572775,  3.33517398, -0.28353876,  2.22228167,
            3.37827965, -0.38039848,  2.30585538,  3.39581202, -0.40315947,
            2.28140439,  3.30744947, -0.34507765,  2.3192153 ,  3.32271189,
           -0.37985939,  2.28781998,  3.25369177, -0.33017489,  2.32175388,
            3.26067362, -0.35074812,  2.28962179,  3.46779196, -0.39984364,
            2.14211826,  3.4457554 , -0.40052427,  2.13728155,  3.39596148,
           -0.36497373,  2.16001236,  3.3636314 , -0.35923228,  2.15903339,
            3.3291937 , -0.32414059,  2.18428759,  3.30650459, -0.32979574,
            2.1748117 ,  3.32143108, -0.45605355,  2.3038486 ,  3.33678197,
           -0.48430847,  2.26971057,  3.24879995, -0.41608446,  2.32082024,
            3.25973346, -0.45113035,  2.28961487,  3.19741299, -0.39557397,
            2.32236279,  3.1971931 , -0.41621798,  2.3052735 ,  3.44021461,
           -0.47273853,  2.22689048,  3.46354428, -0.51426958,  2.18281357,
            3.49363437, -0.5008243 ,  2.13503688,  3.45257068, -0.50686571,
            2.10603664,  3.44270415, -0.44806158,  2.09355419,  3.3961715 ,
           -0.46399712,  2.07018255,  3.36689268, -0.41071836,  2.11293093,
            3.32461982, -0.42956815,  2.08772994,  3.29866943, -0.38373722,
            2.12435652,  3.27252403, -0.40332404,  2.10157392,  3.31920227,
           -0.5049494 ,  2.19213013,  3.33915021, -0.52831958,  2.16856411,
            3.23313056, -0.46891532,  2.21042759,  3.24611249, -0.50114492,
            2.18583826,  3.19224477, -0.4447345 ,  2.21047407,  3.1856924 ,
           -0.47102245,  2.18907544,  3.40336819, -0.51930228,  2.02136741,
            3.38280289, -0.52029279,  2.01653119,  3.33333334, -0.48383743,
            2.04014535,  3.30057189, -0.47722512,  2.0398742 ,  3.26571298,
           -0.44254266,  2.0646762 ,  3.24141431, -0.44880301,  2.054343  ,
            3.26972741, -0.56960534,  2.18886531,  3.28772084, -0.59116116,
            2.17570189,  3.19088038, -0.52806344,  2.19279986,  3.21229579,
           -0.55582271,  2.18450361,  3.13701426, -0.50848498,  2.19643429,
            3.14614673, -0.5187621 ,  2.19014593,  3.37690643, -0.58732176,
            2.10384076,  3.40153911, -0.61981151,  2.07090676,  3.42319638,
           -0.61152296,  2.02411924,  3.38941752, -0.59931715,  2.01868358,
            3.39357139, -0.55442051,  1.98761984,  3.35533851, -0.55346559,
            1.98116648,  3.31180856, -0.51272556,  2.00980668,  3.27450322,
           -0.51458286,  2.00097857,  3.23896789, -0.48534276,  2.02071289,
            3.21898487, -0.49632632,  2.0069175 ,  3.52192528, -0.28999025,
            2.35493796,  3.5271388 , -0.30706495,  2.32594827,  3.33322414,
           -0.20230392,  2.3739012 ,  3.31089434, -0.2259143 ,  2.37789244,
            3.55296511, -0.29701502,  2.29083446,  3.55079984, -0.30531093,
            2.27118903,  3.35745841, -0.19457232,  2.32417735,  3.33999992,
           -0.21695728,  2.31733975,  3.44934802, -0.34779146,  2.39122894,
            3.44831732, -0.36001213,  2.34995726,  3.26499559, -0.26972636,
            2.39806831,  3.24196219, -0.29243581,  2.40635842,  3.577703  ,
           -0.39318218,  2.28336338,  3.5714559 , -0.41372692,  2.24641436,
            3.30055735, -0.26008645,  2.32077599,  3.27720512, -0.27968425,
            2.33787245,  3.51774258, -0.34334245,  2.24679137,  3.50218725,
           -0.36710955,  2.19485931,  3.32505245, -0.25532544,  2.26276997,
            3.30793967, -0.27484423,  2.26221497,  3.4392609 , -0.39732854,
            2.2580115 ,  3.44771345, -0.42396064,  2.20923623,  3.26546364,
           -0.31916586,  2.2713048 ,  3.24207887, -0.34788831,  2.26534563,
            3.46966984, -0.41310082,  2.16993308,  3.48023799, -0.42710778,
            2.14279679,  3.29537209, -0.30871269,  2.21225079,  3.2762483 ,
           -0.33343653,  2.20424416,  3.38474168, -0.4646299 ,  2.27572735,
            3.39037665, -0.4797565 ,  2.23493513,  3.20500539, -0.3866486 ,
            2.27308672,  3.18322778, -0.40561421,  2.28728159,  3.50390653,
           -0.51421283,  2.17001955,  3.50197888, -0.53407969,  2.12592491,
            3.23885869, -0.3765921 ,  2.20173986,  3.21942192, -0.39181208,
            2.21873954,  3.45103913, -0.46871763,  2.11530632,  3.43554185,
           -0.48899081,  2.07090043,  3.2612354 , -0.3741064 ,  2.14402247,
            3.24723056, -0.38853951,  2.14747733,  3.38370316, -0.5220089 ,
            2.13753708,  3.37892318, -0.54745119,  2.08256847,  3.20770363,
           -0.43212314,  2.15000136,  3.18099561, -0.46312392,  2.14777063,
            3.39697016, -0.53450116,  2.04334101,  3.41132915, -0.54887709,
            2.01717963,  3.23422053, -0.4240175 ,  2.09469557,  3.21413038,
           -0.45238115,  2.08021992,  3.32981013, -0.58034382,  2.16594879,
            3.33198693, -0.58591557,  2.12343005,  3.14185005, -0.5020447 ,
            2.16103228,  3.12787003, -0.51359368,  2.17173254,  3.43952143,
           -0.62841171,  2.05960584,  3.44059095, -0.63530392,  2.03100984,
            3.17564464, -0.4950109 ,  2.08218598,  3.16046311, -0.50145599,
            2.10927093,  3.38774319, -0.57350276,  2.01165257,  3.39055073,
           -0.57805567,  1.986929  ,  3.20134114, -0.48815022,  2.02607196,
            3.18956343, -0.49977917,  2.03026875
        };

    BOOST_CHECK( vectorTools::fuzzyEquals( *dc.getBoundaryPoints( ), boundaryPointAnswer ) );

    floatType volumeAnswer =  0.0184456;

    floatVector volumeResult;

    error = dc.performVolumeIntegration( floatVector( points.size( ) / 3, 1 ), 1, volumeResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( volumeResult[ 0 ], volumeAnswer ) );

    floatType areaAnswer =  0.585814;

    floatVector areaResult;

    error = dc.performSurfaceIntegration( floatVector( points.size( ) / 3, 1 ), 1, areaResult );

    BOOST_CHECK( !error );

    BOOST_CHECK( vectorTools::fuzzyEquals( areaResult[ 0 ], areaAnswer ) );

}
