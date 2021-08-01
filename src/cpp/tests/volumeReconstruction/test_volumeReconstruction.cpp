//!The test file for volumeReconstruction.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>
#define USE_EIGEN
#include<vector_tools.h>

#include<volumeReconstruction.h>

typedef volumeReconstruction::errorNode errorNode; //!Redefinition for the error node
typedef volumeReconstruction::errorOut errorOut; //!Redefinition for a pointer to the error node
typedef volumeReconstruction::floatType floatType; //!Define the float values type.
typedef volumeReconstruction::floatVector floatVector; //! Define a vector of floats
typedef volumeReconstruction::floatMatrix floatMatrix; //!Define a matrix of floats
typedef volumeReconstruction::intMatrix intMatrix; //!Define a matrix of ints
typedef volumeReconstruction::uIntType uIntType; //!Define the unsigned int type
typedef volumeReconstruction::uIntVector uIntVector; //!Define a vector of unsigned ints

int test_dualContouring_constructor( std::ofstream &results ){
    /*!
     * Tests for the constructors of dualContouring
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError( ) ){
        
        dc.getError( )->print( );

        results << "test_dualContouring_constructor (test 1) & False\n";
        return 1;
    
    }

    std::shared_ptr< volumeReconstruction::volumeReconstructionBase > yR
        = volumeReconstruction::volumeReconstructionBase( yf ).create( );

    if ( yR->getError( ) ){

        yR->getError( )->print( );

        results << "test_dualContouring_constructor (test 2) & False\n";
        return 1;

    }

   results << "test_dualContouring_constructor & True\n";
   return 0;
}

int test_dualContouring_loadPoints( std::ofstream &results ){
    /*!
     * Test for loading the points into the object
     *
     * :param std::ofstream &results: The output file
     */

    floatVector points = { 1, 2, 3, 4, 5, 6 };

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_loadPoints & False\n";
        return 1;

    }

    std::unique_ptr< errorNode > error;
    error.reset( dc.loadPoints( &points ) );

    if ( error ){

        error->print( );
        results << "test_dualContouring_loadPoints (test 1) & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( *dc.getPoints( ), points ) ){

        results << "test_dualContouring_loadPoints (test 2) & False\n";
        return 1;

    }

    floatVector points2 = { 1, 2, 3, 4, 5 };


    error.reset( dc.loadPoints( &points2 ) );

    if ( !error ){

        error->print( );
        results << "test_dualContouring_loadPoints (test 3) & False\n";
        return 1;

    }

    results << "test_dualContouring_loadPoints & True\n";
    return 0;
}

int test_dualContouring_loadFunction( std::ofstream &results ){
    /*!
     * Test for loading the function into the object
     *
     * :param std::ofstream &results: The output file
     */

    floatVector points = { 1, 2, 3, 4, 5, 6 };

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_loadFunction & False\n";
        return 1;

    }

    std::unique_ptr< errorNode > error;
    error.reset( dc.loadPoints( &points ) );

    if ( error ){

        error->print( );
        results << "test_dualContouring_loadFunction & False\n";
        return 1;

    }

    floatVector function = { -1, 10 };

    error.reset( dc.loadFunction( &function ) );

    if ( error ){

        error->print( );
        results << "test_dualContouring_loadFunction & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( *dc.getFunction( ), function ) ){

        results << "test_dualContouring_loadFunction (test 1) & False\n";
        return 1;

    }

    floatVector function2 = { 2 };

    error.reset( dc.loadFunction( &function2 ) );
    if ( !error ){

        results << "test_dualContouring_loadFunction (test 2) & False\n";
        return 1;

    }

    results << "test_dualContouring_loadFunction & True\n";
    return 0;
}

int test_dualContouring_getFunctionValue( std::ofstream &results ){
    /*!
     * Test for getting a particular value of the function
     *
     * :param std::ofstream &results: The output file
     */

    floatVector points = { 1, 2, 3, 4, 5, 6 };

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_getFunctionValue & False\n";
        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){

        error->print( );
        results << "test_dualContouring_getFunctionValue & False\n";
        return 1;

    }

    floatVector function = { -1, 10 };

    error = dc.loadFunction( &function );

    if ( error ){

        error->print( );
        results << "test_dualContouring_getFunctionValue & False\n";
        return 1;

    }

    floatType result = 0;
    error = dc.getFunctionValue( 0, result );

    if ( error ){

        error->print( );
        results << "test_dualContouring_getFunctionValue & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( function [ 0 ], result ) ){

        results << "test_dualContouring_getFunctionValue (test 1) & False\n";
        return 1;

    }

    volumeReconstruction::dualContouring dc2( yf );

    if ( dc2.getError ( ) ){

        dc2.getError( )->print( );

        results << "test_dualContouring_getFunctionValue & False\n";
        return 1;

    }

    error = dc2.loadPoints( &points );

    if ( error ){

        error->print( );
        results << "test_dualContouring_getFunctionValue & False\n";
        return 1;

    }

    result = 0;
    error = dc2.getFunctionValue( 0, result );

    if ( error ){

        error->print( );
        results << "test_dualContouring_getFunctionValue & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( 1., result ) ){

        results << "test_dualContouring_getFunctionValue (test 2) & False\n";
        return 1;

    }

    results << "test_dualContouring_getFunctionValue & True\n";
    return 0;
}

int test_KDNode_constructor( std::ofstream &results ){
    /*!
     * Test the KD tree creation
     *
     * :param std::ofstream &results: The output file
     */

    floatVector points = { 2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2 };
    uIntVector  ownedIndices = { 0, 2, 4, 6, 8, 10 };
    unsigned int dim = 2;

    volumeReconstruction::KDNode tree( &points, ownedIndices, 0, dim );

    results << "test_KDNode_constructor & True\n";
    return 0;
}

int test_KDNode_getIndex( std::ofstream &results ){
    /*!
     * Test getting the index of the KDNode
     *
     * :param std::ofstream &results: The output file
     */

    floatVector points = { 2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2 };
    uIntVector  ownedIndices = { 0, 2, 4, 6, 8, 10 };
    unsigned int dim = 2;

    volumeReconstruction::KDNode tree( &points, ownedIndices, 0, dim );

    if ( !( tree.getIndex( ), 10 ) ){
        
        results << "test_KDNode_getIndex (test 1) & True\n";
        return 1;

    }

    results << "test_KDNode_getIndex & True\n";
    return 0;
}

int test_KDNode_getMinimumValueDimension( std::ofstream &results ){
    /*!
     * Test getting the minimum of a dimension of the KD tree
     *
     * :param std::ofstream &results: The output file
     */

    floatVector points = { 2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2 };
    uIntVector  ownedIndices = { 0, 2, 4, 6, 8, 10 };
    unsigned int dim = 2;

    volumeReconstruction::KDNode tree( &points, ownedIndices, 0, dim );

    floatType answer1 = 2;
    floatType answer2 = 1;

    floatType result = tree.getMinimumValueDimension( 0 );

    if ( !vectorTools::fuzzyEquals( result, answer1 ) ){
        results << "test_KDNode_getMinimumValueDimension (test 1) & False\n";
        return 1;
    }

    result = tree.getMinimumValueDimension( 1 );

    if ( !vectorTools::fuzzyEquals( result, answer2 ) ){
        results << "test_KDNode_getMinimumValueDimesnion (test 2) & False\n";
        return 1;
    }

    floatVector points2 = { 1, 10, 50, 50, 10, 30, 35, 90, 55, 1, 60, 80, 25, 40, 70, 70, 51, 75 };
    uIntVector ownedIndices2 = { 0, 2, 4, 6, 8, 10, 12, 14, 16 };

    volumeReconstruction::KDNode tree2( &points2, ownedIndices2, 0, dim );

    floatType answer3 = 1;
    floatType answer4 = 1;

    result = tree2.getMinimumValueDimension( 0 );

    if ( !vectorTools::fuzzyEquals( result, answer3 ) ){
        results << "test_KDNode_getMinimumValueDimension (test 3) & False\n";
        return 1;
    }

    result = tree2.getMinimumValueDimension( 1 );

    if ( !vectorTools::fuzzyEquals( result, answer4 ) ){
        results << "test_KDNode_getMinimumValueDimension (test 4) & False\n";
        return 1;
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
            0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782
        };

    uIntVector ownedIndices3 = { 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57 };

    volumeReconstruction::KDNode tree3( &points3, ownedIndices3, 0, 3 );

    floatType answer5 = -0.72313702;
    floatType answer6 = -0.78823811;
    floatType answer7 = -0.99593042;

    if ( !vectorTools::fuzzyEquals( tree3.getMinimumValueDimension( 0 ), answer5 ) ){
        results << "test_KDNode_getMinimumValueDimension (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( tree3.getMinimumValueDimension( 1 ), answer6 ) ){
        results << "test_KDNode_getMinimumValueDimension (test 6) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( tree3.getMinimumValueDimension( 2 ), answer7 ) ){
        results << "test_KDNode_getMinimumValueDimension (test 7) & False\n";
        return 1;
    }

    results << "test_KDNode_getMinimumValueDimension & True\n";
    return 0;
}

int test_KDNode_getMaximumValueDimension( std::ofstream &results ){
    /*!
     * Test getting the maximum of a dimension of the KD tree
     *
     * :param std::ofstream &results: The output file
     */

    floatVector points = { 2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2 };
    uIntVector  ownedIndices = { 0, 2, 4, 6, 8, 10 };
    unsigned int dim = 2;

    volumeReconstruction::KDNode tree( &points, ownedIndices, 0, dim );

    floatType answer1 = 9;
    floatType answer2 = 7;

    floatType result = tree.getMaximumValueDimension( 0 );

    if ( !vectorTools::fuzzyEquals( result, answer1 ) ){
        std::cout << "result: ";
        std::cout << result << "\n";
        results << "test_KDNode_getMaximumValueDimension (test 1) & False\n";
        return 1;
    }

    result = tree.getMaximumValueDimension( 1 );

    if ( !vectorTools::fuzzyEquals( result, answer2 ) ){
        results << "test_KDNode_getMaximumValueDimesnion (test 2) & False\n";
        return 1;
    }

    floatVector points2 = { 1, 10, 50, 50, 10, 30, 35, 90, 55, 1, 60, 80, 25, 40, 70, 70, 51, 75 };
    uIntVector ownedIndices2 = { 0, 2, 4, 6, 8, 10, 12, 14, 16 };

    volumeReconstruction::KDNode tree2( &points2, ownedIndices2, 0, dim );

    floatType answer3 = 70;
    floatType answer4 = 90;

    result = tree2.getMaximumValueDimension( 0 );

    if ( !vectorTools::fuzzyEquals( result, answer3 ) ){
        results << "test_KDNode_getMaximumValueDimension (test 3) & False\n";
        return 1;
    }

    result = tree2.getMaximumValueDimension( 1 );

    if ( !vectorTools::fuzzyEquals( result, answer4 ) ){
        results << "test_KDNode_getMaximumValueDimension (test 4) & False\n";
        return 1;
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
            0.37676147,  0.16014946,  0.85054661,  0.82638332, -0.63426782
        };

    uIntVector ownedIndices3 = { 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57 };

    volumeReconstruction::KDNode tree3( &points3, ownedIndices3, 0, 3 );

    floatType answer5 = 0.92378245;
    floatType answer6 = 0.83431932;
    floatType answer7 = 0.91994519;

    if ( !vectorTools::fuzzyEquals( tree3.getMaximumValueDimension( 0 ), answer5 ) ){
        results << "test_KDNode_getMaximumValueDimension (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( tree3.getMaximumValueDimension( 1 ), answer6 ) ){
        results << "test_KDNode_getMaximumValueDimension (test 6) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( tree3.getMaximumValueDimension( 2 ), answer7 ) ){
        results << "test_KDNode_getMaximumValueDimension (test 7) & False\n";
        return 1;
    }

    results << "test_KDNode_getMaximumValueDimension & True\n";
    return 0;
}

int test_KDNode_getPointsInRange( std::ofstream &results ){
    /*!
     * Get all of the points in a range
     *
     * :param std::ofstream &results: The output file
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

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_KDNode_getPointsInRange (test 1) & False\n";
        return 1;
    }

    floatVector points2 = { 1, 10, 50, 50, 10, 30, 35, 90, 55, 1, 60, 80, 25, 40, 70, 70, 51, 75 };
    uIntVector ownedIndices2 = { 0, 2, 4, 6, 8, 10, 12, 14, 16 };

    volumeReconstruction::KDNode tree2( &points2, ownedIndices2, 0, dim );

    answer = { 2, 12, 16 };

    upperBound = { 52, 76 };
    lowerBound = { 22,  5 };

    result.clear( );

    tree2.getPointsInRange( upperBound, lowerBound, result );

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_KDNode_getPointsInRange (test 2) & False\n";
        return 1;
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

        if ( !std::any_of( result.begin( ), result.end( ),
                           [&]( uIntType r ){ return vectorTools::fuzzyEquals( r, *a ); } ) ){
            vectorTools::print( result );
            vectorTools::print( answer );
            results << "test_KDNode_getPointsInRange (test 3) & False\n";
            return 1;
        }

    }

    results << "test_KDNode_getPointsInRange & True\n";
    return 0;
}

int test_KDNode_getPointsWithinRadiusOfOrigin( std::ofstream &results ){
    /*!
     * Get all of the points within a given radius of the origin
     *
     * :param std::ofstream &results: The output file
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

        if ( !std::any_of( result.begin( ), result.end( ),
                           [&]( uIntType r ){ return vectorTools::fuzzyEquals( r, *a ); } ) ){
            vectorTools::print( result );
            vectorTools::print( answer );
            results << "test_KDNode_getPointsWithinRadiusOfOrigin (test 1) & False\n";
            return 1;
        }

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

        if ( !std::any_of( result.begin( ), result.end( ),
                           [&]( uIntType r ){ return vectorTools::fuzzyEquals( r, *a ); } ) ){
            vectorTools::print( result );
            vectorTools::print( answer );
            results << "test_KDNode_getPointsWithinRadiusOfOrigin (test 2) & False\n";
            return 1;
        }

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

        if ( !std::any_of( result.begin( ), result.end( ),
                           [&]( uIntType r ){ return vectorTools::fuzzyEquals( r, *a ); } ) ){
            vectorTools::print( result );
            vectorTools::print( answer );
            results << "test_KDNode_getPointsWithinRadiusOfOrigin (test 3) & False\n";
            return 1;
        }

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

        if ( !std::any_of( result.begin( ), result.end( ),
                           [&]( uIntType r ){ return vectorTools::fuzzyEquals( r, *a ); } ) ){
            vectorTools::print( result );
            vectorTools::print( answer );
            results << "test_KDNode_getPointsWithinRadiusOfOrigin (test 4) & False\n";
            return 1;
        }

    }

    results << "test_KDNode_getPointsWithinRadiusOfOrigin & True\n";
    return 0;
}

int test_dualContouring_evaluate( std::ofstream &results ){
    /*!
     * Test the dualContouring evaluate function. This prepares the 
     * volume reconstruction object to perform the other functions
     *
     * :param std::ofstream &results: The output file
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

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_evaluate & False\n";
        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){
        error->print( );
        results << "test_dualContouring_evaluate & False\n";
        return 1;
    }

    error = dc.evaluate( );

    if ( error ){
        error->print( );
        results << "test_dualContouring_evaluate & False\n";
        return 1;
    }

    const floatVector upperBoundsAnswer = { 0.98173517, 0.99606496, 0.99247162 };

    const floatVector lowerBoundsAnswer = { -0.95963138, -0.98616276, -0.99593042 };

    const floatType medianNeighborhoodDistanceAnswer = 0.398116;

    const floatVector* upperBoundsResult = dc.getUpperBounds( );
    const floatVector* lowerBoundsResult = dc.getLowerBounds( );
    const floatType*   medianNeighborhoodDistanceResult = dc.getMedianNeighborhoodDistance( );

    if ( !vectorTools::fuzzyEquals( *upperBoundsResult, upperBoundsAnswer ) ){
        vectorTools::print( *upperBoundsResult );
        results << "test_dualContouring_evaluate (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( *lowerBoundsResult, lowerBoundsAnswer ) ){
        vectorTools::print( *lowerBoundsResult );
        results << "test_dualContouring_evaluate (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( *medianNeighborhoodDistanceResult, medianNeighborhoodDistanceAnswer ) ){
        std::cout << *medianNeighborhoodDistanceResult << "\n";
        results << "test__dualContouring_evaluate (test 3) & False\n";
        return 1;
    }

    results << "test_dualContouring_evaluate & True\n";
    return 0;
}

int test_dualContouringInternalPointResidual( std::ofstream &results ){
    /*!
     * Test the dual contouring internal point residual
     *
     * :param std::ofstream &results: The output file
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

    if ( error ){

        error->print( );
        results << "test_dualContouringInternalPointResidual & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( residualResult, residualAnswer ) ){
        results << "test_dualContouringInternalPointResidual (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( vectorTools::appendVectors( jacobian ), jacobianAnswerVec ) ){
        results << "test_dualContouringInternalPointResidual (test 2) & False\n";
        return 1;
    }

    results << "test_dualContouringInternalPointResidual & True\n";
    return 0;
}

int test_dualContouring_performVolumeIntegration( std::ofstream &results ){
    /*!
     * Test volume integration over the reconstructed domain
     *
     * :param std::ofstream &results: The output file
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

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_performVolumeIntegration & False\n";
        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){
        error->print( );
        results << "test_dualContouring_performVolumeIntegration & False\n";
        return 1;
    }

    floatVector functionValues( points.size( ) );
    for ( unsigned int i = 0; i < points.size( ); i+=3 ){

        functionValues[ i + 0 ] = 1;
        functionValues[ i + 1 ] = 2;
        functionValues[ i + 2 ] = 3;

    }

    floatVector integratedVolumeResult;
    floatVector integratedVolumeAnswer = { 7.4237149, 14.8474298, 22.2711447 };

    error = dc.performVolumeIntegration( functionValues, 3, integratedVolumeResult );

    if ( error ){
        error->print( );
        results << "test_dualContouring_performVolumeIntegration & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( integratedVolumeResult, integratedVolumeAnswer ) ){
        vectorTools::print( integratedVolumeResult );
        vectorTools::print( integratedVolumeResult - integratedVolumeAnswer );
        results << "test_dualContouring_performVolumeIntegration (test 1) & False\n";
        return 1;
    }

    results << "test_dualContouring_performVolumeIntegration & True\n";
    return 0;
    
}

int test_dualContouring_performRelativePositionVolumeIntegration( std::ofstream &results ){
    /*!
     * Test volume integration over the reconstructed domain utilizing the
     * relative position to compute a dyadic product.
     *
     * :param std::ofstream &results: The output file
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

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_performRelativePositionVolumeIntegration & False\n";
        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){
        error->print( );
        results << "test_dualContouring_performRelativePositionVolumeIntegration & False\n";
        return 1;
    }

    floatVector functionValues( points.size( ) );
    for ( unsigned int i = 0; i < points.size( ); i+=3 ){

        functionValues[ i + 0 ] = 1;
        functionValues[ i + 1 ] = 2;
        functionValues[ i + 2 ] = 3;

    }

    floatVector integratedVolumeResult;
    floatVector integratedVolumeAnswer = { 0.133381, -0.073798,  0.456491,
                                           0.266761, -0.147596,  0.91298 ,
                                           0.400142, -0.221394,  1.369474 };

    floatVector origin = { 0., 0., 0. };

    error = dc.performRelativePositionVolumeIntegration( functionValues, 3, origin, integratedVolumeResult );

    if ( error ){
        error->print( );
        results << "test_dualContouring_performRelativePositionVolumeIntegration & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( integratedVolumeResult, integratedVolumeAnswer ) ){
        vectorTools::print( integratedVolumeResult );
        vectorTools::print( integratedVolumeResult - integratedVolumeAnswer );
        results << "test_dualContouring_performRelativePositionVolumeIntegration (test 1) & False\n";
        return 1;
    }

    results << "test_dualContouring_performRelativePositionVolumeIntegration & True\n";
    return 0;
    
}

int test_dualContouring_performSurfaceIntegration( std::ofstream &results ){
    /*!
     * Test surface integration over the reconstructed domain
     *
     * :param std::ofstream &results: The output file
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

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_performSurfaceIntegration & False\n";
        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){
        error->print( );
        results << "test_dualContouring_performSurfaceIntegration & False\n";
        return 1;
    }

    floatVector functionValues( points.size( ) );
    for ( unsigned int i = 0; i < points.size( ); i+=3 ){

        functionValues[ i + 0 ] = 1;
        functionValues[ i + 1 ] = 2;
        functionValues[ i + 2 ] = 3;

    }

    floatVector integratedSurfaceResult;
    floatVector integratedSurfaceAnswer = { 23.4326, 46.8652, 70.2978 };

    error = dc.performSurfaceIntegration( functionValues, 3, integratedSurfaceResult );

    if ( error ){
        error->print( );
        results << "test_dualContouring_performSurfaceIntegration & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( integratedSurfaceResult, integratedSurfaceAnswer ) ){
        vectorTools::print( integratedSurfaceResult );
        vectorTools::print( integratedSurfaceResult - integratedSurfaceAnswer );
        results << "test_dualContouring_performSurfaceIntegration (test 1) & False\n";
        return 1;
    }

    results << "test_dualContouring_performSurfaceIntegration & True\n";
    return 0;
    
}

int test_dualContouring_performPositionWeightedSurfaceIntegration( std::ofstream &results ){
    /*!
     * Test position weighted surface integration over the reconstructed domain
     *
     * :param std::ofstream &results: The output file
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

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_performPositionWeightedSurfaceIntegration & False\n";
        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){
        error->print( );
        results << "test_dualContouring_performPositionWeightedSurfaceIntegration & False\n";
        return 1;
    }

    floatVector functionValues( points.size( ) / 3 );
    for ( unsigned int i = 0; i < points.size( ) / 3; i++ ){

        functionValues[ i + 0 ] = 1;

    }

    floatVector integratedSurfaceResult;
    floatVector integratedSurfaceAnswer = { 0.257597, 0.0470696, -0.109122 };

    error = dc.performPositionWeightedSurfaceIntegration( functionValues, 1, integratedSurfaceResult );

    if ( error ){
        error->print( );
        results << "test_dualContouring_performPositionWeightedSurfaceIntegration & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( integratedSurfaceResult, integratedSurfaceAnswer ) ){
        vectorTools::print( integratedSurfaceResult );
        vectorTools::print( integratedSurfaceResult - integratedSurfaceAnswer );
        assert( 1 == 0 );
        results << "test_dualContouring_performPositionWeightedSurfaceIntegration (test 1) & False\n";
        return 1;
    }

    results << "test_dualContouring_performPositionWeightedSurfaceIntegration & True\n";
    return 0;
    
}


int test_dualContouring_performSurfaceFluxIntegration( std::ofstream &results ){
    /*!
     * Test surface flux integration over the reconstructed domain
     *
     * :param std::ofstream &results: The output file
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

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_performSurfaceFluxIntegration & False\n";
        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){
        error->print( );
        results << "test_dualContouring_performSurfaceFluxIntegration & False\n";
        return 1;
    }

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

    if ( error ){
        error->print( );
        results << "test_dualContouring_performSurfaceFluxIntegration & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( integratedSurfaceResult, integratedSurfaceAnswer ) ){
        vectorTools::print( integratedSurfaceResult );
        vectorTools::print( integratedSurfaceResult - integratedSurfaceAnswer );
        results << "test_dualContouring_performSurfaceFluxIntegration (test 1) & False\n";
        return 1;
    }

    results << "test_dualContouring_performSurfaceFluxIntegration & True\n";
    return 0;
    
}

int test_dualContouring_performRelativePositionSurfaceFluxIntegration( std::ofstream &results ){
    /*!
     * Test surface flux integration over the reconstructed domain
     *
     * :param std::ofstream &results: The output file
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

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_performRelativePositionSurfaceFluxIntegration & False\n";
        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){
        error->print( );
        results << "test_dualContouring_performRelativePositionSurfaceFluxIntegration & False\n";
        return 1;
    }

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
    floatVector integratedSurfaceAnswer = { 7.48955819, 22.47199508, 37.44446021, 14.97911639, 29.96155038, 44.93402134 };
    floatVector origin = { 0., 0., 0. };

    error = dc.performRelativePositionSurfaceFluxIntegration( functionValues, 6, origin, integratedSurfaceResult );

    if ( error ){
        error->print( );
        results << "test_dualContouring_performSurfaceFluxIntegration & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( integratedSurfaceResult, integratedSurfaceAnswer ) ){
        vectorTools::print( integratedSurfaceResult );
        vectorTools::print( integratedSurfaceResult - integratedSurfaceAnswer );
        results << "test_dualContouring_performRelativePositionSurfaceFluxIntegration (test 1) & False\n";
        return 1;
    }

    results << "test_dualContouring_performRelativePositionSurfaceFluxIntegration & True\n";
    return 0;
    
}

int test_dualContouring_exportConfiguration( std::ofstream &results ){
    /*!
     * Test the export of the configuration file
     *
     * :param std::ofstream &results: The output file
     */

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError( ) ){
        
        dc.getError( )->print( );

        results << "test_dualContouring_exportConfiguration & False\n";
        return 1;
    
    }

    YAML::Node config = dc.exportConfiguration( );

    if ( !config[ "type" ] ){

        results << "test_dualContouring_exportConfiguration (test 1) & False\n";
        return 1;

    }
    else if ( config[ "type" ].as< std::string >( ).compare( "dual_contouring" ) != 0 ){

        results << "test_dualContouring_exportConfiguration (test 2) & False\n";
        return 1;

    }

    results << "test_dualContouring_exportConfiguration & True\n";
    return 0;
    
}

int test_dualContouring_getSurfaceSubdomains( std::ofstream &results ){
    /*!
     * Test getting sub-domains on the surface of the body
     *
     * :param std::ofstream &results: The output file
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

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_getSurfaceSubdomains & False\n";
        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){
        error->print( );
        results << "test_dualContouring_getSurfaceSubdomains & False\n";
        return 1;
    }

    uIntVector subdomainNodeCountAnswer = { 7, 7, 7, 8, 7, 7, 8, 7 };

    uIntVector subdomainNodesAnswer = { 0,  1,  5,  6, 25,
                                       26, 30, 60, 65, 66,
                                       85, 86, 90, 91, 50,
                                       51, 55, 75, 76, 80,
                                       81, 62, 63, 67, 68,
                                       87, 88, 92, 93, 10,
                                       11, 15, 16, 35, 40,
                                       41,  2,  3,  7,  8,
                                       27, 28, 33, 12, 13,
                                       17, 18, 37, 38, 42,
                                       43, 52, 53, 58, 77,
                                       78, 82, 83 };

    uIntVector subdomainNodeCountResult;
    uIntVector subdomainNodesResult;

    floatType minDistance = 1.0;

    error = dc.getSurfaceSubdomains( minDistance, subdomainNodeCountResult, subdomainNodesResult );

    if ( error ){
        error->print( );
        results << "test_dualContouring_getSurfaceSubdomains & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( subdomainNodeCountResult, subdomainNodeCountAnswer ) ){

        results << "test_dualContouring_getSurfaceSubdomains (test 1) & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( subdomainNodesResult, subdomainNodesAnswer ) ){

        results << "test_dualContouring_getSurfaceSubdomains (test 2) & False\n";
        return 1;

    }

    results << "test_dualContouring_getSurfaceSubdomains & True\n";
    return 0;
    
}

int test_dualContouring_getBoundaryInformation( std::ofstream &results ){
    /*!
     * Check that the boundary points and cells can be extracted
     *
     * :param std::ofstream &results: The output file
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

    YAML::Node yf = YAML::LoadFile( "dualContouring.yaml" );
    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_getBoundaryInformation & False\n";
        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){
        error->print( );
        results << "test_dualContouring_getBoundaryInformation & False\n";
        return 1;
    }

    error = dc.evaluate( );

    if ( error ){
        error->print( );
        results << "test_dualContouring_getBoundaryInformation & False\n";
        return 1;
    }

    const uIntVector boundaryCellsAnswer
        =
        {
             0,  1,  2,  3,  5,  6,  7,  8, 10, 11, 12, 13, 15, 16, 17, 18, 25,
            26, 27, 28, 30, 33, 35, 37, 38, 40, 41, 42, 43, 50, 51, 52, 53, 55,
            58, 60, 62, 63, 65, 66, 67, 68, 75, 76, 77, 78, 80, 81, 82, 83, 85,
            86, 87, 88, 90, 91, 92, 93
        };

    const floatVector boundaryPointsAnswer
        =
        {
            -0.960002, -0.986542, -0.996311, -0.96004 , -0.98658 , -0.49883 ,
            -0.96004 , -0.98658 ,  0.495371, -0.960002, -0.986542,  0.992852,
            -0.96004 , -0.490606, -0.996349, -0.960117, -0.490606, -0.49883 ,
            -0.960117, -0.490606,  0.495371, -0.96004 , -0.490606,  0.99289 ,
            -0.96004 ,  0.500508, -0.996349, -0.960117,  0.500508, -0.49883 ,
            -0.960117,  0.500508,  0.495371, -0.96004 ,  0.500508,  0.99289 ,
            -0.960002,  0.996444, -0.996311, -0.96004 ,  0.996482, -0.49883 ,
            -0.96004 ,  0.996482,  0.495371, -0.960002,  0.996444,  0.992852,
            -0.47429 , -0.98658 , -0.996349, -0.47429 , -0.986658, -0.49883 ,
            -0.47429 , -0.986658,  0.495371, -0.47429 , -0.98658 ,  0.99289 ,
            -0.47429 , -0.490606, -0.996428, -0.47429 , -0.490606,  0.992969,
            -0.47429 ,  0.500508, -0.996428, -0.359348,  0.616812,  0.61188 ,
            -0.474518,  0.50028 ,  0.992869, -0.47429 ,  0.996482, -0.996349,
            -0.47429 ,  0.996561, -0.49883 , -0.474517,  0.996461,  0.495143,
            -0.474541,  0.996424,  0.992832,  0.496394, -0.98658 , -0.996349,
             0.496394, -0.986658, -0.49883 ,  0.496394, -0.986658,  0.495371,
             0.496394, -0.98658 ,  0.99289 ,  0.496394, -0.490606, -0.996428,
             0.496394, -0.490606,  0.992969,  0.496394,  0.500508, -0.996428,
             0.381452,  0.616812,  0.61188 ,  0.496622,  0.50028 ,  0.992869,
             0.496394,  0.996482, -0.996349,  0.496394,  0.996561, -0.49883 ,
             0.496621,  0.996461,  0.495143,  0.496645,  0.996424,  0.992832,
             0.982106, -0.986542, -0.996311,  0.982144, -0.98658 , -0.49883 ,
             0.982144, -0.98658 ,  0.495371,  0.982106, -0.986542,  0.992852,
             0.982144, -0.490606, -0.996349,  0.982221, -0.490606, -0.49883 ,
             0.982221, -0.490606,  0.495371,  0.982144, -0.490606,  0.99289 ,
             0.982144,  0.500508, -0.996349,  0.982221,  0.500508, -0.49883 ,
             0.982221,  0.500508,  0.495371,  0.982144,  0.500508,  0.99289 ,
             0.982106,  0.996444, -0.996311,  0.982144,  0.996482, -0.49883 ,
             0.982144,  0.996482,  0.495371,  0.982106,  0.996444,  0.992852
        };

    const uIntVector *boundaryCellsResult = dc.getBoundaryIDs( );
    const floatVector *boundaryPointsResult = dc.getBoundaryPoints( );

    if ( !boundaryCellsResult ){
        results << "test_dualContouring_getBoundaryInformation (test 1) & False\n";
        return 1;
    }

    if ( !boundaryPointsResult ){
        results << "test_dualContouring_getBoundaryInformation (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( *boundaryCellsResult, boundaryCellsAnswer ) ){

        results << "test_dualContouring_getBoundaryInformation (test 3) & False\n";
        return 1;

    }
    
    if ( !vectorTools::fuzzyEquals( *boundaryPointsResult, boundaryPointsAnswer ) ){

        results << "test_dualContouring_getBoundaryInformation (test 4) & False\n";
        return 1;

    }

    results << "test_dualContouring_getBoundaryInformation & True\n";
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

    test_dualContouring_constructor( results );
    test_dualContouring_loadPoints( results );
    test_dualContouring_loadFunction( results );
    test_dualContouring_getFunctionValue( results );
    test_dualContouring_evaluate( results );
    test_dualContouringInternalPointResidual( results );
    test_dualContouring_performVolumeIntegration( results );
    test_dualContouring_performRelativePositionVolumeIntegration( results );
    test_dualContouring_performSurfaceIntegration( results );
    test_dualContouring_performPositionWeightedSurfaceIntegration( results );
    test_dualContouring_performSurfaceFluxIntegration( results );
    test_dualContouring_performRelativePositionSurfaceFluxIntegration( results );
    test_dualContouring_getSurfaceSubdomains( results );
    test_dualContouring_exportConfiguration( results );
    test_dualContouring_getBoundaryInformation( results );

    test_KDNode_constructor( results );
    test_KDNode_getIndex( results );
    test_KDNode_getMinimumValueDimension( results );
    test_KDNode_getMaximumValueDimension( results );
    test_KDNode_getPointsInRange( results );
    test_KDNode_getPointsWithinRadiusOfOrigin( results );

    //Close the results file
    results.close();
}
