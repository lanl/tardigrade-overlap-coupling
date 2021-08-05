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
    floatVector integratedVolumeAnswer = { 6.52002 , 13.04003 , 19.560051 };

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
    floatVector integratedVolumeAnswer = { 0.48573563, -0.32615985,  0.0715282 ,
                                           0.97147225, -0.6523191 ,  0.14305681,
                                           1.45720788, -0.97847915,  0.21458421 };

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
    floatVector integratedSurfaceAnswer = { 22.3709189 , 44.74184179, 67.11275769 };

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
    floatVector integratedSurfaceAnswer = { 0.903828 , -0.3204276, -0.3170567 };

    error = dc.performPositionWeightedSurfaceIntegration( functionValues, 1, integratedSurfaceResult );

    if ( error ){
        error->print( );
        results << "test_dualContouring_performPositionWeightedSurfaceIntegration & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( integratedSurfaceResult, integratedSurfaceAnswer ) ){
        vectorTools::print( integratedSurfaceResult );
        vectorTools::print( integratedSurfaceResult - integratedSurfaceAnswer );
//        assert( 1 == 0 );
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
    floatVector integratedSurfaceAnswer = { 6.81813769, 20.45440861, 34.09068637, 13.63627627, 27.27255183, 40.90882794 };
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

    uIntVector subdomainNodeCountAnswer = { 31, 50, 26, 31, 33, 37, 22, 42, 24, 32 };

    uIntVector subdomainNodesAnswer = { 244, 241, 240, 225, 224, 221, 220,  13,  12,  11,   3,   2,   0,
                                          1,   8,   9,  10,  44,  45, 112, 113, 114, 115, 120, 121, 140,
                                        141, 146, 147, 148, 149,
                                        315, 314, 307, 306, 305, 304, 291, 290, 287, 286, 275, 271, 270,
                                        267, 266, 211, 210, 207, 206, 205, 204,  87,  86,  79,  78,  77,
                                         76,  75,  74,  73,  72,  66,  65,  64,  63,  62,  61,  60,  57,
                                         56,  98,  99, 168, 169, 184, 185, 186, 187, 190, 191,
                                        327, 326, 325, 324, 319, 318, 299, 298, 295, 294, 279, 219, 218,
                                        215, 214, 199,  81,  80,  67, 106, 107, 110, 111, 171, 197, 198,
                                        311, 310, 303, 302, 301, 300, 283, 282, 281, 280, 261, 260, 203,
                                        202, 201, 200,  91,  90,  89,  88,  85,  84,  83,  82,  71,  70,
                                        161, 180, 181, 182, 183,
                                        263, 262, 247, 246, 243, 242, 227, 226, 223, 222,  55,  54,  51,
                                         50,  17,  16,  15,  14,   4,   5,   6,   7, 116, 117, 118, 119,
                                        142, 143, 144, 145, 160, 162, 163,
                                        272, 257, 256, 253, 252, 249, 248, 245, 237, 236, 233, 232, 229,
                                        228,  47,  22,  21,  20,  19,  18,  28,  29,  30,  34,  35,  36,
                                         37,  46, 122, 123, 126, 127, 130, 131, 132, 150, 151,
                                        312, 309, 308, 289, 288, 285, 284, 269, 268, 265, 264, 209, 208,
                                         92,  53,  52, 164, 165, 166, 167, 188, 189,
                                        274, 259, 258, 255, 254, 251, 250, 239, 238, 235, 234, 231, 230,
                                         59,  58,  49,  48,  23,  24,  25,  26,  27,  31,  32,  33,  39,
                                         40,  41,  42,  43, 124, 125, 128, 129, 136, 137, 138, 139, 157,
                                        158, 159, 170,
                                        278, 277, 276, 273,  69,  38, 133, 134, 135, 152, 153, 154, 155,
                                        156, 172, 173, 174, 175, 176, 177, 178, 179, 194, 196,
                                        323, 322, 321, 320, 317, 316, 313, 297, 296, 293, 292, 217, 216,
                                        213, 212,  96,  95,  94,  93,  68,  97, 100, 101, 102, 103, 104,
                                        105, 108, 109, 192, 193, 195 };

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

        std::cout << "subdomainNodeCountResult: "; vectorTools::print( subdomainNodeCountResult );
        results << "test_dualContouring_getSurfaceSubdomains (test 1) & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( subdomainNodesResult, subdomainNodesAnswer ) ){

        std::cout << "subdomainNodesResult: "; vectorTools::print( subdomainNodesResult );
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

    if ( !boundaryCellsResult ){
        results << "test_dualContouring_getBoundaryInformation (test 1) & False\n";
        return 1;
    }

    if ( !boundaryPointsResult ){
        results << "test_dualContouring_getBoundaryInformation (test 2) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( *boundaryCellsResult, boundaryCellsAnswer ) ){

        std::cerr << "boundaryCellsAnswer: "; vectorTools::print( *boundaryCellsResult );
        results << "test_dualContouring_getBoundaryInformation (test 3) & False\n";
        return 1;

    }
    
    if ( !vectorTools::fuzzyEquals( *boundaryPointsResult, boundaryPointsAnswer ) ){

        std::cerr << "boundaryPointsResult: "; vectorTools::print( *boundaryPointsResult );
        results << "test_dualContouring_getBoundaryInformation (test 4) & False\n";
        return 1;

    }

    results << "test_dualContouring_getBoundaryInformation & True\n";
    return 0;
}

int test_dualContouring_planes( std::ofstream &results ){
    /*!
     * Test of adding planes to the dual-contouring reconstruction
     * 
     * :param std::ofstream &results: The output file
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
    std::string filename = "dualContouring.yaml";

    YAML::Node yf = YAML::LoadFile( filename );

    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){

        error->print( );

        results << "test_dualContouring_planes & False\n";

        return 1;
    }

    error = dc.addBoundingPlanes( planePoints, planeNormals );

    if ( error ){

        error->print( );

        results << "test_dualContouring_planes & False\n";

        return 1;
    }

    error = dc.evaluate( );

    if ( error ){

        error->print( );

        results << "test_dualContouring_planes & False\n";

        return 1;
    }

    floatVector boundaryPointAnswer = { 3.14624245, -0.58973865,  1.96173039,
                                        3.14621245, -0.58976315,  1.9912409 ,
                                        3.1462123 , -0.55710646,  1.96170811,
                                        3.14621234, -0.54642226,  1.97394173,
                                        3.20609047, -0.58976303,  1.96170838,
                                        3.14758837, -0.58976332,  2.02434053,
                                        3.21264493, -0.58976307,  2.10260218,
                                        3.15088781, -0.51528436,  1.96170808,
                                        3.15637448, -0.52816972,  2.01724568,
                                        3.15293863, -0.53259985,  2.14554033,
                                        3.19954193, -0.51989308,  2.22231445,
                                        3.20389118, -0.43909901,  1.9617082 ,
                                        3.21344689, -0.41880961,  2.01710172,
                                        3.18321885, -0.42620054,  2.14074451,
                                        3.17979211, -0.42695534,  2.25749191,
                                        3.20070033, -0.40092358,  2.32778356,
                                        3.196445  , -0.33607496,  2.25098662,
                                        3.19838641, -0.31266659,  2.35286294,
                                        3.32564102, -0.58976348,  1.96170783,
                                        3.35990643, -0.58976352,  2.00945377,
                                        3.34251756, -0.58976311,  2.1340786 ,
                                        3.3661107 , -0.53328505,  1.96170798,
                                        3.31699003, -0.57061446,  2.15746817,
                                        3.32683025, -0.52503339,  2.22321196,
                                        3.30389987, -0.40144022,  1.96170817,
                                        3.31156868, -0.36387611,  1.97545033,
                                        3.2861402 , -0.37643852,  2.11020685,
                                        3.30208444, -0.41851237,  2.23439705,
                                        3.32692949, -0.40637556,  2.32937681,
                                        3.33938063, -0.30007008,  2.01188586,
                                        3.33409423, -0.28720245,  2.10159078,
                                        3.33521977, -0.24908168,  2.24286836,
                                        3.28876121, -0.27880342,  2.3828105 ,
                                        3.33292532, -0.2350507 ,  2.14939876,
                                        3.33650357, -0.2350505 ,  2.2099354 ,
                                        3.42772558, -0.58976316,  1.96170822,
                                        3.45933537, -0.58976317,  2.01006961,
                                        3.43424256, -0.58976321,  2.13356247,
                                        3.46508223, -0.52127788,  1.96170819,
                                        3.49810503, -0.53920553,  2.00947126,
                                        3.47666936, -0.58159057,  2.14942075,
                                        3.44757519, -0.51959154,  2.21833698,
                                        3.43921879, -0.42414819,  1.96170817,
                                        3.47242729, -0.3906117 ,  1.98096535,
                                        3.43707985, -0.45000911,  2.26188838,
                                        3.46046976, -0.39238664,  2.35406438,
                                        3.44405469, -0.30300132,  2.02083616,
                                        3.45079483, -0.26260223,  2.11663748,
                                        3.46232655, -0.26384303,  2.23693969,
                                        3.45210295, -0.30382994,  2.37392987,
                                        3.43051796, -0.23505056,  2.15112966,
                                        3.440795  , -0.23505066,  2.20031471,
                                        3.50504241, -0.53202963,  2.02564324,
                                        3.50504238, -0.54614217,  2.12083118,
                                        3.50504233, -0.49606366,  2.22013826,
                                        3.5050424 , -0.39995012,  2.01278625,
                                        3.50504262, -0.43290066,  2.10702498,
                                        3.50504243, -0.43105794,  2.23243772,
                                        3.50504238, -0.37597857,  2.35189499,
                                        3.50504237, -0.3247448 ,  2.04779941,
                                        3.50504239, -0.26327344,  2.13714537,
                                        3.50504262, -0.29382428,  2.21316026,
                                        3.50504249, -0.3064408 ,  2.35365778 };

    if ( !vectorTools::fuzzyEquals( *dc.getBoundaryPoints( ), boundaryPointAnswer ) ){

        results << "test_dualContouring_planes (test 1) & False\n";

        return 1;

    }

    results << "test_dualContouring_planes & True\n";

    return 0;

}

int test_dualContouring_localBoundary( std::ofstream &results ){
    /*!
     * Test of adding the local boundary to the dual-contouring reconstruction
     * 
     * :param std::ofstream &results: The output file
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
    std::string filename = "dualContouring.yaml";

    YAML::Node yf = YAML::LoadFile( filename );

    volumeReconstruction::dualContouring dc( yf );

    if ( dc.getError ( ) ){

        dc.getError( )->print( );

        results << "test_dualContouring_localBoundary & False\n";

        return 1;

    }

    errorOut error = dc.loadPoints( &points );

    if ( error ){

        error->print( );

        results << "test_dualContouring_localBoundary & False\n";

        return 1;
    }

    error = dc.reconstructInLocalDomain( element );

    if ( error ){

        error->print( );

        results << "test_dualContouring_localBoundary & False\n";

    }

    error = dc.evaluate( );

    if ( error ){

        error->print( );

        results << "test_dualContouring_localBoundary & False\n";

        return 1;
    }

    floatVector boundaryPointAnswer =
        {
            3.1999805344473833, -0.5000046698911095, 1.9999601287613666, 3.224058447116205, -0.5112864982599548, 1.9943332386722705,
            3.2721750268764205, -0.533829250469031, 1.9830865295043565, 3.3202915120185668, -0.556371964459089, 1.9718398449081322,
            3.3684081990555796, -0.5789147744636315, 1.9605931136217274, 3.4165259476550216, -0.6014580782047004, 1.9493461328051747,
            3.4547736898209616, -0.6193772234753362, 1.9404061435532194, 3.188774556277059, -0.5027274242276549, 2.0251068452141694,
            3.220165027829058, -0.5219990275601722, 2.0306699800514507, 3.256568657812959, -0.5409801146832548, 2.027603293036878,
            3.3004454002653434, -0.5601642759866245, 2.01346994617017, 3.338277228846838, -0.5781124920411455, 2.005259901447891,
            3.383101531855535, -0.597632339162113, 1.990599147933953, 3.4310169846493306, -0.619424260927222, 1.9775439823256946,
            3.4973288683373873, -0.6511848413443936, 1.9640033963125967, 3.1663876034351826, -0.5081690879107428, 2.075355078968203,
            3.198783053432488, -0.5198738260847292, 2.057970537544896, 3.223979682948424, -0.5381090791806893, 2.0702528801421196,
            3.28669315937439, -0.5645514068389113, 2.047221131776052, 3.3109471293143073, -0.5758513831305202, 2.0414414440875905,
            3.3637588457348806, -0.6005966053699918, 2.0291051634688464, 3.415968762922292, -0.6251366203077878, 2.0171263916501987,
            3.465698556403386, -0.6482641965550006, 2.0050194641159877, 3.144001511428971, -0.513609585090606, 2.1256008194647835,
            3.1677313663432543, -0.5185683661835107, 2.102649474498326, 3.2112069236078575, -0.5421857484442956, 2.1016686370696873,
            3.257934774869386, -0.5645737290737963, 2.0921477274444356, 3.3039289763682125, -0.5864823753571626, 2.0824150705417446,
            3.3511040440767137, -0.6091570958782887, 2.073007837309153, 3.3997033737129554, -0.6324588948006085, 2.06315411721783,
            3.44778231082752, -0.6557562267361867, 2.0540983674095954, 3.121618498799108, -0.519050675439658, 2.175840439373426,
            3.143707177236447, -0.5265433903204609, 2.162607351621262, 3.191343918836073, -0.5502124621198483, 2.15529087256536,
            3.239451932326336, -0.5730275638466821, 2.1448271196110515, 3.286468661494512, -0.594909179354458, 2.1334253238300867,
            3.3331711588488426, -0.6167972996334962, 2.12253137629022, 3.380584344492135, -0.6392709851817722, 2.1121851384045334,
            3.4263991767923874, -0.6619736655777675, 2.1049756978019083, 3.1034546649768022, -0.5234655694227278, 2.216609611487845,
            3.1210199072424913, -0.5332026290022952, 2.2167645437942594, 3.1695556838093313, -0.5566808387004802, 2.207508327364723,
            3.2181113037529303, -0.5796880671984035, 2.1968902770497967, 3.2659611338653125, -0.6019861067815561, 2.18536754538257,
            3.313124985370782, -0.6238483821511293, 2.1736819304628385, 3.3596034055452755, -0.6455362094593249, 2.1625710204413875,
            3.4131804288567724, -0.6704267998289648, 2.149453428343284, 3.122328226726893, -0.5473783644686429, 2.2547852487090596,
            3.1551525498839066, -0.564416198816963, 2.251802599770206, 3.2014033453962405, -0.5868287663696056, 2.243094262456511,
            3.2492389790671714, -0.6093080520225203, 2.2321058193902066, 3.298851137486882, -0.6318138804583695, 2.2184251592959967,
            3.344867898443173, -0.6542082700123262, 2.2100299967669788, 3.209800869574223, -0.481600290374025, 2.018545130824814,
            3.229456221427316, -0.50283236227285, 2.0030393515900187, 3.274865826111603, -0.5184937973047595, 1.9975196559676762,
            3.3249841336057013, -0.5358189842096766, 1.9913914864540467, 3.378216271863893, -0.5559128782789093, 1.9833466816526792,
            3.432031655910731, -0.5789787828327944, 1.9727162450827669, 3.483389520272626, -0.6004034925504942, 1.9631047869583211,
            3.1942870455155736, -0.4945580101392331, 2.030059863492514, 3.461032962374329, -0.6014974369175489, 1.9962061307656982,
            3.510327490680327, -0.6338992491470723, 1.9798456239601037, 3.177162277310889, -0.4986569826886403, 2.0686577844013962,
            3.484485662764601, -0.6279457967843041, 2.0322081497135835, 3.1593538363804217, -0.4900077581508224, 2.141552863011768,
            3.4638686212584124, -0.6355617146212533, 2.0772926151781363, 3.124738156434781, -0.5050716931494368, 2.2023779309411373,
            3.440630096282687, -0.6410565525180161, 2.126860463283024, 3.1053379261910163, -0.5125025851812145, 2.239033548843242,
            3.122701581631949, -0.5094340370065727, 2.253981025084338, 3.3692871787575984, -0.627558596266531, 2.187402882836444,
            3.421860052876216, -0.6494103796260092, 2.171416448796931, 3.1292209599876983, -0.5315558751009521, 2.2710858348353193,
            3.1619963146120367, -0.5443542848960513, 2.2761589289656947, 3.2041935738009775, -0.5654370113403228, 2.2709505359292987,
            3.2507744831254697, -0.587893301141203, 2.2608504613022498, 3.302526394973658, -0.610099189800276, 2.245879435616255,
            3.355663333319648, -0.6336653221708365, 2.229062201813369, 3.229437723630468, -0.44480997701415403, 2.0556951180825065,
            3.2495955300203483, -0.46910923056476245, 2.0375011273682064, 3.293301613979976, -0.4876215529230281, 2.029067865705054,
            3.3378198448074796, -0.506788065617997, 2.0201963803937124, 3.3859499181745956, -0.5253363223412059, 2.012577513702781,
            3.4376268753251087, -0.5457274911500571, 2.003965163258196, 3.489374701590433, -0.5658246452590928, 1.9956332188821457,
            3.2175036687765632, -0.45233643835990306, 2.0707463909213137, 3.4724157469656802, -0.5698839353905109, 2.026115454180059,
            3.5214071803154265, -0.6029063360719651, 2.0076829412769523, 3.20018682906614, -0.4455687707933215, 2.1374617968907534,
            3.498374730861728, -0.593225617854682, 2.065463848323885, 3.171494794020946, -0.45508768973023755, 2.195405426562143,
            3.479842889054133, -0.5994521542088115, 2.113308914436313, 3.141237177749647, -0.47328680911720766, 2.235807544931972,
            3.4587002691546997, -0.6045340392206372, 2.1639480888135356, 3.1269399349146707, -0.4736450061590119, 2.2758054899578903,
            3.140639496461726, -0.4720977730689735, 2.2930668450016087, 3.3885052452376767, -0.5913196052353498, 2.225947687348648,
            3.4421408665850364, -0.6129411097689624, 2.208111177730339, 3.146347217377218, -0.49287463606340015, 2.307694240617821,
            3.1804377942843187, -0.5074676587578189, 2.313494483220562, 3.222402260506955, -0.5283685147739432, 2.3090312981374335,
            3.269150365346718, -0.5509589330090329, 2.299270955746566, 3.3210165129390883, -0.5733759163633305, 2.284363312001445,
            3.3755049704436306, -0.5964878179287393, 2.265561798037042, 3.249078328923171, -0.4080131843552514, 2.09285092017802,
            3.271369977581782, -0.42981058504320724, 2.0773360807498453, 3.316603479713747, -0.4506426624125557, 2.067089924165023,
            3.360239792749038, -0.47257790183924664, 2.05553678879988, 3.4047148439887525, -0.49228427597467556, 2.0461670726698085,
            3.4527620479992986, -0.5115702460729552, 2.037862815194967, 3.5035008846865487, -0.5296005682557035, 2.03121347218572,
            3.2385598278892416, -0.4158268841490077, 2.103123398684068, 3.4907860426928905, -0.5338037318021146, 2.0617215037367593,
            3.5419376711173705, -0.5639889988210387, 2.0485986725145087, 3.224076652320523, -0.40766869514421294, 2.165258759699081,
            3.5201137363085975, -0.555904956709813, 2.1029501110625772, 3.190479198102804, -0.4209484044427264, 2.22769671754722,
            3.5018996694187603, -0.5624191530408145, 2.150493103869564, 3.1703287253185657, -0.42633053342607663, 2.27169637643435,
            3.4808922560219315, -0.5678684514158878, 2.199831260688791, 3.1497188811515806, -0.4318847353887234, 2.316573880868247,
            3.157191938352076, -0.43350576137347124, 2.327826558372401, 3.4084262833552668, -0.5558527422171164, 2.2629659329078713,
            3.460180075323213, -0.5745788880774929, 2.247704783257525, 3.163007867375596, -0.454599540734302, 2.346722144703588,
            3.1981383813415576, -0.4692161677597629, 2.351615864381882, 3.240647735493597, -0.49078337935402516, 2.3467293447978856,
            3.2874428741512403, -0.513605399459848, 2.337119460613007, 3.3386280669173787, -0.5361334936257728, 2.32301651539484,
            3.3946519104110706, -0.5582475705525858, 2.3034399978709663, 3.26871489971316, -0.37122429162453596, 2.1299983308608765,
            3.289656253562755, -0.39248178384914567, 2.1147154193283777, 3.3345683851601677, -0.41082706955698306, 2.1066646885521005,
            3.3816724704318752, -0.4296748862242608, 2.0985775363858097, 3.4321845181054655, -0.44885240822953176, 2.0908436222674984,
            3.4849169448103434, -0.4703090515215302, 2.081466357888411, 3.536783348823119, -0.49093124514570347, 2.072680601614232,
            3.2535973517794203, -0.38581184835710286, 2.136244517416495, 3.518976766574648, -0.49257752822644046, 2.1038831759479386,
            3.568069436279721, -0.5224302705006497, 2.090850005016926, 3.2380229082632312, -0.3844416367271202, 2.184281606470161,
            3.5434146282228505, -0.5177963641303427, 2.1424189750143223, 3.222451156786926, -0.3784595113236462, 2.2440108692598995,
            3.523979064672467, -0.5251341903001951, 2.18846039081733, 3.193676762334032, -0.3852143534639763, 2.309202541282208,
            3.501086540930684, -0.530743973062475, 2.2375220986017936, 3.167412316583473, -0.39758845743642995, 2.3529573063668523,
            3.1799045275521607, -0.39785632670796606, 2.36312516988477, 3.4284938965314087, -0.5181546620181, 2.299172275039209,
            3.4800748613728256, -0.5380961351172134, 2.285288536400198, 3.1864011197438717, -0.42017292476344054, 2.3840883908300583,
            3.2196302696814, -0.43298943988272387, 2.388073220761789, 3.262034716138119, -0.45443140895291784, 2.382398912669505,
            3.308780699422575, -0.47709292765256345, 2.3724556398099685, 3.360323898914128, -0.49943719563573974, 2.3582987294153077,
            3.414408784036575, -0.5224135568765432, 2.3411849123681114, 3.288346048569361, -0.3344453161880576, 2.1671360473185155,
            3.306972791823767, -0.3575575800308696, 2.149726797229448, 3.3520711830076304, -0.37269832671040415, 2.1446200687589028,
            3.402051097284782, -0.39035858855683714, 2.1381612924371307, 3.4551139899277015, -0.41084451787198767, 2.1297282815881107,
            3.5076501027923306, -0.4339746533190328, 2.1187946238133466, 3.558348791006091, -0.45588562026323, 2.1086156557699316,
            3.274325575503638, -0.3437957601174167, 2.183529473557204, 3.5381841982187745, -0.4574183392474961, 2.1414886425618627,
            3.585270294415278, -0.4893911425555259, 2.1235245327449777, 3.256527692143026, -0.3461479928512843, 2.2284846638858213,
            3.560649461850295, -0.48182798760922374, 2.1798860188567177, 3.2350667287560744, -0.34114065220889456, 2.302591137363466,
            3.540584790617554, -0.48873375660360463, 2.2253844936940212, 3.1991982296482813, -0.36194470296967457, 2.352438387845219,
            3.518029240913693, -0.49394651029830405, 2.27557054479868, 3.1833984554426866, -0.3662592484006274, 2.386698403241355,
            3.2017475518207057, -0.36274940807635603, 2.4042245637442083, 3.4474976389349625, -0.4802956506999781, 2.336427046788089,
            3.5004416628029364, -0.5025847005770998, 2.3190565069957754, 3.2080313369363185, -0.38446686850097905, 2.4194814752090905,
            3.2408238565613092, -0.3975949223912221, 2.4249307888010168, 3.2827949579111833, -0.4184785785840646, 2.4200237405397136,
            3.329495037135696, -0.4409775521626331, 2.4100124064931694, 3.3812821151152344, -0.4631668249576366, 2.394895367639947,
            3.434337263192249, -0.48669386728400343, 2.3772985902171744, 3.3079873969946636, -0.29764647382925397, 2.2042947626904663,
            3.3291016664415474, -0.32189108280057743, 2.1863332175404664, 3.373487920008253, -0.34147767984267463, 2.17705519184208,
            3.4164607406559506, -0.36220886805249075, 2.166467773110755, 3.4611824211518583, -0.3810229399215667, 2.157955096485408,
            3.5099783037680137, -0.40013648914716216, 2.149950665050034, 3.5610290232003847, -0.4188120025303098, 2.1427754853878396,
            3.2982710006154328, -0.30217590883819445, 2.2206042790177793, 3.5481586051361016, -0.42432169550879145, 2.1712876963176386,
            3.598372918062628, -0.45622998859155384, 2.1559107043086234, 3.2787141425039077, -0.29890528837529956, 2.2848569905968366,
            3.5765067864102753, -0.44643197511647364, 2.21332472775557, 3.247325467271933, -0.3103114627589613, 2.3457282527905807,
            3.558845067004131, -0.45231985216379594, 2.2619537336214757, 3.22372546229405, -0.3223614663225544, 2.3826816895769705,
            3.538208655988063, -0.4574541054732487, 2.3121864624555597, 3.207575545055374, -0.3234566092293041, 2.4261108663023174,
            3.2172157775337458, -0.3236782406443122, 2.4391866893182224, 3.4672984112610785, -0.4445960829162207, 2.374165796020323,
            3.520560478806819, -0.4653303168114573, 2.357135045527994, 3.2230917905389176, -0.34461624971023314, 2.455734535003663,
            3.2579077057681194, -0.3594161219678812, 2.4614944886864927, 3.30034411813162, -0.38073540442033865, 2.457180203194386,
            3.3470869985108798, -0.40344337442107553, 2.4475549673029042, 3.3984284064182444, -0.425815179043959, 2.4329367566103572,
            3.454049126564075, -0.44881498677333775, 2.4141339574348213, 3.3276272853491498, -0.2608514072445358, 2.241448343491522,
            3.349547561770517, -0.2820154902602186, 2.2264374054553184, 3.3946402093149377, -0.3021713938976763, 2.216777965772587,
            3.439480959282046, -0.32309434005655285, 2.206374158437049, 3.4863997065774845, -0.34246573361565963, 2.1977763273207085,
            3.536751425157866, -0.3625298345055823, 2.1892070829266372, 3.588435293685315, -0.381442681989335, 2.1819377229357,
            3.3135782759809844, -0.2754314737812399, 2.2446563178498007, 3.57312690723952, -0.38422877096401004, 2.2135491419757187,
            3.6243423089471447, -0.41435870137510394, 2.200442207857703, 3.2993759921510453, -0.26965526998133776, 2.2999450979224054,
            3.600917820246567, -0.40782976583287023, 2.2524685609459967, 3.28236889961363, -0.2672525000341078, 2.3547005739508347,
            3.5823906743762124, -0.4148620173309349, 2.2994067328392043, 3.2600574409301553, -0.2685153527666058, 2.4153333812012088,
            3.560129742324814, -0.42038519335837105, 2.348848819255589, 3.2354545632650717, -0.2769706827485398, 2.46427566568614,
            3.243478352496153, -0.288341821550212, 2.4696710363983807, 3.489010062347211, -0.4074234981543611, 2.410933424518577,
            3.540642446119464, -0.4284770495189538, 2.396282157808702, 3.2483336560435276, -0.31156237762347755, 2.4956797117676057,
            3.2799155423716946, -0.32248394663927804, 2.4981275498807216, 3.3233574594975086, -0.3445008129842154, 2.492478443852093,
            3.370315725071772, -0.3672553834354104, 2.482757262187475, 3.4217654381683698, -0.38936963761197063, 2.46935797708183,
            3.472254673747925, -0.41239191404688913, 2.45316508442892, 3.3472235265451284, -0.22413936124392597, 2.278516470754928,
            3.363069138314193, -0.24194276936348458, 2.265394897515043, 3.4156552308222543, -0.25182849378854855, 2.266491042330745,
            3.4710672253877908, -0.27294316567191945, 2.2579371242446866, 3.5243746061458676, -0.2975206811848064, 2.245837518247907,
            3.5755158130883937, -0.3228015744626086, 2.232684834138764, 3.6241966308063382, -0.3452021202877199, 2.221675199089501,
            3.332163228556783, -0.2337195526233363, 2.2973018764446156, 3.6021211966408977, -0.34473255319327595, 2.259097838963728,
            3.6488440571683514, -0.37809799553324175, 2.2415548008127835, 3.309153000758668, -0.24254195856224608, 2.3407561093194156,
            3.6244106528620406, -0.3714989358162172, 2.2959106886362766, 3.2888020604958452, -0.24664677809446248, 2.3885698036409284,
            3.6022955168981436, -0.37881096447775076, 2.3407058491664148, 3.2686673458981925, -0.2508712093975209, 2.4354613395430063,
            3.544726013757031, -0.365783876878227, 2.3906402845404022, 3.573669501202944, -0.3837430126622774, 2.392096822272117,
            3.2535375068964454, -0.24685319691270194, 2.48901113451211, 3.2787455258083393, -0.255486851746006, 2.4902541608992483,
            3.3188775647766717, -0.2777362918783908, 2.491989284240005, 3.365615884678872, -0.30115674752624444, 2.4778928251733516,
            3.412759608469287, -0.32378289773693225, 2.465473965193198, 3.459062632395011, -0.3445771305451452, 2.455524393062436,
            3.51070467782317, -0.3612144624939231, 2.45085350785132, 3.5560724358265396, -0.3938042011123265, 2.434545014736116,
            3.2729511761680774, -0.2776232748765604, 2.5332536186598458, 3.3046492962391056, -0.28757996093471877, 2.537226928899442,
            3.3487588796176175, -0.3090678749110002, 2.526791073205332, 3.396404171189075, -0.3315180545926889, 2.514924906734217,
            3.445711487594915, -0.3524232418305871, 2.499414422834858, 3.4856313715330587, -0.3765769811722406, 2.494558096785215,
            3.3669162037635614, -0.18724089673701483, 2.315780315362832, 3.3937726230442697, -0.20298941796013625, 2.3066261248612983,
            3.4513338853823874, -0.22531075764878286, 2.2973884618224774, 3.50409178941018, -0.25076718034232814, 2.2843860009580914,
            3.553960271613325, -0.27449631431061083, 2.272397993991964, 3.601551813213038, -0.29668363533237396, 2.261373345368413,
            3.647155551730606, -0.31820145056296006, 2.2505757343494177, 3.3545790748469084, -0.19134316618362676, 2.3406722906155295,
            3.615436750625613, -0.3233353977660824, 2.296004184284071, 3.6472941151407863, -0.3532446920414153, 2.267269405264136,
            3.331633663571382, -0.19767605324088605, 2.3902567542255064, 3.60508492447543, -0.3197755560316544, 2.338107118592952,
            3.625989315839464, -0.3454135056517253, 2.321819652848014, 3.3083007158989095, -0.2048676380726055, 2.4387715916889277,
            3.5823537856658554, -0.3202308900752621, 2.3829468897527666, 3.6021815165444755, -0.34623498297922933, 2.3756457485999714,
            3.2831692915294135, -0.21882998100490408, 2.475255372881354, 3.322481872737792, -0.21658201796379667, 2.483413722256272,
            3.3698900403290106, -0.24039639388683698, 2.474562920195657, 3.417796161377143, -0.2634608042377145, 2.4635038435616625,
            3.4660043760210435, -0.2848926903821183, 2.451782383584823, 3.5121127722012893, -0.3034544543295663, 2.440817257226078,
            3.5501670991930783, -0.3219971927020306, 2.425659018615125, 3.586941218032619, -0.35316909492047216, 2.4228405026260456,
            3.270705036212689, -0.20860455473794973, 2.5368583721592524, 3.2947829328496265, -0.21799826643324655, 2.533928389801373,
            3.34050486894033, -0.2416206115556645, 2.5217111543136617, 3.385988035151569, -0.264656475658833, 2.512738007197333,
            3.433253221692031, -0.2869543754866083, 2.503186319954281, 3.482687078352682, -0.30762591859633276, 2.4905407327761866,
            3.520527792708424, -0.33082678086970185, 2.4810049407484054, 3.3865857337326366, -0.1503926498488736, 2.35298453395121,
            3.407214880492972, -0.15766475883330233, 2.3503344459583198, 3.4572200056354983, -0.18826853762082874, 2.3321334611393683,
            3.507984978880209, -0.21410173585457587, 2.318407556281999, 3.565461828966286, -0.23785034039030287, 2.307858597593824,
            3.6000928950094266, -0.26148822017311757, 2.293036098143868, 3.6403641043117783, -0.28662844438472107, 2.277929955490276,
            3.3748946724187716, -0.1536584675626272, 2.378149383044184, 3.5478654249597352, -0.23354887028297608, 2.342338455575862,
            3.602766071816389, -0.26106203236252323, 2.3217659498454237, 3.6394305080578397, -0.28980972911298375, 2.2953042990256587,
            3.3521188554167423, -0.15953835552990306, 2.428397822753755, 3.3851757435603513, -0.17331662164053374, 2.432398839393265,
            3.4271483815036317, -0.19428173678650545, 2.4168175918071504, 3.5300920065497943, -0.23997974777961165, 2.3974075681385214,
            3.582068959083048, -0.25997356933661, 2.3674942980131197, 3.610061287295044, -0.2876924065301778, 2.3520587767073367,
            3.329154938242791, -0.16642257566000696, 2.4766365464948437, 3.3640064607802507, -0.17652879160449067, 2.4660203256520683,
            3.412800296302094, -0.19870558659071513, 2.4514463289989634, 3.458878357353694, -0.22162753509295688, 2.4468952298797246,
            3.5092043847135983, -0.2427015529518564, 2.4281202225219443, 3.5545506274893803, -0.26340637548823637, 2.4192039208985348,
            3.579332654058262, -0.2866481139753204, 2.4070886496340194, 3.3152372953650113, -0.16764949880360758, 2.513344396829152,
            3.3366177034227436, -0.17680240239924608, 2.5156830297517305, 3.384760115757157, -0.20123126384340567, 2.504311845594198,
            3.4314499500258666, -0.22400382381226985, 2.4938385555938005, 3.4779357910820448, -0.245346633957932, 2.4838698503774683,
            3.522804502355404, -0.26485329996289203, 2.472995021432533, 3.549078232005791, -0.290459266849506, 2.4600084282134076,
            3.401405627462871, -0.12262726617345693, 2.3810210042940176, 3.42516067865342, -0.12005099222410204, 2.3879073376172806,
            3.4725422429735655, -0.14502895217237077, 2.374309806095438, 3.518434293989019, -0.16846651088796985, 2.361825101083268,
            3.559179292196084, -0.19460605486386914, 2.345902720586179, 3.395186665476223, -0.11683292928076927, 2.41351439108969,
            3.426640600379893, -0.13005801189719904, 2.4054375194804263, 3.4748598593383413, -0.15256851977773062, 2.3942660923647354,
            3.5233169732087193, -0.17503503621575203, 2.3821014386399, 3.573214834351267, -0.1978179776733212, 2.3677350012899687,
            3.610757493351051, -0.21795316465761666, 2.3601631452377356, 3.3780856472293057, -0.12148216389275161, 2.4506483368792478,
            3.403931181029558, -0.13325130795230028, 2.4508573834929535, 3.4532715883926794, -0.15771269840947524, 2.44195087224615,
            3.500544909061378, -0.17941460077809143, 2.4293126095954944, 3.5416602089404763, -0.2011610729323164, 2.412143473868147,
            3.580270536255888, -0.22218236894836452, 2.3859922031617606, 3.4202326636400135, -0.15823227177542468, 2.4906569096651547,
            3.4656099318308033, -0.1788637563934917, 2.474676245640877, 3.507140123169552, -0.20051153348444717, 2.45385750093558
        };

    if ( !vectorTools::fuzzyEquals( *dc.getBoundaryPoints( ), boundaryPointAnswer ) ){

        results << "test_dualContouring_localBoundary (test 1) & False\n";

        return 1;

    }

    floatType volumeAnswer = 1;

    floatVector volumeResult;

    error = dc.performVolumeIntegration( floatVector( points.size( ) / 3, 1 ), 1, volumeResult );

    if ( error ){

        error->print( );
        
        results << "test_dualContouring_localBoundary & False\n";

        return 1;

    }

    if ( !vectorTools::fuzzyEquals( volumeResult[ 0 ], volumeAnswer ) ){

        std::cerr << "volumeResult: " << volumeResult[ 0 ] << "\n";
        results << "test_dualContouring_localBoundary (test 2) & False\n";

        return 1;

    }

    results << "test_dualContouring_localBoundary & True\n";

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
//    test_dualContouring_planes( results );
//    test_dualContouring_localBoundary( results );
//
//    test_KDNode_constructor( results );
//    test_KDNode_getIndex( results );
//    test_KDNode_getMinimumValueDimension( results );
//    test_KDNode_getMaximumValueDimension( results );
//    test_KDNode_getPointsInRange( results );
//    test_KDNode_getPointsWithinRadiusOfOrigin( results );

    //Close the results file
    results.close();
}
