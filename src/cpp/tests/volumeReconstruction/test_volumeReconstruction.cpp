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
            3.4460578219915035, -0.25207439818173594, 2.429969888797944, 3.4442857522353925, -0.24739548300586595, 2.421552218777412,
            3.4132731621133923, -0.2368740396124042, 2.438236036066035, 3.4106198851623613, -0.23180884839880722, 2.42916434233434,
            3.381820253924296, -0.22220138915463725, 2.446958557794997, 3.379655711468153, -0.21732230632156962, 2.437936669787781,
            3.5234392766417133, -0.2800894679667475, 2.3773117492896123, 3.523049909521056, -0.2744064454727303, 2.3711537260646787,
            3.5002127108930985, -0.2655387012302839, 2.3890476548949975, 3.4959996290381747, -0.26033265241745307, 2.381425324712378,
            3.4687390250786527, -0.25148260128606226, 2.3972563110706933, 3.461957529032526, -0.24518177552610967, 2.3897474358441904,
            3.4348509792363826, -0.23602946603276515, 2.404971649819214, 3.429182786186924, -0.22998574262163068, 2.3974299323173995,
            3.402787149655732, -0.22088604215487018, 2.41260073682531, 3.3971181567995705, -0.21484629905420313, 2.4050009773383008,
            3.3705361310833193, -0.20531423819432054, 2.420013193564647, 3.3644258252502204, -0.19908235000044308, 2.412405063895113,
            3.3422744838008294, -0.19203093275186875, 2.4255128038743416, 3.339771338966803, -0.18832836029274866, 2.4176217157620723,
            3.543145110407192, -0.2784892595059219, 2.3445709659588543, 3.5401463966511915, -0.2719940413239564, 2.3388537962391682,
            3.515838113910183, -0.2627051784778775, 2.3578142261623736, 3.5093709699043303, -0.25666317205500155, 2.3507025991189088,
            3.482750988985905, -0.24798982023466057, 2.36680482244241, 3.4762922717825107, -0.24174873597593757, 2.359025667139863,
            3.449830618779397, -0.23289696785632027, 2.374426563039828, 3.4434785310551357, -0.22654415580904452, 2.3666091676419794,
            3.417161705606367, -0.2175931874737282, 2.382012777432542, 3.410802498037343, -0.21122483583837043, 2.3741717303153274,
            3.3846486433971834, -0.2022460592250593, 2.3895518005515526, 3.3782981072295932, -0.19584945373555404, 2.3816886066199427,
            3.355572144247802, -0.18898317205200402, 2.3966507818512333, 3.3530854764296842, -0.18467998819774084, 2.387777950719815,
            3.5598935928652433, -0.27602543091119125, 2.310758544474747, 3.5561220564948606, -0.26981603704390955, 2.3044505366292913,
            3.5299077162754875, -0.25974264103360006, 2.325434161350248, 3.52363927775203, -0.2536056422635464, 2.317778084725505,
            3.4966438141773764, -0.24428508902613708, 2.3352673671916175, 3.490333033911066, -0.2379995922760323, 2.3269834481402873,
            3.4636294037919906, -0.2291104271785415, 2.3429243526906833, 3.457389413250867, -0.22273323682094923, 2.3346417833637587,
            3.4308793329450666, -0.21379861381607404, 2.350485551197841, 3.424622939023467, -0.20741270064172115, 2.34219494511198,
            3.3981955679482296, -0.19844231717449762, 2.3580179302978572, 3.391912885777765, -0.19203990256487347, 2.3496964670559555,
            3.3695287572402903, -0.18514385527182994, 2.364494312246624, 3.367612014803457, -0.18089523721799247, 2.354857941796234,
            3.568883702092747, -0.2752883266999627, 2.2830725339676743, 3.557254423149879, -0.26926438961564486, 2.283421737214525,
            3.538224339026788, -0.25765510118984036, 2.2971436632513984, 3.527600366576612, -0.2521967125971917, 2.2944285680702725,
            3.5067045250338476, -0.2418617569210707, 2.306779070493048, 3.4969276632058586, -0.23708177539397116, 2.302416142336889,
            3.4743161815406465, -0.22684242039159871, 2.314337950657204, 3.464885320731068, -0.22214012650124484, 2.309828432383661,
            3.441702158518525, -0.2116450297564794, 2.3218416590043542, 3.4322153957014763, -0.20683623481241495, 2.317431553451831,
            3.4088534940348034, -0.19620149108825766, 2.329463376317555, 3.398816664041547, -0.19069407578133213, 2.3256132190526313,
            3.379959313345873, -0.18228386920440517, 2.3362099393406432, 3.3752509920864644, -0.18000421587916363, 2.3305244826320024,
            3.4809736597789294, -0.30783617296609206, 2.428204689541145, 3.479849977307464, -0.3019359488143939, 2.4178828966864025,
            3.4530713998360976, -0.2948378166436657, 2.4352071683703227, 3.447756700835963, -0.2893839427763853, 2.4265508290540336,
            3.419467279245068, -0.2807233151172316, 2.443959654251351, 3.414667305951551, -0.2747945859566949, 2.4346326388378943,
            3.3867490006000303, -0.2654195662201943, 2.451327319788669, 3.382058692363805, -0.2594917095609858, 2.4418515620713004,
            3.3546424537194826, -0.24964773860650338, 2.458567287649827, 3.350272206175773, -0.24346654260512224, 2.448222960991879,
            3.3246265516929077, -0.2351725750343764, 2.4654884242612347, 3.3211353757971573, -0.22944490815035007, 2.455085607888549,
            3.5228539509987558, -0.3213930584216071, 2.3880195488583857, 3.5227403066513077, -0.3142896530666419, 2.3772821462117455,
            3.5020219601698206, -0.30471540474228526, 2.391957860476706, 3.4992490700500003, -0.30060942303563776, 2.3859491784665288,
            3.4717623881876047, -0.2923760932403331, 2.403217875359936, 3.468780342244251, -0.28852278977153506, 2.3992120548973923,
            3.3406419392868147, -0.23165912729811886, 2.4306968625074035, 3.3363365208067446, -0.22631159367922202, 2.4251888489904116,
            3.3148061757239735, -0.2201928216015774, 2.437472915318331, 3.311457804383092, -0.21713664964073767, 2.4299712280227426,
            3.540022787417503, -0.31989820737890495, 2.353360488156446, 3.5414494010975055, -0.3136033413959247, 2.3448401832658727,
            3.556895402887274, -0.3185641519709995, 2.321243415364253, 3.556756663042087, -0.3118638096584542, 2.3138585192171353,
            3.5712140100889505, -0.3174422550520869, 2.2872015223424484, 3.569620702426979, -0.31251007837748446, 2.2788024662022286,
            3.5017359472892378, -0.34793124595639363, 2.4010222013143654, 3.4995705442060747, -0.34145725419187845, 2.391604984760474,
            3.2826171659564722, -0.24593728873701656, 2.4430457511936394, 3.2800672942046307, -0.24572272257277322, 2.4312406129142694,
            3.4438236799170947, -0.3616417905995129, 2.414050210335913, 3.4362382045393436, -0.35484538966166573, 2.4068282071264226,
            3.315583023551597, -0.3027692878646192, 2.4449413426786104, 3.305732773076087, -0.2945878502404651, 2.440175395015266,
            3.2891352142519406, -0.28934448011742564, 2.4474215179365553, 3.279352492789357, -0.2805921523908987, 2.4409533589643857,
            3.422009196722895, -0.3825655871488697, 2.391483465996079, 3.4334630117704026, -0.38420680121845396, 2.381862382855288,
            3.3892474849042675, -0.3672854668843008, 2.3978754031919682, 3.40133006718554, -0.3706081525391578, 2.39102523573722,
            3.3557587556791053, -0.35189081955220197, 2.4054602210649043, 3.3678659191360274, -0.3553431584941807, 2.398998291643143,
            3.3230356298153625, -0.3364751110374941, 2.4130208308105963, 3.3349397577516333, -0.3398256837549092, 2.4065452626128865,
            3.2896446705972546, -0.32009162935633406, 2.4205488791000835, 3.301532495301922, -0.3236078046713011, 2.4134285898176113,
            3.2611234864680303, -0.3063420434781576, 2.4232006247767752, 3.2690693713979844, -0.307067070061755, 2.417618766734504,
            3.4655450754495503, -0.3914874201293113, 2.3525824563610143, 3.483242103025416, -0.39662885064124376, 2.337715812014453,
            3.5100328297605436, -0.38123285320037953, 2.254591049511847, 3.524697994217847, -0.3839666389140224, 2.2417358843169755,
            3.5201284601537437, -0.3758812698186128, 2.2231694384752267, 3.537608584704591, -0.38133398829868964, 2.208387049354636,
            3.501803666044749, -0.3616729914588916, 2.1992725497131325, 3.512759395663187, -0.36403804506310733, 2.192007662732461,
            3.3871199011404536, -0.38668990365388345, 2.3452382175399573, 3.403979584595572, -0.3907808203135675, 2.334522082103188,
            3.3539651239887562, -0.37153358079168136, 2.3533043155228603, 3.3691487098660797, -0.37570565030126135, 2.3432025421731373,
            3.3215150006500473, -0.356280080417746, 2.3606272635213092, 3.336403038879274, -0.3603684430910786, 2.3508818025814837,
            3.2898908813382275, -0.34106955267477795, 2.3669802646595923, 3.3046951801921436, -0.34518170329577524, 2.357561870760838,
            3.260962184190435, -0.32815296076580835, 2.3709516408840665, 3.2755993592499633, -0.3319697306108886, 2.361926534519157,
            3.2376301336320203, -0.31781020845137353, 2.374193531453732, 3.249182782853652, -0.3207012686560949, 2.364903593263868,
            3.4725143969546126, -0.41252045420269245, 2.2917046355579704, 3.4871889738239994, -0.41678733308389, 2.277915359970577,
            3.2528085082090055, -0.31562341271534455, 2.3476921905541452, 3.262028210111819, -0.314468419359611, 2.3387803149790543,
            3.4816852351956404, -0.40915651324546726, 2.2642747014387727, 3.4970544632884364, -0.41479664342422673, 2.2522167333464367,
            3.4722761402115263, -0.44136239435979857, 2.234333891082725, 3.471277431961096, -0.43304310423273557, 2.226376907291386,
            3.485617751296958, -0.4387416979582558, 2.1983006198097796, 3.485063997351971, -0.4305418797542998, 2.1912776729434498,
            3.4945197058130404, -0.4157054320197084, 2.143437091022918, 3.485332758265044, -0.40937802480282287, 2.1426576432100854,
            3.434886730000531, -0.4685748662409477, 2.2791200405598064, 3.432151935071354, -0.4595297196222442, 2.27101047904367,
            3.409815773584642, -0.45123089168328356, 2.285440910969401, 3.4035954565221544, -0.445761150164981, 2.2779838143114626,
            3.3766794525562966, -0.43977710349367377, 2.293212023191766, 3.373697792745301, -0.43522652422302416, 2.282972555879068,
            3.273323108332781, -0.3926339577379926, 2.3150202596933145, 3.268581809155982, -0.3861276806749439, 2.3052490235464904,
            3.247736035782049, -0.377584031971303, 2.318829277457761, 3.2437173113189717, -0.37136616864326316, 2.308990201610069,
            3.2217402087206377, -0.3660853219397479, 2.3232939530465075, 3.218437363464801, -0.36648122330708477, 2.3103205690607345,
            3.4508716136316204, -0.46348900779165575, 2.2486757778414925, 3.453446872533167, -0.4578149651550237, 2.2378616697789178,
            3.5013257183956066, -0.4573112979264043, 2.148383848083763, 3.499760713484268, -0.4517109217535209, 2.140312809003907,
            3.381122539830029, -0.4820608920792062, 2.2939512822072383, 3.3736273697973638, -0.4749228744921042, 2.286691245179533,
            3.350099993361345, -0.46873346788013176, 2.298764931398119, 3.340630744142351, -0.46179476301871325, 2.2931673048228998,
            3.317383291242145, -0.4541964209686053, 2.3069763556052827, 3.308049554792836, -0.4468398896110066, 2.3008180572810817,
            3.2849027582547583, -0.43888979906726566, 2.3147872382009256, 3.27566714196233, -0.43144410713378684, 2.3086429509578617,
            3.2544322676589243, -0.42418909647904524, 2.3211054474106536, 3.24487632249527, -0.4160191129825102, 2.316358721827514,
            3.4591495771555008, -0.496635783229352, 2.1000677478495513, 3.474340482164226, -0.5008761645347963, 2.0881025929673434,
            3.423222838284699, -0.3234670618004931, 2.4497598989502123, 3.417994147351245, -0.31734215885045297, 2.4417331313996273,
            3.345938597875292, -0.5263789607875691, 2.2411205540876282, 3.357110746323354, -0.5279773788124019, 2.2311849365144196,
            3.395181197953765, -0.3099905227167957, 2.4534088631596718, 3.386165577280111, -0.3028358189497952, 2.447773095662457,
            3.3134554099031877, -0.511958149680374, 2.249649729009933, 3.325502730582714, -0.514104951702753, 2.2387089569753384,
            3.3653990615521017, -0.295759018510726, 2.459415685140817, 3.3564211472813468, -0.28848697178576616, 2.453567935786576,
            3.2800609549768867, -0.4967185481363395, 2.25774357235689, 3.2921111374247407, -0.4997369017906158, 2.2487064616305292,
            3.247178889817539, -0.4812858486358734, 2.265261430576954, 3.25904173398594, -0.48425764521442877, 2.2564153804222853,
            3.2136328236792697, -0.4652176514479397, 2.2720001442097666, 3.224469536478417, -0.46773382230069577, 2.2638865290649353,
            3.3899077298406763, -0.5356735361111893, 2.205372458370357, 3.4071756170203047, -0.5409586907314746, 2.191739866297669,
            3.352748618967473, -0.5203180589048486, 2.2171205688826614, 3.369628130680915, -0.5246513505432882, 2.2055230958994323,
            3.172230334973759, -0.4355066559798206, 2.2527651433987, 3.1809178139946357, -0.4363489656076887, 2.244474426179635,
            3.404149166158367, -0.5331507592360782, 2.1752864943309924, 3.420513746007357, -0.5375120417834575, 2.160072061752978,
            3.4193002816038445, -0.5299364728885131, 2.1410838238322025, 3.4347398937265723, -0.533486949665393, 2.1267984235689488,
            3.432080994012231, -0.5251503259084574, 2.1090473346347376, 3.4481814954551004, -0.5301318506889872, 2.0956587727218796,
            3.4248029208141753, -0.5093884823935739, 2.0505014546037024, 3.4347177106067703, -0.5083769124475992, 2.045351713564541,
            3.3148576823154214, -0.5326146252578108, 2.1959881738251656, 3.3282222098909253, -0.5358294276709145, 2.1896147050107633,
            3.282595747751869, -0.5177604983825858, 2.203740406121658, 3.2948951239975743, -0.5209648814424604, 2.1960766931223534,
            3.2504322258189706, -0.5026631731764668, 2.211209937726092, 3.2624708789440464, -0.5057546708227487, 2.2035176286143097,
            3.218111625484293, -0.48711509454865815, 2.218093755304392, 3.230558306018604, -0.49029400529387335, 2.2110385600924847,
            3.1854637569784474, -0.4709820507735454, 2.2240726266922586, 3.1987331357474424, -0.47402089607122383, 2.217999178161487,
            3.4234867865787884, -0.5405834037357947, 2.016621023574526, 3.4096694856435033, -0.530639222376793, 2.019276176629554,
            3.341700837099875, -0.5733247705726869, 2.1647869496800154, 3.3384956138424395, -0.5665072992529219, 2.1566869373192663,
            3.3127736183728365, -0.5598732204380582, 2.171629463545534, 3.305016158418942, -0.5536198058941094, 2.1637640745320605,
            3.2790159988946113, -0.546056473064957, 2.177916700227871, 3.2707986120084414, -0.5390369472550373, 2.170471510594458,
            3.245950691677984, -0.5305115871981809, 2.184434027487091, 3.237741043612092, -0.5233593570312363, 2.176949306565352,
            3.2136428629133915, -0.5143602436019137, 2.1903530594538867, 3.206445688035602, -0.5070243211069291, 2.1820407784858036,
            3.183829756309355, -0.49960225501534655, 2.1967909105473744, 3.1789566644209537, -0.4927802257337285, 2.187425982499656,
            3.3832671117509427, -0.586044561216901, 2.1245270477348246, 3.3830410477688746, -0.5783381417334218, 2.116944895447199,
            3.4004817580799007, -0.5838238714462687, 2.0916871304800098, 3.402048963491725, -0.5768578140763297, 2.0830853745320774,
            3.416962321800882, -0.582383168928217, 2.0563817490762704, 3.4171489820879217, -0.5751013378053303, 2.0479735795700735,
            3.430896619548436, -0.5815296528756513, 2.020551234023952, 3.429652347692068, -0.576363371842481, 2.01230290445607,
            3.282243987651276, -0.5884266155257303, 2.1815949015352185, 3.2769485697880647, -0.5822358328606956, 2.1733312826037565,
            3.253936469532912, -0.5748794209425736, 2.1855498931751782, 3.244923716023574, -0.5676803721658644, 2.1795426952803725,
            3.2238141337292645, -0.5605308865297478, 2.191639213234175, 3.2146799148397402, -0.5532011218596176, 2.1858576371828784,
            3.4228755576656678, -0.40179983697710414, 2.3357433134645755, 3.4402223511237358, -0.4054296647273062, 2.3227602253028183,
            3.361132017124218, -0.6128030420124913, 2.1338593588796697, 3.3587635274429566, -0.6065031537542409, 2.1241709264599162,
            3.14067380808878, -0.5100065995145919, 2.174894123230866, 3.1384599485300724, -0.509791193751557, 2.1631685054620235,
            3.2668917828680786, -0.5966108701659454, 2.168757881018956, 3.2791692976116433, -0.5989115568306567, 2.156135077317609,
            3.2377453502785376, -0.582410340986619, 2.1740294693683953, 3.2514576136764433, -0.5855679323146294, 2.161593302308314,
            3.2070894308853433, -0.5677954899060274, 2.180481512311928, 3.2182107487906575, -0.5700604249206639, 2.1695586417994424,
            3.3451990062583965, -0.6215060469905215, 2.117131960338011, 3.358870163785911, -0.6250301749107013, 2.105795600677697,
            3.311782072435634, -0.6056848064944585, 2.1244756456537193, 3.3277617675179667, -0.6093608271515745, 2.1099725957496247,
            3.354077666406385, -0.40151270043112597, 2.326998076830839, 3.3475970550864838, -0.3952065941569787, 2.3207177581813414,
            3.279531078828967, -0.5915114759321457, 2.1346599445563768, 3.2947527696202727, -0.5946132123417892, 2.1197147643817624,
            3.324794797692647, -0.3870955543005827, 2.3317473111405955, 3.314710294452759, -0.3795298148942971, 2.3261313879579637,
            3.2482167471643826, -0.5772112262186715, 2.143026587745686, 3.264012697573158, -0.5808540831714766, 2.1287160336950954,
            3.2150471080406335, -0.5617607924342377, 2.15103291772042, 3.2313408430158495, -0.565754002404235, 2.1369370307267443,
            3.182373185726546, -0.5459945711201155, 2.1573747198346713, 3.198133005609872, -0.5500537603440564, 2.144296938978158,
            3.1513596998648126, -0.5300291199273952, 2.1605671316416886, 3.1664895427013575, -0.5340561209618847, 2.1483796818649132,
            3.1302446522985643, -0.5205102341658289, 2.1665584848245576, 3.1399089390573063, -0.5214406658302461, 2.1541340950731898,
            3.356463998079458, -0.6166015310852182, 2.0857254204823796, 3.3721870501326467, -0.6209320721235978, 2.0734717119666186,
            3.3252716990863775, -0.6012167974877466, 2.0908373851303015, 3.3403027721205065, -0.6049333735012102, 2.07792656952337,
            3.294522227868943, -0.5870776083877899, 2.0987793012881664, 3.309061462252795, -0.590005172876631, 2.084404971940799,
            3.2629823984910136, -0.5731081177476242, 2.1084318699540767, 3.2782702379957236, -0.5762299065542935, 2.093440317721761,
            3.2303058793706283, -0.5580477534822933, 2.1167723658428437, 3.2461514844369472, -0.5617590118728775, 2.102577759120767,
            3.1974345187777526, -0.5421914946109799, 2.123167276847813, 3.213313657563009, -0.5461922777792194, 2.109738590706174,
            3.2415462086142792, -0.33555543582643593, 2.32661540090459, 3.2402148376507616, -0.33203109363602756, 2.3177463080931044,
            3.1671789883704187, -0.5261145194413323, 2.1248639013130997, 3.1823641445627846, -0.5301907038664865, 2.112729277152513,
            3.1437557059254604, -0.5169491500934421, 2.1354495406569027, 3.1558959124051893, -0.5175960026700388, 2.1183667966735866,
            3.37176424464947, -0.6134898402729201, 2.053099586346293, 3.387703661859021, -0.617748902142948, 2.040306864693994,
            3.33825964767082, -0.5969434726375595, 2.0585307178189876, 3.3539980681471553, -0.6011693983708901, 2.045957443157672,
            3.3067051256477344, -0.581595545419274, 2.0643109123357766, 3.320872238032134, -0.5848172691667392, 2.0513474738208965,
            3.2763822435798566, -0.5677095997701922, 2.0723039766022446, 3.2902778968484148, -0.5701476741463313, 2.057548900563626,
            3.244566518267604, -0.5531173823534848, 2.0806265624688853, 3.2595094379060323, -0.55594618144413, 2.065344332182678,
            3.213623489584179, -0.5381707153713838, 2.086588161057685, 3.229030745222566, -0.5413845615550652, 2.071670739542474,
            3.184414035720923, -0.5233272359615512, 2.0901410327182557, 3.199794241956794, -0.527306464298209, 2.0774286027650684,
            3.1598791307587435, -0.5133483447489303, 2.100158968777589, 3.173212097884498, -0.5150033936152887, 2.0840678289828456,
            3.3849215820040435, -0.6092678054430452, 2.0206740068152014, 3.39951545676543, -0.6127305483197183, 2.007726901954555,
            3.34936866123408, -0.5922558165879398, 2.0279799536015246, 3.364352757882358, -0.5960951617900672, 2.0154892389113535,
            3.3163120953663636, -0.5769377445474766, 2.036184188535772, 3.32907009809299, -0.5801503662353465, 2.0253899429777165,
            3.2883037925072003, -0.5639403905910715, 2.04308295525824, 3.3001708786395625, -0.5675087372447514, 2.034681681745266,
            3.2584843729824637, -0.5483795869635668, 2.0455588929190887, 3.273383074310539, -0.553935721734369, 2.0380525920332686,
            3.231456427765436, -0.5357634118100445, 2.05200774871404, 3.245967135884607, -0.5399296381279008, 2.041178132376861,
            3.2032450282371, -0.523073481293373, 2.0600916043531594, 3.219499958405662, -0.5283165826215269, 2.0495881718173425,
            3.175996408673293, -0.5107292234735058, 2.0676525831903207, 3.190649132790156, -0.5148200828149323, 2.0563887885138223,
            3.394749423713076, -0.6040527921404887, 1.9906286380038438, 3.4112661997986127, -0.6089933030119868, 1.9788622623345258,
            3.360460774222765, -0.5884602906181162, 1.9999765433140249, 3.3762635652910604, -0.5929551780918644, 1.9880630472671612,
            3.3277584318023647, -0.5744178955706812, 2.0112338952209208, 3.3441650112831742, -0.5788884747513234, 1.9983112102474005,
            3.301407838611527, -0.5625305013961851, 2.018687097384318, 3.3162373890468224, -0.5662270380342451, 2.0060335564180285,
            3.2736622847876475, -0.5481065834887378, 2.0211453304000453, 3.289819962955179, -0.5531330164968739, 2.0101811092640687,
            3.245901274011497, -0.5341858621196371, 2.0250495977617455, 3.2586977215894457, -0.536806672709672, 2.012523096670749,
            3.218346023027837, -0.5226007557059204, 2.035233466015322, 3.2302366305170374, -0.5243038088861738, 2.0215245082969022,
            3.192244337419455, -0.5095363840590683, 2.0389729870210545, 3.2036862974899294, -0.5120739556219025, 2.0283210970583703,
            3.3760601079835935, -0.5876750568526692, 1.973459053543937, 3.3923024225869423, -0.5924764164927186, 1.9617269815017706,
            3.3459276654072854, -0.5744330901868813, 1.982975274465851, 3.362019971039181, -0.5787247848144786, 1.9700366138092518,
            3.3178598709530593, -0.5609630145317124, 1.9886308730021713, 3.3311435678328842, -0.5646984875037067, 1.9784952760235148,
            3.288022341810739, -0.5457920728329211, 1.992236713443808, 3.301232699373, -0.5499255342296995, 1.9833400181482008,
            3.255419898107998, -0.5302395369016873, 1.999071020464143, 3.268586840163435, -0.5341985040168004, 1.989748858373179,
            3.2235242388653904, -0.5157290445386316, 2.0077490556582096, 3.235911058877955, -0.5190439619641098, 1.9978220689164454,
            3.2014736155060146, -0.5049595986325964, 2.0116638366868114, 3.2094209493900636, -0.5069116343690714, 2.004800753005069,
            3.4238445786718565, -0.2556685127077324, 2.436726866605242, 3.4327761607272302, -0.26970205490588084, 2.4269742048780056,
            3.3922460688504237, -0.24071762918141607, 2.4448363823290866, 3.40104541565679, -0.2548331870999961, 2.4349095428477616,
            3.361116560032053, -0.22574390682870718, 2.4526279440922822, 3.3688890980105684, -0.23919495564155094, 2.441434074894677,
            3.507194495341187, -0.28373671787478877, 2.387956717646589, 3.514946404559822, -0.2967104227751343, 2.3776089365500437,
            3.4819226975548503, -0.2712674205182045, 2.3995400995738714, 3.4885151148853493, -0.2837824146673109, 2.3922911391005295,
            3.3511960115686596, -0.2095489653504463, 2.4277743258319497, 3.356280046082364, -0.22288859076653425, 2.4197808912586805,
            3.325092893015557, -0.19862965189808607, 2.4310272039800185, 3.327904114951032, -0.20991579630250193, 2.4256805962680925,
            3.5673546903943416, -0.28357603562962325, 2.2667318972965704, 3.5598455630085146, -0.28661216933585104, 2.2625387386686486,
            3.430883437078004, -0.29938238962125535, 2.443860536180357, 3.4361026711859712, -0.311776979706963, 2.4368333662198314,
            3.397657833734084, -0.28458766436123095, 2.4518799912502214, 3.403759292568995, -0.29779286499772395, 2.4431245338629375,
            3.364977344330567, -0.26911920260803546, 2.4594551476474047, 3.371086204173757, -0.28235857158820077, 2.4506808036954997,
            3.3335099874347645, -0.2528658250250538, 2.465827918585839, 3.3391545698119103, -0.26640943236717834, 2.457240020866048,
            3.304449788594744, -0.23826327414798726, 2.4722999631721576, 3.3089629767955104, -0.24977532813747158, 2.4622951308388137,
            3.510248489164054, -0.3251179544197516, 2.3933282878751942, 3.5151288433065875, -0.3392492060036003, 2.386149703374328,
            3.296384959355409, -0.22841648817001767, 2.4379873294668752, 3.296304295009307, -0.2373653735960362, 2.430821125577508,
            3.579500418272039, -0.32374854039755857, 2.2595784229691978, 3.5706831655783686, -0.3224524674209993, 2.2593582308368383,
            3.5681656293069546, -0.3062285614073203, 2.2463285548261966, 3.553110780425815, -0.30543236631029763, 2.2441687925969127,
            3.4069568154709593, -0.32923422248065964, 2.450737012240235, 3.4126632378218735, -0.34156778704110424, 2.4392933635785656,
            3.3779559240851875, -0.3149213212251228, 2.4561311577546157, 3.384502850671014, -0.32804011803915994, 2.4454325392845937,
            3.3487778284200935, -0.30089604690788574, 2.46210210642527, 3.35231907028009, -0.31290364585570307, 2.4531818115418536,
            3.4812061411944, -0.3513726564390034, 2.406979263194128, 3.487996403033037, -0.36570140019743835, 2.3961854078244147,
            3.4494577427193596, -0.33676135127578866, 2.415158362775885, 3.455367634130511, -0.34903755542851844, 2.4036397384313855,
            3.5720569713904102, -0.33966754097063806, 2.243794698642683, 3.5560060281401316, -0.3410702159811207, 2.23697724260753,
            3.5511196613984057, -0.3221653258286459, 2.2286013299940457, 3.5307681908845137, -0.32531264479946903, 2.2218484812508295,
            3.515352466468677, -0.3069222431771204, 2.2355875249239237, 3.4960153602952775, -0.31164187606267857, 2.2276018988455397,
            3.4283189764419975, -0.36778196942167396, 2.4139041052483114, 3.4351702559243513, -0.38130561223337567, 2.4036852571303045,
            3.396404124211116, -0.35307003638645645, 2.421194324290959, 3.4030753349613017, -0.36571925868125327, 2.4085524815688344,
            3.364274259507925, -0.3384875518596299, 2.42977374166782, 3.370469110697771, -0.351183905646517, 2.4169791753358534,
            3.332251513911179, -0.3233787838562702, 2.4374398069985914, 3.337659776553742, -0.33584389304663687, 2.424730355594579,
            3.3012102849283074, -0.3081989921150656, 2.4449444971830876, 3.305739542468963, -0.32049382653086345, 2.4328505714970805,
            3.2749447079273954, -0.2963986939811597, 2.4454676312451733, 3.2771695281282835, -0.3070717683553507, 2.438224662828363,
            3.467439769722906, -0.37482968984730697, 2.3826590952659603, 3.4772648001708006, -0.389270734802699, 2.3678783840656537,
            3.2571616192911255, -0.27373237667900796, 2.4267735699115796, 3.2536948829849215, -0.28497659075857396, 2.4184037036606307,
            3.5556865162490388, -0.3642730461537741, 2.2154164270256844, 3.536597289936394, -0.3667985296379118, 2.2095238719567414,
            3.527390595063161, -0.3453883744255938, 2.2029817674375383, 3.5117629539480704, -0.35166139352473197, 2.194296564580474,
            3.421195652822518, -0.38512409431677397, 2.366316733553809, 3.430628384397585, -0.398189338300389, 2.350624755938983,
            3.4732548079086727, -0.39805442567451427, 2.322794636597927, 3.4832387116187933, -0.41088225890532576, 2.306640488986448,
            3.51211100711899, -0.3972421684744612, 2.256106147015593, 3.4955102468733656, -0.40044150446465354, 2.251299707293931,
            3.428935988294172, -0.4063369901238218, 2.3078701806143296, 3.4345232558399057, -0.4183603393793595, 2.292663010966526,
            3.3932161302213126, -0.39300612910334587, 2.321387105706639, 3.400825740460857, -0.4049324711643545, 2.3077182923795365,
            3.3577798524400415, -0.37746476399738826, 2.330061084907541, 3.364987589338626, -0.390674326140802, 2.3189617838778123,
            3.3241115872647335, -0.36150704120872756, 2.336920080352659, 3.3303048149464685, -0.3744312926066377, 2.32620113602844,
            3.2943696313984576, -0.34672037530611893, 2.343500364134331, 3.2988980002405595, -0.3589573505469828, 2.3326040346659784,
            3.2667236604595247, -0.3332717290544479, 2.3469934937165347, 3.2706896337530993, -0.34391022312981434, 2.3379857172304646,
            3.255451167504734, -0.3162054080499423, 2.3221043100787964, 3.2524047482506444, -0.3235230771017496, 2.3189902680345043,
            3.335143222871338, -0.40571028093405526, 2.329911132153043, 3.34130139781407, -0.41876571098528054, 2.320337171993528,
            3.3062858522829974, -0.3912719421223186, 2.3341101282823207, 3.3101536661429263, -0.4035805624820523, 2.325993630645271,
            3.414387053120963, -0.4292395172972581, 2.279233543529268, 3.4224610875363113, -0.44088491997413476, 2.2698561621200812,
            3.3845657171792425, -0.41728147134048393, 2.2903993076104165, 3.392539762292514, -0.42956516821748136, 2.2771217284079506,
            3.2818409833357527, -0.3695517127038944, 2.3144766713366605, 3.28517104830309, -0.38170207399045064, 2.3028537455223157,
            3.2546100389825394, -0.3543549417263487, 2.319645589938595, 3.2612330197229964, -0.36775945605962546, 2.3070328735131445,
            3.2304666813141827, -0.34868769596332033, 2.318266180640653, 3.2337914853208134, -0.3583575556769034, 2.3086295340497474,
            3.4632562911768052, -0.44353365762945035, 2.2407999167349417, 3.4678562464037612, -0.4573351970155195, 2.2304896877729727,
            3.495197473728742, -0.44382975146180526, 2.1714302290321226, 3.49050172527061, -0.44332313944714635, 2.1696336726823535,
            3.49696402627518, -0.4243784790950594, 2.1254677327539278, 3.4868571611261667, -0.425942648934872, 2.1221131221983582,
            3.4643054788238175, -0.4058544937942727, 2.1360267141888096, 3.4525620106463335, -0.40895458359715153, 2.130964824414836,
            3.4297582919541045, -0.38943740541627575, 2.1443119616160895, 3.4185798446681046, -0.39163150697298177, 2.1401804956275265,
            3.3961759885527423, -0.3731705718314718, 2.1526455718110036, 3.3849142504955863, -0.37495242812020696, 2.148872299290589,
            3.362579866229229, -0.3573485460134875, 2.1605728422082278, 3.351138325275819, -0.35887516970391364, 2.156996789529316,
            3.3284208526071843, -0.3414340569983749, 2.1684762579756467, 3.316256165640495, -0.34271778285023186, 2.164982229859564,
            3.3018904987821127, -0.32908820198982575, 2.174601493334758, 3.294636069571175, -0.3356107352780443, 2.167293477118018,
            3.417533331000619, -0.4717752985291355, 2.285103937969838, 3.423589528031568, -0.48628496210483696, 2.275788129613061,
            3.387211655062626, -0.457340887100235, 2.2907324391670922, 3.3921944175519934, -0.4688401473273853, 2.2822859046513932,
            3.3548634726426014, -0.4443908673144068, 2.2994768338026255, 3.3591753580881214, -0.45648490940318187, 2.2895401531200124,
            3.3218926135107507, -0.4296443093771258, 2.3068715110951215, 3.326985122246524, -0.44234418957178906, 2.296945687181312,
            3.288140657338497, -0.4136860006384588, 2.314147753342463, 3.293721429763781, -0.42669431748231734, 2.3045995939518282,
            3.2558696907883706, -0.3970202359142508, 2.3220696169097086, 3.260530424372866, -0.41026580679293534, 2.3134443847653556,
            3.2272288499851354, -0.3808594541962993, 2.3251303425771854, 3.2306212058154222, -0.39378938107635586, 2.3185994458823083,
            3.205964117603851, -0.3730609736985813, 2.326342558507279, 3.207434034407037, -0.3816172426277344, 2.3183398721902666,
            3.506790157203222, -0.46373770075556664, 2.1200204085099856, 3.494705453634007, -0.46304293027743826, 2.1156459807156307,
            3.4909077251146474, -0.4461378228807928, 2.1045602125933676, 3.4716706561045063, -0.4464332491872474, 2.100609007590427,
            3.4558766239385426, -0.42933816504335975, 2.1131000274214444, 3.43542023231787, -0.4309153669272478, 2.107752076476268,
            3.4209271717396623, -0.4127258870370844, 2.121485407481098, 3.4012265044447516, -0.41582203037624144, 2.1149035927567397,
            3.387617439316037, -0.3965095066003763, 2.1298254083827537, 3.3684297448149088, -0.4006094240364444, 2.122430810240589,
            3.354285015041047, -0.38057896924643275, 2.1379016484695135, 3.3355511514910683, -0.3853564110457896, 2.1299790380075367,
            3.3202397747591514, -0.36493633735509956, 2.145580150951257, 3.301764644675884, -0.37108468285456225, 2.13646288896477,
            3.290526606346293, -0.35538866163757765, 2.148556951656338, 3.2794365274481523, -0.36322764141538944, 2.139318888914156,
            3.3642148044634332, -0.4885464926074675, 2.2941036635991523, 3.3698571540820477, -0.5016090834203156, 2.2835293444510474,
            3.33306621743765, -0.47402200704243236, 2.3010379020934857, 3.3396858585089446, -0.4866970141834867, 2.289116597266812,
            3.3002010682439455, -0.45911151024546465, 2.309511161512979, 3.3064147452215744, -0.4719983914584358, 2.297975081425962,
            3.2678431435162767, -0.4438079634858804, 2.316966927796278, 3.2738077170573114, -0.456673143307118, 2.30553255278413,
            3.238416635396718, -0.4294531018388064, 2.321622224086195, 3.242591894962287, -0.4416262769129625, 2.3113394172287554,
            3.404953298465634, -0.49626882753883583, 2.262699070308782, 3.413248527322914, -0.5108936668177163, 2.248146808117387,
            3.2091219883810127, -0.40392155161773435, 2.3065685553010304, 3.2100975957370435, -0.41649827500342934, 2.301119471891779,
            3.1927692163708525, -0.3955483773504064, 2.303536907766609, 3.191064814515006, -0.4049227594606738, 2.2982964875722507,
            3.493404364623794, -0.48309899241742316, 2.0950879908482114, 3.474983826387279, -0.48599467757957415, 2.089310811925007,
            3.5391213233135392, -0.26684213473065194, 2.276513483405477, 3.5300592943995404, -0.2707988747501802, 2.2711875015398575,
            3.4674908097023445, -0.46501588796438836, 2.0829437774958484, 3.448959801425009, -0.4704390488368513, 2.074473997949755,
            3.5087614758337584, -0.25137347109114805, 2.2847396580484953, 3.4981712284274593, -0.2551072516682533, 2.2793234399700757,
            3.4326832562065, -0.4534336119443251, 2.0867912770991563, 3.416423115091404, -0.45890450935225763, 2.0787129408520637,
            3.4769065549888736, -0.2361132620821915, 2.292490405984926, 3.4658259687262136, -0.23932929570465766, 2.2874501979054345,
            3.3994594239078637, -0.43948240767499913, 2.093091925504816, 3.3851077869945967, -0.4436771984540959, 2.0865371327891125,
            3.444317971650193, -0.22054736518106666, 2.3003781187516457, 3.4328664370186073, -0.22334303863395186, 2.295648395181606,
            3.367687173030356, -0.42411673859183696, 2.100954201898871, 3.353341196865979, -0.4276230188357593, 2.0950253666266048,
            3.410727031075522, -0.20412449295980045, 2.308851672616869, 3.397948984164764, -0.20625518764090778, 2.304471513983963,
            3.3353499051005517, -0.4088931460799481, 2.108579364506346, 3.320836986674018, -0.4122503949198317, 2.102753825842406,
            3.3813302958701255, -0.19054511713817535, 2.3155476037110847, 3.372148689567745, -0.19534810811270664, 2.30943051589227,
            3.3025631029812663, -0.3948620053850226, 2.1150362373571094, 3.288020644483384, -0.3985722699575862, 2.1088846608188443,
            3.276132107279861, -0.3826339272043512, 2.1210732496843225, 3.267380762634213, -0.38695082497909156, 2.115479845426127,
            3.3528420490388546, -0.5128203596305347, 2.2665678802622296, 3.359094414588002, -0.5261509099263617, 2.256401336859339,
            3.3201250707384977, -0.4978408538343414, 2.2742022673815923, 3.3276389116505345, -0.5109016084684752, 2.2622278695251135,
            3.286800702099811, -0.48271158615631005, 2.2825127468104682, 3.2938913865451807, -0.4958790413911679, 2.270661237592772,
            3.2543281183003208, -0.46742773468588694, 2.289924857615899, 3.261107426365107, -0.48057929632203594, 2.278441131729254,
            3.2237912191594593, -0.45268183819453967, 2.2948389562091385, 3.229562868956014, -0.46561811614600135, 2.2848587716084303,
            3.4594143164129747, -0.313132630881397, 2.4385017912645406, 3.46614164452117, -0.3252235403521873, 2.42876943908763,
            3.391630632695616, -0.5202185211977893, 2.2342717652637667, 3.401382796022582, -0.5341669067297787, 2.219367005275096,
            3.190440282277628, -0.42720035068013473, 2.2825745732129046, 3.1944252748168345, -0.4400179324335823, 2.2726052993988923,
            3.177140477711415, -0.4200222280861496, 2.2819203940618147, 3.1755403381199994, -0.4297230959206614, 2.27164756311172,
            3.463369604087721, -0.5128892170344904, 2.0988344002215436, 3.44740474857736, -0.5161160259235675, 2.094441518954946,
            3.535728685213375, -0.28977271425751816, 2.2550529945737456, 3.51832184177013, -0.2881152721530351, 2.253224575643974,
            3.4472301406028145, -0.4900066950314626, 2.0563839147660308, 3.4341172846543104, -0.496090331144309, 2.0483520690179335,
            3.5020622035499973, -0.2728752132217684, 2.263942864531667, 3.4835561923549796, -0.2712378902477917, 2.2618857417851626,
            3.4195868043614093, -0.4819935823575772, 2.0583638159545643, 3.4054081505130083, -0.48609660196658183, 2.05192544445971,
            3.4685846309697816, -0.2566859623157647, 2.2722261091440785, 3.449425882066942, -0.25551028256726394, 2.269625042229527,
            3.3878928243297346, -0.467202319740412, 2.0657197766415036, 3.3740509334444364, -0.4700921421632252, 2.0604469340204172,
            3.4352156851917885, -0.24094928940812327, 2.2801194074208437, 3.41560290147037, -0.24023439894984797, 2.277013217373366,
            3.3569199199449202, -0.45007330112065713, 2.075335454895758, 3.3426983755753152, -0.45313994578630984, 2.0698294485219977,
            3.4004756063864297, -0.22475884285601874, 2.2881620465012213, 3.379768107690084, -0.22554804513464863, 2.283481216093093,
            3.324606734069841, -0.4340469123412266, 2.0836938215889034, 3.310269654971163, -0.43713255049720684, 2.078148459347726,
            3.3680598791660388, -0.2118473808932283, 2.293673595757019, 3.3558912103979757, -0.21674939814332306, 2.286895265808821,
            3.2919574490774566, -0.4190462976482244, 2.091056887062963, 3.2772926079881386, -0.422614293169067, 2.0850110060238234,
            3.264604458748418, -0.40645716325309406, 2.0972453755215072, 3.25504978550865, -0.41114250486557047, 2.091164210730303,
            3.3146305385628927, -0.5159247809159491, 2.223087743179791, 3.3247002536755, -0.5296376745843181, 2.208359110824794,
            3.2820165105016548, -0.5016562714250529, 2.233058465595612, 3.2924885909576567, -0.5156973485288356, 2.2183805191794685,
            3.318791946047926, -0.2766043921206263, 2.446386150535252, 3.321351672714455, -0.28913549992959814, 2.4379757245182803,
            3.249160420106323, -0.4862476510741709, 2.2404304657498977, 3.260085014172074, -0.5006279236314997, 2.226110661649583,
            3.2892059380720586, -0.260018958062812, 2.449438744643839, 3.2933903249466785, -0.2733044147039981, 2.442117662060069,
            3.2142329117634976, -0.46956195767687653, 2.2471949501986836, 3.227060951220261, -0.4851039898728158, 2.232584937038313,
            3.2712205789855835, -0.2514739653138497, 2.4463013258092574, 3.271191716962157, -0.2590879131473134, 2.4434282686502082,
            3.181795442534511, -0.4533890501538416, 2.2503065173344607, 3.1920142072598163, -0.4685342170538737, 2.238777216995534,
            3.3609483177875092, -0.5280875641454091, 2.193464030862691, 3.3672544424536084, -0.5401992330086028, 2.1800851198392,
            3.1783648665733755, -0.4418152218778302, 2.2281520518165343, 3.1791094119784136, -0.45129707729135254, 2.2223606074122655,
            3.428919034370218, -0.5244077636672754, 2.0413829360123468, 3.4165261970717373, -0.528219020615686, 2.0344759111853095,
            3.480478724401714, -0.2911975820653206, 2.2431818268190864, 3.4630484269709996, -0.2966234190758562, 2.2349203341766852,
            3.4029199861585337, -0.5043112188338965, 2.0349181082524828, 3.3837314080527148, -0.506626537130588, 2.029142986575563,
            3.4469671460214255, -0.27548804751068906, 2.2510231822282103, 3.430734329790566, -0.28105255660017203, 2.242865101305134,
            3.3723776815218907, -0.4913017426855779, 2.0408774801617025, 3.3546167283023753, -0.4936136437939469, 2.0353787835077775,
            3.413580282114665, -0.26035990119950686, 2.2583607712589724, 3.397921551727442, -0.26576785822711196, 2.250454686079716,
            3.342068435855317, -0.47547068793673886, 2.049442230992731, 3.323816566099503, -0.4790662761262528, 2.0426845119510815,
            3.3792366723868006, -0.24643155825702415, 2.26442630178598, 3.3639700867303377, -0.25238000316024256, 2.2561047676126047,
            3.3101094267208238, -0.4598211721314458, 2.0575263672391717, 3.2919214617522625, -0.4640761780568262, 2.050182416051791,
            3.3523328532655583, -0.23639850728804251, 2.268381351971519, 3.344262637531414, -0.2426792910460419, 2.261135921653578,
            3.277416848913284, -0.4449555519431283, 2.064758627121206, 3.2596088393617726, -0.449886065570453, 2.056874338012819,
            3.2499299251968585, -0.4322301484084614, 2.071045552124232, 3.2398510349585856, -0.43680841923993596, 2.0649602580230213,
            3.348732855485503, -0.5508128279770135, 2.1661024366902737, 3.3558155626378343, -0.5620412824190912, 2.151509950199559,
            3.315454950628479, -0.5370554943315469, 2.1766232539290358, 3.323655446732867, -0.5490254302658139, 2.1634140347965083,
            3.281402976371244, -0.5217545336712498, 2.1826080744100897, 3.2885982754258762, -0.5346303441486915, 2.1695281040311754,
            3.2487609437053004, -0.5062134268194098, 2.1893126012565487, 3.2553936648855895, -0.5190089066094885, 2.1760667600516346,
            3.21774716299068, -0.4904000308639502, 2.1959165165263372, 3.2239111475511315, -0.5032279547449312, 2.1818828668021832,
            3.1888412648764812, -0.4741741438931525, 2.2017369665854805, 3.195408611584433, -0.4873743679427066, 2.1867095001622388,
            3.1678399258662395, -0.463821799811523, 2.2041276878079246, 3.170719767141124, -0.4736539582646609, 2.192658550256092,
            3.4958034435960252, -0.33630529678458965, 2.2051776523153817, 3.4814593683280886, -0.34071333102918655, 2.1984307856303293,
            3.4239710567649273, -0.5448390520095577, 2.002166864518847, 3.415109605070069, -0.5473252375999464, 1.9982139057990265,
            3.4642409875966256, -0.3214757473178333, 2.212593539518087, 3.4504001405053577, -0.3241936501537702, 2.207476922500979,
            3.3903896413832606, -0.5276229237649215, 2.0113621793935494, 3.3850923168908804, -0.5338996738034708, 2.0046514262890813,
            3.4329747636488714, -0.3049058310009056, 2.2216456362364996, 3.418546551090256, -0.30762231260158046, 2.216417856846728,
            3.360673130340819, -0.5152151882789521, 2.0169335368171706, 3.354376919203063, -0.5188007206115134, 2.012473994984686,
            3.400647386844727, -0.28917699522806584, 2.2297312367361055, 3.386015112136099, -0.29199355823351875, 2.224373559389465,
            3.3300686474088046, -0.5007947643167482, 2.024161522633644, 3.322719084668099, -0.5031233021259635, 2.020641119232426,
            3.36802269136235, -0.2745328661351545, 2.236775474559698, 3.353266748407568, -0.27780355062297374, 2.2309819770369925,
            3.2982790383902363, -0.4856498198611784, 2.0318201573624837, 3.290228222575182, -0.4873236759529916, 2.0287596683316798,
            3.341143928041413, -0.2619314714056117, 2.243065839771959, 3.331841633430487, -0.26637662985748994, 2.237250582585115,
            3.2660626363856964, -0.47068054777596974, 2.039237637912602, 3.2561206629255244, -0.47140423625678035, 2.0366774415318787,
            3.239587975663328, -0.4551880248826034, 2.0482288475308876, 3.231863058890718, -0.4581653473652545, 2.04404757869206,
            3.388171372047801, -0.37176894862355525, 2.377238845558789, 3.3965831426357282, -0.38464367633282237, 2.36215142651617,
            3.3201156992726553, -0.5786287536631682, 2.172035577821488, 3.326641235377845, -0.5904634721176877, 2.161634067826347,
            3.3553016740459034, -0.3565400019577302, 2.385336688561904, 3.3635245446757707, -0.3696431064336706, 2.3704274099322995,
            3.291564739779427, -0.5650363415617792, 2.177219260500876, 3.2957015974235913, -0.5769744034676061, 2.16952536241608,
            3.3225994622935175, -0.3410801343786274, 2.3927056997215126, 3.330854215130284, -0.3543421136409438, 2.377959055092291,
            3.2590076777018764, -0.5503722682243123, 2.183580811534214, 3.2637989586339913, -0.5630987221777249, 2.175017982793747,
            3.289293109905758, -0.32478419354149013, 2.399039623947675, 3.2984559489598864, -0.33890595053292216, 2.384792302049285,
            3.226158617062527, -0.5347829983253125, 2.1906223535558325, 3.230998216201485, -0.5476241535418914, 2.182323947969542,
            3.2577549534242705, -0.3089300776706915, 2.4017115464272485, 3.2660908672252624, -0.32313102583680336, 2.390252692922892,
            3.1946729205380633, -0.518303175192263, 2.1959358297546148, 3.198979026621988, -0.5315736997490282, 2.1887204029348224,
            3.2405848002243767, -0.3003497491311953, 2.4005527805528715, 3.2410736211696984, -0.3110609523942548, 2.391700432578287,
            3.1651880861177304, -0.5033239394924898, 2.2016480887965355, 3.168834075469796, -0.5146889371752238, 2.1924807270356936,
            3.369872657798564, -0.5901271856636852, 2.128459388596826, 3.374697055697819, -0.6041672658306121, 2.1187342507364355,
            3.156431584100999, -0.49036046103032227, 2.1707626443131693, 3.154235708058694, -0.5008339034395277, 2.1624694959097215,
            3.504058680567762, -0.37571578439649067, 2.1931541381032, 3.4932226255523204, -0.381414668104362, 2.1808068348082505,
            3.438797641014838, -0.588256175089866, 1.9930236173559703, 3.429798981691042, -0.5869207769763174, 1.9925169431563876,
            3.481126050919059, -0.36250211540143074, 2.1785922210010775, 3.465273268384585, -0.36664142955194345, 2.171800368832369,
            3.4268873727476055, -0.5698675076978057, 1.980010302436957, 3.411534104549696, -0.5696406132300961, 1.977276727607232,
            3.4517657387023464, -0.34757054801555104, 2.1865223164445124, 3.436340118546297, -0.3520474870133361, 2.179505824243219,
            3.3944590123797784, -0.5533750931693612, 1.9887695778648338, 3.376942750409766, -0.552322110663577, 1.9863716213228855,
            3.4205463077906932, -0.33138707075236423, 2.1952326509099134, 3.4049014015749157, -0.3368393555404475, 2.1872889810276686,
            3.3612798437831124, -0.5364322236292148, 1.9977939279146495, 3.3423891216215815, -0.535520377841299, 1.9950047255124268,
            3.3882693701658484, -0.31589116039542176, 2.2031165096178413, 3.3727448160497313, -0.32162849368119095, 2.1949371793962564,
            3.328136634300828, -0.5202174126319591, 2.0061643891812455, 3.308351476525093, -0.5198401913794432, 2.002718729678432,
            3.3557346150651983, -0.30100157600737976, 2.2104007316451733, 3.340623887463747, -0.30670785240286125, 2.2023288189740233,
            3.29488311107731, -0.5044199369644934, 2.0141349706988096, 3.2746054491552115, -0.5045941462819369, 2.010094558361649,
            3.3288537261951907, -0.2879549854703812, 2.2170948071789995, 3.319124548921615, -0.29277403931964924, 2.21085830816397,
            3.2600085871929108, -0.4882215737506922, 2.0221590596137973, 3.238662735215936, -0.4900241849506156, 2.016436262512701,
            3.2270498941759396, -0.47618416993899365, 2.0267733755919553, 3.2144669577409, -0.4818717084412728, 2.019202825668567,
            3.2681685819425668, -0.5903911958129093, 2.186994626266351, 3.2733730216247565, -0.5985303979872236, 2.182673115164732,
            3.239038385242612, -0.5759930844576809, 2.192334445020945, 3.247975079791085, -0.585314076663342, 2.1855155591806557,
            3.2114758398234913, -0.5625612646145627, 2.1970831918292055, 3.2176873378280493, -0.5708244137660669, 2.19177122572117,
            3.345256431524829, -0.6136456290962183, 2.139416178249405, 3.353325022485396, -0.6245027898719869, 2.131721015700808,
            3.310616770092136, -0.5973941813964465, 2.1507471194228303, 3.321150720218937, -0.6069728655570762, 2.139330172790703,
            3.180883880406483, -0.5382809842932931, 2.1793006419357237, 3.1894322721114547, -0.5481511935665024, 2.1724150069505512,
            3.15095125233411, -0.5217758382507015, 2.182412462876796, 3.15869495118723, -0.5320973608119105, 2.176447804596619,
            3.131185685975775, -0.5138255256560997, 2.1803585253367133, 3.1339517142567486, -0.5187771657015365, 2.1789561796732313,
            3.496851841352057, -0.4026018286495196, 2.1598200183462937, 3.4861979105418177, -0.4025199852910555, 2.157462856404949,
            3.43084499321327, -0.6009939955891487, 1.9796972694770516, 3.4125399910277316, -0.5985523572785209, 1.9776461766210347,
            3.4682026282167224, -0.38582645144208416, 2.1549495648649644, 3.4562118293161936, -0.38679283829579947, 2.1517767821057974,
            3.410638288817172, -0.5841064853008535, 1.9639764977631857, 3.3939543839682416, -0.583457899858659, 1.9613708658205373,
            3.4378568493280066, -0.37208475540339264, 2.1616110938719078, 3.4231467801434015, -0.37162923272028536, 2.1592081631561735,
            3.3773339377307536, -0.5679390452821882, 1.9722731391226362, 3.361641839406297, -0.5680072898784343, 1.9692068150691844,
            3.4059851431459145, -0.35704424911765453, 2.169159223645238, 3.390192798841256, -0.35639148667156, 2.1667280982886323,
            3.3423599299237794, -0.5514785094722925, 1.980516094867211, 3.327953963984082, -0.5524614017258312, 1.9768658951649423,
            3.373584503154138, -0.34192392637928487, 2.176678525449476, 3.3572333085992674, -0.3410991294319084, 2.174296541889486,
            3.3081979979061495, -0.536126048842883, 1.9879088663073448, 3.2947604301762303, -0.537325654319471, 1.9842473923119472,
            3.34128289389341, -0.32650997312833446, 2.184483283165117, 3.3239119788547367, -0.3253182770037158, 2.182239051930188,
            3.2746362468730745, -0.5212140061740874, 1.995016843345874, 3.261673014423913, -0.5222203172689063, 1.9916216141114986,
            3.313523795151497, -0.31033521012654236, 2.193848088724067, 3.3036670020455605, -0.3128214092606326, 2.1897045582447174,
            3.240296364729076, -0.5074030685770722, 2.000976532026178, 3.2277013512403996, -0.5081856803881997, 1.9978548247401609,
            3.2134847852922213, -0.4964943352636418, 2.005744005303296, 3.2075912034559337, -0.49853520667042966, 2.0027629468264396,
            3.454989404043855, -0.26610794038176383, 2.4202172270707067, 3.462187924911311, -0.27600779562423605, 2.4064610153707338,
            3.363354782977526, -0.22018953539173378, 2.434841127054724, 3.3529536419816117, -0.22867072256033733, 2.43800099411187,
            3.327853139300586, -0.19337900143800507, 2.4120040564595686, 3.3181467976254035, -0.20298213170058702, 2.415413042407908,
            3.3418085707001115, -0.18970019945375155, 2.381408178546708, 3.332377038853407, -0.19905402445374443, 2.3846633503756207,
            3.3567300876539092, -0.18600611246223814, 2.3480869428316984, 3.347369924362723, -0.19527494177277188, 2.3513535411089808,
            3.5741646409204866, -0.2904172453370204, 2.2736196967805697, 3.574035148395789, -0.29897584117770815, 2.2629360501259868,
            3.366069385783489, -0.18480720685347785, 2.324407394813188, 3.3598957798949405, -0.1921977354237334, 2.3233218525572927,
            3.303927131766078, -0.23159973289460514, 2.451438541296758, 3.2932927529522344, -0.2397486162638141, 2.4560019895330947,
            3.527734305139615, -0.3355243100068198, 2.3808409643575192, 3.532911561949785, -0.34175427696258437, 2.3626512743655157,
            3.2981081093725395, -0.22396726686824017, 2.419509265957387, 3.288272230146943, -0.23495775603471966, 2.4197692632099432,
            3.5453136885739696, -0.33323157969033845, 2.3459082008110825, 3.551772854311694, -0.3395601721534377, 2.3286150517241877,
            3.3123841826898075, -0.21972240157759626, 2.3894320134967337, 3.302356433702773, -0.23086791940396675, 2.38984769394357,
            3.5619279549436254, -0.33118761041800454, 2.31247057116793, 3.5668225794087154, -0.33621787373139433, 2.2969315528512806,
            3.3278823251046257, -0.2150528089324073, 2.356935716635524, 3.3175958106162087, -0.22611805579429436, 2.3582953938082967,
            3.576247374397414, -0.32917286067883655, 2.2790474001729164, 3.580993175641437, -0.33410949350302327, 2.263670547023605,
            3.3430376741433805, -0.21107361825942839, 2.323668734322508, 3.332981606731987, -0.22149766041922264, 2.325995738020139,
            3.569163385088081, -0.31176638529892653, 2.2590990251150807, 3.5708029933846053, -0.31598297182274915, 2.2466846884432123,
            3.3538064149111335, -0.20764065393162906, 2.301566643688177, 3.346308472450761, -0.2172477074575349, 2.298646927176111,
            3.428929260633913, -0.3358006263619189, 2.438316250288544, 3.430328313125955, -0.3425691855101243, 2.4183468314709216,
            3.337911156670597, -0.290220060850464, 2.4492065110810204, 3.332175656934985, -0.3005645642591458, 2.449399905744035,
            3.508526209125427, -0.3622599897161765, 2.3902283459446507, 3.510615136599625, -0.36833891167302335, 2.370972282776727,
            3.524560430675355, -0.35803208668111175, 2.3581211667682913, 3.5283727979250417, -0.36535087137553907, 2.337127506536386,
            3.2870521852990926, -0.2469248502745762, 2.392900321830943, 3.281473277334665, -0.25095993425679, 2.3986256307308818,
            3.543183094937875, -0.3554016541197434, 2.324364465046714, 3.544329849614883, -0.36178286553963074, 2.3043572502289464,
            3.3006636336074937, -0.2447660870981636, 2.3594324170245238, 3.295195264429412, -0.2498325288017742, 2.3622250476593396,
            3.557814233454709, -0.35097164688628824, 2.2930688185897696, 3.557021799156096, -0.3566546071631576, 2.2729620983498346,
            3.3149057620506235, -0.24121556597784763, 2.327690770038117, 3.3080013340289427, -0.24871807945514365, 2.3284121128075075,
            3.570775946355546, -0.3468374668316715, 2.2614479794620355, 3.5706569947979387, -0.3536915872757019, 2.241792323213473,
            3.3296157741828964, -0.2374327402592093, 2.2951997751816995, 3.32182644398995, -0.24608539895609938, 2.295535130026906,
            3.5562718306257204, -0.32807172895462894, 2.240126438129817, 3.5516433475209745, -0.33505027805728516, 2.224202841909935,
            3.340323052203459, -0.2326372354264644, 2.276730118978813, 3.333686421346988, -0.2423982260692787, 2.270955963286771,
            3.4506749593978983, -0.3751654334123433, 2.403831362217907, 3.4483713577738135, -0.3798324231349907, 2.386422571212705,
            3.263273301703994, -0.2831517120368759, 2.4327015444101128, 3.258898666269194, -0.29566896910244567, 2.430443593193586,
            3.4999103609543543, -0.3862213671109785, 2.3621486334561843, 3.499729787334884, -0.3928274246612653, 2.3424514259815705,
            3.2638733224607996, -0.266910867594312, 2.4085154780885034, 3.2529004831486343, -0.279127862647503, 2.4089169849432057,
            3.5142469302310726, -0.3822637539174192, 2.3302814150702083, 3.5147259720550696, -0.3895841244813329, 2.3105322218580793,
            3.27895118984767, -0.26336652664222715, 2.374367051315314, 3.267229431124573, -0.2745020916067901, 2.379654791280662,
            3.5299985944260217, -0.37905245335959215, 2.297371098506489, 3.5288443344591283, -0.3857062077297619, 2.2791009013915553,
            3.2925295017075284, -0.2624744757214856, 2.337780475841856, 3.2818643350902166, -0.2719733895793721, 2.3441971705978917,
            3.5423831553492784, -0.37411958018720565, 2.2659803328740775, 3.541097654067738, -0.38024057352855956, 2.246511935349019,
            3.3053188555379065, -0.2617931862253553, 2.302916486233431, 3.295315392394576, -0.270483140787275, 2.3094922186329634,
            3.555617346317593, -0.3696297209483924, 2.233838044576214, 3.556697811020108, -0.37880850481459327, 2.2142796044235786,
            3.3168563743344817, -0.26145578800165364, 2.2707616830422768, 3.3072532139041195, -0.2696610617919812, 2.2774216964794123,
            3.5322883187231517, -0.3505742485678264, 2.215161091318706, 3.5283870367797503, -0.3577650259646182, 2.2006928655895255,
            3.328323712064018, -0.25688544014027737, 2.249546508774693, 3.3192829952444587, -0.2672903433756334, 2.2490168736849907,
            3.4502444305755606, -0.3944859972416671, 2.347111001589077, 3.454133054331134, -0.4025712243388927, 2.3271763032006385,
            3.5024746904700605, -0.40706803467991604, 2.3010962435532876, 3.5033408635191297, -0.4142787101465669, 2.2822389139240165,
            3.2603711381144564, -0.29235243398194016, 2.3539925474510555, 3.250104805607538, -0.30275468010603, 2.356976466883395,
            3.514437885136987, -0.40344855431811055, 2.2734046951138738, 3.5136552235365284, -0.41159730743459916, 2.2570231730680983,
            3.2742642066370924, -0.29075122232226286, 2.318304441888601, 3.264492768384429, -0.29999470176049514, 2.3228120787764177,
            3.2899364508632694, -0.2855724222626368, 2.286601844507419, 3.2795504317230995, -0.29634158994107895, 2.2889973730070206,
            3.303619493321269, -0.2819430672735983, 2.256659853238685, 3.296306193423749, -0.2891644731835839, 2.259264268867157,
            3.3130717393373508, -0.28160058694039325, 2.230484448480384, 3.3070113978445237, -0.288048037048623, 2.2314674145969424,
            3.2589817908587206, -0.321786088410403, 2.3356662729347617, 3.2526657915015496, -0.33169477232008315, 2.3356931553274283,
            3.261539584089884, -0.31486755959731527, 2.2935302968562428, 3.2535446138135593, -0.3218310131181649, 2.298750484301373,
            3.2757562204150332, -0.31184826162528617, 2.260513960440341, 3.267516811034482, -0.321130751843994, 2.260539247401748,
            3.2927756711345237, -0.30418236812337657, 2.2312663930778114, 3.282611717265979, -0.3160969888873483, 2.2301205561513933,
            3.3032074143228516, -0.30170768856051644, 2.20769774380118, 3.2962849875103704, -0.31062061372610816, 2.2048923066535755,
            3.3602358413474036, -0.4145681304834067, 2.317424116671325, 3.365723747628809, -0.4229428273472699, 2.2962501350815336,
            3.2973787790311126, -0.3813229357987727, 2.3185118951666723, 3.2906842025729306, -0.3911722782136358, 2.3191341059895008,
            3.442801127930312, -0.4365124102769159, 2.2652887865378113, 3.44884691730461, -0.44401342576961533, 2.2481718987408867,
            3.2343927676124946, -0.3423411402740297, 2.3015143601109176, 3.224117229853102, -0.3552452596084314, 2.298166532636257,
            3.4768760954368374, -0.4551639337467811, 2.2240236621207563, 3.4810830278513745, -0.4611462415282665, 2.2066805888390153,
            3.250043530056325, -0.3373590810540483, 2.2693741100522136, 3.2392846329289005, -0.34819972354000633, 2.2726551702179094,
            3.490823391959536, -0.4530432958676853, 2.1889145303193773, 3.496830400322915, -0.4587794664015692, 2.172657242340591,
            3.2649867326965563, -0.33461810013665066, 2.233561787320359, 3.255508421994723, -0.3438829758263372, 2.237176462525694,
            3.491039695820283, -0.4288315471648092, 2.1650129884009215, 3.4960354532163804, -0.4371671597301054, 2.1477076082086093,
            3.2791020660236807, -0.3314246109669634, 2.201277192821276, 3.26971649083481, -0.34066323199667936, 2.2046931370288347,
            3.5024534289021756, -0.43093197973479047, 2.1351187479310143, 3.5040479983613926, -0.4395971024203524, 2.1227968257683356,
            3.289490146022143, -0.3282568677138849, 2.1795911048803744, 3.2824711784039993, -0.33714624419052386, 2.177121827974493,
            3.440942927029452, -0.4830845298177863, 2.2698042322030285, 3.4433466023235333, -0.4886156465777005, 2.250607989387939,
            3.207477417626307, -0.3726413542151795, 2.3016430488172794, 3.2007591279984746, -0.37850642116837263, 2.3060010856392275,
            3.456010083853949, -0.47806527117093517, 2.2375890965425396, 3.4619227527546728, -0.4854772542384291, 2.2167749445617386,
            3.2201178846075624, -0.3712298446691266, 2.269057593576848, 3.2129582829958774, -0.379101604270386, 2.269572232412817,
            3.4754005888182267, -0.4764035305014058, 2.2031465705730056, 3.4788491870590446, -0.4827570694235918, 2.1834514286245366,
            3.2351722405219894, -0.3655703274574699, 2.240342358794216, 3.2262810299600417, -0.37569632925173063, 2.240092609627815,
            3.4905386344564326, -0.47323595098397425, 2.170728158911761, 3.4927258884497925, -0.47828009669850957, 2.1521179003368087,
            3.2503469760088968, -0.361047292413149, 2.2083995610191884, 3.240769691663587, -0.3718940015587971, 2.2082843971220147,
            3.5037551919520715, -0.47028506350545046, 2.1381105156189415, 3.5056145803946204, -0.47539594884531783, 2.1199425434358146,
            3.265448883996619, -0.356890202634372, 2.1757367971372163, 3.2555576798725387, -0.36800938736198935, 2.175828545922234,
            3.4949786052001253, -0.4521055753414224, 2.1168713464546474, 3.4919451971020483, -0.45708212596245273, 2.101989343835248,
            3.2765270230548906, -0.3523390915315227, 2.1555861374699545, 3.268890303836381, -0.36284442894669106, 2.15078516242512,
            3.3867648894470617, -0.4951234828929994, 2.2833769630591343, 3.38597240329806, -0.5009044819416859, 2.2658459000455946,
            3.2258373739838415, -0.4175431463818835, 2.3142916323338683, 3.220863800171024, -0.4286259756304911, 2.3108891967830285,
            3.4342884321818623, -0.5064409181449898, 2.2420396246065772, 3.4350324836440533, -0.5127631161782645, 2.2227491693032526,
            3.1976807894995267, -0.3917012293047293, 2.281334399967019, 3.1882209610686183, -0.400186916556477, 2.286872908411414,
            3.4500033226877487, -0.5023498276781059, 2.2105033821038145, 3.4516309820597826, -0.5099091750359112, 2.1903430827042425,
            3.2112600522419186, -0.3895688595611516, 2.2478763646126865, 3.202700848373404, -0.3963034176421588, 2.2552804953449406,
            3.531191185856916, -0.29306317286901484, 2.3669639681930676, 3.5376388115070023, -0.30024155547126513, 2.3521081909540187,
            3.4669990103676778, -0.49939014051857933, 2.176919503213699, 3.4663198234890444, -0.5060378486439091, 2.1571412285836113,
            3.2249129005844956, -0.38706623140119595, 2.2151623819827626, 3.2170515244780566, -0.39306700210838175, 2.2224314981316646,
            3.549982103174955, -0.291350664603946, 2.335873211944739, 3.55365395291844, -0.29801210821728324, 2.3212203081094755,
            3.4803572541742263, -0.49444892326360496, 2.145062895419498, 3.478658509886043, -0.5006406065491477, 2.1254728639496108,
            3.2397762773879246, -0.3834967036413213, 2.1816914962245844, 3.23281069687078, -0.3885429316883437, 2.188819230251732,
            3.5665641563612382, -0.2891820927084907, 2.3018228103591687, 3.56807355882043, -0.2981016174768686, 2.2860619653266094,
            3.492886532604843, -0.4895826653236686, 2.114359896683796, 3.4927610204035333, -0.49798047937289325, 2.0938797718905477,
            3.2550367297799165, -0.37975965220869895, 2.1475095183598265, 3.2488326116406117, -0.38440025934139677, 2.1534876325593495,
            3.473251507299553, -0.4709151660606077, 2.094986799173027, 3.4680341889103707, -0.4785535017552621, 2.079892437836569,
            3.2630484635839956, -0.37804987660454215, 2.128924509611126, 3.259086093651611, -0.3838459151387431, 2.1255571167608145,
            3.3753170334822107, -0.520720632766392, 2.255831738034618, 3.372680034020975, -0.5254098327450998, 2.2388098654441846,
            3.2091090144629986, -0.44039868600535975, 2.2877477625269256, 3.2078611738846887, -0.4522813734953187, 2.281980328810475,
            3.4877009878836467, -0.31992708243885926, 2.4184723373642343, 3.494690190062251, -0.3273260026090885, 2.3987835692613406,
            3.4237447421109923, -0.531164250384141, 2.213951592825392, 3.4231768686470136, -0.5380838918649825, 2.1961808370124714,
            3.18320655606674, -0.4153448587247456, 2.2627500685040345, 3.1738304745656727, -0.42580578814393116, 2.2630379743487947,
            3.436092269553063, -0.5272421050169149, 2.18481811409076, 3.436455597776981, -0.5344781864312637, 2.165162710797943,
            3.197424306335438, -0.411590589828445, 2.2315951073258837, 3.1873898424596856, -0.422146343659709, 2.2335261720405826,
            3.4513151234047004, -0.5237929645758537, 2.1509094047813386, 3.450812883760083, -0.530377296705505, 2.1315333554731195,
            3.211935457624406, -0.4084047860164768, 2.19815854188853, 3.2008674433503033, -0.42019381565615427, 2.199918020350495,
            3.464118839555934, -0.5188138640593989, 2.1181751478055046, 3.464146350967776, -0.5269050418006752, 2.1000516539884777,
            3.229968073092361, -0.40079775575666754, 2.1658628564147255, 3.2169542099524255, -0.41465883937864406, 2.1679330519026214,
            3.4431847216471585, -0.49986259211933337, 2.095674866582955, 3.4387127252425507, -0.5067332040819446, 2.077915563949279,
            3.244995017012623, -0.39792131225702587, 2.1301657063679853, 3.233096058981519, -0.40927258062497995, 2.1354131961962137,
            3.4501036065860737, -0.49626324393139004, 2.068422302837082, 3.4478305665563593, -0.5022932763354442, 2.0533835593126377,
            3.253863719462059, -0.39743220512528193, 2.1060312877846106, 3.245976760962378, -0.40686699511147834, 2.1046617631922593,
            3.3412676287588763, -0.525313768173953, 2.2086566033360238, 3.3431066331123627, -0.5332111308948917, 2.194856448480559,
            3.27404151677789, -0.2490860654194049, 2.4246430128236516, 3.2681633206271217, -0.2539014337922314, 2.4292451202159273,
            3.181662359400788, -0.445830821020756, 2.2386829817753653, 3.175244992256003, -0.45583688387270865, 2.235601927031185,
            3.387710438259003, -0.5358747782575647, 2.1701520716690803, 3.3866962786528925, -0.5425403747972565, 2.149997825918916,
            3.182019459248104, -0.43914344160881913, 2.205771388660743, 3.1736830858424936, -0.4475090807788621, 2.2084000798109997,
            3.4019215421514595, -0.5322407868812504, 2.1361440106600535, 3.4005486695157754, -0.5392354257406389, 2.115723372850765,
            3.196311285960958, -0.4356805748287496, 2.1736651895519943, 3.1874977657345114, -0.4440744469501463, 2.1775874931475436,
            3.4152685092558297, -0.5285944879915453, 2.102912097088704, 3.4138676450098266, -0.5356425904468053, 2.0807832567657614,
            3.2142809040092457, -0.42765288854592154, 2.142616912822674, 3.204705994291879, -0.43688559435660884, 2.146589572010081,
            3.425431386484192, -0.5250532155606386, 2.0679419379111863, 3.4278255254340095, -0.5356861727965311, 2.041298312459536,
            3.227434847015951, -0.4274434295933533, 2.1055127340879918, 3.2216335551836135, -0.4319626538915523, 2.110646255226557,
            3.4124100835164377, -0.5131997393421158, 2.0435944297766646, 3.405306228661312, -0.516306037476944, 2.0314490994669074,
            3.237083936770696, -0.425155565681055, 2.0837093920367424, 3.230440980855375, -0.4343807317965651, 2.0793128731042465,
            3.3733229899988024, -0.5580932839911, 2.1453383875983185, 3.3782360415809127, -0.5643517118725236, 2.1270950367117787,
            3.1697076031403695, -0.463556426033523, 2.1790636926532505, 3.163119897759731, -0.47038728965230425, 2.1805827450544437,
            3.3902364575976014, -0.55579344153814, 2.111333837223087, 3.396608669195898, -0.5628138140298564, 2.092581940010575,
            3.1840416281277584, -0.4588521396303353, 2.1499861669598195, 3.1745878389392286, -0.46996847954438525, 2.1488336753152697,
            3.4055200168661135, -0.5533039006882687, 2.074842635315311, 3.411146971447112, -0.5600660696127538, 2.056937900330616,
            3.2025360275235255, -0.4500030009407602, 2.119520323849455, 3.1906851815858506, -0.4635756740940678, 2.118994615203644,
            3.4218887385393635, -0.5530377166070236, 2.036234896329417, 3.426646477320224, -0.5629067037677632, 2.0190222460009934,
            3.2177438995035357, -0.4476658888259515, 2.0819372397017006, 3.2082212815608138, -0.4564735300638468, 2.086838649137192,
            3.4333238694374058, -0.5558817993058116, 2.0068197976274584, 3.434739141464806, -0.5628801501633274, 1.9970149209473718,
            3.225483162270653, -0.4490109078377473, 2.056382470251699, 3.2185873727511516, -0.4577768405699184, 2.053873315710587,
            3.3482263732021127, -0.58515948902822, 2.1543854396848747, 3.3539391295421317, -0.5924630735880099, 2.133896064320307,
            3.245990259305489, -0.29712133278437836, 2.383040217567703, 3.2371413126876827, -0.307099005187462, 2.3830458794283165,
            3.1623979506670765, -0.49531523974601166, 2.1827220866554295, 3.1519706964798218, -0.5042074968939916, 2.186889474294559,
            3.3880915096486057, -0.6000846413847046, 2.1148019098744344, 3.392527073121864, -0.6063017422312336, 2.0961078234576886,
            3.158131180501745, -0.4876901562067774, 2.1509602307679647, 3.148563025568234, -0.497831169129954, 2.152609267578523,
            3.40514338553374, -0.597747391776854, 2.080575196350889, 3.410921906203553, -0.604170607177329, 2.062104731449239,
            3.170995981039743, -0.48568010573285936, 2.119251322842416, 3.161127027982899, -0.49657091877667997, 2.119858825893918,
            3.4934320695177488, -0.38443218919126687, 2.2497846097901855, 3.489286947082379, -0.3931275483363199, 2.229741197569251,
            3.4211736370573624, -0.5957724416991953, 2.0465130401225524, 3.425654687930756, -0.6008570826173929, 2.029876945109199,
            3.1863354214053348, -0.48042905950650966, 2.0886841934713267, 3.1761000100477563, -0.4914171512305244, 2.090093392754649,
            3.501318173263603, -0.3818428432513188, 2.2159811914491643, 3.5011167000076355, -0.391387092697932, 2.192032407223539,
            3.435191010133242, -0.5938961249170054, 2.012323308516795, 3.439710999144492, -0.5988645384233552, 1.9964405264928,
            3.2015460999877154, -0.4762421231061251, 2.0557859359414152, 3.1915110857699096, -0.4866249793944542, 2.0581571889213617,
            3.4909676110309347, -0.3673718751668159, 2.186925246418183, 3.4856305128879015, -0.3761696164034119, 2.170702882543168,
            3.428499209084511, -0.575687516498915, 1.9928433094488622, 3.429508892061789, -0.5801886884434425, 1.979791649658269,
            3.2123843329361774, -0.4726188185982458, 2.0339678742275598, 3.2048656781154623, -0.48227645989819995, 2.030979085022896,
            3.2874484273323583, -0.5965658177005143, 2.1772733904335992, 3.2905426460328826, -0.598960094509606, 2.1606330385040287,
            3.2021442347397655, -0.5522468661252651, 2.183174916281096, 3.2008779328828085, -0.5595323407540019, 2.185793478419963,
            3.3692006080830845, -0.623660202788668, 2.126164196331072, 3.3710978953944846, -0.6258249843780792, 2.110547963550143,
            3.1326699301339715, -0.5127109326799615, 2.160331884706737, 3.127478624018672, -0.5155585941203563, 2.1679608304880387,
            3.4542070575851738, -0.4145744106160152, 2.288235535623827, 3.454410895260458, -0.4210810788352313, 2.2702386926042983,
            3.383398396165532, -0.6188490365252813, 2.094112363958514, 3.3868937166127493, -0.6223847680966397, 2.0771062925703068,
            3.1455484828949656, -0.5100164394766035, 2.130321915101865, 3.1392827252827233, -0.5131710569885084, 2.1402449846269676,
            3.467214496002156, -0.41139630258373017, 2.260067134006638, 3.4690153942952517, -0.41802857575816776, 2.239788262142867,
            3.4018009142395695, -0.6162606892025775, 2.060631287167831, 3.403033664972109, -0.6190670550997535, 2.044203984911843,
            3.160015072046171, -0.5062423668839051, 2.0985045908384135, 3.154004681891748, -0.5095468746602017, 2.1073186432358026,
            3.477796794638624, -0.41012900660494284, 2.2230353913015617, 3.4795050783584256, -0.4169437168766823, 2.202601089583751,
            3.4166357078103218, -0.6120616139921997, 2.0290413341798983, 3.415899579461432, -0.614206203912159, 2.0121971091014017,
            3.1748898166476733, -0.5018352642219898, 2.067125975555413, 3.168673831702932, -0.5057696668305335, 2.0749296705876756,
            3.4907937575456787, -0.4094865080544312, 2.185358111625816, 3.4973649539061897, -0.4150023602980106, 2.16549302604352,
            3.429438548594685, -0.6079700860430141, 1.9971499380067075, 3.42957120198721, -0.6114349413223586, 1.9809133551905422,
            3.190193456014919, -0.49728892900595556, 2.0348734993145614, 3.1838075982502603, -0.5015111221085934, 2.0424330866652904,
            3.477231596598426, -0.3901305659820917, 2.162227308991889, 3.4756564944616875, -0.39564808686978087, 2.153265843482395,
            3.412503662484033, -0.5889803106096688, 1.977414190509777, 3.4089863274379386, -0.5931250019348235, 1.9643326134444186,
            3.2010547210838984, -0.49357022114751503, 2.0132329779870477, 3.196417513201161, -0.4981998339668314, 2.0147561392987656
        };
    if ( !vectorTools::fuzzyEquals( *dc.getBoundaryPoints( ), boundaryPointAnswer ) ){

        results << "test_dualContouring_localBoundary (test 1) & False\n";

        return 1;

    }

    floatType volumeAnswer = 0.0199049;

    floatVector volumeResult;

    error = dc.performVolumeIntegration( floatVector( points.size( ) / 3, 1 ), 1, volumeResult );

    if ( error ){

        error->print( );
        
        results << "test_dualContouring_localBoundary & False\n";

        return 1;

    }

    if ( !vectorTools::fuzzyEquals( volumeResult[ 0 ], volumeAnswer ) ){

        std::cerr << "areaResult: " << volumeResult[ 0 ] << "\n";
        results << "test_dualContouring_localBoundary (test 2) & False\n";

        return 1;

    }

    floatType areaAnswer = 0.583317;

    floatVector areaResult;

    error = dc.performSurfaceIntegration( floatVector( points.size( ) / 3, 1 ), 1, areaResult );

    if ( error ){

        error->print( );
        
        results << "test_dualContouring_localBoundary & False\n";

        return 1;

    }

    if ( !vectorTools::fuzzyEquals( areaResult[ 0 ], areaAnswer ) ){

        std::cerr << "areaResult: " << areaResult[ 0 ] << "\n";
        results << "test_dualContouring_localBoundary (test 3) & False\n";

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
    test_dualContouring_planes( results );
    test_dualContouring_localBoundary( results );

    test_KDNode_constructor( results );
    test_KDNode_getIndex( results );
    test_KDNode_getMinimumValueDimension( results );
    test_KDNode_getMaximumValueDimension( results );
    test_KDNode_getPointsInRange( results );
    test_KDNode_getPointsWithinRadiusOfOrigin( results );

    //Close the results file
    results.close();
}
