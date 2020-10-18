//!The test file for inputFileProcessor.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>

#include<inputFileProcessor.h>
#include<generateXDMFData.h>

typedef inputFileProcessor::errorNode errorNode; //!Redefinition for the error node
typedef inputFileProcessor::errorOut errorOut; //!Redefinition for a pointer to the error node
typedef inputFileProcessor::floatType floatType; //!Define the float values type.
typedef inputFileProcessor::floatVector floatVector; //! Define a vector of floats
typedef inputFileProcessor::floatMatrix floatMatrix; //!Define a matrix of floats
typedef inputFileProcessor::uIntType uIntType; //!Define an unsigned integer type
typedef inputFileProcessor::uIntVector uIntVector; //!Define a vector of unsigned ints
typedef inputFileProcessor::stringVector stringVector; //!Define a vector of strings
typedef inputFileProcessor::DOFMap DOFMap; //!Define the map between DOF values

errorOut _createXDMFDatafiles( ){

    fileGenerator::fileGenerator fg( "macroscale.yaml" );
    if ( fg.build( ) ){

        fg.getError( )->print( );
        errorOut result = new errorNode( "_createXDMFDatafiles", "Error in creation of the macroscale datafile" );
        return result;

    }

    fg = fileGenerator::fileGenerator( "microscale.yaml" );
    if ( fg.build( ) ){

        fg.getError( )->print( );
        errorOut result = new errorNode( "_createXDMFDatafiles", "Error in creation of the microscale datafile" );
        return result;

    }

    return NULL;

}

int test_openConfigurationFile( std::ofstream &results ){
    /*!
     * Test opening the YAML configuration file.
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print();
        results << "test_openConfigurationFile (test 1) & False\n";
        return 1;
    }

    reader = inputFileProcessor::inputFileProcessor( );
    std::unique_ptr< errorNode > error;
    error.reset( reader.setConfigurationFilename( "" ) );
    if ( !error ){
        results << "test_openConfigurationFile (test 2) & False\n";
        return 1;
    }

    error.reset( reader.setConfigurationFilename( filename ) );
    if ( error ){
        error->print();
        results << "test_openConfigurationFile (test 3) & False\n";
        return 1;
    }

    //Check the variable configuration
    const uIntVector *result = reader.getFreeMacroCellIds( );

    if ( result->size( ) == 0 ){

        results << "test_openConfigurationFile (test 4) & False\n";
        return 1;

    }

    if ( ( *result )[ 0 ] != 1 ){

        results << "test_openConfigurationFile (test 5) & False\n";
        return 1;

    }

    result = reader.getGhostMacroCellIds( );

    if ( result->size( ) == 0 ){

        results << "test_openConfigurationFile (test 6) & False\n";
        return 1;

    }

    if ( ( *result )[ 0 ] != 2 ){

        results << "test_openConfigurationFile (test 7) & False\n";
        return 1;

    }

//    result = reader.getFreeMacroCellMicroDomainCounts( );
//
//    if ( result->size( ) == 0 ){
//
//        results << "test_openConfigurationFile (test 8) & False\n";
//        return 1;
//
//    }
//
//    for ( auto v = result->begin( ); v != result->end( ); v++ ){
//
//        if ( *v != 8 ){
//
//            results << "test_openConfigurationFile (test 9) & False\n";
//            return 1;
//
//        }
//
//    }
//
//    result = reader.getGhostMacroCellMicroDomainCounts( );
//
//    if ( result->size( ) == 0 ){
//
//        results << "test_openConfigurationFile (test 10) & False\n";
//        return 1;
//
//    }
//
//    for ( auto v = result->begin( ); v != result->end( ); v++ ){
//
//        if ( *v != 8 ){
//
//            results << "test_openConfigurationFile (test 11) & False\n";
//            return 1;
//
//        }
//
//    }

    results << "test_openConfigurationFile & True\n";
    return 0;
}

int test_setConfigurationFile( std::ofstream &results ){
    /*!
     * Test setting the YAML configuration file.
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader;

    std::unique_ptr< errorNode > error;
    error.reset( reader.setConfigurationFilename( "" ) );

    if ( !error ){
        results << "test_setConfigurationFile & False\n";
        return 1;
    }

    error.reset( reader.setConfigurationFilename( filename ) );

    if ( error ){
        error->print();
        results << "test_setConfigurationFile (test 1) & False\n";
        return 1;
    }

    results << "test_setConfigurationFile & True\n";
    return 0;
}

int test_initializeFileInterfaces( std::ofstream &results ){
    /*!
     * Test the initialization of the file readers
     *
     * :param std::ofstream &results: The output file.
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print();
        results << "test_initializeFileInterfaces & False\n";
        return 1;
    }

    floatVector answerMacroNodes =
        {
            0,  0, -1,  1,  0, -1,  1,  1, -1,  0,  1, -1,  0,  0,  0,  1,  0,
            0,  1,  1,  0,  0,  1,  0,  0,  0,  1,  1,  0,  1,  1,  1,  1,  0,
            1,  1,  0,  0,  2,  1,  0,  2,  1,  1,  2,  0,  1,  2
        };

    floatVector resultMacroNodes, resultMicroNodes;

    errorOut error = reader._macroscale->readMesh( 1, resultMacroNodes );

    if ( error ){
        error->print();
        results << "test_initializeFileInterfaces & False\n";
        return 1;
    }

    floatVector answerMicroNodes =
        {
            0. , 0. , 0. , 0. , 0. , 0.5, 0. , 0. , 1. , 0. , 0. , 1.5, 0. ,
            0. , 2. , 0. , 0. , 2.5, 0. , 0. , 3. , 0.5, 0. , 0. , 0.5, 0. ,
            0.5, 0.5, 0. , 1. , 0.5, 0. , 1.5, 0.5, 0. , 2. , 0.5, 0. , 2.5,
            0.5, 0. , 3. , 1. , 0. , 0. , 1. , 0. , 0.5, 1. , 0. , 1. , 1. ,
            0. , 1.5, 1. , 0. , 2. , 1. , 0. , 2.5, 1. , 0. , 3. , 0. , 0.5,
            0. , 0. , 0.5, 0.5, 0. , 0.5, 1. , 0. , 0.5, 1.5, 0. , 0.5, 2. ,
            0. , 0.5, 2.5, 0. , 0.5, 3. , 0.5, 0.5, 0. , 0.5, 0.5, 0.5, 0.5,
            0.5, 1. , 0.5, 0.5, 1.5, 0.5, 0.5, 2. , 0.5, 0.5, 2.5, 0.5, 0.5,
            3. , 1. , 0.5, 0. , 1. , 0.5, 0.5, 1. , 0.5, 1. , 1. , 0.5, 1.5,
            1. , 0.5, 2. , 1. , 0.5, 2.5, 1. , 0.5, 3. , 0. , 1. , 0. , 0. ,
            1. , 0.5, 0. , 1. , 1. , 0. , 1. , 1.5, 0. , 1. , 2. , 0. , 1. ,
            2.5, 0. , 1. , 3. , 0.5, 1. , 0. , 0.5, 1. , 0.5, 0.5, 1. , 1. ,
            0.5, 1. , 1.5, 0.5, 1. , 2. , 0.5, 1. , 2.5, 0.5, 1. , 3. , 1. ,
            1. , 0. , 1. , 1. , 0.5, 1. , 1. , 1. , 1. , 1. , 1.5, 1. , 1. ,
            2. , 1. , 1. , 2.5, 1. , 1. , 3.
        };

    if ( !vectorTools::fuzzyEquals( answerMacroNodes, resultMacroNodes ) ){
        results << "test_initializeFileInterfaces (test 1) & False\n";
        return 1;
    }

    error = reader._microscale->readMesh( 1, resultMicroNodes );

    if ( error ){
        error->print();
        results << "test_initializeFileInterfaces & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultMicroNodes, answerMicroNodes ) ){
        results << "test_initializeFileInterfaces (test 2) & False\n";
        return 1;
    }

    results << "test_initializeFileInterfaces & True\n";
    return 0;
}

int test_initializeIncrement( std::ofstream &results ){
    /*!
     * Test the initialization of the processor for the current increment
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_initializeIncrement & False\n";
        return 1;
    }

    errorOut error = reader.initializeIncrement( 1, 1 );
    if ( error ){
        error->print( );
        results << "test_initializeIncrement & False\n";
        return 1;
    }

    //Check that the unique micro-scale nodes have been identified
    const DOFMap microGlobalToLocalMapAnswer
        =
            {
                { 15,  0 },
                { 31,  1 },
                { 13,  2 },
                { 26,  3 },
                { 53,  4 },
                { 21,  5 },
                { 37,  6 },
                { 48,  7 },
                {  5,  8 },
                { 10,  9 },
                {  3, 10 },
                {  4, 11 },
                { 32, 12 },
                { 33, 13 },
                { 34, 14 },
                { 28, 15 },
                { 25, 16 },
                { 50, 17 },
                { 43, 18 },
                { 27, 19 },
                {  1, 20 },
                {  7, 21 },
                { 30, 22 },
                { 16, 23 },
                { 22, 24 },
                {  2, 25 },
                { 46, 26 },
                { 24, 27 },
                { 39, 28 },
                { 40, 29 },
                { 57, 30 },
                { 44, 31 },
                { 58, 32 },
                { 29, 33 },
                { 59, 34 },
                { 11, 35 },
                {  0, 36 },
                { 20, 37 },
                { 60, 38 },
                { 47, 39 },
                { 49, 40 },
                { 17, 41 },
                { 38, 42 },
                { 14, 43 },
                { 55, 44 },
            };

    const DOFMap *microGlobalToLocalResult = reader.getMicroGlobalToLocalDOFMap( );

    for ( auto it = microGlobalToLocalMapAnswer.begin( ); it != microGlobalToLocalMapAnswer.end( ); it++ ){

        auto r = microGlobalToLocalResult->find( it->first );

        if ( r == microGlobalToLocalResult->end( ) ){

            results << "test_initializeIncrement (test 1) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement (test 2) & False\n";
            return 1;

        }

    }

    //Check that the unique macro-scale nodes have been identified
    const DOFMap macroGlobalToLocalMapAnswer
        =
            {
                {  5,  0 },
                {  9,  1 },
                {  8,  2 },
                { 11,  3 },
                {  3,  4 },
                {  1,  5 },
                {  6,  6 },
                { 15,  7 },
                { 12,  8 },
                {  2,  9 },
                { 13, 10 },
                { 14, 11 },
            };

    const DOFMap *macroGlobalToLocalResult = reader.getMacroGlobalToLocalDOFMap( );

    for ( auto it = macroGlobalToLocalMapAnswer.begin( ); it != macroGlobalToLocalMapAnswer.end( ); it++ ){

        auto r = macroGlobalToLocalResult->find( it->first );

        if ( r == macroGlobalToLocalResult->end( ) ){

            results << "test_initializeIncrement (test 3) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement (test 4) & False\n";
            return 1;

        }

    }

    //Check that the micro node weights are initialized properly
    const std::unordered_map< uIntType, floatType > microNodeWeightsAnswer
        =
            {
                { 24, 1.000 },
                { 39, 0.500 },
                { 15, 0.500 },
                { 31, 0.500 },
                { 43, 1.000 },
                { 40, 0.500 },
                { 57, 0.250 },
                { 13, 0.250 },
                { 26, 0.250 },
                { 27, 0.500 },
                { 11, 1.000 },
                {  0, 0.500 },
                {  5, 0.500 },
                { 10, 0.500 },
                { 30, 1.000 },
                { 44, 0.500 },
                { 58, 0.250 },
                { 53, 0.250 },
                { 21, 0.250 },
                {  1, 0.500 },
                { 29, 0.250 },
                { 59, 0.125 },
                { 37, 0.125 },
                { 48, 0.125 },
                {  7, 0.250 },
                { 20, 0.500 },
                { 60, 0.250 },
                {  3, 0.250 },
                {  4, 0.250 },
                { 16, 0.500 },
                { 14, 1.000 },
                { 55, 0.500 },
                { 25, 0.500 },
                { 50, 0.500 },
                { 46, 1.000 },
                { 47, 0.500 },
                { 49, 0.250 },
                { 32, 0.250 },
                { 33, 0.250 },
                { 22, 0.500 },
                { 17, 1.000 },
                { 38, 0.500 },
                { 34, 0.500 },
                { 28, 0.500 },
                {  2, 1.000 },
            };

    const std::unordered_map< uIntType, floatType > *microNodeWeightsResult = reader.getMicroWeights( );

    for ( auto it = microNodeWeightsAnswer.begin( ); it != microNodeWeightsAnswer.end( ); it++ ){

        auto r = microNodeWeightsResult->find( it->first );

        if ( r == microNodeWeightsResult->end( ) ){

            results << "test_initializeIncrement (test 5) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement (test 6) & False\n";
            return 1;

        }

    }

    //Make sure the micro global node id to output index map has been extracted correctly
    const DOFMap microGlobalNodeToOutputMapAnswer
        =
            {
                { 15,  2 },
                { 31,  3 },
                { 13,  9 },
                { 26, 10 },
                { 53, 23 },
                { 21, 24 },
                { 37, 30 },
                { 48, 31 },
                {  5, 16 },
                { 10, 17 },
                {  3, 37 },
                {  4, 38 },
                { 32, 51 },
                { 33, 52 },
                { 34, 58 },
                { 28, 59 },
                { 25, 44 },
                { 50, 45 },
                { 43,  4 },
                { 27, 11 },
                {  1, 25 },
                {  7, 32 },
                { 30, 18 },
                { 16, 39 },
                { 22, 53 },
                {  2, 60 },
                { 46, 46 },
                { 24,  0 },
                { 39,  1 },
                { 40,  7 },
                { 57,  8 },
                { 44, 21 },
                { 58, 22 },
                { 29, 28 },
                { 59, 29 },
                { 11, 14 },
                {  0, 15 },
                { 20, 35 },
                { 60, 36 },
                { 47, 49 },
                { 49, 50 },
                { 17, 56 },
                { 38, 57 },
                { 14, 42 },
                { 55, 43 },
            };

    const DOFMap *microGlobalNodeToOutputMapResult = reader.getMicroNodeIDOutputIndex( );

    for ( auto it = microGlobalNodeToOutputMapAnswer.begin( ); it != microGlobalNodeToOutputMapAnswer.end( ); it++ ){

        auto r = microGlobalNodeToOutputMapResult->find( it->first );

        if ( r == microGlobalNodeToOutputMapResult->end( ) ){

            results << "test_initializeIncrement (test 7) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement (test 8) & False\n";
            return 1;

        }

    }

    //Make sure the macro global node id to output index map has been extracted correctly
    const DOFMap macroGlobalNodeToOutputMapAnswer
        =
            {
                {  5,  4 },
                {  9,  5 },
                {  8,  6 },
                { 11,  7 },
                {  3,  8 },
                {  1,  9 },
                {  6, 10 },
                { 15, 11 },
                { 12, 12 },
                {  2, 13 },
                { 13, 14 },
                { 14, 15 },
            };

    const DOFMap *macroGlobalNodeToOutputMapResult = reader.getMacroNodeIDOutputIndex( );

    for ( auto it = macroGlobalNodeToOutputMapAnswer.begin( ); it != macroGlobalNodeToOutputMapAnswer.end( ); it++ ){

        auto r = macroGlobalNodeToOutputMapResult->find( it->first );

        if ( r == macroGlobalNodeToOutputMapResult->end( ) ){

            results << "test_initializeIncrement (test 9) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement (test 10) & False\n";
            return 1;

        }

    }

    //Make sure the time of the micro increment has been extracted correctly
    const floatType timeAnswer = 1;

    const floatType *timeResult = reader.getMicroTime( );

    if ( !vectorTools::fuzzyEquals( timeAnswer, *timeResult ) ){

        results << "test_initializeIncrement (test 11) & False\n";
        return 1;

    }

    const std::unordered_map< uIntType, floatType > densityAnswer
        =
        {
            { 15, 6.000 },
            { 31, 8.000 },
            { 13, 7.000 },
            { 26, 9.000 },
            { 53, 4.500 },
            { 21, 6.500 },
            { 37, 5.500 },
            { 48, 7.500 },
            {  5, 8.000 },
            { 10, 10.000 },
            {  3, 6.500 },
            {  4, 8.500 },
            { 32, 4.000 },
            { 33, 6.000 },
            { 34, 5.000 },
            { 28, 7.000 },
            { 25, 3.000 },
            { 50, 5.000 },
            { 43, 10.000 },
            { 27, 11.000 },
            {  1, 8.500 },
            {  7, 9.500 },
            { 30, 12.000 },
            { 16, 10.500 },
            { 22, 8.000 },
            {  2, 9.000 },
            { 46, 7.000 },
            { 24, 2.000 },
            { 39, 4.000 },
            { 40, 3.000 },
            { 57, 5.000 },
            { 44, 0.500 },
            { 58, 2.500 },
            { 29, 1.500 },
            { 59, 3.500 },
            { 11, 4.000 },
            {  0, 6.000 },
            { 20, 2.500 },
            { 60, 4.500 },
            { 47, 0.000 },
            { 49, 2.000 },
            { 17, 1.000 },
            { 38, 3.000 },
            { 14, -1.000 },
            { 55, 1.000 },
        };

    const std::unordered_map< uIntType, floatType > *densityResult = reader.getMicroDensities( );

    for ( auto it = densityAnswer.begin( ); it != densityAnswer.end( ); it++ ){

        auto r = densityResult->find( it->first );

        if ( r == densityResult->end( ) ){

            results << "test_initializeIncrement (test 12) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement (test 13) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatType > volumeAnswer
        =
         {
            { 15, -2.000 },
            { 31, -3.500 },
            { 13, 0.500 },
            { 26, -1.000 },
            { 53, -1.800 },
            { 21, -3.300 },
            { 37, 0.700 },
            { 48, -0.800 },
            {  5, 3.000 },
            { 10, 1.500 },
            {  3, 3.200 },
            {  4, 1.700 },
            { 32, 0.900 },
            { 33, -0.600 },
            { 34, 3.400 },
            { 28, 1.900 },
            { 25, -1.600 },
            { 50, -3.100 },
            { 43, -5.000 },
            { 27, -2.500 },
            {  1, -4.800 },
            {  7, -2.300 },
            { 30, 0.000 },
            { 16, 0.200 },
            { 22, -2.100 },
            {  2, 0.400 },
            { 46, -4.600 },
            { 24, 1.000 },
            { 39, -0.500 },
            { 40, 3.500 },
            { 57, 2.000 },
            { 44, 1.200 },
            { 58, -0.300 },
            { 29, 3.700 },
            { 59, 2.200 },
            { 11, 6.000 },
            {  0, 4.500 },
            { 20, 6.200 },
            { 60, 4.700 },
            { 47, 3.900 },
            { 49, 2.400 },
            { 17, 6.400 },
            { 38, 4.900 },
            { 14, 1.400 },
            { 55, -0.100 },
        };

    const std::unordered_map< uIntType, floatType > *volumeResult = reader.getMicroVolumes( );

    for ( auto it = volumeAnswer.begin( ); it != volumeAnswer.end( ); it++ ){

        auto r = volumeResult->find( it->first );

        if ( r == volumeResult->end( ) ){

            results << "test_initializeIncrement (test 14) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement (test 15) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microNodeReferencePositionAnswer
        =
            {
                { 15, { 0.00, 0.00, 1.00 } },
                { 31, { 0.00, 0.00, 1.50 } },
                { 13, { 0.50, 0.00, 1.00 } },
                { 26, { 0.50, 0.00, 1.50 } },
                { 53, { 0.00, 0.50, 1.00 } },
                { 21, { 0.00, 0.50, 1.50 } },
                { 37, { 0.50, 0.50, 1.00 } },
                { 48, { 0.50, 0.50, 1.50 } },
                {  5, { 1.00, 0.00, 1.00 } },
                { 10, { 1.00, 0.00, 1.50 } },
                {  3, { 1.00, 0.50, 1.00 } },
                {  4, { 1.00, 0.50, 1.50 } },
                { 32, { 0.50, 1.00, 1.00 } },
                { 33, { 0.50, 1.00, 1.50 } },
                { 34, { 1.00, 1.00, 1.00 } },
                { 28, { 1.00, 1.00, 1.50 } },
                { 25, { 0.00, 1.00, 1.00 } },
                { 50, { 0.00, 1.00, 1.50 } },
                { 43, { 0.00, 0.00, 2.00 } },
                { 27, { 0.50, 0.00, 2.00 } },
                {  1, { 0.00, 0.50, 2.00 } },
                {  7, { 0.50, 0.50, 2.00 } },
                { 30, { 1.00, 0.00, 2.00 } },
                { 16, { 1.00, 0.50, 2.00 } },
                { 22, { 0.50, 1.00, 2.00 } },
                {  2, { 1.00, 1.00, 2.00 } },
                { 46, { 0.00, 1.00, 2.00 } },
                { 24, { 0.00, 0.00, 0.00 } },
                { 39, { 0.00, 0.00, 0.50 } },
                { 40, { 0.50, 0.00, 0.00 } },
                { 57, { 0.50, 0.00, 0.50 } },
                { 44, { 0.00, 0.50, 0.00 } },
                { 58, { 0.00, 0.50, 0.50 } },
                { 29, { 0.50, 0.50, 0.00 } },
                { 59, { 0.50, 0.50, 0.50 } },
                { 11, { 1.00, 0.00, 0.00 } },
                {  0, { 1.00, 0.00, 0.50 } },
                { 20, { 1.00, 0.50, 0.00 } },
                { 60, { 1.00, 0.50, 0.50 } },
                { 47, { 0.50, 1.00, 0.00 } },
                { 49, { 0.50, 1.00, 0.50 } },
                { 17, { 1.00, 1.00, 0.00 } },
                { 38, { 1.00, 1.00, 0.50 } },
                { 14, { 0.00, 1.00, 0.00 } },
                { 55, { 0.00, 1.00, 0.50 } },
            };

    const std::unordered_map< uIntType, floatVector > *microNodeReferencePositionResult = reader.getMicroNodeReferencePositions( );

    for ( auto it = microNodeReferencePositionAnswer.begin( ); it != microNodeReferencePositionAnswer.end( ); it++ ){

        auto r = microNodeReferencePositionResult->find( it->first );

        if ( r == microNodeReferencePositionResult->end( ) ){

            results << "test_initializeIncrement (test 16) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 17) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroNodeReferencePositionsAnswer
        =
            {
                {  5, { 0.000, 0.000, 0.000 } },
                {  9, { 1.000, 0.000, 0.000 } },
                {  8, { 1.000, 1.000, 0.000 } },
                { 11, { 0.000, 1.000, 0.000 } },
                {  3, { 0.000, 0.000, 1.000 } },
                {  1, { 1.000, 0.000, 1.000 } },
                {  6, { 1.000, 1.000, 1.000 } },
                { 15, { 0.000, 1.000, 1.000 } },
                { 12, { 0.000, 0.000, 2.000 } },
                {  2, { 1.000, 0.000, 2.000 } },
                { 13, { 1.000, 1.000, 2.000 } },
                { 14, { 0.000, 1.000, 2.000 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroNodeReferencePositionsResult = reader.getMacroNodeReferencePositions( );

    for ( auto it = macroNodeReferencePositionsAnswer.begin( ); it != macroNodeReferencePositionsAnswer.end( ); it++ ){

        auto r = macroNodeReferencePositionsResult->find( it->first );

        if ( r == macroNodeReferencePositionsResult->end( ) ){

            results << "test_initializeIncrement (test 18) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 19) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, uIntVector > macroNodeReferenceConnectivityAnswer
        =
        {
            {  1, {  9,  5,  9,  8, 11,  3,  1,  6, 15 } },
            {  2, {  9,  3,  1,  6, 15, 12,  2, 13, 14 } },
        };

    const std::unordered_map< uIntType, uIntVector > *macroNodeReferenceConnectivityResult = reader.getMacroNodeReferenceConnectivity( );

    for ( auto it = macroNodeReferenceConnectivityAnswer.begin( ); it != macroNodeReferenceConnectivityAnswer.end( ); it++ ){

        auto r = macroNodeReferenceConnectivityResult->find( it->first );

        if ( r == macroNodeReferenceConnectivityResult->end( ) ){

            results << "test_initializeIncrement (test 21) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 22) & False\n";
            return 1;

        }

    }

    std::unordered_map< uIntType, floatVector > microDisplacementAnswer
        =
            {
               { 15, { 0.000, -14.400, -30.400 } },
               { 31, { 0.000, -14.400, -32.000 } },
               { 13, { 2.000, -14.400, -30.400 } },
               { 26, { 2.000, -14.400, -32.000 } },
               { 53, { 0.000, -11.200, -30.400 } },
               { 21, { 0.000, -11.200, -32.000 } },
               { 37, { 2.000, -11.200, -30.400 } },
               { 48, { 2.000, -11.200, -32.000 } },
               {  5, { 4.000, -14.400, -30.400 } },
               { 10, { 4.000, -14.400, -32.000 } },
               {  3, { 4.000, -11.200, -30.400 } },
               {  4, { 4.000, -11.200, -32.000 } },
               { 32, { 2.000, -8.000, -30.400 } },
               { 33, { 2.000, -8.000, -32.000 } },
               { 34, { 4.000, -8.000, -30.400 } },
               { 28, { 4.000, -8.000, -32.000 } },
               { 25, { 0.000, -8.000, -30.400 } },
               { 50, { 0.000, -8.000, -32.000 } },
               { 43, { 0.000, -14.400, -33.600 } },
               { 27, { 2.000, -14.400, -33.600 } },
               {  1, { 0.000, -11.200, -33.600 } },
               {  7, { 2.000, -11.200, -33.600 } },
               { 30, { 4.000, -14.400, -33.600 } },
               { 16, { 4.000, -11.200, -33.600 } },
               { 22, { 2.000, -8.000, -33.600 } },
               {  2, { 4.000, -8.000, -33.600 } },
               { 46, { 0.000, -8.000, -33.600 } },
               { 24, { 0.000, -14.400, -27.200 } },
               { 39, { 0.000, -14.400, -28.800 } },
               { 40, { 2.000, -14.400, -27.200 } },
               { 57, { 2.000, -14.400, -28.800 } },
               { 44, { 0.000, -11.200, -27.200 } },
               { 58, { 0.000, -11.200, -28.800 } },
               { 29, { 2.000, -11.200, -27.200 } },
               { 59, { 2.000, -11.200, -28.800 } },
               { 11, { 4.000, -14.400, -27.200 } },
               {  0, { 4.000, -14.400, -28.800 } },
               { 20, { 4.000, -11.200, -27.200 } },
               { 60, { 4.000, -11.200, -28.800 } },
               { 47, { 2.000, -8.000, -27.200 } },
               { 49, { 2.000, -8.000, -28.800 } },
               { 17, { 4.000, -8.000, -27.200 } },
               { 38, { 4.000, -8.000, -28.800 } },
               { 14, { 0.000, -8.000, -27.200 } },
               { 55, { 0.000, -8.000, -28.800 } },
           };


    const std::unordered_map< uIntType, floatVector > *microDisplacementResult = reader.getMicroDisplacements( );

    for ( auto it = microDisplacementAnswer.begin( ); it != microDisplacementAnswer.end( ); it++ ){

        auto r = microDisplacementResult->find( it->first );

        if ( r == microDisplacementResult->end( ) ){

            results << "test_initializeIncrement (test 22) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 23) & False\n";
            return 1;

        }

    }


    const uIntVector freeMacroCellIdsAnswer = { 1 };
    const uIntVector ghostMacroCellIdsAnswer = { 2 };

    const uIntVector *freeMacroCellIdsResult = reader.getFreeMacroCellIds( );
    const uIntVector *ghostMacroCellIdsResult = reader.getGhostMacroCellIds( );

    if ( !vectorTools::fuzzyEquals( freeMacroCellIdsAnswer, *freeMacroCellIdsResult ) ){
        results << "test_initializeIncrement (test 24) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( ghostMacroCellIdsAnswer, *ghostMacroCellIdsResult ) ){
        results << "test_initializeIncrement (test 25) & False\n";
        return 1;
    }

//    const uIntVector freeMacroCellMicroDomainCountsAnswer = { 8 };
//    const uIntVector ghostMacroCellMicroDomainCountsAnswer = { 8 };
//
//    const uIntVector *freeMacroCellMicroDomainCountsResult = reader.getFreeMacroCellMicroDomainCounts( );
//    const uIntVector *ghostMacroCellMicroDomainCountsResult = reader.getGhostMacroCellMicroDomainCounts( );
//
//    if ( !vectorTools::fuzzyEquals( freeMacroCellMicroDomainCountsAnswer, *freeMacroCellMicroDomainCountsResult ) ){
//        results << "test_initializeIncrement (test 26) & False\n";
//        return 1;
//    }
//
//    if ( !vectorTools::fuzzyEquals( ghostMacroCellMicroDomainCountsAnswer, *ghostMacroCellMicroDomainCountsResult ) ){
//        results << "test_initializeIncrement (test 27) & False\n";
//        return 1;
//    }

    const std::unordered_map< uIntType, floatVector > microBodyForcesAnswer
        =
            {
                { 15, { 2.000, 2.000, 9.200 } },
                { 31, { 2.000, 2.000, 12.800 } },
                { 13, { 5.000, 2.000, 9.200 } },
                { 26, { 5.000, 2.000, 12.800 } },
                { 53, { 2.000, 2.410, 9.200 } },
                { 21, { 2.000, 2.410, 12.800 } },
                { 37, { 5.000, 2.410, 9.200 } },
                { 48, { 5.000, 2.410, 12.800 } },
                {  5, { 8.000, 2.000, 9.200 } },
                { 10, { 8.000, 2.000, 12.800 } },
                {  3, { 8.000, 2.410, 9.200 } },
                {  4, { 8.000, 2.410, 12.800 } },
                { 32, { 5.000, 2.820, 9.200 } },
                { 33, { 5.000, 2.820, 12.800 } },
                { 34, { 8.000, 2.820, 9.200 } },
                { 28, { 8.000, 2.820, 12.800 } },
                { 25, { 2.000, 2.820, 9.200 } },
                { 50, { 2.000, 2.820, 12.800 } },
                { 43, { 2.000, 2.000, 16.400 } },
                { 27, { 5.000, 2.000, 16.400 } },
                {  1, { 2.000, 2.410, 16.400 } },
                {  7, { 5.000, 2.410, 16.400 } },
                { 30, { 8.000, 2.000, 16.400 } },
                { 16, { 8.000, 2.410, 16.400 } },
                { 22, { 5.000, 2.820, 16.400 } },
                {  2, { 8.000, 2.820, 16.400 } },
                { 46, { 2.000, 2.820, 16.400 } },
                { 24, { 2.000, 2.000, 2.000 } },
                { 39, { 2.000, 2.000, 5.600 } },
                { 40, { 5.000, 2.000, 2.000 } },
                { 57, { 5.000, 2.000, 5.600 } },
                { 44, { 2.000, 2.410, 2.000 } },
                { 58, { 2.000, 2.410, 5.600 } },
                { 29, { 5.000, 2.410, 2.000 } },
                { 59, { 5.000, 2.410, 5.600 } },
                { 11, { 8.000, 2.000, 2.000 } },
                {  0, { 8.000, 2.000, 5.600 } },
                { 20, { 8.000, 2.410, 2.000 } },
                { 60, { 8.000, 2.410, 5.600 } },
                { 47, { 5.000, 2.820, 2.000 } },
                { 49, { 5.000, 2.820, 5.600 } },
                { 17, { 8.000, 2.820, 2.000 } },
                { 38, { 8.000, 2.820, 5.600 } },
                { 14, { 2.000, 2.820, 2.000 } },
                { 55, { 2.000, 2.820, 5.600 } },
            };

    const std::unordered_map< uIntType, floatVector > *microBodyForcesResult = reader.getMicroBodyForces( );

    for ( auto it = microBodyForcesAnswer.begin( ); it != microBodyForcesAnswer.end( ); it++ ){

        auto r = microBodyForcesResult->find( it->first );

        if ( r == microBodyForcesResult->end( ) ){

            results << "test_initializeIncrement (test 28) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 29) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microSurfaceForcesAnswer
        =
            {
                { 15, { 2.100, 2.100, 3.540 } },
                { 31, { 2.100, 2.100, 4.260 } },
                { 13, { 2.460, 2.100, 3.540 } },
                { 26, { 2.460, 2.100, 4.260 } },
                { 53, { 2.100, 3.705, 3.540 } },
                { 21, { 2.100, 3.705, 4.260 } },
                { 37, { 2.460, 3.705, 3.540 } },
                { 48, { 2.460, 3.705, 4.260 } },
                {  5, { 2.820, 2.100, 3.540 } },
                { 10, { 2.820, 2.100, 4.260 } },
                {  3, { 2.820, 3.705, 3.540 } },
                {  4, { 2.820, 3.705, 4.260 } },
                { 32, { 2.460, 5.310, 3.540 } },
                { 33, { 2.460, 5.310, 4.260 } },
                { 34, { 2.820, 5.310, 3.540 } },
                { 28, { 2.820, 5.310, 4.260 } },
                { 25, { 2.100, 5.310, 3.540 } },
                { 50, { 2.100, 5.310, 4.260 } },
                { 43, { 2.100, 2.100, 4.980 } },
                { 27, { 2.460, 2.100, 4.980 } },
                {  1, { 2.100, 3.705, 4.980 } },
                {  7, { 2.460, 3.705, 4.980 } },
                { 30, { 2.820, 2.100, 4.980 } },
                { 16, { 2.820, 3.705, 4.980 } },
                { 22, { 2.460, 5.310, 4.980 } },
                {  2, { 2.820, 5.310, 4.980 } },
                { 46, { 2.100, 5.310, 4.980 } },
                { 24, { 2.100, 2.100, 2.100 } },
                { 39, { 2.100, 2.100, 2.820 } },
                { 40, { 2.460, 2.100, 2.100 } },
                { 57, { 2.460, 2.100, 2.820 } },
                { 44, { 2.100, 3.705, 2.100 } },
                { 58, { 2.100, 3.705, 2.820 } },
                { 29, { 2.460, 3.705, 2.100 } },
                { 59, { 2.460, 3.705, 2.820 } },
                { 11, { 2.820, 2.100, 2.100 } },
                {  0, { 2.820, 2.100, 2.820 } },
                { 20, { 2.820, 3.705, 2.100 } },
                { 60, { 2.820, 3.705, 2.820 } },
                { 47, { 2.460, 5.310, 2.100 } },
                { 49, { 2.460, 5.310, 2.820 } },
                { 17, { 2.820, 5.310, 2.100 } },
                { 38, { 2.820, 5.310, 2.820 } },
                { 14, { 2.100, 5.310, 2.100 } },
                { 55, { 2.100, 5.310, 2.820 } },
            };

    const std::unordered_map< uIntType, floatVector > *microSurfaceForcesResult = reader.getMicroSurfaceForces( );

    for ( auto it = microSurfaceForcesAnswer.begin( ); it != microSurfaceForcesAnswer.end( ); it++ ){

        auto r = microSurfaceForcesResult->find( it->first );

        if ( r == microSurfaceForcesResult->end( ) ){

            results << "test_initializeIncrement (test 30) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 31) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microExternalForcesAnswer
        =
            {
                { 15, { 4.100, 4.100, 12.740 } },
                { 31, { 4.100, 4.100, 17.060 } },
                { 13, { 7.460, 4.100, 12.740 } },
                { 26, { 7.460, 4.100, 17.060 } },
                { 53, { 4.100, 6.115, 12.740 } },
                { 21, { 4.100, 6.115, 17.060 } },
                { 37, { 7.460, 6.115, 12.740 } },
                { 48, { 7.460, 6.115, 17.060 } },
                {  5, { 10.820, 4.100, 12.740 } },
                { 10, { 10.820, 4.100, 17.060 } },
                {  3, { 10.820, 6.115, 12.740 } },
                {  4, { 10.820, 6.115, 17.060 } },
                { 32, { 7.460, 8.130, 12.740 } },
                { 33, { 7.460, 8.130, 17.060 } },
                { 34, { 10.820, 8.130, 12.740 } },
                { 28, { 10.820, 8.130, 17.060 } },
                { 25, { 4.100, 8.130, 12.740 } },
                { 50, { 4.100, 8.130, 17.060 } },
                { 43, { 4.100, 4.100, 21.380 } },
                { 27, { 7.460, 4.100, 21.380 } },
                {  1, { 4.100, 6.115, 21.380 } },
                {  7, { 7.460, 6.115, 21.380 } },
                { 30, { 10.820, 4.100, 21.380 } },
                { 16, { 10.820, 6.115, 21.380 } },
                { 22, { 7.460, 8.130, 21.380 } },
                {  2, { 10.820, 8.130, 21.380 } },
                { 46, { 4.100, 8.130, 21.380 } },
                { 24, { 4.100, 4.100, 4.100 } },
                { 39, { 4.100, 4.100, 8.420 } },
                { 40, { 7.460, 4.100, 4.100 } },
                { 57, { 7.460, 4.100, 8.420 } },
                { 44, { 4.100, 6.115, 4.100 } },
                { 58, { 4.100, 6.115, 8.420 } },
                { 29, { 7.460, 6.115, 4.100 } },
                { 59, { 7.460, 6.115, 8.420 } },
                { 11, { 10.820, 4.100, 4.100 } },
                {  0, { 10.820, 4.100, 8.420 } },
                { 20, { 10.820, 6.115, 4.100 } },
                { 60, { 10.820, 6.115, 8.420 } },
                { 47, { 7.460, 8.130, 4.100 } },
                { 49, { 7.460, 8.130, 8.420 } },
                { 17, { 10.820, 8.130, 4.100 } },
                { 38, { 10.820, 8.130, 8.420 } },
                { 14, { 4.100, 8.130, 4.100 } },
                { 55, { 4.100, 8.130, 8.420 } },
            };

    const std::unordered_map< uIntType, floatVector > *microExternalForcesResult = reader.getMicroExternalForces( );

    for ( auto it = microExternalForcesAnswer.begin( ); it != microExternalForcesAnswer.end( ); it++ ){

        auto r = microExternalForcesResult->find( it->first );

        if ( r == microExternalForcesResult->end( ) ){

            results << "test_initializeIncrement (test 32) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 33) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microVelocitiesAnswer
        =
            {
               { 15, { 5.000, -8.848, 11.950 } },
               { 31, { 5.000, -8.848, 13.025 } },
               { 13, { 6.250, -8.848, 11.950 } },
               { 26, { 6.250, -8.848, 13.025 } },
               { 53, { 5.000, -10.418, 11.950 } },
               { 21, { 5.000, -10.418, 13.025 } },
               { 37, { 6.250, -10.418, 11.950 } },
               { 48, { 6.250, -10.418, 13.025 } },
               {  5, { 7.500, -8.848, 11.950 } },
               { 10, { 7.500, -8.848, 13.025 } },
               {  3, { 7.500, -10.418, 11.950 } },
               {  4, { 7.500, -10.418, 13.025 } },
               { 32, { 6.250, -11.988, 11.950 } },
               { 33, { 6.250, -11.988, 13.025 } },
               { 34, { 7.500, -11.988, 11.950 } },
               { 28, { 7.500, -11.988, 13.025 } },
               { 25, { 5.000, -11.988, 11.950 } },
               { 50, { 5.000, -11.988, 13.025 } },
               { 43, { 5.000, -8.848, 14.100 } },
               { 27, { 6.250, -8.848, 14.100 } },
               {  1, { 5.000, -10.418, 14.100 } },
               {  7, { 6.250, -10.418, 14.100 } },
               { 30, { 7.500, -8.848, 14.100 } },
               { 16, { 7.500, -10.418, 14.100 } },
               { 22, { 6.250, -11.988, 14.100 } },
               {  2, { 7.500, -11.988, 14.100 } },
               { 46, { 5.000, -11.988, 14.100 } },
               { 24, { 5.000, -8.848, 9.800 } },
               { 39, { 5.000, -8.848, 10.875 } },
               { 40, { 6.250, -8.848, 9.800 } },
               { 57, { 6.250, -8.848, 10.875 } },
               { 44, { 5.000, -10.418, 9.800 } },
               { 58, { 5.000, -10.418, 10.875 } },
               { 29, { 6.250, -10.418, 9.800 } },
               { 59, { 6.250, -10.418, 10.875 } },
               { 11, { 7.500, -8.848, 9.800 } },
               {  0, { 7.500, -8.848, 10.875 } },
               { 20, { 7.500, -10.418, 9.800 } },
               { 60, { 7.500, -10.418, 10.875 } },
               { 47, { 6.250, -11.988, 9.800 } },
               { 49, { 6.250, -11.988, 10.875 } },
               { 17, { 7.500, -11.988, 9.800 } },
               { 38, { 7.500, -11.988, 10.875 } },
               { 14, { 5.000, -11.988, 9.800 } },
               { 55, { 5.000, -11.988, 10.875 } }
            };

    const std::unordered_map< uIntType, floatVector > *microVelocitiesResult = reader.getMicroVelocities( );

    for ( auto it = microVelocitiesAnswer.begin( ); it != microVelocitiesAnswer.end( ); it++ ){

        auto r = microVelocitiesResult->find( it->first );

        if ( r == microVelocitiesResult->end( ) ){

            results << "test_initializeIncrement (test 34) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 35) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microAccelerationsAnswer
        =
            {
               { 15, { 5.7765, 5.9930, 1.1000 } },
               { 31, { 5.7765, 5.9930, 2.1500 } },
               { 13, { 6.1315, 5.9930, 1.1000 } },
               { 26, { 6.1315, 5.9930, 2.1500 } },
               { 53, { 5.7765, 6.4080, 1.1000 } },
               { 21, { 5.7765, 6.4080, 2.1500 } },
               { 37, { 6.1315, 6.4080, 1.1000 } },
               { 48, { 6.1315, 6.4080, 2.1500 } },
               {  5, { 6.4865, 5.9930, 1.1000 } },
               { 10, { 6.4865, 5.9930, 2.1500 } },
               {  3, { 6.4865, 6.4080, 1.1000 } },
               {  4, { 6.4865, 6.4080, 2.1500 } },
               { 32, { 6.1315, 6.8230, 1.1000 } },
               { 33, { 6.1315, 6.8230, 2.1500 } },
               { 34, { 6.4865, 6.8230, 1.1000 } },
               { 28, { 6.4865, 6.8230, 2.1500 } },
               { 25, { 5.7765, 6.8230, 1.1000 } },
               { 50, { 5.7765, 6.8230, 2.1500 } },
               { 43, { 5.7765, 5.9930, 3.2000 } },
               { 27, { 6.1315, 5.9930, 3.2000 } },
               {  1, { 5.7765, 6.4080, 3.2000 } },
               {  7, { 6.1315, 6.4080, 3.2000 } },
               { 30, { 6.4865, 5.9930, 3.2000 } },
               { 16, { 6.4865, 6.4080, 3.2000 } },
               { 22, { 6.1315, 6.8230, 3.2000 } },
               {  2, { 6.4865, 6.8230, 3.2000 } },
               { 46, { 5.7765, 6.8230, 3.2000 } },
               { 24, { 5.7765, 5.9930, -1.0000 } },
               { 39, { 5.7765, 5.9930, 0.0500 } },
               { 40, { 6.1315, 5.9930, -1.0000 } },
               { 57, { 6.1315, 5.9930, 0.0500 } },
               { 44, { 5.7765, 6.4080, -1.0000 } },
               { 58, { 5.7765, 6.4080, 0.0500 } },
               { 29, { 6.1315, 6.4080, -1.0000 } },
               { 59, { 6.1315, 6.4080, 0.0500 } },
               { 11, { 6.4865, 5.9930, -1.0000 } },
               {  0, { 6.4865, 5.9930, 0.0500 } },
               { 20, { 6.4865, 6.4080, -1.0000 } },
               { 60, { 6.4865, 6.4080, 0.0500 } },
               { 47, { 6.1315, 6.8230, -1.0000 } },
               { 49, { 6.1315, 6.8230, 0.0500 } },
               { 17, { 6.4865, 6.8230, -1.0000 } },
               { 38, { 6.4865, 6.8230, 0.0500 } },
               { 14, { 5.7765, 6.8230, -1.0000 } },
               { 55, { 5.7765, 6.8230, 0.0500 } },
           };

    const std::unordered_map< uIntType, floatVector > *microAccelerationsResult = reader.getMicroAccelerations( );

    for ( auto it = microAccelerationsAnswer.begin( ); it != microAccelerationsAnswer.end( ); it++ ){

        auto r = microAccelerationsResult->find( it->first );

        if ( r == microAccelerationsResult->end( ) ){

            results << "test_initializeIncrement (test 36) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 37) & False\n";
            return 1;

        }

    }

    const floatType previousTimeAnswer = 0.;
    const floatType *previousTimeResult = reader.getPreviousMicroTime( );

    if ( !vectorTools::fuzzyEquals( previousTimeAnswer, *previousTimeResult ) ){

        results << "test_initializeIncrement (test 38) & False\n";
        return 1;

    }

    const std::unordered_map< uIntType, floatVector > previousMicroDisplacementAnswer
        =
            {
                { 15, { -8.000, -8.000, -11.200 } },
                { 31, { -8.000, -8.000, -12.800 } },
                { 13, { -6.000, -8.000, -11.200 } },
                { 26, { -6.000, -8.000, -12.800 } },
                { 53, { -8.000, -4.800, -11.200 } },
                { 21, { -8.000, -4.800, -12.800 } },
                { 37, { -6.000, -4.800, -11.200 } },
                { 48, { -6.000, -4.800, -12.800 } },
                {  5, { -4.000, -8.000, -11.200 } },
                { 10, { -4.000, -8.000, -12.800 } },
                {  3, { -4.000, -4.800, -11.200 } },
                {  4, { -4.000, -4.800, -12.800 } },
                { 32, { -6.000, -1.600, -11.200 } },
                { 33, { -6.000, -1.600, -12.800 } },
                { 34, { -4.000, -1.600, -11.200 } },
                { 28, { -4.000, -1.600, -12.800 } },
                { 25, { -8.000, -1.600, -11.200 } },
                { 50, { -8.000, -1.600, -12.800 } },
                { 43, { -8.000, -8.000, -14.400 } },
                { 27, { -6.000, -8.000, -14.400 } },
                {  1, { -8.000, -4.800, -14.400 } },
                {  7, { -6.000, -4.800, -14.400 } },
                { 30, { -4.000, -8.000, -14.400 } },
                { 16, { -4.000, -4.800, -14.400 } },
                { 22, { -6.000, -1.600, -14.400 } },
                {  2, { -4.000, -1.600, -14.400 } },
                { 46, { -8.000, -1.600, -14.400 } },
                { 24, { -8.000, -8.000, -8.000 } },
                { 39, { -8.000, -8.000, -9.600 } },
                { 40, { -6.000, -8.000, -8.000 } },
                { 57, { -6.000, -8.000, -9.600 } },
                { 44, { -8.000, -4.800, -8.000 } },
                { 58, { -8.000, -4.800, -9.600 } },
                { 29, { -6.000, -4.800, -8.000 } },
                { 59, { -6.000, -4.800, -9.600 } },
                { 11, { -4.000, -8.000, -8.000 } },
                {  0, { -4.000, -8.000, -9.600 } },
                { 20, { -4.000, -4.800, -8.000 } },
                { 60, { -4.000, -4.800, -9.600 } },
                { 47, { -6.000, -1.600, -8.000 } },
                { 49, { -6.000, -1.600, -9.600 } },
                { 17, { -4.000, -1.600, -8.000 } },
                { 38, { -4.000, -1.600, -9.600 } },
                { 14, { -8.000, -1.600, -8.000 } },
                { 55, { -8.000, -1.600, -9.600 } },
            };

    const std::unordered_map< uIntType, floatVector > *previousMicroDisplacementResult = reader.getPreviousMicroDisplacements( );

    for ( auto it = previousMicroDisplacementAnswer.begin( ); it != previousMicroDisplacementAnswer.end( ); it++ ){

        auto r = previousMicroDisplacementResult->find( it->first );

        if ( r == previousMicroDisplacementResult->end( ) ){

            results << "test_initializeIncrement (test 39) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 40) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > previousMicroVelocitiesAnswer
        =
            {
               { 15, { 1.200, 1.200, 3.350 } },
               { 31, { 1.200, 1.200, 4.425 } },
               { 13, { 2.450, 1.200, 3.350 } },
               { 26, { 2.450, 1.200, 4.425 } },
               { 53, { 1.200, -0.370, 3.350 } },
               { 21, { 1.200, -0.370, 4.425 } },
               { 37, { 2.450, -0.370, 3.350 } },
               { 48, { 2.450, -0.370, 4.425 } },
               {  5, { 3.700, 1.200, 3.350 } },
               { 10, { 3.700, 1.200, 4.425 } },
               {  3, { 3.700, -0.370, 3.350 } },
               {  4, { 3.700, -0.370, 4.425 } },
               { 32, { 2.450, -1.940, 3.350 } },
               { 33, { 2.450, -1.940, 4.425 } },
               { 34, { 3.700, -1.940, 3.350 } },
               { 28, { 3.700, -1.940, 4.425 } },
               { 25, { 1.200, -1.940, 3.350 } },
               { 50, { 1.200, -1.940, 4.425 } },
               { 43, { 1.200, 1.200, 5.500 } },
               { 27, { 2.450, 1.200, 5.500 } },
               {  1, { 1.200, -0.370, 5.500 } },
               {  7, { 2.450, -0.370, 5.500 } },
               { 30, { 3.700, 1.200, 5.500 } },
               { 16, { 3.700, -0.370, 5.500 } },
               { 22, { 2.450, -1.940, 5.500 } },
               {  2, { 3.700, -1.940, 5.500 } },
               { 46, { 1.200, -1.940, 5.500 } },
               { 24, { 1.200, 1.200, 1.200 } },
               { 39, { 1.200, 1.200, 2.275 } },
               { 40, { 2.450, 1.200, 1.200 } },
               { 57, { 2.450, 1.200, 2.275 } },
               { 44, { 1.200, -0.370, 1.200 } },
               { 58, { 1.200, -0.370, 2.275 } },
               { 29, { 2.450, -0.370, 1.200 } },
               { 59, { 2.450, -0.370, 2.275 } },
               { 11, { 3.700, 1.200, 1.200 } },
               {  0, { 3.700, 1.200, 2.275 } },
               { 20, { 3.700, -0.370, 1.200 } },
               { 60, { 3.700, -0.370, 2.275 } },
               { 47, { 2.450, -1.940, 1.200 } },
               { 49, { 2.450, -1.940, 2.275 } },
               { 17, { 3.700, -1.940, 1.200 } },
               { 38, { 3.700, -1.940, 2.275 } },
               { 14, { 1.200, -1.940, 1.200 } },
               { 55, { 1.200, -1.940, 2.275 } },
           };

    const std::unordered_map< uIntType, floatVector > *previousMicroVelocitiesResult = reader.getPreviousMicroVelocities( );

    for ( auto it = previousMicroVelocitiesAnswer.begin( ); it != previousMicroVelocitiesAnswer.end( ); it++ ){

        auto r = previousMicroVelocitiesResult->find( it->first );

        if ( r == previousMicroVelocitiesResult->end( ) ){

            results << "test_initializeIncrement (test 41) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 42) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > previousMicroAccelerationsAnswer
        =
            {
                { 15, { 4.2500, 4.2500, 6.3500 } },
                { 31, { 4.2500, 4.2500, 7.4000 } },
                { 13, { 4.6050, 4.2500, 6.3500 } },
                { 26, { 4.6050, 4.2500, 7.4000 } },
                { 53, { 4.2500, 4.6650, 6.3500 } },
                { 21, { 4.2500, 4.6650, 7.4000 } },
                { 37, { 4.6050, 4.6650, 6.3500 } },
                { 48, { 4.6050, 4.6650, 7.4000 } },
                {  5, { 4.9600, 4.2500, 6.3500 } },
                { 10, { 4.9600, 4.2500, 7.4000 } },
                {  3, { 4.9600, 4.6650, 6.3500 } },
                {  4, { 4.9600, 4.6650, 7.4000 } },
                { 32, { 4.6050, 5.0800, 6.3500 } },
                { 33, { 4.6050, 5.0800, 7.4000 } },
                { 34, { 4.9600, 5.0800, 6.3500 } },
                { 28, { 4.9600, 5.0800, 7.4000 } },
                { 25, { 4.2500, 5.0800, 6.3500 } },
                { 50, { 4.2500, 5.0800, 7.4000 } },
                { 43, { 4.2500, 4.2500, 8.4500 } },
                { 27, { 4.6050, 4.2500, 8.4500 } },
                {  1, { 4.2500, 4.6650, 8.4500 } },
                {  7, { 4.6050, 4.6650, 8.4500 } },
                { 30, { 4.9600, 4.2500, 8.4500 } },
                { 16, { 4.9600, 4.6650, 8.4500 } },
                { 22, { 4.6050, 5.0800, 8.4500 } },
                {  2, { 4.9600, 5.0800, 8.4500 } },
                { 46, { 4.2500, 5.0800, 8.4500 } },
                { 24, { 4.2500, 4.2500, 4.2500 } },
                { 39, { 4.2500, 4.2500, 5.3000 } },
                { 40, { 4.6050, 4.2500, 4.2500 } },
                { 57, { 4.6050, 4.2500, 5.3000 } },
                { 44, { 4.2500, 4.6650, 4.2500 } },
                { 58, { 4.2500, 4.6650, 5.3000 } },
                { 29, { 4.6050, 4.6650, 4.2500 } },
                { 59, { 4.6050, 4.6650, 5.3000 } },
                { 11, { 4.9600, 4.2500, 4.2500 } },
                {  0, { 4.9600, 4.2500, 5.3000 } },
                { 20, { 4.9600, 4.6650, 4.2500 } },
                { 60, { 4.9600, 4.6650, 5.3000 } },
                { 47, { 4.6050, 5.0800, 4.2500 } },
                { 49, { 4.6050, 5.0800, 5.3000 } },
                { 17, { 4.9600, 5.0800, 4.2500 } },
                { 38, { 4.9600, 5.0800, 5.3000 } },
                { 14, { 4.2500, 5.0800, 4.2500 } },
                { 55, { 4.2500, 5.0800, 5.3000 } },
            };

    const std::unordered_map< uIntType, floatVector > *previousMicroAccelerationsResult = reader.getPreviousMicroAccelerations( );

    for ( auto it = previousMicroAccelerationsAnswer.begin( ); it != previousMicroAccelerationsAnswer.end( ); it++ ){

        auto r = previousMicroAccelerationsResult->find( it->first );

        if ( r == previousMicroAccelerationsResult->end( ) ){

            results << "test_initializeIncrement (test 43) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 44) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microStressesAnswer
        =
            {
                { 15, { 3.090371, 3.765898, 2.979736, 0.679351, 4.038599, 2.932600, 2.403665, 3.732252, 0.398282 } },
                { 31, { 3.726374, 4.279899, 2.727951, 0.186259, 4.194622, 3.003618, 2.646383, 3.839766, -0.401559 } },
                { 13, { 3.606480, 3.522090, 2.946139, 1.672865, 4.876204, 2.621181, 3.337650, 4.342345, 0.722520 } },
                { 26, { 4.242483, 4.036091, 2.694354, 1.179773, 5.032227, 2.692198, 3.580368, 4.449859, -0.077322 } },
                { 53, { 2.629139, 3.951084, 3.104852, 0.011765, 3.878585, 3.875028, 3.382334, 3.573368, 0.176006 } },
                { 21, { 3.265143, 4.465084, 2.853068, -0.481326, 4.034608, 3.946046, 3.625052, 3.680882, -0.623836 } },
                { 37, { 3.145248, 3.707276, 3.071255, 1.005279, 4.716189, 3.563609, 4.316319, 4.183461, 0.500244 } },
                { 48, { 3.781252, 4.221276, 2.819471, 0.512188, 4.872212, 3.634626, 4.559037, 4.290975, -0.299598 } },
                {  5, { 4.122589, 3.278282, 2.912542, 2.666379, 5.713808, 2.309761, 4.271634, 4.952438, 1.046758 } },
                { 10, { 4.758592, 3.792283, 2.660757, 2.173287, 5.869831, 2.380779, 4.514352, 5.059952, 0.246916 } },
                {  3, { 3.661357, 3.463468, 3.037658, 1.998793, 5.553793, 3.252189, 5.250303, 4.793554, 0.824481 } },
                {  4, { 4.297360, 3.977468, 2.785874, 1.505701, 5.709816, 3.323207, 5.493021, 4.901068, 0.024640 } },
                { 32, { 2.684016, 3.892461, 3.196372, 0.337694, 4.556174, 4.506037, 5.294987, 4.024577, 0.277967 } },
                { 33, { 3.320020, 4.406462, 2.944587, -0.155398, 4.712197, 4.577054, 5.537706, 4.132091, -0.521874 } },
                { 34, { 3.200125, 3.648653, 3.162775, 1.331208, 5.393779, 4.194617, 6.228972, 4.634670, 0.602205 } },
                { 28, { 3.836129, 4.162654, 2.910990, 0.838116, 5.549802, 4.265635, 6.471690, 4.742184, -0.197637 } },
                { 25, { 2.167908, 4.136269, 3.229969, -0.655820, 3.718570, 4.817456, 4.361003, 3.414484, -0.046270 } },
                { 50, { 2.803911, 4.650270, 2.978184, -1.148912, 3.874593, 4.888474, 4.603721, 3.521998, -0.846112 } },
                { 43, { 4.362378, 4.793900, 2.476167, -0.306832, 4.350645, 3.074636, 2.889102, 3.947280, -1.201401 } },
                { 27, { 4.878487, 4.550092, 2.442570, 0.686682, 5.188250, 2.763216, 3.823086, 4.557373, -0.877164 } },
                {  1, { 3.901146, 4.979085, 2.601283, -0.974418, 4.190631, 4.017064, 3.867770, 3.788396, -1.423677 } },
                {  7, { 4.417255, 4.735277, 2.567686, 0.019096, 5.028235, 3.705644, 4.801755, 4.398490, -1.099440 } },
                { 30, { 5.394596, 4.306284, 2.408973, 1.680195, 6.025854, 2.451796, 4.757071, 5.167467, -0.552926 } },
                { 16, { 4.933364, 4.491469, 2.534089, 1.012610, 5.865839, 3.394224, 5.735740, 5.008583, -0.775202 } },
                { 22, { 3.956023, 4.920462, 2.692802, -0.648490, 4.868220, 4.648072, 5.780424, 4.239606, -1.321716 } },
                {  2, { 4.472132, 4.676654, 2.659205, 0.345024, 5.705825, 4.336652, 6.714408, 4.849699, -0.997479 } },
                { 46, { 3.439914, 5.164270, 2.726399, -1.642003, 4.030616, 4.959492, 4.846439, 3.629512, -1.645954 } },
                { 24, { 1.818364, 2.737897, 3.483305, 1.665534, 3.726553, 2.790565, 1.918229, 3.517223, 1.997966 } },
                { 39, { 2.454367, 3.251897, 3.231521, 1.172443, 3.882576, 2.861583, 2.160947, 3.624737, 1.198124 } },
                { 40, { 2.334473, 2.494089, 3.449708, 2.659048, 4.564158, 2.479145, 2.852213, 4.127316, 2.322204 } },
                { 57, { 2.970476, 3.008089, 3.197924, 2.165957, 4.720181, 2.550163, 3.094932, 4.234830, 1.522362 } },
                { 44, { 1.357132, 2.923082, 3.608422, 0.997949, 3.566539, 3.732993, 2.896898, 3.358339, 1.775690 } },
                { 58, { 1.993136, 3.437083, 3.356637, 0.504857, 3.722562, 3.804011, 3.139616, 3.465853, 0.975848 } },
                { 29, { 1.873241, 2.679274, 3.574825, 1.991463, 4.404143, 3.421573, 3.830882, 3.968432, 2.099927 } },
                { 59, { 2.509245, 3.193275, 3.323040, 1.498371, 4.560166, 3.492591, 4.073600, 4.075946, 1.300085 } },
                { 11, { 2.850582, 2.250281, 3.416111, 3.652562, 5.401762, 2.167726, 3.786198, 4.737409, 2.646441 } },
                {  0, { 3.486585, 2.764281, 3.164327, 3.159470, 5.557785, 2.238743, 4.028916, 4.844923, 1.846599 } },
                { 20, { 2.389350, 2.435466, 3.541228, 2.984976, 5.241747, 3.110154, 4.764867, 4.578525, 2.424165 } },
                { 60, { 3.025354, 2.949467, 3.289443, 2.491885, 5.397770, 3.181171, 5.007585, 4.686039, 1.624323 } },
                { 47, { 1.412010, 2.864460, 3.699941, 1.323877, 4.244128, 4.364001, 4.809551, 3.809548, 1.877651 } },
                { 49, { 2.048013, 3.378460, 3.448156, 0.830785, 4.400151, 4.435019, 5.052269, 3.917062, 1.077809 } },
                { 17, { 1.928119, 2.620652, 3.666344, 2.317391, 5.081733, 4.052582, 5.743536, 4.419641, 2.201888 } },
                { 38, { 2.564122, 3.134652, 3.414559, 1.824299, 5.237756, 4.123599, 5.986254, 4.527156, 1.402047 } },
                { 14, { 0.895901, 3.108268, 3.733538, 0.330363, 3.406524, 4.675421, 3.875566, 3.199455, 1.553413 } },
                { 55, { 1.531904, 3.622268, 3.481753, -0.162728, 3.562547, 4.746439, 4.118285, 3.306969, 0.753572 } },
            };

    const std::unordered_map< uIntType, floatVector > *microStressesResult = reader.getMicroStresses( );

    for ( auto it = microStressesAnswer.begin( ); it != microStressesAnswer.end( ); it++ ){

        auto r = microStressesResult->find( it->first );

        if ( r == microStressesResult->end( ) ){

            results << "test_initializeIncrement (test 45) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 46) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microInternalForcesAnswer
        =
             {
                { 15, { 2.851562, 4.231133, 1.866341 } },
                { 31, { 3.344066, 4.560423, 0.901855 } },
                { 13, { 3.793505, 3.271966, 2.243144 } },
                { 26, { 4.286010, 3.601256, 1.278659 } },
                { 53, { 2.646243, 5.060558, 1.343330 } },
                { 21, { 3.138748, 5.389848, 0.378845 } },
                { 37, { 3.588187, 4.101391, 1.720134 } },
                { 48, { 4.080691, 4.430681, 0.755648 } },
                {  5, { 4.735449, 2.312798, 2.619948 } },
                { 10, { 5.227953, 2.642088, 1.655463 } },
                {  3, { 4.530131, 3.142223, 2.096937 } },
                {  4, { 5.022635, 3.471513, 1.132452 } },
                { 32, { 3.382869, 4.930815, 1.197123 } },
                { 33, { 3.875373, 5.260105, 0.232638 } },
                { 34, { 4.324812, 3.971648, 1.573927 } },
                { 28, { 4.817317, 4.300938, 0.609442 } },
                { 25, { 2.440925, 5.889983, 0.820320 } },
                { 50, { 2.933429, 6.219273, -0.144166 } },
                { 43, { 3.836570, 4.889713, -0.062630 } },
                { 27, { 4.778514, 3.930546, 0.314174 } },
                {  1, { 3.631252, 5.719138, -0.585641 } },
                {  7, { 4.573196, 4.759971, -0.208837 } },
                { 30, { 5.720458, 2.971378, 0.690977 } },
                { 16, { 5.515139, 3.800803, 0.167967 } },
                { 22, { 4.367877, 5.589395, -0.731847 } },
                {  2, { 5.309821, 4.630228, -0.355044 } },
                { 46, { 3.425933, 6.548563, -1.108651 } },
                { 24, { 1.866553, 3.572553, 3.795311 } },
                { 39, { 2.359057, 3.901843, 2.830826 } },
                { 40, { 2.808497, 2.613386, 4.172115 } },
                { 57, { 3.301001, 2.942676, 3.207630 } },
                { 44, { 1.661235, 4.401978, 3.272301 } },
                { 58, { 2.153739, 4.731268, 2.307815 } },
                { 29, { 2.603179, 3.442811, 3.649104 } },
                { 59, { 3.095683, 3.772101, 2.684619 } },
                { 11, { 3.750441, 1.654218, 4.548919 } },
                {  0, { 4.242945, 1.983508, 3.584433 } },
                { 20, { 3.545122, 2.483643, 4.025908 } },
                { 60, { 4.037627, 2.812933, 3.061423 } },
                { 47, { 2.397860, 4.272235, 3.126094 } },
                { 49, { 2.890364, 4.601525, 2.161609 } },
                { 17, { 3.339804, 3.313068, 3.502898 } },
                { 38, { 3.832308, 3.642358, 2.538412 } },
                { 14, { 1.455917, 5.231403, 2.749290 } },
                { 55, { 1.948421, 5.560693, 1.784805 } },
            };

    const std::unordered_map< uIntType, floatVector > *microInternalForcesResult = reader.getMicroInternalForces( );

    for ( auto it = microInternalForcesAnswer.begin( ); it != microInternalForcesAnswer.end( ); it++ ){

        auto r = microInternalForcesResult->find( it->first );

        if ( r == microInternalForcesResult->end( ) ){

            results << "test_initializeIncrement (test 47) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 48) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microInertialForcesAnswer
        =
            {
                { 15, { 0.915926, -0.776311, -3.061289 } },
                { 31, { 1.713875, 0.043499, -3.050079 } },
                { 13, { 0.277727, -0.625485, -3.510893 } },
                { 26, { 1.075675, 0.194325, -3.499683 } },
                { 53, { 0.331368, -0.491830, -3.876881 } },
                { 21, { 1.129317, 0.327980, -3.865671 } },
                { 37, { -0.306832, -0.341004, -4.326486 } },
                { 48, { 0.491117, 0.478807, -4.315275 } },
                {  5, { -0.360473, -0.474659, -3.960498 } },
                { 10, { 0.437475, 0.345151, -3.949287 } },
                {  3, { -0.945031, -0.190177, -4.776090 } },
                {  4, { -0.147083, 0.629633, -4.764879 } },
                { 32, { -0.891390, -0.056522, -5.142078 } },
                { 33, { -0.093442, 0.763288, -5.130867 } },
                { 34, { -1.529590, 0.094304, -5.591682 } },
                { 28, { -0.731641, 0.914114, -5.580471 } },
                { 25, { -0.253190, -0.207348, -4.692474 } },
                { 50, { 0.544758, 0.612462, -4.681263 } },
                { 43, { 2.511823, 0.863309, -3.038868 } },
                { 27, { 1.873624, 1.014136, -3.488472 } },
                {  1, { 1.927265, 1.147791, -3.854460 } },
                {  7, { 1.289065, 1.298617, -4.304064 } },
                { 30, { 1.235424, 1.164962, -3.938076 } },
                { 16, { 0.650865, 1.449443, -4.753668 } },
                { 22, { 0.704507, 1.583098, -5.119656 } },
                {  2, { 0.066307, 1.733925, -5.569261 } },
                { 46, { 1.342707, 1.432272, -4.670052 } },
                { 24, { -0.679970, -2.415932, -3.083711 } },
                { 39, { 0.117978, -1.596122, -3.072500 } },
                { 40, { -1.318170, -2.265106, -3.533315 } },
                { 57, { -0.520222, -1.445295, -3.522104 } },
                { 44, { -1.264529, -2.131450, -3.899303 } },
                { 58, { -0.466580, -1.311640, -3.888092 } },
                { 29, { -1.902729, -1.980624, -4.348907 } },
                { 59, { -1.104780, -1.160814, -4.337696 } },
                { 11, { -1.956370, -2.114279, -3.982919 } },
                {  0, { -1.158422, -1.294469, -3.971708 } },
                { 20, { -2.540928, -1.829798, -4.798511 } },
                { 60, { -1.742980, -1.009988, -4.787300 } },
                { 47, { -2.487287, -1.696143, -5.164499 } },
                { 49, { -1.689339, -0.876332, -5.153289 } },
                { 17, { -3.125487, -1.545316, -5.614103 } },
                { 38, { -2.327538, -0.725506, -5.602893 } },
                { 14, { -1.849087, -1.846969, -4.714895 } },
                { 55, { -1.051139, -1.027159, -4.703684 } },
            };

    const std::unordered_map< uIntType, floatVector > *microInertialForcesResult = reader.getMicroInertialForces( );

    for ( auto it = microInertialForcesAnswer.begin( ); it != microInertialForcesAnswer.end( ); it++ ){

        auto r = microInertialForcesResult->find( it->first );

        if ( r == microInertialForcesResult->end( ) ){

            results << "test_initializeIncrement (test 49) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 50) & False\n";
            return 1;

        }

    }

    const floatType macroTimeAnswer = 1.;
    const floatType* macroTimeResult = reader.getMacroTime( );

    if ( !vectorTools::fuzzyEquals( macroTimeAnswer, *macroTimeResult ) ){

        results << "test_initializeIncrement (test 51) & False\n";
        return 1;

    }

    const std::unordered_map< uIntType, floatVector > macroDisplacementsAnswer
        =
            {
                {  5, { 0.641275, 0.232390, -2.327989 } },
                {  9, { -0.277488, 1.864821, -0.959118 } },
                {  8, { -1.872777, 2.331052, -2.562691 } },
                { 11, { -0.954015, 0.698621, -3.931561 } },
                {  3, { 0.863789, 1.140577, -2.616417 } },
                {  1, { -0.054974, 2.773008, -1.247547 } },
                {  6, { -1.650263, 3.239239, -2.851120 } },
                { 15, { -0.731501, 1.606808, -4.219990 } },
                { 12, { 1.086303, 2.048764, -2.904846 } },
                {  2, { 0.167540, 3.681195, -1.535975 } },
                { 13, { -1.427749, 4.147426, -3.139548 } },
                { 14, { -0.508987, 2.514995, -4.508419 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroDisplacementsResult = reader.getMacroDisplacements( );

    for ( auto it = macroDisplacementsAnswer.begin( ); it != macroDisplacementsAnswer.end( ); it++ ){

        auto r = macroDisplacementsResult->find( it->first );

        if ( r == macroDisplacementsResult->end( ) ){

            results << "test_initializeIncrement (test 51) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 52) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroDispDOFVectorAnswer
        =
            {
                {  5, { 0.641275, 0.232390, -2.327989, 2.476106, 3.649307, 0.380024, -0.602181, -0.098268, 1.214942, -1.346951, 0.875060, 0.810153 } },
                {  9, { -0.277488, 1.864821, -0.959118, 3.747254, 2.993086, 1.602849, -0.801713, 0.385446, 1.268479, -0.333959, 0.834465, 0.773507 } },
                {  8, { -1.872777, 2.331052, -2.562691, 5.125425, 4.185234, 0.239589, -2.210062, -1.562753, 0.311109, -2.017245, 2.360391, -1.174211 } },
                { 11, { -0.954015, 0.698621, -3.931561, 3.854278, 4.841455, -0.983235, -2.010531, -2.046468, 0.257571, -3.030237, 2.400986, -1.137565 } },
                {  3, { 0.863789, 1.140577, -2.616417, 3.550081, 5.494371, 0.536456, 0.929143, -0.090355, 2.448663, -3.175893, 0.545918, -0.235911 } },
                {  1, { -0.054974, 2.773008, -1.247547, 4.821228, 4.838150, 1.759280, 0.729611, 0.393360, 2.502200, -2.162901, 0.505323, -0.272557 } },
                {  6, { -1.650263, 3.239239, -2.851120, 6.199400, 6.030298, 0.396021, -0.678739, -1.554840, 1.544830, -3.846187, 2.031250, -2.220275 } },
                { 15, { -0.731501, 1.606808, -4.219990, 4.928252, 6.686519, -0.826804, -0.479207, -2.038554, 1.491293, -4.859179, 2.071844, -2.183628 } },
                { 12, { 1.086303, 2.048764, -2.904846, 4.624056, 7.339434, 0.692887, 2.460467, -0.082442, 3.682384, -5.004835, 0.216776, -1.281975 } },
                {  2, { 0.167540, 3.681195, -1.535975, 5.895203, 6.683213, 1.915712, 2.260935, 0.401273, 3.735921, -3.991843, 0.176182, -1.318621 } },
                { 13, { -1.427749, 4.147426, -3.139548, 7.273375, 7.875361, 0.552453, 0.852585, -1.546927, 2.778551, -5.675130, 1.702108, -3.266339 } },
                { 14, { -0.508987, 2.514995, -4.508419, 6.002227, 8.531582, -0.670372, 1.052117, -2.030641, 2.725014, -6.688121, 1.742702, -3.229692 } },
            };
    const std::unordered_map< uIntType, floatVector > *macroDispDOFVectorResult = reader.getMacroDispDOFVector( );

    for ( auto it = macroDispDOFVectorAnswer.begin( ); it != macroDispDOFVectorAnswer.end( ); it++ ){

        auto r = macroDispDOFVectorResult->find( it->first );

        if ( r == macroDispDOFVectorResult->end( ) ){

            results << "test_initializeIncrement (test 53) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 54) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroVelocitiesAnswer
        =
            {
                {  5, { -2.346964, -2.328088, 2.337123, 1.674508, 1.402881, -1.276984, 2.105710, 1.195199, 0.817334, 1.284606, -0.465939, -1.366498 } },
                {  9, { -0.548438, -2.321249, 3.137203, -0.038677, 2.051591, -1.215603, 0.532880, 0.749185, 1.584563, 0.398241, -2.435270, -0.742095 } },
                {  8, { 1.144646, -3.077963, 3.498864, -1.812710, 3.863297, -3.088573, -0.357800, 1.094312, 2.608025, 0.763815, -2.637127, -0.243639 } },
                { 11, { -0.653880, -3.084803, 2.698784, -0.099526, 3.214588, -3.149954, 1.215030, 1.540326, 1.840796, 1.650181, -0.667796, -0.868043 } },
                {  3, { -1.773946, -3.882980, 2.418979, 1.029116, 1.811975, -0.760381, 3.600941, 2.431781, 1.015779, 0.656494, 0.547274, -2.449650 } },
                {  1, { 0.024580, -3.876140, 3.219060, -0.684069, 2.460685, -0.699000, 2.028111, 1.985767, 1.783008, -0.229872, -1.422057, -1.825247 } },
                {  6, { 1.717664, -4.632855, 3.580720, -2.458102, 4.272392, -2.571970, 1.137432, 2.330893, 2.806471, 0.135703, -1.623914, -1.326792 } },
                { 15, { -0.080862, -4.639694, 2.780640, -0.744918, 3.623682, -2.633351, 2.710262, 2.776907, 2.039242, 1.022068, 0.345417, -1.951195 } },
                { 12, { -1.200928, -5.437871, 2.500835, 0.383724, 2.221070, -0.243779, 5.096172, 3.668363, 1.214225, 0.028381, 1.560487, -3.532803 } },
                {  2, { 0.597598, -5.431032, 3.300916, -1.329461, 2.869779, -0.182398, 3.523342, 3.222349, 1.981454, -0.857984, -0.408844, -2.908399 } },
                { 13, { 2.290682, -6.187747, 3.662577, -3.103494, 4.681486, -2.055368, 2.632663, 3.567475, 3.004916, -0.492410, -0.610701, -2.409944 } },
                { 14, { 0.492156, -6.194586, 2.862496, -1.390310, 4.032776, -2.116749, 4.205493, 4.013489, 2.237687, 0.393955, 1.358630, -3.034347 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroVelocitiesResult = reader.getMacroVelocities( );

    for ( auto it = macroVelocitiesAnswer.begin( ); it != macroVelocitiesAnswer.end( ); it++ ){

        auto r = macroVelocitiesResult->find( it->first );

        if ( r == macroVelocitiesResult->end( ) ){

            results << "test_initializeIncrement (test 55) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 56) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroAccelerationsAnswer
        =
            {
               {  5, { 0.427196, 0.928039, 1.624912, -2.356312, 1.564029, -1.719566, 0.672795, -0.997762, -0.440505, 1.065985, 0.969844, 2.956951 } },
               {  9, { 1.700591, 2.076310, 2.157836, -1.571605, 3.403657, -2.637085, 2.643157, 0.496112, -0.696747, 2.936684, -0.903955, 4.777131 } },
               {  8, { 2.377492, 3.010728, 1.465971, -0.839560, 1.625823, -0.674050, 4.129535, -0.565470, -1.308817, 3.939671, 1.079660, 5.432814 } },
               { 11, { 1.104097, 1.862456, 0.933047, -1.624267, -0.213804, 0.243469, 2.159173, -2.059344, -1.052576, 2.068972, 2.953459, 3.612634 } },
               {  3, { -0.120091, 1.778900, 3.245842, -3.033673, 0.251196, -1.879919, 1.501775, -1.509271, 0.832536, 0.079488, 0.817082, 4.057459 } },
               {  1, { 1.153304, 2.927172, 3.778766, -2.248965, 2.090823, -2.797438, 3.472137, -0.015397, 0.576294, 1.950187, -1.056717, 5.877639 } },
               {  6, { 1.830205, 3.861590, 3.086901, -1.516920, 0.312989, -0.834402, 4.958514, -1.076979, -0.035777, 2.953174, 0.926898, 6.533322 } },
               { 15, { 0.556810, 2.713318, 2.553977, -2.301628, -1.526638, 0.083117, 2.988152, -2.570853, 0.220465, 1.082475, 2.800697, 4.713142 } },
               { 12, { -0.667378, 2.629762, 4.866772, -3.711033, -1.061638, -2.040271, 2.330754, -2.020779, 2.105576, -0.907009, 0.664319, 5.157967 } },
               {  2, { 0.606017, 3.778034, 5.399696, -2.926326, 0.777989, -2.957790, 4.301116, -0.526906, 1.849334, 0.963690, -1.209479, 6.978147 } },
               { 13, { 1.282918, 4.712451, 4.707831, -2.194281, -0.999844, -0.994755, 5.787494, -1.588488, 1.237264, 1.966677, 0.774135, 7.633830 } },
               { 14, { 0.009523, 3.564180, 4.174907, -2.978988, -2.839471, -0.077236, 3.817132, -3.082362, 1.493506, 0.095978, 2.647934, 5.813650 } },
           };

    const std::unordered_map< uIntType, floatVector > *macroAccelerationsResult = reader.getMacroAccelerations( );

    for ( auto it = macroAccelerationsAnswer.begin( ); it != macroAccelerationsAnswer.end( ); it++ ){

        auto r = macroAccelerationsResult->find( it->first );

        if ( r == macroAccelerationsResult->end( ) ){

            results << "test_initializeIncrement (test 57) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 58) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > previousMacroDispDOFVectorAnswer
        =
            {
                {  5, { 0.028316, 1.669368, -0.675497, 1.237413, 1.695543, -0.472972, 0.302567, 0.847229, 1.313857, 0.027499, 0.297804, -0.165902 } },
                {  9, { -0.890446, 3.301799, 0.693373, 2.508561, 1.039322, 0.749852, 0.103035, 1.330944, 1.367395, 1.040491, 0.257210, -0.202549 } },
                {  8, { -2.485736, 3.768030, -0.910199, 3.886732, 2.231470, -0.613407, -1.305314, -0.617256, 0.410025, -0.642795, 1.783136, -2.150267 } },
                { 11, { -1.566974, 2.135599, -2.279070, 2.615585, 2.887691, -1.836232, -1.105782, -1.100970, 0.356487, -1.655787, 1.823730, -2.113620 } },
                {  3, { 0.250830, 2.577555, -0.963926, 2.311388, 3.540607, -0.316540, 1.833891, 0.855143, 2.547578, -1.801443, -0.031338, -1.211966 } },
                {  1, { -0.667932, 4.209986, 0.404945, 3.582536, 2.884385, 0.906284, 1.634359, 1.338857, 2.601116, -0.788451, -0.071932, -1.248613 } },
                {  6, { -2.263222, 4.676217, -1.198628, 4.960707, 4.076534, -0.456975, 0.226010, -0.609342, 1.643746, -2.471737, 1.453994, -3.196331 } },
                { 15, { -1.344460, 3.043786, -2.567498, 3.689560, 4.732755, -1.679800, 0.425541, -1.093057, 1.590208, -3.484729, 1.494589, -3.159684 } },
                { 12, { 0.473344, 3.485742, -1.252354, 3.385363, 5.385670, -0.160109, 3.365215, 0.863056, 3.781300, -3.630385, -0.360479, -2.258030 } },
                {  2, { -0.445418, 5.118173, 0.116516, 4.656510, 4.729449, 1.062716, 3.165683, 1.346770, 3.834837, -2.617393, -0.401074, -2.294677 } },
                { 13, { -2.040708, 5.584404, -1.487057, 6.034682, 5.921597, -0.300544, 1.757333, -0.601429, 2.877467, -4.300680, 1.124852, -4.242395 } },
                { 14, { -1.121946, 3.951973, -2.855927, 4.763534, 6.577818, -1.523368, 1.956865, -1.085144, 2.823929, -5.313671, 1.165447, -4.205748 } },
            };

    const std::unordered_map< uIntType, floatVector > *previousMacroDispDOFVectorResult = reader.getPreviousMacroDispDOFVector( );

    for ( auto it = previousMacroDispDOFVectorAnswer.begin( ); it != previousMacroDispDOFVectorAnswer.end( ); it++ ){

        auto r = previousMacroDispDOFVectorResult->find( it->first );

        if ( r == previousMacroDispDOFVectorResult->end( ) ){

            results << "test_initializeIncrement (test 59) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 60) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > previousMacroVelocitiesAnswer
        =
             {
                {  5, { -1.501765, -1.483894, 0.684890, 0.330246, -0.559686, -0.727351, 1.240470, 1.138232, 0.275631, 1.151980, 0.135063, -0.239197 } },
                {  9, { 0.296761, -1.477054, 1.484970, -1.382938, 0.089023, -0.665970, -0.332360, 0.692218, 1.042860, 0.265615, -1.834268, 0.385206 } },
                {  8, { 1.989845, -2.233769, 1.846631, -3.156971, 1.900730, -2.538940, -1.223039, 1.037345, 2.066323, 0.631189, -2.036125, 0.883662 } },
                { 11, { 0.191319, -2.240608, 1.046550, -1.443787, 1.252020, -2.600321, 0.349790, 1.483359, 1.299094, 1.517554, -0.066794, 0.259258 } },
                {  3, { -0.928747, -3.038786, 0.766746, -0.315146, -0.150592, -0.210748, 2.735701, 2.374814, 0.474077, 0.523867, 1.148276, -1.322349 } },
                {  1, { 0.869779, -3.031946, 1.566826, -2.028330, 0.498117, -0.149367, 1.162871, 1.928800, 1.241306, -0.362498, -0.821055, -0.697946 } },
                {  6, { 2.562863, -3.788661, 1.928487, -3.802363, 2.309824, -2.022337, 0.272192, 2.273927, 2.264768, 0.003076, -1.022912, -0.199490 } },
                { 15, { 0.764337, -3.795500, 1.128407, -2.089179, 1.661114, -2.083718, 1.845022, 2.719941, 1.497539, 0.889442, 0.946419, -0.823894 } },
                { 12, { -0.355729, -4.593677, 0.848602, -0.960538, 0.258502, 0.305854, 4.230932, 3.611396, 0.672522, -0.104245, 2.161489, -2.405502 } },
                {  2, { 1.442797, -4.586838, 1.648683, -2.673722, 0.907211, 0.367235, 2.658103, 3.165382, 1.439751, -0.990611, 0.192157, -1.781098 } },
                { 13, { 3.135881, -5.343553, 2.010343, -4.447755, 2.718918, -1.505734, 1.767423, 3.510508, 2.463214, -0.625036, -0.009700, -1.282643 } },
                { 14, { 1.337355, -5.350392, 1.210263, -2.734571, 2.070208, -1.567116, 3.340253, 3.956522, 1.695985, 0.261329, 1.959631, -1.907046 } },
            };

    const std::unordered_map< uIntType, floatVector > *previousMacroVelocitiesResult = reader.getPreviousMacroVelocities( );

    for ( auto it = previousMacroVelocitiesAnswer.begin( ); it != previousMacroVelocitiesAnswer.end( ); it++ ){

        auto r = previousMacroVelocitiesResult->find( it->first );

        if ( r == previousMacroVelocitiesResult->end( ) ){

            results << "test_initializeIncrement (test 61) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 62) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > previousMacroAccelerationsAnswer
        =
            {
               {  5, { 1.590337, -0.180779, 0.211887, -0.897186, 0.426442, -1.677882, -1.198663, 0.035661, 1.004663, 0.476809, -0.315302, 1.383007 } },
               {  9, { 2.863731, 0.967492, 0.744811, -0.112479, 2.266069, -2.595401, 0.771698, 1.529535, 0.748421, 2.347508, -2.189100, 3.203187 } },
               {  8, { 3.540632, 1.901910, 0.052946, 0.619566, 0.488235, -0.632365, 2.258076, 0.467953, 0.136351, 3.350494, -0.205486, 3.858870 } },
               { 11, { 2.267237, 0.753638, -0.479978, -0.165141, -1.351392, 0.285154, 0.287714, -1.025921, 0.392593, 1.479796, 1.668313, 2.038690 } },
               {  3, { 1.043050, 0.670082, 1.832817, -1.574547, -0.886392, -1.838234, -0.369684, -0.475848, 2.277704, -0.509688, -0.468064, 2.483515 } },
               {  1, { 2.316444, 1.818354, 2.365741, -0.789840, 0.953235, -2.755753, 1.600678, 1.018026, 2.021462, 1.361011, -2.341863, 4.303695 } },
               {  6, { 2.993345, 2.752772, 1.673876, -0.057795, -0.824598, -0.792718, 3.087056, -0.043556, 1.409392, 2.363997, -0.358248, 4.959378 } },
               { 15, { 1.719950, 1.604500, 1.140952, -0.842502, -2.664225, 0.124801, 1.116694, -1.537430, 1.665633, 0.493299, 1.515551, 3.139198 } },
               { 12, { 0.495763, 1.520944, 3.453747, -2.251908, -2.199226, -1.998587, 0.459296, -0.987357, 3.550745, -1.496185, -0.620827, 3.584022 } },
               {  2, { 1.769158, 2.669216, 3.986671, -1.467200, -0.359598, -2.916106, 2.429658, 0.506517, 3.294503, 0.374514, -2.494625, 5.404202 } },
               { 13, { 2.446058, 3.603633, 3.294806, -0.735155, -2.137432, -0.953071, 3.916035, -0.555065, 2.682432, 1.377500, -0.511011, 6.059886 } },
               { 14, { 1.172663, 2.455361, 2.761882, -1.519863, -3.977059, -0.035552, 1.945673, -2.048939, 2.938674, -0.493198, 1.362788, 4.239706 } },
           };

    const std::unordered_map< uIntType, floatVector > *previousMacroAccelerationsResult = reader.getPreviousMacroAccelerations( );

    for ( auto it = previousMacroAccelerationsAnswer.begin( ); it != previousMacroAccelerationsAnswer.end( ); it++ ){

        auto r = previousMacroAccelerationsResult->find( it->first );

        if ( r == previousMacroAccelerationsResult->end( ) ){

            results << "test_initializeIncrement (test 63) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 64) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroInternalForcesAnswer
        =
            {
                {  5, { -3.123250, 0.977401, 2.054240, 1.757330, -0.474837, 1.152554, 0.177148, -0.110803, -0.822029, -0.726168, -0.869646, 0.163025 } },
                {  9, { -1.936939, 2.067285, 3.129559, 1.307129, 1.073546, 1.568170, 0.669365, 0.549069, -1.706006, 0.098907, -0.200723, 0.426380 } },
                {  8, { -2.811405, 2.316517, 3.559565, 0.391797, 2.169330, 0.888284, -0.794757, -1.286124, 0.075576, 0.364483, 1.119236, -0.465277 } },
                { 11, { -3.997715, 1.226633, 2.484245, 0.841998, 0.620946, 0.472668, -1.286974, -1.945997, 0.959553, -0.460592, 0.450314, -0.728632 } },
                {  3, { -3.134236, 2.438781, 3.072693, 3.369779, -1.079855, 0.553564, -0.716138, 0.891244, -2.584717, -2.198251, 0.401625, 0.010912 } },
                {  1, { -1.947925, 3.528665, 4.148013, 2.919578, 0.468528, 0.969180, -0.223921, 1.551116, -3.468694, -1.373176, 1.070548, 0.274268 } },
                {  6, { -2.822390, 3.777897, 4.578018, 2.004247, 1.564312, 0.289294, -1.688043, -0.284077, -1.687111, -1.107600, 2.390507, -0.617389 } },
                { 15, { -4.008701, 2.688013, 3.502699, 2.454448, 0.015928, -0.126322, -2.180260, -0.943950, -0.803134, -1.932675, 1.721585, -0.880745 } },
                { 12, { -3.145221, 3.900161, 4.091147, 4.982229, -1.684873, -0.045426, -1.609424, 1.893290, -4.347404, -3.670334, 1.672896, -0.141200 } },
                {  2, { -1.958911, 4.990045, 5.166466, 4.532028, -0.136490, 0.370190, -1.117208, 2.553163, -5.231381, -2.845259, 2.341819, 0.122155 } },
                { 13, { -2.833376, 5.239278, 5.596472, 3.616696, 0.959294, -0.309696, -2.581330, 0.717969, -3.449799, -2.579683, 3.661778, -0.769502 } },
                { 14, { -4.019687, 4.149394, 4.521152, 4.066897, -0.589089, -0.725312, -3.073546, 0.058097, -2.565822, -3.404757, 2.992856, -1.032857 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroInternalForcesResult = reader.getMacroInternalForces( );

    for ( auto it = macroInternalForcesAnswer.begin( ); it != macroInternalForcesAnswer.end( ); it++ ){

        auto r = macroInternalForcesResult->find( it->first );

        if ( r == macroInternalForcesResult->end( ) ){

            results << "test_initializeIncrement (test 65) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 66) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroInertialForcesAnswer
        =
            {
                {  5, { 1.964831, 1.278634, -1.125705, -2.353362, 0.113154, -1.589520, 0.935279, 1.013984, 2.260416, 3.108513, -2.500627, 0.826868 } },
                {  9, { 2.191268, 0.761777, -1.956907, -1.913616, -0.753844, -2.572974, -0.979177, 1.534873, 3.600554, 2.755534, -2.233762, 2.617534 } },
                {  8, { 1.628217, 1.471129, -1.883881, -1.853679, -2.038739, -1.016306, -1.413380, 0.733672, 2.539391, 4.310301, -0.946690, 2.543249 } },
                { 11, { 1.401781, 1.987986, -1.052679, -2.293425, -1.171740, -0.032851, 0.501075, 0.212783, 1.199253, 4.663281, -1.213556, 0.752583 } },
                {  3, { 1.065332, 3.110252, -0.267265, -0.927058, -1.756161, -2.313741, -0.023437, -0.393319, 1.532658, 2.038427, -4.150676, 2.771371 } },
                {  1, { 1.291768, 2.593395, -1.098467, -0.487312, -2.623159, -3.297196, -1.937893, 0.127570, 2.872797, 1.685447, -3.883810, 4.562038 } },
                {  6, { 0.728718, 3.302747, -1.025441, -0.427375, -3.908054, -1.740527, -2.372096, -0.673631, 1.811634, 3.240215, -2.596738, 4.487753 } },
                { 15, { 0.502282, 3.819604, -0.194239, -0.867121, -3.041055, -0.757073, -0.457640, -1.194520, 0.471496, 3.593195, -2.863604, 2.697087 } },
                { 12, { 0.165833, 4.941870, 0.591175, 0.499246, -3.625476, -3.037963, -0.982153, -1.800622, 0.804901, 0.968341, -5.800724, 4.715875 } },
                {  2, { 0.392269, 4.425013, -0.240027, 0.938992, -4.492474, -4.021417, -2.896609, -1.279733, 2.145039, 0.615361, -5.533858, 6.506542 } },
                { 13, { -0.170781, 5.134365, -0.167001, 0.998929, -5.777368, -2.464749, -3.330812, -2.080934, 1.083877, 2.170129, -4.246787, 6.432257 } },
                { 14, { -0.397217, 5.651222, 0.664201, 0.559183, -4.910370, -1.481294, -1.416356, -2.601823, -0.256262, 2.523108, -4.513652, 4.641591 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroInertialForcesResult = reader.getMacroInertialForces( );

    for ( auto it = macroInertialForcesAnswer.begin( ); it != macroInertialForcesAnswer.end( ); it++ ){

        auto r = macroInertialForcesResult->find( it->first );

        if ( r == macroInertialForcesResult->end( ) ){

            results << "test_initializeIncrement (test 67) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 68) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroBodyForcesAnswer
        =
            {
                {  5, { -0.585148, 2.009889, -0.810035, -2.353495, -1.481933, -1.394433, 2.604481, -0.764513, 0.620197, 2.312799, -2.676288, 0.791100 } },
                {  9, { 0.281864, 2.430634, 1.088528, -2.990224, -3.411629, -1.880188, 3.117476, -2.050000, 0.797479, 2.325336, -3.771395, -0.904908 } },
                {  8, { 2.143468, 1.904844, 1.400196, -1.881262, -3.641681, -0.353682, 3.474300, -0.518626, 0.209746, 0.757177, -5.497756, -2.613388 } },
                { 11, { 1.276457, 1.484099, -0.498368, -1.244534, -1.711985, 0.132072, 2.961305, 0.766861, 0.032463, 0.744640, -4.402649, -0.917380 } },
                {  3, { -2.214053, 1.620183, -2.790640, -3.412234, -0.154223, -2.107009, 1.316934, -2.590009, 1.051996, 1.727722, -3.562681, -0.900240 } },
                {  1, { -1.347041, 2.040929, -0.892076, -4.048962, -2.083919, -2.592763, 1.829929, -3.875496, 1.229279, 1.740259, -4.657788, -2.596248 } },
                {  6, { 0.514563, 1.515139, -0.580409, -2.940001, -2.313971, -1.066258, 2.186754, -2.344123, 0.641545, 0.172100, -6.384149, -4.304728 } },
                { 15, { -0.352449, 1.094394, -2.478973, -2.303272, -0.384275, -0.580504, 1.673759, -1.058635, 0.464263, 0.159563, -5.289042, -2.608720 } },
                { 12, { -3.842958, 1.230478, -4.771245, -4.470972, 1.173487, -2.819585, 0.029388, -4.415505, 1.483796, 1.142644, -4.449073, -2.591580 } },
                {  2, { -2.975947, 1.651223, -2.872681, -5.107700, -0.756209, -3.305339, 0.542383, -5.700993, 1.661078, 1.155182, -5.544180, -4.287588 } },
                { 13, { -1.114343, 1.125434, -2.561014, -3.998739, -0.986262, -1.778834, 0.899208, -4.169619, 1.073344, -0.412977, -7.270541, -5.996068 } },
                { 14, { -1.981354, 0.704688, -4.459577, -3.362010, 0.943435, -1.293080, 0.386213, -2.884132, 0.896062, -0.425514, -6.175435, -4.300060 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroBodyForcesResult = reader.getMacroBodyForces( );

    for ( auto it = macroBodyForcesAnswer.begin( ); it != macroBodyForcesAnswer.end( ); it++ ){

        auto r = macroBodyForcesResult->find( it->first );

        if ( r == macroBodyForcesResult->end( ) ){

            results << "test_initializeIncrement (test 69) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 70) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroSurfaceForcesAnswer
        =
            {
                {  5, { 0.700483, 1.063404, -0.561245, 0.196039, 0.442264, -1.455507, 0.778785, -0.781834, 0.010692, -2.689651, -1.646803, -2.731019 } },
                {  9, { 2.343569, 1.285853, -0.326911, 1.727253, 0.486287, -1.879686, 0.780393, -2.137662, 1.301457, -3.572794, -1.020227, -2.294672 } },
                {  8, { 3.754636, 3.023779, -0.181710, 1.511737, -1.141960, -0.566068, 0.325507, -0.975232, -0.594964, -4.529299, -0.496431, -3.091869 } },
                { 11, { 2.111550, 2.801331, -0.416044, -0.019477, -1.185983, -0.141889, 0.323899, 0.380597, -1.885729, -3.646155, -1.123007, -3.528216 } },
                {  3, { -0.684575, 2.891510, -1.646963, 1.544569, 1.248986, -0.666638, 0.967332, 1.124401, 0.193638, -0.784134, -2.692686, -3.146953 } },
                {  1, { 0.958511, 3.113959, -1.412630, 3.075784, 1.293009, -1.090817, 0.968940, -0.231428, 1.484403, -1.667278, -2.066109, -2.710605 } },
                {  6, { 2.369578, 4.851886, -1.267429, 2.860268, -0.335238, 0.222800, 0.514055, 0.931002, -0.412019, -2.623782, -1.542313, -3.507802 } },
                { 15, { 0.726492, 4.629437, -1.501763, 1.329054, -0.379261, 0.646979, 0.512447, 2.286831, -1.702784, -1.740639, -2.168890, -3.944150 } },
                { 12, { -2.069633, 4.719617, -2.732682, 2.893100, 2.055709, 0.122230, 1.155880, 3.030635, 0.376583, 1.121382, -3.738568, -3.562886 } },
                {  2, { -0.426547, 4.942065, -2.498348, 4.424315, 2.099731, -0.301949, 1.157488, 1.674806, 1.667349, 0.238239, -3.111991, -3.126539 } },
                { 13, { 0.984520, 6.679992, -2.353148, 4.208799, 0.471484, 1.011669, 0.702603, 2.837237, -0.229073, -0.718266, -2.588195, -3.923736 } },
                { 14, { -0.658566, 6.457543, -2.587482, 2.677584, 0.427461, 1.435848, 0.700994, 4.193065, -1.519838, 0.164878, -3.214772, -4.360083 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroSurfaceForcesResult = reader.getMacroSurfaceForces( );

    for ( auto it = macroSurfaceForcesAnswer.begin( ); it != macroSurfaceForcesAnswer.end( ); it++ ){

        auto r = macroSurfaceForcesResult->find( it->first );

        if ( r == macroSurfaceForcesResult->end( ) ){

            results << "test_initializeIncrement (test 71) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 72) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroExternalForcesAnswer
        =
            {
                {  5, { 0.115335, 3.073293, -1.371280, -2.157457, -1.039669, -2.849940, 3.383265, -1.546346, 0.630889, -0.376852, -4.323092, -1.939920 } },
                {  9, { 2.625433, 3.716487, 0.761617, -1.262971, -2.925342, -3.759873, 3.897868, -4.187662, 2.098937, -1.247459, -4.791622, -3.199580 } },
                {  8, { 5.898105, 4.928624, 1.218485, -0.369525, -4.783642, -0.919750, 3.799808, -1.493858, -0.385218, -3.772122, -5.994187, -5.705257 } },
                { 11, { 3.388007, 4.285430, -0.914412, -1.264011, -2.897968, -0.009817, 3.285205, 1.147458, -1.853266, -2.901515, -5.525657, -4.445597 } },
                {  3, { -2.898628, 4.511694, -4.437603, -1.867664, 1.094763, -2.773648, 2.284266, -1.465608, 1.245634, 0.943587, -6.255366, -4.047193 } },
                {  1, { -0.388530, 5.154888, -2.304706, -0.973178, -0.790910, -3.683581, 2.798869, -4.106924, 2.713682, 0.072981, -6.723897, -5.306853 } },
                {  6, { 2.884141, 6.367025, -1.847838, -0.079733, -2.649210, -0.843458, 2.700809, -1.413120, 0.229527, -2.451682, -7.926462, -7.812531 } },
                { 15, { 0.374043, 5.723830, -3.980735, -0.974219, -0.763536, 0.066476, 2.186206, 1.228196, -1.238521, -1.581076, -7.457932, -6.552870 } },
                { 12, { -5.912592, 5.950094, -7.503927, -1.577872, 3.229195, -2.697355, 1.185267, -1.384870, 1.860379, 2.264027, -8.187641, -6.154466 } },
                {  2, { -3.402494, 6.593289, -5.371029, -0.683386, 1.343522, -3.607288, 1.699871, -4.026186, 3.328427, 1.393420, -8.656171, -7.414127 } },
                { 13, { -0.129822, 7.805425, -4.914161, 0.210060, -0.514778, -0.767165, 1.601810, -1.332382, 0.844272, -1.131243, -9.858737, -9.919804 } },
                { 14, { -2.639920, 7.162231, -7.047059, -0.684426, 1.370896, 0.142768, 1.087207, 1.308934, -0.623776, -0.260636, -9.390206, -8.660143 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroExternalForcesResult = reader.getMacroExternalForces( );

    for ( auto it = macroExternalForcesAnswer.begin( ); it != macroExternalForcesAnswer.end( ); it++ ){

        auto r = macroExternalForcesResult->find( it->first );

        if ( r == macroExternalForcesResult->end( ) ){

            results << "test_initializeIncrement (test 73) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement (test 74) & False\n";
            return 1;

        }

    }

    const uIntVector *freeMicroNodeIds = reader.getFreeMicroNodeIds( );
    const uIntVector *ghostMicroNodeIds = reader.getGhostMicroNodeIds( );

    for ( auto n  = ghostMicroNodeIds->begin( );
               n != ghostMicroNodeIds->end( );
               n++ ){

        if ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) != freeMicroNodeIds->end( ) ){
            std::cout << "*n: " << *n << "\n";
            results << "test_initializeIncrement (test 75) & False\n";
            return 1;
        }

    }

    uIntVector nodes;
    const stringVector *freeMicroDomainNames = reader.getFreeMicroDomainNames( );
    for ( auto domain  = freeMicroDomainNames->begin( );
               domain != freeMicroDomainNames->end( );
               domain++ ){

        reader._microscale->getSubDomainNodes( 0, *domain, nodes );

        for ( auto n = nodes.begin( ); n != nodes.end( ); n++ ){

            if ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) == freeMicroNodeIds->end( ) ){

                results << "test_initializeIncrement (test 76) & False\n";
                return 1;

            }

        }

    }

    const stringVector *ghostMicroDomainNames = reader.getGhostMicroDomainNames( );

    for ( auto domain  = ghostMicroDomainNames->begin( );
               domain != ghostMicroDomainNames->end( );
               domain++ ){

        reader._microscale->getSubDomainNodes( 0, *domain, nodes );

        for ( auto n = nodes.begin( ); n != nodes.end( ); n++ ){

            if ( ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) == freeMicroNodeIds->end( ) ) &&
                 ( std::find( ghostMicroNodeIds->begin( ), ghostMicroNodeIds->end( ), *n ) == ghostMicroNodeIds->end( ) ) ){

                results << "test_initializeIncrement (test 77) & False\n";
                return 1;

            }

        }

    }

    const uIntVector *freeMacroNodeIds = reader.getFreeMacroNodeIds( );
    const uIntVector *ghostMacroNodeIds = reader.getGhostMacroNodeIds( );

    const stringVector *ghostMacroDomainNames = reader.getGhostMacroDomainNames( );
    for ( auto domain  = ghostMacroDomainNames->begin( );
               domain != ghostMacroDomainNames->end( );
               domain++ ){

        reader._macroscale->getSubDomainNodes( 0, *domain, nodes );

        for ( auto n = nodes.begin( ); n != nodes.end( ); n++ ){

            if ( ( std::find( ghostMacroNodeIds->begin( ), ghostMacroNodeIds->end( ), *n ) == ghostMacroNodeIds->end( ) ) ){

                results << "test_initializeIncrement (test 78) & False\n";
                return 1;

            }

        }

    }

    const stringVector *freeMacroDomainNames = reader.getFreeMacroDomainNames( );
    for ( auto domain  = freeMacroDomainNames->begin( );
               domain != freeMacroDomainNames->end( );
               domain++ ){

        reader._macroscale->getSubDomainNodes( 0, *domain, nodes );

        for ( auto n = nodes.begin( ); n != nodes.end( ); n++ ){

            if ( ( std::find( ghostMacroNodeIds->begin( ), ghostMacroNodeIds->end( ), *n ) == ghostMacroNodeIds->end( ) ) &&
                 ( std::find( freeMacroNodeIds->begin( ), freeMacroNodeIds->end( ), *n ) == freeMacroNodeIds->end( ) ) ){

                results << "test_initializeIncrement (test 79) & False\n";
                return 1;

            }

        }

    }

    const DOFMap *microGlobalToLocalDOFMap = reader.getMicroGlobalToLocalDOFMap( );

    if ( microGlobalToLocalDOFMap->size( ) != ( freeMicroNodeIds->size( ) + ghostMicroNodeIds->size( ) ) ){

        results << "test_initializeIncrement (test 80) & False\n";
        return 1;

    }

    for ( auto n  = freeMicroNodeIds->begin( );
               n != freeMicroNodeIds->end( );
               n++ ){

        if ( microGlobalToLocalDOFMap->find( *n ) == microGlobalToLocalDOFMap->end( ) ){

            results << "test_initializeIncrement (test 81) & False\n";
            return 1;

        }

    }

    for ( auto n  = ghostMicroNodeIds->begin( );
               n != ghostMicroNodeIds->end( );
               n++ ){

        if ( microGlobalToLocalDOFMap->find( *n ) == microGlobalToLocalDOFMap->end( ) ){

            results << "test_initializeIncrement (test 82) & False\n";
            return 1;

        }

    }

    const DOFMap *macroGlobalToLocalDOFMap = reader.getMacroGlobalToLocalDOFMap( );

    if ( macroGlobalToLocalDOFMap->size( ) != ( freeMacroNodeIds->size( ) + ghostMacroNodeIds->size( ) ) ){

        results << "test_initializeIncrement (test 83) & False\n";
        return 1;

    }

    for ( auto n  = freeMacroNodeIds->begin( );
               n != freeMacroNodeIds->end( );
               n++ ){

        if ( macroGlobalToLocalDOFMap->find( *n ) == macroGlobalToLocalDOFMap->end( ) ){

            results << "test_initializeIncrement (test 84) & False\n";
            return 1;

        }

    }

    for ( auto n  = ghostMacroNodeIds->begin( );
               n != ghostMacroNodeIds->end( );
               n++ ){

        if ( macroGlobalToLocalDOFMap->find( *n ) == macroGlobalToLocalDOFMap->end( ) ){

            results << "test_initializeIncrement (test 85) & False\n";
            return 1;

        }

    }

    if ( !reader.microBodyForceDefined( ) ){

        results << "test_initializeIncrement (test 86) & False\n";
        return 1;

    }

    if ( !reader.microSurfaceForceDefined( ) ){

        results << "test_initializeIncrement (test 87) & False\n";
        return 1;

    }

    if ( !reader.microAccelerationDefined( ) ){

        results << "test_initializeIncrement (test 88) & False\n";
        return 1;

    }

    if ( reader.useReconstructedMassCenters( ) ){

        results << "test_initializeIncrement (test 89) & False\n";
        return 1;

    }

    if ( !reader.microVelocitiesDefined( ) ){

        results << "test_initializeIncrement (test 90) & False\n";
        return 1;

    }

    if ( !reader.macroAccelerationDefined( ) ){

        results << "test_initializeIncrement (test 91) & False\n";
        return 1;

    }

    if ( !reader.macroVelocitiesDefined( ) ){

        results << "test_initializeIncrement (test 92) & False\n";
        return 1;

    }


    if ( !reader.microInternalForceDefined( ) ){

        results << "test_initializeIncrement (test 93) & False\n";
        return 1;

    }

    if ( !reader.macroInternalForceDefined( ) ){

        results << "test_initializeIncrement (test 94) & False\n";
        return 1;

    }

    if ( !reader.macroInertialForceDefined( ) ){

        results << "test_initializeIncrement (test 95) & False\n";
        return 1;

    }

    if ( !reader.macroExternalForceDefined( ) ){

        results << "test_initializeIncrement (test 96) & False\n";
        return 1;

    }

    const std::string macroReferenceDensityTypesAnswer = "constant";
    const floatVector macroReferenceDensitiesAnswer = { 2. };
    const std::unordered_map< unsigned int, floatVector > *macroReferenceDensitiesResult = reader.getMacroReferenceDensities( );
    const std::unordered_map< unsigned int, std::string > *macroReferenceDensityTypesResult = reader.getMacroReferenceDensityTypes( );

    for ( auto mRDR = macroReferenceDensitiesResult->begin( ); mRDR != macroReferenceDensitiesResult->end( ); mRDR++ ){

        if ( !vectorTools::fuzzyEquals( macroReferenceDensitiesAnswer, mRDR->second ) ){

            results << "test_initializeIncrement (test 97) & False\n";
            return 1;

        }

    }

    for ( auto mRDTR = macroReferenceDensityTypesResult->begin( ); mRDTR != macroReferenceDensityTypesResult->end( ); mRDTR++ ){

        if ( macroReferenceDensityTypesAnswer.compare( mRDTR->second ) != 0 ){

            results << "test_initializeIncrement (test 98) & False\n";
            return 1;

        }

    }

    const std::string macroReferenceMomentOfInertiaTypesAnswer = "constant";
    const floatVector macroReferenceMomentsOfInertiaAnswer = { 1e-5, 2e-5, 3e-5,
                                                               2e-5, 4e-5, 5e-5,
                                                               3e-5, 5e-5, 6e-5 };
    const std::unordered_map< unsigned int, floatVector > *macroReferenceMomentsOfInertiaResult
        = reader.getMacroReferenceMomentsOfInertia( );
    const std::unordered_map< unsigned int, std::string > *macroReferenceMomentOfInertiaTypesResult
        = reader.getMacroReferenceMomentOfInertiaTypes( );

    for ( auto mRMIR  = macroReferenceMomentsOfInertiaResult->begin( );
               mRMIR != macroReferenceMomentsOfInertiaResult->end( );
               mRMIR++ ){

        if ( !vectorTools::fuzzyEquals( macroReferenceMomentsOfInertiaAnswer, mRMIR->second ) ){

            results << "test_initializeIncrement (test 99) & False\n";
            return 1;

        }

    }

    for ( auto mRMITR  = macroReferenceMomentOfInertiaTypesResult->begin( );
               mRMITR != macroReferenceMomentOfInertiaTypesResult->end( ); mRMITR++ ){

        if ( macroReferenceMomentOfInertiaTypesAnswer.compare( mRMITR->second ) != 0 ){

            results << "test_initializeIncrement (test 100) & False\n";
            return 1;

        }

    }

    if ( !reader.microSurfaceForceDefined( ) ){

        results << "test_initializeIncrement (test 101) & False\n";
        return 1;

    }

    if ( !reader.microExternalForceDefined( ) ){

        results << "test_initializeIncrement (test 102) & False\n";
        return 1;

    }

    if ( !reader.extractPreviousDOFValues( ) ){

        results << "test_initializeIncrement (test 103) & False\n";
        return 1;

    }

    const floatType DtAnswer = 1.;
    const floatType* DtResult = reader.getDt( );

    if ( !vectorTools::fuzzyEquals( DtAnswer, *DtResult ) ){

        results << "test_initializeIncrement (test 104) & False\n";
        return 1;

    }

    const floatType newmarkGammaAnswer = 0.50;
    const floatType newmarkBetaAnswer  = 0.25;

    if ( !vectorTools::fuzzyEquals( newmarkGammaAnswer, *reader.getNewmarkGamma( ) ) ){

        results << "test_initializeIncrement (test 105) & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( newmarkBetaAnswer, *reader.getNewmarkBeta( ) ) ){

        results << "test_initializeIncrement (test 106) & False\n";
        return 1;

    }

    const std::unordered_map< uIntType, stringVector > macroCellToDomainMapAnswer
        =
        {
            { 1, { "ghost_nodeset_volume_1", "ghost_nodeset_volume_2", "ghost_nodeset_volume_3", "ghost_nodeset_volume_4",
                   "ghost_nodeset_volume_5", "ghost_nodeset_volume_6", "ghost_nodeset_volume_7", "ghost_nodeset_volume_8" } },
            { 2, { "free_nodeset_volume_1", "free_nodeset_volume_2", "free_nodeset_volume_3", "free_nodeset_volume_4",
                   "free_nodeset_volume_5", "free_nodeset_volume_6", "free_nodeset_volume_7", "free_nodeset_volume_8" } },
        };

    const std::unordered_map< uIntType, stringVector > *macroCellToDomainMapResult = reader.getMacroCellToDomainMap( );

    for ( auto a = macroCellToDomainMapAnswer.begin( ); a != macroCellToDomainMapAnswer.end( ); a++ ){

        auto r = macroCellToDomainMapResult->find( a->first );

        if ( r == macroCellToDomainMapResult->end( ) ){

            results << "test_initializeIncrement (test 107) & False\n";
            return 1;

        }

        if ( a->second.size( ) != r->second.size( ) ){

            results << "test_initializeIncrement (test 108) & False\n";
            return 1;

        }

        for ( unsigned int i = 0; i < a->second.size( ); i++ ){

            if ( a->second[ i ].compare( r->second[ i ] ) != 0 ){

                results << "test_initializeIncrement (test 109) & False\n";
                return 1;

            }

        }        

    }

    const std::unordered_map< std::string, uIntType > microDomainIDMapAnswer
        =
        {
            { "free_nodeset_volume_1",   0 },
            { "free_nodeset_volume_2",   1 },
            { "free_nodeset_volume_3",   2 },
            { "free_nodeset_volume_4",   3 },
            { "free_nodeset_volume_5",   4 },
            { "free_nodeset_volume_6",   5 },
            { "free_nodeset_volume_7",   6 },
            { "free_nodeset_volume_8",   7 },
            { "ghost_nodeset_volume_1",  8 },
            { "ghost_nodeset_volume_2",  9 },
            { "ghost_nodeset_volume_3", 10 },
            { "ghost_nodeset_volume_4", 11 },
            { "ghost_nodeset_volume_5", 12 },
            { "ghost_nodeset_volume_6", 13 },
            { "ghost_nodeset_volume_7", 14 },
            { "ghost_nodeset_volume_8", 15 },
        };

    const std::unordered_map< std::string, uIntType > *microDomainIDMapResult = reader.getMicroDomainIDMap( );

    for ( auto a = microDomainIDMapAnswer.begin( ); a != microDomainIDMapAnswer.end( ); a++ ){

        auto r = microDomainIDMapResult->find( a->first );

        if ( r == microDomainIDMapResult->end( ) ){

            results << "test_initializeIncrement (test 110) & False\n";
            return 1;

        }

        if ( r->second != a->second ){

            results << "test_initializeIncrement (test 111) & False\n";
            return 1;

        }

    }

    results << "test_initializeIncrement & True\n";
    return 0;
}

int test_initializeIncrement_Arlequin( std::ofstream &results ){
    /*!
     * Test the initialization of the processor for the current increment
     * if Arlequin mode is defined
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig_Arlequin.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_initializeIncrement_Arlequin & False\n";
        return 1;
    }

    errorOut error = reader.initializeIncrement( 1, 1 );
    if ( error ){
        error->print( );
        results << "test_initializeIncrement_Arlequin & False\n";
        return 1;
    }

    //Check that the unique micro-scale nodes have been identified
    const DOFMap microGlobalToLocalMapAnswer
        =
            {
                { 15,  0 },
                { 31,  1 },
                { 13,  2 },
                { 26,  3 },
                { 53,  4 },
                { 21,  5 },
                { 37,  6 },
                { 48,  7 },
                {  5,  8 },
                { 10,  9 },
                {  3, 10 },
                {  4, 11 },
                { 32, 12 },
                { 33, 13 },
                { 34, 14 },
                { 28, 15 },
                { 25, 16 },
                { 50, 17 },
                { 43, 18 },
                { 27, 19 },
                {  1, 20 },
                {  7, 21 },
                { 30, 22 },
                { 16, 23 },
                { 22, 24 },
                {  2, 25 },
                { 46, 26 },
                { 24, 27 },
                { 39, 28 },
                { 40, 29 },
                { 57, 30 },
                { 44, 31 },
                { 58, 32 },
                { 29, 33 },
                { 59, 34 },
                { 11, 35 },
                {  0, 36 },
                { 20, 37 },
                { 60, 38 },
                { 47, 39 },
                { 49, 40 },
                { 17, 41 },
                { 38, 42 },
                { 14, 43 },
                { 55, 44 },
            };

    const DOFMap *microGlobalToLocalResult = reader.getMicroGlobalToLocalDOFMap( );

    for ( auto it = microGlobalToLocalMapAnswer.begin( ); it != microGlobalToLocalMapAnswer.end( ); it++ ){

        auto r = microGlobalToLocalResult->find( it->first );

        if ( r == microGlobalToLocalResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 1) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement_Arlequin (test 2) & False\n";
            return 1;

        }

    }

    //Check that the unique macro-scale nodes have been identified
    const DOFMap macroGlobalToLocalMapAnswer
        =
            {
                {  5,  0 },
                {  9,  1 },
                {  8,  2 },
                { 11,  3 },
                {  3,  4 },
                {  1,  5 },
                {  6,  6 },
                { 15,  7 },
                { 12,  8 },
                {  2,  9 },
                { 13, 10 },
                { 14, 11 },
            };

    const DOFMap *macroGlobalToLocalResult = reader.getMacroGlobalToLocalDOFMap( );

    for ( auto it = macroGlobalToLocalMapAnswer.begin( ); it != macroGlobalToLocalMapAnswer.end( ); it++ ){

        auto r = macroGlobalToLocalResult->find( it->first );

        if ( r == macroGlobalToLocalResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 3) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement_Arlequin (test 4) & False\n";
            return 1;

        }

    }

    //Check that the micro node weights are initialized properly
    const std::unordered_map< uIntType, floatType > microNodeWeightsAnswer
        =
            {
                { 24, 1.000 },
                { 39, 0.500 },
                { 15, 0.500 },
                { 31, 0.500 },
                { 43, 1.000 },
                { 40, 0.500 },
                { 57, 0.250 },
                { 13, 0.250 },
                { 26, 0.250 },
                { 27, 0.500 },
                { 11, 1.000 },
                {  0, 0.500 },
                {  5, 0.500 },
                { 10, 0.500 },
                { 30, 1.000 },
                { 44, 0.500 },
                { 58, 0.250 },
                { 53, 0.250 },
                { 21, 0.250 },
                {  1, 0.500 },
                { 29, 0.250 },
                { 59, 0.125 },
                { 37, 0.125 },
                { 48, 0.125 },
                {  7, 0.250 },
                { 20, 0.500 },
                { 60, 0.250 },
                {  3, 0.250 },
                {  4, 0.250 },
                { 16, 0.500 },
                { 14, 1.000 },
                { 55, 0.500 },
                { 25, 0.500 },
                { 50, 0.500 },
                { 46, 1.000 },
                { 47, 0.500 },
                { 49, 0.250 },
                { 32, 0.250 },
                { 33, 0.250 },
                { 22, 0.500 },
                { 17, 1.000 },
                { 38, 0.500 },
                { 34, 0.500 },
                { 28, 0.500 },
                {  2, 1.000 },
            };

    const std::unordered_map< uIntType, floatType > *microNodeWeightsResult = reader.getMicroWeights( );

    for ( auto it = microNodeWeightsAnswer.begin( ); it != microNodeWeightsAnswer.end( ); it++ ){

        auto r = microNodeWeightsResult->find( it->first );

        if ( r == microNodeWeightsResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 5) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement_Arlequin (test 6) & False\n";
            return 1;

        }

    }

    //Make sure the micro global node id to output index map has been extracted correctly
    const DOFMap microGlobalNodeToOutputMapAnswer
        =
            {
                { 15,  2 },
                { 31,  3 },
                { 13,  9 },
                { 26, 10 },
                { 53, 23 },
                { 21, 24 },
                { 37, 30 },
                { 48, 31 },
                {  5, 16 },
                { 10, 17 },
                {  3, 37 },
                {  4, 38 },
                { 32, 51 },
                { 33, 52 },
                { 34, 58 },
                { 28, 59 },
                { 25, 44 },
                { 50, 45 },
                { 43,  4 },
                { 27, 11 },
                {  1, 25 },
                {  7, 32 },
                { 30, 18 },
                { 16, 39 },
                { 22, 53 },
                {  2, 60 },
                { 46, 46 },
                { 24,  0 },
                { 39,  1 },
                { 40,  7 },
                { 57,  8 },
                { 44, 21 },
                { 58, 22 },
                { 29, 28 },
                { 59, 29 },
                { 11, 14 },
                {  0, 15 },
                { 20, 35 },
                { 60, 36 },
                { 47, 49 },
                { 49, 50 },
                { 17, 56 },
                { 38, 57 },
                { 14, 42 },
                { 55, 43 },
            };

    const DOFMap *microGlobalNodeToOutputMapResult = reader.getMicroNodeIDOutputIndex( );

    for ( auto it = microGlobalNodeToOutputMapAnswer.begin( ); it != microGlobalNodeToOutputMapAnswer.end( ); it++ ){

        auto r = microGlobalNodeToOutputMapResult->find( it->first );

        if ( r == microGlobalNodeToOutputMapResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 7) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement_Arlequin (test 8) & False\n";
            return 1;

        }

    }

    //Make sure the macro global node id to output index map has been extracted correctly
    const DOFMap macroGlobalNodeToOutputMapAnswer
        =
            {
                {  5,  4 },
                {  9,  5 },
                {  8,  6 },
                { 11,  7 },
                {  3,  8 },
                {  1,  9 },
                {  6, 10 },
                { 15, 11 },
                { 12, 12 },
                {  2, 13 },
                { 13, 14 },
                { 14, 15 },
            };

    const DOFMap *macroGlobalNodeToOutputMapResult = reader.getMacroNodeIDOutputIndex( );

    for ( auto it = macroGlobalNodeToOutputMapAnswer.begin( ); it != macroGlobalNodeToOutputMapAnswer.end( ); it++ ){

        auto r = macroGlobalNodeToOutputMapResult->find( it->first );

        if ( r == macroGlobalNodeToOutputMapResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 9) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement_Arlequin (test 10) & False\n";
            return 1;

        }

    }

    //Make sure the time of the micro increment has been extracted correctly
    const floatType timeAnswer = 1;

    const floatType *timeResult = reader.getMicroTime( );

    if ( !vectorTools::fuzzyEquals( timeAnswer, *timeResult ) ){

        results << "test_initializeIncrement_Arlequin (test 11) & False\n";
        return 1;

    }

    const std::unordered_map< uIntType, floatType > densityAnswer
        =
        {
            { 15, 6.000 },
            { 31, 8.000 },
            { 13, 7.000 },
            { 26, 9.000 },
            { 53, 4.500 },
            { 21, 6.500 },
            { 37, 5.500 },
            { 48, 7.500 },
            {  5, 8.000 },
            { 10, 10.000 },
            {  3, 6.500 },
            {  4, 8.500 },
            { 32, 4.000 },
            { 33, 6.000 },
            { 34, 5.000 },
            { 28, 7.000 },
            { 25, 3.000 },
            { 50, 5.000 },
            { 43, 10.000 },
            { 27, 11.000 },
            {  1, 8.500 },
            {  7, 9.500 },
            { 30, 12.000 },
            { 16, 10.500 },
            { 22, 8.000 },
            {  2, 9.000 },
            { 46, 7.000 },
            { 24, 2.000 },
            { 39, 4.000 },
            { 40, 3.000 },
            { 57, 5.000 },
            { 44, 0.500 },
            { 58, 2.500 },
            { 29, 1.500 },
            { 59, 3.500 },
            { 11, 4.000 },
            {  0, 6.000 },
            { 20, 2.500 },
            { 60, 4.500 },
            { 47, 0.000 },
            { 49, 2.000 },
            { 17, 1.000 },
            { 38, 3.000 },
            { 14, -1.000 },
            { 55, 1.000 },
        };

    const std::unordered_map< uIntType, floatType > *densityResult = reader.getMicroDensities( );

    for ( auto it = densityAnswer.begin( ); it != densityAnswer.end( ); it++ ){

        auto r = densityResult->find( it->first );

        if ( r == densityResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 12) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement_Arlequin (test 13) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatType > volumeAnswer
        =
         {
            { 15, -2.000 },
            { 31, -3.500 },
            { 13, 0.500 },
            { 26, -1.000 },
            { 53, -1.800 },
            { 21, -3.300 },
            { 37, 0.700 },
            { 48, -0.800 },
            {  5, 3.000 },
            { 10, 1.500 },
            {  3, 3.200 },
            {  4, 1.700 },
            { 32, 0.900 },
            { 33, -0.600 },
            { 34, 3.400 },
            { 28, 1.900 },
            { 25, -1.600 },
            { 50, -3.100 },
            { 43, -5.000 },
            { 27, -2.500 },
            {  1, -4.800 },
            {  7, -2.300 },
            { 30, 0.000 },
            { 16, 0.200 },
            { 22, -2.100 },
            {  2, 0.400 },
            { 46, -4.600 },
            { 24, 1.000 },
            { 39, -0.500 },
            { 40, 3.500 },
            { 57, 2.000 },
            { 44, 1.200 },
            { 58, -0.300 },
            { 29, 3.700 },
            { 59, 2.200 },
            { 11, 6.000 },
            {  0, 4.500 },
            { 20, 6.200 },
            { 60, 4.700 },
            { 47, 3.900 },
            { 49, 2.400 },
            { 17, 6.400 },
            { 38, 4.900 },
            { 14, 1.400 },
            { 55, -0.100 },
        };

    const std::unordered_map< uIntType, floatType > *volumeResult = reader.getMicroVolumes( );

    for ( auto it = volumeAnswer.begin( ); it != volumeAnswer.end( ); it++ ){

        auto r = volumeResult->find( it->first );

        if ( r == volumeResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 14) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": " << r->second << "\n";
            std::cout << it->first << ": " << it->second << "\n";
            results << "test_initializeIncrement_Arlequin (test 15) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microNodeReferencePositionAnswer
        =
            {
                { 15, { 0.00, 0.00, 1.00 } },
                { 31, { 0.00, 0.00, 1.50 } },
                { 13, { 0.50, 0.00, 1.00 } },
                { 26, { 0.50, 0.00, 1.50 } },
                { 53, { 0.00, 0.50, 1.00 } },
                { 21, { 0.00, 0.50, 1.50 } },
                { 37, { 0.50, 0.50, 1.00 } },
                { 48, { 0.50, 0.50, 1.50 } },
                {  5, { 1.00, 0.00, 1.00 } },
                { 10, { 1.00, 0.00, 1.50 } },
                {  3, { 1.00, 0.50, 1.00 } },
                {  4, { 1.00, 0.50, 1.50 } },
                { 32, { 0.50, 1.00, 1.00 } },
                { 33, { 0.50, 1.00, 1.50 } },
                { 34, { 1.00, 1.00, 1.00 } },
                { 28, { 1.00, 1.00, 1.50 } },
                { 25, { 0.00, 1.00, 1.00 } },
                { 50, { 0.00, 1.00, 1.50 } },
                { 43, { 0.00, 0.00, 2.00 } },
                { 27, { 0.50, 0.00, 2.00 } },
                {  1, { 0.00, 0.50, 2.00 } },
                {  7, { 0.50, 0.50, 2.00 } },
                { 30, { 1.00, 0.00, 2.00 } },
                { 16, { 1.00, 0.50, 2.00 } },
                { 22, { 0.50, 1.00, 2.00 } },
                {  2, { 1.00, 1.00, 2.00 } },
                { 46, { 0.00, 1.00, 2.00 } },
                { 24, { 0.00, 0.00, 0.00 } },
                { 39, { 0.00, 0.00, 0.50 } },
                { 40, { 0.50, 0.00, 0.00 } },
                { 57, { 0.50, 0.00, 0.50 } },
                { 44, { 0.00, 0.50, 0.00 } },
                { 58, { 0.00, 0.50, 0.50 } },
                { 29, { 0.50, 0.50, 0.00 } },
                { 59, { 0.50, 0.50, 0.50 } },
                { 11, { 1.00, 0.00, 0.00 } },
                {  0, { 1.00, 0.00, 0.50 } },
                { 20, { 1.00, 0.50, 0.00 } },
                { 60, { 1.00, 0.50, 0.50 } },
                { 47, { 0.50, 1.00, 0.00 } },
                { 49, { 0.50, 1.00, 0.50 } },
                { 17, { 1.00, 1.00, 0.00 } },
                { 38, { 1.00, 1.00, 0.50 } },
                { 14, { 0.00, 1.00, 0.00 } },
                { 55, { 0.00, 1.00, 0.50 } },
            };

    const std::unordered_map< uIntType, floatVector > *microNodeReferencePositionResult = reader.getMicroNodeReferencePositions( );

    for ( auto it = microNodeReferencePositionAnswer.begin( ); it != microNodeReferencePositionAnswer.end( ); it++ ){

        auto r = microNodeReferencePositionResult->find( it->first );

        if ( r == microNodeReferencePositionResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 16) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 17) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroNodeReferencePositionsAnswer
        =
            {
                {  5, { 0.000, 0.000, 0.000 } },
                {  9, { 1.000, 0.000, 0.000 } },
                {  8, { 1.000, 1.000, 0.000 } },
                { 11, { 0.000, 1.000, 0.000 } },
                {  3, { 0.000, 0.000, 1.000 } },
                {  1, { 1.000, 0.000, 1.000 } },
                {  6, { 1.000, 1.000, 1.000 } },
                { 15, { 0.000, 1.000, 1.000 } },
                { 12, { 0.000, 0.000, 2.000 } },
                {  2, { 1.000, 0.000, 2.000 } },
                { 13, { 1.000, 1.000, 2.000 } },
                { 14, { 0.000, 1.000, 2.000 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroNodeReferencePositionsResult = reader.getMacroNodeReferencePositions( );

    for ( auto it = macroNodeReferencePositionsAnswer.begin( ); it != macroNodeReferencePositionsAnswer.end( ); it++ ){

        auto r = macroNodeReferencePositionsResult->find( it->first );

        if ( r == macroNodeReferencePositionsResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 18) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 19) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, uIntVector > macroNodeReferenceConnectivityAnswer
        =
        {
            {  1, {  9,  5,  9,  8, 11,  3,  1,  6, 15 } },
            {  2, {  9,  3,  1,  6, 15, 12,  2, 13, 14 } },
        };

    const std::unordered_map< uIntType, uIntVector > *macroNodeReferenceConnectivityResult = reader.getMacroNodeReferenceConnectivity( );

    for ( auto it = macroNodeReferenceConnectivityAnswer.begin( ); it != macroNodeReferenceConnectivityAnswer.end( ); it++ ){

        auto r = macroNodeReferenceConnectivityResult->find( it->first );

        if ( r == macroNodeReferenceConnectivityResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 21) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 22) & False\n";
            return 1;

        }

    }

    std::unordered_map< uIntType, floatVector > microDisplacementAnswer
        =
            {
               { 15, { 0.000, -14.400, -30.400 } },
               { 31, { 0.000, -14.400, -32.000 } },
               { 13, { 2.000, -14.400, -30.400 } },
               { 26, { 2.000, -14.400, -32.000 } },
               { 53, { 0.000, -11.200, -30.400 } },
               { 21, { 0.000, -11.200, -32.000 } },
               { 37, { 2.000, -11.200, -30.400 } },
               { 48, { 2.000, -11.200, -32.000 } },
               {  5, { 4.000, -14.400, -30.400 } },
               { 10, { 4.000, -14.400, -32.000 } },
               {  3, { 4.000, -11.200, -30.400 } },
               {  4, { 4.000, -11.200, -32.000 } },
               { 32, { 2.000, -8.000, -30.400 } },
               { 33, { 2.000, -8.000, -32.000 } },
               { 34, { 4.000, -8.000, -30.400 } },
               { 28, { 4.000, -8.000, -32.000 } },
               { 25, { 0.000, -8.000, -30.400 } },
               { 50, { 0.000, -8.000, -32.000 } },
               { 43, { 0.000, -14.400, -33.600 } },
               { 27, { 2.000, -14.400, -33.600 } },
               {  1, { 0.000, -11.200, -33.600 } },
               {  7, { 2.000, -11.200, -33.600 } },
               { 30, { 4.000, -14.400, -33.600 } },
               { 16, { 4.000, -11.200, -33.600 } },
               { 22, { 2.000, -8.000, -33.600 } },
               {  2, { 4.000, -8.000, -33.600 } },
               { 46, { 0.000, -8.000, -33.600 } },
               { 24, { 0.000, -14.400, -27.200 } },
               { 39, { 0.000, -14.400, -28.800 } },
               { 40, { 2.000, -14.400, -27.200 } },
               { 57, { 2.000, -14.400, -28.800 } },
               { 44, { 0.000, -11.200, -27.200 } },
               { 58, { 0.000, -11.200, -28.800 } },
               { 29, { 2.000, -11.200, -27.200 } },
               { 59, { 2.000, -11.200, -28.800 } },
               { 11, { 4.000, -14.400, -27.200 } },
               {  0, { 4.000, -14.400, -28.800 } },
               { 20, { 4.000, -11.200, -27.200 } },
               { 60, { 4.000, -11.200, -28.800 } },
               { 47, { 2.000, -8.000, -27.200 } },
               { 49, { 2.000, -8.000, -28.800 } },
               { 17, { 4.000, -8.000, -27.200 } },
               { 38, { 4.000, -8.000, -28.800 } },
               { 14, { 0.000, -8.000, -27.200 } },
               { 55, { 0.000, -8.000, -28.800 } },
           };


    const std::unordered_map< uIntType, floatVector > *microDisplacementResult = reader.getMicroDisplacements( );

    for ( auto it = microDisplacementAnswer.begin( ); it != microDisplacementAnswer.end( ); it++ ){

        auto r = microDisplacementResult->find( it->first );

        if ( r == microDisplacementResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 22) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 23) & False\n";
            return 1;

        }

    }


    const uIntVector freeMacroCellIdsAnswer = { 1 };
    const uIntVector ghostMacroCellIdsAnswer = { 2 };

    const uIntVector *freeMacroCellIdsResult = reader.getFreeMacroCellIds( );
    const uIntVector *ghostMacroCellIdsResult = reader.getGhostMacroCellIds( );

    if ( !vectorTools::fuzzyEquals( freeMacroCellIdsAnswer, *freeMacroCellIdsResult ) ){
        results << "test_initializeIncrement_Arlequin (test 24) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( ghostMacroCellIdsAnswer, *ghostMacroCellIdsResult ) ){
        results << "test_initializeIncrement_Arlequin (test 25) & False\n";
        return 1;
    }

//    const uIntVector freeMacroCellMicroDomainCountsAnswer = { 8 };
//    const uIntVector ghostMacroCellMicroDomainCountsAnswer = { 8 };
//
//    const uIntVector *freeMacroCellMicroDomainCountsResult = reader.getFreeMacroCellMicroDomainCounts( );
//    const uIntVector *ghostMacroCellMicroDomainCountsResult = reader.getGhostMacroCellMicroDomainCounts( );
//
//    if ( !vectorTools::fuzzyEquals( freeMacroCellMicroDomainCountsAnswer, *freeMacroCellMicroDomainCountsResult ) ){
//        results << "test_initializeIncrement_Arlequin (test 26) & False\n";
//        return 1;
//    }
//
//    if ( !vectorTools::fuzzyEquals( ghostMacroCellMicroDomainCountsAnswer, *ghostMacroCellMicroDomainCountsResult ) ){
//        results << "test_initializeIncrement_Arlequin (test 27) & False\n";
//        return 1;
//    }

    const std::unordered_map< uIntType, floatVector > microBodyForcesAnswer
        =
            {
                { 15, { 2.000, 2.000, 9.200 } },
                { 31, { 2.000, 2.000, 12.800 } },
                { 13, { 5.000, 2.000, 9.200 } },
                { 26, { 5.000, 2.000, 12.800 } },
                { 53, { 2.000, 2.410, 9.200 } },
                { 21, { 2.000, 2.410, 12.800 } },
                { 37, { 5.000, 2.410, 9.200 } },
                { 48, { 5.000, 2.410, 12.800 } },
                {  5, { 8.000, 2.000, 9.200 } },
                { 10, { 8.000, 2.000, 12.800 } },
                {  3, { 8.000, 2.410, 9.200 } },
                {  4, { 8.000, 2.410, 12.800 } },
                { 32, { 5.000, 2.820, 9.200 } },
                { 33, { 5.000, 2.820, 12.800 } },
                { 34, { 8.000, 2.820, 9.200 } },
                { 28, { 8.000, 2.820, 12.800 } },
                { 25, { 2.000, 2.820, 9.200 } },
                { 50, { 2.000, 2.820, 12.800 } },
                { 43, { 2.000, 2.000, 16.400 } },
                { 27, { 5.000, 2.000, 16.400 } },
                {  1, { 2.000, 2.410, 16.400 } },
                {  7, { 5.000, 2.410, 16.400 } },
                { 30, { 8.000, 2.000, 16.400 } },
                { 16, { 8.000, 2.410, 16.400 } },
                { 22, { 5.000, 2.820, 16.400 } },
                {  2, { 8.000, 2.820, 16.400 } },
                { 46, { 2.000, 2.820, 16.400 } },
                { 24, { 2.000, 2.000, 2.000 } },
                { 39, { 2.000, 2.000, 5.600 } },
                { 40, { 5.000, 2.000, 2.000 } },
                { 57, { 5.000, 2.000, 5.600 } },
                { 44, { 2.000, 2.410, 2.000 } },
                { 58, { 2.000, 2.410, 5.600 } },
                { 29, { 5.000, 2.410, 2.000 } },
                { 59, { 5.000, 2.410, 5.600 } },
                { 11, { 8.000, 2.000, 2.000 } },
                {  0, { 8.000, 2.000, 5.600 } },
                { 20, { 8.000, 2.410, 2.000 } },
                { 60, { 8.000, 2.410, 5.600 } },
                { 47, { 5.000, 2.820, 2.000 } },
                { 49, { 5.000, 2.820, 5.600 } },
                { 17, { 8.000, 2.820, 2.000 } },
                { 38, { 8.000, 2.820, 5.600 } },
                { 14, { 2.000, 2.820, 2.000 } },
                { 55, { 2.000, 2.820, 5.600 } },
            };

    const std::unordered_map< uIntType, floatVector > *microBodyForcesResult = reader.getMicroBodyForces( );

    for ( auto it = microBodyForcesAnswer.begin( ); it != microBodyForcesAnswer.end( ); it++ ){

        auto r = microBodyForcesResult->find( it->first );

        if ( r == microBodyForcesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 28) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 29) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microSurfaceForcesAnswer
        =
            {
                { 15, { 2.100, 2.100, 3.540 } },
                { 31, { 2.100, 2.100, 4.260 } },
                { 13, { 2.460, 2.100, 3.540 } },
                { 26, { 2.460, 2.100, 4.260 } },
                { 53, { 2.100, 3.705, 3.540 } },
                { 21, { 2.100, 3.705, 4.260 } },
                { 37, { 2.460, 3.705, 3.540 } },
                { 48, { 2.460, 3.705, 4.260 } },
                {  5, { 2.820, 2.100, 3.540 } },
                { 10, { 2.820, 2.100, 4.260 } },
                {  3, { 2.820, 3.705, 3.540 } },
                {  4, { 2.820, 3.705, 4.260 } },
                { 32, { 2.460, 5.310, 3.540 } },
                { 33, { 2.460, 5.310, 4.260 } },
                { 34, { 2.820, 5.310, 3.540 } },
                { 28, { 2.820, 5.310, 4.260 } },
                { 25, { 2.100, 5.310, 3.540 } },
                { 50, { 2.100, 5.310, 4.260 } },
                { 43, { 2.100, 2.100, 4.980 } },
                { 27, { 2.460, 2.100, 4.980 } },
                {  1, { 2.100, 3.705, 4.980 } },
                {  7, { 2.460, 3.705, 4.980 } },
                { 30, { 2.820, 2.100, 4.980 } },
                { 16, { 2.820, 3.705, 4.980 } },
                { 22, { 2.460, 5.310, 4.980 } },
                {  2, { 2.820, 5.310, 4.980 } },
                { 46, { 2.100, 5.310, 4.980 } },
                { 24, { 2.100, 2.100, 2.100 } },
                { 39, { 2.100, 2.100, 2.820 } },
                { 40, { 2.460, 2.100, 2.100 } },
                { 57, { 2.460, 2.100, 2.820 } },
                { 44, { 2.100, 3.705, 2.100 } },
                { 58, { 2.100, 3.705, 2.820 } },
                { 29, { 2.460, 3.705, 2.100 } },
                { 59, { 2.460, 3.705, 2.820 } },
                { 11, { 2.820, 2.100, 2.100 } },
                {  0, { 2.820, 2.100, 2.820 } },
                { 20, { 2.820, 3.705, 2.100 } },
                { 60, { 2.820, 3.705, 2.820 } },
                { 47, { 2.460, 5.310, 2.100 } },
                { 49, { 2.460, 5.310, 2.820 } },
                { 17, { 2.820, 5.310, 2.100 } },
                { 38, { 2.820, 5.310, 2.820 } },
                { 14, { 2.100, 5.310, 2.100 } },
                { 55, { 2.100, 5.310, 2.820 } },
            };

    const std::unordered_map< uIntType, floatVector > *microSurfaceForcesResult = reader.getMicroSurfaceForces( );

    for ( auto it = microSurfaceForcesAnswer.begin( ); it != microSurfaceForcesAnswer.end( ); it++ ){

        auto r = microSurfaceForcesResult->find( it->first );

        if ( r == microSurfaceForcesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 30) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 31) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microExternalForcesAnswer
        =
            {
                { 15, { 4.100, 4.100, 12.740 } },
                { 31, { 4.100, 4.100, 17.060 } },
                { 13, { 7.460, 4.100, 12.740 } },
                { 26, { 7.460, 4.100, 17.060 } },
                { 53, { 4.100, 6.115, 12.740 } },
                { 21, { 4.100, 6.115, 17.060 } },
                { 37, { 7.460, 6.115, 12.740 } },
                { 48, { 7.460, 6.115, 17.060 } },
                {  5, { 10.820, 4.100, 12.740 } },
                { 10, { 10.820, 4.100, 17.060 } },
                {  3, { 10.820, 6.115, 12.740 } },
                {  4, { 10.820, 6.115, 17.060 } },
                { 32, { 7.460, 8.130, 12.740 } },
                { 33, { 7.460, 8.130, 17.060 } },
                { 34, { 10.820, 8.130, 12.740 } },
                { 28, { 10.820, 8.130, 17.060 } },
                { 25, { 4.100, 8.130, 12.740 } },
                { 50, { 4.100, 8.130, 17.060 } },
                { 43, { 4.100, 4.100, 21.380 } },
                { 27, { 7.460, 4.100, 21.380 } },
                {  1, { 4.100, 6.115, 21.380 } },
                {  7, { 7.460, 6.115, 21.380 } },
                { 30, { 10.820, 4.100, 21.380 } },
                { 16, { 10.820, 6.115, 21.380 } },
                { 22, { 7.460, 8.130, 21.380 } },
                {  2, { 10.820, 8.130, 21.380 } },
                { 46, { 4.100, 8.130, 21.380 } },
                { 24, { 4.100, 4.100, 4.100 } },
                { 39, { 4.100, 4.100, 8.420 } },
                { 40, { 7.460, 4.100, 4.100 } },
                { 57, { 7.460, 4.100, 8.420 } },
                { 44, { 4.100, 6.115, 4.100 } },
                { 58, { 4.100, 6.115, 8.420 } },
                { 29, { 7.460, 6.115, 4.100 } },
                { 59, { 7.460, 6.115, 8.420 } },
                { 11, { 10.820, 4.100, 4.100 } },
                {  0, { 10.820, 4.100, 8.420 } },
                { 20, { 10.820, 6.115, 4.100 } },
                { 60, { 10.820, 6.115, 8.420 } },
                { 47, { 7.460, 8.130, 4.100 } },
                { 49, { 7.460, 8.130, 8.420 } },
                { 17, { 10.820, 8.130, 4.100 } },
                { 38, { 10.820, 8.130, 8.420 } },
                { 14, { 4.100, 8.130, 4.100 } },
                { 55, { 4.100, 8.130, 8.420 } },
            };

    const std::unordered_map< uIntType, floatVector > *microExternalForcesResult = reader.getMicroExternalForces( );

    for ( auto it = microExternalForcesAnswer.begin( ); it != microExternalForcesAnswer.end( ); it++ ){

        auto r = microExternalForcesResult->find( it->first );

        if ( r == microExternalForcesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 32) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 33) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microVelocitiesAnswer
        =
            {
               { 15, { 5.000, -8.848, 11.950 } },
               { 31, { 5.000, -8.848, 13.025 } },
               { 13, { 6.250, -8.848, 11.950 } },
               { 26, { 6.250, -8.848, 13.025 } },
               { 53, { 5.000, -10.418, 11.950 } },
               { 21, { 5.000, -10.418, 13.025 } },
               { 37, { 6.250, -10.418, 11.950 } },
               { 48, { 6.250, -10.418, 13.025 } },
               {  5, { 7.500, -8.848, 11.950 } },
               { 10, { 7.500, -8.848, 13.025 } },
               {  3, { 7.500, -10.418, 11.950 } },
               {  4, { 7.500, -10.418, 13.025 } },
               { 32, { 6.250, -11.988, 11.950 } },
               { 33, { 6.250, -11.988, 13.025 } },
               { 34, { 7.500, -11.988, 11.950 } },
               { 28, { 7.500, -11.988, 13.025 } },
               { 25, { 5.000, -11.988, 11.950 } },
               { 50, { 5.000, -11.988, 13.025 } },
               { 43, { 5.000, -8.848, 14.100 } },
               { 27, { 6.250, -8.848, 14.100 } },
               {  1, { 5.000, -10.418, 14.100 } },
               {  7, { 6.250, -10.418, 14.100 } },
               { 30, { 7.500, -8.848, 14.100 } },
               { 16, { 7.500, -10.418, 14.100 } },
               { 22, { 6.250, -11.988, 14.100 } },
               {  2, { 7.500, -11.988, 14.100 } },
               { 46, { 5.000, -11.988, 14.100 } },
               { 24, { 5.000, -8.848, 9.800 } },
               { 39, { 5.000, -8.848, 10.875 } },
               { 40, { 6.250, -8.848, 9.800 } },
               { 57, { 6.250, -8.848, 10.875 } },
               { 44, { 5.000, -10.418, 9.800 } },
               { 58, { 5.000, -10.418, 10.875 } },
               { 29, { 6.250, -10.418, 9.800 } },
               { 59, { 6.250, -10.418, 10.875 } },
               { 11, { 7.500, -8.848, 9.800 } },
               {  0, { 7.500, -8.848, 10.875 } },
               { 20, { 7.500, -10.418, 9.800 } },
               { 60, { 7.500, -10.418, 10.875 } },
               { 47, { 6.250, -11.988, 9.800 } },
               { 49, { 6.250, -11.988, 10.875 } },
               { 17, { 7.500, -11.988, 9.800 } },
               { 38, { 7.500, -11.988, 10.875 } },
               { 14, { 5.000, -11.988, 9.800 } },
               { 55, { 5.000, -11.988, 10.875 } }
            };

    const std::unordered_map< uIntType, floatVector > *microVelocitiesResult = reader.getMicroVelocities( );

    for ( auto it = microVelocitiesAnswer.begin( ); it != microVelocitiesAnswer.end( ); it++ ){

        auto r = microVelocitiesResult->find( it->first );

        if ( r == microVelocitiesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 34) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 35) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microAccelerationsAnswer
        =
            {
               { 15, { 5.7765, 5.9930, 1.1000 } },
               { 31, { 5.7765, 5.9930, 2.1500 } },
               { 13, { 6.1315, 5.9930, 1.1000 } },
               { 26, { 6.1315, 5.9930, 2.1500 } },
               { 53, { 5.7765, 6.4080, 1.1000 } },
               { 21, { 5.7765, 6.4080, 2.1500 } },
               { 37, { 6.1315, 6.4080, 1.1000 } },
               { 48, { 6.1315, 6.4080, 2.1500 } },
               {  5, { 6.4865, 5.9930, 1.1000 } },
               { 10, { 6.4865, 5.9930, 2.1500 } },
               {  3, { 6.4865, 6.4080, 1.1000 } },
               {  4, { 6.4865, 6.4080, 2.1500 } },
               { 32, { 6.1315, 6.8230, 1.1000 } },
               { 33, { 6.1315, 6.8230, 2.1500 } },
               { 34, { 6.4865, 6.8230, 1.1000 } },
               { 28, { 6.4865, 6.8230, 2.1500 } },
               { 25, { 5.7765, 6.8230, 1.1000 } },
               { 50, { 5.7765, 6.8230, 2.1500 } },
               { 43, { 5.7765, 5.9930, 3.2000 } },
               { 27, { 6.1315, 5.9930, 3.2000 } },
               {  1, { 5.7765, 6.4080, 3.2000 } },
               {  7, { 6.1315, 6.4080, 3.2000 } },
               { 30, { 6.4865, 5.9930, 3.2000 } },
               { 16, { 6.4865, 6.4080, 3.2000 } },
               { 22, { 6.1315, 6.8230, 3.2000 } },
               {  2, { 6.4865, 6.8230, 3.2000 } },
               { 46, { 5.7765, 6.8230, 3.2000 } },
               { 24, { 5.7765, 5.9930, -1.0000 } },
               { 39, { 5.7765, 5.9930, 0.0500 } },
               { 40, { 6.1315, 5.9930, -1.0000 } },
               { 57, { 6.1315, 5.9930, 0.0500 } },
               { 44, { 5.7765, 6.4080, -1.0000 } },
               { 58, { 5.7765, 6.4080, 0.0500 } },
               { 29, { 6.1315, 6.4080, -1.0000 } },
               { 59, { 6.1315, 6.4080, 0.0500 } },
               { 11, { 6.4865, 5.9930, -1.0000 } },
               {  0, { 6.4865, 5.9930, 0.0500 } },
               { 20, { 6.4865, 6.4080, -1.0000 } },
               { 60, { 6.4865, 6.4080, 0.0500 } },
               { 47, { 6.1315, 6.8230, -1.0000 } },
               { 49, { 6.1315, 6.8230, 0.0500 } },
               { 17, { 6.4865, 6.8230, -1.0000 } },
               { 38, { 6.4865, 6.8230, 0.0500 } },
               { 14, { 5.7765, 6.8230, -1.0000 } },
               { 55, { 5.7765, 6.8230, 0.0500 } },
           };

    const std::unordered_map< uIntType, floatVector > *microAccelerationsResult = reader.getMicroAccelerations( );

    for ( auto it = microAccelerationsAnswer.begin( ); it != microAccelerationsAnswer.end( ); it++ ){

        auto r = microAccelerationsResult->find( it->first );

        if ( r == microAccelerationsResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 36) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 37) & False\n";
            return 1;

        }

    }

    const floatType previousTimeAnswer = 0.;
    const floatType *previousTimeResult = reader.getPreviousMicroTime( );

    if ( !vectorTools::fuzzyEquals( previousTimeAnswer, *previousTimeResult ) ){

        results << "test_initializeIncrement_Arlequin (test 38) & False\n";
        return 1;

    }

    const std::unordered_map< uIntType, floatVector > previousMicroDisplacementAnswer
        =
            {
                { 15, { -8.000, -8.000, -11.200 } },
                { 31, { -8.000, -8.000, -12.800 } },
                { 13, { -6.000, -8.000, -11.200 } },
                { 26, { -6.000, -8.000, -12.800 } },
                { 53, { -8.000, -4.800, -11.200 } },
                { 21, { -8.000, -4.800, -12.800 } },
                { 37, { -6.000, -4.800, -11.200 } },
                { 48, { -6.000, -4.800, -12.800 } },
                {  5, { -4.000, -8.000, -11.200 } },
                { 10, { -4.000, -8.000, -12.800 } },
                {  3, { -4.000, -4.800, -11.200 } },
                {  4, { -4.000, -4.800, -12.800 } },
                { 32, { -6.000, -1.600, -11.200 } },
                { 33, { -6.000, -1.600, -12.800 } },
                { 34, { -4.000, -1.600, -11.200 } },
                { 28, { -4.000, -1.600, -12.800 } },
                { 25, { -8.000, -1.600, -11.200 } },
                { 50, { -8.000, -1.600, -12.800 } },
                { 43, { -8.000, -8.000, -14.400 } },
                { 27, { -6.000, -8.000, -14.400 } },
                {  1, { -8.000, -4.800, -14.400 } },
                {  7, { -6.000, -4.800, -14.400 } },
                { 30, { -4.000, -8.000, -14.400 } },
                { 16, { -4.000, -4.800, -14.400 } },
                { 22, { -6.000, -1.600, -14.400 } },
                {  2, { -4.000, -1.600, -14.400 } },
                { 46, { -8.000, -1.600, -14.400 } },
                { 24, { -8.000, -8.000, -8.000 } },
                { 39, { -8.000, -8.000, -9.600 } },
                { 40, { -6.000, -8.000, -8.000 } },
                { 57, { -6.000, -8.000, -9.600 } },
                { 44, { -8.000, -4.800, -8.000 } },
                { 58, { -8.000, -4.800, -9.600 } },
                { 29, { -6.000, -4.800, -8.000 } },
                { 59, { -6.000, -4.800, -9.600 } },
                { 11, { -4.000, -8.000, -8.000 } },
                {  0, { -4.000, -8.000, -9.600 } },
                { 20, { -4.000, -4.800, -8.000 } },
                { 60, { -4.000, -4.800, -9.600 } },
                { 47, { -6.000, -1.600, -8.000 } },
                { 49, { -6.000, -1.600, -9.600 } },
                { 17, { -4.000, -1.600, -8.000 } },
                { 38, { -4.000, -1.600, -9.600 } },
                { 14, { -8.000, -1.600, -8.000 } },
                { 55, { -8.000, -1.600, -9.600 } },
            };

    const std::unordered_map< uIntType, floatVector > *previousMicroDisplacementResult = reader.getPreviousMicroDisplacements( );

    for ( auto it = previousMicroDisplacementAnswer.begin( ); it != previousMicroDisplacementAnswer.end( ); it++ ){

        auto r = previousMicroDisplacementResult->find( it->first );

        if ( r == previousMicroDisplacementResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 39) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 40) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > previousMicroVelocitiesAnswer
        =
            {
               { 15, { 1.200, 1.200, 3.350 } },
               { 31, { 1.200, 1.200, 4.425 } },
               { 13, { 2.450, 1.200, 3.350 } },
               { 26, { 2.450, 1.200, 4.425 } },
               { 53, { 1.200, -0.370, 3.350 } },
               { 21, { 1.200, -0.370, 4.425 } },
               { 37, { 2.450, -0.370, 3.350 } },
               { 48, { 2.450, -0.370, 4.425 } },
               {  5, { 3.700, 1.200, 3.350 } },
               { 10, { 3.700, 1.200, 4.425 } },
               {  3, { 3.700, -0.370, 3.350 } },
               {  4, { 3.700, -0.370, 4.425 } },
               { 32, { 2.450, -1.940, 3.350 } },
               { 33, { 2.450, -1.940, 4.425 } },
               { 34, { 3.700, -1.940, 3.350 } },
               { 28, { 3.700, -1.940, 4.425 } },
               { 25, { 1.200, -1.940, 3.350 } },
               { 50, { 1.200, -1.940, 4.425 } },
               { 43, { 1.200, 1.200, 5.500 } },
               { 27, { 2.450, 1.200, 5.500 } },
               {  1, { 1.200, -0.370, 5.500 } },
               {  7, { 2.450, -0.370, 5.500 } },
               { 30, { 3.700, 1.200, 5.500 } },
               { 16, { 3.700, -0.370, 5.500 } },
               { 22, { 2.450, -1.940, 5.500 } },
               {  2, { 3.700, -1.940, 5.500 } },
               { 46, { 1.200, -1.940, 5.500 } },
               { 24, { 1.200, 1.200, 1.200 } },
               { 39, { 1.200, 1.200, 2.275 } },
               { 40, { 2.450, 1.200, 1.200 } },
               { 57, { 2.450, 1.200, 2.275 } },
               { 44, { 1.200, -0.370, 1.200 } },
               { 58, { 1.200, -0.370, 2.275 } },
               { 29, { 2.450, -0.370, 1.200 } },
               { 59, { 2.450, -0.370, 2.275 } },
               { 11, { 3.700, 1.200, 1.200 } },
               {  0, { 3.700, 1.200, 2.275 } },
               { 20, { 3.700, -0.370, 1.200 } },
               { 60, { 3.700, -0.370, 2.275 } },
               { 47, { 2.450, -1.940, 1.200 } },
               { 49, { 2.450, -1.940, 2.275 } },
               { 17, { 3.700, -1.940, 1.200 } },
               { 38, { 3.700, -1.940, 2.275 } },
               { 14, { 1.200, -1.940, 1.200 } },
               { 55, { 1.200, -1.940, 2.275 } },
           };

    const std::unordered_map< uIntType, floatVector > *previousMicroVelocitiesResult = reader.getPreviousMicroVelocities( );

    for ( auto it = previousMicroVelocitiesAnswer.begin( ); it != previousMicroVelocitiesAnswer.end( ); it++ ){

        auto r = previousMicroVelocitiesResult->find( it->first );

        if ( r == previousMicroVelocitiesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 41) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 42) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > previousMicroAccelerationsAnswer
        =
            {
                { 15, { 4.2500, 4.2500, 6.3500 } },
                { 31, { 4.2500, 4.2500, 7.4000 } },
                { 13, { 4.6050, 4.2500, 6.3500 } },
                { 26, { 4.6050, 4.2500, 7.4000 } },
                { 53, { 4.2500, 4.6650, 6.3500 } },
                { 21, { 4.2500, 4.6650, 7.4000 } },
                { 37, { 4.6050, 4.6650, 6.3500 } },
                { 48, { 4.6050, 4.6650, 7.4000 } },
                {  5, { 4.9600, 4.2500, 6.3500 } },
                { 10, { 4.9600, 4.2500, 7.4000 } },
                {  3, { 4.9600, 4.6650, 6.3500 } },
                {  4, { 4.9600, 4.6650, 7.4000 } },
                { 32, { 4.6050, 5.0800, 6.3500 } },
                { 33, { 4.6050, 5.0800, 7.4000 } },
                { 34, { 4.9600, 5.0800, 6.3500 } },
                { 28, { 4.9600, 5.0800, 7.4000 } },
                { 25, { 4.2500, 5.0800, 6.3500 } },
                { 50, { 4.2500, 5.0800, 7.4000 } },
                { 43, { 4.2500, 4.2500, 8.4500 } },
                { 27, { 4.6050, 4.2500, 8.4500 } },
                {  1, { 4.2500, 4.6650, 8.4500 } },
                {  7, { 4.6050, 4.6650, 8.4500 } },
                { 30, { 4.9600, 4.2500, 8.4500 } },
                { 16, { 4.9600, 4.6650, 8.4500 } },
                { 22, { 4.6050, 5.0800, 8.4500 } },
                {  2, { 4.9600, 5.0800, 8.4500 } },
                { 46, { 4.2500, 5.0800, 8.4500 } },
                { 24, { 4.2500, 4.2500, 4.2500 } },
                { 39, { 4.2500, 4.2500, 5.3000 } },
                { 40, { 4.6050, 4.2500, 4.2500 } },
                { 57, { 4.6050, 4.2500, 5.3000 } },
                { 44, { 4.2500, 4.6650, 4.2500 } },
                { 58, { 4.2500, 4.6650, 5.3000 } },
                { 29, { 4.6050, 4.6650, 4.2500 } },
                { 59, { 4.6050, 4.6650, 5.3000 } },
                { 11, { 4.9600, 4.2500, 4.2500 } },
                {  0, { 4.9600, 4.2500, 5.3000 } },
                { 20, { 4.9600, 4.6650, 4.2500 } },
                { 60, { 4.9600, 4.6650, 5.3000 } },
                { 47, { 4.6050, 5.0800, 4.2500 } },
                { 49, { 4.6050, 5.0800, 5.3000 } },
                { 17, { 4.9600, 5.0800, 4.2500 } },
                { 38, { 4.9600, 5.0800, 5.3000 } },
                { 14, { 4.2500, 5.0800, 4.2500 } },
                { 55, { 4.2500, 5.0800, 5.3000 } },
            };

    const std::unordered_map< uIntType, floatVector > *previousMicroAccelerationsResult = reader.getPreviousMicroAccelerations( );

    for ( auto it = previousMicroAccelerationsAnswer.begin( ); it != previousMicroAccelerationsAnswer.end( ); it++ ){

        auto r = previousMicroAccelerationsResult->find( it->first );

        if ( r == previousMicroAccelerationsResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 43) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 44) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microStressesAnswer
        =
            {
                { 15, { 3.090371, 3.765898, 2.979736, 0.679351, 4.038599, 2.932600, 2.403665, 3.732252, 0.398282 } },
                { 31, { 3.726374, 4.279899, 2.727951, 0.186259, 4.194622, 3.003618, 2.646383, 3.839766, -0.401559 } },
                { 13, { 3.606480, 3.522090, 2.946139, 1.672865, 4.876204, 2.621181, 3.337650, 4.342345, 0.722520 } },
                { 26, { 4.242483, 4.036091, 2.694354, 1.179773, 5.032227, 2.692198, 3.580368, 4.449859, -0.077322 } },
                { 53, { 2.629139, 3.951084, 3.104852, 0.011765, 3.878585, 3.875028, 3.382334, 3.573368, 0.176006 } },
                { 21, { 3.265143, 4.465084, 2.853068, -0.481326, 4.034608, 3.946046, 3.625052, 3.680882, -0.623836 } },
                { 37, { 3.145248, 3.707276, 3.071255, 1.005279, 4.716189, 3.563609, 4.316319, 4.183461, 0.500244 } },
                { 48, { 3.781252, 4.221276, 2.819471, 0.512188, 4.872212, 3.634626, 4.559037, 4.290975, -0.299598 } },
                {  5, { 4.122589, 3.278282, 2.912542, 2.666379, 5.713808, 2.309761, 4.271634, 4.952438, 1.046758 } },
                { 10, { 4.758592, 3.792283, 2.660757, 2.173287, 5.869831, 2.380779, 4.514352, 5.059952, 0.246916 } },
                {  3, { 3.661357, 3.463468, 3.037658, 1.998793, 5.553793, 3.252189, 5.250303, 4.793554, 0.824481 } },
                {  4, { 4.297360, 3.977468, 2.785874, 1.505701, 5.709816, 3.323207, 5.493021, 4.901068, 0.024640 } },
                { 32, { 2.684016, 3.892461, 3.196372, 0.337694, 4.556174, 4.506037, 5.294987, 4.024577, 0.277967 } },
                { 33, { 3.320020, 4.406462, 2.944587, -0.155398, 4.712197, 4.577054, 5.537706, 4.132091, -0.521874 } },
                { 34, { 3.200125, 3.648653, 3.162775, 1.331208, 5.393779, 4.194617, 6.228972, 4.634670, 0.602205 } },
                { 28, { 3.836129, 4.162654, 2.910990, 0.838116, 5.549802, 4.265635, 6.471690, 4.742184, -0.197637 } },
                { 25, { 2.167908, 4.136269, 3.229969, -0.655820, 3.718570, 4.817456, 4.361003, 3.414484, -0.046270 } },
                { 50, { 2.803911, 4.650270, 2.978184, -1.148912, 3.874593, 4.888474, 4.603721, 3.521998, -0.846112 } },
                { 43, { 4.362378, 4.793900, 2.476167, -0.306832, 4.350645, 3.074636, 2.889102, 3.947280, -1.201401 } },
                { 27, { 4.878487, 4.550092, 2.442570, 0.686682, 5.188250, 2.763216, 3.823086, 4.557373, -0.877164 } },
                {  1, { 3.901146, 4.979085, 2.601283, -0.974418, 4.190631, 4.017064, 3.867770, 3.788396, -1.423677 } },
                {  7, { 4.417255, 4.735277, 2.567686, 0.019096, 5.028235, 3.705644, 4.801755, 4.398490, -1.099440 } },
                { 30, { 5.394596, 4.306284, 2.408973, 1.680195, 6.025854, 2.451796, 4.757071, 5.167467, -0.552926 } },
                { 16, { 4.933364, 4.491469, 2.534089, 1.012610, 5.865839, 3.394224, 5.735740, 5.008583, -0.775202 } },
                { 22, { 3.956023, 4.920462, 2.692802, -0.648490, 4.868220, 4.648072, 5.780424, 4.239606, -1.321716 } },
                {  2, { 4.472132, 4.676654, 2.659205, 0.345024, 5.705825, 4.336652, 6.714408, 4.849699, -0.997479 } },
                { 46, { 3.439914, 5.164270, 2.726399, -1.642003, 4.030616, 4.959492, 4.846439, 3.629512, -1.645954 } },
                { 24, { 1.818364, 2.737897, 3.483305, 1.665534, 3.726553, 2.790565, 1.918229, 3.517223, 1.997966 } },
                { 39, { 2.454367, 3.251897, 3.231521, 1.172443, 3.882576, 2.861583, 2.160947, 3.624737, 1.198124 } },
                { 40, { 2.334473, 2.494089, 3.449708, 2.659048, 4.564158, 2.479145, 2.852213, 4.127316, 2.322204 } },
                { 57, { 2.970476, 3.008089, 3.197924, 2.165957, 4.720181, 2.550163, 3.094932, 4.234830, 1.522362 } },
                { 44, { 1.357132, 2.923082, 3.608422, 0.997949, 3.566539, 3.732993, 2.896898, 3.358339, 1.775690 } },
                { 58, { 1.993136, 3.437083, 3.356637, 0.504857, 3.722562, 3.804011, 3.139616, 3.465853, 0.975848 } },
                { 29, { 1.873241, 2.679274, 3.574825, 1.991463, 4.404143, 3.421573, 3.830882, 3.968432, 2.099927 } },
                { 59, { 2.509245, 3.193275, 3.323040, 1.498371, 4.560166, 3.492591, 4.073600, 4.075946, 1.300085 } },
                { 11, { 2.850582, 2.250281, 3.416111, 3.652562, 5.401762, 2.167726, 3.786198, 4.737409, 2.646441 } },
                {  0, { 3.486585, 2.764281, 3.164327, 3.159470, 5.557785, 2.238743, 4.028916, 4.844923, 1.846599 } },
                { 20, { 2.389350, 2.435466, 3.541228, 2.984976, 5.241747, 3.110154, 4.764867, 4.578525, 2.424165 } },
                { 60, { 3.025354, 2.949467, 3.289443, 2.491885, 5.397770, 3.181171, 5.007585, 4.686039, 1.624323 } },
                { 47, { 1.412010, 2.864460, 3.699941, 1.323877, 4.244128, 4.364001, 4.809551, 3.809548, 1.877651 } },
                { 49, { 2.048013, 3.378460, 3.448156, 0.830785, 4.400151, 4.435019, 5.052269, 3.917062, 1.077809 } },
                { 17, { 1.928119, 2.620652, 3.666344, 2.317391, 5.081733, 4.052582, 5.743536, 4.419641, 2.201888 } },
                { 38, { 2.564122, 3.134652, 3.414559, 1.824299, 5.237756, 4.123599, 5.986254, 4.527156, 1.402047 } },
                { 14, { 0.895901, 3.108268, 3.733538, 0.330363, 3.406524, 4.675421, 3.875566, 3.199455, 1.553413 } },
                { 55, { 1.531904, 3.622268, 3.481753, -0.162728, 3.562547, 4.746439, 4.118285, 3.306969, 0.753572 } },
            };

    const std::unordered_map< uIntType, floatVector > *microStressesResult = reader.getMicroStresses( );

    for ( auto it = microStressesAnswer.begin( ); it != microStressesAnswer.end( ); it++ ){

        auto r = microStressesResult->find( it->first );

        if ( r == microStressesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 45) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 46) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microInternalForcesAnswer
        =
             {
                { 15, { 2.851562, 4.231133, 1.866341 } },
                { 31, { 3.344066, 4.560423, 0.901855 } },
                { 13, { 3.793505, 3.271966, 2.243144 } },
                { 26, { 4.286010, 3.601256, 1.278659 } },
                { 53, { 2.646243, 5.060558, 1.343330 } },
                { 21, { 3.138748, 5.389848, 0.378845 } },
                { 37, { 3.588187, 4.101391, 1.720134 } },
                { 48, { 4.080691, 4.430681, 0.755648 } },
                {  5, { 4.735449, 2.312798, 2.619948 } },
                { 10, { 5.227953, 2.642088, 1.655463 } },
                {  3, { 4.530131, 3.142223, 2.096937 } },
                {  4, { 5.022635, 3.471513, 1.132452 } },
                { 32, { 3.382869, 4.930815, 1.197123 } },
                { 33, { 3.875373, 5.260105, 0.232638 } },
                { 34, { 4.324812, 3.971648, 1.573927 } },
                { 28, { 4.817317, 4.300938, 0.609442 } },
                { 25, { 2.440925, 5.889983, 0.820320 } },
                { 50, { 2.933429, 6.219273, -0.144166 } },
                { 43, { 3.836570, 4.889713, -0.062630 } },
                { 27, { 4.778514, 3.930546, 0.314174 } },
                {  1, { 3.631252, 5.719138, -0.585641 } },
                {  7, { 4.573196, 4.759971, -0.208837 } },
                { 30, { 5.720458, 2.971378, 0.690977 } },
                { 16, { 5.515139, 3.800803, 0.167967 } },
                { 22, { 4.367877, 5.589395, -0.731847 } },
                {  2, { 5.309821, 4.630228, -0.355044 } },
                { 46, { 3.425933, 6.548563, -1.108651 } },
                { 24, { 1.866553, 3.572553, 3.795311 } },
                { 39, { 2.359057, 3.901843, 2.830826 } },
                { 40, { 2.808497, 2.613386, 4.172115 } },
                { 57, { 3.301001, 2.942676, 3.207630 } },
                { 44, { 1.661235, 4.401978, 3.272301 } },
                { 58, { 2.153739, 4.731268, 2.307815 } },
                { 29, { 2.603179, 3.442811, 3.649104 } },
                { 59, { 3.095683, 3.772101, 2.684619 } },
                { 11, { 3.750441, 1.654218, 4.548919 } },
                {  0, { 4.242945, 1.983508, 3.584433 } },
                { 20, { 3.545122, 2.483643, 4.025908 } },
                { 60, { 4.037627, 2.812933, 3.061423 } },
                { 47, { 2.397860, 4.272235, 3.126094 } },
                { 49, { 2.890364, 4.601525, 2.161609 } },
                { 17, { 3.339804, 3.313068, 3.502898 } },
                { 38, { 3.832308, 3.642358, 2.538412 } },
                { 14, { 1.455917, 5.231403, 2.749290 } },
                { 55, { 1.948421, 5.560693, 1.784805 } },
            };

    const std::unordered_map< uIntType, floatVector > *microInternalForcesResult = reader.getMicroInternalForces( );

    for ( auto it = microInternalForcesAnswer.begin( ); it != microInternalForcesAnswer.end( ); it++ ){

        auto r = microInternalForcesResult->find( it->first );

        if ( r == microInternalForcesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 47) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 48) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > microInertialForcesAnswer
        =
            {
                { 15, { 0.915926, -0.776311, -3.061289 } },
                { 31, { 1.713875, 0.043499, -3.050079 } },
                { 13, { 0.277727, -0.625485, -3.510893 } },
                { 26, { 1.075675, 0.194325, -3.499683 } },
                { 53, { 0.331368, -0.491830, -3.876881 } },
                { 21, { 1.129317, 0.327980, -3.865671 } },
                { 37, { -0.306832, -0.341004, -4.326486 } },
                { 48, { 0.491117, 0.478807, -4.315275 } },
                {  5, { -0.360473, -0.474659, -3.960498 } },
                { 10, { 0.437475, 0.345151, -3.949287 } },
                {  3, { -0.945031, -0.190177, -4.776090 } },
                {  4, { -0.147083, 0.629633, -4.764879 } },
                { 32, { -0.891390, -0.056522, -5.142078 } },
                { 33, { -0.093442, 0.763288, -5.130867 } },
                { 34, { -1.529590, 0.094304, -5.591682 } },
                { 28, { -0.731641, 0.914114, -5.580471 } },
                { 25, { -0.253190, -0.207348, -4.692474 } },
                { 50, { 0.544758, 0.612462, -4.681263 } },
                { 43, { 2.511823, 0.863309, -3.038868 } },
                { 27, { 1.873624, 1.014136, -3.488472 } },
                {  1, { 1.927265, 1.147791, -3.854460 } },
                {  7, { 1.289065, 1.298617, -4.304064 } },
                { 30, { 1.235424, 1.164962, -3.938076 } },
                { 16, { 0.650865, 1.449443, -4.753668 } },
                { 22, { 0.704507, 1.583098, -5.119656 } },
                {  2, { 0.066307, 1.733925, -5.569261 } },
                { 46, { 1.342707, 1.432272, -4.670052 } },
                { 24, { -0.679970, -2.415932, -3.083711 } },
                { 39, { 0.117978, -1.596122, -3.072500 } },
                { 40, { -1.318170, -2.265106, -3.533315 } },
                { 57, { -0.520222, -1.445295, -3.522104 } },
                { 44, { -1.264529, -2.131450, -3.899303 } },
                { 58, { -0.466580, -1.311640, -3.888092 } },
                { 29, { -1.902729, -1.980624, -4.348907 } },
                { 59, { -1.104780, -1.160814, -4.337696 } },
                { 11, { -1.956370, -2.114279, -3.982919 } },
                {  0, { -1.158422, -1.294469, -3.971708 } },
                { 20, { -2.540928, -1.829798, -4.798511 } },
                { 60, { -1.742980, -1.009988, -4.787300 } },
                { 47, { -2.487287, -1.696143, -5.164499 } },
                { 49, { -1.689339, -0.876332, -5.153289 } },
                { 17, { -3.125487, -1.545316, -5.614103 } },
                { 38, { -2.327538, -0.725506, -5.602893 } },
                { 14, { -1.849087, -1.846969, -4.714895 } },
                { 55, { -1.051139, -1.027159, -4.703684 } },
            };

    const std::unordered_map< uIntType, floatVector > *microInertialForcesResult = reader.getMicroInertialForces( );

    for ( auto it = microInertialForcesAnswer.begin( ); it != microInertialForcesAnswer.end( ); it++ ){

        auto r = microInertialForcesResult->find( it->first );

        if ( r == microInertialForcesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 49) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 50) & False\n";
            return 1;

        }

    }

    const floatType macroTimeAnswer = 1.;
    const floatType* macroTimeResult = reader.getMacroTime( );

    if ( !vectorTools::fuzzyEquals( macroTimeAnswer, *macroTimeResult ) ){

        results << "test_initializeIncrement_Arlequin (test 51) & False\n";
        return 1;

    }

    const std::unordered_map< uIntType, floatVector > macroDisplacementsAnswer
        =
            {
                {  5, { 0.641275, 0.232390, -2.327989 } },
                {  9, { -0.277488, 1.864821, -0.959118 } },
                {  8, { -1.872777, 2.331052, -2.562691 } },
                { 11, { -0.954015, 0.698621, -3.931561 } },
                {  3, { 0.863789, 1.140577, -2.616417 } },
                {  1, { -0.054974, 2.773008, -1.247547 } },
                {  6, { -1.650263, 3.239239, -2.851120 } },
                { 15, { -0.731501, 1.606808, -4.219990 } },
                { 12, { 1.086303, 2.048764, -2.904846 } },
                {  2, { 0.167540, 3.681195, -1.535975 } },
                { 13, { -1.427749, 4.147426, -3.139548 } },
                { 14, { -0.508987, 2.514995, -4.508419 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroDisplacementsResult = reader.getMacroDisplacements( );

    for ( auto it = macroDisplacementsAnswer.begin( ); it != macroDisplacementsAnswer.end( ); it++ ){

        auto r = macroDisplacementsResult->find( it->first );

        if ( r == macroDisplacementsResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 51) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 52) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroDispDOFVectorAnswer
        =
            {
                {  5, { 0.641275, 0.232390, -2.327989, 2.476106, 3.649307, 0.380024, -0.602181, -0.098268, 1.214942, -1.346951, 0.875060, 0.810153 } },
                {  9, { -0.277488, 1.864821, -0.959118, 3.747254, 2.993086, 1.602849, -0.801713, 0.385446, 1.268479, -0.333959, 0.834465, 0.773507 } },
                {  8, { -1.872777, 2.331052, -2.562691, 5.125425, 4.185234, 0.239589, -2.210062, -1.562753, 0.311109, -2.017245, 2.360391, -1.174211 } },
                { 11, { -0.954015, 0.698621, -3.931561, 3.854278, 4.841455, -0.983235, -2.010531, -2.046468, 0.257571, -3.030237, 2.400986, -1.137565 } },
                {  3, { 0.863789, 1.140577, -2.616417, 3.550081, 5.494371, 0.536456, 0.929143, -0.090355, 2.448663, -3.175893, 0.545918, -0.235911 } },
                {  1, { -0.054974, 2.773008, -1.247547, 4.821228, 4.838150, 1.759280, 0.729611, 0.393360, 2.502200, -2.162901, 0.505323, -0.272557 } },
                {  6, { -1.650263, 3.239239, -2.851120, 6.199400, 6.030298, 0.396021, -0.678739, -1.554840, 1.544830, -3.846187, 2.031250, -2.220275 } },
                { 15, { -0.731501, 1.606808, -4.219990, 4.928252, 6.686519, -0.826804, -0.479207, -2.038554, 1.491293, -4.859179, 2.071844, -2.183628 } },
                { 12, { 1.086303, 2.048764, -2.904846, 4.624056, 7.339434, 0.692887, 2.460467, -0.082442, 3.682384, -5.004835, 0.216776, -1.281975 } },
                {  2, { 0.167540, 3.681195, -1.535975, 5.895203, 6.683213, 1.915712, 2.260935, 0.401273, 3.735921, -3.991843, 0.176182, -1.318621 } },
                { 13, { -1.427749, 4.147426, -3.139548, 7.273375, 7.875361, 0.552453, 0.852585, -1.546927, 2.778551, -5.675130, 1.702108, -3.266339 } },
                { 14, { -0.508987, 2.514995, -4.508419, 6.002227, 8.531582, -0.670372, 1.052117, -2.030641, 2.725014, -6.688121, 1.742702, -3.229692 } },
            };
    const std::unordered_map< uIntType, floatVector > *macroDispDOFVectorResult = reader.getMacroDispDOFVector( );

    for ( auto it = macroDispDOFVectorAnswer.begin( ); it != macroDispDOFVectorAnswer.end( ); it++ ){

        auto r = macroDispDOFVectorResult->find( it->first );

        if ( r == macroDispDOFVectorResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 53) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 54) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroVelocitiesAnswer
        =
            {
                {  5, { -2.346964, -2.328088, 2.337123, 1.674508, 1.402881, -1.276984, 2.105710, 1.195199, 0.817334, 1.284606, -0.465939, -1.366498 } },
                {  9, { -0.548438, -2.321249, 3.137203, -0.038677, 2.051591, -1.215603, 0.532880, 0.749185, 1.584563, 0.398241, -2.435270, -0.742095 } },
                {  8, { 1.144646, -3.077963, 3.498864, -1.812710, 3.863297, -3.088573, -0.357800, 1.094312, 2.608025, 0.763815, -2.637127, -0.243639 } },
                { 11, { -0.653880, -3.084803, 2.698784, -0.099526, 3.214588, -3.149954, 1.215030, 1.540326, 1.840796, 1.650181, -0.667796, -0.868043 } },
                {  3, { -1.773946, -3.882980, 2.418979, 1.029116, 1.811975, -0.760381, 3.600941, 2.431781, 1.015779, 0.656494, 0.547274, -2.449650 } },
                {  1, { 0.024580, -3.876140, 3.219060, -0.684069, 2.460685, -0.699000, 2.028111, 1.985767, 1.783008, -0.229872, -1.422057, -1.825247 } },
                {  6, { 1.717664, -4.632855, 3.580720, -2.458102, 4.272392, -2.571970, 1.137432, 2.330893, 2.806471, 0.135703, -1.623914, -1.326792 } },
                { 15, { -0.080862, -4.639694, 2.780640, -0.744918, 3.623682, -2.633351, 2.710262, 2.776907, 2.039242, 1.022068, 0.345417, -1.951195 } },
                { 12, { -1.200928, -5.437871, 2.500835, 0.383724, 2.221070, -0.243779, 5.096172, 3.668363, 1.214225, 0.028381, 1.560487, -3.532803 } },
                {  2, { 0.597598, -5.431032, 3.300916, -1.329461, 2.869779, -0.182398, 3.523342, 3.222349, 1.981454, -0.857984, -0.408844, -2.908399 } },
                { 13, { 2.290682, -6.187747, 3.662577, -3.103494, 4.681486, -2.055368, 2.632663, 3.567475, 3.004916, -0.492410, -0.610701, -2.409944 } },
                { 14, { 0.492156, -6.194586, 2.862496, -1.390310, 4.032776, -2.116749, 4.205493, 4.013489, 2.237687, 0.393955, 1.358630, -3.034347 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroVelocitiesResult = reader.getMacroVelocities( );

    for ( auto it = macroVelocitiesAnswer.begin( ); it != macroVelocitiesAnswer.end( ); it++ ){

        auto r = macroVelocitiesResult->find( it->first );

        if ( r == macroVelocitiesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 55) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 56) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroAccelerationsAnswer
        =
            {
               {  5, { 0.427196, 0.928039, 1.624912, -2.356312, 1.564029, -1.719566, 0.672795, -0.997762, -0.440505, 1.065985, 0.969844, 2.956951 } },
               {  9, { 1.700591, 2.076310, 2.157836, -1.571605, 3.403657, -2.637085, 2.643157, 0.496112, -0.696747, 2.936684, -0.903955, 4.777131 } },
               {  8, { 2.377492, 3.010728, 1.465971, -0.839560, 1.625823, -0.674050, 4.129535, -0.565470, -1.308817, 3.939671, 1.079660, 5.432814 } },
               { 11, { 1.104097, 1.862456, 0.933047, -1.624267, -0.213804, 0.243469, 2.159173, -2.059344, -1.052576, 2.068972, 2.953459, 3.612634 } },
               {  3, { -0.120091, 1.778900, 3.245842, -3.033673, 0.251196, -1.879919, 1.501775, -1.509271, 0.832536, 0.079488, 0.817082, 4.057459 } },
               {  1, { 1.153304, 2.927172, 3.778766, -2.248965, 2.090823, -2.797438, 3.472137, -0.015397, 0.576294, 1.950187, -1.056717, 5.877639 } },
               {  6, { 1.830205, 3.861590, 3.086901, -1.516920, 0.312989, -0.834402, 4.958514, -1.076979, -0.035777, 2.953174, 0.926898, 6.533322 } },
               { 15, { 0.556810, 2.713318, 2.553977, -2.301628, -1.526638, 0.083117, 2.988152, -2.570853, 0.220465, 1.082475, 2.800697, 4.713142 } },
               { 12, { -0.667378, 2.629762, 4.866772, -3.711033, -1.061638, -2.040271, 2.330754, -2.020779, 2.105576, -0.907009, 0.664319, 5.157967 } },
               {  2, { 0.606017, 3.778034, 5.399696, -2.926326, 0.777989, -2.957790, 4.301116, -0.526906, 1.849334, 0.963690, -1.209479, 6.978147 } },
               { 13, { 1.282918, 4.712451, 4.707831, -2.194281, -0.999844, -0.994755, 5.787494, -1.588488, 1.237264, 1.966677, 0.774135, 7.633830 } },
               { 14, { 0.009523, 3.564180, 4.174907, -2.978988, -2.839471, -0.077236, 3.817132, -3.082362, 1.493506, 0.095978, 2.647934, 5.813650 } },
           };

    const std::unordered_map< uIntType, floatVector > *macroAccelerationsResult = reader.getMacroAccelerations( );

    for ( auto it = macroAccelerationsAnswer.begin( ); it != macroAccelerationsAnswer.end( ); it++ ){

        auto r = macroAccelerationsResult->find( it->first );

        if ( r == macroAccelerationsResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 57) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 58) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > previousMacroDispDOFVectorAnswer
        =
            {
                {  5, { 0.028316, 1.669368, -0.675497, 1.237413, 1.695543, -0.472972, 0.302567, 0.847229, 1.313857, 0.027499, 0.297804, -0.165902 } },
                {  9, { -0.890446, 3.301799, 0.693373, 2.508561, 1.039322, 0.749852, 0.103035, 1.330944, 1.367395, 1.040491, 0.257210, -0.202549 } },
                {  8, { -2.485736, 3.768030, -0.910199, 3.886732, 2.231470, -0.613407, -1.305314, -0.617256, 0.410025, -0.642795, 1.783136, -2.150267 } },
                { 11, { -1.566974, 2.135599, -2.279070, 2.615585, 2.887691, -1.836232, -1.105782, -1.100970, 0.356487, -1.655787, 1.823730, -2.113620 } },
                {  3, { 0.250830, 2.577555, -0.963926, 2.311388, 3.540607, -0.316540, 1.833891, 0.855143, 2.547578, -1.801443, -0.031338, -1.211966 } },
                {  1, { -0.667932, 4.209986, 0.404945, 3.582536, 2.884385, 0.906284, 1.634359, 1.338857, 2.601116, -0.788451, -0.071932, -1.248613 } },
                {  6, { -2.263222, 4.676217, -1.198628, 4.960707, 4.076534, -0.456975, 0.226010, -0.609342, 1.643746, -2.471737, 1.453994, -3.196331 } },
                { 15, { -1.344460, 3.043786, -2.567498, 3.689560, 4.732755, -1.679800, 0.425541, -1.093057, 1.590208, -3.484729, 1.494589, -3.159684 } },
                { 12, { 0.473344, 3.485742, -1.252354, 3.385363, 5.385670, -0.160109, 3.365215, 0.863056, 3.781300, -3.630385, -0.360479, -2.258030 } },
                {  2, { -0.445418, 5.118173, 0.116516, 4.656510, 4.729449, 1.062716, 3.165683, 1.346770, 3.834837, -2.617393, -0.401074, -2.294677 } },
                { 13, { -2.040708, 5.584404, -1.487057, 6.034682, 5.921597, -0.300544, 1.757333, -0.601429, 2.877467, -4.300680, 1.124852, -4.242395 } },
                { 14, { -1.121946, 3.951973, -2.855927, 4.763534, 6.577818, -1.523368, 1.956865, -1.085144, 2.823929, -5.313671, 1.165447, -4.205748 } },
            };

    const std::unordered_map< uIntType, floatVector > *previousMacroDispDOFVectorResult = reader.getPreviousMacroDispDOFVector( );

    for ( auto it = previousMacroDispDOFVectorAnswer.begin( ); it != previousMacroDispDOFVectorAnswer.end( ); it++ ){

        auto r = previousMacroDispDOFVectorResult->find( it->first );

        if ( r == previousMacroDispDOFVectorResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 59) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 60) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > previousMacroVelocitiesAnswer
        =
             {
                {  5, { -1.501765, -1.483894, 0.684890, 0.330246, -0.559686, -0.727351, 1.240470, 1.138232, 0.275631, 1.151980, 0.135063, -0.239197 } },
                {  9, { 0.296761, -1.477054, 1.484970, -1.382938, 0.089023, -0.665970, -0.332360, 0.692218, 1.042860, 0.265615, -1.834268, 0.385206 } },
                {  8, { 1.989845, -2.233769, 1.846631, -3.156971, 1.900730, -2.538940, -1.223039, 1.037345, 2.066323, 0.631189, -2.036125, 0.883662 } },
                { 11, { 0.191319, -2.240608, 1.046550, -1.443787, 1.252020, -2.600321, 0.349790, 1.483359, 1.299094, 1.517554, -0.066794, 0.259258 } },
                {  3, { -0.928747, -3.038786, 0.766746, -0.315146, -0.150592, -0.210748, 2.735701, 2.374814, 0.474077, 0.523867, 1.148276, -1.322349 } },
                {  1, { 0.869779, -3.031946, 1.566826, -2.028330, 0.498117, -0.149367, 1.162871, 1.928800, 1.241306, -0.362498, -0.821055, -0.697946 } },
                {  6, { 2.562863, -3.788661, 1.928487, -3.802363, 2.309824, -2.022337, 0.272192, 2.273927, 2.264768, 0.003076, -1.022912, -0.199490 } },
                { 15, { 0.764337, -3.795500, 1.128407, -2.089179, 1.661114, -2.083718, 1.845022, 2.719941, 1.497539, 0.889442, 0.946419, -0.823894 } },
                { 12, { -0.355729, -4.593677, 0.848602, -0.960538, 0.258502, 0.305854, 4.230932, 3.611396, 0.672522, -0.104245, 2.161489, -2.405502 } },
                {  2, { 1.442797, -4.586838, 1.648683, -2.673722, 0.907211, 0.367235, 2.658103, 3.165382, 1.439751, -0.990611, 0.192157, -1.781098 } },
                { 13, { 3.135881, -5.343553, 2.010343, -4.447755, 2.718918, -1.505734, 1.767423, 3.510508, 2.463214, -0.625036, -0.009700, -1.282643 } },
                { 14, { 1.337355, -5.350392, 1.210263, -2.734571, 2.070208, -1.567116, 3.340253, 3.956522, 1.695985, 0.261329, 1.959631, -1.907046 } },
            };

    const std::unordered_map< uIntType, floatVector > *previousMacroVelocitiesResult = reader.getPreviousMacroVelocities( );

    for ( auto it = previousMacroVelocitiesAnswer.begin( ); it != previousMacroVelocitiesAnswer.end( ); it++ ){

        auto r = previousMacroVelocitiesResult->find( it->first );

        if ( r == previousMacroVelocitiesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 61) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 62) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > previousMacroAccelerationsAnswer
        =
            {
               {  5, { 1.590337, -0.180779, 0.211887, -0.897186, 0.426442, -1.677882, -1.198663, 0.035661, 1.004663, 0.476809, -0.315302, 1.383007 } },
               {  9, { 2.863731, 0.967492, 0.744811, -0.112479, 2.266069, -2.595401, 0.771698, 1.529535, 0.748421, 2.347508, -2.189100, 3.203187 } },
               {  8, { 3.540632, 1.901910, 0.052946, 0.619566, 0.488235, -0.632365, 2.258076, 0.467953, 0.136351, 3.350494, -0.205486, 3.858870 } },
               { 11, { 2.267237, 0.753638, -0.479978, -0.165141, -1.351392, 0.285154, 0.287714, -1.025921, 0.392593, 1.479796, 1.668313, 2.038690 } },
               {  3, { 1.043050, 0.670082, 1.832817, -1.574547, -0.886392, -1.838234, -0.369684, -0.475848, 2.277704, -0.509688, -0.468064, 2.483515 } },
               {  1, { 2.316444, 1.818354, 2.365741, -0.789840, 0.953235, -2.755753, 1.600678, 1.018026, 2.021462, 1.361011, -2.341863, 4.303695 } },
               {  6, { 2.993345, 2.752772, 1.673876, -0.057795, -0.824598, -0.792718, 3.087056, -0.043556, 1.409392, 2.363997, -0.358248, 4.959378 } },
               { 15, { 1.719950, 1.604500, 1.140952, -0.842502, -2.664225, 0.124801, 1.116694, -1.537430, 1.665633, 0.493299, 1.515551, 3.139198 } },
               { 12, { 0.495763, 1.520944, 3.453747, -2.251908, -2.199226, -1.998587, 0.459296, -0.987357, 3.550745, -1.496185, -0.620827, 3.584022 } },
               {  2, { 1.769158, 2.669216, 3.986671, -1.467200, -0.359598, -2.916106, 2.429658, 0.506517, 3.294503, 0.374514, -2.494625, 5.404202 } },
               { 13, { 2.446058, 3.603633, 3.294806, -0.735155, -2.137432, -0.953071, 3.916035, -0.555065, 2.682432, 1.377500, -0.511011, 6.059886 } },
               { 14, { 1.172663, 2.455361, 2.761882, -1.519863, -3.977059, -0.035552, 1.945673, -2.048939, 2.938674, -0.493198, 1.362788, 4.239706 } },
           };

    const std::unordered_map< uIntType, floatVector > *previousMacroAccelerationsResult = reader.getPreviousMacroAccelerations( );

    for ( auto it = previousMacroAccelerationsAnswer.begin( ); it != previousMacroAccelerationsAnswer.end( ); it++ ){

        auto r = previousMacroAccelerationsResult->find( it->first );

        if ( r == previousMacroAccelerationsResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 63) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 64) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroInternalForcesAnswer
        =
            {
                {  5, { -3.123250, 0.977401, 2.054240, 1.757330, -0.474837, 1.152554, 0.177148, -0.110803, -0.822029, -0.726168, -0.869646, 0.163025 } },
                {  9, { -1.936939, 2.067285, 3.129559, 1.307129, 1.073546, 1.568170, 0.669365, 0.549069, -1.706006, 0.098907, -0.200723, 0.426380 } },
                {  8, { -2.811405, 2.316517, 3.559565, 0.391797, 2.169330, 0.888284, -0.794757, -1.286124, 0.075576, 0.364483, 1.119236, -0.465277 } },
                { 11, { -3.997715, 1.226633, 2.484245, 0.841998, 0.620946, 0.472668, -1.286974, -1.945997, 0.959553, -0.460592, 0.450314, -0.728632 } },
                {  3, { -3.134236, 2.438781, 3.072693, 3.369779, -1.079855, 0.553564, -0.716138, 0.891244, -2.584717, -2.198251, 0.401625, 0.010912 } },
                {  1, { -1.947925, 3.528665, 4.148013, 2.919578, 0.468528, 0.969180, -0.223921, 1.551116, -3.468694, -1.373176, 1.070548, 0.274268 } },
                {  6, { -2.822390, 3.777897, 4.578018, 2.004247, 1.564312, 0.289294, -1.688043, -0.284077, -1.687111, -1.107600, 2.390507, -0.617389 } },
                { 15, { -4.008701, 2.688013, 3.502699, 2.454448, 0.015928, -0.126322, -2.180260, -0.943950, -0.803134, -1.932675, 1.721585, -0.880745 } },
                { 12, { -3.145221, 3.900161, 4.091147, 4.982229, -1.684873, -0.045426, -1.609424, 1.893290, -4.347404, -3.670334, 1.672896, -0.141200 } },
                {  2, { -1.958911, 4.990045, 5.166466, 4.532028, -0.136490, 0.370190, -1.117208, 2.553163, -5.231381, -2.845259, 2.341819, 0.122155 } },
                { 13, { -2.833376, 5.239278, 5.596472, 3.616696, 0.959294, -0.309696, -2.581330, 0.717969, -3.449799, -2.579683, 3.661778, -0.769502 } },
                { 14, { -4.019687, 4.149394, 4.521152, 4.066897, -0.589089, -0.725312, -3.073546, 0.058097, -2.565822, -3.404757, 2.992856, -1.032857 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroInternalForcesResult = reader.getMacroInternalForces( );

    for ( auto it = macroInternalForcesAnswer.begin( ); it != macroInternalForcesAnswer.end( ); it++ ){

        auto r = macroInternalForcesResult->find( it->first );

        if ( r == macroInternalForcesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 65) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 66) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroInertialForcesAnswer
        =
            {
                {  5, { 1.964831, 1.278634, -1.125705, -2.353362, 0.113154, -1.589520, 0.935279, 1.013984, 2.260416, 3.108513, -2.500627, 0.826868 } },
                {  9, { 2.191268, 0.761777, -1.956907, -1.913616, -0.753844, -2.572974, -0.979177, 1.534873, 3.600554, 2.755534, -2.233762, 2.617534 } },
                {  8, { 1.628217, 1.471129, -1.883881, -1.853679, -2.038739, -1.016306, -1.413380, 0.733672, 2.539391, 4.310301, -0.946690, 2.543249 } },
                { 11, { 1.401781, 1.987986, -1.052679, -2.293425, -1.171740, -0.032851, 0.501075, 0.212783, 1.199253, 4.663281, -1.213556, 0.752583 } },
                {  3, { 1.065332, 3.110252, -0.267265, -0.927058, -1.756161, -2.313741, -0.023437, -0.393319, 1.532658, 2.038427, -4.150676, 2.771371 } },
                {  1, { 1.291768, 2.593395, -1.098467, -0.487312, -2.623159, -3.297196, -1.937893, 0.127570, 2.872797, 1.685447, -3.883810, 4.562038 } },
                {  6, { 0.728718, 3.302747, -1.025441, -0.427375, -3.908054, -1.740527, -2.372096, -0.673631, 1.811634, 3.240215, -2.596738, 4.487753 } },
                { 15, { 0.502282, 3.819604, -0.194239, -0.867121, -3.041055, -0.757073, -0.457640, -1.194520, 0.471496, 3.593195, -2.863604, 2.697087 } },
                { 12, { 0.165833, 4.941870, 0.591175, 0.499246, -3.625476, -3.037963, -0.982153, -1.800622, 0.804901, 0.968341, -5.800724, 4.715875 } },
                {  2, { 0.392269, 4.425013, -0.240027, 0.938992, -4.492474, -4.021417, -2.896609, -1.279733, 2.145039, 0.615361, -5.533858, 6.506542 } },
                { 13, { -0.170781, 5.134365, -0.167001, 0.998929, -5.777368, -2.464749, -3.330812, -2.080934, 1.083877, 2.170129, -4.246787, 6.432257 } },
                { 14, { -0.397217, 5.651222, 0.664201, 0.559183, -4.910370, -1.481294, -1.416356, -2.601823, -0.256262, 2.523108, -4.513652, 4.641591 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroInertialForcesResult = reader.getMacroInertialForces( );

    for ( auto it = macroInertialForcesAnswer.begin( ); it != macroInertialForcesAnswer.end( ); it++ ){

        auto r = macroInertialForcesResult->find( it->first );

        if ( r == macroInertialForcesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 67) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 68) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroBodyForcesAnswer
        =
            {
                {  5, { -0.585148, 2.009889, -0.810035, -2.353495, -1.481933, -1.394433, 2.604481, -0.764513, 0.620197, 2.312799, -2.676288, 0.791100 } },
                {  9, { 0.281864, 2.430634, 1.088528, -2.990224, -3.411629, -1.880188, 3.117476, -2.050000, 0.797479, 2.325336, -3.771395, -0.904908 } },
                {  8, { 2.143468, 1.904844, 1.400196, -1.881262, -3.641681, -0.353682, 3.474300, -0.518626, 0.209746, 0.757177, -5.497756, -2.613388 } },
                { 11, { 1.276457, 1.484099, -0.498368, -1.244534, -1.711985, 0.132072, 2.961305, 0.766861, 0.032463, 0.744640, -4.402649, -0.917380 } },
                {  3, { -2.214053, 1.620183, -2.790640, -3.412234, -0.154223, -2.107009, 1.316934, -2.590009, 1.051996, 1.727722, -3.562681, -0.900240 } },
                {  1, { -1.347041, 2.040929, -0.892076, -4.048962, -2.083919, -2.592763, 1.829929, -3.875496, 1.229279, 1.740259, -4.657788, -2.596248 } },
                {  6, { 0.514563, 1.515139, -0.580409, -2.940001, -2.313971, -1.066258, 2.186754, -2.344123, 0.641545, 0.172100, -6.384149, -4.304728 } },
                { 15, { -0.352449, 1.094394, -2.478973, -2.303272, -0.384275, -0.580504, 1.673759, -1.058635, 0.464263, 0.159563, -5.289042, -2.608720 } },
                { 12, { -3.842958, 1.230478, -4.771245, -4.470972, 1.173487, -2.819585, 0.029388, -4.415505, 1.483796, 1.142644, -4.449073, -2.591580 } },
                {  2, { -2.975947, 1.651223, -2.872681, -5.107700, -0.756209, -3.305339, 0.542383, -5.700993, 1.661078, 1.155182, -5.544180, -4.287588 } },
                { 13, { -1.114343, 1.125434, -2.561014, -3.998739, -0.986262, -1.778834, 0.899208, -4.169619, 1.073344, -0.412977, -7.270541, -5.996068 } },
                { 14, { -1.981354, 0.704688, -4.459577, -3.362010, 0.943435, -1.293080, 0.386213, -2.884132, 0.896062, -0.425514, -6.175435, -4.300060 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroBodyForcesResult = reader.getMacroBodyForces( );

    for ( auto it = macroBodyForcesAnswer.begin( ); it != macroBodyForcesAnswer.end( ); it++ ){

        auto r = macroBodyForcesResult->find( it->first );

        if ( r == macroBodyForcesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 69) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 70) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroSurfaceForcesAnswer
        =
            {
                {  5, { 0.700483, 1.063404, -0.561245, 0.196039, 0.442264, -1.455507, 0.778785, -0.781834, 0.010692, -2.689651, -1.646803, -2.731019 } },
                {  9, { 2.343569, 1.285853, -0.326911, 1.727253, 0.486287, -1.879686, 0.780393, -2.137662, 1.301457, -3.572794, -1.020227, -2.294672 } },
                {  8, { 3.754636, 3.023779, -0.181710, 1.511737, -1.141960, -0.566068, 0.325507, -0.975232, -0.594964, -4.529299, -0.496431, -3.091869 } },
                { 11, { 2.111550, 2.801331, -0.416044, -0.019477, -1.185983, -0.141889, 0.323899, 0.380597, -1.885729, -3.646155, -1.123007, -3.528216 } },
                {  3, { -0.684575, 2.891510, -1.646963, 1.544569, 1.248986, -0.666638, 0.967332, 1.124401, 0.193638, -0.784134, -2.692686, -3.146953 } },
                {  1, { 0.958511, 3.113959, -1.412630, 3.075784, 1.293009, -1.090817, 0.968940, -0.231428, 1.484403, -1.667278, -2.066109, -2.710605 } },
                {  6, { 2.369578, 4.851886, -1.267429, 2.860268, -0.335238, 0.222800, 0.514055, 0.931002, -0.412019, -2.623782, -1.542313, -3.507802 } },
                { 15, { 0.726492, 4.629437, -1.501763, 1.329054, -0.379261, 0.646979, 0.512447, 2.286831, -1.702784, -1.740639, -2.168890, -3.944150 } },
                { 12, { -2.069633, 4.719617, -2.732682, 2.893100, 2.055709, 0.122230, 1.155880, 3.030635, 0.376583, 1.121382, -3.738568, -3.562886 } },
                {  2, { -0.426547, 4.942065, -2.498348, 4.424315, 2.099731, -0.301949, 1.157488, 1.674806, 1.667349, 0.238239, -3.111991, -3.126539 } },
                { 13, { 0.984520, 6.679992, -2.353148, 4.208799, 0.471484, 1.011669, 0.702603, 2.837237, -0.229073, -0.718266, -2.588195, -3.923736 } },
                { 14, { -0.658566, 6.457543, -2.587482, 2.677584, 0.427461, 1.435848, 0.700994, 4.193065, -1.519838, 0.164878, -3.214772, -4.360083 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroSurfaceForcesResult = reader.getMacroSurfaceForces( );

    for ( auto it = macroSurfaceForcesAnswer.begin( ); it != macroSurfaceForcesAnswer.end( ); it++ ){

        auto r = macroSurfaceForcesResult->find( it->first );

        if ( r == macroSurfaceForcesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 71) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 72) & False\n";
            return 1;

        }

    }

    const std::unordered_map< uIntType, floatVector > macroExternalForcesAnswer
        =
            {
                {  5, { 0.115335, 3.073293, -1.371280, -2.157457, -1.039669, -2.849940, 3.383265, -1.546346, 0.630889, -0.376852, -4.323092, -1.939920 } },
                {  9, { 2.625433, 3.716487, 0.761617, -1.262971, -2.925342, -3.759873, 3.897868, -4.187662, 2.098937, -1.247459, -4.791622, -3.199580 } },
                {  8, { 5.898105, 4.928624, 1.218485, -0.369525, -4.783642, -0.919750, 3.799808, -1.493858, -0.385218, -3.772122, -5.994187, -5.705257 } },
                { 11, { 3.388007, 4.285430, -0.914412, -1.264011, -2.897968, -0.009817, 3.285205, 1.147458, -1.853266, -2.901515, -5.525657, -4.445597 } },
                {  3, { -2.898628, 4.511694, -4.437603, -1.867664, 1.094763, -2.773648, 2.284266, -1.465608, 1.245634, 0.943587, -6.255366, -4.047193 } },
                {  1, { -0.388530, 5.154888, -2.304706, -0.973178, -0.790910, -3.683581, 2.798869, -4.106924, 2.713682, 0.072981, -6.723897, -5.306853 } },
                {  6, { 2.884141, 6.367025, -1.847838, -0.079733, -2.649210, -0.843458, 2.700809, -1.413120, 0.229527, -2.451682, -7.926462, -7.812531 } },
                { 15, { 0.374043, 5.723830, -3.980735, -0.974219, -0.763536, 0.066476, 2.186206, 1.228196, -1.238521, -1.581076, -7.457932, -6.552870 } },
                { 12, { -5.912592, 5.950094, -7.503927, -1.577872, 3.229195, -2.697355, 1.185267, -1.384870, 1.860379, 2.264027, -8.187641, -6.154466 } },
                {  2, { -3.402494, 6.593289, -5.371029, -0.683386, 1.343522, -3.607288, 1.699871, -4.026186, 3.328427, 1.393420, -8.656171, -7.414127 } },
                { 13, { -0.129822, 7.805425, -4.914161, 0.210060, -0.514778, -0.767165, 1.601810, -1.332382, 0.844272, -1.131243, -9.858737, -9.919804 } },
                { 14, { -2.639920, 7.162231, -7.047059, -0.684426, 1.370896, 0.142768, 1.087207, 1.308934, -0.623776, -0.260636, -9.390206, -8.660143 } },
            };

    const std::unordered_map< uIntType, floatVector > *macroExternalForcesResult = reader.getMacroExternalForces( );

    for ( auto it = macroExternalForcesAnswer.begin( ); it != macroExternalForcesAnswer.end( ); it++ ){

        auto r = macroExternalForcesResult->find( it->first );

        if ( r == macroExternalForcesResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 73) & False\n";
            return 1;

        }
        else if ( !vectorTools::fuzzyEquals( r->second, it->second ) ){

            std::cout << r->first << ": "; vectorTools::print( r->second );
            std::cout << it->first << ": "; vectorTools::print( it->second );
            results << "test_initializeIncrement_Arlequin (test 74) & False\n";
            return 1;

        }

    }

    const uIntVector *freeMicroNodeIds = reader.getFreeMicroNodeIds( );
    const uIntVector *ghostMicroNodeIds = reader.getGhostMicroNodeIds( );

    for ( auto n  = ghostMicroNodeIds->begin( );
               n != ghostMicroNodeIds->end( );
               n++ ){

        if ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) != freeMicroNodeIds->end( ) ){
            std::cout << "*n: " << *n << "\n";
            results << "test_initializeIncrement_Arlequin (test 75) & False\n";
            return 1;
        }

    }

    uIntVector nodes;
    const stringVector *freeMicroDomainNames = reader.getFreeMicroDomainNames( );
    for ( auto domain  = freeMicroDomainNames->begin( );
               domain != freeMicroDomainNames->end( );
               domain++ ){

        reader._microscale->getSubDomainNodes( 0, *domain, nodes );

        for ( auto n = nodes.begin( ); n != nodes.end( ); n++ ){

            if ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) == freeMicroNodeIds->end( ) ){

                results << "test_initializeIncrement_Arlequin (test 76) & False\n";
                return 1;

            }

        }

    }

    const stringVector *ghostMicroDomainNames = reader.getGhostMicroDomainNames( );

    for ( auto domain  = ghostMicroDomainNames->begin( );
               domain != ghostMicroDomainNames->end( );
               domain++ ){

        reader._microscale->getSubDomainNodes( 0, *domain, nodes );

        for ( auto n = nodes.begin( ); n != nodes.end( ); n++ ){

            if ( ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) == freeMicroNodeIds->end( ) ) &&
                 ( std::find( ghostMicroNodeIds->begin( ), ghostMicroNodeIds->end( ), *n ) == ghostMicroNodeIds->end( ) ) ){

                results << "test_initializeIncrement_Arlequin (test 77) & False\n";
                return 1;

            }

        }

    }

    const uIntVector *freeMacroNodeIds = reader.getFreeMacroNodeIds( );
    const uIntVector *ghostMacroNodeIds = reader.getGhostMacroNodeIds( );

    const stringVector *ghostMacroDomainNames = reader.getGhostMacroDomainNames( );
    for ( auto domain  = ghostMacroDomainNames->begin( );
               domain != ghostMacroDomainNames->end( );
               domain++ ){

        reader._macroscale->getSubDomainNodes( 0, *domain, nodes );

        for ( auto n = nodes.begin( ); n != nodes.end( ); n++ ){

            if ( ( std::find( ghostMacroNodeIds->begin( ), ghostMacroNodeIds->end( ), *n ) == ghostMacroNodeIds->end( ) ) ){

                results << "test_initializeIncrement_Arlequin (test 78) & False\n";
                return 1;

            }

        }

    }

    const stringVector *freeMacroDomainNames = reader.getFreeMacroDomainNames( );
    for ( auto domain  = freeMacroDomainNames->begin( );
               domain != freeMacroDomainNames->end( );
               domain++ ){

        reader._macroscale->getSubDomainNodes( 0, *domain, nodes );

        for ( auto n = nodes.begin( ); n != nodes.end( ); n++ ){

            if ( ( std::find( ghostMacroNodeIds->begin( ), ghostMacroNodeIds->end( ), *n ) == ghostMacroNodeIds->end( ) ) &&
                 ( std::find( freeMacroNodeIds->begin( ), freeMacroNodeIds->end( ), *n ) == freeMacroNodeIds->end( ) ) ){

                results << "test_initializeIncrement_Arlequin (test 79) & False\n";
                return 1;

            }

        }

    }

    const DOFMap *microGlobalToLocalDOFMap = reader.getMicroGlobalToLocalDOFMap( );

    if ( microGlobalToLocalDOFMap->size( ) != ( freeMicroNodeIds->size( ) + ghostMicroNodeIds->size( ) ) ){

        results << "test_initializeIncrement_Arlequin (test 80) & False\n";
        return 1;

    }

    for ( auto n  = freeMicroNodeIds->begin( );
               n != freeMicroNodeIds->end( );
               n++ ){

        if ( microGlobalToLocalDOFMap->find( *n ) == microGlobalToLocalDOFMap->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 81) & False\n";
            return 1;

        }

    }

    for ( auto n  = ghostMicroNodeIds->begin( );
               n != ghostMicroNodeIds->end( );
               n++ ){

        if ( microGlobalToLocalDOFMap->find( *n ) == microGlobalToLocalDOFMap->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 82) & False\n";
            return 1;

        }

    }

    const DOFMap *macroGlobalToLocalDOFMap = reader.getMacroGlobalToLocalDOFMap( );

    if ( macroGlobalToLocalDOFMap->size( ) != ( freeMacroNodeIds->size( ) + ghostMacroNodeIds->size( ) ) ){

        results << "test_initializeIncrement_Arlequin (test 83) & False\n";
        return 1;

    }

    for ( auto n  = freeMacroNodeIds->begin( );
               n != freeMacroNodeIds->end( );
               n++ ){

        if ( macroGlobalToLocalDOFMap->find( *n ) == macroGlobalToLocalDOFMap->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 84) & False\n";
            return 1;

        }

    }

    for ( auto n  = ghostMacroNodeIds->begin( );
               n != ghostMacroNodeIds->end( );
               n++ ){

        if ( macroGlobalToLocalDOFMap->find( *n ) == macroGlobalToLocalDOFMap->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 85) & False\n";
            return 1;

        }

    }

    if ( !reader.microBodyForceDefined( ) ){

        results << "test_initializeIncrement_Arlequin (test 86) & False\n";
        return 1;

    }

    if ( !reader.microSurfaceForceDefined( ) ){

        results << "test_initializeIncrement_Arlequin (test 87) & False\n";
        return 1;

    }

    if ( !reader.microAccelerationDefined( ) ){

        results << "test_initializeIncrement_Arlequin (test 88) & False\n";
        return 1;

    }

    if ( reader.useReconstructedMassCenters( ) ){

        results << "test_initializeIncrement_Arlequin (test 89) & False\n";
        return 1;

    }

    if ( !reader.microVelocitiesDefined( ) ){

        results << "test_initializeIncrement_Arlequin (test 90) & False\n";
        return 1;

    }

    if ( !reader.macroAccelerationDefined( ) ){

        results << "test_initializeIncrement_Arlequin (test 91) & False\n";
        return 1;

    }

    if ( !reader.macroVelocitiesDefined( ) ){

        results << "test_initializeIncrement_Arlequin (test 92) & False\n";
        return 1;

    }


    if ( !reader.microInternalForceDefined( ) ){

        results << "test_initializeIncrement_Arlequin (test 93) & False\n";
        return 1;

    }

    if ( !reader.macroInternalForceDefined( ) ){

        results << "test_initializeIncrement_Arlequin (test 94) & False\n";
        return 1;

    }

    if ( !reader.macroInertialForceDefined( ) ){

        results << "test_initializeIncrement_Arlequin (test 95) & False\n";
        return 1;

    }

    if ( !reader.macroExternalForceDefined( ) ){

        results << "test_initializeIncrement_Arlequin (test 96) & False\n";
        return 1;

    }

    const std::string macroReferenceDensityTypesAnswer = "constant";
    const floatVector macroReferenceDensitiesAnswer = { 2. };
    const std::unordered_map< unsigned int, floatVector > *macroReferenceDensitiesResult = reader.getMacroReferenceDensities( );
    const std::unordered_map< unsigned int, std::string > *macroReferenceDensityTypesResult = reader.getMacroReferenceDensityTypes( );

    for ( auto mRDR = macroReferenceDensitiesResult->begin( ); mRDR != macroReferenceDensitiesResult->end( ); mRDR++ ){

        if ( !vectorTools::fuzzyEquals( macroReferenceDensitiesAnswer, mRDR->second ) ){

            results << "test_initializeIncrement_Arlequin (test 97) & False\n";
            return 1;

        }

    }

    for ( auto mRDTR = macroReferenceDensityTypesResult->begin( ); mRDTR != macroReferenceDensityTypesResult->end( ); mRDTR++ ){

        if ( macroReferenceDensityTypesAnswer.compare( mRDTR->second ) != 0 ){

            results << "test_initializeIncrement_Arlequin (test 98) & False\n";
            return 1;

        }

    }

    const std::string macroReferenceMomentOfInertiaTypesAnswer = "constant";
    const floatVector macroReferenceMomentsOfInertiaAnswer = { 1e-5, 2e-5, 3e-5,
                                                               2e-5, 4e-5, 5e-5,
                                                               3e-5, 5e-5, 6e-5 };
    const std::unordered_map< unsigned int, floatVector > *macroReferenceMomentsOfInertiaResult
        = reader.getMacroReferenceMomentsOfInertia( );
    const std::unordered_map< unsigned int, std::string > *macroReferenceMomentOfInertiaTypesResult
        = reader.getMacroReferenceMomentOfInertiaTypes( );

    for ( auto mRMIR  = macroReferenceMomentsOfInertiaResult->begin( );
               mRMIR != macroReferenceMomentsOfInertiaResult->end( );
               mRMIR++ ){

        if ( !vectorTools::fuzzyEquals( macroReferenceMomentsOfInertiaAnswer, mRMIR->second ) ){

            results << "test_initializeIncrement_Arlequin (test 99) & False\n";
            return 1;

        }

    }

    for ( auto mRMITR  = macroReferenceMomentOfInertiaTypesResult->begin( );
               mRMITR != macroReferenceMomentOfInertiaTypesResult->end( ); mRMITR++ ){

        if ( macroReferenceMomentOfInertiaTypesAnswer.compare( mRMITR->second ) != 0 ){

            results << "test_initializeIncrement_Arlequin (test 100) & False\n";
            return 1;

        }

    }

    if ( !reader.microSurfaceForceDefined( ) ){

        results << "test_initializeIncrement_Arlequin (test 101) & False\n";
        return 1;

    }

    if ( !reader.microExternalForceDefined( ) ){

        results << "test_initializeIncrement_Arlequin (test 102) & False\n";
        return 1;

    }

    if ( !reader.extractPreviousDOFValues( ) ){

        results << "test_initializeIncrement_Arlequin (test 103) & False\n";
        return 1;

    }

    const floatType DtAnswer = 1.;
    const floatType* DtResult = reader.getDt( );

    if ( !vectorTools::fuzzyEquals( DtAnswer, *DtResult ) ){

        results << "test_initializeIncrement_Arlequin (test 104) & False\n";
        return 1;

    }

    const floatType newmarkGammaAnswer = 0.50;
    const floatType newmarkBetaAnswer  = 0.25;

    if ( !vectorTools::fuzzyEquals( newmarkGammaAnswer, *reader.getNewmarkGamma( ) ) ){

        results << "test_initializeIncrement_Arlequin (test 105) & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( newmarkBetaAnswer, *reader.getNewmarkBeta( ) ) ){

        results << "test_initializeIncrement_Arlequin (test 106) & False\n";
        return 1;

    }

    const std::unordered_map< uIntType, stringVector > macroCellToDomainMapAnswer
        =
        {
            { 1, { "ghost_nodeset_volume_1", "ghost_nodeset_volume_2", "ghost_nodeset_volume_3", "ghost_nodeset_volume_4",
                   "ghost_nodeset_volume_5", "ghost_nodeset_volume_6", "ghost_nodeset_volume_7", "ghost_nodeset_volume_8" } },
            { 2, { "free_nodeset_volume_1", "free_nodeset_volume_2", "free_nodeset_volume_3", "free_nodeset_volume_4",
                   "free_nodeset_volume_5", "free_nodeset_volume_6", "free_nodeset_volume_7", "free_nodeset_volume_8" } },
        };

    const std::unordered_map< uIntType, stringVector > *macroCellToDomainMapResult = reader.getMacroCellToDomainMap( );

    for ( auto a = macroCellToDomainMapAnswer.begin( ); a != macroCellToDomainMapAnswer.end( ); a++ ){

        auto r = macroCellToDomainMapResult->find( a->first );

        if ( r == macroCellToDomainMapResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 107) & False\n";
            return 1;

        }

        if ( a->second.size( ) != r->second.size( ) ){

            results << "test_initializeIncrement_Arlequin (test 108) & False\n";
            return 1;

        }

        for ( unsigned int i = 0; i < a->second.size( ); i++ ){

            if ( a->second[ i ].compare( r->second[ i ] ) != 0 ){

                results << "test_initializeIncrement_Arlequin (test 109) & False\n";
                return 1;

            }

        }        

    }

    const std::unordered_map< std::string, uIntType > microDomainIDMapAnswer
        =
        {
            { "free_nodeset_volume_1",   0 },
            { "free_nodeset_volume_2",   1 },
            { "free_nodeset_volume_3",   2 },
            { "free_nodeset_volume_4",   3 },
            { "free_nodeset_volume_5",   4 },
            { "free_nodeset_volume_6",   5 },
            { "free_nodeset_volume_7",   6 },
            { "free_nodeset_volume_8",   7 },
            { "ghost_nodeset_volume_1",  8 },
            { "ghost_nodeset_volume_2",  9 },
            { "ghost_nodeset_volume_3", 10 },
            { "ghost_nodeset_volume_4", 11 },
            { "ghost_nodeset_volume_5", 12 },
            { "ghost_nodeset_volume_6", 13 },
            { "ghost_nodeset_volume_7", 14 },
            { "ghost_nodeset_volume_8", 15 },
        };

    const std::unordered_map< std::string, uIntType > *microDomainIDMapResult = reader.getMicroDomainIDMap( );

    for ( auto a = microDomainIDMapAnswer.begin( ); a != microDomainIDMapAnswer.end( ); a++ ){

        auto r = microDomainIDMapResult->find( a->first );

        if ( r == microDomainIDMapResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 110) & False\n";
            return 1;

        }

        if ( r->second != a->second ){

            results << "test_initializeIncrement_Arlequin (test 111) & False\n";
            return 1;

        }

    }

    //Test the Arlequin weights
    const std::unordered_map< uIntType, floatType > arlequinWeightsAnswer
        =
        {
            {  5, 0.000 },
            {  9, 0.000 },
            {  8, 0.000 },
            { 11, 0.000 },
            {  3, 0.500 },
            {  1, 0.500 },
            {  6, 0.500 },
            { 15, 0.500 },
            { 12, 1.000 },
            {  2, 1.000 },
            { 13, 1.000 },
            { 14, 1.000 },
        };

    if ( !reader.useArlequinCoupling( ) ){

        results << "test_initializeIncrement_Arlequin (test 112) & False\n";
        return 1;

    }

    const std::unordered_map< uIntType, floatType > *arlequinWeightsResult = reader.getMacroArlequinWeights( );

    for ( auto a = arlequinWeightsAnswer.begin( ); a != arlequinWeightsAnswer.end( ); a++ ){

        auto r = arlequinWeightsResult->find( a->first );

        if ( r == arlequinWeightsResult->end( ) ){

            results << "test_initializeIncrement_Arlequin (test 113) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( a->second, r->second ) ){

            results << "test_initializeIncrement_Arlequin (test 114) & False\n";
            return 1;

        }

    }

    results << "test_initializeIncrement_Arlequin & True\n";
    return 0;
}


int test_getFreeMicroDomainNames( std::ostream &results ){
    /*!
     * Test getting a pointer to the free micro domain names
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_getFreeMicroDomainNames & False\n";
        return 1;
    }

    stringVector answer = { "free_nodeset_volume_1",
                            "free_nodeset_volume_2",
                            "free_nodeset_volume_3",
                            "free_nodeset_volume_4",
                            "free_nodeset_volume_5",
                            "free_nodeset_volume_6",
                            "free_nodeset_volume_7",
                            "free_nodeset_volume_8" };

    const stringVector *result = reader.getFreeMicroDomainNames( );

    unsigned int indx = 0;
    for ( auto it = result->begin( ); it != result->end( ); it++ ){

        if ( it->compare( answer[ indx ] ) != 0 ){

            results << "test_getFreeMicroDomainNames (test 1) & False\n";
            return 1;

        }

        indx++;

    }

    results << "test_getFreeMicroDomainNames & True\n";
    return 0;
}

int test_getGhostMicroDomainNames( std::ostream &results ){
    /*!
     * Test getting a pointer to the ghost micro domain names
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_getGhostMicroDomainNames & False\n";
        return 1;
    }

    stringVector answer = { "ghost_nodeset_volume_1",
                            "ghost_nodeset_volume_2",
                            "ghost_nodeset_volume_3",
                            "ghost_nodeset_volume_4",
                            "ghost_nodeset_volume_5",
                            "ghost_nodeset_volume_6",
                            "ghost_nodeset_volume_7",
                            "ghost_nodeset_volume_8" };

    const stringVector *result = reader.getGhostMicroDomainNames( );

    unsigned int indx = 0;
    for ( auto it = result->begin( ); it != result->end( ); it++ ){

        if ( it->compare( answer[ indx ] ) != 0 ){

            results << "test_getGhostMicroDomainNames (test 1) & False\n";
            return 1;

        }

        indx++;

    }

    results << "test_getGhostMicroDomainNames & True\n";
    return 0;
}

int test_getCouplingInitialization( std::ofstream &results ){
    /*!
     * Test getting the coupling initialization from the configuration file
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_getCouplingInitialization & False\n";
        return 1;
    }

    YAML::Node couplingInitialization = reader.getCouplingInitialization( );

    if ( !couplingInitialization ){
        results << "test_getCouplingInitialization (test 1) & False\n";
        return 1;
    }

    std::string typeAnswer = "use_first_increment";
    if ( couplingInitialization[ "type" ].as<std::string>( ).compare( typeAnswer ) ){
        results << "test_getCouplingInitialization (test 2) & False\n";
        return 1;
    }

    std::string projectionTypeAnswer = "averaged_l2_projection";
    bool useReconstructedMassCentersAnswer = false;
    floatType potentialWeightingFactorAnswer = 0.5;
    floatType kineticWeightingFactorAnswer = 0.5;
    std::string potentialPartitioningTypeAnswer = "volume_fraction";
    std::string kineticPartitioningTypeAnswer = "volume_fraction";

    if ( couplingInitialization[ "projection_type" ] ){

        if( couplingInitialization[ "projection_type" ].as< std::string >( ).compare( projectionTypeAnswer ) != 0 ){

            results << "test_getCouplingInitialization (test 3) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 4) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "use_reconstructed_mass_centers" ] ){

        if ( couplingInitialization[ "use_reconstructed_mass_centers" ].as< bool >( ) != useReconstructedMassCentersAnswer ){

            results << "test_getCouplingInitialization (test 5) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 6) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "potential_energy_weighting_factor" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "potential_energy_weighting_factor" ].as< floatType >( ),
                                        potentialWeightingFactorAnswer ) ){
    
            results << "test_getCouplingInitialization (test 7) & False\n";
            return 1;
    
        }

    }
    else{

        results << "test_getCouplingInitialization (test 8) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "kinetic_energy_weighting_factor" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "kinetic_energy_weighting_factor" ].as< floatType >( ),
                                        kineticWeightingFactorAnswer ) ){
    
            results << "test_getCouplingInitialization (test 9) & False\n";
            return 1;
    
        }

    }
    else{

        results << "test_getCouplingInitialization (test 10) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "potential_energy_partitioning_coefficient" ][ "type" ] ){

        if ( couplingInitialization[ "potential_energy_partitioning_coefficient" ][ "type" ].as< std::string >( ).compare( potentialPartitioningTypeAnswer ) != 0 ){
    
            results << "test_getCouplingInitialization (test 11) & False\n";
            return 1;
    
        }

    }
    else{

        results << "test_getCouplingInitialization (test 12) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "kinetic_energy_partitioning_coefficient" ][ "type" ] ){

        if ( couplingInitialization[ "kinetic_energy_partitioning_coefficient" ][ "type" ].as< std::string >( ).compare( kineticPartitioningTypeAnswer ) != 0 ){
    
            results << "test_getCouplingInitialization (test 13) & False\n";
            return 1;
    
        }

    }
    else{

        results << "test_getCouplingInitialization (test 14) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "macro_proportionality_coefficient" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "macro_proportionality_coefficient" ].as< floatType >( ),
                                        1e-3 ) ){

            results << "test_getCouplingInitialization (test 15) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 16) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "micro_proportionality_coefficient" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "micro_proportionality_coefficient" ].as< floatType >( ),
                                        1e-3 ) ){

            results << "test_getCouplingInitialization (test 17) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 18) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "macro_internal_force_sign" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "macro_internal_force_sign" ].as< floatType >( ), -1. ) ){

            results << "test_getCouplingInitialization (test 19) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 20) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "macro_external_force_sign" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "macro_external_force_sign" ].as< floatType >( ), 1. ) ){

            results << "test_getCouplingInitialization (test 21) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 22) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "micro_internal_force_sign" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "micro_internal_force_sign" ].as< floatType >( ), 1. ) ){

            results << "test_getCouplingInitialization (test 23) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 24) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "micro_external_force_sign" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "micro_external_force_sign" ].as< floatType >( ), 1. ) ){

            results << "test_getCouplingInitialization (test 25) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 26) & False\n";
        return 1;

    }

    if ( !couplingInitialization[ "extract_previous_dof_values" ].as< bool >( ) ){

        results << "test_getCouplingInitialization (test 27) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "previous_micro_increment" ].as< uIntType >( ) != 0 ){

        results << "test_getCouplingInitialization (test 28) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "previous_macro_increment" ].as< uIntType >( ) != 0 ){

        results << "test_getCouplingInitialization (test 29) & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( couplingInitialization[ "update_displacement" ][ "Newmark-beta_parameters" ][ "beta" ].as< floatType >( ), 0.25 ) ){

        results << "test_getCouplingInitialization (test 30) & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( couplingInitialization[ "update_displacement" ][ "Newmark-beta_parameters" ][ "gamma" ].as< floatType >( ), 0.5 ) ){

        results << "test_getCouplingInitialization (test 31) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "output_reference_information" ][ "filename" ].as< std::string >( ).compare( "reference_information" ) != 0 ){

        results << "test_getCouplingInitialization (test 32) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "output_homogenized_response" ][ "filename" ].as< std::string >( ).compare( "homogenized_response" ) != 0 ){

        results << "test_getCouplingInitialization (test 33) & False\n";
        return 1;

    }

    if ( ( couplingInitialization[ "output_updated_dof" ] ) && ( !couplingInitialization[ "output_updated_dof" ].IsScalar( ) ) ){

        if ( couplingInitialization[ "output_updated_dof" ][ "macroscale_filename" ].as< std::string >( ).compare( "macroscale_dof" ) != 0 ){
    
            results << "test_getCouplingInitialization (test 34) & False\n";
            return 1;
    
        }
    }
    else{

        results << "test_getCouplingInitialization (test 34) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "output_updated_dof" ][ "microscale_filename" ].as< std::string >( ).compare( "microscale_dof" ) != 0 ){

        results << "test_getCouplingInitialization (test 35) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "reference_filename" ].as< std::string >( ).compare( "reference_information.xdmf" ) != 0 ){

        std::cout << couplingInitialization[ "reference_filename" ].as< std::string >( ) << "\n";
        results << "test_getCouplingInitialization (test 36) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "output_homogenized_response" ][ "filetype" ].as< std::string >( ).compare( "XDMF" ) != 0 ){

        results << "test_getCouplingInitialization (test 37) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "output_homogenized_response" ][ "mode" ].as< std::string >( ).compare( "write" ) != 0 ){

        results << "test_getCouplingInitialization (test 38) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "output_updated_dof" ][ "macroscale_filetype" ].as< std::string >( ).compare( "XDMF" ) != 0 ){

        results << "test_getCouplingInitialization (test 39) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "output_updated_dof" ][ "microscale_filetype" ].as< std::string >( ).compare( "XDMF" ) != 0 ){

        results << "test_getCouplingInitialization (test 40) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "macro_inertial_force_sign" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "macro_inertial_force_sign" ].as< floatType >( ), 1. ) ){

            results << "test_getCouplingInitialization (test 41) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 42) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "macro_body_force_sign" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "macro_body_force_sign" ].as< floatType >( ), 1. ) ){

            results << "test_getCouplingInitialization (test 43) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 44) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "macro_surface_force_sign" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "macro_surface_force_sign" ].as< floatType >( ), 1. ) ){

            results << "test_getCouplingInitialization (test 45) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 46) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "micro_inertial_force_sign" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "micro_inertial_force_sign" ].as< floatType >( ), 1. ) ){

            results << "test_getCouplingInitialization (test 47) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 48) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "micro_body_force_sign" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "micro_body_force_sign" ].as< floatType >( ), 1. ) ){

            results << "test_getCouplingInitialization (test 49) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 50) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "micro_surface_force_sign" ] ){

        if ( !vectorTools::fuzzyEquals( couplingInitialization[ "micro_surface_force_sign" ].as< floatType >( ), 1. ) ){

            results << "test_getCouplingInitialization (test 51) & False\n";
            return 1;

        }

    }
    else{

        results << "test_getCouplingInitialization (test 52) & False\n";
        return 1;

    }

    if ( couplingInitialization[ "solve_coupling_odes_at_microdomains" ].as< bool >( ) ){

        results << "test_getCouplingInitialization (test 53) & False\n";
        return 1;

    }

    if ( reader.solveCouplingODEsAtMicroDomains( ) ){

        results << "test_getCouplingInitialization (test 54) & False\n";
        return 1;

    }

    results << "test_getCouplingInitialization & True\n";
    return 0;

}

int test_getVolumeReconstructionConfig( std::ofstream &results ){
    /*!
     * Test getting the volume reconstruction configuration from the configuration file
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_getVolumeReconstructionConfig & False\n";
        return 1;
    }

    const YAML::Node vRInitialization = reader.getVolumeReconstructionConfig( );

    if ( !vRInitialization ){
        results << "test_getVolumeReconstructionConfig (test 1) & False\n";
        return 1;
    }

    std::string typeAnswer = "dual_contouring";
    if ( vRInitialization[ "type" ].as<std::string>( ).compare( typeAnswer ) ){
        results << "test_getVolumeReconstructionConfig (test 2) & False\n";
        return 1;
    }

    floatType toleranceAnswer = 1e-2;
    if ( !vectorTools::fuzzyEquals( vRInitialization[ "element_contain_tolerance" ].as< floatType >( ), toleranceAnswer ) ){
        results << "test_getVolumeReconstructionConfig (test 3) & False\n";
        return 1;
    }

    floatType useMacroNormalsAnswer = true;
    if ( !( vRInitialization[ "use_macro_normals" ].as< bool >( ) == useMacroNormalsAnswer ) ){
        results << "test_getVolumeReconstruction (test 4) & False\n";
        return 1;
    }

    results << "test_getVolumeReconstructionConfig & True\n";
    return 0;

}

int test_getFreeMacroDomainNames( std::ofstream &results ){
    /*!
     * Test getting the free macro volume sets from the configuration file
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_getFreeMacroDomainNames & False\n";
        return 1;
    }

    stringVector answer = { "free_nodes" };

    const stringVector *result = reader.getFreeMacroDomainNames( );

    unsigned int indx = 0;
    for ( auto it = result->begin( ); it != result->end( ); it++ ){

        if ( it->compare( answer[ indx ] ) != 0 ){

            results << "test_getFreeMacroDomainNames (test 1) & False\n";
            return 1;

        }

        indx++;

    }

    results << "test_getFreeMacroDomainNames & True\n";
    return 0;

}

int test_getGhostMacroDomainNames( std::ofstream &results ){
    /*!
     * Test getting the ghost macro volume sets from the configuration file
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_getGhostMacroDomainNames & False\n";
        return 1;
    }

    stringVector answer = { "ghost_nodes" };

    const stringVector *result = reader.getGhostMacroDomainNames( );

    unsigned int indx = 0;
    for ( auto it = result->begin( ); it != result->end( ); it++ ){

        if ( it->compare( answer[ indx ] ) != 0 ){

            results << "test_getGhostMacroDomainNames (test 1) & False\n";
            return 1;

        }

        indx++;

    }

    results << "test_getGhostMacroDomainNames & True\n";
    return 0;

}

int test_getMicroDomainSurfaceApproximateSplitCount( std::ostream &results ){
    /*!
     * Test getting a pointer to the approximate number of surfaces to split the micro
     * domains into
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_getMicroDomainSurfaceApproximateSplitCount & False\n";
        return 1;
    }

    std::unordered_map< std::string, uIntType > answer
        =
        {
            { "free_nodeset_volume_1", 6 },
            { "free_nodeset_volume_2", 6 },
            { "free_nodeset_volume_3", 6 },
            { "free_nodeset_volume_4", 6 },
            { "free_nodeset_volume_5", 6 },
            { "free_nodeset_volume_6", 6 },
            { "free_nodeset_volume_7", 6 },
            { "free_nodeset_volume_8", 6 },
            { "ghost_nodeset_volume_1", 6 },
            { "ghost_nodeset_volume_2", 6 },
            { "ghost_nodeset_volume_3", 6 },
            { "ghost_nodeset_volume_4", 6 },
            { "ghost_nodeset_volume_5", 6 },
            { "ghost_nodeset_volume_6", 6 },
            { "ghost_nodeset_volume_7", 6 },
            { "ghost_nodeset_volume_8", 6 },
        };

    const std::unordered_map< std::string, uIntType > *result = reader.getMicroDomainSurfaceApproximateSplitCount( );

    for ( auto a = answer.begin( ); a != answer.end( ); a++ ){

        auto r = result->find( a->first );

        if ( r == result->end( ) ){

            results << "test_getMicroDomainSurfaceApproximateSplitCount (test 1) & False\n";
            return 1;

        }

        if ( r->second != a->second ){

            results << "test_getMicroDomainSurfaceApproximateSplitCount (test 2) & False\n";
            return 1;

        }

    }

    results << "test_getMicroDomainSurfaceApproximateSplitCount & True\n";
    return 0;

}

//int test_getFreeMicroSurfaceApproximateSplitCount( std::ostream &results ){
//    /*!
//     * Test getting a pointer to the approximate number of surfaces to split a micro
//     * domain into.
//     *
//     * :param std::ofstream &results: The output file
//     */
//
//    std::string filename = "testConfig.yaml";
//    inputFileProcessor::inputFileProcessor reader( filename );
//
//    if ( reader.getError( ) ){
//        reader.getError( )->print( );
//        results << "test_getFreeMicroSurfaceApproximateSplitCount & False\n";
//        return 1;
//    }
//
//    uIntVector answer( 8, 6 );
//
//    const uIntVector *result = reader.getFreeMicroSurfaceApproximateSplitCount( );
//
//    unsigned int indx = 0;
//    for ( auto it = result->begin( ); it != result->end( ); it++ ){
//
//        if ( !vectorTools::fuzzyEquals( *it, answer[ indx ] ) ){
//
//            results << "test_getFreeMicroSurfaceApproximateSplitCount (test 1) & False\n";
//            return 1;
//
//        }
//
//        indx++;
//
//    }
//
//    results << "test_getFreeMicroSurfaceApproximateSplitCount & True\n";
//    return 0;
//}
//
//int test_getGhostMicroSurfaceApproximateSplitCount( std::ostream &results ){
//    /*!
//     * Test getting a pointer to the approximate number of surfaces to split a micro
//     * domain into.
//     *
//     * :param std::ofstream &results: The output file
//     */
//
//    std::string filename = "testConfig.yaml";
//    inputFileProcessor::inputFileProcessor reader( filename );
//
//    if ( reader.getError( ) ){
//        reader.getError( )->print( );
//        results << "test_getGhostMicroSurfaceApproximateSplitCount & False\n";
//        return 1;
//    }
//
//    uIntVector answer( 8, 6 );
//
//    const uIntVector *result = reader.getGhostMicroSurfaceApproximateSplitCount( );
//
//    unsigned int indx = 0;
//    for ( auto it = result->begin( ); it != result->end( ); it++ ){
//
//        if ( !vectorTools::fuzzyEquals( *it, answer[ indx ] ) ){
//
//            results << "test_getGhostMicroSurfaceApproximateSplitCount (test 1) & False\n";
//            return 1;
//
//        }
//
//        indx++;
//
//    }
//
//    results << "test_getGhostMicroSurfaceApproximateSplitCount & True\n";
//    return 0;
//}

int test_outputReferenceInformation( std::ofstream &results ){
    /*!
     * Test whether the reference information should be output
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_outputReferenceInformation & False\n";
        return 1;
    }

    if ( !reader.outputReferenceInformation( ) ){

        results << "test_outputReferenceInformation (test 1) & False\n";
        return 1;

    }

    results << "test_outputReferenceInformation & True\n";
    return 1;

}

int test_outputHomogenizedInformation( std::ofstream &results ){
    /*!
     * Test whether the homogenized information should be output
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_outputHomogenizedInformation & False\n";
        return 1;
    }

    if ( !reader.outputHomogenizedInformation( ) ){

        results << "test_outputHomogenizedInformation (test 1) & False\n";
        return 1;

    }

    results << "test_outputHomogenizedInformation & True\n";
    return 1;

}

int test_outputUpdatedDOF( std::ofstream &results ){
    /*!
     * Test whether the updated degree of freedom information should be output
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_outputUpdatedDOF & False\n";
        return 1;
    }

    if ( !reader.outputUpdatedDOF( ) ){

        results << "test_outputUpdatedDOF (test 1) & False\n";
        return 1;

    }

    results << "test_outputUpdatedDOF & True\n";
    return 1;

}

int main(){
    /*!
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or False 
    if the test passes or fails respectively.
    */

    //Run the setup
    std::unique_ptr< errorNode > error;
    error.reset( _createXDMFDatafiles( ) );
    if ( error ){

        error->print( );
        return 1;

    }

    //Open the results file
    std::ofstream results;
    results.open("results.tex");

    test_openConfigurationFile( results );
    test_setConfigurationFile( results );
    test_initializeFileInterfaces( results );
    test_initializeIncrement( results );
    test_initializeIncrement_Arlequin( results );
    test_getFreeMicroDomainNames( results );
    test_getGhostMicroDomainNames( results );
    test_getFreeMacroDomainNames( results );
    test_getGhostMacroDomainNames( results );
    test_getCouplingInitialization( results );
    test_getVolumeReconstructionConfig( results );
    test_getMicroDomainSurfaceApproximateSplitCount( results );
    test_outputReferenceInformation( results );
    test_outputHomogenizedInformation( results );
    test_outputUpdatedDOF( results );

    //Close the results file
    results.close();

}
