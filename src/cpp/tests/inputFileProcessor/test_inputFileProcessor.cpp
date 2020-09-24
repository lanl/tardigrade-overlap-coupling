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

    result = reader.getFreeMacroCellMicroDomainCounts( );

    if ( result->size( ) == 0 ){

        results << "test_openConfigurationFile (test 8) & False\n";
        return 1;

    }

    for ( auto v = result->begin( ); v != result->end( ); v++ ){

        if ( *v != 8 ){

            results << "test_openConfigurationFile (test 9) & False\n";
            return 1;

        }

    }

    result = reader.getGhostMacroCellMicroDomainCounts( );

    if ( result->size( ) == 0 ){

        results << "test_openConfigurationFile (test 10) & False\n";
        return 1;

    }

    for ( auto v = result->begin( ); v != result->end( ); v++ ){

        if ( *v != 8 ){

            results << "test_openConfigurationFile (test 11) & False\n";
            return 1;

        }

    }

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

    const uIntVector freeMacroCellMicroDomainCountsAnswer = { 8 };
    const uIntVector ghostMacroCellMicroDomainCountsAnswer = { 8 };

    const uIntVector *freeMacroCellMicroDomainCountsResult = reader.getFreeMacroCellMicroDomainCounts( );
    const uIntVector *ghostMacroCellMicroDomainCountsResult = reader.getGhostMacroCellMicroDomainCounts( );

    if ( !vectorTools::fuzzyEquals( freeMacroCellMicroDomainCountsAnswer, *freeMacroCellMicroDomainCountsResult ) ){
        results << "test_initializeIncrement (test 26) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( ghostMacroCellMicroDomainCountsAnswer, *ghostMacroCellMicroDomainCountsResult ) ){
        results << "test_initializeIncrement (test 27) & False\n";
        return 1;
    }

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

    const std::unordered_map< uIntType, floatVector > microSurfaceTractionsAnswer
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

    const std::unordered_map< uIntType, floatVector > *microSurfaceTractionsResult = reader.getMicroSurfaceTractions( );

    for ( auto it = microSurfaceTractionsAnswer.begin( ); it != microSurfaceTractionsAnswer.end( ); it++ ){

        auto r = microSurfaceTractionsResult->find( it->first );

        if ( r == microSurfaceTractionsResult->end( ) ){

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
















    const floatVector macroDisplacementsAnswer =
        {
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
            0., 0., -0.001,
            0., 0., -0.001,
        };

    const floatVector *macroDisplacementsResult = reader.getMacroDisplacements( );

    if ( !vectorTools::fuzzyEquals( macroDisplacementsAnswer, *macroDisplacementsResult ) ){
        vectorTools::print( *macroDisplacementsResult );
        results << "test_initializeIncrement (test 12) & False\n";
        return 1;
    }

//    const uIntVector *nonOverlappedMicroNodeIds = reader.getNonOverlappedMicroNodeIds( );
    const uIntVector *freeMicroNodeIds = reader.getFreeMicroNodeIds( );
    const uIntVector *ghostMicroNodeIds = reader.getGhostMicroNodeIds( );

//    for ( auto n  = freeMicroNodeIds->begin( );
//               n != freeMicroNodeIds->end( );
//               n++ ){
//
//        if ( std::find( nonOverlappedMicroNodeIds->begin( ), nonOverlappedMicroNodeIds->end( ), *n ) != nonOverlappedMicroNodeIds->end( ) ){
//            std::cout << "*n: " << *n << "\n";
//            results << "test_initializeIncrement (test 13) & False\n";
//            return 1;
//        }
//    }

    for ( auto n  = ghostMicroNodeIds->begin( );
               n != ghostMicroNodeIds->end( );
               n++ ){

        if ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) != freeMicroNodeIds->end( ) ){
            std::cout << "*n: " << *n << "\n";
            results << "test_initializeIncrement (test 14) & False\n";
            return 1;
        }

//        if ( std::find( nonOverlappedMicroNodeIds->begin( ), nonOverlappedMicroNodeIds->end( ), *n ) != nonOverlappedMicroNodeIds->end( ) ){
//            std::cout << "*n: " << *n << "\n";
//            results << "test_initializeIncrement (test 15) & False\n";
//            return 1;
//        }

    }

//    const stringVector *nonOverlappedMicroSurfaceNames = reader.getNonOverlappedMicroSurfaceNames( );
//    uIntVector nodes;
//
//    for ( auto surface  = nonOverlappedMicroSurfaceNames->begin( );
//               surface != nonOverlappedMicroSurfaceNames->end( );
//               surface++ ){
//
//        reader._microscale->getSubDomainNodes( 0, *surface, nodes );
//
//        for ( auto n = nodes.begin( ); n != nodes.end( ); n++ ){
//
//            if ( std::find( nonOverlappedMicroNodeIds->begin( ), nonOverlappedMicroNodeIds->end( ), *n ) == nonOverlappedMicroNodeIds->end( ) ){
//                results << "test_initializeIncrement (test 16) & False\n";
//                return 1;
//            }
//
//        }
//
//    }

    uIntVector nodes;
    const stringVector *freeMicroDomainNames = reader.getFreeMicroDomainNames( );
    for ( auto domain  = freeMicroDomainNames->begin( );
               domain != freeMicroDomainNames->end( );
               domain++ ){

        reader._microscale->getSubDomainNodes( 0, *domain, nodes );

        for ( auto n = nodes.begin( ); n != nodes.end( ); n++ ){

//            if ( ( std::find( nonOverlappedMicroNodeIds->begin( ), nonOverlappedMicroNodeIds->end( ), *n ) == nonOverlappedMicroNodeIds->end( ) ) &&
//                 ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) == freeMicroNodeIds->end( ) ) ){
            if ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) == freeMicroNodeIds->end( ) ){

                results << "test_initializeIncrement (test 17) & False\n";
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

//            if ( ( std::find( nonOverlappedMicroNodeIds->begin( ), nonOverlappedMicroNodeIds->end( ), *n ) == nonOverlappedMicroNodeIds->end( ) ) &&
//                 ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) == freeMicroNodeIds->end( ) ) &&
//                 ( std::find( ghostMicroNodeIds->begin( ), ghostMicroNodeIds->end( ), *n ) == ghostMicroNodeIds->end( ) ) ){
            if ( ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) == freeMicroNodeIds->end( ) ) &&
                 ( std::find( ghostMicroNodeIds->begin( ), ghostMicroNodeIds->end( ), *n ) == ghostMicroNodeIds->end( ) ) ){

                results << "test_initializeIncrement (test 18) & False\n";
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

                results << "test_initializeIncrement (test 19) & False\n";
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

                results << "test_initializeIncrement (test 20) & False\n";
                return 1;

            }

        }

    }

    const DOFMap *microGlobalToLocalDOFMap = reader.getMicroGlobalToLocalDOFMap( );

    if ( microGlobalToLocalDOFMap->size( ) != ( freeMicroNodeIds->size( ) + ghostMicroNodeIds->size( ) ) ){

        results << "test_initializeIncrement (test 21) & False\n";
        return 1;

    }

    for ( auto n  = freeMicroNodeIds->begin( );
               n != freeMicroNodeIds->end( );
               n++ ){

        if ( microGlobalToLocalDOFMap->find( *n ) == microGlobalToLocalDOFMap->end( ) ){

            results << "test_initializeIncrement (test 22) & False\n";
            return 1;

        }

    }

    for ( auto n  = ghostMicroNodeIds->begin( );
               n != ghostMicroNodeIds->end( );
               n++ ){

        if ( microGlobalToLocalDOFMap->find( *n ) == microGlobalToLocalDOFMap->end( ) ){

            results << "test_initializeIncrement (test 23) & False\n";
            return 1;

        }

    }

    const DOFMap *macroGlobalToLocalDOFMap = reader.getMacroGlobalToLocalDOFMap( );

    if ( macroGlobalToLocalDOFMap->size( ) != ( freeMacroNodeIds->size( ) + ghostMacroNodeIds->size( ) ) ){

        results << "test_initializeIncrement (test 24) & False\n";
        return 1;

    }

    for ( auto n  = freeMacroNodeIds->begin( );
               n != freeMacroNodeIds->end( );
               n++ ){

        if ( macroGlobalToLocalDOFMap->find( *n ) == macroGlobalToLocalDOFMap->end( ) ){

            results << "test_initializeIncrement (test 25) & False\n";
            return 1;

        }

    }

    for ( auto n  = ghostMacroNodeIds->begin( );
               n != ghostMacroNodeIds->end( );
               n++ ){

        if ( macroGlobalToLocalDOFMap->find( *n ) == macroGlobalToLocalDOFMap->end( ) ){

            results << "test_initializeIncrement (test 26) & False\n";
            return 1;

        }

    }

    const floatVector *macroDispDOFVectorResult = reader.getMacroDispDOFVector( );

    const floatVector macroDispDOFVectorAnswer =
        {
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0., 0., -0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        };
   
    if ( !vectorTools::fuzzyEquals( macroDispDOFVectorAnswer, *macroDispDOFVectorResult ) ){
        results << "test_initializeIncrement (test 27) & False\n";
        return 1;
    }


    if ( reader.microBodyForceDefined( ) ){

        results << "test_initializeIncrement (test 29) & False\n";
        return 1;

    }

    if ( !reader.microAccelerationDefined( ) ){

        results << "test_initializeIncrement (test 31) & False\n";
        return 1;

    }

    if ( reader.useReconstructedMassCenters( ) ){

        results << "test_initializeIncrement (test 32) & False\n";
        return 1;

    }

    if ( !reader.microVelocitiesDefined( ) ){

        results << "test_initializeIncrement (test 34) & False\n";
        return 1;

    }

    const floatVector macroAccelerationsAnswer = { 0., 0., -0.004, 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    const floatVector *macroAccelerationsResult = reader.getMacroAccelerations( );

    for ( auto mAR = macroAccelerationsResult->begin( ); mAR != macroAccelerationsResult->end( ); mAR++ ){

        if ( !vectorTools::fuzzyEquals( *mAR, macroAccelerationsAnswer[ ( mAR - macroAccelerationsResult->begin( ) ) % 12 ], 1e-5, 1e-5 ) ){
    
            results << "test_initializeIncrement (test 36) & False\n";
            return 1;
    
        }

    }

    if ( !reader.macroAccelerationDefined( ) ){

        results << "test_initializeIncrement (test 37) & False\n";
        return 1;

    }

    const floatVector macroVelocitiesAnswer = { 0., 0., -0.002, 0., 0., 0., 0., 0., 0., 0., 0., 0. };

    const floatVector *macroVelocitiesResult = reader.getMacroVelocities( );

    for ( auto mVR = macroVelocitiesResult->begin( ); mVR != macroVelocitiesResult->end( ); mVR++ ){

        if ( !vectorTools::fuzzyEquals( *mVR, macroVelocitiesAnswer[ ( mVR - macroVelocitiesResult->begin( ) ) % 12 ], 1e-5, 1e-5 ) ){

            std::cout << *mVR << "\n"; 
            results << "test_initializeIncrement (test 38) & False\n";
            return 1;
    
        }

    }

    if ( !reader.macroVelocitiesDefined( ) ){

        results << "test_initializeIncrement (test 39) & False\n";
        return 1;

    }


    if ( !reader.microInternalForceDefined( ) ){

        results << "test_initializeIncrement (test 41) & False\n";
        return 1;

    }

    const floatVector macroInternalForcesAnswer = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    const floatVector *macroInternalForcesResult = reader.getMacroInternalForces( );

    for ( auto mIFR = macroInternalForcesResult->begin( ); mIFR != macroInternalForcesResult->end( ); mIFR++ ){

        if ( !vectorTools::fuzzyEquals( *mIFR, macroInternalForcesAnswer[ ( mIFR - macroInternalForcesResult->begin( ) ) % 12 ], 1e-5, 1e-4 ) ){

            std::cout << mIFR - macroInternalForcesResult->begin( ) << ": " << *mIFR << "\n";
            results << "test_initializeIncrement (test 44) & False\n";
            return 1;

        }

    }

    if ( !reader.macroInternalForceDefined( ) ){

        results << "test_initializeIncrement (test 45) & False\n";
        return 1;

    }

    const floatVector macroInertialForcesAnswer = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    const floatVector *macroInertialForcesResult = reader.getMacroInertialForces( );

    for ( auto mIFR = macroInertialForcesResult->begin( ); mIFR != macroInertialForcesResult->end( ); mIFR++ ){

        if ( !vectorTools::fuzzyEquals( *mIFR, macroInertialForcesAnswer[ ( mIFR - macroInertialForcesResult->begin( ) ) % 12 ], 1e-5, 1e-4 ) ){

            std::cout << mIFR - macroInertialForcesResult->begin( ) << ": " << *mIFR << "\n";
            results << "test_initializeIncrement (test 46) & False\n";
            return 1;

        }

    }

    if ( !reader.macroInertialForceDefined( ) ){

        results << "test_initializeIncrement (test 47) & False\n";
        return 1;

    }

    const std::string macroReferenceDensityTypesAnswer = "constant";
    const floatVector macroReferenceDensitiesAnswer = { 2. };
    const std::unordered_map< unsigned int, floatVector > *macroReferenceDensitiesResult = reader.getMacroReferenceDensities( );
    const std::unordered_map< unsigned int, std::string > *macroReferenceDensityTypesResult = reader.getMacroReferenceDensityTypes( );

    for ( auto mRDR = macroReferenceDensitiesResult->begin( ); mRDR != macroReferenceDensitiesResult->end( ); mRDR++ ){

        if ( !vectorTools::fuzzyEquals( macroReferenceDensitiesAnswer, mRDR->second ) ){

            results << "test_initializeIncrement (test 48) & False\n";
            return 1;

        }

    }

    for ( auto mRDTR = macroReferenceDensityTypesResult->begin( ); mRDTR != macroReferenceDensityTypesResult->end( ); mRDTR++ ){

        if ( macroReferenceDensityTypesAnswer.compare( mRDTR->second ) != 0 ){

            results << "test_initializeIncrement (test 49) & False\n";
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

            results << "test_initializeIncrement (test 50) & False\n";
            return 1;

        }

    }

    for ( auto mRMITR  = macroReferenceMomentOfInertiaTypesResult->begin( );
               mRMITR != macroReferenceMomentOfInertiaTypesResult->end( ); mRMITR++ ){

        if ( macroReferenceMomentOfInertiaTypesAnswer.compare( mRMITR->second ) != 0 ){

            results << "test_initializeIncrement (test 51) & False\n";
            return 1;

        }

    }

    const floatVector macroExternalForcesAnswer = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    const floatVector *macroExternalForcesResult = reader.getMacroExternalForces( );

    for ( auto mIFR = macroExternalForcesResult->begin( ); mIFR != macroExternalForcesResult->end( ); mIFR++ ){

        if ( !vectorTools::fuzzyEquals( *mIFR, macroExternalForcesAnswer[ ( mIFR - macroExternalForcesResult->begin( ) ) % 12 ], 1e-5, 1e-4 ) ){

            std::cout << mIFR - macroExternalForcesResult->begin( ) << ": " << *mIFR << "\n";
            results << "test_initializeIncrement (test 52) & False\n";
            return 1;

        }

    }

    if ( reader.macroExternalForceDefined( ) ){

        results << "test_initializeIncrement (test 53) & False\n";
        return 1;

    }

    if ( reader.microSurfaceTractionDefined( ) ){

        results << "test_initializeIncrement (test 55) & False\n";
        return 1;

    }

    if ( reader.microExternalForceDefined( ) ){

        results << "test_initializeIncrement (test 57) & False\n";
        return 1;

    }


    const floatVector previousMacroDispDOFVectorAnswer = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    const floatVector *previousMacroDispDOFVectorResult = reader.getPreviousMacroDispDOFVector( );

    for ( auto v  = previousMacroDispDOFVectorResult->begin( );
               v != previousMacroDispDOFVectorResult->end( );
               v++ ){

        uIntType indx = v - previousMacroDispDOFVectorResult->begin( );

        if ( !vectorTools::fuzzyEquals( previousMacroDispDOFVectorAnswer[ indx % 12 ], *v, 1e-5, 1e-5 ) ){

            std::cout << *v << "\n";
            results << "test_initializeIncrement (test 61) & False\n";
            return 1;

        }

    }

    const floatVector previousMacroVelocitiesAnswer = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    const floatVector *previousMacroVelocitiesResult = reader.getPreviousMacroVelocities( );

    for ( auto v  = previousMacroVelocitiesResult->begin( );
               v != previousMacroVelocitiesResult->end( );
               v++ ){

        uIntType indx = v - previousMacroVelocitiesResult->begin( );

        if ( !vectorTools::fuzzyEquals( previousMacroVelocitiesAnswer[ indx % 12 ], *v, 1e-5, 1e-5 ) ){

            std::cout << *v << "\n";
            results << "test_initializeIncrement (test 62) & False\n";
            return 1;

        }

    }

    const floatVector previousMacroAccelerationsAnswer = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    const floatVector *previousMacroAccelerationsResult = reader.getPreviousMacroAccelerations( );

    for ( auto v  = previousMacroAccelerationsResult->begin( );
               v != previousMacroAccelerationsResult->end( );
               v++ ){

        uIntType indx = v - previousMacroAccelerationsResult->begin( );

        if ( !vectorTools::fuzzyEquals( previousMacroAccelerationsAnswer[ indx % 12 ], *v, 1e-5, 1e-5 ) ){

            std::cout << *v << "\n";
            results << "test_initializeIncrement (test 63) & False\n";
            return 1;

        }

    }

    if ( !reader.extractPreviousDOFValues( ) ){

        results << "test_initializeIncrement (test 64) & False\n";
        return 1;

    }

    const floatType microTimeAnswer = 1.;
    const floatType* microTimeResult = reader.getMicroTime( );

    if ( !vectorTools::fuzzyEquals( microTimeAnswer, *microTimeResult ) ){

        results << "test_initializeIncrement (test 64) & False\n";
        return 1;

    }

    const floatType macroTimeAnswer = 1.;
    const floatType* macroTimeResult = reader.getMacroTime( );

    if ( !vectorTools::fuzzyEquals( macroTimeAnswer, *macroTimeResult ) ){

        results << "test_initializeIncrement (test 65) & False\n";
        return 1;

    }

    const floatType DtAnswer = 1.;
    const floatType* DtResult = reader.getDt( );

    if ( !vectorTools::fuzzyEquals( DtAnswer, *DtResult ) ){

        results << "test_initializeIncrement (test 66) & False\n";
        return 1;

    }

    const floatType newmarkGammaAnswer = 0.50;
    const floatType newmarkBetaAnswer  = 0.25;

    if ( !vectorTools::fuzzyEquals( newmarkGammaAnswer, *reader.getNewmarkGamma( ) ) ){

        results << "test_initializeIncrement (test 67) & False\n";
        return 1;

    }

    if ( !vectorTools::fuzzyEquals( newmarkBetaAnswer, *reader.getNewmarkBeta( ) ) ){

        results << "test_initializeIncrement (test 68) & False\n";
        return 1;

    }

    results << "test_initializeIncrement & True\n";
    return 0;
}

int test_getFreeMicroDomainNames( std::ostream &results ){
    /*!
     * Test getting a pointer to the free micro domain names
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "../testFiles/testConfig.yaml";
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

    std::string filename = "../testFiles/testConfig.yaml";
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

//int test_getFreeMicroSurfaceNames( std::ostream &results ){
//    /*!
//     * Test getting a pointer to the free micro surface names
//     *
//     * :param std::ofstream &results: The output file
//     */
//
//    std::string filename = "../testFiles/testConfig.yaml";
//    inputFileProcessor::inputFileProcessor reader( filename );
//
//    if ( reader.getError( ) ){
//        reader.getError( )->print( );
//        results << "test_getFreeMicroSurfaceNames & False\n";
//        return 1;
//    }
//
//    stringVector answer = { "free_nodeset_surface_1",
//                            "free_nodeset_surface_2",
//                            "free_nodeset_surface_3",
//                            "free_nodeset_surface_4",
//                            "free_nodeset_surface_5",
//                            "free_nodeset_surface_6",
//                            "free_nodeset_surface_7",
//                            "free_nodeset_surface_8" };
//
//    const stringVector *result = reader.getFreeMicroSurfaceNames( );
//
//    unsigned int indx = 0;
//    for ( auto it = result->begin( ); it != result->end( ); it++ ){
//
//        if ( it->compare( answer[ indx ] ) != 0 ){
//
//            results << "test_getFreeMicroSurfaceNames (test 1) & False\n";
//            return 1;
//
//        }
//
//        indx++;
//
//    }
//
//    results << "test_getFreeMicroSurfaceNames & True\n";
//    return 0;
//}
//
//int test_getGhostMicroSurfaceNames( std::ostream &results ){
//    /*!
//     * Test getting a pointer to the ghost micro surface names
//     *
//     * :param std::ofstream &results: The output file
//     */
//
//    std::string filename = "../testFiles/testConfig.yaml";
//    inputFileProcessor::inputFileProcessor reader( filename );
//
//    if ( reader.getError( ) ){
//        reader.getError( )->print( );
//        results << "test_getGhostMicroSurfaceNames & False\n";
//        return 1;
//    }
//
//    stringVector answer = { "ghost_nodeset_surface_1",
//                            "ghost_nodeset_surface_2",
//                            "ghost_nodeset_surface_3",
//                            "ghost_nodeset_surface_4",
//                            "ghost_nodeset_surface_5",
//                            "ghost_nodeset_surface_6",
//                            "ghost_nodeset_surface_7",
//                            "ghost_nodeset_surface_8" };
//
//    const stringVector *result = reader.getGhostMicroSurfaceNames( );
//
//    unsigned int indx = 0;
//    for ( auto it = result->begin( ); it != result->end( ); it++ ){
//
//        if ( it->compare( answer[ indx ] ) != 0 ){
//
//            results << "test_getGhostMicroSurfaceNames (test 1) & False\n";
//            return 1;
//
//        }
//
//        indx++;
//
//    }
//
//    results << "test_getGhostMicroDomainNames & True\n";
//    return 0;
//}
//
//int test_getNonOverlappedMicroSurfaceNames( std::ostream &results ){
//    /*!
//     * Test getting a pointer to the non-overlapped micro surface names
//     *
//     * :param std::ofstream &results: The output file
//     */
//
//    std::string filename = "../testFiles/testConfig.yaml";
//    inputFileProcessor::inputFileProcessor reader( filename );
//
//    if ( reader.getError( ) ){
//        reader.getError( )->print( );
//        results << "test_getNonOverlappedMicroDomainNames & False\n";
//        return 1;
//    }
//
//    stringVector answer = { "non_overlapped_nodeset_surface_1",
//                            "non_overlapped_nodeset_surface_2",
//                            "non_overlapped_nodeset_surface_3",
//                            "non_overlapped_nodeset_surface_4",
//                            "non_overlapped_nodeset_surface_5",
//                            "non_overlapped_nodeset_surface_6",
//                            "non_overlapped_nodeset_surface_7",
//                            "non_overlapped_nodeset_surface_8" };
//
//    const stringVector *result = reader.getNonOverlappedMicroSurfaceNames( );
//
//    unsigned int indx = 0;
//    for ( auto it = result->begin( ); it != result->end( ); it++ ){
//
//        if ( it->compare( answer[ indx ] ) != 0 ){
//
//            results << "test_getNonOverlappedMicroSurfaceNames (test 1) & False\n";
//            return 1;
//
//        }
//
//        indx++;
//
//    }
//
//    results << "test_getNonOverlappedMicroDomainNames & True\n";
//    return 0;
//
//}

int test_getCouplingInitialization( std::ofstream &results ){
    /*!
     * Test getting the coupling initialization from the configuration file
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "../testFiles/testConfig.yaml";
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

    std::string projectionTypeAnswer = "direct_projection";
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

    results << "test_getCouplingInitialization & True\n";
    return 0;

}

int test_getVolumeReconstructionConfig( std::ofstream &results ){
    /*!
     * Test getting the volume reconstruction configuration from the configuration file
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "../testFiles/testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_getVolumeReconstructionConfig & False\n";
        return 1;
    }

    const YAML::Node couplingInitialization = reader.getVolumeReconstructionConfig( );

    if ( !couplingInitialization ){
        results << "test_getVolumeReconstructionConfig (test 1) & False\n";
        return 1;
    }

    std::string typeAnswer = "dual_contouring";
    if ( couplingInitialization[ "type" ].as<std::string>( ).compare( typeAnswer ) ){
        results << "test_getVolumeReconstructionConfig (test 2) & False\n";
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

    std::string filename = "../testFiles/testConfig.yaml";
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

    std::string filename = "../testFiles/testConfig.yaml";
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

int test_getFreeMicroSurfaceApproximateSplitCount( std::ostream &results ){
    /*!
     * Test getting a pointer to the approximate number of surfaces to split a micro
     * domain into.
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "../testFiles/testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_getFreeMicroSurfaceApproximateSplitCount & False\n";
        return 1;
    }

    uIntVector answer( 8, 6 );

    const uIntVector *result = reader.getFreeMicroSurfaceApproximateSplitCount( );

    unsigned int indx = 0;
    for ( auto it = result->begin( ); it != result->end( ); it++ ){

        if ( !vectorTools::fuzzyEquals( *it, answer[ indx ] ) ){

            results << "test_getFreeMicroSurfaceApproximateSplitCount (test 1) & False\n";
            return 1;

        }

        indx++;

    }

    results << "test_getFreeMicroSurfaceApproximateSplitCount & True\n";
    return 0;
}

int test_getGhostMicroSurfaceApproximateSplitCount( std::ostream &results ){
    /*!
     * Test getting a pointer to the approximate number of surfaces to split a micro
     * domain into.
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "../testFiles/testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_getGhostMicroSurfaceApproximateSplitCount & False\n";
        return 1;
    }

    uIntVector answer( 8, 6 );

    const uIntVector *result = reader.getGhostMicroSurfaceApproximateSplitCount( );

    unsigned int indx = 0;
    for ( auto it = result->begin( ); it != result->end( ); it++ ){

        if ( !vectorTools::fuzzyEquals( *it, answer[ indx ] ) ){

            results << "test_getGhostMicroSurfaceApproximateSplitCount (test 1) & False\n";
            return 1;

        }

        indx++;

    }

    results << "test_getGhostMicroSurfaceApproximateSplitCount & True\n";
    return 0;
}

int test_outputReferenceInformation( std::ofstream &results ){
    /*!
     * Test whether the reference information should be output
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "../testFiles/testConfig.yaml";
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

    std::string filename = "../testFiles/testConfig.yaml";
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

    std::string filename = "../testFiles/testConfig.yaml";
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
//    test_getFreeMicroDomainNames( results );
//    test_getGhostMicroDomainNames( results );
////    test_getFreeMicroSurfaceNames( results );
////    test_getGhostMicroSurfaceNames( results );
//    test_getFreeMacroDomainNames( results );
//    test_getGhostMacroDomainNames( results );
////    test_getNonOverlappedMicroSurfaceNames( results );
//    test_getCouplingInitialization( results );
//    test_getVolumeReconstructionConfig( results );
//    test_getFreeMicroSurfaceApproximateSplitCount( results );
//    test_getGhostMicroSurfaceApproximateSplitCount( results );
//    test_outputReferenceInformation( results );
//    test_outputHomogenizedInformation( results );
//    test_outputUpdatedDOF( results );
//
    //Close the results file
    results.close();

}
