//!The test file for inputFileProcessor.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>

#include<inputFileProcessor.h>

typedef inputFileProcessor::errorNode errorNode; //!Redefinition for the error node
typedef inputFileProcessor::errorOut errorOut; //!Redefinition for a pointer to the error node
typedef inputFileProcessor::floatType floatType; //!Define the float values type.
typedef inputFileProcessor::floatVector floatVector; //! Define a vector of floats
typedef inputFileProcessor::floatMatrix floatMatrix; //!Define a matrix of floats
typedef inputFileProcessor::uIntVector uIntVector; //!Define a vector of unsigned ints
typedef inputFileProcessor::stringVector stringVector; //!Define a vector of strings

int test_openConfigurationFile( std::ofstream &results ){
    /*!
     * Test opening the YAML configuration file.
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "../testFiles/testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print();
        results << "test_openConfigurationFile (test 1) & False\n";
        return 1;
    }

    reader = inputFileProcessor::inputFileProcessor( );
    errorOut error = reader.setConfigurationFilename( "" );
    if ( !error ){
        results << "test_openConfigurationFile (test 2) & False\n";
        return 1;
    }

    error = reader.setConfigurationFilename( filename );
    if ( error ){
        error->print();
        results << "test_openConfigurationFile (test 3) & False\n";
        return 1;
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

    std::string filename = "../testFiles/testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader;

    errorOut error = reader.setConfigurationFilename( "" );

    if ( !error ){
        results << "test_setConfigurationFile & False\n";
        return 1;
    }

    error = reader.setConfigurationFilename( filename );

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

    std::string filename = "../testFiles/testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print();
        results << "test_initializeFileInterfaces & False\n";
        return 1;
    }

    floatVector answerMacroNodes =
        {
            1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 
            1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 2, 1, 1, 2, 0, 0, 
            2, 1, 0, 2, 0, 0, 3, 0, 1, 3, 1, 1, 3, 1, 0, 3
        };

    floatVector answerMicroNodes;

    floatVector resultMacroNodes, resultMicroNodes;

    errorOut error = reader._macroscale->readMesh( 10, resultMacroNodes );

    if ( error ){
        error->print();
        results << "test_initializeFileInterfaces & False\n";
        return 1;
    }

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

    if ( resultMicroNodes.size() != 80703 ){
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

    std::string filename = "../testFiles/testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_initializeIncrement & False\n";
        return 1;
    }

    errorOut error = reader.initializeIncrement( 1 );
    if ( error ){
        error->print( );
        results << "test_initializeIncrement & False\n";
        return 1;
    }

    floatType resultAnswer = 2000.;

    const floatVector *densityResult = reader.getMicroDensities( );

    for ( auto it = densityResult->begin( ); it != densityResult->end( ); it++ ){

        if ( !vectorTools::fuzzyEquals( *it, resultAnswer ) ){
            results << "test_initializeIncrement (test 1) & False\n";
            return 1;
        }

    }

    floatType volumeTotalAnswer = 3;

    const floatVector *volumeResult = reader.getMicroVolumes( );
    floatType volumeResultTotal = 0;
    
    for ( auto it = volumeResult->begin( ); it != volumeResult->end( ); it++ ){

        volumeResultTotal += *it;

    }

    if ( !vectorTools::fuzzyEquals( volumeTotalAnswer, volumeResultTotal ) ){
        results << "test_initializeIncrement (test 2) & False\n";
        return 1;
    }

    floatVector answer = { 0., 0., 0.01 };

    const floatVector *displacementResult = reader.getMicroDisplacements( );

    unsigned int itmp = 0;

    for ( auto it = displacementResult->begin( ); it != displacementResult->end( ); it++ ){

        if ( !vectorTools::fuzzyEquals( *it, answer[ itmp ] ) ){

            std::cout << it - displacementResult->begin( ) << "\n";
            std::cout << *it << "\n";
            std::cout << answer[ itmp ] << "\n";
            
            results << "test_initializeIncrement (test 3) & False\n";
            return 1;

        }
        
        itmp++;
        if ( itmp >= 3 ){

            itmp = 0;

        }

    }

    const floatVector *microNodeReferencePositionResult = reader.getMicroNodeReferencePositions( );

    if ( !vectorTools::fuzzyEquals( microNodeReferencePositionResult->size( ), displacementResult->size( ) ) ){
        results << "test_initializeIncrement (test 4) & False\n";
        return 1;
    }

    const floatVector macroNodeReferencePositionsAnswer = { 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1,
                                                            1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1,
                                                            0, 1, 2, 1, 1, 2, 0, 0, 2, 1, 0, 2,
                                                            0, 0, 3, 0, 1, 3, 1, 1, 3, 1, 0, 3, };

    const floatVector *macroNodeReferencePositionsResult = reader.getMacroNodeReferencePositions( );

    if ( !vectorTools::fuzzyEquals( macroNodeReferencePositionsAnswer, *macroNodeReferencePositionsResult ) ){
        results << "test_initializeIncrement (test 5) & False\n";
        return 1;
    }

    const uIntVector macroNodeReferenceConnectivityAnswer = { 9,  0,  1,  2,  3,  4, 5, 6,  7,
                                                              9,  8,  7,  4,  9, 10, 3, 0, 11,
                                                              9, 12, 13, 14, 15, 10, 8, 9, 11 };

    const uIntVector *macroNodeReferenceConnectivityResult = reader.getMacroNodeReferenceConnectivity( );

    if ( !vectorTools::fuzzyEquals( macroNodeReferenceConnectivityAnswer, *macroNodeReferenceConnectivityResult ) ){
        results << "test_initializeIncrement (test 6) & False\n";
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

int test_getFreeMicroSurfaceNames( std::ostream &results ){
    /*!
     * Test getting a pointer to the free micro surface names
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "../testFiles/testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_getFreeMicroSurfaceNames & False\n";
        return 1;
    }

    stringVector answer = { "free_nodeset_surface_1",
                            "free_nodeset_surface_2",
                            "free_nodeset_surface_3",
                            "free_nodeset_surface_4",
                            "free_nodeset_surface_5",
                            "free_nodeset_surface_6",
                            "free_nodeset_surface_7",
                            "free_nodeset_surface_8" };

    const stringVector *result = reader.getFreeMicroSurfaceNames( );

    unsigned int indx = 0;
    for ( auto it = result->begin( ); it != result->end( ); it++ ){

        if ( it->compare( answer[ indx ] ) != 0 ){

            results << "test_getFreeMicroSurfaceNames (test 1) & False\n";
            return 1;

        }

        indx++;

    }

    results << "test_getFreeMicroSurfaceNames & True\n";
    return 0;
}

int test_getGhostMicroSurfaceNames( std::ostream &results ){
    /*!
     * Test getting a pointer to the ghost micro surface names
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "../testFiles/testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_getGhostMicroSurfaceNames & False\n";
        return 1;
    }

    stringVector answer = { "ghost_nodeset_surface_1",
                            "ghost_nodeset_surface_2",
                            "ghost_nodeset_surface_3",
                            "ghost_nodeset_surface_4",
                            "ghost_nodeset_surface_5",
                            "ghost_nodeset_surface_6",
                            "ghost_nodeset_surface_7",
                            "ghost_nodeset_surface_8" };

    const stringVector *result = reader.getGhostMicroSurfaceNames( );

    unsigned int indx = 0;
    for ( auto it = result->begin( ); it != result->end( ); it++ ){

        if ( it->compare( answer[ indx ] ) != 0 ){

            results << "test_getGhostMicroSurfaceNames (test 1) & False\n";
            return 1;

        }

        indx++;

    }

    results << "test_getGhostMicroDomainNames & True\n";
    return 0;
}

int test_getNonOverlappedMicroSurfaceNames( std::ostream &results ){
    /*!
     * Test getting a pointer to the non-overlapped micro surface names
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "../testFiles/testConfig.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    if ( reader.getError( ) ){
        reader.getError( )->print( );
        results << "test_getNonOverlappedMicroDomainNames & False\n";
        return 1;
    }

    stringVector answer = { "non_overlapped_nodeset_surface_1",
                            "non_overlapped_nodeset_surface_2",
                            "non_overlapped_nodeset_surface_3",
                            "non_overlapped_nodeset_surface_4",
                            "non_overlapped_nodeset_surface_5",
                            "non_overlapped_nodeset_surface_6",
                            "non_overlapped_nodeset_surface_7",
                            "non_overlapped_nodeset_surface_8" };

    const stringVector *result = reader.getNonOverlappedMicroSurfaceNames( );

    unsigned int indx = 0;
    for ( auto it = result->begin( ); it != result->end( ); it++ ){

        if ( it->compare( answer[ indx ] ) != 0 ){

            results << "test_getNonOverlappedMicroSurfaceNames (test 1) & False\n";
            return 1;

        }

        indx++;

    }

    results << "test_getNonOverlappedMicroDomainNames & True\n";
    return 0;

}

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

    results << "test_getCouplingInitialization & True\n";
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

    test_openConfigurationFile( results );
    test_setConfigurationFile( results );
    test_initializeFileInterfaces( results );
    test_initializeIncrement( results );
    test_getFreeMicroDomainNames( results );
    test_getGhostMicroDomainNames( results );
    test_getFreeMicroSurfaceNames( results );
    test_getGhostMicroSurfaceNames( results );
    test_getNonOverlappedMicroSurfaceNames( results );
    test_getCouplingInitialization( results );

    //Close the results file
    results.close();
}
