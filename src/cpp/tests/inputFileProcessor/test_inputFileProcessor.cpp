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
typedef inputFileProcessor::uIntType uIntType; //!Define an unsigned integer type
typedef inputFileProcessor::uIntVector uIntVector; //!Define a vector of unsigned ints
typedef inputFileProcessor::stringVector stringVector; //!Define a vector of strings
typedef inputFileProcessor::DOFMap DOFMap; //!Define the map between DOF values

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

    errorOut error = reader._macroscale->readMesh( 1, resultMacroNodes );

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

    errorOut error = reader.initializeIncrement( 1, 1 );
    if ( error ){
        error->print( );
        results << "test_initializeIncrement & False\n";
        return 1;
    }

    floatType resultAnswer = 2.;

    const floatVector *densityResult = reader.getMicroDensities( );

    for ( auto it = densityResult->begin( ); it != densityResult->end( ); it++ ){

        if ( !vectorTools::fuzzyEquals( *it, resultAnswer ) ){ //NOTE: Tolerance due to inertial effects
            std::cout << *it << "\n";
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

    floatVector answer = { 0., 0., 0.001 };

    const floatVector *displacementResult = reader.getMicroDisplacements( );

    unsigned int itmp = 0;

    for ( auto it = displacementResult->begin( ); it != displacementResult->end( ); it++ ){

        if ( !vectorTools::fuzzyEquals( *it, answer[ itmp ] ) ){ //Note tolerance due to inertia term

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

    const uIntVector macroNodeReferenceConnectivityCellIndicesAnswer = { 0, 9, 18 };

    const uIntVector *macroNodeReferenceConnectivityCellIndicesResult = reader.getMacroNodeReferenceConnectivityCellIndices( );

    if ( !vectorTools::fuzzyEquals( macroNodeReferenceConnectivityCellIndicesAnswer, *macroNodeReferenceConnectivityCellIndicesResult ) ){
        results << "test_initializeIncrement (test 7) & False\n";
        return 1;
    }

    const uIntVector freeMacroCellIdsAnswer = { 1 };
    const uIntVector ghostMacroCellIdsAnswer = { 2 };

    const uIntVector *freeMacroCellIdsResult = reader.getFreeMacroCellIds( );
    const uIntVector *ghostMacroCellIdsResult = reader.getGhostMacroCellIds( );

    if ( !vectorTools::fuzzyEquals( freeMacroCellIdsAnswer, *freeMacroCellIdsResult ) ){
        results << "test_initializeIncrement (test 8) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( ghostMacroCellIdsAnswer, *ghostMacroCellIdsResult ) ){
        results << "test_initializeIncrement (test 9) & False\n";
        return 1;
    }

    const uIntVector freeMacroCellMicroDomainCountsAnswer = { 8 };
    const uIntVector ghostMacroCellMicroDomainCountsAnswer = { 8 };

    const uIntVector *freeMacroCellMicroDomainCountsResult = reader.getFreeMacroCellMicroDomainCounts( );
    const uIntVector *ghostMacroCellMicroDomainCountsResult = reader.getGhostMacroCellMicroDomainCounts( );

    if ( !vectorTools::fuzzyEquals( freeMacroCellMicroDomainCountsAnswer, *freeMacroCellMicroDomainCountsResult ) ){
        results << "test_initializeIncrement (test 10) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( ghostMacroCellMicroDomainCountsAnswer, *ghostMacroCellMicroDomainCountsResult ) ){
        results << "test_initializeIncrement (test 11) & False\n";
        return 1;
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

    const uIntVector *nonOverlappedMicroNodeIds = reader.getNonOverlappedMicroNodeIds( );
    const uIntVector *freeMicroNodeIds = reader.getFreeMicroNodeIds( );
    const uIntVector *ghostMicroNodeIds = reader.getGhostMicroNodeIds( );

    for ( auto n  = freeMicroNodeIds->begin( );
               n != freeMicroNodeIds->end( );
               n++ ){

        if ( std::find( nonOverlappedMicroNodeIds->begin( ), nonOverlappedMicroNodeIds->end( ), *n ) != nonOverlappedMicroNodeIds->end( ) ){
            std::cout << "*n: " << *n << "\n";
            results << "test_initializeIncrement (test 13) & False\n";
            return 1;
        }
    }

    for ( auto n  = ghostMicroNodeIds->begin( );
               n != ghostMicroNodeIds->end( );
               n++ ){

        if ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) != freeMicroNodeIds->end( ) ){
            std::cout << "*n: " << *n << "\n";
            results << "test_initializeIncrement (test 14) & False\n";
            return 1;
        }

        if ( std::find( nonOverlappedMicroNodeIds->begin( ), nonOverlappedMicroNodeIds->end( ), *n ) != nonOverlappedMicroNodeIds->end( ) ){
            std::cout << "*n: " << *n << "\n";
            results << "test_initializeIncrement (test 15) & False\n";
            return 1;
        }

    }

    const stringVector *nonOverlappedMicroSurfaceNames = reader.getNonOverlappedMicroSurfaceNames( );
    uIntVector nodes;

    for ( auto surface  = nonOverlappedMicroSurfaceNames->begin( );
               surface != nonOverlappedMicroSurfaceNames->end( );
               surface++ ){

        reader._microscale->getSubDomainNodes( 0, *surface, nodes );

        for ( auto n = nodes.begin( ); n != nodes.end( ); n++ ){

            if ( std::find( nonOverlappedMicroNodeIds->begin( ), nonOverlappedMicroNodeIds->end( ), *n ) == nonOverlappedMicroNodeIds->end( ) ){
                results << "test_initializeIncrement (test 16) & False\n";
                return 1;
            }

        }

    }

    const stringVector *freeMicroDomainNames = reader.getFreeMicroDomainNames( );
    for ( auto domain  = freeMicroDomainNames->begin( );
               domain != freeMicroDomainNames->end( );
               domain++ ){

        reader._microscale->getSubDomainNodes( 0, *domain, nodes );

        for ( auto n = nodes.begin( ); n != nodes.end( ); n++ ){

            if ( ( std::find( nonOverlappedMicroNodeIds->begin( ), nonOverlappedMicroNodeIds->end( ), *n ) == nonOverlappedMicroNodeIds->end( ) ) &&
                 ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) == freeMicroNodeIds->end( ) ) ){

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

            if ( ( std::find( nonOverlappedMicroNodeIds->begin( ), nonOverlappedMicroNodeIds->end( ), *n ) == nonOverlappedMicroNodeIds->end( ) ) &&
                 ( std::find( freeMicroNodeIds->begin( ), freeMicroNodeIds->end( ), *n ) == freeMicroNodeIds->end( ) ) &&
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

            if ( ( std::find( ghostMacroNodeIds->begin( ), ghostMacroNodeIds->end( ), *n ) == ghostMicroNodeIds->end( ) ) &&
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

    const floatVector microBodyForcesAnswer = { 0., 0., 0. };
    const floatVector *microBodyForcesResult = reader.getMicroBodyForces( );

    if ( !vectorTools::fuzzyEquals( *microBodyForcesResult, microBodyForcesAnswer ) ){

        results << "test_initializeIncrement (test 28) & False\n";
        return 1;

    }

    if ( reader.microBodyForceDefined( ) ){

        results << "test_initializeIncrement (test 29) & False\n";
        return 1;

    }

    const floatVector microAccelerationsAnswer = { 0., 0., 0.004 };
    const floatVector *microAccelerationsResult = reader.getMicroAccelerations( );

    for ( auto mAR = microAccelerationsResult->begin( ); mAR != microAccelerationsResult->end( ); mAR++ ){

        if ( !vectorTools::fuzzyEquals( *mAR, microAccelerationsAnswer[ ( mAR - microAccelerationsResult->begin( ) ) % 3 ], 1e-5, 1e-5 ) ){
    
            results << "test_initializeIncrement (test 30) & False\n";
            return 1;
    
        }

    }

    if ( !reader.microAccelerationDefined( ) ){

        results << "test_initializeIncrement (test 31) & False\n";
        return 1;

    }

    if ( reader.useReconstructedMassCenters( ) ){

        results << "test_initializeIncrement (test 32) & False\n";
        return 1;

    }

    const floatVector microVelocitiesAnswer = { 0., 0., 0.002 };

    const floatVector *microVelocitiesResult = reader.getMicroVelocities( );

    for ( auto mVR = microVelocitiesResult->begin( ); mVR != microVelocitiesResult->end( ); mVR++ ){

        if ( !vectorTools::fuzzyEquals( *mVR, microVelocitiesAnswer[ ( mVR - microVelocitiesResult->begin( ) ) % 3 ], 1e-5, 1e-5 ) ){
    
            results << "test_initializeIncrement (test 33) & False\n";
            return 1;
    
        }

    }

    if ( !reader.microVelocitiesDefined( ) ){

        results << "test_initializeIncrement (test 34) & False\n";
        return 1;

    }

    unsigned int dim = 3;
    const floatVector microStressesAnswer( dim * dim * reader.getMicroDensities( )->size( ) );
    const floatVector *microStressesResult = reader.getMicroStresses( );

    if ( !vectorTools::fuzzyEquals( *microStressesResult, microStressesAnswer, 1e-5, 1e-1 ) ){ //NOTE: The tolerance is high because stresses can be large

        vectorTools::print( *microStressesResult );
        results << "test_initializeIncrement (test 35) & False\n";
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

    const floatVector microInternalForcesAnswer = { 0., 0., 0. };
    const floatVector *microInternalForcesResult = reader.getMicroInternalForces( );

    for ( auto mIFR = microInternalForcesResult->begin( ); mIFR != microInternalForcesResult->end( ); mIFR++ ){

        if ( !vectorTools::fuzzyEquals( *mIFR, microInternalForcesAnswer[ ( mIFR - microInternalForcesResult->begin( ) ) % 3 ], 1e-5, 1e-4 ) ){

            std::cout << mIFR - microInternalForcesResult->begin( ) << ": " << *mIFR << "\n";
            results << "test_initializeIncrement (test 40) & False\n";
            return 1;

        }

    }

    if ( !reader.microInternalForceDefined( ) ){

        results << "test_initializeIncrement (test 41) & False\n";
        return 1;

    }

    const floatVector microInertialForcesAnswer = { 0., 0., 0. };
    const floatVector *microInertialForcesResult = reader.getMicroInertialForces( );

    for ( auto mIFR = microInertialForcesResult->begin( ); mIFR != microInertialForcesResult->end( ); mIFR++ ){

        if ( !vectorTools::fuzzyEquals( *mIFR, microInertialForcesAnswer[ ( mIFR - microInertialForcesResult->begin( ) ) % 3 ], 1e-5, 1e-4 ) ){

            std::cout << mIFR - microInertialForcesResult->begin( ) << ": " << *mIFR << "\n";
            results << "test_initializeIncrement (test 42) & False\n";
            return 1;

        }

    }

    if ( !reader.microInertialForceDefined( ) ){

        results << "test_initializeIncrement (test 43) & False\n";
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

    const floatVector microSurfaceTractionsAnswer = { 0., 0., 0. };
    const floatVector *microSurfaceTractionsResult = reader.getMicroSurfaceTractions( );

    if ( !vectorTools::fuzzyEquals( *microSurfaceTractionsResult, microSurfaceTractionsAnswer ) ){

        vectorTools::print( *microSurfaceTractionsResult );
        results << "test_initializeIncrement (test 54) & False\n";
        return 1;

    }

    if ( reader.microSurfaceTractionDefined( ) ){

        results << "test_initializeIncrement (test 55) & False\n";
        return 1;

    }

    const floatVector microExternalForcesAnswer = { 0., 0., 0. };
    const floatVector *microExternalForcesResult = reader.getMicroExternalForces( );

    if ( !vectorTools::fuzzyEquals( *microExternalForcesResult, microExternalForcesAnswer ) ){

        vectorTools::print( *microExternalForcesResult );
        results << "test_initializeIncrement (test 56) & False\n";
        return 1;

    }

    if ( reader.microExternalForceDefined( ) ){

        results << "test_initializeIncrement (test 57) & False\n";
        return 1;

    }

    const floatVector previousMicroDisplacementsAnswer = { 0., 0., 0. };
    const floatVector *previousMicroDisplacementsResult = reader.getPreviousMicroDisplacements( );

    for ( auto v  = previousMicroDisplacementsResult->begin( );
               v != previousMicroDisplacementsResult->end( );
               v++ ){

        uIntType indx = v - previousMicroDisplacementsResult->begin( );

        if ( !vectorTools::fuzzyEquals( previousMicroDisplacementsAnswer[ indx % 3 ], *v, 1e-5, 1e-5 ) ){

            std::cout << *v << "\n";
            results << "test_initializeIncrement (test 58) & False\n";
            return 1;

        }

    }

    const floatVector previousMicroVelocitiesAnswer = { 0., 0., 0. };
    const floatVector *previousMicroVelocitiesResult = reader.getPreviousMicroVelocities( );

    for ( auto v  = previousMicroVelocitiesResult->begin( );
               v != previousMicroVelocitiesResult->end( );
               v++ ){

        uIntType indx = v - previousMicroVelocitiesResult->begin( );

        if ( !vectorTools::fuzzyEquals( previousMicroVelocitiesAnswer[ indx % 3 ], *v, 1e-5, 1e-5 ) ){

            std::cout << *v << "\n";
            results << "test_initializeIncrement (test 59) & False\n";
            return 1;

        }

    }

    const floatVector previousMicroAccelerationsAnswer = { 0., 0., 0.0 };
    const floatVector *previousMicroAccelerationsResult = reader.getPreviousMicroAccelerations( );

    for ( auto v  = previousMicroAccelerationsResult->begin( );
               v != previousMicroAccelerationsResult->end( );
               v++ ){

        uIntType indx = v - previousMicroAccelerationsResult->begin( );

        if ( !vectorTools::fuzzyEquals( previousMicroAccelerationsAnswer[ indx % 3 ], *v, 1e-5, 1e-5 ) ){

            std::cout << *v << "\n";
            results << "test_initializeIncrement (test 60) & False\n";
            return 1;

        }

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

    if ( couplingInitialization[ "output_updated_dof" ][ "macroscale_filename" ].as< std::string >( ).compare( "macroscale_dof" ) != 0 ){

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
    test_getFreeMacroDomainNames( results );
    test_getGhostMacroDomainNames( results );
    test_getNonOverlappedMicroSurfaceNames( results );
    test_getCouplingInitialization( results );
    test_getVolumeReconstructionConfig( results );
    test_getFreeMicroSurfaceApproximateSplitCount( results );
    test_getGhostMicroSurfaceApproximateSplitCount( results );
    test_outputReferenceInformation( results );
    test_outputHomogenizedInformation( results );
    test_outputUpdatedDOF( results );

    //Close the results file
    results.close();
}
