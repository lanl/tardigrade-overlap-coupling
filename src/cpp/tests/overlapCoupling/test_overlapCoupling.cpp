//!The test file for dataFileInterface.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>
#define USE_EIGEN
#include<vector_tools.h>

#include<overlapCoupling.h>

typedef overlapCoupling::errorNode errorNode; //!Redefinition for the error node
typedef overlapCoupling::errorOut errorOut; //!Redefinition for a pointer to the error node
typedef overlapCoupling::floatType floatType; //!Define the float values type.
typedef overlapCoupling::floatVector floatVector; //! Define a vector of floats
typedef overlapCoupling::floatMatrix floatMatrix; //!Define a matrix of floats
typedef overlapCoupling::uIntVector uIntVector; //!Define a vector of unsigned ints

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

    //Close the results file
    results.close();
}
