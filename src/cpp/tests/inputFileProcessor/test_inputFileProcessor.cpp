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

int test_openConfigurationFile( std::ofstream &results ){
    /*!
     * Test opening the YAML configuration file.
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testFile.yaml";
    inputFileProcessor::inputFileProcessor reader( filename );

    errorOut error = reader.openConfigurationFile( );

    if ( error ){
        error->print();
        results << "test_openConfigurationFile (test 1) & False\n";
        return 1;
    }

    reader = inputFileProcessor::inputFileProcessor( );
    error = reader.openConfigurationFile( );
    if ( !error ){
        results << "test_openConfigurationFile (test 2) & False\n";
        return 1;
    }

    error = reader.openConfigurationFile( filename );
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

    std::string filename = "testFile.yaml";
    inputFileProcessor::inputFileProcessor reader;

    errorOut error = reader.openConfigurationFile( );

    if ( !error ){
        results << "test_setConfigurationFile & False\n";
        return 1;
    }

    reader.setConfigurationFilename( filename );
    error = reader.openConfigurationFile( );

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

    errorOut error = reader.openConfigurationFile( );
    if ( error ){
        error->print();
        results << "test_initializeFileInterfaces & False\n";
        return 1;
    }

    error = reader.initializeFileInterfaces( );

    if( error ){
        error->print();
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

    error = reader._macroscale->readMesh( 10, resultMacroNodes );

    if ( error ){
        error->print();
        results << "test_initializeFileInterfaces & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answerMacroNodes, resultMacroNodes ) ){
        results << "test_initializeFileInterfaces (test 1) & False\n";
        return 1;
    }

    error = reader._microscale->readMesh( 10, resultMicroNodes );

    if ( error ){
        error->print();
        results << "test_initializeFileInterfaces & False\n";
        return 1;
    }

    unsigned int rs = 0;
    for ( unsigned int i = 0; i < resultMicroNodes.size(); i++ ){
        std::cout << resultMicroNodes[ i ] << ", ";

        rs++;

        if ( rs > 3 * 5 ){
            
            std::cout << "\n";

            rs = 0;

        }
    }

    results << "test_initializeFileInterfaces & True\n";
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

    //Close the results file
    results.close();
}
