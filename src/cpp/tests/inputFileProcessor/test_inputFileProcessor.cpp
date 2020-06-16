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
    inputFileProcessor::inputFileProcessorBase reader( filename );

    errorOut error = reader.openConfigurationFile( );

    if ( error ){
        error->print();
        results << "test_openConfigurationFile (test 1) & False\n";
        return 1;
    }

    reader = inputFileProcessor::inputFileProcessorBase( );
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
    inputFileProcessor::inputFileProcessorBase reader;

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

    //Close the results file
    results.close();
}
