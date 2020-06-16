/*=============================================================================
|                                inputFileProcessor                           |
===============================================================================
| Process input files to put them in a format that can be read by the overlap |
| coupling toolchain.                                                         |
=============================================================================*/

#ifndef INPUTFILEPROCESSOR_H
#define INPUTFILEPROCESSOR

#define USE_EIGEN
#include<vector_tools.h>
#include<error_tools.h>
#include<yaml-cpp/yaml.h>

namespace inputFileProcessor{

    typedef errorTools::Node errorNode; //!Redefinition for the error node
    typedef errorNode* errorOut; //!Redefinition for a pointer to the error node
    typedef double floatType; //!Define the float values type.
    typedef std::vector< floatType > floatVector; //! Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!Define a matrix of floats
    typedef std::vector< unsigned int > uIntVector; //!Define a vector of unsigned ints

    class inputFileProcessorBase {
        /*!
         * The base class for the input file processor reader
         */


        public:

            //Constructors
            inputFileProcessorBase( );
            inputFileProcessorBase( const std::string &yamlConfigurationFilename );

            //Functions
            errorOut setConfigurationFilename( const std::string &configurationFilename );
            errorOut openConfigurationFile( );
            errorOut openConfigurationFile( const std::string &configurationFilename );
            
        private:
            std::string _configFilename = "";
            YAML::Node config;

    };
}

#endif
