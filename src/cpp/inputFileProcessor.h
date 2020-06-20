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

#include<dataFileInterface.h>

namespace inputFileProcessor{

    typedef errorTools::Node errorNode; //!Redefinition for the error node
    typedef errorNode* errorOut; //!Redefinition for a pointer to the error node
    typedef double floatType; //!Define the float values type.
    typedef std::vector< floatType > floatVector; //! Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!Define a matrix of floats
    typedef std::vector< unsigned int > uIntVector; //!Define a vector of unsigned ints

    class dataFileReaderBase;

    class inputFileProcessor {
        /*!
         * The class from the file processor which 
         * enables access to output files for use in 
         * coupling a micromorphic macro-scale to a 
         * subscale. The configuration file is in the
         * following format. Parentheses indicate the
         * defaults.
         *
         * macroscale_definition:
         *     filename: path_to_macroscale_file
         *     mode: ( "read" )
         *     filetype: ("XDMF")
         *
         * microscale_definition:
         *     filename: path_to_microscale_file
         *     mode: ( "read" )
         *     filetype: ( "XDMF" )
         *
         * Use of the function entails the following:
         * reader = inputFileProcessor( YAML_filename, macroDataReader, microDataReader );
         */


        public:

            //Constructors
            inputFileProcessor( );
            inputFileProcessor( const std::string &yamlConfigurationFilename );


            //Destructor
            ~inputFileProcessor( );

            //Functions
            errorOut setConfigurationFilename( const std::string &configurationFilename );
            errorOut openConfigurationFile( );
            errorOut openConfigurationFile( const std::string &configurationFilename );

            errorOut initializeFileReaders( );

        private:
            //Private functions

            //Private Attributes
            std::string _configFilename = "";
            YAML::Node _config;

            std::shared_ptr< dataFileInterface::dataFileBase > _macroscale;
            std::shared_ptr< dataFileInterface::dataFileBase > _microscale;

    };

}

#endif
