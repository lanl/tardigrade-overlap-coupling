/*=============================================================================
|                                inputFileProcessor                           |
===============================================================================
| Process input files to put them in a format that can be read by the overlap |
| coupling toolchain.                                                         |
=============================================================================*/

#include<inputFileProcessor.h>

namespace inputFileProcessor{

    inputFileProcessor::inputFileProcessor( ){
        /*!
         * The input file processor constructor
         */

        return;
    }

    inputFileProcessor::inputFileProcessor( const std::string &configurationFilename ){
        /*!
         * The input file processor constructor
         *
         * :param const std::string &configurationFilename: The filename for the YAML configuration file.
         */

        //Set the configuration filename
        _configFilename = configurationFilename;

        return;
    }

    inputFileProcessor::~inputFileProcessor( ){
        /*!
         * The destructor
         */

        //Write the configuration filename out
        if ( !(_configFilename.compare( "" ) == 0 ) ){
            std::ofstream yamlOut( _configFilename + ".as_evaluated" );
            yamlOut << _config;

        }
    }

    errorOut inputFileProcessor::setConfigurationFilename( const std::string &configurationFilename ){
        /*!
         * Set the configuration filename
         *
         * :param std::string &configurationFilename: The configuration filename.
         */

        _configFilename = configurationFilename;
        return NULL;
    }

    errorOut inputFileProcessor::openConfigurationFile( ){
        /*!
         * Open the configuration file
         */

        if ( _configFilename.compare( "" ) == 0 ){
            return new errorNode( "inputFileProcessor",
                                  "The configuration filename has not been set" );
        }

        //Open the YAML configuration file
        try {

            _config = YAML::LoadFile( _configFilename );

        }
        catch ( YAML::BadFile ){
            return new errorNode( "inputFileProcessor",
                                 "Bad file passed to YAML file" );
        }
        catch ( ... ){
            return new errorNode( "inputFileProcessor",
                                  "Invalid YAML file" );
        }

        return NULL;
    }

    errorOut inputFileProcessor::openConfigurationFile( const std::string &configurationFilename ){
        /*!
         * Open the configuration file
         *
         * :param const std::string &configurationFilename: The configuration filename.
         */

        setConfigurationFilename( configurationFilename );
        return openConfigurationFile( );
    }

    errorOut inputFileProcessor::initializeFileReaders( ){
        /*!
         * Initialize the file readers
         */

        if ( _config[ "macroscale_definition" ] ){

            //Set the Default values
            if ( !_config[ "macroscale_definition" ]["mode"] ){

                _config[ "macroscale_definition" ][ "mode" ] = "read";

            }
            if ( !_config[ "macroscale_definition" ]["filetype"] ){

                _config[ "macroscale_definition" ][ "filetype" ] = "XDMF";

            }

            _macroscale = dataFileInterface::dataFileBase.Create( _config[ "macroscale_definition" ] );
//
//            _macroscale = std::make_shared<macroReader>( _config[ "macroscale_definition" ] );
//
        }
        else{

            return new errorNode( "inputFileProcessor",
                                  "There is no macroscale_definition in the YAML configuration file" );

        }
        if ( _config[ "microscale_definition" ] ){

            //Set the Default values
            if ( !_config[ "microscale_definition" ]["mode"] ){

                _config[ "microscale_definition" ][ "mode" ] = "read";

            }
            if ( !_config[ "microscale_definition" ]["filetype"] ){

                _config[ "microscale_definition" ][ "filetype" ] = "XDMF";

            }
//
//            _microscale = std::make_shared<microReader>( _config[ "microscale_definition" ] );
//
        }
        else{

            return new errorNode( "inputFileProcessor",
                                  "There is no microscale_definition in the YAML configuration file" );

        }

        return NULL;
    }

}
