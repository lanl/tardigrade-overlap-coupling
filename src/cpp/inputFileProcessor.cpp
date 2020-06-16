/*=============================================================================
|                                inputFileProcessor                           |
===============================================================================
| Process input files to put them in a format that can be read by the overlap |
| coupling toolchain.                                                         |
=============================================================================*/

#include<inputFileProcessor.h>

namespace inputFileProcessor{

    inputFileProcessorBase::inputFileProcessorBase( ){
        /*!
         * The base input file processor constructor
         */

        return;
    }

    inputFileProcessorBase::inputFileProcessorBase( const std::string &configurationFilename ){
        /*!
         * The base input file processor constructor
         *
         * :param const std::string &configurationFilename: The filename for the YAML configuration file.
         */

        //Set the configuration filename
        _configFilename = configurationFilename;

        return;
    }

    errorOut inputFileProcessorBase::setConfigurationFilename( const std::string &configurationFilename ){
        /*!
         * Set the configuration filename
         *
         * :param std::string &configurationFilename: The configuration filename.
         */

        _configFilename = configurationFilename;
        return NULL;
    }

    errorOut inputFileProcessorBase::openConfigurationFile( ){
        /*!
         * Open the configuration file
         */

        if ( _configFilename.compare( "" ) == 0 ){
            return new errorNode( "inputFileProcessorBase",
                                  "The configuration filename has not been set" );
        }

        //Open the YAML configuration file
        try {

            config = YAML::LoadFile( _configFilename );

        }
        catch ( YAML::BadFile ){
            return new errorNode( "inputFileProcessorBase",
                                 "Bad file passed to YAML file" );
        }
        catch ( ... ){
            return new errorNode( "inputFileProcessorBase",
                                  "Invalid YAML file" );
        }

        return NULL;
    }

    errorOut inputFileProcessorBase::openConfigurationFile( const std::string &configurationFilename ){
        /*!
         * Open the configuration file
         *
         * :param const std::string &configurationFilename: The configuration filename.
         */

        setConfigurationFilename( configurationFilename );
        return openConfigurationFile( );
    }

}
