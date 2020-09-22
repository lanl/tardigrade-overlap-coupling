/*=============================================================================
|                             generateXDMFData.cpp                            |
===============================================================================
| Functions used to generate XDMF data which can either be used to test the   |
| overlapCoupling class or to create filtering domains.                       |
=============================================================================*/

#include<generateXDMFData.h>

namespace fileGenerator{

    fileGenerator::fileGenerator( ){
        /*!
         * The default constructor
         */

        return;
    }

    fileGenerator::fileGenerator( const std::string &yamlFilename ){
        /*!
         * The constructor which will read in the YAML configuration file
         *
         * :param const std::string &yamlFilename: The YAML configuration file
         *
         */


        try{

            _config = YAML::LoadFile( yamlFilename );

        }
        catch( std::exception &e ){

            _error.reset( new errorNode( "fileGenerator", e.what( ) ) );

        }

        if ( !_config ){
            _error.reset( new errorNode( "fileGenerator", "File not found" ) );
        }

    }

    const std::unique_ptr< errorNode > &fileGenerator::getError( ){
        /*!
         * Get a reference to the error object
         */

        return _error;
    }

}
