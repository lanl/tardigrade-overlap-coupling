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

        //Initialize the output writer
        if ( !_config[ "output_configuration" ] ){
            std::cerr << "'output_configuration' not specified in the YAML file. Using defaults\n";
            _config[ "output_configuration" ][ "filename" ] = "xdmf_out";

        }

        //Force overwrite of 'mode' to be 'write'        
        _config[ "output_configuration" ][ "mode" ] = "write";

        //Force the filetype to be XDMF
        _config[ "output_configuration" ][ "filetype" ] = "XDMF";

        _writer = dataFileInterface::dataFileBase( _config[ "output_configuration" ] ).create( );

        if ( _writer->_error ){

            _error.reset( new errorNode( "fileGenerator", "Error when forming dataFileInterface writer" ) );
            _error->addNext( _writer->_error );
            return;

        }

        return;

    }

    const std::unique_ptr< errorNode > &fileGenerator::getError( ){
        /*!
         * Get a reference to the error object
         */

        return _error;
    }

}
