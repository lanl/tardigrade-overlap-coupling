/*=============================================================================
|                              dataFileInterface                              |
===============================================================================
| A collection of tools that allow for data-files to be interacted wtih both  |
| for reading and writing of those files.                                     |
===============================================================================
| Currently supported file standards:                                         |
| XDMF                                                                        |
=============================================================================*/

#include<dataFileInterface.h>
#include<boost/algorithm/string.hpp>

namespace dataFileInterface{

    dataFileBase::dataFileBase( ){
        /*!
         * The base dataFile constructor
         */

        _error = NULL;

        return;
    }

    dataFileBase::dataFileBase( const YAML::Node &config ){
        /*!
         * The base dataFile constructor
         *
         * filename and mode must both be defined in the file. filename indicates
         * the data-base filename and mode indicates if it should be "read" or "write"
         * or some other user-defined mode.
         *
         * :param const YAML::Node &config: The YAML configuration file associated with the datafile
         */

        //Set the configuration file
        _config = config;

        if ( _config["filename"] ){
            _filename = _config["filename"].as<std::string>();
        }
        else{
            _error = new errorNode( "dataFileBase", "The filename must be defined in the YAML configuration file" );
            return;
        }

        if ( _config["mode"] ){
            _mode = _config["mode"].as<std::string>();
            boost::algorithm::to_lower( _mode );
        }
        else{
            _error = new errorNode( "dataFileBase", "The mode for the data file must be defined in the YAML configuration file" );
            return;
        }

        _error = NULL;

        return;
    }

    std::unique_ptr<dataFileBase> dataFileBase::Create( const std::string &type ){
        /*!
         * Create a new dataFile object
         *
         * :param const std::string &type: The name of the object to be created.
         */

        auto T = registryMap.find( type );

        if ( T == registryMap.end() ){
            return NULL; //The requested class is not defined
        }
        else{ //Register new dataFile objects below
            if ( T->second == XDMF ){
                return make_unique<XDMFDataFile>( _config );
            }
            else{
                return NULL;
            }
        }
        return NULL;
    }


    /*=========================================================================
    |                               XDMFDataFile                              |
    =========================================================================*/

    XDMFDataFile::XDMFDataFile( ) : dataFileBase( ){
        /*!
         * The XDMFDataFile default constructor
         */

        return;
    }

    XDMFDataFile::XDMFDataFile( const YAML::Node &configuration ) : dataFileBase( configuration ){
        /*!
         * The XDMFDataFile constructor with a YAML node
         *
         * :param const YAML::Node &configuration: The YAML configuration file for the XDMFDataFile object
         */

        //Initialize Readmode if required
        if ( _mode.compare( "read" ) );
        _initializeReadMode( );
        if ( _error ){
            return;
        }
        
        //Initialize Writemode if required
        //TODO: Implement this

    }

    void XDMFDataFile::_initializeReadMode( ){
        /*
         * Initialize the read mode for the XDMF file-type
         */

        //Initialize the reader
        _reader = XdmfReader::New( );

        //Read in the XDMF domain
        try{
            _domain = shared_dynamic_cast< XdmfDomain >( _reader->read( _filename ) );
        }
        catch( XdmfError &e ){
            _error = new errorNode( "XDMFDataFile", e.what() );
        }
    }
}
