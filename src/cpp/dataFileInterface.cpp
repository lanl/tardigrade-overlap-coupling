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

    errorOut dataFileBase::readMesh( const unsigned int increment, floatVector &nodalPositions ){
        /*!
         * Read a mesh from the datafile.
         *
         * :param const unsigned int increment: The increment at which to get the mesh
         * :param floatVector &nodalPositions: The positions of the nodes in the mesh organized
         *     in row-major format [ node, coordinates ]
         */

        return new errorNode( "readMesh", "The readMesh function is not defined" );
    }


    /*=========================================================================
    |                               XDMFDataFile                              |
    =========================================================================*/

    XDMFDataFile::XDMFDataFile( ) : dataFileBase( ){
        /*!
         * The XDMFDataFile default constructor
         *
         * The format of the XDMF file is assumed to have one unstrutured grid collection
         * of type temporal.
         *
         * TODO: Add capability for other mesh types
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

    errorOut XDMFDataFile::readMesh( const unsigned int increment, floatVector &nodalPositions ){
        /*!
         * Read the mesh from the XDMF datafile
         *
         * Only unstructured grids are current considered
         *
         * :param const unsigned int increment: The increment at which to get the mesh
         * :param floatVector &nodalPositions: The positions of the nodes in the mesh organized
         *     in row-major format [ node, coordinates ]
         */

        unsigned int nGridCollections = _domain->getNumberGridCollections( );

        //Return an error if no grid collections are defined
        if ( nGridCollections == 0 ){
            return new errorNode( "readMesh", "No grid collections are defined in the XDMF file" );
        }

        //Print warning message if the number of grid collections is greater than 1
        if ( nGridCollections > 1 ){
                std::cerr << "WARNING: The number of grid collections is greater than 1. Only the first one will be used.\n";
        }

        //Extract the grid holder
        shared_ptr< XdmfGridCollection > _gridHolder = _domain->getGridCollection( 0 );

        //Get the data for the unstructured grids
        unsigned int nUnstructuredGrids = _gridHolder->getNumberUnstructuredGrids( );
        if ( nUnstructuredGrids == 0 ){
            return new errorNode( "readMesh", "There are no unstructured grids defined in the output file" );
        }

        if ( increment >= nUnstructuredGrids ){
            return new errorNode( "readMesh", "The requested increment is higher than the number of grids" );
        }

        //Get the unstructured grid
        shared_ptr< XdmfUnstructuredGrid > _grid = _gridHolder->getUnstructuredGrid( increment );

        //Get the geometry of the grid
        shared_ptr< XdmfGeometry > _geom = _grid->getGeometry( );
        _geom->read( ); //Required

        if ( _geom->getType() != XdmfGeometryType::XYZ()){
            return new errorNode( "readMesh", "The geometry type must be XYZ" );
        }

        //Set the nodal positions
        nodalPositions.resize( _geom->getSize( ) );
        _geom->getValues( 0, nodalPositions.data(), _geom->getSize( ), 1, 1 );

        return NULL;
    }
}
