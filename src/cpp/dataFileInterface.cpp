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

    dataFileBase::dataFileBase( const YAML::Node &config, errorOut error ) : dataFileBase( config ){
        /*!
         * The base dataFile constructor
         *
         * filename and mode must both be defined in the file. filename indicates
         * the data-base filename and mode indicates if it should be "read" or "write"
         * or some other user-defined mode.
         *
         * :param const YAML::Node &config: The YAML configuration file associated with the datafile
         * :param errorOut error: The current value of the error
         */

        //Copy the error
        _error = error;

        return;
    }

    std::shared_ptr<dataFileBase> dataFileBase::create( ){
        /*!
         * Create a new dataFile object using the configuration file
         */

        if ( _config[ "filetype" ] ){
            return dataFileBase::create( _config[ "filetype" ].as< std::string >( ) );
        }
        
        _error = new errorNode( "create", "The filetype is not defined" );
        return std::make_shared< dataFileBase >( _config, _error );
    }
    std::shared_ptr<dataFileBase> dataFileBase::create( const std::string &type ){
        /*!
         * Create a new dataFile object
         *
         * :param const std::string &type: The name of the object to be created.
         */

        auto T = registryMap.find( type );

        if ( T == registryMap.end() ){
            _error = new errorNode( "create", "The filetype ( " + type + " ) is not recognized" );
            return std::make_shared< dataFileBase >( _config, _error ); //The requested class is not defined
        }
        else{ //Register new dataFile objects below
            if ( T->second == XDMF ){
                return std::make_shared<XDMFDataFile>( _config );
            }
            else{
                _error = new errorNode( "create", "The filetype ( " + type + " ) is not recognized" );
                return std::make_shared< dataFileBase >( _config, _error ); //The requested class is not defined
            }
        }

        _error = new errorNode( "create", "You should never get here..." );
        return std::make_shared< dataFileBase >( _config, _error ); //The requested class is not defined
    }

    //Virtual functions to be overloaded
    
    errorOut dataFileBase::getNumIncrements( unsigned int &numIncrements ){
        /*!
         * Get the number of increments from the datafile.
         *
         * :param unsigned int &numIncrements: The number of increments
         */

        return new errorNode( "getNumIncrements", "The getNumIncrements function is not defined" ); 
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

    errorOut XDMFDataFile::getNumIncrements( unsigned int &numIncrements ){
        /*!
         * Get the number of increments from the XDMF datafile
         *
         * Currently only considers unstructured grids
         *
         * :param unsigned int &numIncrements: The number of increments
         */

        unsigned int nGridCollections = _domain->getNumberGridCollections( );
        //Print warning message if the number of root-level grid collections is greater than 1
        if ( nGridCollections > 1 ){
                std::cerr << "WARNING: The number of root-level grid collections is greater than 1. Only the first one will be used.\n";
        }

        shared_ptr< XdmfGridCollection > _gridHolder;
        errorOut error = getXDMFGridCollection( 0, _gridHolder );

        if ( error ){ 
            errorOut result = new errorNode( "getNumIncrements",
                                             "Error in getting the XDMF grid collection" );
            result->addNext( error );
            return result;
        }

        numIncrements = _gridHolder->getNumberUnstructuredGrids( );
        return NULL;
    }

    errorOut XDMFDataFile::getXDMFGridCollection( const unsigned int gridCollectionNum, shared_ptr< XdmfGridCollection > &gridHolder ){
        /*!
         * Get a grid collection from the XDMF datafile
         *
         * :param const unsigned int gridCollectionNum: The number of the grid collection to retrieve.
         * :param shared_ptr< XdmfGridCollection > &gridHolder: The returned grid collection.
         */

        unsigned int nGridCollections = _domain->getNumberGridCollections( );

        //Return an error if no grid collections are defined
        if ( nGridCollections <= gridCollectionNum ){
            return new errorNode( "readMesh", "No grid collections are defined in the XDMF file" );
        }

        //Extract the grid holder
        gridHolder = _domain->getGridCollection( gridCollectionNum );

        return NULL;
    }

    errorOut XDMFDataFile::getUnstructuredGrid( const unsigned int increment, shared_ptr< XdmfUnstructuredGrid > &unstructuredGrid ){
        /*!
         * Read an unstructured grid from the datafile
         *
         * :param const unsigned int increment: The increment of the unstructured grid to extract
         * :param shared_ptr< XdmfUnstructuredGrid > &unstructuredGrid: The pointer to the unstructured grid
         */

        //Check the number of grid collections defined on the root level. It should be 1.
        unsigned int numIncrements;
        errorOut error = XDMFDataFile::getNumIncrements( numIncrements );

        if ( error ){
            errorOut result = new errorNode( "readMesh", "Error in getting the number of increments" );
            result->addNext( error );
            return result;
        }

        //Get the grid holder
        shared_ptr< XdmfGridCollection > gridHolder;
        error = getXDMFGridCollection( 0, gridHolder );

        if ( error ){
            errorOut result = new errorNode( "readMesh", "Error in getting the grid collection" );
            result->addNext( error );
            return result;
        }

        //Get the data for the unstructured grids
        unsigned int nUnstructuredGrids = gridHolder->getNumberUnstructuredGrids( );
        if ( nUnstructuredGrids == 0 ){
            return new errorNode( "readMesh", "There are no unstructured grids defined in the output file" );
        }

        if ( increment >= nUnstructuredGrids ){
            return new errorNode( "readMesh", "The requested increment is higher than the number of grids" );
        }

        //Get the unstructured grid
        unstructuredGrid = gridHolder->getUnstructuredGrid( increment );

        return NULL;
    }

    errorOut XDMFDataFile::readMesh( const unsigned int increment, floatVector &nodalPositions ){
        /*!
         * Read the mesh from the XDMF datafile
         *
         * Only unstructured grids are currently considered
         *
         * :param const unsigned int increment: The increment at which to get the mesh
         * :param floatVector &nodalPositions: The positions of the nodes in the mesh organized
         *     in row-major format [ node, coordinates ]
         */

        shared_ptr< XdmfUnstructuredGrid > grid;
        getUnstructuredGrid( increment, grid );

        //Get the geometry of the grid
        shared_ptr< XdmfGeometry > geom = grid->getGeometry( );
        geom->read( ); //Required to read from the HDF5 file

        if ( geom->getType() != XdmfGeometryType::XYZ()){
            return new errorNode( "readMesh", "The geometry type must be XYZ" );
        }

        //Set the nodal positions
        nodalPositions.resize( geom->getSize( ) );
        geom->getValues( 0, nodalPositions.data(), geom->getSize( ), 1, 1 );

        return NULL;
    }
}