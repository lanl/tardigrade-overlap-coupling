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

    errorOut dataFileBase::getIncrementTime( const uIntType &increment, floatType &time ){
        /*!
         * Get the increment time from the datafile
         *
         * :param const uIntType &increment: The increment at which to extract the time
         * :param floatType &time: The time associated with the increment
         */

        ( void ) increment;
        ( void ) time;
        return new errorNode( "getIncrementTime", "The getIncrementTime function is not defined" );
            
    }

    errorOut dataFileBase::getNumIncrements( uIntType &numIncrements ){
        /*!
         * Get the number of increments from the datafile.
         *
         * :param uIntType &numIncrements: The number of increments
         */

        ( void ) numIncrements;

        return new errorNode( "getNumIncrements", "The getNumIncrements function is not defined" ); 
    }

    errorOut dataFileBase::getNumNodes( const uIntType increment, uIntType &numNodes ){
        /*!
         * Get the number of nodes from the datafile.
         *
         * :param const uIntType increment: The increment at which to get the nodes
         * :param uIntType &numIncrements: The number of increments
         */

        ( void ) increment;
        ( void ) numNodes;

        return new errorNode( "getNumIncrements", "The getNumNodes function is not defined" ); 
    }

    errorOut dataFileBase::readMesh( const uIntType increment, floatVector &nodalPositions ){
        /*!
         * Read a mesh from the datafile.
         *
         * :param const uIntType increment: The increment at which to get the mesh
         * :param floatVector &nodalPositions: The positions of the nodes in the mesh organized
         *     in row-major format [ node, coordinates ]
         */

        ( void ) increment;
        ( void ) nodalPositions;

        return new errorNode( "readMesh", "The readMesh function is not defined" );
    }

    errorOut dataFileBase::getDomainNodes( const uIntType increment, const std::string domainName, uIntVector &domainNodes ){
        /*!
         * Get the nodes in the provided domain.
         *
         * :param const uIntType increment: The increment at which to get the domain values
         * :param const std::string domainName: The name of the domain to be extracted.
         * :param uIntVector &domainNodes: The nodes located in the domain.
         */

        ( void ) increment;
        ( void ) domainName;
        ( void ) domainNodes;

        return new errorNode( "getDomainNodes", "The getDomainNodes function is not defined" );
    }

    errorOut dataFileBase::getNumDomainNodes( const uIntType increment, const std::string domainName,
                                              uIntType &numDomainNodes ){
        /*!
         * Get the number of nodes in the provided domain.
         *
         * :param const uIntType increment: The increment at which to get the domain's node count
         * :param const std::string domainName: The name of the domain to be interrogated.
         * :param uIntType numDomainNodes: The number of nodes in the domain
         */

        ( void ) increment;
        ( void ) domainName;
        ( void ) numDomainNodes;

        return new errorNode( "getDomainNodes", "The getNumDomainNodes function is not defined" );
    }

    errorOut dataFileBase::getSetNames( const uIntType increment, stringVector &setNames ){
        /*!
         * Get the set names from the provided domain at the current increment
         *
         * :param const uIntType increment: The increment at which to get the sets
         * :param stringVector &setNames: The names of the sets
         */

        ( void ) increment;
        ( void ) setNames;

        return new errorNode( "getSetNames", "The getSetNames function is not defined" );
    }

    errorOut dataFileBase::getSolutionData( const uIntType increment, const std::string &dataName, const std::string &dataCenter,
                                        floatVector &data ){
        /*!
         * Get the values of the mesh data named dataName. This should work for
         * nodal values as well as cell values.
         *
         * :param const uIntType increment: The increment at which to get the data
         * :param const std::string &dataName: The name of the data
         * :param const std::string &dataCenter: The type of the data. This will either be "Node" or "Cell"
         * :param floatVector &data: The output data vector
         *
         */

        ( void ) increment;
        ( void ) dataName;
        ( void ) dataCenter;
        ( void ) data;

        return new errorNode( "getSolutionData", "The getSolutionData function is not defined" );
    }

    errorOut dataFileBase::getSolutionVectorDataFromComponents( const uIntType increment,
                                                                const stringVector &componentNames,
                                                                const std::string &dataCenter, floatVector &data ){
        /*!
         * Extract solution vector ( and tensor ) data from a file where the components are saved individually
         *
         * :param const uIntType increment: The increment at which to get the data
         * :param const stringVector &componentNames: The name of the data's components
         * :param const std::string &dataCenter: The type of the data. This will either be "Node" or "Cell"
         * :param floatVector &data: The output data vector
         */

        unsigned int nComponents = componentNames.size( );

        //Check the component names vector
        if ( nComponents == 0 ){

            data.clear( );
            return NULL;

        }

        unsigned int nDataPoints = 0;
        errorOut error = NULL;

        //Loop over the component names
        floatVector componentData;

        for ( auto cN = componentNames.begin( ); cN != componentNames.end( ); cN++ ){

            error = getSolutionData( increment, *cN, dataCenter, componentData );

            if ( error ){

                errorOut result = new errorNode( "getSolutionVectorDataFromComponents",
                                                 "Error in the extraction of component " + *cN );
                result->addNext( error );
                return result;

            }

            if ( ( cN - componentNames.begin( ) ) == 0 ){

                data.clear( );
                nDataPoints = componentData.size( );
                data = floatVector( nComponents * nDataPoints );

            }

            if ( componentData.size( ) != nDataPoints ){

                return new errorNode( "getSolutionVectorDataFromComponents",
                                      "The component " + *cN + " does not have a consistent size with preceeding components" );

            }

            for ( auto c = componentData.begin( ); c != componentData.end( ); c++ ){

                data[ ( c - componentData.begin( ) ) * nComponents + cN - componentNames.begin( ) ] = *c;

            } 

        }
       
       return NULL;

    }

    errorOut dataFileBase::getMeshData( const uIntType increment,
                                        floatVector &nodePositions, uIntVector &connectivity, uIntVector &connectivityCellIndices,
                                        uIntType &cellCounts) {
        /*!
         * Get the mesh data from the datafile.
         *
         * :param const uIntType increment: The increment at which to get the data
         * :param floatVector &nodePositions: The node positions in the format [ x1, y1, z1, x2, y2, z2, ... ]
         * :param uIntVector &connectivity: The connectivity description in XDMF format [ element_type_1, ..., element_type_2, ..., ]
         * :param uIntVector &connectivityCellIndices: The indicices at which a new cell is defined in the connectivity vector.
         * :param uIntType &cellCounts: The number of cells present.
         */

        ( void ) increment;
        ( void ) nodePositions;
        ( void ) connectivity;
        ( void ) connectivityCellIndices;
        ( void ) cellCounts;

        return new errorNode( "getMeshData", "The getMeshData function is not defined" );
    }

    errorOut dataFileBase::connectivityToCellIndices( const uIntType &nCells, const uIntVector &connectivity,
                                                      uIntVector &connectivityCellIndices ){
        /*!
         * Compute the connectivity to cell indices
         *
         * :param const uIntType &nCells: The number of cells
         * :param const uIntVector &connectivity: The connectivity vector
         * :param uIntVector &connectivityCellIndices: The indices in the connectivity vector at which the cells are defined.
         */

        //Set the connectivity cell index vector size
        connectivityCellIndices = uIntVector( nCells, 0 );

        unsigned int index;
        uIntType index_connectivity = 0;
        uIntType cellDataCount;

        uIntType elementType;
        uIntType nFaces;

        for ( auto it = connectivityCellIndices.begin( ) + 1; it != connectivityCellIndices.end( ); it++ ){

            //Set the index of the connectivityCellIndices vector
            index = it - connectivityCellIndices.begin( );

            //Check if the connectivity vector is long enough
            if ( connectivity.size( ) <= index_connectivity ){

                return new errorNode( "connectivityToCellIndices",
                                      "The connectivity vector is to short for cell " + std::to_string( index ) );

            }

            //Get the element type
            elementType = connectivity[ index_connectivity ];

            //Perform different operations based on the element type
            if ( ( elementType == 2 ) && ( elementType == 3 ) ){ //Polyline and Polygon
                
                cellDataCount = connectivity[ index_connectivity + 1 ] + 1;

            }
            else if ( elementType == 16 ){ //Polyhedron

                if ( connectivity.size( ) <= index_connectivity + 1 ){

                    return new errorNode( "connectivityToCellIndices",
                                          "The connectivity vector is too short for the polyhedron definition of cell "
                                          + std::to_string( index ) );

                }

                cellDataCount = 1;
                nFaces = connectivity[ index_connectivity + cellDataCount ];

                for ( uIntType n = 0; n < nFaces; n++ ){

                    cellDataCount += 1; //Move to the new face

                    if ( connectivity.size( ) <= cellDataCount ){

                        return new errorNode( "connectivityToCellIndices",
                                              "The connectivity vector is too short for the polyhedron definition of cell "
                                              + std::to_string( index ) );

                    }

                    cellDataCount += connectivity[ index_connectivity + cellDataCount ]; //Add the number of nodes on the face

                }

            }
            else{

                auto elType = cellNodeCount.find( elementType );

                if ( elType == cellNodeCount.end( ) ){

                    return new errorNode( "connectivityToCellIndices",
                                          "The cell type " + std::to_string( elementType ) + " is not recognized" );

                }

                cellDataCount = elType->second;

            }

            //Update the cell indices
            index_connectivity += cellDataCount + 1;
            connectivityCellIndices[ index ] = index_connectivity;

        }

        return NULL;
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
        if ( _mode.compare( "read" ) == 0 ){
            _initializeReadMode( );
        }
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

    errorOut XDMFDataFile::getNumIncrements( uIntType &numIncrements ){
        /*!
         * Get the number of increments from the XDMF datafile
         *
         * Currently only considers unstructured grids
         *
         * :param uIntType &numIncrements: The number of increments
         */

        uIntType nGridCollections = _domain->getNumberGridCollections( );
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

    errorOut XDMFDataFile::getXDMFGridCollection( const uIntType gridCollectionNum, shared_ptr< XdmfGridCollection > &gridHolder ){
        /*!
         * Get a grid collection from the XDMF datafile
         *
         * :param const uIntType gridCollectionNum: The number of the grid collection to retrieve.
         * :param shared_ptr< XdmfGridCollection > &gridHolder: The returned grid collection.
         */

        uIntType nGridCollections = _domain->getNumberGridCollections( );

        //Return an error if no grid collections are defined
        if ( nGridCollections <= gridCollectionNum ){
            return new errorNode( "readMesh", "No grid collections are defined in the XDMF file" );
        }

        //Extract the grid holder
        gridHolder = _domain->getGridCollection( gridCollectionNum );

        return NULL;
    }

    errorOut XDMFDataFile::getUnstructuredGrid( const uIntType increment, shared_ptr< XdmfUnstructuredGrid > &unstructuredGrid ){
        /*!
         * Read an unstructured grid from the datafile
         *
         * :param const uIntType increment: The increment of the unstructured grid to extract
         * :param shared_ptr< XdmfUnstructuredGrid > &unstructuredGrid: The pointer to the unstructured grid
         */

        //Check the number of grid collections defined on the root level. It should be 1.
        uIntType numIncrements;
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
        uIntType nUnstructuredGrids = gridHolder->getNumberUnstructuredGrids( );
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

    errorOut XDMFDataFile::readMesh( const uIntType increment, floatVector &nodalPositions ){
        /*!
         * Read the mesh from the XDMF datafile
         *
         * Only unstructured grids are currently considered
         *
         * :param const uIntType increment: The increment at which to get the mesh
         * :param floatVector &nodalPositions: The positions of the nodes in the mesh organized
         *     in row-major format [ node, coordinates ]
         */

        //Get the grid
        shared_ptr< XdmfUnstructuredGrid > grid;
        errorOut error = getUnstructuredGrid( increment, grid );
        if ( error ){
            errorOut result = new errorNode( "readMesh", "Error in extraction of the grid" );
            result->addNext( error );
            return result;
        }

        //Get the geometry of the grid
        shared_ptr< XdmfGeometry > geom = grid->getGeometry( );
        geom->read( ); //Required to read from the HDF5 file

        if ( geom->getType() != XdmfGeometryType::XYZ()){
            return new errorNode( "readMesh", "The geometry type must be XYZ" );
        }

        //Get the nodal positions
        nodalPositions.resize( geom->getSize( ) );
        geom->getValues( 0, nodalPositions.data(), geom->getSize( ), 1, 1 );

//        //Get the set names
//        std::cout << "the number of sets: " << grid->getNumberSets( ) << "\n";
//        for ( uIntType s = 0; s < grid->getNumberSets( ); s++ ){
//            auto set = grid->getSet( s );
//            set->read( );
//            std::cout << "  name: " << set->getName( ) << "\n";
//            std::cout << "    type: ";
//            if ( set->getType() == XdmfSetType::Node( ) ){
//                std::cout << "Node\n";
//            }
//            if ( set->getType() == XdmfSetType::Cell( ) ){
//                std::cout << "Element\n";
//            }
//            std::cout << "    number of values: " << set->getSize( ) << "\n";
//
//        }

        return NULL;
    }

    errorOut XDMFDataFile::getNumNodes( const uIntType increment, uIntType &numNodes ){
        /*!
         * Get the number of nodes in the datafile
         *
         * :param const uIntType increment: The increment at which to get the nodes
         * :param uIntType &numNodes: The number of nodes
         */

        shared_ptr< XdmfUnstructuredGrid > grid;
        errorOut error = getUnstructuredGrid( increment, grid );

        if ( error ){
            errorOut result = new errorNode( "getNumNodes", "Error in extraction of grid" );
            result->addNext( error );
            return result;
        }

        if ( grid->getGeometry()->getType() != XdmfGeometryType::XYZ( ) ){
            return new errorNode( "getNumNodes", "The geometry type must be XYZ" );
        }

        numNodes = grid->getGeometry()->getSize() / 3;

        return NULL;
        
    }

    errorOut XDMFDataFile::getDomainNodes( const uIntType increment, const std::string domainName, uIntVector &domainNodes ){
        /*!
         * Get the nodes in the given domain
         *
         * :param const uIntType increment: The time increment at which to get the domain's nodes
         * :param const std::string domainName: The name of the domain
         * :param uIntVector &domainNodes: The nodes located in the domain
         */

        //Get the grid
        shared_ptr< XdmfUnstructuredGrid > grid;
        errorOut error = getUnstructuredGrid( increment, grid );
        if ( error ){
            errorOut result = new errorNode( "getDomainNodes", "Error in extraction of the grid" );
            result->addNext( error );
            return result;
        }

        //Get the set
        shared_ptr< XdmfSet > set = grid->getSet( domainName );
        if ( !set ){
            return new errorNode( "getDomainNodes", "No domain of name " + domainName + " found" );
        }

        //Required for HDF5 files
        set->read( );

        //If the set-type is Node, we can just get the values directly
        domainNodes.clear();
        if ( set->getType( ) == XdmfSetType::Node( ) ){

            domainNodes.resize( set->getSize( ) );
            set->getValues( 0, domainNodes.data(), set->getSize( ), 1, 1 );
            return NULL;

        }
        //TODO: Add more ability to extract the nodes from non-cell entities
        else{

            return new errorNode( "getDomainNodes", "The set type is not recognized. It must be Node" );

        }

        return NULL;
    }

    errorOut XDMFDataFile::getNumDomainNodes( const uIntType increment, const std::string domainName,
                                              uIntType &numDomainNodes ){
        /*!
         * Get the number of nodes in the provided domain.
         *
         * :param const uIntType increment: The increment at which to get the domain's node count
         * :param const std::string domainName: The name of the domain to be interrogated.
         * :param uIntType numDomainNodes: The number of nodes in the domain
         */

        //Get the grid
        shared_ptr< XdmfUnstructuredGrid > grid;
        errorOut error = getUnstructuredGrid( increment, grid );
        if ( error ){
            errorOut result = new errorNode( "getDomainNodes", "Error in extraction of the grid" );
            result->addNext( error );
            return result;
        }

        //Get the set
        shared_ptr< XdmfSet > set = grid->getSet( domainName );
        if ( !set ){
            return new errorNode( "getDomainNodes", "No domain of name " + domainName + " found" );
        }

        numDomainNodes = set->getSize( );

        return NULL;

    }

    errorOut XDMFDataFile::getSetNames( const uIntType increment, stringVector &setNames ){
        /*!
         * Get the set names from the provided domain at the current increment
         *
         * :param const uIntType increment: The increment at which to get the sets
         * :param stringVector &setNames: The names of the sets
         */

        //Get the grid
        shared_ptr< XdmfUnstructuredGrid > grid;
        errorOut error = getUnstructuredGrid( increment, grid );
        if ( error ){
            errorOut result = new errorNode( "getSetNames", "Error in extraction of the grid" );
            result->addNext( error );
            return result;
        }

        //Get the sets from the grid
        setNames = stringVector( grid->getNumberSets( ) );

        for ( uIntType i = 0; i < setNames.size(); i++ ){

            setNames[ i ] = grid->getSet( i )->getName( );

        }

        return NULL;
    }

    errorOut XDMFDataFile::getSolutionData( const uIntType increment, const std::string &dataName, const std::string &dataCenter,
                                            floatVector &data ){
        /*!
         * Get the values of the solution data named dataName. This should work for
         * nodal values as well as cell values.
         *
         * :param const uIntType increment: The increment at which to get the data
         * :param const std::string &dataName: The name of the data
         * :param const std::string &dataCenter: The type of the data. This will either be "Node" or "Cell"
         * :param floatVector &data: The output data vector
         *
         */

        //Get the grid
        shared_ptr< XdmfUnstructuredGrid > grid;
        errorOut error = getUnstructuredGrid( increment, grid );
        if ( error ){
            errorOut result = new errorNode( "getSolutionData", "Error in the extraction of the grid" );
            result->addNext( error );
            return result;
        }

        //Set the name and center
        std::string name;
        shared_ptr< const XdmfAttributeCenter > center;


        //Set the center
        std::string centerName = dataCenter;
        boost::algorithm::to_lower( centerName );

        if ( centerName.compare( "node" ) == 0 ){
            center = XdmfAttributeCenter::Node( );
        }
        else if ( centerName.compare( "cell" ) == 0 ){
            center = XdmfAttributeCenter::Cell( );
        }
        else{
            return new errorNode( "getSolutionData", "The dataCenter must either be 'Node' or 'Cell'" ); 
        }


        //Find the attribute name and type that matches up with the requested values
        shared_ptr< XdmfAttribute > attribute;
        for ( uIntType a = 0; a < grid->getNumberAttributes( ); a++ ){

            //Get the current attribute
            attribute = grid->getAttribute( a );

            //Check if the data is right
            if ( ( dataName.compare( attribute->getName( ) ) == 0 ) &&
                 ( attribute->getCenter( ) == center ) ){

                data = floatVector( attribute->getSize( ) );
                attribute->read( );
                grid->getAttribute( a )->getValues( 0, data.data( ), grid->getAttribute( a )->getSize( ), 1, 1 );

                return NULL;

            }

        }

        return new errorNode( "getSolutionData",
                              "Attribute with dataName '" + dataName + "' and center '" + dataCenter + "' was not found" );

    }

    errorOut XDMFDataFile::getMeshData( const uIntType increment,
                                        floatVector &nodePositions, uIntVector &connectivity, uIntVector &connectivityCellIndices,
                                        uIntType &cellCounts ){
        /*!
         * Get the mesh data from the datafile.
         *
         * :param const uIntType increment: The increment at which to get the data
         * :param floatVector &nodePositions: The node positions in the format [ x1, y1, z1, x2, y2, z2, ... ]
         * :param uIntVector &connectivity: The connectivity description in XDMF format [ element_type_1, ..., element_type_2, ..., ]
         *     Note the order that the elements appear in the vector are assumed to be their cell Ids starting at zero
         * :param uIntVector &connectivityCellIndices: The indicices at which a new cell is defined in the connectivity vector.
         * :param uIntType &cellCounts: The number of cells present.
         */

        //Get the grid
        shared_ptr< XdmfUnstructuredGrid > grid;
        errorOut error = getUnstructuredGrid( increment, grid );

        if ( error ){
            errorOut result = new errorNode( "getMeshData", "Error in the extraction of the mesh's grid" );
            result->addNext( error );
            return result;
        }

        //Get the geometry
        shared_ptr< XdmfGeometry > geom = grid->getGeometry( );

        //Read the geometry into memory
        geom->read( );

        //Set the node positions
        nodePositions = floatVector( geom->getSize( ) );
        geom->getValues( 0, nodePositions.data( ), geom->getSize( ), 1, 1 );

        //Get the topology
        shared_ptr< XdmfTopology > topology = grid->getTopology( );

        //Read the topology into memory
        topology->read( );

        //Set the topology values
        connectivity = uIntVector( topology->getSize( ) );
        topology->getValues( 0, connectivity.data( ), topology->getSize( ), 1, 1 );

        if ( !_config[ "cell_id_variable_name" ] ){
            return new errorNode( "getMeshData", "The key 'cell_id_variable_name' is not defined" );
        }
        else if ( !_config[ "cell_id_variable_name" ].IsScalar( ) ){
            return new errorNode( "getMeshData", "The key 'cell_id_variable_name' must be a scalar value" );
        }
        else if ( !grid->getAttribute( _config[ "cell_id_variable_name" ].as< std::string >( ) ) ){
            return new errorNode( "getMeshData", "The 'cell_id_variable_name' specified does not exist in the output file" );
        }
        else{
            shared_ptr< XdmfAttribute > attribute = grid->getAttribute( _config[ "cell_id_variable_name" ].as< std::string >( ) );
            cellCounts = attribute->getSize( );
        }

        error = connectivityToCellIndices( cellCounts, connectivity, connectivityCellIndices );

        return NULL;
    }

    errorOut XDMFDataFile::getIncrementTime( const uIntType &increment, floatType &time ){
        /*!
         * Get the increment time from the datafile
         *
         * :param const uIntType &increment: The increment at which to extract the time
         * :param floatType &time: The time associated with the increment
         */

        //Check the number of increments
        uIntType numIncrements;
        errorOut error = getNumIncrements( numIncrements );

        if ( error ){

            errorOut result = new errorNode( "getIncrementTime", "Error in the extraction of the increment count" );
            result->addNext( error );
            return result;

        }

        shared_ptr< XdmfUnstructuredGrid > grid;
        error = getUnstructuredGrid( increment, grid );

        if ( error ){

            errorOut result = new errorNode( "getIncrementTime", "Error in the extraction of the grid" );
            result->addNext( error );
            return result;

        }

        time = grid->getTime( )->getValue( );

        return NULL;
            
    }

}
