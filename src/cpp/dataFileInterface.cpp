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

        //Make sure that the configuration file actually has something in it
        if( !config[ "filename" ] ){

            _error.reset( new errorNode( "dataFileBase", "The configuration file is empty" ) );
            return;

        }

        //Set the configuration file
        _config = YAML::Clone( config );

        if ( _config["filename"] ){
            _filename = _config["filename"].as<std::string>();
        }
        else{
            _error.reset( new errorNode( "dataFileBase", "The filename must be defined in the YAML configuration file" ) );
            return;
        }

        if ( _config["mode"] ){
            _mode = _config["mode"].as<std::string>();
            boost::algorithm::to_lower( _mode );
        }
        else{
            _error.reset( new errorNode( "dataFileBase", "The mode for the data file must be defined in the YAML configuration file" ) );
            return;
        }

        _error = NULL;

        return;
    }

    dataFileBase::dataFileBase( const YAML::Node &config, std::unique_ptr< errorNode > &error ) : dataFileBase( config ){
        /*!
         * The base dataFile constructor
         *
         * filename and mode must both be defined in the file. filename indicates
         * the data-base filename and mode indicates if it should be "read" or "write"
         * or some other user-defined mode.
         *
         * :param const YAML::Node &config: The YAML configuration file associated with the datafile
         * :param std::unique_ptr< errorNode > error: The current value of the error
         */

        //Copy the error
        _error = std::move( error );

        return;
    }

    std::shared_ptr<dataFileBase> dataFileBase::create( ){
        /*!
         * Create a new dataFile object using the configuration file
         */

        if ( _config[ "filetype" ] ){
            return dataFileBase::create( _config[ "filetype" ].as< std::string >( ) );
        }
        
        _error.reset( new errorNode( "create", "The filetype is not defined" ) );
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
            _error.reset( new errorNode( "create", "The filetype ( " + type + " ) is not recognized" ) );
            return std::make_shared< dataFileBase >( _config, _error ); //The requested class is not defined
        }
        else{ //Register new dataFile objects below
            if ( T->second == XDMF ){
                return std::make_shared<XDMFDataFile>( _config );
            }
            else{
                _error.reset( new errorNode( "create", "The filetype ( " + type + " ) is not recognized" ) );
                return std::make_shared< dataFileBase >( _config, _error ); //The requested class is not defined
            }
        }

        _error.reset( new errorNode( "create", "You should never get here..." ) );
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

    errorOut dataFileBase::getNodeIds( const uIntType increment, const std::string &nodeIdAttributeName, uIntVector &domainNodes ){
        /*!
         * Get the global node ids from the domain.
         *
         * TODO: Add specification of collection here ( and in other file reads )
         *
         * :param const uIntType increment: The increment at which to extract the data
         * :param const std::string &nodeIdAttributeName: The name of the node-id attribute
         * :param uIntVector &domainNodes: The nodes in the domain
         */

        ( void ) increment;
        ( void ) domainNodes;
        ( void ) nodeIdAttributeName;

        return new errorNode( "getNodeIds", "Not implemented" );
    }

    errorOut dataFileBase::getCellIds( const uIntType increment, const std::string &cellIdAttributeName, uIntVector &domainCells ){
        /*!
         * Get the global cell ids from the domain.
         *
         * TODO: Add specification of collection here ( and in other file reads )
         *
         * :param const uIntType increment: The increment at which to extract the data
         * :param const std::string &cellIdAttributeName: The name of the cell-id attribute
         * :param uIntVector &domainNodes: The nodes in the domain
         */

        ( void ) increment;
        ( void ) domainCells;
        ( void ) cellIdAttributeName;

        return new errorNode( "getCellIds", "Not implemented" );
    }

    errorOut dataFileBase::getSubDomainNodes( const uIntType increment, const std::string subDomainName, uIntVector &domainNodes ){
        /*!
         * Get the nodes in the provided sub-domain ( node set ).
         *
         * :param const uIntType increment: The increment at which to get the domain values
         * :param const std::string subDomainName: The name of the sub-domain to be extracted.
         * :param uIntVector &domainNodes: The nodes located in the domain.
         */

        ( void ) increment;
        ( void ) subDomainName;
        ( void ) domainNodes;

        return new errorNode( "getSubDomainNodes", "The getSubDomainNodes function is not defined" );
    }

    errorOut dataFileBase::getNumSubDomainNodes( const uIntType increment, const std::string subDomainName,
                                              uIntType &numSubDomainNodes ){
        /*!
         * Get the number of nodes in the provided sub-domain.
         *
         * :param const uIntType increment: The increment at which to get the domain's node count
         * :param const std::string subDomainName: The name of the sub-domain to be interrogated.
         * :param uIntType numSubDomainNodes: The number of nodes in the domain
         */

        ( void ) increment;
        ( void ) subDomainName;
        ( void ) numSubDomainNodes;

        return new errorNode( "getNumSubDomainNodes", "The getNumSubDomainNodes function is not defined" );
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

            uIntType componentIndex = cN - componentNames.begin( );

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
                uIntType dataIndex = ( c - componentData.begin( ) ) * nComponents;

                uIntType index = dataIndex + componentIndex;

                data[ index ] = *c;

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

    errorOut dataFileBase::initializeIncrement( const floatType time, const uIntType &reference_increment,
                                                const uIntType &collectionNumber, uIntType &increment ){
        /*!
         * Initialize a new increment
         *
         * :param const floatType time: The increment's time
         * :param const uIntType &reference_increment: The reference increment for the current increment
         * :param const uIntType &collectionNumber: The collection number that the increment exists in
         * :param uIntType &increment: The initialized increment's number
         */

        ( void ) time;
        ( void ) reference_increment;
        ( void ) collectionNumber;
        ( void ) increment;

        return new errorNode( "initializeIncrement", "Not implemented" );

    }

    errorOut dataFileBase::writeIncrementMeshData( const uIntType increment, const uIntType collectionNumber,
                                                   const uIntVector &nodeIds, const uIntMatrix &nodeSets,
                                                   const stringVector &nodeSetNames, const floatVector &nodePositions,
                                                   const uIntVector &elementIds, const uIntMatrix &elementSets,
                                                   const stringVector &elementSetNames, const uIntVector &connectivity ){
        /*!
         * Write the given increment's mesh data to the file.
         *
         * :param const uIntType increment: The increment at which to write the mesh data
         * :param const uIntType collectionNumber: The collection number within which to write the mesh data
         *     This is zero if only one "collection" of increments is available. If there are multiple 
         *     collections of data which contain meshes which can be evolving over time there will be more.
         * :param const uIntVector &nodeIds: The global node ids
         * :param const uIntMatrix &nodeSets: The global node ids in the nodesets ( node set, global node ids )
         * :param const stringVector &nodeSetNames: The names of the nodesets
         * :param const floatVector &nodePositions: The node positions in the format [ x1, y1, z1, x2, y2, z2, ... ]
         * :param const uIntVector &elementIds: The global Id numbers of the element
         * :param const uIntMatrix &elementSets: The global element ids in the element sets ( element set, global element ids )
         * :param const stringVector &elementSetNames: The names of the element sets
         * :param const uIntVector &connectivity: The connectivity description in XDMF format [ element_type_1, ..., element_type_2, ..., ]
         */

        ( void ) increment;
        ( void ) collectionNumber;
        ( void ) nodeIds;
        ( void ) nodeSets;
        ( void ) nodeSetNames;
        ( void ) nodePositions;
        ( void ) elementIds;
        ( void ) elementSets;
        ( void ) elementSetNames;
        ( void ) connectivity;

        return new errorNode( "writeIncrementMeshData", "Not implemented" );
    }

    errorOut dataFileBase::writeScalarSolutionData( const uIntType increment, const uIntType collectionNumber,
                                                    const std::string &dataName, const std::string &dataType,
                                                    const floatVector &data ){
        /*!
         * Write the given increment's solution data to the file
         *
         * :param const uIntType increment: The increment to write the solution data to
         * :param const uInttype collectionNumber: The temporal collection to write the data to
         * :param const std::string &dataName: The name of the data
         * :param const std::string &dataType: The type of data ( "Node" or "Cell" )
         * :param const floatVector &data: The data to write
         */

        ( void ) increment;
        ( void ) collectionNumber;
        ( void ) dataName;
        ( void ) dataType;
        ( void ) data;

        return new errorNode( "writeScalarSolutionData", "Not implemented" );
    }

    errorOut dataFileBase::addRootCollection( const std::string &collectionName, const std::string &collectionDescription,
                                              uIntType &collectionNumber ){
        /*!
         * Add another collection to the root
         *
         * :param const std::string &collectionName: The name of the collection
         * :param const std::string &collectionDescription: The description of the collection
         * :param uIntType &collectionNumber: The number of the new collection
         */

        ( void ) collectionName;
        ( void ) collectionDescription;
        ( void ) collectionNumber;

        return new errorNode( "addRootCollection", "Not implemented" );
    }

    errorOut dataFileBase::writeSolutionData( const uIntType increment, const uIntType collectionNumber,
                                              const stringVector &dataNames, const std::string &dataType,
                                              const floatVector &data ){
        /*!
         * Write solution data to the file where the solution data is a vector of values.
         *
         * :param const uIntType increment: The increment to write to the file
         * :param const uIntType collectionNumber: The collection to write to
         * :param const stringVector &dataNames: The names of the data i.e. [ "v_1", "v_2", "v_3" ]
         * :param const std::string &dataType: The type of data ( "Node" or "Cell" )
         * :param const floatVector &data: The data to write, stored as [ v_11, v_12, v_13, v_21, v_22, v_23, v_31, v_32, v_33, v_41, ... ]
         *     where the example is given to mimic the dataNames example
         */

        if ( ( data.size( ) % dataNames.size( ) ) != 0 ){

            return new errorNode( "writeSolutionData", "The data and dataNames vectors don't have consistent sizes" );

        }

        errorOut error;
        floatVector componentValues( data.size( ) / dataNames.size( ), 0 );

        for ( uIntType i = 0; i < dataNames.size( ); i++ ){

            //Extract the data values
            uIntType indx = 0;
            for ( uIntType j = i; j < data.size( ); j += dataNames.size( ), indx++ ){

                componentValues[ indx ] = data[ j ];

            }

            error = writeScalarSolutionData( increment, collectionNumber, dataNames[ i ], dataType, componentValues );

            if ( error ){

                std::string outstr = "Error when writing ";
                outstr += dataNames[ i ];
                outstr += " to the output file\n";
                errorOut result = new errorNode( "writeSolutionData", outstr );
                result->addNext( error );
                return result;

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
         *
         * The format of the XDMF file is assumed to have one unstrutured grid collection
         * of type temporal.
         *
         * TODO: Add capability for other mesh types
         */

        return;
    }

    XDMFDataFile::~XDMFDataFile( ){
        /*!
         * The XDMFDataFile destructor
         */

        if ( ( _mode.compare( "write" ) == 0 ) && ( _domain != nullptr ) ){

            //Form the writer
            shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( _filename + ".h5", false );
            heavyWriter->setReleaseData( true );
            shared_ptr< XdmfWriter > writer = XdmfWriter::New( _filename + ".xdmf", heavyWriter );
    
            //Write the data to file
            _domain->accept( writer );

        }

        return;

    }

    XDMFDataFile::XDMFDataFile( const YAML::Node &configuration ) : dataFileBase( configuration ){
        /*!
         * The XDMFDataFile constructor with a YAML node
         *
         * :param const YAML::Node &configuration: The YAML configuration file for the XDMFDataFile object
         */

        //Initialize read or write mode if required
        if ( _mode.compare( "read" ) == 0 ){
            _initializeReadMode( );
            return;
        }
        else if ( _mode.compare( "write" ) == 0 ){
            _initializeWriteMode( _config );
            return;
        }
        else{
            _error.reset( new errorNode( "XDMFDataFile", "The data file mode " + _mode + " is not recognized" ) );
            return;
        }

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
            _error.reset( new errorNode( "XDMFDataFile", e.what() ) );
        }
    }

    void XDMFDataFile::_initializeWriteMode( YAML::Node &configuration ){
        /*!
         * Initialize the write mode for the XDMF file-type
         *
         * :param YAML::Node &configuration: The configuration file
         */

        //Try to extract configuration information

        std::string domainName   = "DOMAIN";
        std::string domainString = "Default domain name for micromorphic overlap coupling / filter output";


        if ( configuration [ "domain_name" ] ){

            try{

                domainName = configuration[ "domain_name" ].as< std::string >( );

            }
            catch( std::exception &e ){

                std::string output = "Error in YAML file for the output domain name: ";
                output += e.what( );
                _error.reset( new errorNode( "XDMFDataFile", output ) );
                return;

            }

        }
        else{

            configuration [ "domain_name" ] = domainName;

        }

        if ( configuration [ "domain_string" ] ){

            try{

                domainName = configuration[ "domain_string" ].as< std::string >( );

            }
            catch( std::exception &e ){

                std::string output = "Error in YAML file for the output domain string: ";
                output += e.what( );
                _error.reset( new errorNode( "XDMFDataFile", output ) );
                return;

            }

        }
        else{

            configuration[ "domain_string" ] = domainString;

        }

        if ( configuration [ "append_to_existing_file" ] ){

            try{

                _append = configuration[ "append_to_existing_file" ].as< bool >( );

            }
            catch( std::exception &e ){

                std::string output = "Error in YAML file for the append to existing file flag: ";
                output += e.what( );
                _error.reset( new errorNode( "XDMFDataFile", output ) );
                return;

            }

        }

        //Initialize the writer
        shared_ptr< XdmfWriter > writer;
        try{

            shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( _filename + ".h5", true );
            heavyWriter->setReleaseData( true );
            writer = XdmfWriter::New( _filename + ".xdmf", heavyWriter );

        }
        catch( XdmfError &e ){

            std::string output = "Error in forming the XDMF writer: ";
            output += e.what( );
            _error.reset( new errorNode( "XDMFDataFile", output ) );
            return;

        }

        if ( _append ){

            try{
                //Try to read the domain in from the file ( if it exists )
                std::cout << "Checking for existing XDMF file ( may throw XDMF warning )\n";
                _reader = XdmfReader::New( );
                _domain = shared_dynamic_cast<XdmfDomain>( _reader->read( _filename + ".xdmf" ) );

                uIntType numIncrements;
                _error.reset( getNumIncrements( numIncrements ) );

                if ( _error ){
                    return;
                }

                _increment_reference_grids = uIntVector( numIncrements, 0 );
            }
            catch( XdmfError &e ){

                //Initialize a new XDMF domain
                _domain = XdmfDomain::New( );
                shared_ptr< XdmfInformation > domaininfo = XdmfInformation::New( domainName, domainString );
                _domain->insert( domaininfo );
        
                //Create the main temporal collection
                shared_ptr< XdmfGridCollection > _gridHolder = XdmfGridCollection::New( );
                _gridHolder->setType( XdmfGridCollectionType::Temporal( ) );
                shared_ptr< XdmfInformation > _holderInfo
                    = XdmfInformation::New( "Main_Temporal_Collection", "The main temporal ( or iteration ) collection" );
                _gridHolder->insert( _holderInfo );
                _domain->insert( _gridHolder );
                _domain->accept( writer );
            }

        }
        else{

            //Initialize a new XDMF domain
            _domain = XdmfDomain::New( );
            shared_ptr< XdmfInformation > domaininfo = XdmfInformation::New( domainName, domainString );
            _domain->insert( domaininfo );
        
            //Create the main temporal collection
            shared_ptr< XdmfGridCollection > _gridHolder = XdmfGridCollection::New( );
            _gridHolder->setType( XdmfGridCollectionType::Temporal( ) );
            shared_ptr< XdmfInformation > _holderInfo
                = XdmfInformation::New( "Main_Temporal_Collection", "The main temporal ( or iteration ) collection" );
            _gridHolder->insert( _holderInfo );
            _domain->insert( _gridHolder );
            _domain->accept( writer );

        }

        return;

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

    errorOut XDMFDataFile::getSubDomainNodes( const uIntType increment, const std::string subDomainName, uIntVector &domainNodes ){
        /*!
         * Get the nodes in the given sub-domain
         *
         * :param const uIntType increment: The time increment at which to get the domain's nodes
         * :param const std::string subDomainName: The name of the sub-domain
         * :param uIntVector &domainNodes: The nodes located in the domain
         */

        //Get the grid
        shared_ptr< XdmfUnstructuredGrid > grid;
        errorOut error = getUnstructuredGrid( increment, grid );
        if ( error ){
            errorOut result = new errorNode( "getSubDomainNodes", "Error in extraction of the grid" );
            result->addNext( error );
            return result;
        }

        //Get the set
        shared_ptr< XdmfSet > set = grid->getSet( subDomainName );
        if ( !set ){
            return new errorNode( "getSubDomainNodes", "No domain of name " + subDomainName + " found" );
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

            return new errorNode( "getSubDomainNodes", "The set type is not recognized. It must be Node" );

        }

        return NULL;
    }

    errorOut XDMFDataFile::getNumSubDomainNodes( const uIntType increment, const std::string subDomainName,
                                              uIntType &numSubDomainNodes ){
        /*!
         * Get the number of nodes in the provided domain.
         *
         * :param const uIntType increment: The increment at which to get the domain's node count
         * :param const std::string subDomainName: The name of the domain to be interrogated.
         * :param uIntType numSubDomainNodes: The number of nodes in the domain
         */

        //Get the grid
        shared_ptr< XdmfUnstructuredGrid > grid;
        errorOut error = getUnstructuredGrid( increment, grid );
        if ( error ){
            errorOut result = new errorNode( "getSubDomainNodes", "Error in extraction of the grid" );
            result->addNext( error );
            return result;
        }

        //Get the set
        shared_ptr< XdmfSet > set = grid->getSet( subDomainName );
        if ( !set ){
            return new errorNode( "getSubDomainNodes", "No domain of name " + subDomainName + " found" );
        }

        numSubDomainNodes = set->getSize( );

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
         * :param uIntVector &connectivityCellIndices: The indices at which a new cell is defined in the connectivity vector.
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
        if ( topology->getType( ) == XdmfTopologyType::Mixed( ) ){

            connectivity = uIntVector( topology->getSize( ) );
            topology->getValues( 0, connectivity.data( ), topology->getSize( ), 1, 1 );

        }
        else{

            // Initialize the connectivity vector size
            cellCounts   = topology->getSize( ) / topology->getType( )->getNodesPerElement( );
            connectivity = uIntVector( topology->getSize( ) + cellCounts, 0 );
            uIntVector elementPoints( topology->getType( )->getNodesPerElement( ) );

            for ( unsigned int i = 0; i < cellCounts; i++ ){

                // Set the element type
                connectivity[ ( topology->getType( )->getNodesPerElement( ) + 1 ) * i ] = topology->getType( )->getID( );

                // Extract the element nodes
                topology->getValues( topology->getType( )->getNodesPerElement( ) * i, elementPoints.data( ), topology->getType( )->getNodesPerElement( ), 1, 1 );

                // Set the element nodes to the connectivity vector
                std::copy_n( elementPoints.begin( ), elementPoints.size( ), connectivity.begin( ) + ( topology->getType( )->getNodesPerElement( ) + 1 ) * i + 1 );

            }

        }

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

        if ( error ){

            errorOut result = new errorNode( "getMeshData", "Error when getting the cell indices for the connectivity vector" );
            result->addNext( error );
            return result;

        }

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

    errorOut XDMFDataFile::initializeIncrement( const floatType time, const uIntType &reference_increment, const uIntType &collectionNumber, uIntType &increment ){
        /*!
         * Initialize a new increment
         *
         * :param const floatType time: The increment's time
         * :param const uIntType &reference_increment: The reference increment for the current increment
         * :param const uInttype &collectionNumber: The collection which contains the new increment
         * :param uIntType &increment: The initialized increment's number
         */

        //Initialize the grid
        shared_ptr< XdmfUnstructuredGrid > grid = XdmfUnstructuredGrid::New( );

        //Write the timestamp information
        shared_ptr< XdmfTime > untime = XdmfTime::New( time );
        shared_ptr< XdmfInformation > untimeinfo = XdmfInformation::New( "Time", "This is the current value of the timestep" );
        untime->insert( untimeinfo );
        grid->setTime( untime );

        //Store the grid in the collection
        if ( _domain->getNumberGridCollections( ) <= collectionNumber ){
            return new errorNode( "initializeIncrement",
                                  "The collection number " + std::to_string( collectionNumber ) + " is out of range" );
        }
        shared_ptr< XdmfGridCollection > collection = _domain->getGridCollection( collectionNumber ); //TODO: This should be verified to be general enough
        collection->insert( grid );

//        //Form the writer
//        shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( _filename + ".h5", false );
//        heavyWriter->setReleaseData( true );
//        shared_ptr< XdmfWriter > writer = XdmfWriter::New( _filename + ".xdmf", heavyWriter );
//
//        //Write the data to file
//        _domain->accept( writer );

        //Set the reference increment
        _increment_reference_grids.push_back( reference_increment );

        increment  = _increment_reference_grids.size( ) - 1;

        return NULL;
    }

    errorOut XDMFDataFile::writeIncrementMeshData( const uIntType increment, const uIntType collectionNumber,
                                                   const uIntVector &nodeIds, const uIntMatrix &nodeSets,
                                                   const stringVector &nodeSetNames, const floatVector &nodePositions,
                                                   const uIntVector &elementIds, const uIntMatrix &elementSets,
                                                   const stringVector &elementSetNames, const uIntVector &connectivity ){
        /*!
         * Write the given increment's mesh data to the file. ( only works for 3D data )
         *
         * :param const uIntType increment: The increment at which to write the mesh data
         * :param const uIntType collectionNumber: The collection number to store the increment in
         * :param const uIntVector &nodeIds: The global node ids
         * :param const uIntMatrix &nodeSets: The global node ids of the nodesets ( n sets, n nodes in set )
         * :param const stringVector &nodeSetNames: The names of the nodesets
         * :param const floatVector &nodePositions: The node positions in the format [ x1, y1, z1, x2, y2, z2, ... ]
         * :param const uIntVector &elementIds: The global element id numbers
         * :param const uIntMatrix &elementSets: The global element ids in the element sets ( element set, global element ids )
         * :param const stringVector &elementSetNames: The names of the element sets
         * :param const uIntVector &connectivity: The connectivity description in XDMF format [ element_type_1, ..., element_type_2, ..., ]
         */

        if ( _domain->getGridCollection( collectionNumber )->getNumberUnstructuredGrids( ) < increment ){

            return new errorNode( "writeIncrementMeshData",
                                  "The increment to write increment to ( " + std::to_string( increment ) +
                                  " ) is not defined in the XDMF data file" );

        }

        if ( _increment_reference_grids.size( ) <= increment ){

            return new errorNode( "writeIncrementMeshData",
                                  "The increment to write increment to ( " + std::to_string( increment ) +
                                  " ) is larger than the increment reference grids vector can allow" );

        }

        const uIntType _reference_increment = _increment_reference_grids[ increment ];

        //Construct the output grid
        shared_ptr< XdmfUnstructuredGrid > grid = _domain->getGridCollection( collectionNumber )->getUnstructuredGrid( increment );

        //If the reference increment is not the current increment, we don't need to write out the mesh
        //but can point to the old mesh
        
        if ( _reference_increment != increment ){

            //Point the new grid to the current reference gitd's geometry, topology, and ids
            shared_ptr< XdmfUnstructuredGrid > reference_grid
                = _domain->getGridCollection( 0 )->getUnstructuredGrid( _reference_increment );

            //Reference the geometry and topology to the previous grid
            grid->setGeometry( reference_grid->getGeometry( ) );
            grid->setTopology( reference_grid->getTopology( ) );

            //Get the Node and Element ID references
            grid->insert( reference_grid->getAttribute( "NODEID" ) );
            grid->insert( reference_grid->getAttribute( "ELEMID" ) );

            //Get the nodesets
            for ( unsigned int i = 0; i < reference_grid->getNumberSets( ); i++ ){
                grid->insert( reference_grid->getSet( i ) );
            }

//            //Build the writer
//            shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( _filename + ".h5", false );
//            heavyWriter->setReleaseData( true );
//            shared_ptr< XdmfWriter > writer = XdmfWriter::New( _filename + ".xdmf", heavyWriter );
//
//            //Write the grid to the file
//            shared_ptr< XdmfGridCollection > collection = _domain->getGridCollection( collectionNumber ); //TODO: This should be verified to be general enough
//            _domain->accept( writer );

            return NULL;

        }

        if ( nodeSets.size( ) != nodeSetNames.size( ) ){

            return new errorNode( "writeIncrementMeshData",
                                  "The size of the node sets vector and the node set names vector are not the same size" );

        }

        //Save the node id information
        shared_ptr< XdmfAttribute > _nodeIds = XdmfAttribute::New( );
        _nodeIds->setType( XdmfAttributeType::GlobalId( ) );
        _nodeIds->setCenter( XdmfAttributeCenter::Node( ) );
        _nodeIds->setName( "NODEID" );
        _nodeIds->insert( 0, nodeIds.data( ), nodeIds.size( ), 1, 1 );
        shared_ptr< XdmfInformation > nodeIdsInfo = XdmfInformation::New( "ID", "The nodal IDs" );
        _nodeIds->insert( nodeIdsInfo );
        grid->insert( _nodeIds );

        //Save the nodesets
        std::vector< shared_ptr< XdmfSet > > XdmfNodeSets( nodeSets.size( ) );

        for ( auto it = nodeSets.begin( ); it != nodeSets.end( ); it++ ){

            unsigned int index = it - nodeSets.begin( );

            XdmfNodeSets[ index ] = XdmfSet::New( );
            XdmfNodeSets[ index ]->setType( XdmfSetType::Node( ) );
            XdmfNodeSets[ index ]->setName( nodeSetNames[ index ] );
            XdmfNodeSets[ index ]->insert( 0, nodeSets[ index ].data( ), nodeSets[ index ].size( ), 1, 1 );
            grid->insert( XdmfNodeSets[ index ] );

        }

        //Save the nodal coordinates
        shared_ptr< XdmfGeometry > nodeGeometry = XdmfGeometry::New( );
        nodeGeometry->setType( XdmfGeometryType::XYZ( ) );
        nodeGeometry->setName( "Coordinates" );
        nodeGeometry->insert( 0, nodePositions.data( ), nodePositions.size( ), 1, 1 );
        shared_ptr< XdmfInformation > nodeGeometryInfo = XdmfInformation::New( "Coordinates", "Coordinates of the nodes in x1, y1, z1, x2, ... format " );
        nodeGeometry->insert( nodeGeometryInfo );
        grid->setGeometry( nodeGeometry );

        //Set the topology
        shared_ptr< XdmfTopology > topology = XdmfTopology::New( );
        topology->setType( XdmfTopologyType::Mixed( ) );
        topology->setName( "Topology" );
        topology->insert( 0, connectivity.data( ), connectivity.size( ), 1, 1 );
        grid->setTopology( topology );

        //Set the element ID numbers
        shared_ptr< XdmfAttribute > _elementIds = XdmfAttribute::New( );
        _elementIds->setType( XdmfAttributeType::GlobalId( ) );
        _elementIds->setCenter( XdmfAttributeCenter::Cell( ) );
        _elementIds->setName( "ELEMID" );
        _elementIds->insert( 0, elementIds.data( ), elementIds.size( ), 1, 1 );
        shared_ptr< XdmfInformation > elementIdsInfo = XdmfInformation::New( "ID", "The element IDs" );
        _elementIds->insert( elementIdsInfo );
        grid->insert( _elementIds );

        //Save the element sets
        std::vector< shared_ptr< XdmfSet > > XdmfElementSets( elementSets.size( ) );

        for ( auto it = elementSets.begin( ); it != elementSets.end( ); it++ ){

            unsigned int index = it - elementSets.begin( );

            XdmfElementSets[ index ] = XdmfSet::New( );
            XdmfElementSets[ index ]->setType( XdmfSetType::Cell( ) );
            XdmfElementSets[ index ]->setName( elementSetNames[ index ] );
            XdmfElementSets[ index ]->insert( 0, elementSets[ index ].data( ), elementSets[ index ].size( ), 1, 1 );
            grid->insert( XdmfElementSets[ index ] );

        }

//        //Build the writer
//        shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( _filename + ".h5", false );
//        heavyWriter->setReleaseData( true );
//        shared_ptr< XdmfWriter > writer = XdmfWriter::New( _filename + ".xdmf", heavyWriter );
//
//        //Write to the output file
//        _domain->accept( writer );

        return NULL;

    }

    errorOut XDMFDataFile::addRootCollection( const std::string &collectionName, const std::string &collectionDescription,
                                              uIntType &collectionNumber ){
        /*!
         * Add another collection to the root. This is currently assumed to be a temporal collection
         *
         * :param const std::string &collectionName: The name of the collection
         * :param const std::string &collectionDescription: The description of the collection
         * :param uIntType &collectionNumber: The number of the new collection
         */

        //Create the main temporal collection
        shared_ptr< XdmfGridCollection > _gridHolder = XdmfGridCollection::New( );
        _gridHolder->setType( XdmfGridCollectionType::Temporal( ) );
        shared_ptr< XdmfInformation > _holderInfo
            = XdmfInformation::New( collectionName, collectionDescription );
        _gridHolder->insert( _holderInfo );
        _domain->insert( _gridHolder );

//        //Build the writer
//        shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( _filename + ".h5", false );
//        heavyWriter->setReleaseData( true );
//        shared_ptr< XdmfWriter > writer = XdmfWriter::New( _filename + ".xdmf", heavyWriter );
//
//        _domain->accept( writer );

        collectionNumber = _domain->getNumberGridCollections( ) - 1;

        return NULL;
    }

    errorOut XDMFDataFile::getNodeIds( const uIntType increment, const std::string &nodeIdAttributeName, uIntVector &domainNodes ){
        /*!
         * Get the global node ids from the domain.
         *
         * TODO: Add specification of collection here ( and in other file reads )
         *
         * :param const uIntType increment: The increment at which to extract the data
         * :param const std::string &nodeIdAttributeName: The name of the node id attribute
         * :param uIntVector &domainNodes: The nodes in the domain
         */

        shared_ptr< XdmfUnstructuredGrid > grid;
        errorOut error = getUnstructuredGrid( increment, grid );

        if ( error ){

            errorOut result = new errorNode( "getNodeIds", "Error in getting the unstructured grid" );
            result->addNext( error );
            return result;

        }

        try{

            shared_ptr< XdmfAttribute > _attribute = grid->getAttribute( nodeIdAttributeName );
            _attribute->read( );
            domainNodes.resize( _attribute->getSize( ) );
            _attribute->getValues( 0, domainNodes.data(), _attribute->getSize( ), 1, 1 );

        }
        catch( XdmfError &e ){

            std::string outstr = "Error in extraction of the node ids: ";
            outstr += e.what( );

            return new errorNode( "getNodeIds", outstr );
        }

        return NULL;

    }

    errorOut XDMFDataFile::getCellIds( const uIntType increment, const std::string &cellIdAttributeName, uIntVector &domainCells ){
        /*!
         * Get the global cell ids from the domain.
         *
         * TODO: Add specification of collection here ( and in other file reads )
         *
         * :param const uIntType increment: The increment at which to extract the data
         * :param const std::string &cellIdAttributeName: The name of the cell-id attribute
         * :param uIntVector &domainNodes: The nodes in the domain
         */

        shared_ptr< XdmfUnstructuredGrid > grid;
        errorOut error = getUnstructuredGrid( increment, grid );

        if ( error ){

            errorOut result = new errorNode( "getCellIds", "Error in getting the unstructured grid" );
            result->addNext( error );
            return result;

        }

        try{

            shared_ptr< XdmfAttribute > _attribute = grid->getAttribute( cellIdAttributeName );
            _attribute->read( );
            domainCells.resize( _attribute->getSize( ) );
            _attribute->getValues( 0, domainCells.data(), _attribute->getSize( ), 1, 1 );

        }
        catch( XdmfError &e ){

            std::string outstr = "Error in extraction of the cell ids: ";
            outstr += e.what( );

            return new errorNode( "getCellIds", outstr );
        }

        return NULL;
    }

    errorOut XDMFDataFile::writeScalarSolutionData( const uIntType increment, const uIntType collectionNumber,
                                                    const std::string &dataName, const std::string &dataType,
                                                    const floatVector &data ){
        /*!
         * Write the given increment's solution data to the file
         *
         * :param const uIntType increment: The increment to write the solution data to
         * :param const uInttype collectionNumber: The temporal collection to write the data to
         * :param const std::string &dataName: The name of the data
         * :param const std::string &dataType: The type of data ( "Node" or "Cell" )
         * :param const floatVector &data: The data to write
         */

        if ( _domain->getNumberGridCollections( ) <= collectionNumber ){

            return new errorNode( "writeScalarSolutionData",
                                  "The collection number " + std::to_string( collectionNumber ) +
                                  " is larger than the number of collections" );

        }

        if ( _domain->getGridCollection( collectionNumber )->getNumberUnstructuredGrids( ) <= increment ){

            return new errorNode( "writeScalarSolutionData",
                                  "The increment number " + std::to_string( increment ) +
                                  " is larger than the number of increments defined in collection " +
                                  std::to_string( collectionNumber ) );

        }

        shared_ptr< XdmfUnstructuredGrid > grid
            = _domain->getGridCollection( collectionNumber )->getUnstructuredGrid( increment );

        shared_ptr< XdmfAttribute > solutionVector = XdmfAttribute::New( );
        solutionVector->setType( XdmfAttributeType::Scalar( ) );
        solutionVector->setName( dataName );

        std::string _dataType = dataType;
        boost::algorithm::to_lower( _dataType );

        if ( _dataType.compare( "node" ) == 0 ){

            solutionVector->setCenter( XdmfAttributeCenter::Node( ) );

        }
        else if ( _dataType.compare( "cell" ) == 0 ){

            solutionVector->setCenter( XdmfAttributeCenter::Cell( ) );

        }
        else{

            std::string outstr = "data type ";
            outstr += dataType;
            outstr += " is not recognized";
            return new errorNode( "writeScalarSolutionData", outstr );

        }

        solutionVector->insert( 0, data.data( ), data.size( ), 1, 1 );

        shared_ptr< XdmfInformation > solutionVectorInformation = XdmfInformation::New( dataName, "Quantity " + dataName );

        solutionVector->insert( solutionVectorInformation );

        grid->insert( solutionVector );

//        //Build the writer
//        shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( _filename + ".h5", false );
//        heavyWriter->setReleaseData( true );
//        shared_ptr< XdmfWriter > writer = XdmfWriter::New( _filename + ".xdmf", heavyWriter );
//
//        _domain->accept( writer );

        return NULL;

    }

}
