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

        initialize( );

        return;
    }

    void inputFileProcessor::initialize( ){
        /*!
         * Initialize the input file processor
         */

        //Clear the error
        _error = NULL;

        //Open the configuration file
        errorOut error = openConfigurationFile( );
        if ( error ){
            _error = new errorNode( "initialize", "Error in opening the configuration file" );
            _error->addNext( error );
            return;
        }

        //Initialize the file interfaces
        error = initializeFileInterfaces( );
        if( error ){
            _error = new errorNode( "initialize", "Error in data-file interface initialization" );
            _error->addNext( error );
            return;
        }

        //Initialize the coupling domains
        error = initializeCouplingDomains( );
        if ( error ){
            _error = new errorNode( "initialize", "Error in initialization of the coupling domains" );
            _error->addNext( error );
            return;
        }

        //Check the coupling initialization
        error = checkCouplingInitialization( );
        if ( error ){
            _error = new errorNode( "initialize", "Error in the coupling initialization configuration" );
            _error->addNext( error );
            return;
        }

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

        //Set the new configuration filename
        _configFilename = configurationFilename;

        //Reset whether the increment has been initialized
        _increment_initialized = false;

        initialize( );

        return _error;
    }

    errorOut inputFileProcessor::openConfigurationFile( ){
        /*!
         * Open the configuration file
         */

        if ( _configFilename.compare( "" ) == 0 ){
            return new errorNode( "openConfigurationFile",
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

    errorOut inputFileProcessor::initializeFileInterfaces( ){
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

            _macroscale = dataFileInterface::dataFileBase( _config[ "macroscale_definition" ] ).create( _config[ "macroscale_definition" ][ "filetype" ].as< std::string >( ) );

            if ( _macroscale->_error ){
                return _macroscale->_error;
            }

        }
        else{

            return new errorNode( "initializeFileInterfaces",
                                  "There is no 'macroscale_definition' in the YAML configuration file" );

        }
        if ( _config[ "microscale_definition" ] ){

            //Set the Default values
            if ( !_config[ "microscale_definition" ]["mode"] ){

                _config[ "microscale_definition" ][ "mode" ] = "read";

            }
            if ( !_config[ "microscale_definition" ]["filetype"] ){

                _config[ "microscale_definition" ][ "filetype" ] = "XDMF";

            }

            //Open the micro-scale data file
            _microscale = dataFileInterface::dataFileBase( _config[ "microscale_definition" ] ).create( _config[ "microscale_definition" ][ "filetype" ].as< std::string >( ) );

            if ( _microscale->_error ){
                return _microscale->_error;
            }

        }
        else{

            return new errorNode( "initializeFileInterfaces",
                                  "There is no 'microscale_definition' in the YAML configuration file" );

        }

        //Check that all required sets in the configuration file are defined correctly in the datasets
        
        //Check that micro-scale volume-surface node pairs are defined correctly
        if ( _config[ "microscale_definition" ][ "volume_surface_node_pairs" ] ){

            YAML::Node VSPairs = _config[ "microscale_definition" ][ "volume_surface_node_pairs" ];

            if ( !VSPairs.IsSequence( ) ){
                return new errorNode( "initializeFileInterfaces",
                                      "'volume_surface_node_pairs' is required to be a sequence" );
            }

            unsigned int indx = 1;
            for ( auto it = VSPairs.begin( ); it != VSPairs.end(); it++ ){

                if ( !( *it )[ "volume" ] ){
                    return new errorNode( "initializeFileInterfaces",
                                          "'volume' must be defined for each entry of 'volume_surface_node_pairs'" );
                }
                else if( !( *it )[ "volume" ].IsScalar( ) ){
                    return new errorNode( "initializeFileInterface",
                                          "'volume' should only define a single nodeset ( pair " + std::to_string( indx ) + " )" );
                }
                if ( !( *it )[ "surface" ] ){
                    return new errorNode( "initializeFileInterfaces",
                                          "'surface' must be defined for each entry of 'volume_surface_node_pairs'" );
                }
                else if( !( *it )[ "surface" ].IsScalar( ) ){
                    return new errorNode( "initializeFileInterface",
                                          "'surface' should only define a single nodeset ( pair " + std::to_string( indx ) + " )" );
                }

                indx++;

            }

        }

        return NULL;
    }

    errorOut inputFileProcessor::setMicroNodeWeights( const unsigned int increment ){
        /*!
         * Compute the weights of the micro-nodes
         *
         * :param const unsigned int increment: The increment at which to compute the weights
         */

        //Get the number of nodes in the microscale
        unsigned int numNodes;
        errorOut error = _microscale->getNumNodes( increment, numNodes );

        if ( error ){
            errorOut result = new errorNode( "setMicroNodeWeights", "Error in the extraction of the number of micro-scale nodes" );
            result->addNext( error );
            return result;
        }

        //Initialize the weight vector
        _microDomainWeights = floatVector( numNodes, 0 );

        //Loop through the free micro-surface sets
        uIntVector setNodes;
        for ( auto setName = _free_micro_surface_sets.begin(); setName != _free_micro_surface_sets.end(); setName++ ){

            //Get the surface set
            error = _microscale->getDomainNodes( increment, *setName, setNodes );

            if ( error ){
                errorOut result = new errorNode( "setMicroNodeWeights",
                                                 "Error in the extraction of the free micro surface set " + *setName );
                result->addNext( error );
                return result;
            }

            //Loop over the nodes
            for ( auto n = setNodes.begin(); n != setNodes.end(); n++ ){
                
                _microDomainWeights[ *n ] += 1;

            }

        }

        //Loop through the ghost micro-surface sets
        for ( auto setName = _ghost_micro_surface_sets.begin(); setName != _ghost_micro_surface_sets.end(); setName++ ){

            //Get the surface set
            error = _microscale->getDomainNodes( increment, *setName, setNodes );

            if ( error ){
                errorOut result = new errorNode( "setMicroNodeWeights",
                                                 "Error in the extraction of the ghost micro surface set " + *setName );
                result->addNext( error );
                return result;
            }

            //Loop over the nodes
            for ( auto n = setNodes.begin(); n != setNodes.end(); n++ ){
                
                _microDomainWeights[ *n ] += 1;

            }

        }

        //Loop through the non-overlapped sets
        for ( auto setName = _non_overlapped_micro_surface_sets.begin();
              setName != _non_overlapped_micro_surface_sets.end(); setName++ ){

            //Get the surface set
            error = _microscale->getDomainNodes( increment, *setName, setNodes );

            if ( error ){
                errorOut result = new errorNode( "setMicroNodeWeights",
                                                 "Error in the extraction of the non-overlapped micro surface set " + *setName );
                result->addNext( error );
                return result;
            }

            //Loop over the nodes
            for ( auto n = setNodes.begin(); n != setNodes.end(); n++ ){
                
                _microDomainWeights[ *n ] += 1;

            }

        }

        //Compute the weights
        for ( unsigned int i = 0; i < numNodes; i++ ){

            if ( _microDomainWeights[ i ] > 0 ){

                _microDomainWeights[ i ] = 1. / _microDomainWeights[ i ];

            }
            else{

                _microDomainWeights[ i ] = 1.;

            }

        }

        return NULL;
    }

    errorOut inputFileProcessor::setSurfaceSets( const unsigned int increment ){
        /*!
         * Set the surface sets
         *
         * :param const unsigned int increment: The increment to extract the surface sets
         */

        //Extract the set names from the microscale simulation
        stringVector setNames;
        errorOut error = _microscale->getSetNames( increment, setNames );

        if ( error ){
            errorOut result = new errorNode( "setSurfaceSets",
                                             "Error in extraction of the current increment's set names" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut inputFileProcessor::initializeIncrement( const unsigned int increment ){
        /*!
         * Initialize the processor for the indicated increment.
         *
         * :param const unsigned int increment: The increment to prepare for.
         */

        //Check if the requested increment is the currently initialized increment
        //If so, we don't need to run the initialization
        if ( ( increment == _current_increment ) && ( _increment_initialized ) ){
            return NULL;
        }

        errorOut error;
        //Collect the sets
        error = setSurfaceSets( increment );

        //Set the weights of the micro-nodes
        error = setMicroNodeWeights( increment );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in computation of the micro-node weights" );
            result->addNext( error );
            return result;
        }

        //Extract the densities of the micro-nodes
        error = extractMicroNodeDensities( increment );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extraction of the micro-node densities" );
            result->addNext( error );
            return result;
        }

        //Extract the volumes of the micro-nodes
        error = extractMicroNodeVolumes( increment );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extraction of the micro-node volumes" );
            result->addNext( error );
            return result;
        }

        //Set the current increment
        _current_increment = increment;
        _increment_initialized = true;

        return NULL;
    }

    errorOut inputFileProcessor::initializeCouplingDomains( ){
        /*!
         * Initialize the coupling domains
         */

        if ( _config[ "free_macroscale_domains" ] ){

            errorOut error = checkCommonDomainConfiguration( _config[ "free_macroscale_domains" ], _ghost_micro_surface_sets );
            if ( error ){
                errorOut result = new errorNode( "initializeCouplingDomains",
                                                 "Error in input-file check of the free macroscale domains" );
                result->addNext( error );
                return result;
            }

        }

        if ( _config[ "ghost_macroscale_domains" ] ){

            errorOut error = checkCommonDomainConfiguration( _config[ "ghost_macroscale_domains" ], _free_micro_surface_sets );
            if ( error ){
                errorOut result = new errorNode( "initializeCouplingDomains",
                                                 "Error in input-file check of the ghost macroscale domains" );
                result->addNext( error );
                return result;
            }

        }

        if ( _config[ "non_overlapped_microscale_domains" ] ){

            if ( _config[ "non_overlapped_microscale_domains" ].IsSequence( ) ){

                unsigned int nNonOverlappedMicroscaleDomains = 0;
                unsigned int indx = 1;
                _non_overlapped_micro_surface_sets.clear();

                for ( auto domain = _config[ "non_overlapped_microscale_domains" ].begin( );
                      domain != _config[ "non_overlapped_microscale_domains" ].end( );
                      domain++ ){

                    if ( !domain->IsScalar( ) ){

                        return new errorNode( "initializeCouplingDomains",
                                              "Entry " + std::to_string( indx ) + " of non_overlapped_microscale_domains is not a scalar" );

                    }

                    nNonOverlappedMicroscaleDomains++;
                    indx++;

                }

                indx = 0;
                _non_overlapped_micro_surface_sets = stringVector( nNonOverlappedMicroscaleDomains );

                for ( auto domain = _config[ "non_overlapped_microscale_domains" ].begin( );
                      domain != _config[ "non_overlapped_microscale_domains" ].end( );
                      domain++ ){

                    _non_overlapped_micro_surface_sets[ indx ] = domain->as< std::string >( );
                    indx++;

                }

            }

        }

        //Make sure that no nodeset that appears in the ghost micro-scale also appears in the free or non-overlapped micro-scales
        for ( auto nodeset = _ghost_micro_surface_sets.begin( ); nodeset != _ghost_micro_surface_sets.end( ); nodeset++ ){

            if ( std::find( _free_micro_surface_sets.begin( ), _free_micro_surface_sets.end( ), *nodeset ) != _free_micro_surface_sets.end( ) ){
                return new errorNode( "initializeCouplingDomains",
                                      *nodeset + " appears in the ghost and free micro-surface nodeset definitions" );
            }

        }

        return NULL;

    }

    errorOut inputFileProcessor::checkCommonDomainConfiguration( const YAML::Node &domainConfig,
                                                                 stringVector &volumeNodesets ){
        /*!
         * Extract common values in the configuration of a domain
         *
         * :param YAML::Node &domainConfig: The configuration of a particular domain.
         * :param stringVector &volumeNodesets: The nodeset names for the surfaces of the
         *     micro domains
         */

        if ( ( !domainConfig.IsSequence() ) && ( !domainConfig.IsNull( ) ) ){

            std::cout << domainConfig[ "macro_free_1" ][ "macro_nodeset" ].as<std::string>( ) <<"\n";

            return new errorNode( "checkCommonDomainConfiguration",
                                  "The definition of the domains must either be empty or a sequence" );

        }

        unsigned int indx = 1;
        unsigned int indx2 = 1;
        unsigned int nVolumeNodesets = 0;
        volumeNodesets.clear();
        for ( auto domain = domainConfig.begin( ); domain != domainConfig.end(); domain++ ){

            if ( !( *domain )[ "macro_nodeset" ] ){

                return new errorNode( "checkCommonDomainConfiguration",
                                      "The macro-nodeset is not defined in entry " + std::to_string( indx ) );

            }
            if ( ( *domain )[ "macro_nodeset" ].IsNull( ) ){
                return new errorNode( "checkCommonDomainConfiguration",
                                      "'macro_nodeset' cannot be empty in entry " + std::to_string( indx ) );
            }

            if ( !( *domain )[ "micro_nodesets" ] ){

                return new errorNode( "checkCommonDomainConfiguration",
                                      "The micro-nodeset is not defined in entry " + std::to_string( indx ) );

            }
            if ( !( *domain )[ "micro_nodesets" ].IsSequence( ) ){

                return new errorNode( "checkCommonDomainConfiguration",
                                      "The micro-nodesets are not defined as a sequence in entry " + std::to_string( indx ) ); 

            }

            indx2 = 1;
            for ( auto nodeset = ( *domain )[ "micro_nodesets" ].begin( ); nodeset != ( *domain )[ "micro_nodesets" ].end( ); nodeset++ ){

                if ( !nodeset->IsScalar( ) ){

                    return new errorNode( "checkCommonDomainConfiguration",
                                          "Micro-nodeset entry " + std::to_string( indx2 ) + " of domain entry " + std::to_string( indx ) + " is not a Scalar" );

                }

                nVolumeNodesets++;
                indx2++;

            }

            indx++;

        }

        //Extract the surface nodesets
        volumeNodesets.reserve( nVolumeNodesets );
        for ( auto domain = domainConfig.begin( ); domain != domainConfig.end( ); domain++ ){
            
            for ( auto nodeset = ( *domain )[ "micro_nodesets" ].begin( ); nodeset != ( *domain )[ "micro_nodesets" ].end( ); nodeset++ ){

                if ( std::find( volumeNodesets.begin( ), volumeNodesets.end( ), nodeset->as< std::string >( ) ) != volumeNodesets.end( ) ){

                    return new errorNode( "checkCommonDomainConfiguration",
                                          nodeset->as< std::string >( ) + " appears more than once in the coupling definition" );

                }

                volumeNodesets.push_back( nodeset->as< std::string >( ) );

            }

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroNodeDensities( const unsigned int &increment ){
        /*!
         * Extract the node densities for the micro domain at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        //Check if the density name has been defined
        if ( !_config[ "microscale_definition" ][ "density_variable_name" ] ){
        
            return new errorNode( "extractMicroNodeDensities", "The density variable name is not defined" );

        }

        //Re-size the micro node density vector
        errorOut error = _microscale->getMeshData( increment, 
                                                   _config[ "microscale_definition" ][ "density_variable_name" ].as< std::string >( ),
                                                   "Node", _microDensity );

        if ( error ){

            errorOut result = new errorNode( "extractMicroNodeDensities", "Error in extraction of the micro densitys" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroNodeVolumes( const unsigned int &increment ){
        /*!
         * Extract the node volumes for the micro domain at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        //Check if the volume name has been defined
        if ( !_config[ "microscale_definition" ][ "volume_variable_name" ] ){
        
            return new errorNode( "extractMicroNodeVolumes", "The volume variable name is not defined" );

        }

        //Re-size the micro node density vector
        errorOut error = _microscale->getMeshData( increment, 
                                                   _config[ "microscale_definition" ][ "volume_variable_name" ].as< std::string >( ),
                                                   "Node", _microVolume );

        if ( error ){

            errorOut result = new errorNode( "extractMicroNodeVolumes", "Error in extraction of the micro volumes" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    const floatVector* inputFileProcessor::getMicroDensities( ){
        /*!
         * Get a pointer to the density
         */

        return &_microDensity;
    }

    const floatVector* inputFileProcessor::getMicroVolumes( ){
        /*!
         * Get a pointer to the volumes
         */

        return &_microVolume;
    }

    const floatVector* inputFileProcessor::getMicroWeights( ){
        /*!
         * Get the micro weights
         */

        return &_microDomainWeights;
    }

    const stringVector* inputFileProcessor::getFreeMicroDomainNames( ){
        /*!
         * Get the free domain names
         */

        return &_free_micro_surface_sets;
    }

    const stringVector* inputFileProcessor::getGhostMicroDomainNames( ){
        /*!
         * Get the ghost domain names
         */

        return &_ghost_micro_surface_sets;
    }

    const stringVector* inputFileProcessor::getNonOverlappedMicroDomainNames( ){
        /*!
         * Get the ghost domain names
         */

        return &_non_overlapped_micro_surface_sets;
    }

    const YAML::Node inputFileProcessor::getCouplingInitialization( ){
        /*!
         * Get the coupling initialization from the configuration file.
         */

        return _config[ "coupling_initialization" ];
    }

    errorOut inputFileProcessor::checkCouplingInitialization( ){
        /*!
         * Check the coupling initialization
         */

        if ( !_config[ "coupling_initialization" ][ "type" ] ){

            _config[ "coupling_initialization" ][ "type" ] = "use_initial_state";

        }

        return NULL;

    }

}
