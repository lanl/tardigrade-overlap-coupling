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

        //Check the volume reconstruction initialization
        error = checkVolumeReconstructionInitialization( );
        if ( error ){
            _error = new errorNode( "initialize", "Error in the volume reconstruction initialization" );
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

    errorOut inputFileProcessor::setSurfaceSets( const unsigned int microIncrement ){
        /*!
         * Set the surface sets
         *
         * :param const unsigned int microIncrement: The increment to extract the micro surface sets
         */

        //Extract the set names from the microscale simulation
        stringVector setNames;
        errorOut error = _microscale->getSetNames( microIncrement, setNames );

        if ( error ){
            errorOut result = new errorNode( "setSurfaceSets",
                                             "Error in extraction of the current micro increment's set names" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut inputFileProcessor::initializeIncrement( const unsigned int microIncrement, const unsigned int macroIncrement ){
        /*!
         * Initialize the processor for the indicated increment.
         *
         * :param const unsigned int microIncrement: The micro increment to prepare for.
         * :param const unsigned int macroIncrement: The micro increment to prepare for.
         */

        //Check if the requested increment is the currently initialized increment
        //If so, we don't need to re-run the initialization
        if ( ( macroIncrement == _current_macroIncrement ) &&
             ( microIncrement == _current_microIncrement ) &&
             ( _increment_initialized ) ){
            return NULL;
        }

        errorOut error = NULL;
        //Collect the sets
        error = setSurfaceSets( microIncrement );

        //Set the weights of the micro-nodes
        error = setMicroNodeWeights( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in computation of the micro-node weights" );
            result->addNext( error );
            return result;
        }

        //Extract the densities of the micro-nodes
        error = extractMicroNodeDensities( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extraction of the micro-node densities" );
            result->addNext( error );
            return result;
        }

        //Extract the volumes of the micro-nodes
        error = extractMicroNodeVolumes( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extraction of the micro-node volumes" );
            result->addNext( error );
            return result;
        }

        //Extract the reference positions of the micro-nodes
        error = extractReferenceMicroMeshData( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extraction of the micro-node mesh information" );
            result->addNext( error );
            return result;
        }

        //Extract the reference positions of the macro-nodes
        error = extractReferenceMacroMeshData( macroIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extraction of the macro-node mesh information" );
            result->addNext( error );
            return result;
        }

        //Extract the micro displacements
        error = extractMicroDisplacements( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extraction of the micro displacements" );
            result->addNext( error );
            return result;
        }

        //Extract the micro body forces
        error = extractMicroBodyForces( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the micro body forces" );
            result->addNext( error );
            return result;
        }

        //Extract the micro accelerations
        error = extractMicroAccelerations( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the micro accelerations" );
            result->addNext( error );
            return result;
        }

        //Extract the micro stresses
        error = extractMicroStresses( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the micro stresses" );
            result->addNext( error );
            return result;
        }

        //Extract the macro displacements
        error = extractMacroDisplacements( macroIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extraction of the macro displacements" );
            result->addNext( error );
            return result;
        }

        //Extract the macro displacement DOF vector
        error = extractMacroDispDOFVector( macroIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extraction of the macro displacement DOF vector" );
            result->addNext( error );
            return result;
        }

        //Set the unique macro and micro nodes
        error = setMicroNodeIndexMappings( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in setting the unique micro node index mappings" );
            result->addNext( error );
            return result;
        }

        error = setMacroNodeIndexMappings( macroIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in setting the unique macro node index mappings" );
            result->addNext( error );
            return result;
        }

        //Set the current increment
        _current_macroIncrement = macroIncrement;
        _current_microIncrement = microIncrement;
        _increment_initialized = true;

        return NULL;
    }

    errorOut inputFileProcessor::initializeCouplingDomains( ){
        /*!
         * Initialize the coupling domains
         */

        if ( _config[ "free_macroscale_domains" ] ){

            errorOut error = checkCommonDomainConfiguration( _config[ "free_macroscale_domains" ],
                                                             _free_macro_cell_ids, _free_macro_cell_micro_domain_counts,
                                                             _free_macro_volume_sets,
                                                             _ghost_micro_volume_sets );
            if ( error ){
                errorOut result = new errorNode( "initializeCouplingDomains",
                                                 "Error in input-file check of the free macroscale domains" );
                result->addNext( error );
                return result;
            }

            error = checkCommonVolumeToSurfaceMapping( _ghost_micro_volume_sets, _ghost_micro_surface_sets );
            if ( error ){
                errorOut result = new errorNode( "initializeCouplingDomains",
                                                 "Error in input-file check of the micro volume to surface mapping for the free macroscale domains" );
                result->addNext( error );
                return result;
            }

        }

        if ( _config[ "ghost_macroscale_domains" ] ){

            errorOut error = checkCommonDomainConfiguration( _config[ "ghost_macroscale_domains" ],
                                                             _ghost_macro_cell_ids, _ghost_macro_cell_micro_domain_counts,
                                                             _ghost_macro_volume_sets,
                                                             _free_micro_volume_sets );
            if ( error ){
                errorOut result = new errorNode( "initializeCouplingDomains",
                                                 "Error in input-file check of the ghost macroscale domains" );
                result->addNext( error );
                return result;
            }

            error = checkCommonVolumeToSurfaceMapping( _free_micro_volume_sets, _free_micro_surface_sets );
            if ( error ){
                errorOut result = new errorNode( "initializeCouplingDomains",
                                                 "Error in input-file check of the micro volume to surface mapping for the free macroscale domains" );
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

            if ( std::find( _non_overlapped_micro_surface_sets.begin( ), _non_overlapped_micro_surface_sets.end( ), *nodeset ) != _non_overlapped_micro_surface_sets.end( ) ){
                return new errorNode( "initializeCouplingDomains",
                                      *nodeset + " appears in the ghost and free micro-surface nodeset definitions" );
            }
        }

        for ( auto nodeset = _free_micro_surface_sets.begin( ); nodeset != _free_micro_surface_sets.end( ); nodeset++ ){

            if ( std::find( _non_overlapped_micro_surface_sets.begin( ), _non_overlapped_micro_surface_sets.end( ), *nodeset ) != _non_overlapped_micro_surface_sets.end( ) ){
                return new errorNode( "initializeCouplingDomains",
                                      *nodeset + " appears in the free and non-overlapped micro-surface nodeset definitions" );
            }

        }

        //Make sure that no volume nodeset that appears in the ghost micro-scale also appears in the free micro-scale
        for ( auto nodeset = _ghost_micro_volume_sets.begin( ); nodeset != _ghost_micro_volume_sets.end( ); nodeset++ ){

            if ( std::find( _free_micro_volume_sets.begin( ), _free_micro_volume_sets.end( ), *nodeset ) != _free_micro_volume_sets.end( ) ){
                return new errorNode( "initializeCouplingDomains",
                                      *nodeset + " appears in the ghost and free micro-volume nodeset definitions" );
            }

        }

        return NULL;

    }

    errorOut inputFileProcessor::checkCommonDomainConfiguration( const YAML::Node &domainConfig,
                                                                 uIntVector &macroCellIds,
                                                                 uIntVector &macroCellMicroDomainCounts,
                                                                 stringVector &macroVolumeNodesets,
                                                                 stringVector &microVolumeNodesets ){
        /*!
         * Extract common values in the configuration of a domain
         *
         * :param YAML::Node &domainConfig: The configuration of a particular domain.
         * :param uIntVector &macroCellIds: The macro-cell Ids corresponding to the domain
         * :param uIntVector &macroCellMicroDomainCounts: The number of micro domains in each macro domain
         * :param stringVector &macroVolumeNodesets: The nodeset names for the nodes in the
         *     macro domains
         * :param stringVector &microVolumeNodesets: The nodeset names for the nodes in the
         *     micro domains
         */

        if ( ( !domainConfig.IsSequence() ) && ( !domainConfig.IsNull( ) ) ){

            return new errorNode( "checkCommonDomainConfiguration",
                                  "The definition of the domains must either be empty or a sequence" );

        }

        unsigned int indx = 1;
        unsigned int indx2 = 1;
        unsigned int nVolumeNodesets = 0;
        microVolumeNodesets.clear();
        for ( auto domain = domainConfig.begin( ); domain != domainConfig.end(); domain++ ){

            if ( !( *domain )[ "macro_nodeset" ] ){

                return new errorNode( "checkCommonDomainConfiguration",
                                      "The macro-nodeset is not defined in entry " + std::to_string( indx ) );

            }
            if ( !( *domain )[ "macro_cell" ].IsScalar( ) ){
                return new errorNode( "checkCommonDomainConfiguration",
                                      "'macro_cell' must be defined as the cell ( element ) corresponding with the nodeset. It is empty in entry " + std::to_string( indx ) );
            }
            if ( ( *domain )[ "macro_nodeset" ].IsNull( ) ){
                return new errorNode( "checkCommonDomainConfiguration",
                                      "'macro_nodeset' cannot be empty in entry " + std::to_string( indx ) );
            }

            if ( !( *domain )[ "macro_nodeset" ] ){
                return new errorNode( "checkCommonDomainConfiguration",
                                      "The macro-nodeset is not defined in entry " + std::to_string( indx ) );
            }

            if ( !( *domain )[ "macro_nodeset" ].IsScalar( ) ){
                return new errorNode( "checkCommonDomainConfiguration",
                                      "The macro-nodeset must be a scalar string value " + std::to_string( indx ) );
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

        //Extract the volume nodesets
        macroCellIds.reserve( domainConfig.size( ) );
        macroCellMicroDomainCounts.reserve( domainConfig.size( ) );
        macroVolumeNodesets.reserve( domainConfig.size( ) );
        microVolumeNodesets.reserve( nVolumeNodesets );
        for ( auto domain = domainConfig.begin( ); domain != domainConfig.end( ); domain++ ){

            macroCellIds.push_back( ( *domain )[ "macro_cell" ].as< unsigned int >( ) );
            macroVolumeNodesets.push_back( ( *domain )[ "macro_nodeset" ].as< std::string >( ) );
            macroCellMicroDomainCounts.push_back( ( *domain )[ "micro_nodesets" ].size( ) );

            for ( auto nodeset = ( *domain )[ "micro_nodesets" ].begin( ); nodeset != ( *domain )[ "micro_nodesets" ].end( ); nodeset++ ){

                if ( std::find( microVolumeNodesets.begin( ), microVolumeNodesets.end( ), nodeset->as< std::string >( ) ) != microVolumeNodesets.end( ) ){

                    return new errorNode( "checkCommonDomainConfiguration",
                                          nodeset->as< std::string >( ) + " appears more than once in the coupling definition" );

                }

                microVolumeNodesets.push_back( nodeset->as< std::string >( ) );

            }

        }

        return NULL;

    }

    errorOut inputFileProcessor::checkCommonVolumeToSurfaceMapping( const stringVector &microVolumeNodesets, 
                                                                    stringVector &microSurfaceNodesets ){
        /*!
         * Check the common volume to surface node mapping
         *
         * :param const stringVector &microVolumeNodesets: The nodesets which represent the volume of the micro-domains
         * :param const stringVector &microSurfaceNodesets: The nodesets which represent the surfaces of the micro-domains
         */

        if ( !_config[ "microscale_definition" ][ "volume_surface_node_pairs" ] ){
            
            return new errorNode( "checkCommonVolumeToSurfaceMapping",
                                  "volume_surface_node_pairs must be defined as a keyname attribute of 'microscale_definition'" );

        }

        if ( ( !_config[ "microscale_definition" ][ "volume_surface_node_pairs" ].IsSequence( ) ) && ( microVolumeNodesets.size( ) > 0 ) ){
            return new errorNode( "checkCommonVolumeToSurfaceMapping",
                                  "volume_surface_node_pairs must be a sequence" );

        }

        //Loop through the micro-volume nodesets
        microSurfaceNodesets = stringVector( microVolumeNodesets.size( ) );
        unsigned int indx = 0;
        unsigned int p;
        bool volumeFound = false;
        for ( auto volumeName = microVolumeNodesets.begin( ); volumeName != microVolumeNodesets.end( ); volumeName++ ){

            p = 1;
            volumeFound = false;
            for ( auto pair = _config[ "microscale_definition" ][ "volume_surface_node_pairs" ].begin( );
                       pair != _config[ "microscale_definition" ][ "volume_surface_node_pairs" ].end( );
                  pair++ ){

                if ( !( *pair )[ "volume" ] ){

                    return new errorNode( "checkCommonVolumeToSurfaceMapping",
                                          "'volume' is not a key of 'volume_surface_node_pairs' entry " + std::to_string( p ) );

                }
                if ( !( *pair )[ "surface" ] ){

                    return new errorNode( "checkCommonVolumeToSurfaceMapping",
                                          "'surface' is not a key of 'volume_surface_node_pairs' entry " + std::to_string( p ) );

                }

                if ( ( *pair )[ "volume" ].as< std::string >( ).compare( *volumeName ) == 0 ){

                    microSurfaceNodesets[ indx ] = ( *pair )[ "surface" ].as< std::string >( );
                    volumeFound = true;
                    break;

                }

                p++;

            }

            if ( !volumeFound ){

                return new errorNode( "checkCommonVolumeToSurfaceMapping"
                                      "'" + *volumeName + "' not found in 'volume_surface_node_pairs'" );

            }

            indx++;

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
        errorOut error = _microscale->getSolutionData( increment, 
                                                       _config[ "microscale_definition" ][ "density_variable_name" ].as< std::string >( ),
                                                       "Node", _microDensities );

        if ( error ){

            errorOut result = new errorNode( "extractMicroNodeDensities", "Error in extraction of the micro densities" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroBodyForces( const unsigned int &increment ){
        /*!
         * Extract the node micro-body forces at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        //Check if the body force name has been defined
        _microBodyForceFlag = false;
        if ( ( !_config[ "microscale_definition" ][ "body_force_variable_name" ] ) ||
             ( _config[ "microscale_definition" ][ "body_force_variable_name" ].as< std::string >( ).compare( "NULL" ) == 0 ) ){

            _config [ "microscale_definition" ][ "body_force_variable_name" ] = "NULL"; //Indicate that the body force is assumed to be zero
            _microBodyForces = { 0., 0., 0. }; //Set the body force to zero

            return NULL;

        }

        //Extract the micro body force vector
        errorOut error =
            _microscale->getSolutionData( increment,
                                          _config[ "microscale_definition" ][ "body_force_variable_name" ].as< std::string > ( ),
                                          "Node", _microBodyForces );

        if ( error ){

            errorOut result = new errorNode( "extractMicroBodyForce", "Error in extraction of the micro body forces" );
            result->addNext( error );
            return result;

        }

        _microBodyForceFlag = true;

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroAccelerations( const unsigned int &increment ){
        /*!
         * Extract the node micro-body forces at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        //Check if the body force name has been defined
        _microAccelerationFlag = false;
        if ( ( !_config[ "microscale_definition" ][ "acceleration_variable_name" ] ) ||
             ( _config[ "microscale_definition" ][ "acceleration_variable_name" ].as< std::string >( ).compare( "NULL" ) == 0 ) ){

            _config [ "microscale_definition" ][ "acceleration_variable_name" ] = "NULL"; //Indicate that the acceleration is assumed to be zero
            _microAccelerations = { 0., 0., 0. }; //Set the acceleration to zero

            return NULL;

        }

        //Extract the micro acceleration vector
        errorOut error =
            _microscale->getSolutionData( increment,
                                          _config[ "microscale_definition" ][ "acceleration_variable_name" ].as< std::string > ( ),
                                          "Node", _microAccelerations );

        if ( error ){

            errorOut result = new errorNode( "extractMicroAccelerations", "Error in extraction of the micro accelerations" );
            result->addNext( error );
            return result;

        }

        _microAccelerationFlag = true;

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroStresses( const unsigned int &increment ){
        /*!
         * Extract the micro-stresses at the indicated increment
         *
         * :param const unsigne int &increment: The current increment
         */

        if ( !_config[ "microscale_definition" ][ "stress_variable_names" ] ){

            return new errorNode( "extractMicroStresses",
                                  "The stress variable names have not been defined under microscale_definition -> stress_variable_names" );

        }

        uIntType dim = ( uIntType )std::sqrt( ( floatType )_config[ "microscale_definition" ][ "stress_variable_names" ].size( ) );

        if ( _config[ "microscale_definition" ][ "stress_variable_names" ].size( ) != dim * dim ){

            return new errorNode( "extractMicroStresses",
                                  "The dimensionality of the stresses must be a perfect square ( 1 value, 4 values, or 9 values ) for 1d, 2d, or 3d" );

        }

        stringVector stressKeys( dim * dim );

        for ( uIntType i = 0; i < dim; i++ ){

            for ( uIntType j = 0; j < dim; j++ ){

                stressKeys[ dim * i + j ] = "s" + std::to_string( i + 1 ) + std::to_string( j + 1 );

            }

        }

        uIntType numNodes;
        errorOut error = _microscale->getNumNodes( increment, numNodes );

        if ( error ){

            errorOut result = new errorNode( "extractMicroStresses",
                                             "Error in getting the number of nodes in the microscale" );
            result->addNext( error );
            return result;

        }

        _microStresses = floatVector( dim * dim * numNodes, 0 );

        for ( auto sK = stressKeys.begin( ); sK != stressKeys.end( ); sK++ ){

            if ( !_config[ "microscale_definition" ][ "stress_variable_names" ][ *sK ] ){

                return new errorNode( "extractMicroStresses",
                                      "The micro stress component " + *sK + " is not defined in the input file" );

            }

            floatVector stressComponent;

            error = _microscale->getSolutionData( increment,
                                                  _config[ "microscale_definition" ][ "stress_variable_names" ][ *sK ].as< std::string >( ),
                                                  "Node", stressComponent );

            if ( error ){

                errorOut result = new errorNode( "extractMicroStresses",
                                                 "Error in extracting the stress component " + *sK + " from the input file" );
                result->addNext( error );
                return result;

            }

            for ( auto sC = stressComponent.begin( ); sC != stressComponent.end( ); sC++ ){

                _microStresses[ dim * dim * ( sC - stressComponent.begin( ) ) + sK - stressKeys.begin( ) ] = *sC;

            }

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
        errorOut error = _microscale->getSolutionData( increment, 
                                                       _config[ "microscale_definition" ][ "volume_variable_name" ].as< std::string >( ),
                                                       "Node", _microVolumes );

        if ( error ){

            errorOut result = new errorNode( "extractMicroNodeVolumes", "Error in extraction of the micro volumes" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroDisplacements( const unsigned int &increment ){
        /*!
         * Extract the positions of the nodes in the micro domain.
         *
         * :param const unsigned int &increment: The current increment
         */

        //Variable definitions
        long unsigned int indx = 0;

        floatVector disp;

        //Expected displacement keys
        stringVector displacementKeys = { "u1", "u2", "u3" }; //The dimension is assumed to be 3

        //Initialize the micro displacements vector
        _microDisplacements.clear( );
        unsigned int numNodes;

        errorOut error = _microscale->getNumNodes( increment, numNodes );

        if ( error ){
            
            errorOut result = new errorNode( "extractMicroDisplacements", "Error in the determination of the number of micro-nodes" );
            result->addNext( error );
            return result;

        }

        _microDisplacements.resize( displacementKeys.size( ) * numNodes );

        //Check if the displacement names have been defined
        if ( !_config[ "microscale_definition" ][ "displacement_variable_names" ] ){

            return new errorNode( "extractMicroDisplacements", "The names of the displacements of the micro-positions are not defined" );

        }
        else{

            for ( auto it = displacementKeys.begin( ); it != displacementKeys.end( ); it++ ){

                if ( !_config[ "microscale_definition" ][ "displacement_variable_names" ][ *it ] ){

                    return new errorNode( "extractMicroDisplacements", "'" + *it + " is not defined in 'displacement_variable_names'");

                }
                else if ( !_config[ "microscale_definition" ][ "displacement_variable_names" ][ *it ].IsScalar( ) ){

                    return new errorNode( "extractMicroDisplacements", "'" + *it + "' is not a scalar" );

                }
                else{

                    //Get the displacement
                    errorOut error =
                        _microscale->getSolutionData( increment, 
                                                      _config[ "microscale_definition" ][ "displacement_variable_names" ][ *it ].as< std::string >( ),
                                                      "Node", disp );

                    if ( error ){
                        
                        errorOut result = new errorNode( "extractMicroDisplacements", "Error in extraction of '" + *it + "'" );
                        result->addNext( error );
                        return result;

                    }

                    //Set the displacements
                    for ( unsigned int i = 0; i < disp.size( ); i++ ){

                        _microDisplacements[ displacementKeys.size( ) * i + indx ] = disp[ i ];

                    }

                }

                indx++;

            }

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroDisplacements( const unsigned int &increment ){
        /*!
         * Extract the positions of the nodes in the macro domain.
         *
         * :param const unsigned int &increment: The current increment
         */

        //Variable definitions
        long unsigned int indx = 0;

        floatVector disp;

        //Expected displacement keys
        stringVector displacementKeys = { "u1", "u2", "u3" }; //The dimension is assumed to be 3

        //Initialize the micro displacements vector
        _macroDisplacements.clear( );
        unsigned int numNodes;

        errorOut error = _macroscale->getNumNodes( increment, numNodes );

        if ( error ){
            
            errorOut result = new errorNode( "extractMacroDisplacements", "Error in the determination of the number of macro-nodes" );
            result->addNext( error );
            return result;

        }

        _macroDisplacements.resize( displacementKeys.size( ) * numNodes );

        //Check if the displacement names have been defined
        if ( !_config[ "macroscale_definition" ][ "displacement_variable_names" ] ){

            return new errorNode( "extractMacroDisplacements", "The names of the displacements of the macro-positions are not defined" );

        }
        else{

            for ( auto it = displacementKeys.begin( ); it != displacementKeys.end( ); it++ ){

                if ( !_config[ "macroscale_definition" ][ "displacement_variable_names" ][ *it ] ){

                    return new errorNode( "extractMacroDisplacements", "'" + *it + " is not defined in 'displacement_variable_names'");

                }
                else if ( !_config[ "macroscale_definition" ][ "displacement_variable_names" ][ *it ].IsScalar( ) ){

                    return new errorNode( "extractMacroDisplacements", "'" + *it + "' is not a scalar" );

                }
                else{

                    //Get the displacement
                    errorOut error
                        = _macroscale->getSolutionData( increment, 
                                                        _config[ "macroscale_definition" ][ "displacement_variable_names" ][ *it ].as< std::string >( ),
                                                        "Node", disp );

                    if ( error ){
                        
                        errorOut result = new errorNode( "extractMacroDisplacements", "Error in extraction of '" + *it + "'" );
                        result->addNext( error );
                        return result;

                    }

                    //Set the displacements
                    for ( unsigned int i = 0; i < disp.size( ); i++ ){

                        _macroDisplacements[ displacementKeys.size( ) * i + indx ] = disp[ i ];

                    }

                }

                indx++;

            }

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroDispDOFVector( const unsigned int &increment ){
        /*!
         * Extract the displacement degrees of freedom of the nodes in the macro domain.
         *
         * :param const unsigned int &increment: The current increment
         */

        //Variable definitions
        long unsigned int indx = 0;

        floatVector disp;

        //Expected displacement keys
        stringVector displacementKeys = { "u1", "u2", "u3",
                                          "phi11", "phi12", "phi13", 
                                          "phi21", "phi22", "phi23", 
                                          "phi31", "phi32", "phi33" }; //The dimension is assumed to be 3

        //Initialize the micro displacement degree of freedom vector
        _macroDispDOFVector.clear( );
        unsigned int numNodes;

        errorOut error = _macroscale->getNumNodes( increment, numNodes );

        if ( error ){
            
            errorOut result = new errorNode( "extractMacroDispDOFVector", "Error in the determination of the number of macro-nodes" );
            result->addNext( error );
            return result;

        }

        _macroDispDOFVector.resize( displacementKeys.size( ) * numNodes );

        //Check if the displacement names have been defined
        if ( !_config[ "macroscale_definition" ][ "displacement_variable_names" ] ){

            return new errorNode( "extractMacroDispDOFVector", "The names of the displacements of the macro-positions are not defined" );

        }
        else{

            for ( auto it = displacementKeys.begin( ); it != displacementKeys.end( ); it++ ){

                if ( !_config[ "macroscale_definition" ][ "displacement_variable_names" ][ *it ] ){

                    return new errorNode( "extractMacroDispDOFVector", "'" + *it + " is not defined in 'displacement_variable_names'");

                }
                else if ( !_config[ "macroscale_definition" ][ "displacement_variable_names" ][ *it ].IsScalar( ) ){

                    return new errorNode( "extractMacroDispDOFVector", "'" + *it + "' is not a scalar" );

                }
                else{

                    //Get the displacement
                    errorOut error
                        = _macroscale->getSolutionData( increment, 
                                                        _config[ "macroscale_definition" ][ "displacement_variable_names" ][ *it ].as< std::string >( ),
                                                        "Node", disp );

                    if ( error ){
                        
                        errorOut result = new errorNode( "extractMacroDispDOFVector", "Error in extraction of '" + *it + "'" );
                        result->addNext( error );
                        return result;

                    }

                    //Set the displacements
                    for ( unsigned int i = 0; i < disp.size( ); i++ ){

                        _macroDispDOFVector[ displacementKeys.size( ) * i + indx ] = disp[ i ];

                    }

                }

                indx++;

            }

        }

        return NULL;

    }


    errorOut inputFileProcessor::extractReferenceMicroMeshData( const unsigned int &increment ){
        /*!
         * Extract the mesh data for the micro-scale domain in the reference configuration
         *
         * :param const unsigned int &increment: The increment at which to extract the micro-mesh data
         */

        errorOut error = _microscale->getMeshData( increment, _microNodeReferencePositions,
                                                   _microNodeReferenceConnectivity,
                                                   _microNodeReferenceConnectivityCellIndices, _microCellCounts );

        if ( error ){
            errorOut result = new errorNode( "extractMicroMeshData", "Error in the extraction of the micro-mesh information" );
            result->addNext( error );
            return result;
        }

        return NULL;

    }

    errorOut inputFileProcessor::extractReferenceMacroMeshData( const unsigned int &increment ){
        /*!
         * Extract the mesh data for the macro-scale domain in the reference configuration
         *
         * :param const unsigned int &increment: The increment at which to extract the macro-mesh data
         */

        errorOut error = _macroscale->getMeshData( increment, _macroNodeReferencePositions,
                                                   _macroNodeReferenceConnectivity, _macroNodeReferenceConnectivityCellIndices,
                                                   _macroCellCounts );

        if ( error ){
            errorOut result = new errorNode( "extractMacroMeshData", "Error in the extraction of the macro-mesh information" );
            result->addNext( error );
            return result;
        }

        return NULL;

    }

    const floatVector* inputFileProcessor::getMicroDensities( ){
        /*!
         * Get a pointer to the density
         */

        return &_microDensities;
    }

    const floatVector* inputFileProcessor::getMicroBodyForces( ){
        /*!
         * Get a pointer to the micro body forces
         */

        return &_microBodyForces;
    }

    const floatVector* inputFileProcessor::getMicroAccelerations( ){
        /*!
         * Get a pointer to the micro accelerations
         */

        return &_microAccelerations;
    }

    const floatVector* inputFileProcessor::getMicroStresses( ){
        /*!
         * Get a pointer to the micro stresses
         */

        return &_microStresses;
    }

    const floatVector* inputFileProcessor::getMicroVolumes( ){
        /*!
         * Get a pointer to the volumes
         */

        return &_microVolumes;
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

        return &_free_micro_volume_sets;
    }

    const stringVector* inputFileProcessor::getGhostMicroDomainNames( ){
        /*!
         * Get the ghost domain names
         */

        return &_ghost_micro_volume_sets;
    }

    const stringVector* inputFileProcessor::getFreeMicroSurfaceNames( ){
        /*!
         * Get the free domain names
         */

        return &_free_micro_surface_sets;
    }

    const stringVector* inputFileProcessor::getGhostMicroSurfaceNames( ){
        /*!
         * Get the ghost domain names
         */

        return &_ghost_micro_surface_sets;
    }

    const stringVector* inputFileProcessor::getNonOverlappedMicroSurfaceNames( ){
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

            _config[ "coupling_initialization" ][ "type" ] = "use_first_increment";

        }

        if ( !_config[ "coupling_initialization" ][ "projection_type" ] ){

            _config[ "coupling_initialization" ][ "projection_type" ] = "direct_projection";

        }

        if ( _config[ "coupling_initialization" ][ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){

            _computeMicroShapeFunctions = true;

        }

        if ( !_config[ "coupling_initialization" ][ "use_reconstructed_mass_centers" ] ){

            _config[ "coupling_initialization" ][ "use_reconstructed_mass_centers" ] = true;

        }

        return NULL;

    }

    errorOut inputFileProcessor::getUniqueNodesInDomains( const unsigned int &increment,
                                                          const std::shared_ptr< dataFileInterface::dataFileBase > &dataFile,
                                                          const stringVector &domainNames, uIntVector &uniqueIds ){
        /*!
         *
         * Determine the unique nodes in a collection of domains
         *
         * :param const std::shared_ptr< dataFileInterface::dataFileBase > &dataFile: The datafile to read
         * :param const stringVector &domainNames: The names of the domains to search
         * :param uIntVector &uniqueIds: The resulting unique Ids.
         *
         */

        //Loop over the free micro domains to determine the approximate
        //size of the micro-domains
        
        unsigned int approximateSize = 0;
        unsigned int n;

        for ( auto domain = domainNames.begin( );
                   domain != domainNames.end( );
                   domain++ ){

            dataFile->getNumDomainNodes( increment, *domain, n );
            approximateSize += n;

        }

        //Loop through the domains and store the unique node ids
        uniqueIds.clear( );
        uniqueIds.reserve( approximateSize );

        uIntVector nodes;
        errorOut error = NULL;

        for ( auto domain =  domainNames.begin( );
                   domain != domainNames.end( );
                   domain++ ){

            error = dataFile->getDomainNodes( increment, *domain, nodes );

            if ( error ){

                errorOut result = new errorNode( "getUniqueNodesInDomains", "Error in getting the nodes of '" + *domain + "'" );
                result->addNext( error );
                return result;

            }

            //Loop through the nodes
            for ( auto node =  nodes.begin( );
                       node != nodes.end( );
                       node++ ){
                  
                //If the node is not found add it 
                if ( std::find( uniqueIds.begin( ), uniqueIds.end( ), *node ) == uniqueIds.end( ) ){

                    uniqueIds.push_back( *node );

                }

            }

        }

        return NULL;

    }

    errorOut inputFileProcessor::setMicroNodeIndexMappings( const unsigned int &increment ){
        /*
         * Set the micro node index mappings for the output matrices.
         * The order that the node appears in the unique node vector
         * is the index
         *
         * If a micro-scale node is found in free then it cannot be ghost.
         * We give preference to free nodes since the micro-scale is assumed
         * to be more accurate than the macroscale.
         *
         * If a micro-scale node is found in the non-overlapped domains, then it
         * cannot be free. We give preference to non-overlapped nodes over 
         * free nodes since a micro-domain which is non-overlapped is assumed
         * to be more accurate than the intermediate domain
         */

        //Get the unique nodes in the non-overlapped, free and ghost domains
        errorOut error = getUniqueNodesInDomains( increment, _microscale, _non_overlapped_micro_surface_sets,
                                                  _unique_non_overlapped_micro_nodes );

        if ( error ){

            errorOut result = new errorNode( "setMicroNodeIndexMappings",
                                             "Error in determining the unique non-overlapped microscale nodes" );
            result->addNext( error );
            return result;

        }

        error = getUniqueNodesInDomains( increment, _microscale, _free_micro_volume_sets, _unique_free_micro_nodes );

        if ( error ){

            errorOut result = new errorNode( "setMicroNodeIndexMappings",
                                             "Error in determining the unique free microscale nodes" );
            result->addNext( error );
            return result;

        }

        error = getUniqueNodesInDomains( increment, _microscale, _ghost_micro_volume_sets, _unique_ghost_micro_nodes );

        if ( error ){

            errorOut result = new errorNode( "setMicroNodeIndexMappings",
                                             "Error in determining the unique ghost microscale nodes" );
            result->addNext( error );
            return result;

        }

        //Remove nodes found in the non-overlapped nodes from the free and ghost nodes
        unsigned int numNodes = 0;
        unsigned int n;

        //Approximate the size of the duplicate nodes. At worst, this will be the size of the
        //nodes on the surfaces of the non-overlapped domains
        for ( auto domain =  _non_overlapped_micro_surface_sets.begin( );
                   domain != _non_overlapped_micro_surface_sets.end( );
                   domain++ ){

            _microscale->getNumDomainNodes( increment, *domain, n );
            numNodes += n;

        }

        //Loop through the free nodes to find duplicates
        uIntVector duplicateNodes;
        duplicateNodes.reserve( numNodes );
        
        for ( auto node =  _unique_free_micro_nodes.begin( );
                   node != _unique_free_micro_nodes.end( );
                   node++ ){

            //If the free node is found in the non-overlapped nodes add it to the duplicates
            if ( std::find( _unique_non_overlapped_micro_nodes.begin( ), _unique_non_overlapped_micro_nodes.end( ),  *node )
                 != _unique_non_overlapped_micro_nodes.end( ) ){

                duplicateNodes.push_back( node - _unique_free_micro_nodes.begin( ) );

            }

        }

        //Sort the duplicate nodes
        std::sort( duplicateNodes.begin( ), duplicateNodes.end( ) );

        //Remove the duplicate nodes
        error = removeIndicesFromVector( _unique_free_micro_nodes, duplicateNodes.begin( ), duplicateNodes.end( ) );

        if ( error ){

            errorOut result = new errorNode( "setMicroNodeIndexMappings",
                                             "Error in the removal of the non-overlapped duplicate values from the free micro-node vector" );
            result->addNext( error );
            return result;

        }

        //Loop through the ghost nodes to find duplicates
        duplicateNodes.clear( );
        duplicateNodes.reserve( numNodes );
        
        for ( auto node =  _unique_ghost_micro_nodes.begin( );
                   node != _unique_ghost_micro_nodes.end( );
                   node++ ){

            //If the ghost node is found in the free nodes add it to the duplicates
            if ( std::find( _unique_non_overlapped_micro_nodes.begin( ), _unique_non_overlapped_micro_nodes.end( ),  *node )
                 != _unique_non_overlapped_micro_nodes.end( ) ){

                duplicateNodes.push_back( node - _unique_ghost_micro_nodes.begin( ) );

            }

        }

        //Sort the duplicate nodes
        std::sort( duplicateNodes.begin( ), duplicateNodes.end( ) );

        //Remove the duplicate nodes
        error = removeIndicesFromVector( _unique_ghost_micro_nodes, duplicateNodes.begin( ), duplicateNodes.end( ) );

        if ( error ){

            errorOut result = new errorNode( "setMicroNodeIndexMappings",
                                             "Error in the removal of non-overlapped duplicate values from the ghost vector" );
            result->addNext( error );
            return result;

        }

        //Find any duplicated between the free nodes and the ghost nodes
        duplicateNodes.clear( );
        duplicateNodes.reserve( numNodes );
        
        for ( auto node =  _unique_ghost_micro_nodes.begin( );
                   node != _unique_ghost_micro_nodes.end( );
                   node++ ){

            //If the ghost node is found in the free nodes add it to the duplicates
            if ( std::find( _unique_free_micro_nodes.begin( ), _unique_free_micro_nodes.end( ),  *node )
                 != _unique_free_micro_nodes.end( ) ){

                duplicateNodes.push_back( node - _unique_ghost_micro_nodes.begin( ) );

            }

        }

        //Sort the duplicate nodes
        std::sort( duplicateNodes.begin( ), duplicateNodes.end( ) );

        //Remove the duplicate nodes
        error = removeIndicesFromVector( _unique_ghost_micro_nodes, duplicateNodes.begin( ), duplicateNodes.end( ) );

        if ( error ){

            errorOut result = new errorNode( "setMicroNodeIndexMappings",
                                             "Error in the removal of free duplicate values from the ghost vector" );
            result->addNext( error );
            return result;

        }

        //Set the local to global maps
        _global_to_local_micro_node_map.clear();
        _global_to_local_micro_node_map.reserve( _unique_free_micro_nodes.size( ) + _unique_ghost_micro_nodes.size( ) );

        for ( auto node  = _unique_free_micro_nodes.begin( );
                   node != _unique_free_micro_nodes.end( );
                   node++ ){

            _global_to_local_micro_node_map[ *node ] = node - _unique_free_micro_nodes.begin( );

        }

        for ( auto node  = _unique_ghost_micro_nodes.begin( );
                   node != _unique_ghost_micro_nodes.end( );
                   node++ ){

            _global_to_local_micro_node_map[ *node ] = node - _unique_ghost_micro_nodes.begin( )
                                                     + _unique_free_micro_nodes.size( );

        }

        return NULL;
    }

    errorOut inputFileProcessor::setMacroNodeIndexMappings( const unsigned int &increment ){
        /*
         * Set the macro node index mappings for the output matrices.
         * The order that the node appears in the unique node vectors
         * is the index.
         *
         * If a macro-scale node is found in ghost then it cannot be free.
         * We give preference to ghost nodes since the micro-scale is assumed
         * to be more accurate than the macroscale.
         */

        //Get the unique nodes in the free and ghost domains
        errorOut error = getUniqueNodesInDomains( increment, _macroscale, _free_macro_volume_sets, _unique_free_macro_nodes );

        if ( error ){

            errorOut result = new errorNode( "setMacroNodeIndexMappings",
                                             "Error in determining the unique free macroscale nodes" );
            result->addNext( error );
            return result;

        }

        error = getUniqueNodesInDomains( increment, _macroscale, _ghost_macro_volume_sets, _unique_ghost_macro_nodes );

        if ( error ){

            errorOut result = new errorNode( "setMacroNodeIndexMappings",
                                             "Error in determining the unique ghost macroscale nodes" );
            result->addNext( error );
            return result;

        }

        //Remove nodes found in the ghost nodes from the free nodes
        unsigned int numNodes = 0;
        unsigned int n;

        //Approximate the size of the duplicate nodes. At worst, this will be the size of the
        //nodes on the surfaces of the ghost domains. Since there should be many less nodes 
        //in the macroscale than the microscale we will just loop through the existing macro-scale
        //nodesets of the volumes.

        for ( auto domain =  _ghost_macro_volume_sets.begin( );
                   domain != _ghost_macro_volume_sets.end( );
                   domain++ ){

            _macroscale->getNumDomainNodes( increment, *domain, n );
            numNodes += n;

        }

        //Loop through the free nodes to find duplicates
        uIntVector duplicateNodes;
        duplicateNodes.reserve( numNodes );
        
        for ( auto node =  _unique_free_macro_nodes.begin( );
                   node != _unique_free_macro_nodes.end( );
                   node++ ){

            //If the free node is found in the ghost nodes add it to the duplicates
            if ( std::find( _unique_ghost_macro_nodes.begin( ), _unique_ghost_macro_nodes.end( ),  *node )
                 != _unique_ghost_macro_nodes.end( ) ){

                duplicateNodes.push_back( node - _unique_free_macro_nodes.begin( ) );

            }

        }

        //Sort the duplicate nodes
        std::sort( duplicateNodes.begin( ), duplicateNodes.end( ) );

        //Remove the duplicate nodes
        error = removeIndicesFromVector( _unique_free_macro_nodes, duplicateNodes.begin( ), duplicateNodes.end( ) );

        if ( error ){

            errorOut result = new errorNode( "setMacroNodeIndexMappings",
                                             "Error in the removal of the duplicate values from the vector" );
            result->addNext( error );
            return result;

        }

        //Set the global to local maps
        _global_to_local_macro_node_map.clear();
        _global_to_local_macro_node_map.reserve( _unique_free_macro_nodes.size( ) + _unique_ghost_macro_nodes.size( ) );

        for ( auto node  = _unique_free_macro_nodes.begin( );
                   node != _unique_free_macro_nodes.end( );
                   node++ ){

            _global_to_local_macro_node_map[ *node ] = node - _unique_free_macro_nodes.begin( );

        }

        for ( auto node  = _unique_ghost_macro_nodes.begin( );
                   node != _unique_ghost_macro_nodes.end( );
                   node++ ){

            _global_to_local_macro_node_map[ *node ] = node - _unique_ghost_macro_nodes.begin( )
                                                     + _unique_free_macro_nodes.size( );

        }

        return NULL;
    }

    template< typename T, typename Iter >
    errorOut inputFileProcessor::removeIndicesFromVector( std::vector< T > & v, Iter begin, Iter end ){
        /*!
         * Remove the specified indices from the vector
         *
         * Method from https://codereview.stackexchange.com/questions/206686/removing-by-indices-several-elements-from-a-vector?rq=1
         *
         * :param std::vector< T > & v: The vector to be modified
         * :param Iter begin: The beginning iterator of the index vector
         * :param Iter end: The end iterator of the index vector
         */

        if ( !std::is_sorted( begin, end ) ){

            return new errorNode( "removeIndicesFromVector", "The index vector is not sorted" );

        }

        auto rm_iter = begin;
        std::size_t current_index = 0;

        const auto pred = [&]( const T& ){
            
            // any more to remove?
            if ( rm_iter == end ) { return false; }

            // is this one specified?
            if ( *rm_iter == current_index++ ){ return ++rm_iter, true; }

            return false;

        }; 

        v.erase( std::remove_if( v.begin( ), v.end( ), pred ), v.end( ) );

        return NULL;
    }

    errorOut inputFileProcessor::checkVolumeReconstructionInitialization( ){
        /*!
         * Check the initialization of the volume reconstruction
         */

        if ( !_config[ "volume_reconstruction" ] ){

            _config[ "volume_reconstruction" ][ "type" ] = "dual_contouring";

        }

        _volumeReconstructionConfig = _config[ "volume_reconstruction" ];

        return NULL;
    }

    YAML::Node inputFileProcessor::getVolumeReconstructionConfig( ){
        /*!
         * Return the volume reconstruction configuration
         */

        return _volumeReconstructionConfig;
    }

    bool inputFileProcessor::microBodyForceDefined( ){
        /*!
         * Get whether the micro-body force has been defined
         */

        return _microBodyForceFlag;
    }

    bool inputFileProcessor::microAccelerationDefined( ){
        /*!
         * Get whether the micro-acceleration has been defined
         */

        return _microAccelerationFlag;
    }

    const floatVector* inputFileProcessor::getMicroDisplacements( ){
        /*!
         * Get the micro-displacements
         */

        return &_microDisplacements;
    }

    const floatVector* inputFileProcessor::getMacroDisplacements( ){
        /*!
         * Get the macro-displacements
         */

        return &_macroDisplacements;
    }

    const floatVector* inputFileProcessor::getMacroDispDOFVector( ){
        /*!
         * Get the macro-displacement DOF vector
         */

        return &_macroDispDOFVector;
    }

    const floatVector* inputFileProcessor::getMicroNodeReferencePositions( ){
        /*!
         * Get the nodal positions from which the displacements are
         * referenced.
         */

        return &_microNodeReferencePositions;
    }

    const floatVector* inputFileProcessor::getMacroNodeReferencePositions( ){
        /*!
         * Get the nodal positions from which the displacements are
         * referenced.
         */

        return &_macroNodeReferencePositions;
    }

    const uIntVector* inputFileProcessor::getMacroNodeReferenceConnectivity( ){
        /*!
         * Get the nodal positions from which the displacements are
         * referenced.
         */

        return &_macroNodeReferenceConnectivity;
    }

    const uIntVector* inputFileProcessor::getMacroNodeReferenceConnectivityCellIndices( ){
        /*!
         * Get the nodal positions from which the displacements are
         * referenced.
         */

        return &_macroNodeReferenceConnectivityCellIndices;
    }

    const uIntVector* inputFileProcessor::getFreeMacroCellIds( ){
        /*!
         * Get the free macro cell ids
         */

        return &_free_macro_cell_ids;
    }

    const uIntVector* inputFileProcessor::getGhostMacroCellIds( ){
        /*!
         * Get the ghost macro cell ids
         */

        return &_ghost_macro_cell_ids;
    }

    const uIntVector* inputFileProcessor::getFreeMacroCellMicroDomainCounts( ){
        /*!
         * Get the free macro cell ids
         */

        return &_free_macro_cell_micro_domain_counts;
    }

    const uIntVector* inputFileProcessor::getGhostMacroCellMicroDomainCounts( ){
        /*!
         * Get the ghost macro cell ids
         */

        return &_ghost_macro_cell_micro_domain_counts;
    }

    bool inputFileProcessor::computeMicroShapeFunctions( ){
        /*!
         * Return whether the micro-shape functions should be computed
         */

        return _computeMicroShapeFunctions;
    }

    const uIntVector *inputFileProcessor::getNonOverlappedMicroNodeIds( ){
        /*!
         * Get the free micro-node ids
         */

        return &_unique_non_overlapped_micro_nodes;
    }

    const uIntVector *inputFileProcessor::getFreeMicroNodeIds( ){
        /*!
         * Get the free micro-node ids
         */

        return &_unique_free_micro_nodes;
    }

    const uIntVector *inputFileProcessor::getGhostMicroNodeIds( ){
        /*!
         * Get the ghost micro-node ids
         */

        return &_unique_ghost_micro_nodes;
    }

    const stringVector *inputFileProcessor::getFreeMacroDomainNames( ){
        /*!
         * Get the free macro volume sets
         */

        return &_free_macro_volume_sets;
    }

    const stringVector *inputFileProcessor::getGhostMacroDomainNames( ){
        /*!
         * Get the ghost macro volume sets
         */

        return &_ghost_macro_volume_sets;
    }

    const uIntVector *inputFileProcessor::getFreeMacroNodeIds( ){
        /*!
         * Get the free micro-node ids
         */

        return &_unique_free_macro_nodes;
    }

    const uIntVector *inputFileProcessor::getGhostMacroNodeIds( ){
        /*!
         * Get the ghost micro-node ids
         */

        return &_unique_ghost_macro_nodes;
    }

    const DOFMap *inputFileProcessor::getMicroGlobalToLocalDOFMap( ){
        /*!
         * Get the DOF map from the global to local DOF Id for the micro nodes
         */

        return &_global_to_local_micro_node_map;
    }

    const DOFMap *inputFileProcessor::getMacroGlobalToLocalDOFMap( ){
        /*!
         * Get the DOF map from the global to local DOF Id for the macro nodes
         */

        return &_global_to_local_macro_node_map;
    }

    bool inputFileProcessor::useReconstructedMassCenters( ){
        /*!
         * Query whether the reconstructed mass centers should be used or not
         * in the homogenization.
         */

        return _config[ "coupling_initialization" ][ "use_reconstructed_mass_centers" ].as< bool >( );

    }

}
