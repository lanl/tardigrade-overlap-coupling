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

        //Check the coupling initialization
        error = checkCouplingInitialization( );
        if ( error ){
            _error = new errorNode( "initialize", "Error in the coupling initialization configuration" );
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
                errorOut result = new errorNode( "initializeFileInterfaces", "Error in construction of the macroscale interface" );
                result->addNext( _macroscale->_error );
                return result;
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
                errorOut result = new errorNode( "initializeFileInterfaces", "Error in construction of the microscale interface" );
                result->addNext( _microscale->_error );
                return result;
            }

        }
        else{

            return new errorNode( "initializeFileInterfaces",
                                  "There is no 'microscale_definition' in the YAML configuration file" );

        }

        return NULL;
    }

    errorOut inputFileProcessor::setMicroNodeWeights( const unsigned int increment ){
        /*!
         * Compute the weights of the micro-nodes and stores them in a map from the global
         * node id to the weight.
         *
         * :param const unsigned int increment: The increment at which to compute the weights
         */

        //Initialize the weight map
        _microDomainWeights.clear( );
        for ( auto n = _unique_free_micro_nodes.begin( ); n != _unique_free_micro_nodes.end( ); n++ ){
            _microDomainWeights.emplace( *n, 0 );
        }
        for ( auto n = _unique_ghost_micro_nodes.begin( ); n != _unique_ghost_micro_nodes.end( ); n++ ){
            _microDomainWeights.emplace( *n, 0 );
        }

        //Loop through the free micro-volume sets
        uIntVector setNodes;
        errorOut error;
//        for ( auto setName = _free_micro_surface_sets.begin(); setName != _free_micro_surface_sets.end(); setName++ ){
        for ( auto setName = _free_micro_volume_sets.begin(); setName != _free_micro_volume_sets.end(); setName++ ){

            //Get the volume set
            error = _microscale->getSubDomainNodes( increment, *setName, setNodes );

            if ( error ){
                errorOut result = new errorNode( "setMicroNodeWeights",
                                                 "Error in the extraction of the free micro volume set " + *setName );
                result->addNext( error );
                return result;
            }

            //Loop over the nodes
            for ( auto n = setNodes.begin(); n != setNodes.end(); n++ ){
                
                _microDomainWeights[ *n ] += 1;

            }

        }

        //Loop through the ghost micro-volume sets
//        for ( auto setName = _ghost_micro_surface_sets.begin(); setName != _ghost_micro_surface_sets.end(); setName++ ){
        for ( auto setName = _ghost_micro_volume_sets.begin(); setName != _ghost_micro_volume_sets.end(); setName++ ){

            //Get the volume set
            error = _microscale->getSubDomainNodes( increment, *setName, setNodes );

            if ( error ){
                errorOut result = new errorNode( "setMicroNodeWeights",
                                                 "Error in the extraction of the ghost micro volume set " + *setName );
                result->addNext( error );
                return result;
            }

            //Loop over the nodes
            for ( auto n = setNodes.begin(); n != setNodes.end(); n++ ){
                
                _microDomainWeights[ *n ] += 1;

            }

        }

        //Compute the weights
        for ( auto w = _microDomainWeights.begin( ); w != _microDomainWeights.end( ); w++ ){

            _microDomainWeights[ w->first ] = 1. / w->second;

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

        //Set the weights of the micro-nodes
        error = setMicroNodeWeights( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in computation of the micro-node weights" );
            result->addNext( error );
            return result;
        }

        //Set the output index maps
        error = setMicroNodeOutputIndexMappings( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in setting the micro node to output index map" );
            result->addNext( error );
            return result;
        }

        error = setMacroNodeOutputIndexMappings( macroIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in setting the macro node to output index map" );
            result->addNext( error );
            return result;
        }

        //Extract the timestamp of the micro domain at the indicated increment
        error = extractMicroTime( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extraction of the micro timestamp" );
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

        //Extract the micro surface forces
        error = extractMicroSurfaceForces( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the micro surface forces" );
            result->addNext( error );
            return result;
        }

        //Extract the micro external forces
        error = extractMicroExternalForces( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the micro external forces" );
            result->addNext( error );
            return result;
        }

        //Extract the micro velocities
        error = extractMicroVelocities( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the micro velocities" );
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

        if ( _extractPreviousDOFValues ){

            bool tmpFlag;
            uIntType previousMicroIncrement = _config[ "coupling_initialization" ][ "previous_micro_increment" ].as< uIntType >( );

            //Extract the previous time
            error = extractMicroTime( previousMicroIncrement, _previousMicroTime );

            if ( error ){
                errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the previous micro time" );
                result->addNext( error );
                return result;
            }

            //Extract the micro displacements
            error = extractMicroDisplacements( previousMicroIncrement, tmpFlag, _previousMicroDisplacements );
    
            if ( error ){
                errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the previous micro velocities" );
                result->addNext( error );
                return result;
            }
    
            //Extract the micro velocities
            error = extractMicroVelocities( previousMicroIncrement, tmpFlag, _previousMicroVelocities );
    
            if ( error ){
                errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the previous micro velocities" );
                result->addNext( error );
                return result;
            }
    
            //Extract the micro accelerations
            error = extractMicroAccelerations( previousMicroIncrement, tmpFlag, _previousMicroAccelerations );
    
            if ( error ){
                errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the previous micro accelerations" );
                result->addNext( error );
                return result;
            }

        }

        //Extract the micro stresses
        error = extractMicroStresses( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the micro stresses" );
            result->addNext( error );
            return result;
        }

        //Extract the micro internal forces
        error = extractMicroInternalForces( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the micro internal forces" );
            result->addNext( error );
            return result;
        }

        //Extract the micro inertial forces
        error = extractMicroInertialForces( microIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the micro inertial forces" );
            result->addNext( error );
            return result;
        }

        //Extract the timestamp of the macro domain at the indicated increment
        error = extractMacroTime( macroIncrement );
    
        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extraction of the macro timestamp" );
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
    
        //Extract the macro velocities
        error = extractMacroVelocities( macroIncrement );
    
        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the macro velocities" );
            result->addNext( error );
            return result;
        }
    
        //Extract the macro accelerations
        error = extractMacroAccelerations( macroIncrement );
    
        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the macro accelerations" );
            result->addNext( error );
            return result;
        }
    
        if ( _extractPreviousDOFValues ){
    
            bool tmpFlag;
            uIntType previousMacroIncrement = _config[ "coupling_initialization" ][ "previous_macro_increment" ].as< uIntType >( );
    
            //Extract the previous time
            error = extractMacroTime( previousMacroIncrement, _previousMacroTime );
    
            if ( error ){
                errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the previous micro time" );
                result->addNext( error );
                return result;
            }
    
            if ( !vectorTools::fuzzyEquals( _microTime - _previousMicroTime, _macroTime - _previousMacroTime ) && !_isFiltering ){
    
                return new errorNode( "initializeIncrement",
                                      "The change in time between increments for the macro-scale and micro-scale is not consistent" );
    
            }
    
            _Dt = _microTime - _previousMicroTime;
    
            //Extract the macro displacements
            error = extractMacroDispDOFVector( previousMacroIncrement, tmpFlag, _previousMacroDispDOFVector );
        
            if ( error ){
                errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the previous macro displacements" );
                result->addNext( error );
                return result;
            }
        
            //Extract the macro velocities
            error = extractMacroVelocities( previousMacroIncrement, tmpFlag, _previousMacroVelocities );
        
            if ( error ){
                errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the previous macro velocities" );
                result->addNext( error );
                return result;
            }
        
            //Extract the macro accelerations
            error = extractMacroAccelerations( previousMacroIncrement, tmpFlag, _previousMacroAccelerations );
        
            if ( error ){
                errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the previous macro accelerations" );
                result->addNext( error );
                return result;
            }
    
        }
    
        //Extract the macro internal forces
        error = extractMacroInternalForces( macroIncrement );
    
        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the macro internal forces" );
            result->addNext( error );
            return result;
        }
    
        //Extract the macro inertial forces
        error = extractMacroInertialForces( macroIncrement );
    
        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the macro inertial forces" );
            result->addNext( error );
            return result;
        }
    
        //Extract the macro body forces / couples
        error = extractMacroBodyForces( microIncrement );
    
        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the macro body forces and couples" );
            result->addNext( error );
            return result;
        }
    
        //Extract the macro surface forces / couples
        error = extractMacroSurfaceForces( microIncrement );
    
        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the macro surface forces and couples" );
            result->addNext( error );
            return result;
        }
    
    
        //Extract the macro external forces
        error = extractMacroExternalForces( macroIncrement );
    
        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the macro external forces" );
            result->addNext( error );
            return result;
        }

        //Extract the macro lumped-mass matrix
        error = extractMacroLumpedMassMatrix( macroIncrement );

        if ( error ){
            errorOut result = new errorNode( "initializeIncrement", "Error in the extract of the macro lumped-mass matrix" );
            result->addNext( error );
            return result;
        }
    
        if ( _useArlequinMethod ){
    
            error = extractMacroArlequinWeights( macroIncrement );
    
            if ( error ){
    
                errorOut result = new errorNode( "initializeIncrement", "Error in the extraction of the Arlequin macro node weights\n" );
                result->addNext( error );
                return result;
    
            }
    
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

        errorOut error;
        if ( _config[ "free_macroscale_domains" ] ){

            _microDomainSurfaceCount.clear( );
            _macroCellDomainMap.clear( );
            _macroReferenceDensityTypes.clear( );
            _macroReferenceMomentOfInertiaTypes.clear( );
            _macroReferenceDensities.clear( );
            _macroReferenceMomentsOfInertia.clear( );
            error = checkCommonDomainConfiguration( _config[ "free_macroscale_domains" ],
                                                    _free_macro_cell_ids,
                                                    _free_macro_volume_sets,
                                                    _ghost_micro_volume_sets,
                                                    _microDomainSurfaceCount,
                                                    _macroCellDomainMap,
                                                    _freeMacroMassPropertiesRequired,
                                                    _macroReferenceDensityTypes,
                                                    _macroReferenceMomentOfInertiaTypes,
                                                    _macroReferenceDensities,
                                                    _macroReferenceMomentsOfInertia );
            if ( error ){
                errorOut result = new errorNode( "initializeCouplingDomains",
                                                 "Error in input-file check of the free macroscale domains" );
                result->addNext( error );
                return result;
            }

        }

        if ( _config[ "ghost_macroscale_domains" ] ){

            if ( _useArlequinMethod ){
                error = checkCommonDomainConfiguration( _config[ "ghost_macroscale_domains" ],
                                                        _ghost_macro_cell_ids,
                                                        _ghost_macro_volume_sets,
                                                        _free_micro_volume_sets,
                                                        _microDomainSurfaceCount,
                                                        _macroCellDomainMap,
                                                        _ghostMacroMassPropertiesRequired,
                                                        _macroReferenceDensityTypes,
                                                        _macroReferenceMomentOfInertiaTypes,
                                                        _macroReferenceDensities,
                                                        _macroReferenceMomentsOfInertia );
            }
            else{
                error = checkCommonDomainConfiguration( _config[ "ghost_macroscale_domains" ],
                                                        _ghost_macro_cell_ids,
                                                        _ghost_macro_volume_sets,
                                                        _free_micro_volume_sets,
                                                        _microDomainSurfaceCount,
                                                        _macroCellDomainMap );
            }
            if ( error ){
                errorOut result = new errorNode( "initializeCouplingDomains",
                                                 "Error in input-file check of the ghost macroscale domains" );
                result->addNext( error );
                return result;
            }

        }

        //Make sure that no volume nodeset that appears in the ghost micro-scale also appears in the free micro-scale
        for ( auto nodeset = _ghost_micro_volume_sets.begin( ); nodeset != _ghost_micro_volume_sets.end( ); nodeset++ ){

            if ( std::find( _free_micro_volume_sets.begin( ), _free_micro_volume_sets.end( ), *nodeset ) != _free_micro_volume_sets.end( ) ){
                return new errorNode( "initializeCouplingDomains",
                                      *nodeset + " appears in the ghost and free micro-volume nodeset definitions" );
            }

        }

        //Initialize the internal ID number of the micro domains
        _microDomainIDMap.clear( );
        _microDomainIDMap.reserve( _free_micro_volume_sets.size( ) + _ghost_micro_volume_sets.size( ) );

        uIntType index = 0;
        for ( auto nodeset = _free_micro_volume_sets.begin( ); nodeset != _free_micro_volume_sets.end( ); nodeset++, index++ ){

            _microDomainIDMap.emplace( *nodeset, index );

        }
        for ( auto nodeset = _ghost_micro_volume_sets.begin( ); nodeset != _ghost_micro_volume_sets.end( ); nodeset++, index++ ){

            _microDomainIDMap.emplace( *nodeset, index );

        }

        return NULL;

    }

    errorOut inputFileProcessor::checkCommonDomainConfiguration( YAML::Node domainConfig,
                                                                 uIntVector &macroCellIds,
//                                                                 uIntVector &macroCellMicroDomainCounts,
                                                                 stringVector &macroVolumeNodesets,
                                                                 stringVector &microVolumeNodesets,
                                                                 std::unordered_map< std::string, uIntType > &microSurfaceDomainCount,
                                                                 std::unordered_map< uIntType, stringVector > &macroCellDomainMap ){
        /*!
         * Extract common values in the configuration of a domain
         *
         * :param YAML::Node domainConfig: The configuration of a particular domain.
         * :param uIntVector &macroCellIds: The macro-cell Ids corresponding to the domain
         * :param stringVector &macroVolumeNodesets: The nodeset names for the nodes in the
         *     macro domains
         * :param stringVector &microVolumeNodesets: The nodeset names for the nodes in the
         *     micro domains
         * :param std::unordered_map< std::string, uIntType > &microSurfaceDomainCount: The approximate number of subdomains to split the
         *     micro domain's surface into
         */

        bool massPropertyDefinitionRequired = false;
        std::unordered_map< unsigned int, floatVector > density;
        std::unordered_map< unsigned int, floatVector > microInertia;
        std::unordered_map< unsigned int, std::string > densityTypes;
        std::unordered_map< unsigned int, std::string > microInertiaTypes;

        return checkCommonDomainConfiguration( domainConfig, macroCellIds,// macroCellMicroDomainCounts,
                                               macroVolumeNodesets, microVolumeNodesets,
                                               microSurfaceDomainCount, macroCellDomainMap, massPropertyDefinitionRequired,
                                               densityTypes, microInertiaTypes, density, microInertia );

    }

    errorOut inputFileProcessor::checkCommonDomainConfiguration( YAML::Node domainConfig,
                                                                 uIntVector &macroCellIds,
//                                                                 uIntVector &macroCellMicroDomainCounts,
                                                                 stringVector &macroVolumeNodesets,
                                                                 stringVector &microVolumeNodesets,
                                                                 std::unordered_map< std::string, uIntType > &microSurfaceDomainCount,
                                                                 std::unordered_map< unsigned int, stringVector > &macroCellToDomainMap,
                                                                 const bool &massPropertyDefinitionRequired,
                                                                 std::unordered_map< unsigned int, std::string > &densityTypes,
                                                                 std::unordered_map< unsigned int, std::string > &microInertiaTypes,
                                                                 std::unordered_map< unsigned int, floatVector > &density,
                                                                 std::unordered_map< unsigned int, floatVector > &microInertia ){
        /*!
         * Extract common values in the configuration of a domain
         *
         * :param YAML::Node domainConfig: The configuration of a particular domain.
         * :param uIntVector &macroCellIds: The macro-cell Ids corresponding to the domain
         * :param stringVector &macroVolumeNodesets: The nodeset names for the nodes in the
         *     macro domains
         * :param stringVector &microVolumeNodesets: The nodeset names for the nodes in the
         *     micro domains
         * :param std::unordered_map< std::string, uIntType > &microSurfaceDomainCount: The approximate number of subdomains to split the
         *     micro domain's surface into
         * :param const bool &massPropertyDefinitionRequired: Flag which indicates if the mass properties
         *     must be defined for the element
         * :param std::unordered_map< unsigned int, std::string > &densityTypes: The types of the densities
         *     Currently only constant is accepted
         * :param std::unordered_map< unsigned int, std::string > &microInertiaTypes: The types of the micro inertias
         *     Currently only constant is accepted
         * :param std::unordered_map< unsigned int, floatVector > &density: The density of the element. Either a single value if it is constant
         *     or multiple defined at each gauss point.
         * :param std::unordered_map< unsigned int, floatVector > &microInertia: The micro inertia of the element. Either a single symmetric matrix if it is 
         *     constant or multiple if defined for each gauss point.
         * :param std::unordered_map< unsigned int, &macroCellToDomainMap: The map from the macro cell to all of the micro 
         *     domains it contains 
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

            try {

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
                for ( auto nodeset  = ( *domain )[ "micro_nodesets" ].begin( );
                           nodeset != ( *domain )[ "micro_nodesets" ].end( );
                           nodeset++ ){
    
                    if ( !( *nodeset )[ "name" ] ){
    
                        return new errorNode( "checkCommonDomainConfiguration",
                                              "The keyword 'name' is not defined in micro-nodeset entry " + std::to_string( indx2 ) +
                                              " of domain entry " + std::to_string( indx ) + " is not defined" );
    
                    }
    
                    if ( !( *nodeset )[ "number_of_surface_microdomains" ] ){
                           
                           ( *nodeset )[ "number_of_surface_microdomains" ] = _defaultNumberOfMicroDomainSurfaceRegions;
    
                    }
                    else if ( !( *nodeset )[ "number_of_surface_microdomains" ].IsScalar( ) ){
    
                        return new errorNode( "checkCommonDomainConfiguration",
                                              "Micro-nodeset 'number_of_surface_microdomains' in entry " + std::to_string( indx2 ) +
                                              " of domain entry " + std::to_string( indx ) + " must be a scalar integer" );
    
                    }
    
                    if ( !( *nodeset )[ "name" ].IsScalar( ) ){
    
                        return new errorNode( "checkCommonDomainConfiguration",
                                              "Micro-nodeset entry " + std::to_string( indx2 ) + " of domain entry " + 
                                              std::to_string( indx ) + " is not a Scalar" );
    
                    }
    
                    nVolumeNodesets++;
                    indx2++;
    
                }
                
                if ( massPropertyDefinitionRequired ){
        
                    //Extract the reference Density
                    try{

                        if ( !( *domain )[ "reference_density" ] ){
        
                            return new errorNode( "checkCommonDomainConfiguration",
                                                  "The reference density is required for the macro-domain in entry " + std::to_string( indx ) +
                                                  " but is not defined." );
        
                        }
                        if ( !( *domain )[ "reference_density" ][ "type" ] ){
        
                            std::string outstr  = "The type of the reference density must be defined.";
                                        outstr += "  Acceptable types are:\n";
                                        outstr += "    constant";
                            return new errorNode( "checkCommonDomainConfiguration", outstr );
        
                        }
    
                        if ( ( *domain )[ "reference_density" ][ "type" ].as< std::string >( ).compare( "constant" ) == 0 ){
        
                            if ( !( *domain )[ "reference_density" ][ "value" ] ){
        
                                std::string  outstr = "The value of the reference density for macro domain " + std::to_string( indx ) + " is not defined";
                                            outstr += "The format is:\n  value: floating_point_value";
                                return new errorNode( "checkCommonDomainConfiguration", outstr );
        
                            }
                            else{
        
                                densityTypes.emplace( ( *domain )[ "macro_cell" ].as< unsigned int >( ),
                                                      ( *domain )[ "reference_density" ][ "type" ].as< std::string >( ) );
                                
                                floatVector tmp = { ( *domain )[ "reference_density" ][ "value" ].as< floatType >( ) };
                                density.emplace( ( *domain )[ "macro_cell" ].as< unsigned int >( ), tmp );
        
                            }
        
                        }
                        else{
        
                            std::string outstr  = "The reference density type for macro-domain " + std::to_string( indx ) + " is not recognized.\n";
                                        outstr += "  type: " + ( *domain )[ "reference_density" ][ "type" ].as< std::string >( );
        
                            return new errorNode( "checkCommonDomainConfiguration", outstr );
        
                        }

                    }

                    catch ( std::exception &e ){

                        std::string outstr  = "Unexpected error encountered in the reference density definition of macro-domain " + std::to_string( indx ) + ".\n";
                                    outstr += "This is likely due to a problem in the YAML configuration file.\n";
                                    outstr += "The original error message was:\n";
                                    outstr += e.what( );

                        return new errorNode( "checkCommonDomainConfiguration", outstr );

                    }
   
                    try{

                        if ( !( *domain )[ "reference_moment_of_inertia" ] ){
        
                            return new errorNode( "checkCommonDomainConfiguration",
                                                  "The reference moment of inertia is required for the macro-domain in entry " +
                                                  std::to_string( indx ) + " but is not defined." );
        
                        }

                        if ( !( *domain )[ "reference_moment_of_inertia" ][ "type" ] ){

                            return new errorNode( "checkCommonDomainConfiguration",
                                                  "The reference moment of inertia type is required for the macro-domain in entry " +
                                                  std::to_string( indx ) + " but is not defined." );

                        }
                        if ( ( *domain )[ "reference_moment_of_inertia" ][ "type" ].as< std::string >( ).compare( "constant" ) == 0 ){
        
                            if ( !( *domain )[ "reference_moment_of_inertia" ][ "value" ] ){
        
                                std::string  outstr = "The values of the reference moment of inertia for macro domain " + std::to_string( indx ) + " are not defined";
                                            outstr += "The format is:\n  value: [ I11, I12, I13, I22, I23, I33 ]";
                                return new errorNode( "checkCommonDomainConfiguration", outstr );
        
                            }
                            else if ( !( *domain )[ "reference_moment_of_inertia" ][ "value" ].IsSequence( ) ){
        
                                std::string  outstr = "The values of the reference moment of inertia for macro domain " + std::to_string( indx ) + " are not defined as a sequence";
                                            outstr += "The format is:\n  value: [ I11, I12, I13, I22, I23, I33 ]";
                                return new errorNode( "checkCommonDomainConfiguration", outstr );
        
                            }
                            else{

                                uIntType vindex = 0;
                                uIntType ncomponents =  ( *domain )[ "reference_moment_of_inertia" ][ "value" ].size( );
                                if ( ncomponents != ( _dim * ( _dim + 1 ) / 2 ) ){

                                    return new errorNode( "checkCommonDomainConfiguration",
                                                          "Six terms are required for the definition of a constant reference micro moment of inertia of macro-domain " +
                                                          std::to_string( indx ) + " and " + std::to_string( ncomponents ) );

                                }
                                floatVector tmp;
                                for ( auto v  = ( *domain )[ "reference_moment_of_inertia" ][ "value" ].begin( );
                                           v != ( *domain )[ "reference_moment_of_inertia" ][ "value" ].end( );
                                           v++, vindex++ ){
                
                                    if ( v->IsScalar( ) ){
                
                                        tmp.push_back( v->as< floatType >( ) );
                
                                    }
                                    else{
                
                                        return new errorNode( "checkCommonDomainConfiguration",
                                                              "The micro-inertia entry " + std::to_string( vindex ) +
                                                              " is not a scalar value" );
                
                                    }
                
                                }

                                floatVector domainMicroInertia = { tmp[ 0 ], tmp[ 1 ], tmp[ 2 ],
                                                                   tmp[ 1 ], tmp[ 3 ], tmp[ 4 ],
                                                                   tmp[ 2 ], tmp[ 4 ], tmp[ 5 ] }; //Assuming 3D
                                microInertia.emplace( ( *domain )[ "macro_cell" ].as< unsigned int >( ),
                                                      domainMicroInertia ); 
                                microInertiaTypes.emplace( ( *domain )[ "macro_cell" ].as< unsigned int >( ),
                                                           ( *domain )[ "reference_moment_of_inertia" ][ "type" ].as< std::string >( ) );
        
                            }
        
                        }
                        else{
        
                            std::string outstr  = "The reference density type for macro-domain " + std::to_string( indx ) + " is not recognized.\n";
                                        outstr += "  type: " + ( *domain )[ "reference_density" ][ "type" ].as< std::string >( );
        
                            return new errorNode( "checkCommonDomainConfiguration", outstr );
        
                        }
        
//                        for ( auto v  = ( *domain )[ "reference_moment_of_inertia" ].begin( );
//                                   v != ( *domain )[ "reference_moment_of_inertia" ].end( );
//                                   v++ ){
//        
//                            if ( v->IsScalar( ) ){
//        
//                                microInertia.push_back( v->as< floatType >( ) );
//        
//                            }
//                            else{
//        
//                                return new errorNode( "checkCommonDomainConfiguration",
//                                                      "The micro-inertia must be constant over the free micromorphic element" );
//        
//                            }
//        
//                        }

                    }

                    catch ( std::exception &e ){

                        std::string outstr  = "Unexpected error encountered in the reference moment of inertia definition of macro-domain " + std::to_string( indx ) + ".\n";
                                    outstr += "This is likely due to a problem in the YAML configuration file.\n";
                                    outstr += "The original error message was:\n";
                                    outstr += e.what( );

                        return new errorNode( "checkCommonDomainConfiguration", outstr );

                    }
        
                }

            }
            catch( std::exception &e ){

                std::string outstr  = "Unexpected error encountered when processing macro-domain " + std::to_string( indx ) + ".\nLikely a YAML formatting error.\n";
                            outstr += "The original error message is:\n";
                            outstr += e.what( );
                return new errorNode( "checkCommonDomainConfiguration", outstr );

            }

            indx++;

        }

        //Extract the volume nodesets
        macroCellIds.reserve( domainConfig.size( ) );
//        macroCellMicroDomainCounts.reserve( domainConfig.size( ) );
        macroVolumeNodesets.reserve( domainConfig.size( ) );
        microVolumeNodesets.reserve( nVolumeNodesets );
        for ( auto domain = domainConfig.begin( ); domain != domainConfig.end( ); domain++ ){

            macroCellIds.push_back( ( *domain )[ "macro_cell" ].as< unsigned int >( ) );
            macroVolumeNodesets.push_back( ( *domain )[ "macro_nodeset" ].as< std::string >( ) );
//            macroCellMicroDomainCounts.push_back( ( *domain )[ "micro_nodesets" ].size( ) );

            for ( auto nodeset = ( *domain )[ "micro_nodesets" ].begin( ); nodeset != ( *domain )[ "micro_nodesets" ].end( ); nodeset++ ){

                std::string nodesetName = ( *nodeset )[ "name" ].as< std::string >( );
                uIntType numberOfSurfaceMicroDomains = ( *nodeset )[ "number_of_surface_microdomains" ].as< uIntType >( );

                if ( std::find( microVolumeNodesets.begin( ), microVolumeNodesets.end( ), nodesetName ) != microVolumeNodesets.end( ) ){

                    return new errorNode( "checkCommonDomainConfiguration",
                                          nodeset->as< std::string >( ) + " appears more than once in the coupling definition" );

                }

                microVolumeNodesets.push_back( nodesetName );
                microSurfaceDomainCount.emplace( nodesetName, numberOfSurfaceMicroDomains );

            }

            if ( macroCellToDomainMap.find( macroCellIds.back( ) ) == macroCellToDomainMap.end( ) ){
    
                macroCellToDomainMap.emplace( macroCellIds.back( ), microVolumeNodesets );
    
            }
            else{
                return new errorNode( "checkCommonDomainConfiguration",
                                      "Macro cell " + std::to_string( macroCellIds.back( ) ) + " appears more than once in the coupling definition" );
            }

        }
            

        return NULL;

    }

//    errorOut inputFileProcessor::checkCommonVolumeToSurfaceMapping( const stringVector &microVolumeNodesets, 
//                                                                    stringVector &microSurfaceNodesets ){
//        /*!
//         * Check the common volume to surface node mapping
//         *
//         * :param const stringVector &microVolumeNodesets: The nodesets which represent the volume of the micro-domains
//         * :param const stringVector &microSurfaceNodesets: The nodesets which represent the surfaces of the micro-domains
//         */
//
//        if ( !_config[ "microscale_definition" ][ "volume_surface_node_pairs" ] ){
//            
//            return new errorNode( "checkCommonVolumeToSurfaceMapping",
//                                  "volume_surface_node_pairs must be defined as a keyname attribute of 'microscale_definition'" );
//
//        }
//
//        if ( ( !_config[ "microscale_definition" ][ "volume_surface_node_pairs" ].IsSequence( ) ) && ( microVolumeNodesets.size( ) > 0 ) ){
//            return new errorNode( "checkCommonVolumeToSurfaceMapping",
//                                  "volume_surface_node_pairs must be a sequence" );
//
//        }
//
//        //Loop through the micro-volume nodesets
//        microSurfaceNodesets = stringVector( microVolumeNodesets.size( ) );
//        unsigned int indx = 0;
//        unsigned int p;
//        bool volumeFound = false;
//        for ( auto volumeName = microVolumeNodesets.begin( ); volumeName != microVolumeNodesets.end( ); volumeName++ ){
//
//            p = 1;
//            volumeFound = false;
//            for ( auto pair = _config[ "microscale_definition" ][ "volume_surface_node_pairs" ].begin( );
//                       pair != _config[ "microscale_definition" ][ "volume_surface_node_pairs" ].end( );
//                  pair++ ){
//
//                if ( !( *pair )[ "volume" ] ){
//
//                    return new errorNode( "checkCommonVolumeToSurfaceMapping",
//                                          "'volume' is not a key of 'volume_surface_node_pairs' entry " + std::to_string( p ) );
//
//                }
//                if ( !( *pair )[ "surface" ] ){
//
//                    return new errorNode( "checkCommonVolumeToSurfaceMapping",
//                                          "'surface' is not a key of 'volume_surface_node_pairs' entry " + std::to_string( p ) );
//
//                }
//
//                if ( ( *pair )[ "volume" ].as< std::string >( ).compare( *volumeName ) == 0 ){
//
//                    microSurfaceNodesets[ indx ] = ( *pair )[ "surface" ].as< std::string >( );
//                    volumeFound = true;
//                    break;
//
//                }
//
//                p++;
//
//            }
//
//            if ( !volumeFound ){
//
//                return new errorNode( "checkCommonVolumeToSurfaceMapping"
//                                      "'" + *volumeName + "' not found in 'volume_surface_node_pairs'" );
//
//            }
//
//            indx++;
//
//        }
//
//        return NULL;
//
//    } 

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

        //Initialize the size of the micro densities map
        _microDensities.clear( );
        _microDensities.reserve( _microGlobalNodeIDOutputIndex.size( ) );

        //Get the values of the micro densities from the output file
        floatVector values;
        errorOut error = _microscale->getSolutionData( increment, 
                                                       _config[ "microscale_definition" ][ "density_variable_name" ].as< std::string >( ),
                                                       "Node", values );

        if ( error ){

            errorOut result = new errorNode( "extractMicroNodeDensities", "Error in extraction of the micro densities" );
            result->addNext( error );
            return result;

        }

        for ( auto n = _microGlobalNodeIDOutputIndex.begin( ); n != _microGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMicroNodeDensities", "The density vector is too short for the required index" );

            }

            _microDensities.emplace( n->first, values[ n->second ] );

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroBodyForces( const unsigned int &increment ){
        /*!
         * Extract the node micro-body forces at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        stringVector variableKeys =
            {
                "F1", "F2", "F3",
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "body_force_variable_names";
        YAML::Node configuration = _config[ "microscale_definition" ][ configurationName.c_str( ) ];

        floatVector values;
        errorOut error = inputFileProcessor::extractDataFileProperties( _microscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, _microBodyForceFlag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMicroBodyForces",
                                             "Error in the extraction of the micro body forces" );
            result->addNext( error );
            return result;

        }

        if ( !_microBodyForceFlag ){

            _microBodyForces.clear( );
            return NULL;

        }

        floatType _sign = _config[ "coupling_initialization" ][ "micro_body_force_sign" ].as< floatType >( );

        _microBodyForces.clear( );
        _microBodyForces.reserve( _microGlobalNodeIDOutputIndex.size( ) );
        for ( auto n = _microGlobalNodeIDOutputIndex.begin( ); n != _microGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMicroBodyForces", "The micro body force vector is too short for the required index" );

            }

            floatVector bf( values.begin( ) + variableKeys.size( ) * n->second,
                            values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) );
            bf *= _sign;

            _microBodyForces.emplace( n->first, bf );

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroSurfaceForces( const unsigned int &increment ){
        /*!
         * Extract the node micro-surface forces at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        stringVector variableKeys =
            {
                "F1", "F2", "F3",
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "surface_force_variable_names";
        YAML::Node configuration = _config[ "microscale_definition" ][ configurationName.c_str( ) ];

        floatVector values;
        errorOut error = inputFileProcessor::extractDataFileProperties( _microscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, _microSurfaceForceFlag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMicroSurfaceForces",
                                             "Error in the extraction of the micro surface forces" );
            result->addNext( error );
            return result;

        }

        if ( !_microSurfaceForceFlag ){

            _microSurfaceForces.clear( );
            return NULL;

        }

        floatType _sign = _config[ "coupling_initialization" ][ "micro_surface_force_sign" ].as< floatType >( );

        _microSurfaceForces.clear( );
        _microSurfaceForces.reserve( _microGlobalNodeIDOutputIndex.size( ) );
        for ( auto n = _microGlobalNodeIDOutputIndex.begin( ); n != _microGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMicroSurfaceForces",
                                      "The micro surface forces vector is too short for the required index" );

            }

            floatVector st( values.begin( ) + variableKeys.size( ) * n->second,
                            values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) );
            st *= _sign;

            _microSurfaceForces.emplace( n->first, st );

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroExternalForces( const unsigned int &increment ){
        /*!
         * Extract the node micro external forces at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        stringVector variableKeys =
            {
                "F1", "F2", "F3",
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "external_force_variable_names";
        YAML::Node configuration = _config[ "microscale_definition" ][ configurationName.c_str( ) ];

        floatVector values;
        errorOut error = inputFileProcessor::extractDataFileProperties( _microscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, _microExternalForceFlag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMicroSurfaceForces",
                                             "Error in the extraction of the macro velocities" );
            result->addNext( error );
            return result;

        }

        if ( !_microExternalForceFlag ){

            if( _microSurfaceForceFlag && _microBodyForceFlag ){

                _microExternalForces.clear( );
                _microExternalForces.reserve( _microGlobalNodeIDOutputIndex.size( ) );
                for ( auto n = _microGlobalNodeIDOutputIndex.begin( ); n != _microGlobalNodeIDOutputIndex.end( ); n++ ){

                    _microExternalForces.emplace( n->first, _microSurfaceForces[ n->first ] + _microBodyForces[ n->first ] );

                }
                _microExternalForceFlag = true;

            }
            else if ( _microSurfaceForceFlag ){

                _microExternalForces = _microSurfaceForces;
                _microExternalForceFlag = true;

            }
            else if ( _microBodyForceFlag ){

                _microExternalForces = _microBodyForces;
                _microExternalForceFlag = true;

            }
            else{
                _microExternalForces.clear( );
            }

        }
        else{

            _microExternalForces.clear( );
            _microExternalForces.reserve( _microGlobalNodeIDOutputIndex.size( ) );
            floatType _sign = _config[ "coupling_initialization" ][ "micro_external_force_sign" ].as< floatType >( );

            for ( auto n = _microGlobalNodeIDOutputIndex.begin( ); n != _microGlobalNodeIDOutputIndex.end( ); n++ ){
    
                if ( n->second >= values.size( ) ){
    
                    return new errorNode( "extractMicroExternalForces",
                                          "The micro external force vector is too short for the required index" );
    
                }
    
                floatVector ef( values.begin( ) + variableKeys.size( ) * n->second,
                                values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) );
                ef *= _sign;
    
                _microExternalForces.emplace( n->first, ef );
    
            }

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroAccelerations( const unsigned int &increment ){
        /*!
         * Extract the node micro-accelerations at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        errorOut error = extractMicroAccelerations( increment, _microAccelerationFlag, _microAccelerations );

        if ( error ){

            errorOut result = new errorNode( "extractMicroAccelerations",
                                             "Error in the extraction of the micro accelerations" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroAccelerations( const unsigned int &increment, bool &flag,
                                                            std::unordered_map< uIntType, floatVector > &microAccelerations ){
        /*!
         * Extract the node micro-accelerations at the indicated increment
         *
         * :param const unsigned int &increment: The increment at which to make the extraction
         * :param floatVector &flag: The flag to indicate if the accelerations were defined in the file
         * :param std::unordered_map< uIntType, floatVector > &microAccelerations: The map to store the micro-accelerations in
         */

        stringVector variableKeys =
            {
                "a1", "a2", "a3",
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "acceleration_variable_names";
        YAML::Node configuration = _config[ "microscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        errorOut error = inputFileProcessor::extractDataFileProperties( _microscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, flag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMicroAccelerations",
                                             "Error in the extraction of the micro accelerations" );
            result->addNext( error );
            return result;

        }

        if ( !flag ){

            microAccelerations.clear( );
            return NULL;

        }

        microAccelerations.clear( );
        microAccelerations.reserve( _microGlobalNodeIDOutputIndex.size( ) );
        for ( auto n = _microGlobalNodeIDOutputIndex.begin( ); n != _microGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMicroVelocities",
                                      "The micro velocity vector is too short for the required index" );

            }

            microAccelerations.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                               values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroVelocities( const unsigned int &increment ){
        /*!
         * Extract the node micro-velocities at the indicated increment
         *
         * TODO: Move this and other vector / tensor extraction routines to a common function
         *
         * :param const unsigned int &increment: The current increment
         */

        errorOut error = extractMicroVelocities( increment, _microVelocityFlag, _microVelocities );

        if ( error ){

            errorOut result = new errorNode( "extractMicroVelocities",
                                             "Error in the extraction of the micro velocities" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroVelocities( const unsigned int &increment, bool &flag,
                                                         std::unordered_map< uIntType, floatVector > &microVelocities ){
        /*!
         * Extract the node micro-velocities at the indicated increment
         *
         * TODO: Move this and other vector / tensor extraction routines to a common function
         *
         * :param const unsigned int &increment: The current increment
         * :param bool &flag: The flag to indicate if the values are defined in the input file
         * :param std::unordered_map< uIntType, floatVector > &microVelocities: The map to store the velocities in
         */

        stringVector variableKeys =
            {
                "v1", "v2", "v3",
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "velocity_variable_names";
        YAML::Node configuration = _config[ "microscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        errorOut error = inputFileProcessor::extractDataFileProperties( _microscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, flag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMicroVelocities",
                                             "Error in the extraction of the micro velocities" );
            result->addNext( error );
            return result;

        }

        if ( !flag ){

            microVelocities.clear( );
            return NULL;

        }

        microVelocities.clear( );
        microVelocities.reserve( _microGlobalNodeIDOutputIndex.size( ) );
        for ( auto n = _microGlobalNodeIDOutputIndex.begin( ); n != _microGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMicroVelocities",
                                      "The micro velocity vector is too short for the required index" );

            }

            microVelocities.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                            values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroVelocities( const unsigned int &increment ){
        /*!
         * Extract the node macro-velocities at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        errorOut error = extractMacroVelocities( increment, _macroVelocityFlag, _macroVelocities );

        if ( error ){

            errorOut result = new errorNode( "extractMacroVelocities",
                                             "Error in the extraction of the macro velocities" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroVelocities( const unsigned int &increment, bool &flag,
                                                         std::unordered_map< uIntType, floatVector > &macroVelocities ){
        /*!
         * Extract the node macro-velocities at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         * :param bool &flag: The flag indicating if the values are defined in the file
         * :param std::unordered_map< uIntType, floatVector > &macroVelocities: The map to store the values in
         */

        stringVector variableKeys =
            {
                "v1", "v2", "v3",
                "phiDot11", "phiDot12", "phiDot13",
                "phiDot21", "phiDot22", "phiDot23",
                "phiDot31", "phiDot32", "phiDot33",
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "velocity_variable_names";
        YAML::Node configuration = _config[ "macroscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        errorOut error = inputFileProcessor::extractDataFileProperties( _macroscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, flag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMacroVelocities",
                                             "Error in the extraction of the macro velocities" );
            result->addNext( error );
            return result;

        }

        if ( !flag ){

            _macroVelocities.clear( );
            return NULL;

        }
        else{

            macroVelocities.clear( );
            macroVelocities.reserve( _macroGlobalNodeIDOutputIndex.size( ) );
            for ( auto n = _macroGlobalNodeIDOutputIndex.begin( ); n != _macroGlobalNodeIDOutputIndex.end( ); n++ ){

                if ( n->second > values.size( ) ){

                    return new errorNode( "extractMacroVelocities",
                                          "The values vector is too small for the index required by micro node " +
                                          std::to_string( n->first ) );

                }

                macroVelocities.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                                values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );

            }

        }

        return NULL;

    }


    errorOut inputFileProcessor::extractMacroAccelerations( const unsigned int &increment ){
        /*!
         * Extract the node macro-accelerations at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        errorOut error = extractMacroAccelerations( increment, _macroAccelerationFlag, _macroAccelerations );

        if ( error ){

            errorOut result = new errorNode( "extractMacroAccelerations",
                                             "Error in the extraction of the macro accelerations" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroAccelerations( const unsigned int &increment, bool &flag,
                                                            std::unordered_map< uIntType, floatVector > &macroAccelerations ){
        /*!
         * Extract the node macro-accelerations at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         * :param bool &flag: The flag indicating if the values are defined in the file
         * :param std::unordered_map< uIntType, floatVector > &macroAccelerations: The vector to store the values in
         */

        stringVector variableKeys =
            {
                "a1", "a2", "a3",
                "phiDotDot11", "phiDotDot12", "phiDotDot13",
                "phiDotDot21", "phiDotDot22", "phiDotDot23",
                "phiDotDot31", "phiDotDot32", "phiDotDot33",
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "acceleration_variable_names";
        YAML::Node configuration = _config[ "macroscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        errorOut error = inputFileProcessor::extractDataFileProperties( _macroscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, flag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMacroAccelerations",
                                             "Error in the extraction of the macro accelerations" );
            result->addNext( error );
            return result;

        }

        if ( !flag ){

            macroAccelerations.clear( );
            return NULL;

        }

        macroAccelerations.clear( );
        macroAccelerations.reserve( _macroGlobalNodeIDOutputIndex.size( ) );
        for ( auto n = _macroGlobalNodeIDOutputIndex.begin( ); n != _macroGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMacroAccelerations",
                                      "The values vector is too small for the index required by micro-node " +
                                      std::to_string( n->first ) );

            }

            macroAccelerations.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                               values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroInternalForces( const unsigned int &increment ){
        /*!
         * Extract the node macro internal forces at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        stringVector variableKeys =
            {
                "F1", "F2", "F3",
                "C11", "C12", "C13",
                "C21", "C22", "C23",
                "C31", "C32", "C33"
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "internal_force_variable_names";
        YAML::Node configuration = _config[ "macroscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        errorOut error = inputFileProcessor::extractDataFileProperties( _macroscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, _macroInternalForceFlag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMacroInternalForces",
                                             "Error in the extraction of the macro internal forces" );
            result->addNext( error );
            return result;

        }

        if ( !_macroInternalForceFlag ){

            _macroInternalForces.clear( );
            return NULL;

        }

        _macroInternalForces.clear( );
        _macroInternalForces.reserve( _macroGlobalNodeIDOutputIndex.size( ) );
        floatType _sign = _config[ "coupling_initialization" ][ "macro_internal_force_sign" ].as< floatType >( );
        for( auto n = _macroGlobalNodeIDOutputIndex.begin( ); n != _macroGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMacroInternalForces",
                                      "The index required by macro node " + std::to_string( n->first ) + " is too large for the values vector" );

            }

            _macroInternalForces.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                                 values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );
            _macroInternalForces[ n->first ] *= _sign;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroBodyForces( const unsigned int &increment ){
        /*!
         * Extract the node macro body forces and couples at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        stringVector variableKeys =
            {
                "F1", "F2", "F3",
                "C11", "C12", "C13",
                "C21", "C22", "C23",
                "C31", "C32", "C33"
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "body_force_variable_names";
        YAML::Node configuration = _config[ "macroscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        errorOut error = inputFileProcessor::extractDataFileProperties( _macroscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, _macroBodyForceFlag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMacroBodyForces",
                                             "Error in the extraction of the macro body forces" );
            result->addNext( error );
            return result;

        }

        if ( !_macroBodyForceFlag ){

            _macroBodyForces.clear( );
            return NULL;

        }

        _macroBodyForces.clear( );
        _macroBodyForces.reserve( _macroGlobalNodeIDOutputIndex.size( ) );
        floatType _sign = _config[ "coupling_initialization" ][ "macro_body_force_sign" ].as< floatType >( );
        for( auto n = _macroGlobalNodeIDOutputIndex.begin( ); n != _macroGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMacroBodyForces",
                                      "The index required by macro node " + std::to_string( n->first ) + " is too large for the values vector" );

            }

            _macroBodyForces.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                             values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );
            _macroBodyForces[ n->first ] *= _sign;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroSurfaceForces( const unsigned int &increment ){
        /*!
         * Extract the node macro surface forces and couples at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        stringVector variableKeys =
            {
                "F1", "F2", "F3",
                "C11", "C12", "C13",
                "C21", "C22", "C23",
                "C31", "C32", "C33"
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "surface_force_variable_names";
        YAML::Node configuration = _config[ "macroscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        errorOut error = inputFileProcessor::extractDataFileProperties( _macroscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, _macroSurfaceForceFlag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMacroSurfaceForces",
                                             "Error in the extraction of the macro surface forces" );
            result->addNext( error );
            return result;

        }

        if ( !_macroSurfaceForceFlag ){

            _macroSurfaceForces.clear( );
            return NULL;

        }

        _macroSurfaceForces.clear( );
        _macroSurfaceForces.reserve( _macroGlobalNodeIDOutputIndex.size( ) );
        floatType _sign = _config[ "coupling_initialization" ][ "macro_surface_force_sign" ].as< floatType >( );
        for( auto n = _macroGlobalNodeIDOutputIndex.begin( ); n != _macroGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMacroSurfaceForces",
                                      "The index required by macro node " + std::to_string( n->first ) + " is too large for the values vector" );

            }

            _macroSurfaceForces.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                                values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );
            _macroSurfaceForces[ n->first ] *= _sign;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroExternalForces( const unsigned int &increment ){
        /*!
         * Extract the node macro internal forces at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        stringVector variableKeys =
            {
                "F1", "F2", "F3",
                "C11", "C12", "C13",
                "C21", "C22", "C23",
                "C31", "C32", "C33"
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "external_force_variable_names";
        YAML::Node configuration = _config[ "macroscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        errorOut error = inputFileProcessor::extractDataFileProperties( _macroscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, _macroExternalForceFlag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMacroExternalForces",
                                             "Error in the extraction of the macro external forces" );
            result->addNext( error );
            return result;

        }

        if ( !_macroExternalForceFlag ){

            if( _macroSurfaceForceFlag && _macroBodyForceFlag ){

                _macroExternalForces.clear( );
                _macroExternalForces.reserve( _macroGlobalNodeIDOutputIndex.size( ) );
                for ( auto n = _macroGlobalNodeIDOutputIndex.begin( ); n != _macroGlobalNodeIDOutputIndex.end( ); n++ ){

                    _macroExternalForces.emplace( n->first, _macroSurfaceForces[ n->first ] + _macroBodyForces[ n->first ] );

                }
                _macroExternalForceFlag = true;

            }
            else if ( _macroSurfaceForceFlag ){

                _macroExternalForces = _macroSurfaceForces;
                _macroExternalForceFlag = true;

            }
            else if ( _macroBodyForceFlag ){

                _macroExternalForces = _macroBodyForces;
                _macroExternalForceFlag = true;

            }
            else{
                _macroExternalForces.clear( );
            }

        }
        else{

            _macroExternalForces.clear( );
            _macroExternalForces.reserve( _macroGlobalNodeIDOutputIndex.size( ) );
            floatType _sign = _config[ "coupling_initialization" ][ "macro_external_force_sign" ].as< floatType >( );

            for ( auto n = _macroGlobalNodeIDOutputIndex.begin( ); n != _macroGlobalNodeIDOutputIndex.end( ); n++ ){
    
                if ( n->second >= values.size( ) ){
    
                    return new errorNode( "extractMacroExternalForces",
                                          "The macro external force vector is too short for the required index" );
    
                }
    
                floatVector ef( values.begin( ) + variableKeys.size( ) * n->second,
                                values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) );
                ef *= _sign;
    
                _macroExternalForces.emplace( n->first, ef );
    
            }

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroInertialForces( const unsigned int &increment ){
        /*!
         * Extract the node macro inertial forces at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        stringVector variableKeys =
            {
                "F1", "F2", "F3",
                "C11", "C12", "C13",
                "C21", "C22", "C23",
                "C31", "C32", "C33"
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "inertial_force_variable_names";
        YAML::Node configuration = _config[ "macroscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        errorOut error = inputFileProcessor::extractDataFileProperties( _macroscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, _macroInertialForceFlag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMacroInertialForces",
                                             "Error in the extraction of the macro inertial forces" );
            result->addNext( error );
            return result;

        }

        if ( !_macroInertialForceFlag ){

            _macroInertialForces.clear( );
            return NULL;

        }

        floatType _sign = _config[ "coupling_initialization" ][ "macro_inertial_force_sign" ].as< floatType >( );
        _macroInertialForces.clear( );
        _macroInertialForces.reserve( _macroGlobalNodeIDOutputIndex.size( ) );
        for ( auto n = _macroGlobalNodeIDOutputIndex.begin( ); n != _macroGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMacroInertialForces",
                                      "The index required for macro node " + std::to_string( n->first ) +
                                      " is too large for the values vector" );

            }

            _macroInertialForces.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                                 values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );
            _macroInertialForces[ n->first ] *= _sign;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroLumpedMassMatrix( const unsigned int &increment ){
        /*!
         * Extract the macro lumped-mass matrix
         *
         * :param const unsigned int &increment: The current increment
         */

        stringVector variableKeys =
            {
                "M1", "M2", "M3",
                "I11", "I12", "I13",
                "I21", "I22", "I23",
                "I31", "I32", "I33"
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "lumped_mass_matrix_variable_names";
        YAML::Node configuration = _config[ "macroscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        errorOut error = inputFileProcessor::extractDataFileProperties( _macroscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, _macroLumpedMassMatrixFlag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMacroLumpedMassMatrix",
                                             "Error in the extraction of the macro lumped mass matrix" );
            result->addNext( error );
            return result;

        }

        if ( !_macroLumpedMassMatrixFlag ){

            _macroLumpedMassMatrix.clear( );
            return NULL;

        }

        _macroLumpedMassMatrix.clear( );
        _macroLumpedMassMatrix.reserve( _macroGlobalNodeIDOutputIndex.size( ) );
        floatType _sign = _config[ "coupling_initialization" ][ "macro_body_force_sign" ].as< floatType >( );
        for( auto n = _macroGlobalNodeIDOutputIndex.begin( ); n != _macroGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMacroLumpedMassMatrix",
                                      "The index required by macro node " + std::to_string( n->first ) + " is too large for the values vector" );

            }

            _macroLumpedMassMatrix.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                                   values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );
            _macroLumpedMassMatrix[ n->first ] *= _sign;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractDataFileProperties( std::shared_ptr< dataFileInterface::dataFileBase > &dataFile,
                                                            const unsigned int &increment, const stringVector &variableKeys,
                                                            const std::string &dataType,
                                                            const bool &populateWithNullOnUndefined,
                                                            const std::string &configurationName,
                                                            YAML::Node &configuration, bool &populatedFlag, floatVector &properties ){
        /*!
         * Extract properties from the datafile. The variableKeys are assumed to be scalar coefficients of some property.
         *
         * :param std::shared_ptr< dataFileInterface::dataFileBase > &dataFile: The datafile to read the information from.
         * :param const unsigned int &increment: The increment at which to extract the data
         * :param const stringVector &variableKeys: The names of the coefficients to be extracted
         * :param const std::string &dataType: The type of data "Node" and "Cell" are currently supported
         * :param const bool &populateWithNullOnUndefined: If true and the variable definition is not found in the datafile
         *     the variableKeys will be set equal to "NULL". If false, an error will be raised
         * :param const std::string &configurationName: The name of the configuration node
         * :param YAML::Node &configuration: The YAML Node which discribes the configuration. This node should be something like
         *     _config[ "macroscale_definition" ][ "variable_names" ] where the map from variableKeys to variable names
         *     is defined as key: name
         * :param bool &populatedFlag: A flag which is set to true if the output is populated from the datafile. False if it is
         *     set to a zero vector the length of variableKeys ( the default )
         * :param floatVector &properties: The properties vector extracted from the datafile
         */

        //Check if the variable property names have been defined
        populatedFlag = false;

        if ( variableKeys.size( ) == 0 ){

            return new errorNode( "extractDataFileProperties", "No variable keys have been defined" );

        }

        bool missingKey = false;
        if ( !configuration ){

            missingKey = true;

        }
        else{

            if ( !configuration.IsScalar( ) ){

                for ( auto vK = variableKeys.begin( ); vK != variableKeys.end( ); vK++ ){
        
                    if ( configuration[ *vK ].IsNull( ) ){
        
                        missingKey = true;
                        break;
        
                    }
        
                }

            }
            else if ( configuration.IsSequence( ) ){

                std::string output = "The configuration for " + configurationName + " is set as a sequence when it shouldn't be.\n";
                output            += "The formation should be:\n";
                output            += "  " + configurationName + ":\n";
                for ( auto vK = variableKeys.begin( ); vK != variableKeys.end( ); vK++ ){

                    output += "          " + *vK + ": " + *vK + "_variable_name\n";

                }
                

            }
            else{

                missingKey = true;

            }

        }

        if ( ( !configuration ) ||
             ( missingKey ) ||
             ( configuration[ variableKeys[ 0 ] ].as< std::string >( ).compare( "NULL" ) == 0 ) ){

            if ( populateWithNullOnUndefined ){

                configuration = YAML::Node( );

                for ( auto vK = variableKeys.begin( ); vK != variableKeys.end( ); vK++ ){
   
                    configuration[ *vK ] = "NULL";
    
                }
   
                properties = floatVector( variableKeys.size( ), 0 ); //Set the properties to zero
                return NULL;

            }
            else{

                std::string output  = "The configuration is not fully defined for " + configurationName + ".\n";
                output             += "  The definition of the variable components should be performed as:\n";
                output             += configurationName + ":\n";
                for ( auto vK = variableKeys.begin( ); vK != variableKeys.end( ); vK++ ){

                    output += "          " + *vK + ": " + *vK + "_variable_name\n";

                }
                return new errorNode( "extractDataFileProperties", output );

            }

            return NULL;

        }

        //Get the variable names
        stringVector variableNames( variableKeys.size( ) );
        for ( auto vK  = variableKeys.begin( );  vK != variableKeys.end( ); vK++ ){

            //Get the variable name associated with this component
            if ( !configuration[ *vK ].IsScalar( ) ){

                std::string output  = "The definition of " + configurationName + " variable key " + *vK +
                                      " is either not defined in the input file or incorrectly defined.\n";
                output             += "The definition of the variable components should be performed as:\n";
                output             += configurationName + ":\n";
                for ( auto vK = variableKeys.begin( ); vK != variableKeys.end( ); vK++ ){

                    output += "          " + *vK + ": " + *vK + "_variable_name\n";

                }

                return new errorNode( "extractDataFileProperties", output );

            }
            
            variableNames[ vK - variableKeys.begin( ) ]
                = configuration[ *vK ].as< std::string >( );

        }

        //Extract the velocity vector
        errorOut error = dataFile->getSolutionVectorDataFromComponents( increment, variableNames,
                                                                        dataType, properties );

        if ( error ){

            errorOut result = new errorNode( "extractDataFileProperties",
                                             "Error in the extraction of the datafile properties" );
            result->addNext( error );
            return result;

        }

        populatedFlag = true;

        return NULL;
    }

    errorOut inputFileProcessor::extractMicroStresses( const unsigned int &increment ){
        /*!
         * Extract the micro-stresses at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        stringVector variableKeys =
            {
                "s11", "s12", "s13", "s21", "s22", "s23", "s31", "s32", "s33"
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = false;
        bool _tmpFlag;

        std::string configurationName = "stress_variable_names";
        YAML::Node configuration = _config[ "microscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        errorOut error = inputFileProcessor::extractDataFileProperties( _microscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, _tmpFlag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMicroStresses",
                                             "Error in the extraction of the micro stresses" );
            result->addNext( error );
            return result;

        }

        _microStresses.clear( );
        _microStresses.reserve( _microGlobalNodeIDOutputIndex.size( ) );
        for ( auto n = _microGlobalNodeIDOutputIndex.begin( ); n != _microGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMicroVelocities",
                                      "The micro stress vector is too short for the required index" );

            }

            _microStresses.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                           values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroInternalForces( const unsigned int &increment ){
        /*!
         * Extract the internal force vector for the micro nodes
         *
         * :param const unsigned int &increment: The increment at which to extract the vector
         */

        stringVector variableKeys =
            {
                "F1", "F2", "F3"
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "internal_force_variable_names";
        YAML::Node configuration = _config[ "microscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        errorOut error = inputFileProcessor::extractDataFileProperties( _microscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, _microInternalForceFlag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMicroInternalForces",
                                             "Error in the extraction of the micro internal forces" );
            result->addNext( error );
            return result;

        }

        if ( !_microInternalForceFlag ){

            _microInternalForces.clear( );
            return NULL;

        }


        _microInternalForces.clear( );
        _microInternalForces.reserve( _microGlobalNodeIDOutputIndex.size( ) );
        floatType _sign = _config[ "coupling_initialization" ][ "micro_internal_force_sign" ].as< floatType >( );
        for ( auto n = _microGlobalNodeIDOutputIndex.begin( ); n != _microGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMicroVelocities",
                                      "The micro internal force vector is too short for the required index" );

            }

            _microInternalForces.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                                 values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );
            _microInternalForces[ n->first ] *= _sign;

        }

        return NULL;
    }

    errorOut inputFileProcessor::extractMicroInertialForces( const unsigned int &increment ){
        /*!
         * Extract the inertial force vector for the micro nodes
         *
         * :param const unsigned int &increment: The increment at which to extract the vector
         */

        stringVector variableKeys =
            {
                "F1", "F2", "F3"
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = true;

        std::string configurationName = "inertial_force_variable_names";
        YAML::Node configuration = _config[ "microscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        errorOut error = inputFileProcessor::extractDataFileProperties( _microscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, _microInertialForceFlag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMicroInertialForces",
                                             "Error in the extraction of the micro inertial forces" );
            result->addNext( error );
            return result;

        }

        if ( !_microInertialForceFlag ){

            _microInertialForces.clear( );
            return NULL;

        }


        _microInertialForces.clear( );
        _microInertialForces.reserve( _microGlobalNodeIDOutputIndex.size( ) );
        floatType _sign = _config[ "coupling_initialization" ][ "micro_inertial_force_sign" ].as< floatType >( );
        for ( auto n = _microGlobalNodeIDOutputIndex.begin( ); n != _microGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMicroVelocities",
                                      "The micro inertial force vector is too short for the required index" );

            }


            _microInertialForces.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                                 values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );
            _microInertialForces[ n->first ] *= _sign;

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

        //Initialize the size of the micro volumes map
        _microVolumes.clear( );
        _microVolumes.reserve( _microGlobalNodeIDOutputIndex.size( ) );

        //Get the values of the micro volumes from the output file
        floatVector values;
        errorOut error = _microscale->getSolutionData( increment,
                                                       _config[ "microscale_definition" ][ "volume_variable_name" ].as< std::string >( ),
                                                       "Node", values );

        if ( error ){

            errorOut result = new errorNode( "extractMicroNodeVolumes", "Error in extraction of the micro volumes" );
            result->addNext( error );
            return result;

        }

        for ( auto n = _microGlobalNodeIDOutputIndex.begin( ); n != _microGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMicroNodeVolumes", "The volume vector is too short for the required index" );

            }

            _microVolumes.emplace( n->first, values[ n->second ] );

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroArlequinWeights( const unsigned int &increment ){
        /*!
         * Extract the macro Arlequin weights for the macro nodes at the indicated increment
         *
         * :param const unsigned int &increment: The current increment
         */

        //Check if the volume name has been defined
        if ( !_config[ "coupling_initialization" ][ "arlequin_weighting_variable_name" ] ){
        
            return new errorNode( "extractMacroArlequinWeights", "The Arlequin weight variable name is not defined" );

        }

        //Initialize the size of the macro Arlequin weights map
        _macroArlequinWeights.clear( );
        _macroArlequinWeights.reserve( _macroGlobalNodeIDOutputIndex.size( ) );

        //Get the values of the micro volumes from the output file
        floatVector values;

        if ( !_isFiltering ){

            errorOut error = _macroscale->getSolutionData( increment,
                                                           _config[ "coupling_initialization" ][ "arlequin_weighting_variable_name" ].as< std::string >( ),
                                                           "Node", values );
    
            if ( error ){
    
                errorOut result = new errorNode( "extractMacroArlequinWeights", "Error in extraction of the macro Arlequin nodal weights" );
                result->addNext( error );
                return result;
    
            }

        }

        for ( auto n = _macroGlobalNodeIDOutputIndex.begin( ); n != _macroGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( _isFiltering ){

                _macroArlequinWeights.emplace( n->first, 0 );

            }
            else{

                if ( n->second >= values.size( ) ){
    
                    return new errorNode( "extractMacroArlequinWeights", "The Arlequin weights vector is too short for the required index" );
    
                }
    
                _macroArlequinWeights.emplace( n->first, values[ n->second ] );

            }

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroDisplacements( const unsigned int &increment ){
        /*!
         * Extract the node micro-displacements at the indicated increment
         *
         * TODO: Move this and other vector / tensor extraction routines to a common function
         *
         * :param const unsigned int &increment: The current increment
         */

        bool flag;

        errorOut error = extractMicroDisplacements( increment, flag, _microDisplacements );

        if ( error ){

            errorOut result = new errorNode( "extractMicroDisplacements",
                                             "Error in the extraction of the micro displacements" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroDisplacements( const unsigned int &increment, bool &flag,
                                                            std::unordered_map< uIntType, floatVector > &microDisplacements ){
        /*!
         * Extract the node micro-displacements at the indicated increment
         *
         * TODO: Move this and other vector / tensor extraction routines to a common function
         *
         * :param const unsigned int &increment: The current increment
         * :param bool &flag: The flag to indicate if the values are defined in the input file
         * :param std::unordered_map< uIntType, floatVector > &microDisplacements: The map to store the displacements in
         */

        stringVector variableKeys =
            {
                "u1", "u2", "u3",
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = false;

        std::string configurationName = "displacement_variable_names";
        YAML::Node configuration = _config[ "microscale_definition" ][ configurationName.c_str( ) ];

        floatVector values;
        errorOut error = inputFileProcessor::extractDataFileProperties( _microscale, increment, variableKeys, dataType,
                                                                        populateWithNullOnUndefined, configurationName,
                                                                        configuration, flag, values );

        if ( error ){

            errorOut result = new errorNode( "extractMicroDisplacements",
                                             "Error in the extraction of the micro displacements" );
            result->addNext( error );
            return result;

        }

        microDisplacements.clear( );
        microDisplacements.reserve( _microGlobalNodeIDOutputIndex.size( ) );

        for ( auto n = _microGlobalNodeIDOutputIndex.begin( ); n != _microGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= values.size( ) ){

                return new errorNode( "extractMicroDisplacements",
                                      "The displacement vector is too short for the required index" );

            }

            microDisplacements.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                               values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroDisplacements( const unsigned int &increment ){
        /*!
         * Extract the positions of the nodes in the macro domain.
         *
         * TODO: Replace with generalized form
         *
         * :param const unsigned int &increment: The current increment
         */

        stringVector variableKeys =
            {
                "u1", "u2", "u3",
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = false;
        bool flag;

        std::string configurationName = "displacement_variable_names";
        YAML::Node configuration = _config[ "macroscale_definition" ][ configurationName.c_str( ) ];

        floatVector values;
        if ( !_isFiltering ){

            errorOut error = inputFileProcessor::extractDataFileProperties( _macroscale, increment, variableKeys, dataType,
                                                                            populateWithNullOnUndefined, configurationName,
                                                                            configuration, flag, values );
    
            if ( error ){
    
                errorOut result = new errorNode( "extractMacroDisplacements",
                                                 "Error in the extraction of the macro displacements" );
                result->addNext( error );
                return result;
    
            }
        }

        _macroDisplacements.clear( );
        _macroDisplacements.reserve( _macroGlobalNodeIDOutputIndex.size( ) );

        for ( auto n = _macroGlobalNodeIDOutputIndex.begin( ); n != _macroGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( _isFiltering ){

                _macroDisplacements.emplace( n->first, floatVector( variableKeys.size( ), 0 ) );

            }
            else{

                if ( n->second >= values.size( ) ){
    
                    return new errorNode( "extractMacroDisplacements",
                                          "The outputDOF vector is too short for the required index" );
    
                }

                _macroDisplacements.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                                    values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );
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

        bool flag;
        errorOut error = extractMacroDispDOFVector( increment, flag, _macroDispDOFVector );

        if ( error ){

            errorOut result = new errorNode( "extractMicroDispDOFVector",
                                             "Error in the extraction of the micro displacement degree of freedom vector" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroDispDOFVector( const unsigned int &increment, bool &flag, std::unordered_map< uIntType, floatVector > &macroDispDOFVector ){
        /*!
         * Extract the displacement degrees of freedom of the nodes in the macro domain.
         *
         * TODO: Replace with generalized form
         *
         * :param const unsigned int &increment: The current increment
         * :param bool &flag: The flag to indicate if the values are defined in the input file
         * :param std::unordered_map< uIntType, floatVector > &microDispDOFVector: The vector to store the displacements DOF vector in
         */
        stringVector variableKeys =
            {
                 "u1", "u2", "u3",
                 "phi11", "phi12", "phi13",
                 "phi21", "phi22", "phi23",
                 "phi31", "phi32", "phi33"
            };

        std::string dataType = "Node";

        bool populateWithNullOnUndefined = false;

        std::string configurationName = "displacement_variable_names";
        YAML::Node configuration = _config[ "macroscale_definition" ][ configurationName.c_str( ) ];
        floatVector values;

        if ( !_isFiltering ){

            errorOut error = inputFileProcessor::extractDataFileProperties( _macroscale, increment, variableKeys, dataType,
                                                                            populateWithNullOnUndefined, configurationName,
                                                                            configuration, flag, values );
    
            if ( error ){
    
                errorOut result = new errorNode( "extractMicroDispDOFVector",
                                                 "Error in the extraction of the micro displacement degree of freedom vector" );
                result->addNext( error );
                return result;
    
            }

        }

        macroDispDOFVector.clear( );
        macroDispDOFVector.reserve( _macroGlobalNodeIDOutputIndex.size( ) );

        for ( auto n = _macroGlobalNodeIDOutputIndex.begin( ); n != _macroGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( _isFiltering ){

                macroDispDOFVector.emplace( n->first, floatVector( variableKeys.size( ), 0 ) );

            }
            else{

                if ( n->second >= values.size( ) ){
    
                    return new errorNode( "extractMacroDispDOFVector",
                                          "The DOF vector is too short for the required index" );
    
                }

                macroDispDOFVector.emplace( n->first, floatVector( values.begin( ) + variableKeys.size( ) * n->second,
                                                                   values.begin( ) + variableKeys.size( ) * ( n->second + 1 ) ) );

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

        //Initialize the size of the micro reference positions map
        _microNodeReferencePositions.clear( );
        _microNodeReferencePositions.reserve( _microGlobalNodeIDOutputIndex.size( ) );
        floatVector referencePositions;
        uIntVector referenceConnectivity;
        uIntVector connectivityCellIndices;
        uIntType cellCounts;

        uIntVector elementConnectivity;

        errorOut error = _microscale->getMeshData( increment, referencePositions, referenceConnectivity,
                                                   connectivityCellIndices, cellCounts );

        if ( error ){
            errorOut result = new errorNode( "extractMicroMeshData", "Error in the extraction of the micro-mesh information" );
            result->addNext( error );
            return result;
        }

        for ( auto n = _microGlobalNodeIDOutputIndex.begin( ); n != _microGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= referencePositions.size( ) ){

                return new errorNode( "extractReferenceMicroMeshData",
                                      "The reference positions vector is too short for the required index" );

            }

            _microNodeReferencePositions.emplace( n->first, floatVector( referencePositions.begin( ) + _dim * n->second,
                                                                         referencePositions.begin( ) + _dim * ( n->second + 1 )
                                                                       ) );

        }

        return NULL;

    }

    errorOut inputFileProcessor::extractReferenceMacroMeshData( const unsigned int &increment ){
        /*!
         * Extract the mesh data for the macro-scale domain in the reference configuration
         *
         * :param const unsigned int &increment: The increment at which to extract the macro-mesh data
         */

        //Initialize the size of the micro reference positions map
        _macroNodeReferencePositions.clear( );
        _macroNodeReferencePositions.reserve( _macroGlobalNodeIDOutputIndex.size( ) );
        floatVector referencePositions;
        uIntVector referenceConnectivity;
        uIntVector connectivityCellIndices;
        uIntType cellCounts;

        uIntVector elementConnectivity;

        errorOut error = _macroscale->getMeshData( increment, referencePositions, referenceConnectivity,
                                                   connectivityCellIndices, cellCounts );

        if ( error ){
            errorOut result = new errorNode( "extractReferenceMacroMeshData", "Error in the extraction of the macro-mesh information" );
            result->addNext( error );
            return result;
        }

        for ( auto n = _macroGlobalNodeIDOutputIndex.begin( ); n != _macroGlobalNodeIDOutputIndex.end( ); n++ ){

            if ( n->second >= referencePositions.size( ) ){

                return new errorNode( "extractReferenceMacroMeshData",
                                      "The reference positions vector is too short for the required index" );

            }

            _macroNodeReferencePositions.emplace( n->first, floatVector( referencePositions.begin( ) + _dim * n->second,
                                                                         referencePositions.begin( ) + _dim * ( n->second + 1 )
                                                                       ) );

        }

        uIntVector elementIds;
        std::string elementIDName = _config[ "macroscale_definition" ][ "cell_id_variable_name" ].as< std::string >( );
        error = _macroscale->getCellIds( increment, elementIDName, elementIds );

        if ( cellCounts > elementIds.size( ) ){

            return new errorNode( "extractReferenceMacroMeshData",
                                  "The connectivity and the number of elements in the macroscale aren't consistent" );

        }

        _macroNodeReferenceConnectivity.clear( );
        _macroNodeReferenceConnectivity.reserve( _free_macro_cell_ids.size( ) + _ghost_macro_cell_ids.size( ) );

        uIntType startIndex, stopIndex;
        for ( auto cell = _free_macro_cell_ids.begin( ); cell != _free_macro_cell_ids.end( ); cell++ ){

            auto index = std::find( elementIds.begin( ), elementIds.end( ), *cell );

            if ( index == elementIds.end( ) ){

                return new errorNode( "Free macro cell " + std::to_string( *cell ) + " was not found in the element ids" );

            }

            startIndex = connectivityCellIndices[ *index ];
            if ( ( *index + 1 ) >= connectivityCellIndices.size( ) ){

                stopIndex = referenceConnectivity.size( );

            }
            else{

                stopIndex = connectivityCellIndices[ *index + 1 ];

            }

            _macroNodeReferenceConnectivity.emplace( *cell, uIntVector( referenceConnectivity.begin( ) + startIndex,
                                                                        referenceConnectivity.begin( ) + stopIndex ) );


        }

        for ( auto cell = _ghost_macro_cell_ids.begin( ); cell != _ghost_macro_cell_ids.end( ); cell++ ){

            auto index = std::find( elementIds.begin( ), elementIds.end( ), *cell );

            if ( index == elementIds.end( ) ){

                return new errorNode( "Ghost macro cell " + std::to_string( *cell ) + " was not found in the element ids" );

            }

            startIndex = connectivityCellIndices[ *index ];
            if ( ( *index + 1 ) >= connectivityCellIndices.size( ) ){

                stopIndex = referenceConnectivity.size( );

            }
            else{

                stopIndex = connectivityCellIndices[ *index + 1 ];

            }

            _macroNodeReferenceConnectivity.emplace( *cell, uIntVector( referenceConnectivity.begin( ) + startIndex,
                                                                        referenceConnectivity.begin( ) + stopIndex ) );

        }

        return NULL;

    }

    const floatType* inputFileProcessor::getMicroTime( ){
        /*!
         * Get the timestamp of the micro increment
         */

        return &_microTime;
    }

    const floatType* inputFileProcessor::getPreviousMicroTime( ){
        /*!
         * Get the timestamp of the micro increment
         */

        if ( _extractPreviousDOFValues ){

            return &_previousMicroTime;

        }
        else{

            return &_microTime;
            
        }

    }

    const std::unordered_map< uIntType, floatType >* inputFileProcessor::getMicroDensities( ){
        /*!
         * Get a pointer to the density
         */

        return &_microDensities;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMicroBodyForces( ){
        /*!
         * Get a pointer to the micro body forces
         */

        return &_microBodyForces;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMicroSurfaceForces( ){
        /*!
         * Get a pointer to the micro surface forces
         */

        return &_microSurfaceForces;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMicroExternalForces( ){
        /*!
         * Get a pointer to the micro body forces
         */

        return &_microExternalForces;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMicroVelocities( ){
        /*!
         * Get a pointer to the micro velocities
         */

        return &_microVelocities;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMicroAccelerations( ){
        /*!
         * Get a pointer to the micro accelerations
         */

        return &_microAccelerations;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getPreviousMicroDisplacements( ){
        /*!
         * Get a pointer to the previous micro displacements
         */

        if ( _extractPreviousDOFValues ){

            return &_previousMicroDisplacements;

        }
        else{

            return &_microDisplacements;

        }
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getPreviousMicroVelocities( ){
        /*!
         * Get a pointer to the previous micro velocities
         */

        if ( _extractPreviousDOFValues ){

            return &_previousMicroVelocities;

        }
        else{

            return &_microVelocities;

        }
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getPreviousMicroAccelerations( ){
        /*!
         * Get a pointer to the micro accelerations
         */

        if ( _extractPreviousDOFValues ){

            return &_previousMicroAccelerations;

        }
        else{

            return &_microAccelerations;

        }
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMicroStresses( ){
        /*!
         * Get a pointer to the micro stresses
         */

        return &_microStresses;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMicroInternalForces( ){
        /*!
         * Get a pointer to the micro internal forces
         */

        return &_microInternalForces;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMicroInertialForces( ){
        /*!
         * Get a pointer to the micro inertial forces
         */

        return &_microInertialForces;
    }

    const std::unordered_map< uIntType, floatType >* inputFileProcessor::getMicroVolumes( ){
        /*!
         * Get a pointer to the volumes
         */

        return &_microVolumes;
    }

    const std::unordered_map< uIntType, floatType >* inputFileProcessor::getMicroWeights( ){
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

    const std::unordered_map< std::string, uIntType > *inputFileProcessor::getMicroDomainSurfaceApproximateSplitCount( ){
        /*!
         * Get a constant reference to the micro-domain to approximate surface split count map
         */
        return &_microDomainSurfaceCount;
    }

//    const uIntVector* inputFileProcessor::getFreeMicroSurfaceApproximateSplitCount( ){
//        /*!
//         * Get the free micro-surface approximate split counts. I.e. the number of 
//         * surfaces a given micro domain should be split into ( approximately )
//         */
//
//        return &_free_micro_surface_approximate_split_count;
//    }
//
//    const uIntVector* inputFileProcessor::getGhostMicroSurfaceApproximateSplitCount( ){
//        /*!
//         * Get the ghost micro-surface approximate split counts. I.e. the number of 
//         * surfaces a given micro domain should be split into ( approximately )
//         */
//
//        return &_ghost_micro_surface_approximate_split_count;
//    }
//
//    const stringVector* inputFileProcessor::getFreeMicroSurfaceNames( ){
//        /*!
//         * Get the free domain names
//         */
//
//        return &_free_micro_surface_sets;
//    }
//
//    const stringVector* inputFileProcessor::getGhostMicroSurfaceNames( ){
//        /*!
//         * Get the ghost domain names
//         */
//
//        return &_ghost_micro_surface_sets;
//    }
//
//    const stringVector* inputFileProcessor::getNonOverlappedMicroSurfaceNames( ){
//        /*!
//         * Get the ghost domain names
//         */
//
//        return &_non_overlapped_micro_surface_sets;
//    }

    const YAML::Node inputFileProcessor::getCouplingInitialization( ){
        /*!
         * Get the coupling initialization from the configuration file.
         */

        return _config[ "coupling_initialization" ];
    }

    bool inputFileProcessor::extractPreviousDOFValues( ){
        /*!
         * Get the flag for if the previous accelerations and velocities were
         * supposed to have been extracted.
         */

        return _extractPreviousDOFValues;

    }

    errorOut inputFileProcessor::checkCouplingInitialization( ){
        /*!
         * Check the coupling initialization
         */

        if ( !_config[ "coupling_initialization" ][ "type" ] ){

            _config[ "coupling_initialization" ][ "type" ] = "use_first_increment";

        }

        if ( !_config[ "coupling_initialization" ][ "projection_type" ] ){

            _config[ "coupling_initialization" ][ "projection_type" ] = "averaged_l2_projection";

        }

        if ( _config[ "coupling_initialization" ][ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){

            _computeMicroShapeFunctions = true;

        }

        if ( _config[ "coupling_initialization" ][ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){

            std::cerr << "WARNING: direct_projection can give unexpected results.\n";
            std::cerr << "         It is suggested that either l2_projection or\n"; 
            std::cerr << "         averaged_l2_projection ( recommended ) are used\n";
            std::cerr << "         for small problems and arlequin is used for large\n";
            std::cerr << "         problems.\n";

        }

        _useArlequinMethod = false;
        if ( _config[ "coupling_initialization" ][ "projection_type" ].as< std::string >( ).compare( "arlequin" )  == 0 ){

            if ( !_config[ "coupling_initialization" ][ "arlequin_weighting_variable_name" ] ){

                return new errorNode( "checkCouplingInitialization",
                                      "For 'arlequin' projection the weighting variable name at the macro nodes must be defined in 'arlequin_weighting_variable_name' in 'coupling_initialization'" );

            }

            if ( !_config[ "coupling_initialization" ][ "arlequin_penalty_parameter" ] ){

                _config[ "coupling_initialization" ][ "arlequin_penalty_parameter" ] = 1;
                _arlequinPenaltyParameter = 1;

            }
            else{

                _arlequinPenaltyParameter = _config[ "coupling_initialization" ][ "arlequin_penalty_parameter" ].as< floatType >( );

            }

            if ( !_config[ "coupling_initialization" ][ "arlequin_update_penalty_parameter" ] ){

                _config[ "coupling_initialization" ][ "arlequin_update_penalty_parameter" ] = 1;
                _arlequinUpdatePenaltyParameter = 1;

            }
            else{

                _arlequinUpdatePenaltyParameter = _config[ "coupling_initialization" ][ "arlequin_update_penalty_parameter" ].as< floatType >( );
            }

            _useArlequinMethod = true;

        }

        if ( !_config[ "coupling_initialization" ][ "use_reconstructed_mass_centers" ] ){

            _config[ "coupling_initialization" ][ "use_reconstructed_mass_centers" ] = true;

        }

        if ( !_config [ "coupling_initialization" ][ "potential_energy_weighting_factor" ] ){

            _config[ "coupling_initialization" ][ "potential_energy_weighting_factor" ] = 0.5; //Default to 0.5

        }

        if ( !_config[ "coupling_initialization" ][ "potential_energy_partitioning_coefficient" ] ){

            _config[ "coupling_initialization" ][ "potential_energy_partitioning_coefficient" ][ "type" ] = "volume_fraction"; //Default to volume fraction

        }

        if ( !_config [ "coupling_initialization" ][ "kinetic_energy_weighting_factor" ] ){

            _config[ "coupling_initialization" ][ "kinetic_energy_weighting_factor" ] = 0.5; //Default to 0.5

        }

        if ( !_config[ "coupling_initialization" ][ "kinetic_energy_partitioning_coefficient" ] ){

            _config[ "coupling_initialization" ][ "kinetic_energy_partitioning_coefficient" ][ "type" ] = "volume_fraction"; //Default to volume fraction

        }

        if ( !_config[ "coupling_initialization" ][ "macro_proportionality_coefficient" ] ){

            _config[ "coupling_initialization" ][ "macro_proportionality_coefficient" ] = 1e-3; //Default to 1e-3

        }

        if ( !_config[ "coupling_initialization" ][ "micro_proportionality_coefficient" ] ){

            _config[ "coupling_initialization" ][ "micro_proportionality_coefficient" ] = 1e-3; //Default to 1e-3

        }

        if ( !_config[ "coupling_initialization" ][ "macro_internal_force_sign" ] ){

            _config[ "coupling_initialization" ][ "macro_internal_force_sign" ] = -1; //Default to -1 to be consistent with the micromorphic implementation

        }

        if ( !_config[ "coupling_initialization" ][ "macro_inertial_force_sign" ] ){

            _config[ "coupling_initialization" ][ "macro_inertial_force_sign" ] = 1; //Default to 1

        }

        if ( !_config[ "coupling_initialization" ][ "macro_body_force_sign" ] ){

            _config[ "coupling_initialization" ][ "macro_body_force_sign" ] = 1; //Default to 1

        }

        if ( !_config[ "coupling_initialization" ][ "macro_surface_force_sign" ] ){

            _config[ "coupling_initialization" ][ "macro_surface_force_sign" ] = 1; //Default to 1

        }

        if ( !_config[ "coupling_initialization" ][ "macro_external_force_sign" ] ){

            _config[ "coupling_initialization" ][ "macro_external_force_sign" ] = 1; //Default to 1

        }

        if ( !_config[ "coupling_initialization" ][ "micro_internal_force_sign" ] ){

            _config[ "coupling_initialization" ][ "micro_internal_force_sign" ] = 1; //Default to 1

        }

        if ( !_config[ "coupling_initialization" ][ "micro_inertial_force_sign" ] ){

            _config[ "coupling_initialization" ][ "micro_inertial_force_sign" ] = 1; //Default to 1

        }

        if ( !_config[ "coupling_initialization" ][ "micro_body_force_sign" ] ){

            _config[ "coupling_initialization" ][ "micro_body_force_sign" ] = 1; //Default to 1

        }

        if ( !_config[ "coupling_initialization" ][ "micro_surface_force_sign" ] ){

            _config[ "coupling_initialization" ][ "micro_surface_force_sign" ] = 1; //Default to 1

        }

        if ( !_config[ "coupling_initialization" ][ "micro_external_force_sign" ] ){

            _config[ "coupling_initialization" ][ "micro_external_force_sign" ] = 1; //Default to 1

        }

        if ( _config[ "coupling_initialization" ][ "extract_previous_dof_values" ] ){

            if ( !_config[ "coupling_initialization" ][ "extract_previous_dof_values" ].IsScalar( ) ){

                return new errorNode( "'extract_previous_dof_values' must be a boolean value" );

            }
            else if ( _config[ "coupling_initialization" ][ "extract_previous_dof_values" ].as< bool >( ) ){

                _extractPreviousDOFValues = true;

            }

            if ( _extractPreviousDOFValues ){

                if ( !_config[ "coupling_initialization" ][ "previous_micro_increment" ] ){

                    return new errorNode( "checkCouplingInitialization",
                                          "'previous_micro_increment' is not defined in 'coupling_initialization' when the user has requested that the previous values of the degrees of freedom are extracted" );

                }
                else if(  !_config[ "coupling_initialization" ][ "previous_micro_increment" ].IsScalar( ) ){

                    return new errorNode( "checkCouplingInitialization",
                                          "'previous_micro_increment' must be defined as a scalar integer value indicating the increment ( i.e. timestep number ) which defines the previous dof values at the micro scale" );

                }

                if ( !_config[ "coupling_initialization" ][ "previous_macro_increment" ] ){

                    return new errorNode( "checkCouplingInitialization",
                                          "'previous_macro_increment' is not defined in 'coupling_initialization' when the user has requested that the previous dof values are extracted" );

                }
                else if(  !_config[ "coupling_initialization" ][ "previous_macro_increment" ].IsScalar( ) ){

                    return new errorNode( "checkCouplingInitialization",
                                          "'previous_macro_increment' must be defined as a scalar integer value indicating the increment ( i.e. timestep number ) which defines the last converged dof values at the macro scale" );

                }

            }

        }
        else{

            _config[ "coupling_initialization" ][ "extract_previous_dof_values" ] = false;

        }

        if ( _config[ "coupling_initialization" ][ "update_displacement" ] ){

            if ( !_config[ "coupling_initialization" ][ "extract_previous_dof_values" ].as< bool > ( ) ){

                if ( !_config[ "coupling_initialization" ][ "update_displacement" ][ "Dt" ] ){

                    return new errorNode( "checkCouplingInitialization",
                                          "If the previous DOF values are not to be extracted and the displacement is to be updated, 'Dt' must be defined under 'update_displacement'" );

                }
                else{

                    _Dt = _config[ "coupling_initialization" ][ "update_displacement" ][ "Dt" ].as< floatType >( );

                }

            }
            else if ( _config[ "coupling_initialization" ][ "update_displacement" ][ "Dt" ] ){

                std::cerr << "WARNING: Dt is specified when the previous increment has been indicated.\n";
                std::cerr << "         The Dt in the input file will be ignored\n";
                _config[ "coupling_initialization" ][ "update_displacement" ][ "Dt" ] = "NULL";

            }

            if ( !_config[ "coupling_initialization" ][ "update_displacement" ][ "Newmark-beta_parameters" ][ "gamma" ] ){

                _config[ "coupling_initialization" ][ "update_displacement" ][ "Newmark-beta_parameters" ][ "gamma" ] = 0.5;
                _newmarkGamma = 0.5;

            }

            if ( !_config[ "coupling_initialization" ][ "update_displacement" ][ "Newmark-beta_parameters" ][ "beta" ] ){

                _config[ "coupling_initialization" ][ "update_displacement" ][ "Newmark-beta_parameters" ][ "beta" ] = 0.25;
                _newmarkBeta = 0.25;

            }

        }
        else{

            _config[ "coupling_initialization" ][ "update_displacement" ] = false;

        }

        if ( _config[ "coupling_initialization" ][ "output_reference_information" ] ){

            if ( !_config[ "coupling_initialization" ][ "output_reference_information" ][ "filename" ] ){

                _config[ "coupling_initialization" ][ "output_reference_information" ][ "filename" ] = "reference_information";

            }
            _outputReferenceInformation = true;

            //Check if the file is writable
            std::string tmp = _config[ "coupling_initialization" ][ "output_reference_information" ][ "filename" ].as< std::string >( );
            std::ofstream myfile;
            myfile.open( tmp );

            if ( !myfile.is_open( ) ){

                std::remove( tmp.c_str( ) );
                return new errorNode( "checkCouplingInitialization",
                                      "The 'output_reference_information' filename cannot be opened" );

            }

            std::remove( tmp.c_str( ) );

            _config[ "coupling_initialization" ][ "reference_filename" ] = "reference_information.xdmf";

        }

        if ( _config[ "coupling_initialization" ][ "output_homogenized_response" ] ){

            _outputHomogenizedInformation = true;

            if ( !_config[ "coupling_initialization" ][ "output_homogenized_response" ][ "filename" ] ){

                _config[ "coupling_initialization" ][ "output_homogenized_response" ][ "filename" ] = "homogenized_response";

            }

            std::string tmp = _config[ "coupling_initialization" ][ "output_homogenized_response" ][ "filename" ].as< std::string >( );
            std::ofstream myfile;
            myfile.open( tmp );

            if ( !myfile.is_open( ) ){

                std::remove( tmp.c_str( ) );
                return new errorNode( "checkCouplingInitialization",
                                      "The 'output_homogenized_response' filename cannot be opened" );

            }

            std::remove( tmp.c_str( ) );

            if ( !_config[ "coupling_initialization" ][ "output_homogenized_response" ][ "filetype" ] ){

                _config [ "coupling_initialization" ][ "output_homogenized_response" ][ "filetype" ] = "XDMF";

            }

            if ( !_config[ "coupling_initialization" ][ "output_homogenized_response" ][ "mode" ] ){

                _config [ "coupling_initialization" ][ "output_homogenized_response" ][ "mode" ] = "write";

            }

        }
        else{

            _config[ "coupling_initialization" ][ "output_homogenized_response" ] = false;
            _outputHomogenizedInformation = false;

        }

        if ( _config[ "coupling_initialization" ][ "output_updated_dof" ] ){

            if ( _config[ "coupling_initialization" ][ "update_displacement" ].IsScalar( ) &&
                 !_config[ "coupling_initialization" ][ "update_displacement" ].as< bool >( ) ){

                return new errorNode( "checkCouplingInitialization",
                                      "If 'output_updated_dof' is enabled, then 'update_displacement' must be as well" );

            }

            _outputUpdatedDOF = true;

            if ( !_config[ "coupling_initialization" ][ "output_updated_dof" ][ "macroscale_filename" ] ){

                _config[ "coupling_initialization" ][ "output_updated_dof" ][ "macroscale_filename" ] = "macroscale_dof";

            }

            if ( !_config[ "coupling_initialization" ][ "output_updated_dof" ][ "microscale_filename" ] ){

                _config[ "coupling_initialization" ][ "output_updated_dof" ][ "microscale_filename" ] = "microscale_dof";

            }

            std::string tmp = _config[ "coupling_initialization" ][ "output_updated_dof" ][ "macroscale_filename" ].as< std::string >( );
            std::ofstream myfile;
            myfile.open( tmp );

            if ( !myfile.is_open( ) ){

                std::remove( tmp.c_str( ) );
                return new errorNode( "checkCouplingInitialization",
                                      "The 'output_updated_dof' macroscale filename cannot be opened" );

            }

            std::remove( tmp.c_str( ) );

            tmp = _config[ "coupling_initialization" ][ "output_updated_dof" ][ "microscale_filename" ].as< std::string >( );
            myfile.open( tmp );

            if ( !myfile.is_open( ) ){

                std::remove( tmp.c_str( ) );
                return new errorNode( "checkCouplingInitialization",
                                      "The 'output_updated_dof' microscale filename cannot be opened" );

            }

            std::remove( tmp.c_str( ) );

            if ( !_config[ "coupling_initialization" ][ "output_updated_dof" ][ "macroscale_filetype" ] ){

                _config[ "coupling_initialization" ][ "output_updated_dof" ][ "macroscale_filetype" ] = "XDMF";

            }

            if ( !_config[ "coupling_initialization" ][ "output_updated_dof" ][ "microscale_filetype" ] ){

                _config[ "coupling_initialization" ][ "output_updated_dof" ][ "microscale_filetype" ] = "XDMF";

            }

        }
        else{

            _config[ "coupling_initialization" ][ "output_updated_dof" ] = false;

        }

        if ( !_config[ "coupling_initialization" ][ "solve_coupling_odes_at_microdomains" ] ){

            _config[ "coupling_initialization" ][ "solve_coupling_odes_at_microdomains" ] = false; //TODO: investigate this

        }
        else{

            try{
                _config[ "coupling_initialization" ][ "solve_coupling_odes_at_microdomains" ].as< bool >( );
            }
            catch( ... ){
                return new errorNode( "checkCouplingInitialization", "'solve_coupling_odes_at_micronodes' must be a boolean" ); 
            }

            if ( _config[ "coupling_initialization" ][ "solve_coupling_odes_at_microdomains" ].as< bool >( ) ){

                return new errorNode( "checkCouplingInitialization", "'solve_coupling_odes_at_micronodes' has not been implemented" );

            }


            if ( _config[ "coupling_initialization" ][ "solve_coupling_odes_at_microdomains" ].as< bool >( ) ){

                if ( _config[ "coupling_initialization" ][ "projection_type" ].as< std::string >( ).compare( "averaged_l2_projection" ) != 0 ){
                    return new errorNode( "checkCouplingInitialization", "'solve_coupling_odes_at_micronodes' can only be used with 'averaged_l2_projection'" );
                }

            }

        }
        _couplingODESolutionLocationFlag = _config[ "coupling_initialization" ][ "solve_coupling_odes_at_microdomains" ].as< bool >( );

        if ( !_config[ "coupling_initialization" ][ "use_reconstructed_volume_for_mass_matrix" ] ){

            _config[ "coupling_initialization" ][ "use_reconstructed_volume_for_mass_matrix" ] = true;
            _useReconstructedVolumeForMassMatrix = true;

        }

        _useReconstructedVolumeForMassMatrix
            = _config[ "coupling_initialization" ][ "use_reconstructed_volume_for_mass_matrix" ].as< bool >( );

        _isFiltering = false;
        if ( _config[ "coupling_initialization" ][ "apply_micro_to_macro_filter" ] ){

            _isFiltering = _config[ "coupling_initialization" ][ "apply_micro_to_macro_filter" ].as< bool >( );
            _ghostMacroMassPropertiesRequired = false;

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

            dataFile->getNumSubDomainNodes( increment, *domain, n );
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

            error = dataFile->getSubDomainNodes( increment, *domain, nodes );

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
//        errorOut error = getUniqueNodesInDomains( increment, _microscale, _non_overlapped_micro_surface_sets,
//                                                  _unique_non_overlapped_micro_nodes );
//
//        if ( error ){
//
//            errorOut result = new errorNode( "setMicroNodeIndexMappings",
//                                             "Error in determining the unique non-overlapped microscale nodes" );
//            result->addNext( error );
//            return result;
//
//        }

        errorOut error = getUniqueNodesInDomains( increment, _microscale, _free_micro_volume_sets, _unique_free_micro_nodes );

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

//        //Remove nodes found in the non-overlapped nodes from the free and ghost nodes
//        unsigned int numNodes = 0;
//        unsigned int n;
//
//        //Approximate the size of the duplicate nodes. At worst, this will be the size of the
//        //nodes on the surfaces of the non-overlapped domains
//        for ( auto domain =  _non_overlapped_micro_surface_sets.begin( );
//                   domain != _non_overlapped_micro_surface_sets.end( );
//                   domain++ ){
//
//            _microscale->getNumSubDomainNodes( increment, *domain, n );
//            numNodes += n;
//
//        }
//
//        //Loop through the free nodes to find duplicates
//        uIntVector duplicateNodes;
//        duplicateNodes.reserve( numNodes );
//        
//        for ( auto node =  _unique_free_micro_nodes.begin( );
//                   node != _unique_free_micro_nodes.end( );
//                   node++ ){
//
//            //If the free node is found in the non-overlapped nodes add it to the duplicates
//            if ( std::find( _unique_non_overlapped_micro_nodes.begin( ), _unique_non_overlapped_micro_nodes.end( ),  *node )
//                 != _unique_non_overlapped_micro_nodes.end( ) ){
//
//                duplicateNodes.push_back( node - _unique_free_micro_nodes.begin( ) );
//
//            }
//
//        }
//
//        //Sort the duplicate nodes
//        std::sort( duplicateNodes.begin( ), duplicateNodes.end( ) );
//
//        //Remove the duplicate nodes
//        error = removeIndicesFromVector( _unique_free_micro_nodes, duplicateNodes.begin( ), duplicateNodes.end( ) );
//
//        if ( error ){
//
//            errorOut result = new errorNode( "setMicroNodeIndexMappings",
//                                             "Error in the removal of the non-overlapped duplicate values from the free micro-node vector" );
//            result->addNext( error );
//            return result;
//
//        }
//
//        //Loop through the ghost nodes to find duplicates
//        duplicateNodes.clear( );
//        duplicateNodes.reserve( numNodes );
//        
//        for ( auto node =  _unique_ghost_micro_nodes.begin( );
//                   node != _unique_ghost_micro_nodes.end( );
//                   node++ ){
//
//            //If the ghost node is found in the free nodes add it to the duplicates
//            if ( std::find( _unique_non_overlapped_micro_nodes.begin( ), _unique_non_overlapped_micro_nodes.end( ),  *node )
//                 != _unique_non_overlapped_micro_nodes.end( ) ){
//
//                duplicateNodes.push_back( node - _unique_ghost_micro_nodes.begin( ) );
//
//            }
//
//        }
//
//        //Sort the duplicate nodes
//        std::sort( duplicateNodes.begin( ), duplicateNodes.end( ) );
//
//        //Remove the duplicate nodes
//        error = removeIndicesFromVector( _unique_ghost_micro_nodes, duplicateNodes.begin( ), duplicateNodes.end( ) );
//
//        if ( error ){
//
//            errorOut result = new errorNode( "setMicroNodeIndexMappings",
//                                             "Error in the removal of non-overlapped duplicate values from the ghost vector" );
//            result->addNext( error );
//            return result;
//
//        }

        //Find any duplicated between the free nodes and the ghost nodes
        uIntVector duplicateNodes;
        duplicateNodes.clear( );
        duplicateNodes.reserve( _unique_ghost_micro_nodes.size( ) );
        
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

            _macroscale->getNumSubDomainNodes( increment, *domain, n );
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

    errorOut inputFileProcessor::extractMicroTime( const unsigned int &increment ){
        /*!
         * Extract the timestamp of the micro increment
         *
         * :param const unsigned int &increment: The increment at which to extract the time
         */

        errorOut error = _microscale->getIncrementTime( increment, _microTime );

        if ( error ){
            errorOut result = new errorNode( "extractMicroTime",
                                             "Error in the extraction of the micro domain's time at increment " +
                                             std::to_string( increment ) );
            result->addNext( error );
            return result;
        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMicroTime( const unsigned int &increment, floatType &microTime ){
        /*!
         * Extract the timestamp of the micro increment
         *
         * :param const unsigned int &increment: The increment at which to extract the time
         * :param floatType m&icroTime: The variable to store the micro time in
         */

        errorOut error = _microscale->getIncrementTime( increment, microTime );

        if ( error ){
            errorOut result = new errorNode( "extractMicroTime",
                                             "Error in the extraction of the micro domain's time at increment " +
                                             std::to_string( increment ) );
            result->addNext( error );
            return result;
        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroTime( const unsigned int &increment ){
        /*!
         * Extract the timestamp of the macro increment
         *
         * :param const unsigned int &increment: The increment at which to extract the time
         */

        errorOut error = _macroscale->getIncrementTime( increment, _macroTime );

        if ( error ){
            errorOut result = new errorNode( "extractMacroTime",
                                             "Error in the extraction of the macro domain's time at increment " +
                                             std::to_string( increment ) );
            result->addNext( error );
            return result;
        }

        return NULL;

    }

    errorOut inputFileProcessor::extractMacroTime( const unsigned int &increment, floatType &macroTime ){
        /*!
         * Extract the timestamp of the macro increment
         *
         * :param const unsigned int &increment: The increment at which to extract the time
         * :param floatType &macroTime: The variable to store the macro time in
         */

        errorOut error = _macroscale->getIncrementTime( increment, macroTime );

        if ( error ){
            errorOut result = new errorNode( "extractMacroTime",
                                             "Error in the extraction of the macro domain's time at increment " +
                                             std::to_string( increment ) );
            result->addNext( error );
            return result;
        }

        return NULL;

    }

    errorOut inputFileProcessor::checkVolumeReconstructionInitialization( ){
        /*!
         * Check the initialization of the volume reconstruction
         */

        if ( !_config[ "volume_reconstruction" ] ){

            _config[ "volume_reconstruction" ][ "type" ] = "dual_contouring";

        }

        if ( !_config[ "volume_reconstruction" ][ "element_contain_tolerance" ] ){

            _config[ "volume_reconstruction" ][ "element_contain_tolerance" ] = 1e-2;

        }

        if ( !_config[ "volume_reconstruction" ][ "use_macro_normals" ] ){

            _config[ "volume_reconstruction" ][ "use_macro_normals" ] = true;

        }

        _volumeReconstructionConfig = _config[ "volume_reconstruction" ];

        return NULL;
    }

    errorOut inputFileProcessor::setMicroNodeOutputIndexMappings( const unsigned int &increment ){
        /*!
         * Set the map between the micro node ID number and it's index in the output file.
         *
         * :param const unsigned int &increment: The increment to do the extraction
         */

        uIntVector nodeIDs;
        std::string attributeName = _config[ "microscale_definition" ][ "node_id_variable_name" ].as< std::string >( );
        errorOut error = _microscale->getNodeIds( increment, attributeName, nodeIDs );

        if ( error ){

            errorOut result = new errorNode( "_setMicroNodeOutputIndexMappings", "Error when getting the node ids" );
            result->addNext( error );
            return result;

        }

        _microGlobalNodeIDOutputIndex.reserve( _unique_free_micro_nodes.size( ) + _unique_ghost_micro_nodes.size( ) );

        for ( auto n = _unique_free_micro_nodes.begin( ); n != _unique_free_micro_nodes.end( ); n++ ){

            auto it = std::find( nodeIDs.begin( ), nodeIDs.end( ), *n );

            if ( it == nodeIDs.end( ) ){

                std::string outstr = "Free micro node ";
                outstr += *n;
                outstr += " not found in nodeIds";

                return new errorNode( "_setMicroNodeOutputIndexMappings", outstr );

            }

            uIntType index = it - nodeIDs.begin( );

            _microGlobalNodeIDOutputIndex.emplace( *n, index );

        } 

        for ( auto n = _unique_ghost_micro_nodes.begin( ); n != _unique_ghost_micro_nodes.end( ); n++ ){

            auto it = std::find( nodeIDs.begin( ), nodeIDs.end( ), *n );

            if ( it == nodeIDs.end( ) ){

                std::string outstr = "Ghost micro node ";
                outstr += *n;
                outstr += " not found in nodeIds";

                return new errorNode( "_setMicroNodeOutputIndexMappings", outstr );

            }

            uIntType index = it - nodeIDs.begin( );

            _microGlobalNodeIDOutputIndex.emplace( *n, index );

        } 

        return NULL;

    }

    errorOut inputFileProcessor::setMacroNodeOutputIndexMappings( const unsigned int &increment ){
        /*!
         * Set the map between the macro node ID number and it's index in the output file.
         *
         * :param const unsigned int &increment: The increment to do the extraction
         */

        uIntVector nodeIDs;
        if ( !_config[ "macroscale_definition" ][ "node_id_variable_name" ] ){

            return new errorNode( "setMacroNodeOutputIndexMappings", "'node_id_variable_name' not defined in the input file" );

        }
        std::string attributeName = _config[ "macroscale_definition" ][ "node_id_variable_name" ].as< std::string >( );
        errorOut error = _macroscale->getNodeIds( increment, attributeName, nodeIDs );

        if ( error ){

            errorOut result = new errorNode( "setMacroNodeOutputIndexMappings", "Error when getting the node ids" );
            result->addNext( error );
            return result;

        }

        _macroGlobalNodeIDOutputIndex.reserve( _unique_free_macro_nodes.size( ) + _unique_ghost_macro_nodes.size( ) );

        for ( auto n = _unique_free_macro_nodes.begin( ); n != _unique_free_macro_nodes.end( ); n++ ){

            auto it = std::find( nodeIDs.begin( ), nodeIDs.end( ), *n );

            if ( it == nodeIDs.end( ) ){

                std::string outstr = "Free macro node ";
                outstr += *n;
                outstr += " not found in nodeIds";

                return new errorNode( "setMacroNodeOutputIndexMappings", outstr );

            }

            uIntType index = it - nodeIDs.begin( );

            _macroGlobalNodeIDOutputIndex.emplace( *n, index );

        } 

        for ( auto n = _unique_ghost_macro_nodes.begin( ); n != _unique_ghost_macro_nodes.end( ); n++ ){

            auto it = std::find( nodeIDs.begin( ), nodeIDs.end( ), *n );

            if ( it == nodeIDs.end( ) ){

                std::string outstr = "Ghost macro node ";
                outstr += *n;
                outstr += " not found in nodeIds";

                return new errorNode( "setMacroNodeOutputIndexMappings", outstr );

            }

            uIntType index = it - nodeIDs.begin( );

            _macroGlobalNodeIDOutputIndex.emplace( *n, index );

        } 

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

    bool inputFileProcessor::microSurfaceForceDefined( ){
        /*!
         * Get whether the micro-surface force has been defined
         */

        return _microSurfaceForceFlag;
    }

    bool inputFileProcessor::microExternalForceDefined( ){
        /*!
         * Get whether the micro external force has been defined
         */

        return _microExternalForceFlag;
    }

    bool inputFileProcessor::microInternalForceDefined( ){
        /*!
         * Get whether the micro internal force has been defined
         */

        return _microInternalForceFlag;
    }

    bool inputFileProcessor::microInertialForceDefined( ){
        /*!
         * Get whether the micro inertial force has been defined
         */

        return _microInertialForceFlag;
    }

    bool inputFileProcessor::microVelocitiesDefined( ){
        /*!
         * Get whether the micro-velocities has been defined
         */

        return _microVelocityFlag;
    }

    bool inputFileProcessor::microAccelerationDefined( ){
        /*!
         * Get whether the micro-acceleration has been defined
         */

        return _microAccelerationFlag;
    }

    bool inputFileProcessor::macroVelocitiesDefined( ){
        /*!
         * Get whether the macro-velocities has been defined
         */

        return _macroVelocityFlag;
    }

    bool inputFileProcessor::macroAccelerationDefined( ){
        /*!
         * Get whether the macro-acceleration has been defined
         */

        return _macroAccelerationFlag;
    }

    bool inputFileProcessor::macroInternalForceDefined( ){
        /*!
         * Get whether the macro internal force has been defined
         */

        return _macroInternalForceFlag;
    }

    bool inputFileProcessor::macroExternalForceDefined( ){
        /*!
         * Get whether the macro external force has been defined
         */

        return _macroExternalForceFlag;
    }

    bool inputFileProcessor::macroInertialForceDefined( ){
        /*!
         * Get whether the macro inertial force has been defined
         */

        return _macroInertialForceFlag;
    }

    bool inputFileProcessor::macroLumpedMassMatrixDefined( ){
        /*!
         * Get whether the macro lumped mass matrix has been defined
         */

        return _macroLumpedMassMatrixFlag;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMicroDisplacements( ){
        /*!
         * Get the micro-displacements
         */

        return &_microDisplacements;
    }

    const floatType* inputFileProcessor::getMacroTime( ){
        /*!
         * Get the timestamp of the macro increment
         */

        return &_macroTime;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMacroDisplacements( ){
        /*!
         * Get the macro-displacements
         */

        return &_macroDisplacements;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMacroDispDOFVector( ){
        /*!
         * Get the macro-displacement DOF vector
         */

        return &_macroDispDOFVector;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMacroVelocities( ){
        /*!
         * Get a pointer to the macro velocities
         */

        return &_macroVelocities;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMacroAccelerations( ){
        /*!
         * Get a pointer to the macro accelerations
         */

        return &_macroAccelerations;

    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getPreviousMacroDispDOFVector( ){
        /*!
         * Get a pointer to the previous macro values of the displacement degrees of freedom
         *
         * If extraction of those values as not requested, the current macro displacement degrees of freedom will be returned
         */

        if ( _extractPreviousDOFValues ){

            return &_previousMacroDispDOFVector;

        }
        else{

            return &_macroDispDOFVector;

        }

    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getPreviousMacroVelocities( ){
        /*!
         * Get a pointer to the previous macro velocities
         *
         * If extraction of those values as not requested, the current macro velocities will be returned
         */

        if ( _extractPreviousDOFValues ){

            return &_previousMacroVelocities;

        }
        else{

            return &_macroVelocities;

        }

    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getPreviousMacroAccelerations( ){
        /*!
         * Get a pointer to the previous macro accelerations
         *
         * If extraction of those values as not requested, the current macro accelerations will be returned
         */

        if ( _extractPreviousDOFValues ){

            return &_previousMacroAccelerations;

        }
        else{

            return &_macroAccelerations;

        }

    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMacroInternalForces( ){
        /*!
         * Get a pointer to the macro internal forces
         */

        return &_macroInternalForces;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMacroExternalForces( ){
        /*!
         * Get a pointer to the macro internal forces
         */

        return &_macroExternalForces;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMacroInertialForces( ){
        /*!
         * Get a pointer to the macro inertial forces
         */

        return &_macroInertialForces;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMacroBodyForces( ){
        /*!
         * Get a pointer to the macro body forces
         */

        return &_macroBodyForces;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMacroSurfaceForces( ){
        /*!
         * Get a pointer to the macro surface forces
         */

        return &_macroSurfaceForces;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMicroNodeReferencePositions( ){
        /*!
         * Get the nodal positions from which the displacements are
         * referenced.
         */

        return &_microNodeReferencePositions;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMacroNodeReferencePositions( ){
        /*!
         * Get the nodal positions from which the displacements are
         * referenced.
         */

        return &_macroNodeReferencePositions;
    }

    const std::unordered_map< uIntType, uIntVector >* inputFileProcessor::getMacroNodeReferenceConnectivity( ){
        /*!
         * Get the nodal positions from which the displacements are
         * referenced.
         */

        return &_macroNodeReferenceConnectivity;
    }

//    const uIntVector* inputFileProcessor::getMacroNodeReferenceConnectivityCellIndices( ){
//        /*!
//         * Get the nodal positions from which the displacements are
//         * referenced.
//         */
//
//        return &_macroNodeReferenceConnectivityCellIndices;
//    }

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

//    const uIntVector* inputFileProcessor::getFreeMacroCellMicroDomainCounts( ){
//        /*!
//         * Get the free macro cell ids
//         */
//
//        return &_free_macro_cell_micro_domain_counts;
//    }
//
//    const uIntVector* inputFileProcessor::getGhostMacroCellMicroDomainCounts( ){
//        /*!
//         * Get the ghost macro cell ids
//         */
//
//        return &_ghost_macro_cell_micro_domain_counts;
//    }

    bool inputFileProcessor::computeMicroShapeFunctions( ){
        /*!
         * Return whether the micro-shape functions should be computed
         */

        return _computeMicroShapeFunctions;
    }

//    const uIntVector *inputFileProcessor::getNonOverlappedMicroNodeIds( ){
//        /*!
//         * Get the free micro-node ids
//         */
//
//        return &_unique_non_overlapped_micro_nodes;
//    }

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

    const std::unordered_map< unsigned int, floatVector > *inputFileProcessor::getMacroReferenceDensities( ){
        /*!
         * Get the macro reference densities for the macro domains
         */

        return &_macroReferenceDensities;

    }

    const std::unordered_map< unsigned int, floatVector > *inputFileProcessor::getMacroReferenceMomentsOfInertia( ){
        /*!
         * Get the macro reference moments of inertia
         */

        return &_macroReferenceMomentsOfInertia;

    }

    const std::unordered_map< unsigned int, std::string > *inputFileProcessor::getMacroReferenceDensityTypes( ){
        /*!
         * Get the macro reference density types for the macro domains
         */

        return &_macroReferenceDensityTypes;

    }

    const std::unordered_map< unsigned int, std::string > *inputFileProcessor::getMacroReferenceMomentOfInertiaTypes( ){
        /*!
         * Get the macro reference moment of inertia types
         */

        return &_macroReferenceMomentOfInertiaTypes;

    }

    const floatType* inputFileProcessor::getDt( ){
        /*!
         * Get the change in time for the DOF update
         */

        return &_Dt;
    }

    const floatType* inputFileProcessor::getNewmarkGamma( ){
        /*!
         * Get the gamma parameter for the Newmark beta method
         */

        return &_newmarkGamma;
    }

    const floatType* inputFileProcessor::getNewmarkBeta( ){
        /*!
         * Get the beta parameter for the Newmark beta method
         */

        return &_newmarkBeta;
    }

    bool inputFileProcessor::outputReferenceInformation( ){
        /*!
         * Get a flag for whether the reference information should be output to a file
         */

        return _outputReferenceInformation;

    }

    bool inputFileProcessor::outputHomogenizedInformation( ){
        /*!
         * Get a flag for whether the homogenized information should be output to a file
         */

        return _outputHomogenizedInformation;

    }

    bool inputFileProcessor::outputUpdatedDOF( ){
        /*!
         * Get a flag for whether the updated DOF should be output to a file
         */

        return _outputUpdatedDOF;

    }

    const DOFMap* inputFileProcessor::getMacroNodeIDOutputIndex( ){
        /*!
         * Returns a constant reference to the macro node id to output index map
         */

        return &_macroGlobalNodeIDOutputIndex;
    }

    const DOFMap* inputFileProcessor::getMicroNodeIDOutputIndex( ){
        /*!
         * Returns a constant reference to the micro node id to output index map
         */

        return &_microGlobalNodeIDOutputIndex;

    }

    const std::unordered_map< uIntType, stringVector >* inputFileProcessor::getMacroCellToDomainMap( ){
        /*!
         * Returns a constant reference to the macro cell to micro domain map
         */

        return &_macroCellDomainMap;
    }

    const std::unordered_map< std::string, uIntType >* inputFileProcessor::getMicroDomainIDMap( ){
        /*!
         * Returns a constant reference to the micro domain id map
         */

        return &_microDomainIDMap;
    }

    bool inputFileProcessor::solveCouplingODEsAtMicroDomains( ){
        /*!
         * Return whether the coupling ODEs should be solved at the micro domains or not
         */

        return _couplingODESolutionLocationFlag;
    }

    bool inputFileProcessor::useArlequinCoupling( ){
        /*!
         * Return whether Arlequin coupling should be used
         */

        return _useArlequinMethod;
    }

    const std::unordered_map< uIntType, floatType >* inputFileProcessor::getMacroArlequinWeights( ){
        /*!
         * Get the macro Arlequin weights at the micro nodes
         */

        return &_macroArlequinWeights;
    }

    const std::unordered_map< uIntType, floatVector >* inputFileProcessor::getMacroLumpedMassMatrix( ){
        /*!
         * Get the macro lumped-mass matrix
         */

        return &_macroLumpedMassMatrix;
    }

    const floatType* inputFileProcessor::getArlequinPenaltyParameter( ){
        /*!
         * Get the Arlequin penalty parameter
         */

        return &_arlequinPenaltyParameter;
    }

    const floatType* inputFileProcessor::getArlequinUpdatePenaltyParameter( ){
        /*!
         * Get the Arlequin penalty parameter
         */

        return &_arlequinUpdatePenaltyParameter;
    }

    bool inputFileProcessor::useReconstructedVolumeForMassMatrix( ){
        /*!
         * Use the reconstructed volume properties for the mass matrix rather than the 
         * discrete points.
         */

        return _useReconstructedVolumeForMassMatrix;
    }

    bool inputFileProcessor::isFiltering( ){
        /*!
         * Flag to indicate if the coupling is a filtering operation meaning no
         * macroscale information is returned to the microscale
         */

        return _isFiltering;
    }

}
