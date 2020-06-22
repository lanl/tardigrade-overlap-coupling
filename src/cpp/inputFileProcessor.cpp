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
        std::vector< std::string > setNames;
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

        return NULL;
    }

    errorOut inputFileProcessor::initializeCouplingDomains( ){
        /*!
         * Initialize the coupling domains
         */

        if ( _config[ "free_macroscale_domains" ] ){

            errorOut error = checkCommonDomainConfiguration( _config[ "free_macroscale_domains" ] );
            if ( error ){
                errorOut result = new errorNode( "initializeCouplingDomains",
                                                 "Error in input-file check of the free macroscale domains" );
                result->addNext( error );
                return result;
            }

        }

        if ( _config[ "ghost_macroscale_domains" ] ){

            //In the configuration of the micro-scale domains
            errorOut error = checkCommonDomainConfiguration( _config[ "ghost_macroscale_domains" ] );
            if ( error ){
                errorOut result = new errorNode( "initializeCouplingDomains",
                                                 "Error in input-file check of the ghost macroscale domains" );
                result->addNext( error );
                return result;
            }

        }

        return NULL;

    }

    errorOut inputFileProcessor::checkCommonDomainConfiguration( const YAML::Node &domainConfig ){
        /*!
         * Check the common values in the configuration of a domain
         *
         * :param YAML::Node &domainConfig: The configuration of a particular domain.
         */

        if ( ( !domainConfig.IsSequence() ) && ( !domainConfig.IsNull( ) ) ){

            std::cout << domainConfig[ "macro_free_1" ][ "macro_nodeset" ].as<std::string>( ) <<"\n";

            return new errorNode( "checkCommonDomainConfiguration",
                                  "The definition of the domains must either be empty or a sequence" );

        }

        unsigned int indx = 1;
        unsigned int indx2 = 1;
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

                indx2++;

            }

            indx++;

        }

        return NULL;
    }

}
