/*!============================================================================
|                               overlapCoupling                               |
===============================================================================
| The implementation of the Overlap coupling method                           |
=============================================================================*/

#include<overlapCoupling.h>
#include<Eigen/SparseQR>

namespace overlapCoupling{

    overlapCoupling::overlapCoupling( ){
        /*!
         * The default constructor
         */

        return;
    }

    overlapCoupling::overlapCoupling( const std::string &configurationFilename ){
        /*!
         * The constructor where the configuration filename is provided
         *
         * :param const std::string &configurationFilename: The configuration filename
         */

        errorOut error = setConfigurationFilename( configurationFilename );

        if ( error ){
            _error = new errorNode( "overlapCoupling", "Error when setting the configuration filename" );
            _error->addNext( error );
        }

        return;
    }

    errorOut overlapCoupling::setConfigurationFilename( const std::string &configurationFilename ){
        /*!
         * Set the configuration filename.
         *
         * :param const std::string &configurationFilename: The configuration filename.
         */

        _error = NULL;

        errorOut error = _inputProcessor.setConfigurationFilename( configurationFilename );

        if ( error ){
            errorOut result = new errorNode( "setConfigurationFilename",
                                             "Error in setting the configuration filename of the input processor" );
            result->addNext( error );
            return result;
        }
        return NULL;
    }

    const errorOut overlapCoupling::getConstructorError( ){
        /*!
         * Return the current value of the error during the construction.
         */

        return _error;
    }

    errorOut overlapCoupling::processIncrement( const unsigned int &microIncrement,
                                                const unsigned int &macroIncrement ){
        /*!
         * Process the indicated increment
         *
         * :param const unsigned int &microIncrement: The micro increment to process
         * :param const unsigned int &macroIncrement: The macro increment to process
         */

        //Initialize the input processor
        errorOut error = _inputProcessor.initializeIncrement( microIncrement, macroIncrement );
        if ( error ){
            errorOut result = new errorNode( "processIncrement", "Error in initialization of the input processor" );
            result->addNext( error );
            return result;
        }

        //Compute the centers of mass of the free and ghost domains
        error = computeIncrementCentersOfMass( microIncrement, macroIncrement,
                                               _freeMicroDomainMasses, _ghostMicroDomainMasses,
                                               _freeMicroDomainCentersOfMass, _ghostMicroDomainCentersOfMass );
        if ( error ){
            errorOut result = new errorNode( "processIncrement", "Error in computation of the domain centers of mass" );
            result->addNext( error );
            return result;
        }

        //Project the degrees of freedom
        error = projectDegreesOfFreedom( );

        if ( error ){
            errorOut result = new errorNode( "processIncrement", "Error in the projection of the ghost degrees of freedom" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut overlapCoupling::initializeCoupling( ){
        /*!
         * Initialize the coupling between the domains
         *
         * Configuration for this process is located in the YAML file under
         * the root level key "coupling_initialization". If this is not 
         * defined, a default strategy will be employed. This strategy will
         * be written out to configurationFilename.as_evaluated.
         */

        //Get the coupling initialization from the configuration file
        const YAML::Node couplingInitialization = _inputProcessor.getCouplingInitialization( );

        if ( !couplingInitialization ){
            return new errorNode ( "initializeCoupling", "The coupling initialization configuration is not defined" );
        }

        errorOut error;
        if ( couplingInitialization[ "type" ].as< std::string >( ).compare( "use_first_increment" ) == 0 ){
            error = setReferenceStateFromIncrement( 0, 0 );
        }
        else{
            return new errorNode( "initializeCoupling",
                                  "The coupling initialization type '" + couplingInitialization[ "type" ].as< std::string >( )
                                + "' is not recognized" );
        }
        
        if ( error ){
            errorOut result = new errorNode( "initializeCoupling", "Error in initialization of the coupling" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut overlapCoupling::setReferenceStateFromIncrement( const unsigned int &microIncrement,
                                                              const unsigned int &macroIncrement ){
        /*!
         * Set the reference state from the indicated increment
         *
         * :param const unsigned int &microIncrement: The micro increment at which to set the reference state
         * :param const unsigned int &macroIncrement: The macro increment at which to set the reference state
         */

        //Initialize the input processor
        errorOut error = _inputProcessor.initializeIncrement( microIncrement, macroIncrement );
        if ( error ){
            errorOut result = new errorNode( "processIncrement", "Error in initialization of the input processor" );
            result->addNext( error );
            return result;
        }

        //Get the macro cell ids
        const uIntVector *freeMacroCellIDs = _inputProcessor.getFreeMacroCellIds( );
        const uIntVector *ghostMacroCellIDs = _inputProcessor.getGhostMacroCellIds( );

        //Get the macro domain names
        const stringVector *freeMacroDomainNames = _inputProcessor.getFreeMacroDomainNames( );
        const stringVector *ghostMacroDomainNames = _inputProcessor.getGhostMacroDomainNames( );

        //Get the micro cell names
        const stringVector *freeMicroDomainNames = _inputProcessor.getFreeMicroDomainNames( );
        const stringVector *ghostMicroDomainNames = _inputProcessor.getGhostMicroDomainNames( );

        //Get the macro cell counts
        const uIntVector *freeMacroCellMicroDomainCounts = _inputProcessor.getFreeMacroCellMicroDomainCounts( );
        const uIntVector *ghostMacroCellMicroDomainCounts = _inputProcessor.getGhostMacroCellMicroDomainCounts( );

        //Set the output vector sizes
        unsigned int numFreeMicroDomains  = _inputProcessor.getFreeMicroDomainNames( )->size( );
        unsigned int numGhostMicroDomains = _inputProcessor.getGhostMicroDomainNames( )->size( );

        //Set the reference micro domain mass vector sizes
        _referenceFreeMicroDomainMasses = floatVector( numFreeMicroDomains, 0 );
        _referenceGhostMicroDomainMasses = floatVector( numGhostMicroDomains, 0 );

        //Set the reference micro domain center of mass vector sizes
        _referenceFreeMicroDomainCentersOfMass = floatVector( _dim * numFreeMicroDomains, 0 );
        _referenceGhostMicroDomainCentersOfMass = floatVector( _dim * numGhostMicroDomains, 0 );

        //Loop over the free macro-scale cells
        unsigned int cellIndex;
        unsigned int nMicroDomains;
        unsigned int domainIndex = 0;

        uIntVector macroNodes;

        floatVector domainReferenceXiVectors;
        floatVector domainCenterOfMassShapeFunctionValues;
        floatVector domainMicroPositionShapeFunctionValues;

        for ( auto cellID  = freeMacroCellIDs->begin( );
                   cellID != freeMacroCellIDs->end( );
                   cellID++ ){

                //Set the index
                cellIndex = cellID - freeMacroCellIDs->begin( );

                //Set the number of micro-domains encompassed by the cell
                nMicroDomains = ( *freeMacroCellMicroDomainCounts )[ cellIndex ];

                //Get the macro-node set
                error = _inputProcessor._macroscale->getDomainNodes( macroIncrement, ( *freeMacroDomainNames )[ cellIndex ], macroNodes );

                if ( error ){

                    errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                                     "Error in extracting the free macro-node set" );
                    result->addNext( error );
                    return result;

                }

                //Loop over the micro-domains
                for ( auto domain  = ghostMicroDomainNames->begin( ) + domainIndex;
                           domain != ghostMicroDomainNames->begin( ) + domainIndex + nMicroDomains;
                           domain++ ){

                    error = processDomainReference( microIncrement,
                                                    domain - ghostMicroDomainNames->begin( ), *domain,
                                                    *cellID, macroNodes,
                                                    _referenceGhostMicroDomainMasses[ domain - ghostMicroDomainNames->begin( ) ],
                                                    _referenceGhostMicroDomainCentersOfMass,
                                                    domainReferenceXiVectors,
                                                    domainCenterOfMassShapeFunctionValues,
                                                    domainMicroPositionShapeFunctionValues );

                    if ( error ){

                        errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                                         "Error in processing '" + *domain + "' for a free reference state" );
                        result->addNext( error );
                        return result;

                    }

                }

                //TODO: Add direct projection projector construction here

                domainIndex += nMicroDomains;

        }

        //Loop over the ghost macro-scale cells

        domainIndex = 0;

        for ( auto cellID  = ghostMacroCellIDs->begin( );
                   cellID != ghostMacroCellIDs->end( );
                   cellID++ ){

                //Set the index
                cellIndex = cellID - ghostMacroCellIDs->begin( );

                //Set the number of micro-domains encompassed by the cell
                nMicroDomains = ( *ghostMacroCellMicroDomainCounts )[ cellIndex ];

                //Get the macro-node set
                error = _inputProcessor._macroscale->getDomainNodes( macroIncrement, ( *ghostMacroDomainNames )[ cellIndex ], macroNodes );

                if ( error ){

                    errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                                     "Error in extracting the ghost macro-node set" );
                    result->addNext( error );
                    return result;

                }

                //Loop over the micro-domains
                for ( auto domain  = freeMicroDomainNames->begin( ) + domainIndex;
                           domain != freeMicroDomainNames->begin( ) + domainIndex + nMicroDomains;
                           domain++ ){

                    error = processDomainReference( microIncrement,
                                                    domain - freeMicroDomainNames->begin( ), *domain,
                                                    *cellID, macroNodes,
                                                    _referenceFreeMicroDomainMasses[ domain - freeMicroDomainNames->begin( ) ],
                                                    _referenceFreeMicroDomainCentersOfMass,
                                                    domainReferenceXiVectors,
                                                    domainCenterOfMassShapeFunctionValues,
                                                    domainMicroPositionShapeFunctionValues );

                    if ( error ){

                        errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                                         "Error in processing '" + *domain + "' for a ghost reference state" );
                        result->addNext( error );
                        return result;

                    }

                }

                //TODO: Add direct projection projector construction here

                domainIndex += nMicroDomains;

        }

        //Form the projectors
        error = formTheProjectors( );

        if ( error ){

            errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                             "Error in the formation of the projectors" );
            result->addNext( error );
            return result;

        }

        return NULL;
    }

    errorOut overlapCoupling::formTheProjectors( ){
        /*!
         * Form the projection operators
         */

        const YAML::Node config = _inputProcessor.getCouplingInitialization( );

        if ( config[ "projection_type" ].as< std::string >( ).compare( "l2_projection" ) == 0 ){

            errorOut error = formL2Projectors( );

            if ( error ){

                errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                                 "Error in the formation of the L2 projectors" );
                result->addNext( error );
                return result;

            }

        }
        else if ( config[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){

            errorOut error = formDirectProjectionProjectors( );

            if ( error ){

                errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                                 "Error in the formation of the direct projection projectors" );
                result->addNext( error );
                return result;

            }

        }
        else{

            return new errorNode( "formTheProjectors",
                                  "'projection_type' '" + config[ "projection_type" ].as< std::string >( ) + "' not recognized" );

        }

        return NULL;
    }

    errorOut overlapCoupling::formL2Projectors( ){
        /*!
         * Form the projectors if the L2 projection is to be used
         */
            
        //Set the dimension of the displacement DOF
        unsigned int nDispMicroDOF = _dim;

        unsigned int nDispMacroDOF = _dim + _dim * _dim;

        //Get the number of micro degrees of freedom which are free and ghost
        unsigned int nFreeMicroDOF = nDispMicroDOF * _inputProcessor.getFreeMicroNodeIds( )->size( );
        unsigned int nGhostMicroDOF = nDispMicroDOF * _inputProcessor.getGhostMicroNodeIds( )->size( );

        //Get the number of macro degrees of freedom which are free and ghost
        unsigned int nFreeMacroDOF = nDispMacroDOF * _inputProcessor.getFreeMacroNodeIds( )->size( );
        unsigned int nGhostMacroDOF = nDispMacroDOF * _inputProcessor.getGhostMacroNodeIds( )->size( );

        //Extract the part of the shapefunction matrix that interpolates between ghost micromorphic DOF and free classical DOF
        Eigen::SparseQR< SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
        std::cout << "Performing QR decomposition of NQDhat\n";
        solver.compute( _N.block( 0, nFreeMacroDOF, nFreeMicroDOF, nGhostMacroDOF ) );

        if ( solver.info( ) != Eigen::Success ){

            return new errorNode( "formL2Projectors",
                                  "The QR decomposition of the ghost macro to free micro interpolation matrix failed" );

        }

        //Form the identity matrix
        SparseMatrix I( nFreeMicroDOF, nFreeMicroDOF );
        I.reserve( nFreeMicroDOF );
        for ( unsigned int i = 0; i < nFreeMicroDOF; i++ ){
            I.insert( i, i ) = 1;
        }

        std::cout << "Performing linear solve for BDhatQ\n";
        _L2_BDhatQ = solver.solve( I.toDense( ) );
        _L2_BDhatD = -_L2_BDhatQ * _N.block( 0, 0, nFreeMicroDOF, nFreeMacroDOF );

        _L2_BQhatQ = _N.block( nFreeMicroDOF, nFreeMacroDOF, nGhostMicroDOF, nGhostMacroDOF ) * _L2_BDhatQ;
        _L2_BQhatD = _N.block( nFreeMicroDOF, 0, nGhostMicroDOF, nFreeMacroDOF )
                   + _N.block( nFreeMicroDOF, nFreeMacroDOF, nGhostMicroDOF, nGhostMacroDOF ) * _L2_BDhatD;

        std::cout << "Exiting L2 projector formation\n";

        return NULL;
    }

    errorOut overlapCoupling::formDirectProjectionProjectors( ){
        /*!
         * Form the projectors if the direct projection is to be used
         */

        return NULL;
    }

    errorOut overlapCoupling::processDomainReference( const unsigned int &microIncrement,
                                                      const unsigned int &domainIndex, const std::string &domainName,
                                                      const unsigned int cellID, const uIntVector &macroNodes,
                                                      floatType   &referenceMicroDomainMass,
                                                      floatVector &referenceMicroDomainCentersOfMass,
                                                      floatVector &domainReferenceXiVectors,
                                                      floatVector &domainCenterOfMassShapeFunctionValues,
                                                      floatVector &domainMicroPositionShapeFunctionValues ){
        /*!
         * Process the domain for use with preparing the reference configuration
         *
         * :param const unsigned int &microIncrement: The micro-scale increment at which to compute the reference
         * :param const unsigned int &domainIndex: The index of the domain
         * :param const std::string &domainName: The name of the domain
         * :param const unsigned int cellID: The global cell ID number
         * :param const uIntVector &macroNodes: The nodes of the macro domain
         * :param floatType   &referenceMicroDomainMass: The reference mass of the micro domain. This
         *     should be a reference to the micro-domain mass vector
         * :param floatVector &referenceMicroDomainCentersOfMass: The reference micro-domain center of mass
         *     vector for all of the micro domains
         * :param floatVector &domainReferenceXiVectors: The reference Xi vectors for the domain.
         * :param floatVector &domainCenterOfMassShapeFunctionValues: The shape function values at the 
         *     center of mass
         * :param floatVector &domainMicroPositionShapeFunctionValues: The shape function values at the
         *     micro node positions inside of the domain. This is only computed if indicated by the 
         *     configuration file.
         */

        //Process the domain mass data
        floatVector domainCenterOfMass;

        errorOut error = processDomainMassData( microIncrement, domainName, referenceMicroDomainMass,
                                                domainCenterOfMass, domainReferenceXiVectors );

        if ( error ){
            
            errorOut result = new errorNode( "processDomainReference",
                                             "Error in processing the mass data for the micro domain '" + domainName + "'" );
            result->addNext( error );
            return result;

        }

        //Save the center of mass
        for ( unsigned int i = 0; i < domainCenterOfMass.size( ); i++ ){

            referenceMicroDomainCentersOfMass[ _dim * domainIndex + i ] = domainCenterOfMass[ i ];

        }

        //Compute the domain's shape function information
        error = computeDomainShapeFunctionInformation( cellID, domainName, microIncrement, domainCenterOfMass,
                                                       domainCenterOfMassShapeFunctionValues,
                                                       domainMicroPositionShapeFunctionValues );

        if ( error ){
            errorOut result = new errorNode( "processDomainReference",
                                             "Error in computing the shape function values for the domain center of mass or the domainMicroPositionShapeFunctionValues" );
            result->addNext( error );
            return result;
        }

        //Get the domain node ids
        uIntVector domainNodes;
        error = _inputProcessor._microscale->getDomainNodes( microIncrement, domainName, domainNodes );

        if ( error ){
            errorOut result = new errorNode( "processDomainReference",
                                             "Error in extracting the micro-node set" );
            result->addNext( error );
            return result;
        }

        //Add the domain's contribution to the shape function matrix
        error = addDomainContributionToInterpolationMatrix( domainNodes, macroNodes, domainReferenceXiVectors,
                                                            domainCenterOfMassShapeFunctionValues );

        if ( error ){
            errorOut result = new errorNode( "processDomainReference",
                                             "Error in adding part of the shapefunction matrix determined from '" + domainName + "'" );
            result->addNext( error );
            return result;
        }

        return NULL;

    }

    errorOut overlapCoupling::processDomainMassData( const unsigned int &microIncrement, const std::string &domainName,
                                                     floatType &domainMass, floatVector &domainCenterOfMass,
                                                     floatVector &domainXiVectors ){
        /*!
         * Process a micro-scale domain
         *
         * :param const unsigned int microIncrement: The micro increment to process
         * :param const std::string &domainName: The name of the domain
         * :param floatType &domainMass: The mass of the domain
         * :param floatVector &domainCenterOfMass
         * :param floatVector &domainXiVectors
         */

        //Get the domain's nodes
        uIntVector domainNodes;

        errorOut error = _inputProcessor._microscale->getDomainNodes( microIncrement, domainName, domainNodes );

        if ( error ){

            errorOut result = new errorNode( "processDomain",
                                             "Error in getting the nodes from the micro domain '" + domainName + "'" );
            result->addNext( error );
            return result;

        }

        //Compute the center of mass of the domain
        error = DOFProjection::computeDomainCenterOfMass( _dim, domainNodes, *_inputProcessor.getMicroVolumes( ),
                                                          *_inputProcessor.getMicroDensities( ),
                                                          *_inputProcessor.getMicroNodeReferencePositions( ),
                                                          *_inputProcessor.getMicroDisplacements( ),
                                                          *_inputProcessor.getMicroWeights( ),
                                                          domainMass, domainCenterOfMass );

        if ( error ){

            errorOut result = new errorNode( "processDomain", "Error in calculation of '" + domainName + "' center of mass" );
            result->addNext( error );
            return result;

        }

        //Compute the relative position vectors
        error = DOFProjection::computeDomainXis( _dim, domainNodes,
                                                 *_inputProcessor.getMicroNodeReferencePositions( ),
                                                 *_inputProcessor.getMicroDisplacements( ),
                                                 domainCenterOfMass, domainXiVectors );

        if ( error ){
            
            errorOut result = new errorNode( "processDomain", "Error in calculation of '" + domainName + "' xi vectors" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut overlapCoupling::computeDomainShapeFunctionInformation( const unsigned int &cellID,
                                                                     const std::string &domainName,
                                                                     const unsigned int &microIncrement,
                                                                     const floatVector &domainCenterOfMass,
                                                                     floatVector &domainCenterOfMassShapeFunctionValues,
                                                                     floatVector &domainMicroPositionShapeFunctionValues ){
        /*!
         * Compute the shape function values at the required locations
         *
         * :param const unsigned int &cellID: The ID number of the cell which contains the domain
         * :param const std::string &domainName: The name of the nodeset which defines the micro-scale domain
         * :param const unsigned int &microIncrement: The micro-scale increment to analyze
         * :param const floatVector &domainCenterOfMass: The center of mass of the domain
         * :param floatVector &domainCenterOfMassShapeFunctionValues: The shapefunction values at the center of mass
         * :param floatVector &domainMicroPositionShapeFunctionValues: The shapefunction values at all of the micro nodes
         *     contained within the domain.
         */

        //Compute the shape functions of the domain's center of mass
        errorOut error = computeShapeFunctionsAtPoints( cellID,
                                                        *_inputProcessor.getMacroNodeReferencePositions( ),
                                                        *_inputProcessor.getMacroDisplacements( ),
                                                        *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                        *_inputProcessor.getMacroNodeReferenceConnectivityCellIndices( ),
                                                        domainCenterOfMass, domainCenterOfMassShapeFunctionValues );

        if ( error ){

            errorOut result = new errorNode( "computeDomainShapeFunctionInformation",
                                             "Error in the computation of the center of mass for a micro domain" );
            result->addNext( error );
            return result;

        }

        //Get the domain's nodes
        uIntVector domainNodes;

        error = _inputProcessor._microscale->getDomainNodes( microIncrement, domainName, domainNodes );

        if ( error ){

            errorOut result = new errorNode( "computeDomainShapeFunctionInformation",
                                             "Error in the extraction of the nodes in the micro domain" );
            result->addNext( error );
            return result;

        }

        if ( _inputProcessor.computeMicroShapeFunctions( ) ){

            //Get the micro-node positions
            floatVector microNodePositions( _dim * domainNodes.size( ) );
    
            unsigned int index;
            const floatVector *microReferencePositions = _inputProcessor.getMicroNodeReferencePositions( );
            const floatVector *microDisplacements      = _inputProcessor.getMicroDisplacements( );

            for ( auto it = domainNodes.begin( ); it != domainNodes.end( ); it++ ){
    
                index = it - domainNodes.begin( );
                for ( unsigned int i = 0; i < _dim; i++ ){
                    microNodePositions[ index ] = ( *microReferencePositions )[ index + i ] + ( *microDisplacements )[ index + i ];
                }
    
            }
    
            //Compute the shape function values at the micro positions
            error = computeShapeFunctionsAtPoints( cellID,
                                                   *_inputProcessor.getMacroNodeReferencePositions( ),
                                                   *_inputProcessor.getMacroDisplacements( ),
                                                   *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                   *_inputProcessor.getMacroNodeReferenceConnectivityCellIndices( ),
                                                   microNodePositions,
                                                   domainMicroPositionShapeFunctionValues );
    
            if ( error ){
    
                errorOut result = new errorNode( "computeDomainShapeFunctionInformation",
                                                 "Error in the computation of the center of mass for a micro domain" );
                result->addNext( error );
                return result;
    
            }

        }

        return NULL;

    }

    errorOut overlapCoupling::computeIncrementCentersOfMass( const unsigned int microIncrement, const unsigned int macroIncrement,
                                                             floatVector &freeDomainMass, floatVector &ghostDomainMass,
                                                             floatVector &freeDomainCM, floatVector &ghostDomainCM ){
        /*!
         * Compute the centers of mass for an micro increment. Also computes the mass of the micro-scale domains
         *
         * :param const unsigned int microIncrement: The micro increment at which to compute the centers of mass
         * :param const unsigned int macroIncrement: The macro increment at which to compute the centers of mass
         * :param floatVector &freeDomainMass: The mass of the free domains
         * :param floatVector &ghostDomainMass: The mass of the ghost domains
         * :param floatVector &freeDomainCM: The center of mass of the free domains
         * :param floatVector &ghostDomainCM: The center of mass of the ghost domains
         */

        //Compute the centers of mass of each domain
        errorOut error = _inputProcessor.initializeIncrement( microIncrement, macroIncrement );
        if ( error ){
            errorOut result = new errorNode( "computeInitialCentersOfMass", "Error in initialization of the initial increment" );
            result->addNext( error );
            return result;
        }

        //Loop over the free micro-domains
        const stringVector* freeDomains = _inputProcessor.getFreeMicroDomainNames( );
        uIntVector domainNodes;

        floatVector domainCM( _dim );

        freeDomainMass = floatVector( freeDomains->size( ) );
        freeDomainCM = floatVector( _dim * freeDomains->size( ) );

        unsigned int indx = 0;

        for ( auto name = freeDomains->begin( ); name != freeDomains->end( ); name++ ){

            error = _inputProcessor._microscale->getDomainNodes( microIncrement, *name, domainNodes );

            if ( error ){

                errorOut result = new errorNode( "computeIncrementCentersOfMass", "Error in extraction of the free domain's nodes" );
                result->addNext( error );
                return result;

            }

            error = DOFProjection::computeDomainCenterOfMass( _dim, domainNodes, *_inputProcessor.getMicroVolumes( ),
                                                              *_inputProcessor.getMicroDensities( ),
                                                              *_inputProcessor.getMicroNodeReferencePositions( ),
                                                              *_inputProcessor.getMicroDisplacements( ),
                                                              *_inputProcessor.getMicroWeights( ),
                                                              freeDomainMass[ indx ], domainCM );

            if ( error ){

                errorOut result = new errorNode( "computeIncrementCentersOfMass", "Error in calculation of '" + *name + "' center of mass" );
                result->addNext( error );
                return result;

            }

            for ( unsigned int i = 0; i < _dim; i++ ){
                freeDomainCM[ _dim * indx + i ] = domainCM[ i ];
            }

            indx++;

        }

        //Loop over the ghost micro-domains
        const stringVector* ghostDomains = _inputProcessor.getGhostMicroDomainNames( );

        ghostDomainMass = floatVector( freeDomains->size( ) );
        ghostDomainCM = floatVector( _dim * freeDomains->size( ) );

        indx = 0;

        for ( auto name = ghostDomains->begin( ); name != ghostDomains->end( ); name++ ){

            error = _inputProcessor._microscale->getDomainNodes( microIncrement, *name, domainNodes );

            if ( error ){

                errorOut result = new errorNode( "computeIncrementCentersOfMass", "Error in extraction of the ghost domain's nodes" );
                result->addNext( error );
                return result;

            }

            error = DOFProjection::computeDomainCenterOfMass( _dim, domainNodes, *_inputProcessor.getMicroVolumes( ),
                                                              *_inputProcessor.getMicroDensities( ),
                                                              *_inputProcessor.getMicroNodeReferencePositions( ),
                                                              *_inputProcessor.getMicroWeights( ),
                                                              ghostDomainMass[ indx ], domainCM );

            if ( error ){

                errorOut result = new errorNode( "computeIncrementCentersOfMass", "Error in calculation of '" + *name + "' center of mass" );
                result->addNext( error );
                return result;

            }

            for ( unsigned int i = 0; i < _dim; i++ ){
                ghostDomainCM[ _dim * indx + i ] = domainCM[ i ];
            }

            indx++;

        }

        return NULL;
    }

    errorOut overlapCoupling::computeShapeFunctionsAtPoints( const unsigned int cellID,
                                                             const floatVector &nodeLocations,
                                                             const uIntVector &connectivity,
                                                             const uIntVector &connectivityCellIndices,
                                                             const floatVector &points,
                                                             floatVector &shapeFunctions ){
        /*!
         * Compute the shape functions of a given macro-scale domain at the given points
         *
         * :param const unsigned int &cellID: The cell ID at which to compute the shape functions
         * :param const floatVector &nodeLocations: The nodal location vector
         * :param const uIntVector &connectivity: The connectivity vector
         * :param const uIntVector &connectivityCellIndices: The indices of the different cells
         *     in the connectivity vector.
         * :param const floatVector &points: The points at which to compute the shape functions
         * :param floatVector &shapeFunctions: The shapefunctions at the points
         */

        //Make sure the cellID is allowable
        if ( cellID >= connectivityCellIndices.size( ) ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "The cellID is too large for the connectivity cell indices vector" );

        }

        //Make sure the number of points and the size of the output are consistent
        unsigned nPoints = points.size( ) / _dim;
        if ( ( points.size( ) % _dim ) > 0 ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "The points vector is inconsistent with the dimension\n"
                                  "    points.size( ): " + std::to_string( points.size( ) ) + "\n" +
                                  "    nPoints: " + std::to_string( nPoints ) );

        }

        //Get the XDMF cell type
        unsigned int index0 = connectivityCellIndices[ cellID ];
        unsigned int cellType = connectivity[ index0 ];

        //Get the element name
        auto it = elib::XDMFTypeToElementName.find( cellType );

        if ( it == elib::XDMFTypeToElementName.end( ) ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "The cell type " + std::to_string(cellType) + " is not supported" );

        }
        
        //Get the number of nodes
        auto it2 = elib::XDMFTypeToNodeCount.find( cellType );

        if ( it2 == elib::XDMFTypeToNodeCount.end( ) ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "The cell type " + std::to_string( cellType ) + " is not found in the node count map" );

        }

        shapeFunctions.clear();
        shapeFunctions.reserve( it2->second * nPoints );

        //Get the nodes from the file
        elib::vecOfvec nodes( it2->second, elib::vec( _dim, 0 ) );
        for ( unsigned int n = 0; n < it2->second; n++ ){

            for ( unsigned int i = 0; i < _dim; i++ ){

                nodes[ n ][ i ] = nodeLocations[ _dim * connectivity[ index0 + 1 + n ] + i ];

            }

        }
        
        //Get the element
        auto qrule = elib::default_qrules.find( it->second );
        if ( qrule == elib::default_qrules.end( ) ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "The element type " + it->second + " is not found in the default quadrature rules map" );

        }

        std::unique_ptr< elib::Element > element = elib::build_element_from_string( it->second, { }, nodes, qrule->second );

        //Compute the shape functions at each point
        floatVector point;
        floatVector localPosition, pointShapeFunctions;
        errorOut error;

        //Loop over the output vector
        for ( unsigned int p = 0; p < nPoints; p++ ){

            point = floatVector( points.begin( ) + _dim * p, points.begin( ) + _dim * ( p + 1 ) );

            error = element->compute_local_coordinates( point, localPosition );

            if ( error ) {

                return new errorNode( "computeShapeFunctionsAtPoints",
                                      "Error in computing the local coordinates for point " + std::to_string( p ) );

            }

            error = element->get_shape_functions( localPosition, pointShapeFunctions );

            if ( error ) {

                return new errorNode( "computeShapeFunctionsAtPoints",
                                      "Error in the computation of the shape functions for point " + std::to_string( p ) );

            }

            for ( unsigned int i = 0; i < pointShapeFunctions.size( ); i++ ){

                shapeFunctions.push_back( pointShapeFunctions[ i ] );

            }

        }

        return NULL;
            
    }

    errorOut overlapCoupling::computeShapeFunctionsAtPoints( const unsigned int cellID,
                                                             const floatVector &nodeReferenceLocations,
                                                             const floatVector &nodeDisplacements,
                                                             const uIntVector &connectivity,
                                                             const uIntVector &connectivityCellIndices,
                                                             const floatVector &points,
                                                             floatVector &shapeFunctions ){
        /*!
         * Compute the shape functions of a given macro-scale domain at the given points
         *
         * :param const unsigned int &cellID: The cell ID at which to compute the shape functions
         * :param const floatVector &nodeReferenceLocations: The nodal reference location vector
         * :param const floatVector &nodeDisplacements: The nodal reference location vector
         * :param const uIntVector &connectivity: The connectivity vector
         * :param const uIntVector &connectivityCellIndices: The indices of the different cells
         *     in the connectivity vector.
         * :param const floatVector &points: The points at which to compute the shape functions
         * :param floatVector &shapeFunctions: The shapefunctions at the points
         */

        //Make sure the cellID is allowable
        if ( cellID >= connectivityCellIndices.size( ) ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "The cellID is too large for the connectivity cell indices vector" );

        }

        //Make sure the number of points and the size of the output are consistent
        unsigned nPoints = points.size( ) / _dim;
        if ( ( points.size( ) % _dim ) > 0 ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "The points vector is inconsistent with the dimension\n"
                                  "    points.size( ): " + std::to_string( points.size( ) ) + "\n" +
                                  "    nPoints: " + std::to_string( nPoints ) );

        }

        //Get the XDMF cell type
        unsigned int index0 = connectivityCellIndices[ cellID ];
        unsigned int cellType = connectivity[ index0 ];

        //Get the element name
        auto it = elib::XDMFTypeToElementName.find( cellType );

        if ( it == elib::XDMFTypeToElementName.end( ) ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "The cell type " + std::to_string(cellType) + " is not supported" );

        }
        
        //Get the number of nodes
        auto it2 = elib::XDMFTypeToNodeCount.find( cellType );

        if ( it2 == elib::XDMFTypeToNodeCount.end( ) ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "The cell type " + std::to_string( cellType ) + " is not found in the node count map" );

        }

        shapeFunctions.clear();
        shapeFunctions.reserve( it2->second * nPoints );

        //Get the nodes from the file
        elib::vecOfvec nodes( it2->second, elib::vec( _dim, 0 ) );
        for ( unsigned int n = 0; n < it2->second; n++ ){

            for ( unsigned int i = 0; i < _dim; i++ ){

                nodes[ n ][ i ] = nodeReferenceLocations[ _dim * connectivity[ index0 + 1 + n ] + i ]
                                + nodeDisplacements[ _dim * connectivity[ index0 + 1 + n ] + i ];

            }

        }
        
        //Get the element
        auto qrule = elib::default_qrules.find( it->second );
        if ( qrule == elib::default_qrules.end( ) ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "The element type " + it->second + " is not found in the default quadrature rules map" );

        }

        std::unique_ptr< elib::Element > element = elib::build_element_from_string( it->second, { }, nodes, qrule->second );

        //Compute the shape functions at each point
        floatVector point;
        floatVector localPosition, pointShapeFunctions;
        errorOut error;

        //Loop over the output vector
        for ( unsigned int p = 0; p < nPoints; p++ ){

            point = floatVector( points.begin( ) + _dim * p, points.begin( ) + _dim * ( p + 1 ) );

            error = element->compute_local_coordinates( point, localPosition );

            if ( error ) {

                return new errorNode( "computeShapeFunctionsAtPoints",
                                      "Error in computing the local coordinates for point " + std::to_string( p ) );

            }

            error = element->get_shape_functions( localPosition, pointShapeFunctions );

            if ( error ) {

                return new errorNode( "computeShapeFunctionsAtPoints",
                                      "Error in the computation of the shape functions for point " + std::to_string( p ) );

            }

            for ( unsigned int i = 0; i < pointShapeFunctions.size( ); i++ ){

                shapeFunctions.push_back( pointShapeFunctions[ i ] );

            }

        }

        return NULL;
            
    }

    errorOut overlapCoupling::computeShapeFunctionsAtReferenceCentersOfMass( ){
        /*!
         * Compute the shape functions at the reference centers of mass
         */

        //Loop over the free domains
        unsigned int index;
        unsigned int comStart = 0;
        unsigned int nMicroDomains;
        floatVector microDomainCOMs;
        floatVector macroDomainShapeFunctions;
        errorOut error;

        const uIntVector *freeMacroCellIds = _inputProcessor.getFreeMacroCellIds( );
        const uIntVector *freeMacroCellMicroDomainCounts = _inputProcessor.getFreeMacroCellMicroDomainCounts( );

        _referenceGhostMicroDomainCenterOfMassShapeFunctions.clear( );

        for ( auto domain = freeMacroCellIds->begin( ); domain != freeMacroCellIds->end( ); domain++ ){

            //Set the index
            index = domain - freeMacroCellIds->begin( );

            //Get the number of free micro-domains
            nMicroDomains = ( *freeMacroCellMicroDomainCounts )[ index ];

            //Get the centers of mass of the macro-domain
            microDomainCOMs = floatVector( _referenceGhostMicroDomainCentersOfMass.begin( ) + _dim * comStart,
                                           _referenceGhostMicroDomainCentersOfMass.begin( ) + _dim * ( comStart + nMicroDomains ) );

            error = computeShapeFunctionsAtPoints( *domain,
                                                   *_inputProcessor.getMacroNodeReferencePositions( ),
                                                   *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                   *_inputProcessor.getMacroNodeReferenceConnectivityCellIndices( ),
                                                   microDomainCOMs, macroDomainShapeFunctions );

            if ( error ){

                errorOut result = new errorNode( "computeShapeFunctionsAtReferenceCentersOfMass",
                                                 "Error in computation of the shape functions at the points" );
                result->addNext( error );
                return result;

            }

            //Append the shape functions to the storage vector
            _referenceGhostMicroDomainCenterOfMassShapeFunctions.insert( _referenceGhostMicroDomainCenterOfMassShapeFunctions.end( ),
                                                                         macroDomainShapeFunctions.begin( ),
                                                                         macroDomainShapeFunctions.end( ) );

            //Update comStart
            comStart += nMicroDomains;
        }

        //Loop over the ghost domains
        
        comStart = 0;

        const uIntVector *ghostMacroCellIds = _inputProcessor.getGhostMacroCellIds( );
        const uIntVector *ghostMacroCellMicroDomainCounts = _inputProcessor.getGhostMacroCellMicroDomainCounts( );

        _referenceFreeMicroDomainCenterOfMassShapeFunctions.clear( );

        for ( auto domain = ghostMacroCellIds->begin( ); domain != ghostMacroCellIds->end( ); domain++ ){

            //Set the index
            index = domain - ghostMacroCellIds->begin( );

            //Get the number of free micro-domains
            nMicroDomains = ( *ghostMacroCellMicroDomainCounts )[ index ];

            //Get the centers of mass of the macro-domain
            microDomainCOMs = floatVector( _referenceFreeMicroDomainCentersOfMass.begin( ) + _dim * comStart,
                                           _referenceFreeMicroDomainCentersOfMass.begin( ) + _dim * ( comStart + nMicroDomains ) );

            error = computeShapeFunctionsAtPoints( *domain,
                                                   *_inputProcessor.getMacroNodeReferencePositions( ),
                                                   *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                   *_inputProcessor.getMacroNodeReferenceConnectivityCellIndices( ),
                                                   microDomainCOMs, macroDomainShapeFunctions );

            if ( error ){

                errorOut result = new errorNode( "computeShapeFunctionsAtReferenceCentersOfMass",
                                                 "Error in computation of the shape functions at the points" );
                result->addNext( error );
                return result;

            }

            //Append the shape functions to the storage vector
            _referenceFreeMicroDomainCenterOfMassShapeFunctions.insert( _referenceFreeMicroDomainCenterOfMassShapeFunctions.end( ),
                                                                        macroDomainShapeFunctions.begin( ),
                                                                        macroDomainShapeFunctions.end( ) );

            //Update comStart
            comStart += nMicroDomains;
        }

        return NULL;
        
    }

    errorOut overlapCoupling::addDomainContributionToInterpolationMatrix( const uIntVector  &domainNodes,
                                                                          const uIntVector  &macroNodes,
                                                                          const floatVector &domainReferenceXis,
                                                                          const floatVector &domainCenterOfMassShapeFunctionValues ){
        /*!
         * Add the contribution of a domain to the interpolation matrices
         *
         * :param const std::string &domainName: The name of the domain
         * :param const floatVector &domainReferenceXis: The micro-position vectors of the nodes within the domain
         * :param const floatVector &domainCenterOfMassShapeFunctionValues: The shape function values at the centers of 
         *     mass of the micro-domains.
         */

        //Initialize the sparse matrix
        SparseMatrix domainN;

        //Extract the DOF maps
        const DOFMap *microGlobalToLocalDOFMap = _inputProcessor.getMicroGlobalToLocalDOFMap( );
        const DOFMap *macroGlobalToLocalDOFMap = _inputProcessor.getMacroGlobalToLocalDOFMap( );

        //Form the interpolation matrix contributions from the current domain
        errorOut error = DOFProjection::formMacroDomainToMicroInterpolationMatrix( _dim,
                                                                                   microGlobalToLocalDOFMap->size( ),
                                                                                   macroGlobalToLocalDOFMap->size( ),
                                                                                   domainNodes, macroNodes, domainReferenceXis,
                                                                                   domainCenterOfMassShapeFunctionValues,
                                                                                   *_inputProcessor.getMicroWeights( ),
                                                                                   domainN,
                                                                                   microGlobalToLocalDOFMap,
                                                                                   macroGlobalToLocalDOFMap );

        if ( error ){

            errorOut result = new errorNode( "addDomainContributionToInterpolationMatrix",
                                             "Error in computation of the contribution of the domain to the interpolation matrix" );
            result->addNext( error );
            return result;

        }

        //Add the contribution to the total shapefunction matrix
        if ( _N.nonZeros( ) > 0 ){
            _N += domainN;
        }
        else{
            _N = domainN;
        }

        return NULL;

    }

    errorOut overlapCoupling::projectDegreesOfFreedom( ){
        /*!
         * Project the degrees of freedom of the ghost nodes
         * for the current increment
         */

        //Get the displacement vectors
        const floatVector *macroDispDOFVector = _inputProcessor.getMacroDispDOFVector( );
        const floatVector *microDisplacements = _inputProcessor.getMicroDisplacements( );

        //Get the free and ghost node ids
        const uIntVector *freeMacroNodeIds = _inputProcessor.getFreeMacroNodeIds( );
        const uIntVector *ghostMacroNodeIds = _inputProcessor.getGhostMacroNodeIds( );

        const uIntVector *freeMicroNodeIds = _inputProcessor.getFreeMicroNodeIds( );
        const uIntVector *ghostMicroNodeIds = _inputProcessor.getGhostMicroNodeIds( );

        //Set the number of displacement degrees of freedom
        unsigned int nMacroDispDOF = _dim + _dim * _dim;
        unsigned int nMicroDispDOF = _dim;

        //Assemble the free displacements
        floatVector freeMacroDisplacements( nMacroDispDOF * freeMacroNodeIds->size( ) );
        floatVector freeMicroDisplacements( nMicroDispDOF * freeMicroNodeIds->size( ) );

        for ( auto it = freeMacroNodeIds->begin( ); it != freeMacroNodeIds->end( ); it++ ){

            auto map = _inputProcessor.getMacroGlobalToLocalDOFMap( )->find( *it );

            if ( map == _inputProcessor.getMacroGlobalToLocalDOFMap( )->end( ) ){

                return new errorNode( "projectDegreesOfFreedom",
                                      "Global degree of freedom '" + std::to_string( *it ) + "' not found in degree of freedom map" );

            }

            //Set the macro displacements
            for ( unsigned int i = 0; i < nMacroDispDOF; i++ ){

                freeMacroDisplacements[ nMacroDispDOF * ( map->second ) + i ]
                    = ( *macroDispDOFVector )[ nMacroDispDOF * map->first + i ];

            }

            //Set the micro deformation phi

        }

        for ( auto it = freeMicroNodeIds->begin( ); it != freeMicroNodeIds->end( ); it++ ){

            for ( unsigned int i = 0; i < nMicroDispDOF; i++ ){

                freeMicroDisplacements[ nMicroDispDOF * ( it - freeMicroNodeIds->begin( ) ) + i ]
                    = ( *microDisplacements )[ nMicroDispDOF * ( *it ) + i ];

            }

        }

        //Map the macro and micro free displacements to Eigen matrices
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > Q( freeMicroDisplacements.data(), freeMicroDisplacements.size( ), 1 );
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > D( freeMacroDisplacements.data(), freeMacroDisplacements.size( ), 1 );

        //Map the output vectors to Eigen matrices
        _projected_ghost_macro_displacement.clear( );
        _projected_ghost_macro_displacement.resize( nMacroDispDOF * ghostMacroNodeIds->size( ) );

        _projected_ghost_micro_displacement.clear( );
        _projected_ghost_micro_displacement.resize( nMicroDispDOF * ghostMicroNodeIds->size( ) );

        Eigen::Map< Eigen::Matrix< floatType, -1,  1 > > Qhat( _projected_ghost_micro_displacement.data(),
                                                               _projected_ghost_micro_displacement.size( ), 1 );

        Eigen::Map< Eigen::Matrix< floatType, -1,  1 > > Dhat( _projected_ghost_macro_displacement.data(),
                                                               _projected_ghost_macro_displacement.size( ), 1 );

        YAML::Node config = _inputProcessor.getCouplingInitialization( );

        if ( config[ "projection_type" ].as< std::string >( ).compare( "l2_projection" ) == 0 ){

            Dhat = _L2_BDhatQ * Q + _L2_BDhatD * D;
            Qhat = _L2_BQhatQ * Q + _L2_BQhatD * D;

        }
        else if ( config[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){

            Dhat = _DP_BDhatQ * Q + _DP_BDhatD * D;
            Qhat = _DP_BQhatQ * Q + _DP_BQhatD * D;

        }
        else{

            return new errorNode( "projectDegreesOfFreedom",
                                  "'projection_type' '" + config[ "projection_type" ].as< std::string >( ) + "' is not recognized" );

        }

        std::cout << "exiting projection\n";
        return NULL;

    }

    const floatVector* overlapCoupling::getReferenceFreeMicroDomainMasses( ){
        /*!
         * Get access to the reference free micro-domain mass
         */

        return &_referenceFreeMicroDomainMasses;
    }

    const floatVector* overlapCoupling::getReferenceGhostMicroDomainMasses( ){
        /*!
         * Get access to the reference ghost micro-domain masses
         */

        return &_referenceGhostMicroDomainMasses;
    }

    const floatVector* overlapCoupling::getReferenceFreeMicroDomainCentersOfMass( ){
        /*!
         * Get access to the reference free micro-domain centers of mass
         */

        return &_referenceFreeMicroDomainCentersOfMass;
    }

    const floatVector* overlapCoupling::getReferenceGhostMicroDomainCentersOfMass( ){
        /*!
         * Get access to the reference ghost micro-domain centers of mass
         */

        return &_referenceGhostMicroDomainCentersOfMass;
    }

    const floatVector* overlapCoupling::getFreeMicroDomainMasses( ){
        /*!
         * Get access to the free micro-domain mass
         */

        return &_freeMicroDomainMasses;
    }

    const floatVector* overlapCoupling::getGhostMicroDomainMasses( ){
        /*!
         * Get access to the ghost micro-domain masses
         */

        return &_ghostMicroDomainMasses;
    }

    const floatVector* overlapCoupling::getFreeMicroDomainCentersOfMass( ){
        /*!
         * Get access to the free micro-domain centers of mass
         */

        return &_freeMicroDomainCentersOfMass;
    }

    const floatVector* overlapCoupling::getGhostMicroDomainCentersOfMass( ){
        /*!
         * Get access to the ghost micro-domain centers of mass
         */

        return &_ghostMicroDomainCentersOfMass;
    }

    const floatVector* overlapCoupling::getReferenceFreeMicroDomainCenterOfMassShapeFunctions( ){
        /*!
         * Get access to the shapefunction values of the reference free micro domain centers of mass
         */

        return &_referenceFreeMicroDomainCenterOfMassShapeFunctions;
    }

    const floatVector* overlapCoupling::getReferenceGhostMicroDomainCenterOfMassShapeFunctions( ){
        /*!
         * Get access to the shapefunction values of the reference ghost micro domain centers of mass
         */

        return &_referenceGhostMicroDomainCenterOfMassShapeFunctions; 
    }

    const floatVector* overlapCoupling::getProjectedGhostMacroDisplacement( ){
        /*!
         * Get access to the projected ghost macro displacements
         */

        return &_projected_ghost_macro_displacement;
    }

    const floatVector* overlapCoupling::getProjectedGhostMicroDisplacement( ){
        /*!
         * Get access to the projected ghost macro displacements
         */

        return &_projected_ghost_micro_displacement;
    }

}
