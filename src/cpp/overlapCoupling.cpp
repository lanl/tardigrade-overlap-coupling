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

    errorOut overlapCoupling::getConstructorError( ){
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

        //Homogenize the material properties at the micro-scale to the macro-scale
        error = homogenizeMicroScale( microIncrement );

        if ( error ){

            errorOut result = new errorNode( "processIncrement", "Error in the homogenization of the micro-scale to the macro-scale" );
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

        errorOut error = NULL;
        if ( couplingInitialization[ "type" ].as< std::string >( ).compare( "use_first_increment" ) == 0 ){
            error = setReferenceStateFromIncrement( 0, 0 );

            bool save_reference_positions = false;

            if ( couplingInitialization[ "projection_type" ].as< std::string >( ).compare(  "direct_projection" ) == 0 ){

                save_reference_positions = true;

            }

            if ( save_reference_positions ){

                _macroReferencePositions = *_inputProcessor.getMacroNodeReferencePositions( )
                                         + *_inputProcessor.getMacroDisplacements( );

                _microReferencePositions = *_inputProcessor.getMicroNodeReferencePositions( )
                                         + *_inputProcessor.getMicroDisplacements( );

            }
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

        if ( _inputProcessor.getCouplingInitialization( )[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){

            _macroNodeProjectedMass
                = floatVector( _inputProcessor.getMacroGlobalToLocalDOFMap( )->size( ), 0 );

            _macroNodeProjectedMassMomentOfInertia
                = floatVector( _dim * _dim * _inputProcessor.getMacroGlobalToLocalDOFMap( )->size( ), 0 );

            _macroNodeMassRelativePositionConstant
                = floatVector( _dim * _inputProcessor.getMacroGlobalToLocalDOFMap( )->size( ), 0 );

        }

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

                domainIndex += nMicroDomains;

        }

        //Compress the shape-function matrix
        _N.makeCompressed( );

        //Form the projectors
        error = formTheProjectors( microIncrement, macroIncrement );

        if ( error ){

            errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                             "Error in the formation of the projectors" );
            result->addNext( error );
            return result;

        }

        return NULL;
    }

    errorOut overlapCoupling::formTheProjectors( const unsigned int &microIncrement, const unsigned int &macroIncrement ){
        /*!
         * Form the projection operators
         *
         * :param const unsigned int &microIncrement: The increment in the micro-scale
         * :param const unsigned int &macroIncrement: The increment in the macro-scale
         */

        const YAML::Node config = _inputProcessor.getCouplingInitialization( );

        if ( config[ "projection_type" ].as< std::string >( ).compare( "l2_projection" ) == 0 ){

            errorOut error = formL2Projectors( );

            if ( error ){

                errorOut result = new errorNode( "formTheProjectors",
                                                 "Error in the formation of the L2 projectors" );
                result->addNext( error );
                return result;

            }

        }
        else if ( config[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){

            errorOut error = formDirectProjectionProjectors( microIncrement, macroIncrement );

            if ( error ){

                errorOut result = new errorNode( "formTheProjectors",
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
        SparseMatrix NQDhat = _N.block( 0, nFreeMacroDOF, nFreeMicroDOF, nGhostMacroDOF );
        NQDhat.makeCompressed( );
        solver.compute( NQDhat );//_N.block( 0, nFreeMacroDOF, nFreeMicroDOF, nGhostMacroDOF ) );

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

        return NULL;
    }

    errorOut overlapCoupling::formDirectProjectionProjectors( const unsigned int &microIncrement, const unsigned int &macroIncrement ){
        /*!
         * Form the projectors if the direct projection is to be used
         *
         * :param const unsigned int &microIncrement: The micro increment at which to form the projectors
         * :param const unsigned int &macroIncrement: The macro increment at which to form the projectors
         */

        //Get the ghost macro cell IDs
        const uIntVector *ghostMacroCellIDs = _inputProcessor.getGhostMacroCellIds( );
        const uIntVector *ghostMacroCellMicroDomainCounts = _inputProcessor.getGhostMacroCellMicroDomainCounts( );

        const stringVector *ghostMacroDomainNames = _inputProcessor.getGhostMacroDomainNames( );
        const stringVector *freeMicroDomainNames = _inputProcessor.getFreeMicroDomainNames( );

        unsigned int cellIndex;
        unsigned int nMicroDomains;

        errorOut error = NULL;

        uIntVector macroNodes;

        unsigned int domainIndex = 0;

        //Form the projector from the free micro-scale to the ghost macro-scale
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

                errorOut result = new errorNode( "formDirectProjectionProjectors",
                                                 "Error in extracting the ghost macro-node set" );
                result->addNext( error );
                return result;

            }

            //Loop over the free micro-domains
            for ( auto domain  = freeMicroDomainNames->begin( ) + domainIndex;
                       domain != freeMicroDomainNames->begin( ) + domainIndex + nMicroDomains;
                       domain++ ){

                error = addDomainContributionToDirectFreeMicroToGhostMacroProjector( cellIndex, *cellID, microIncrement, 
                                                                                     *domain, macroNodes );

                if ( error ){

                    errorOut result = new errorNode( "formDirectProjectionProjectors",
                                                     "Error in processing free micro-scale domain '" + *domain + "' for a ghost macro domain reference state" );
                    result->addNext( error );
                    return result;

                }

            }

            domainIndex += nMicroDomains;

        }

        //Set the dimension of the displacement DOF
        unsigned int nDispMicroDOF = _dim;

        unsigned int nDispMacroDOF = _dim + _dim * _dim;

        //Get the number of micro degrees of freedom which are free and ghost
        unsigned int nFreeMicroDOF = nDispMicroDOF * _inputProcessor.getFreeMicroNodeIds( )->size( );
        unsigned int nGhostMicroDOF = nDispMicroDOF * _inputProcessor.getGhostMicroNodeIds( )->size( );

        //Get the number of macro degrees of freedom which are free and ghost
        unsigned int nFreeMacroDOF = nDispMacroDOF * _inputProcessor.getFreeMacroNodeIds( )->size( );
        unsigned int nGhostMacroDOF = nDispMacroDOF * _inputProcessor.getGhostMacroNodeIds( )->size( );

        //Assemble the remaining projectors

        _DP_BDhatD = -_DP_BDhatQ * _N.block( 0, 0, nFreeMicroDOF, nFreeMacroDOF );

        _DP_BQhatQ = _N.block( nFreeMicroDOF, nFreeMacroDOF, nGhostMicroDOF, nGhostMacroDOF ) * _DP_BDhatQ;

        _DP_BQhatD = _N.block( nFreeMicroDOF, 0, nGhostMicroDOF, nFreeMacroDOF )
                   + _N.block( nFreeMicroDOF, nFreeMacroDOF, nGhostMicroDOF, nGhostMacroDOF ) * _DP_BDhatD;

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

        //If the projection time is the direct projection method, we need to save some values
        if ( _inputProcessor.getCouplingInitialization( )[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){

            //Save the contributions of the domain to the direct projection values
            error = addDomainToDirectProjectionReferenceValues( domainNodes, macroNodes, domainReferenceXiVectors,
                                                                domainMicroPositionShapeFunctionValues );

            if ( error ){

                errorOut result = new errorNode( "processDomainReference",
                                                 "Error in saving the direct projection reference values" );
                result->addNext( error );
                return result;

            }

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
                                             "Error in the computation of the shape function at the center of mass for a micro domain" );
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
                    microNodePositions[ _dim * index + i ] = ( *microReferencePositions )[ _dim * ( *it ) + i ]
                                                           + ( *microDisplacements )[ _dim * ( *it ) + i ];
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
                                                 "Error in the computation of the shape function at the center of mass for a micro domain" );
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
        errorOut error = NULL;

        //Loop over the output vector
        for ( unsigned int p = 0; p < nPoints; p++ ){

            point = floatVector( points.begin( ) + _dim * p, points.begin( ) + _dim * ( p + 1 ) );

            error = element->compute_local_coordinates( point, localPosition );

            if ( !element->local_point_inside( localPosition ) ){

                shapeFunctions.push_back( 0. );
                continue;

            }

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
        errorOut error = NULL;

        //Loop over the output vector
        for ( unsigned int p = 0; p < nPoints; p++ ){

            point = floatVector( points.begin( ) + _dim * p, points.begin( ) + _dim * ( p + 1 ) );

            error = element->compute_local_coordinates( point, localPosition );

            if ( !element->local_point_inside( localPosition ) ){

                shapeFunctions.push_back( 0. );
                continue;

            }

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
        errorOut error = NULL;

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

        return NULL;

    }

    errorOut overlapCoupling::addDomainToDirectProjectionReferenceValues( const uIntVector &domainNodes,
                                                                          const uIntVector &macroNodes,
                                                                          const floatVector &domainReferenceXiVectors,
                                                                          const floatVector &domainMicroPositionShapeFunctionValues ){
        /*!
         * Add the current domain information to the direct projection reference values
         *
         * :param const uIntVector &domainNodes: The nodes associated with the current domain
         * :param const uIntVector &macroNodes: The macro nodes associated with the current domain
         * :param const floatVector &domainReferenceXiVectors: The relative position vectors of the domain
         * :param const floatVector &domainMicroPositionShapeFunctionValues: The shape function values at the micro
         *     node positions.
         */

        const floatVector *microDensities = _inputProcessor.getMicroDensities( );
        const floatVector *microVolumes   = _inputProcessor.getMicroVolumes( );
        const floatVector *microWeights   = _inputProcessor.getMicroWeights( );

        //Additional values
        unsigned int m, n, p;

        floatType microMass, weight, sf;
        floatVector Xi;

        //Loop through the micro nodes
        for ( auto microNode  = domainNodes.begin( );
                   microNode != domainNodes.end( );
                   microNode++ ){

            //Set the index
            m = microNode - domainNodes.begin( );

            //Compute the micro-mass
            microMass = ( *microDensities )[ *microNode ] * ( *microVolumes )[ *microNode ];

            //Extract the micro Xi vector
            Xi = floatVector( domainReferenceXiVectors.begin( ) + _dim * m,
                              domainReferenceXiVectors.begin( ) + _dim * ( m + 1 ) );

            //Extract the weighting value
            weight = ( *microWeights )[ *microNode ];

            //Loop through the macro nodes
            for ( auto macroNode  = macroNodes.begin( );
                       macroNode != macroNodes.end( );
                       macroNode++ ){

                //Set the index
                n = macroNode - macroNodes.begin( );

                auto indx = _inputProcessor.getMacroGlobalToLocalDOFMap( )->find( *macroNode );

                if ( indx == _inputProcessor.getMacroGlobalToLocalDOFMap( )->end( ) ){

                    return new errorNode( "addDomainToDirectProjectionReferenceValues",
                                          "Macro node '" + std::to_string( n ) + "' not found in global to local macro node map" );

                }
                else{

                    p = indx->second;

                }

                //Get the shape function
                sf = domainMicroPositionShapeFunctionValues[ macroNodes.size( ) * m + n ];

                //Add the contribution to the nodal mass
                _macroNodeProjectedMass[ p ] += microMass * sf * weight;

                //Add the contribution to the mass-moment of inertia
                for ( unsigned int I = 0; I < _dim; I++ ){

                    for ( unsigned int J = 0; J < _dim; J++ ){

                        _macroNodeProjectedMassMomentOfInertia[ _dim * _dim * p + _dim * I + J ]
                            += microMass * sf * weight * Xi[ I ] * Xi[ J ];

                    }
                }

                //Add the contribution to the mass relative position constant
                for ( unsigned int I = 0; I < _dim; I++ ){

                    _macroNodeMassRelativePositionConstant[ _dim * p + I ] += microMass * sf * weight * Xi[ I ];

                }

            }

        }

        return NULL;
    }

    errorOut overlapCoupling::addDomainContributionToDirectFreeMicroToGhostMacroProjector( const unsigned int &cellIndex,
                                                                                           const unsigned int &cellID,
                                                                                           const unsigned int &microIncrement,
                                                                                           const std::string &domainName,
                                                                                           const uIntVector &macroNodes ){
        /*!
         * Compute the current domain's contribution to the direct free micro to ghost macro
         * projection matrix
         */

        //Compute the shape functions at the micro-nodes for the domain

        //Get the domain node ids
        uIntVector domainNodes;
        errorOut error = _inputProcessor._microscale->getDomainNodes( microIncrement, domainName, domainNodes );

        if ( error ){

            errorOut result = new errorNode( "addDomainContributionToDirectFreeMicroToGhostMacroProjector",
                                             "Error in extracting the domain ( " + domainName + " ) nodes" );
            result->addNext( error );
            return result;

        }

        //Get the micro-node positions and relative position vectors
        floatVector microNodePositions( _dim * domainNodes.size( ) );
 
        unsigned int index;
        const floatVector *microReferencePositions = _inputProcessor.getMicroNodeReferencePositions( );
        const floatVector *microDisplacements      = _inputProcessor.getMicroDisplacements( );
        floatVector domainReferenceXiVectors( _dim * domainNodes.size( ) );

        const floatVector domainCenterOfMass = floatVector( _referenceFreeMicroDomainCentersOfMass.begin( ) + _dim * cellIndex,
                                                            _referenceFreeMicroDomainCentersOfMass.begin( ) + _dim * ( cellIndex + 1 ) );
 
        for ( auto it = domainNodes.begin( ); it != domainNodes.end( ); it++ ){
 
            index = it - domainNodes.begin( );
            for ( unsigned int i = 0; i < _dim; i++ ){
                microNodePositions[ _dim * index + i ] = ( *microReferencePositions )[ _dim * ( *it ) + i ]
                                                       + ( *microDisplacements )[ _dim * ( *it ) + i ];

                domainReferenceXiVectors[ _dim * index + i ] = microNodePositions[ _dim * index + i ]
                                                             - domainCenterOfMass[ i ];
            }
 
        }
 
        //Compute the shape function values at the micro positions
        floatVector domainMicroPositionShapeFunctionValues;
        error = computeShapeFunctionsAtPoints( cellID,
                                               *_inputProcessor.getMacroNodeReferencePositions( ),
                                               *_inputProcessor.getMacroDisplacements( ),
                                               *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                               *_inputProcessor.getMacroNodeReferenceConnectivityCellIndices( ),
                                               microNodePositions,
                                               domainMicroPositionShapeFunctionValues );

        if ( error ){
 
            errorOut result = new errorNode( "addDomainContributionToDirectFreeMicroToGhostMacroProjector",
                                             "Error in the computation of the shape functions at the center of mass for a micro domain" );
            result->addNext( error );
            return result;
 
        }

        //Construct the contribution of the domain
        
        //Get the domain's mass properties
        floatVector domainMacroNodeProjectedMass( macroNodes.size( ) );
        floatVector domainMacroNodeProjectedMassMomentOfInertia( _dim * _dim * macroNodes.size( ) );
        floatVector domainMacroNodeProjectedMassRelativePositionConstant( _dim * macroNodes.size( ) );
        
        auto microGlobalToLocalDOFMap = _inputProcessor.getMicroGlobalToLocalDOFMap( );
        auto macroGlobalToLocalDOFMap = _inputProcessor.getMacroGlobalToLocalDOFMap( );

        //We also specify the global to local map for the projector.
        //This is not the same as the existing global map for the macro nodes
        //because we store the free nodes first, so we have to shift them by 
        //the number of free macro-scale nodes
        DOFMap projectorMacroGlobalToLocalDOFMap;
        projectorMacroGlobalToLocalDOFMap.reserve( _inputProcessor.getGhostMacroNodeIds( )->size( ) );

        for ( auto md  = macroNodes.begin( );
                   md != macroNodes.end( );
                   md++ ){

            auto indx = macroGlobalToLocalDOFMap->find( *md );

            if ( indx == microGlobalToLocalDOFMap->end( ) ){

                return new errorNode( "addDomainContributionToDirectFreeMicroToGhostMacroProjector",
                                      "'" + std::to_string( *md ) + "' not found in the DOF map" );

            }

            projectorMacroGlobalToLocalDOFMap.emplace( indx->first, indx->second - _inputProcessor.getFreeMacroNodeIds( )->size( ) );

            unsigned int mdi = md - macroNodes.begin( );

            domainMacroNodeProjectedMass[ mdi ] = _macroNodeProjectedMass[ indx->second ];
            
            for ( unsigned int i = 0; i < _dim; i++ ){

                domainMacroNodeProjectedMassRelativePositionConstant[ _dim * mdi + i ]
                    = _macroNodeMassRelativePositionConstant[ _dim * indx->second + i ];

                for ( unsigned int j = 0; j < _dim; j++ ){

                    domainMacroNodeProjectedMassMomentOfInertia[ _dim * _dim * mdi + _dim * i + j ]
                        = _macroNodeProjectedMassMomentOfInertia[ _dim * _dim * indx->second + _dim * i + j ];

                }

            }

        }

        //Construct the projector contribution for this domain
        SparseMatrix domainProjector;

        error = DOFProjection::formMicroDomainToMacroProjectionMatrix( _dim,
                                                                       _inputProcessor.getFreeMicroNodeIds( )->size( ),
                                                                       _inputProcessor.getGhostMacroNodeIds( )->size( ),
                                                                       domainNodes, macroNodes,
                                                                       *_inputProcessor.getMicroVolumes( ),
                                                                       *_inputProcessor.getMicroDensities( ),
                                                                       *_inputProcessor.getMicroWeights( ),
                                                                       domainReferenceXiVectors,
                                                                       domainMicroPositionShapeFunctionValues,
                                                                       domainMacroNodeProjectedMass,
                                                                       domainMacroNodeProjectedMassMomentOfInertia,
                                                                       domainMacroNodeProjectedMassRelativePositionConstant,
                                                                       domainProjector,
                                                                       microGlobalToLocalDOFMap,
                                                                       &projectorMacroGlobalToLocalDOFMap );

        if ( error ){
            errorOut result = new errorNode( "addDomainContributionToDirectFreeMicroToGhostMacroProjector",
                                             "Error in the computation of the domain's contribution to the micro to macro projection matrix" );
            result->addNext( error );
            return result;
        }

        if ( _DP_BQhatQ.nonZeros( ) == 0 ){

            _DP_BDhatQ = domainProjector;

        }
        else{

            _DP_BDhatQ += domainProjector;

        }

        return NULL;

    }

    errorOut overlapCoupling::homogenizeMicroScale( const unsigned int &microIncrement ){
        /*!
         * Homogenize the micro-scale properties to the macro scale.
         *
         * :param const unsigned int &microIncrement: The increment at the micro-scale to homogenize
         */

        //Loop through the free macro-scale cells
        unsigned int microDomainStartIndex = 0;
        errorOut error = NULL;

        uIntVector microDomainNodeIds;
        floatVector microNodePositions;
        std::shared_ptr< volumeReconstruction::volumeReconstructionBase > reconstructedVolume;

        for ( auto macroCell  = _inputProcessor.getFreeMacroCellIds( )->begin( );
                   macroCell != _inputProcessor.getFreeMacroCellIds( )->end( );
                   macroCell++ ){

            //Set the macro index
            unsigned int macroIndex = macroCell - _inputProcessor.getFreeMacroCellIds( )->begin( );

            //Get the number of micro domain in this macro cell
            unsigned int nCellMicroDomains = ( *_inputProcessor.getFreeMacroCellMicroDomainCounts( ) )[ macroIndex ];

            //Domain micro centers of mass
            floatVector microDomainCentersOfMass( _ghostMicroDomainCentersOfMass.begin( ) + _dim * microDomainStartIndex,
                                                  _ghostMicroDomainCentersOfMass.begin( ) + _dim * ( microDomainStartIndex + nCellMicroDomains ) );

            //Domain surface appproximate number of decompositions
            const uIntVector *microDomainSurfaceDecompositions = _inputProcessor.getGhostMicroSurfaceApproximateSplitCount( );

            unsigned int microIndex = microDomainStartIndex;

            for ( auto microDomain  = _inputProcessor.getGhostMicroDomainNames( )->begin( ) + microDomainStartIndex;
                       microDomain != _inputProcessor.getGhostMicroDomainNames( )->begin( ) + microDomainStartIndex + nCellMicroDomains;
                       microDomain++, microIndex++ ){

                microNodePositions.clear( );
                reconstructedVolume.reset( );

                //Reconstruct the micro-domain's volume
                error = reconstructDomain( microIncrement, *microDomain, microDomainNodeIds, microNodePositions, reconstructedVolume );

                if ( error ){

                    errorOut result = new errorNode( "homogenizeMicroScale",
                                                     "Error in the reconstruction of the microscale domain" );
                    result->addNext( error );
                    return result;

                }

                //Compute the volume averages
                floatVector domainCenterOfMass( microDomainCentersOfMass.begin( ) + _dim * microIndex,
                                                microDomainCentersOfMass.begin( ) + _dim * ( microIndex + 1 ) );

                error = computeDomainVolumeAverages( *macroCell, microDomainNodeIds,
                                                     reconstructedVolume, &domainCenterOfMass );

                if ( error ){

                    errorOut result = new errorNode( "computeDomainVolumeAverages",
                                                     "Error in the computation of the volume averages of the microscale domain" );
                    result->addNext( error );
                    return result;

                }
                
                //Compute the surface averages
                error = computeDomainSurfaceAverages( *macroCell, microDomainNodeIds,
                                                      ( *microDomainSurfaceDecompositions )[ microIndex ],
                                                      reconstructedVolume );

                if ( error ){

                    errorOut result = new errorNode( "homogenizeMicroScale",
                                                     "Error in the computation of the surface averages of the microscale domain" );
                    result->addNext( error );
                    return result;

                }
                
            }

            //Compute the approximate stresses
            error = computeHomogenizedStresses( *macroCell );

            if ( error ){

                errorOut result = new errorNode( "homogenizeMicroScale",
                                                 "Error in the computation of the homogenized stresses" );
                result->addNext( error );
                return result;

            }

            //Increment the start index of the micro domain
            microDomainStartIndex += nCellMicroDomains;

        }

//        //Loop through the ghost macro-scale cells
//        microDomainStartIndex = 0;
//        for ( auto macroCell  = _inputProcessor.getGhostMacroDomainNames( )->begin( );
//                   macroCell != _inputProcessor.getGhostMacroDomainNames( )->end( );
//                   macroCell++ ){
//
//            //Set the macro index
//            unsigned int macroIndex = macroCell - _inputProcessor.getGhostMacroDomainNames( )->begin( );
//
//            //Get the number of micro domain in this macro cell
//            unsigned int nCellMicroDomains = ( *_inputProcessor.getGhostMacroCellMicroDomainCounts( ) )[ macroIndex ];
//
//            for ( auto microDomain  = _inputProcessor.getFreeMicroDomainNames( )->begin( ) + microDomainStartIndex;
//                       microDomain != _inputProcessor.getFreeMicroDomainNames( )->begin( ) + microDomainStartIndex + nCellMicroDomains;
//                       microDomain++ ){
//
//                floatVector microNodePositions;
//                std::shared_ptr< volumeReconstruction::volumeReconstructionBase > reconstructedVolume;
//
//                //Compute the volume averages
//                error = computeDomainVolumeAverages( *macroCell, *microDomain, reconstructedVolume );
//                
//                //Compute the surface averages
//                
//            }
//
//            //Compute the approximate stresses
//
//            //Increment the start index of the micro domain
//            microDomainStartIndex += nCellMicroDomains;
//
//        }

        return NULL;
    }

    errorOut overlapCoupling::reconstructDomain( const unsigned int &microIncrement, const std::string &microDomainName,
                                                 uIntVector &microDomainNodes, floatVector &microNodePositions,
                                                 std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume ){
        /*!
         * Reconstruct the micro-domain's volume to perform volume and surface integrals over that
         * domain.
         *
         * :param const unsigned int &microIncrement: The increment at which to extract the micro-positions
         * :param const std::string &microDomainName: The name of the micro-domain to be re-constructed.
         * :param uIntVector &microDomainNodes: The nodes associated with the micro domain
         * :param floatVector &microNodePositions: The positions of the micro nodes for the current domain
         * :param std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume: The reconstructed
         *     volume ready for additional processing.
         */

        //Get the domain node ids
        errorOut error = _inputProcessor._microscale->getDomainNodes( microIncrement, microDomainName, microDomainNodes );

        if ( error ){

            errorOut result = new errorNode( "reconstructDomain",
                                             "Error in getting the node ids for the domain ( " + microDomainName + " )" );
            result->addNext( error );
            return result;

        }

        //Get the micro-node positions
        microNodePositions.clear( );
        microNodePositions.resize( _dim * microDomainNodes.size( ) );
 
        unsigned int index;
        const floatVector *microReferencePositions = _inputProcessor.getMicroNodeReferencePositions( );
        const floatVector *microDisplacements      = _inputProcessor.getMicroDisplacements( );
 
        for ( auto it = microDomainNodes.begin( ); it != microDomainNodes.end( ); it++ ){
 
            index = it - microDomainNodes.begin( );

            for ( unsigned int i = 0; i < _dim; i++ ){

                microNodePositions[ _dim * index + i ] = ( *microReferencePositions )[ _dim * ( *it ) + i ]
                                                       + ( *microDisplacements )[ _dim * ( *it ) + i ];

            }
 
        }

        //Pass the base name of the output file to the volume reconstruction configuration to be used if output has been requested
        YAML::Node volumeReconstructionConfig = _inputProcessor.getVolumeReconstructionConfig( );
        volumeReconstructionConfig[ "baseOutputFilename" ] = microDomainName + "_" + std::to_string( microIncrement );

        //Get the volume reconstruction object
        reconstructedVolume
            = volumeReconstruction::volumeReconstructionBase( volumeReconstructionConfig ).create( );

        if ( reconstructedVolume->getError( ) ){

            errorOut result = new errorNode( "reconstructDomain",
                                             "Error in creating the volume reconstruction object for " + microDomainName );

            result->addNext( reconstructedVolume->getError( ) );
            return result;

        }

        //Load the micro points
        error = reconstructedVolume->loadPoints( &microNodePositions );

        if ( error ){

            errorOut result = new errorNode( "reconstructDomain",
                                             "Error in loading the micro-scale points for " + microDomainName );
            result->addNext( error );
            return result;

        }

        //Reconstruct the volume
        error = reconstructedVolume->evaluate( );

        if ( error ){

            errorOut result = new errorNode( "reconstructDomain",
                                             "Error in loading the micro-scale points for " + microDomainName );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut overlapCoupling::computeDomainVolumeAverages( const uIntType &macroCellID, const uIntVector &microDomainNodeIDs,
                                                           std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume,
                                                           const floatVector *microDomainCenterOfMass ){
        /*!
         * Compute the required volume averages over the micro-domain.
         *
         * :param const uIntType &macroCellID: The ID number of the macro-cell associated with the micro domain.
         * :param const uIntVector &microDomainNodeIDs: The micro domain's node ids
         * :param volumeReconstruction::volumeReconstructionBase &reconstructedVolume: The reconstructed volume
         *     ready to have volume integrals computed over.
         * :param const floatVector *microDomainCenterOfMass: The center of mass for the micro-domain computed
         *     directly from the particles.
         */

        //Assemble the averaging vector at the nodes for non-volume weighted averaging quantities
        unsigned int dataCountAtPoint = 1  //Volume calculation
                                      + 1  //Density
                                      + 9; //Stress

        if ( _inputProcessor.useReconstructedMassCenters( ) ){
            dataCountAtPoint += _dim; //Domain center of mass
        }

        //Add the micro body force if it is defined
        if ( _inputProcessor.microBodyForceDefined( ) ){

            dataCountAtPoint += _dim; //Add the micro body force

        }

        //Add the micro inertia term if the acceleration is defined
        if ( _inputProcessor.microAccelerationDefined( ) ){

            dataCountAtPoint += _dim; //Add the micro inertia term

        }

        floatVector dataAtMicroPoints( dataCountAtPoint * microDomainNodeIDs.size( ), 0 );

        const floatVector *microDensities = _inputProcessor.getMicroDensities( );

        const floatVector *microBodyForces = _inputProcessor.getMicroBodyForces( );

        const floatVector *microAccelerations = _inputProcessor.getMicroAccelerations( );

        const floatVector *microReferencePositions = _inputProcessor.getMicroNodeReferencePositions( );

        const floatVector *microDisplacements = _inputProcessor.getMicroDisplacements( );

        const floatVector *microStresses = _inputProcessor.getMicroStresses( );

        if ( microStresses->size( ) != _dim * _dim * microDensities->size( ) ){

            return new errorNode( "computeDomainVolumeAverages",
                                  "The micro stress vector size is not consistent wtih the number of points and the dimension" );

        }

        unsigned int index = 0;
        unsigned int localIndex = 0;

        for ( auto node = microDomainNodeIDs.begin( ); node != microDomainNodeIDs.end( ); node++, index++ ){

            dataAtMicroPoints[ dataCountAtPoint * index + 0 ] = 1.;                           //Integrate the volume of the domain
            dataAtMicroPoints[ dataCountAtPoint * index + 1 ] = ( *microDensities )[ *node ]; //Integrate the density of the domain

            //Integrate the micro stresses
            for ( unsigned int i = 0; i < _dim * _dim; i++ ){
                dataAtMicroPoints[ dataCountAtPoint * index + 2 + i ] = ( *microStresses )[ _dim * _dim * ( *node ) + i ];
            }

            localIndex = 11;

            if ( _inputProcessor.useReconstructedMassCenters( ) ){

                //Integrate for the domain's center of mass
                for ( unsigned int i = 0; i < _dim; i++ ){
    
                    dataAtMicroPoints[ dataCountAtPoint * index + localIndex + i ] =
                        ( *microDensities )[ *node ] * ( ( *microReferencePositions )[ _dim * ( *node )  + i ]
                                                       + ( *microDisplacements )[ _dim * ( *node ) + i ] );
    
                }

                //Local index
                localIndex += 3;

            }

            //Add the micro body forces
            if ( _inputProcessor.microBodyForceDefined( ) ){

                for ( unsigned int i = 0; i < _dim; i++ ){

                    dataAtMicroPoints[ dataCountAtPoint * index + localIndex + i ]
                        = ( *microBodyForces )[ _dim * ( *node ) + i ]; //Integrate the body forces of the domain

                }

                localIndex += _dim;

            }

            //Add the micro acclerations
            if ( _inputProcessor.microAccelerationDefined( ) ){

                for ( unsigned int i = 0; i < _dim; i++ ){

                    dataAtMicroPoints[ dataCountAtPoint * index + localIndex + i ]
                        = ( *microDensities )[ *node ] * ( *microAccelerations )[ _dim * ( *node ) + i ]; //Integrate the accelerations of the domain

                }

                localIndex += _dim;

            }

        }

        //Compute the initial volume averages
        floatVector integratedValues;
        errorOut error = reconstructedVolume->performVolumeIntegration( dataAtMicroPoints, dataCountAtPoint, integratedValues );

        if ( error ){

            errorOut result = new errorNode( "computeDomainVolumeAverages",
                                             "Error in computing the initial volume averages" );
            result->addNext( error );
            return result;

        }

        if ( homogenizedVolumes.find( macroCellID ) == homogenizedVolumes.end( ) ){

            //Save the values
            floatVector tmp = { integratedValues[ 0 ] };
            homogenizedVolumes.emplace( macroCellID, tmp ); //Save the volume

            tmp = { integratedValues[ 1 ] / integratedValues[ 0 ] };
            homogenizedDensities.emplace( macroCellID, tmp ); //Save the density

            tmp = floatVector( integratedValues.begin( ) + 2,
                               integratedValues.begin( ) + 2 + _dim * _dim ) / integratedValues[ 0 ];
            homogenizedSymmetricMicroStresses.emplace( macroCellID, tmp ); //Save the symmetric micro stress

            localIndex = 11;

            if ( _inputProcessor.useReconstructedMassCenters( ) ){

                tmp = floatVector( integratedValues.begin( ) + localIndex, integratedValues.begin( ) + localIndex + _dim ) / integratedValues[ 1 ];
                homogenizedCentersOfMass.emplace( macroCellID, tmp ); //Save the center of mass

                localIndex += 3;

            }
            else{

                homogenizedCentersOfMass.emplace( macroCellID, *microDomainCenterOfMass );

            }

            //Save the body force
            if ( _inputProcessor.microBodyForceDefined( ) ){

                tmp.resize( _dim );

                for ( unsigned int i = 0; i < _dim; i++ ){

                    tmp[ i ] = integratedValues[ localIndex + i ] / integratedValues[ 0 ];

                }

                localIndex += _dim;

            }
            else{

                tmp = *microBodyForces;

            }

            homogenizedBodyForces.emplace( macroCellID, tmp );

            //Save the acceleration
            if ( _inputProcessor.microAccelerationDefined( ) ){

                tmp.resize( _dim );

                for ( unsigned int i = 0; i < _dim; i++ ){

                    tmp[ i ] = integratedValues[ localIndex + i ] / integratedValues[ 0 ];

                }

                localIndex += _dim;

            }
            else{

                tmp = *microAccelerations;

            }

            homogenizedAccelerations.emplace( macroCellID, tmp );

        }
        else{

            homogenizedVolumes[ macroCellID ].push_back( integratedValues[ 0 ] ); //Save the volume

            homogenizedDensities[ macroCellID ].push_back( integratedValues[ 1 ] / integratedValues[ 0 ] ); //Save the density

            homogenizedSymmetricMicroStresses[ macroCellID ]
                = vectorTools::appendVectors( { homogenizedSymmetricMicroStresses[ macroCellID ],
                                                floatVector( integratedValues.begin( ) + 2,
                                                             integratedValues.begin( ) + 2 + _dim * _dim ) / integratedValues[ 0 ] } );


            localIndex = 11;

            //Save the centers of mass
            if ( _inputProcessor.useReconstructedMassCenters( ) ){

                homogenizedCentersOfMass[ macroCellID ]
                    = vectorTools::appendVectors( { homogenizedCentersOfMass[ macroCellID ],
                                                    floatVector( integratedValues.begin( ) + localIndex,
                                                                 integratedValues.begin( ) + localIndex + _dim ) / integratedValues[ 1 ] } );

                localIndex += 3;
            }
            else{

                for ( unsigned int i = 0; i < _dim; i++ ){

                    homogenizedCentersOfMass[ macroCellID ].push_back( ( *microDomainCenterOfMass )[ i ] );

                }

            }

            //Save the body forces
            if ( _inputProcessor.microBodyForceDefined( ) ){

                homogenizedBodyForces[ macroCellID ]
                    = vectorTools::appendVectors( { homogenizedBodyForces[ macroCellID ],
                                                    floatVector( integratedValues.begin( ) + localIndex,
                                                                 integratedValues.begin( ) + localIndex + _dim ) / integratedValues[ 1 ] } );

                localIndex += _dim;

            }
            else{

                for ( unsigned int i = 0; i < _dim; i++ ){

                    homogenizedBodyForces[ macroCellID ].push_back( ( *microBodyForces )[ i ] );

                }

            }

            //Save the accelerations
            if ( _inputProcessor.microAccelerationDefined( ) ){

                homogenizedAccelerations[ macroCellID ]
                    = vectorTools::appendVectors( { homogenizedAccelerations[ macroCellID ],
                                                    floatVector( integratedValues.begin( ) + localIndex,
                                                                 integratedValues.begin( ) + localIndex + _dim ) / integratedValues[ 1 ] } );

                localIndex += _dim;

            }
            else{

                for ( unsigned int i = 0; i < _dim; i++ ){

                    homogenizedAccelerations[ macroCellID ].push_back( ( *microAccelerations )[ i ] );

                }

            }

        }

        //Perform the relative position volume integrations
        dataCountAtPoint = 0;

        //Add the micro body force couple if it is defined
        if ( _inputProcessor.microBodyForceDefined( ) ){

            dataCountAtPoint += _dim * _dim; //Add the micro body force couple

        }

        //Add the micro inertia term if the acceleration is defined
        if ( _inputProcessor.microAccelerationDefined( ) ){

            dataCountAtPoint += _dim * _dim; //Add the micro inertia term

        }

        integratedValues.clear( );
        integratedValues.resize( 0 );
        floatVector microPosition( _dim );

        if ( dataCountAtPoint > 0 ){

            index = 0;
            floatVector dataAtMicroPoints( dataCountAtPoint * microDomainNodeIDs.size( ), 0 );
    
            for ( auto node = microDomainNodeIDs.begin( ); node != microDomainNodeIDs.end( ); node++, index++ ){

                //Extract the micro relative position
                for ( unsigned int i = 0; i < _dim; i++ ){

                    microPosition[ i ] = ( ( *microReferencePositions )[ _dim * ( *node )  + i ]
                                         + ( *microDisplacements )[ _dim * ( *node ) + i ] );

                }

                localIndex = 0;

                //Add the contributions to the micro body couple
                if ( _inputProcessor.microBodyForceDefined( ) ){

                    floatVector microBodyForce( microBodyForces->begin( ) + _dim * ( *node ),
                                                microBodyForces->begin( ) + _dim * ( *node + 1 ) );

                    floatVector integrand
                        = ( *microDensities )[ *node ] * vectorTools::appendVectors( vectorTools::dyadic( microBodyForce, microPosition ) );
    
                    for ( unsigned int i = 0; i < _dim * _dim; i++ ){
    
                        dataAtMicroPoints[ dataCountAtPoint * index + localIndex + i ]
                            = integrand[ i ]; //Integrate the body force couple over the domain
    
                    }
    
                    localIndex += _dim * _dim;
    
                }
    
                //Add the contributions to the micro inertia
                if ( _inputProcessor.microAccelerationDefined( ) ){

                    floatVector microRelativeAcceleration( microAccelerations->begin( ) + _dim * ( *node ),
                                                           microAccelerations->begin( ) + _dim * ( *node + 1 ) );
                    microRelativeAcceleration -= floatVector( homogenizedAccelerations[ macroCellID ].end( ) - _dim,
                                                              homogenizedAccelerations[ macroCellID ].end( ) );

                    floatVector integrand
                        = ( *microDensities )[ *node ] * vectorTools::appendVectors( vectorTools::dyadic( microRelativeAcceleration, microPosition ) );
    
                    for ( unsigned int i = 0; i < _dim * _dim; i++ ){
    
                        dataAtMicroPoints[ dataCountAtPoint * index + localIndex + i ]
                            = integrand[ i ]; //Integrate the accelerations of the domain
    
                    }
    
                    localIndex += _dim * _dim;
    
                }

            }

            floatVector centerOfMass( homogenizedCentersOfMass[ macroCellID ].end( ) - _dim,
                                      homogenizedCentersOfMass[ macroCellID ].end( ) );

            error = reconstructedVolume->performRelativePositionVolumeIntegration( dataAtMicroPoints, dataCountAtPoint,
                                                                                   centerOfMass, integratedValues );

            if ( error ){

                errorOut result = new errorNode( "computeDomainVolumeAverages",
                                                 "Error in the computation of the relative position volume integration" );
                result->addNext( error );
                return result;

            }

        }

        if ( !_inputProcessor.microBodyForceDefined( ) ){

            integratedValues = vectorTools::appendVectors( { floatVector( _dim * _dim, 0 ), integratedValues } );

        }
        if ( !_inputProcessor.microAccelerationDefined( ) ){

            integratedValues = vectorTools::appendVectors( { integratedValues, floatVector( _dim * _dim, 0 ) } );

        }

        if ( homogenizedBodyForceCouples.find( macroCellID ) == homogenizedBodyForceCouples.end( ) ){

            floatVector tmp( integratedValues.begin( ),
                             integratedValues.begin( ) + _dim * _dim );
            tmp /= ( homogenizedVolumes[ macroCellID ].back( ) * homogenizedDensities[ macroCellID ].back( ) );
            homogenizedBodyForceCouples.emplace( macroCellID, tmp );

            tmp = floatVector( integratedValues.begin( ) + _dim * _dim,
                               integratedValues.begin( ) + 2 * _dim * _dim );
            tmp /= ( homogenizedVolumes[ macroCellID ].back( ) * homogenizedDensities[ macroCellID ].back( ) );
            homogenizedMicroInertias.emplace( macroCellID, tmp );

        }
        else{
            floatVector tmp( integratedValues.begin( ),
                             integratedValues.begin( ) + _dim * _dim );
            tmp /= ( homogenizedVolumes[ macroCellID ].back( ) * homogenizedDensities[ macroCellID ].back( ) );
            homogenizedBodyForceCouples[ macroCellID ] = vectorTools::appendVectors( { homogenizedBodyForceCouples[ macroCellID ],
                                                                                       tmp } );
            tmp = floatVector( integratedValues.begin( ) + _dim * _dim,
                               integratedValues.begin( ) + 2 * _dim * _dim );
            tmp /= ( homogenizedVolumes[ macroCellID ].back( ) * homogenizedDensities[ macroCellID ].back( ) );
            homogenizedMicroInertias[ macroCellID ] = vectorTools::appendVectors( { homogenizedMicroInertias[ macroCellID ],
                                                                                    tmp } );

        }

        return NULL;

    }

    errorOut overlapCoupling::computeDomainSurfaceAverages( const uIntType &macroCellID, const uIntVector &microDomainNodeIDs,
                                                            const uIntType &microDomainSurfaceDecompositionCount,
                                                            std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume ){
        /*!
         * Compute the required surface averages over the micro-domain.
         *
         * :param const uIntType &macroCellID: The ID of the macro-cell associated with the micro domain.
         * :param const std::string &microDomainNodeIDs: The IDs of the nodes in the micro-domain to have the surface averages
         *     computed over.
         * :param const uIntType &microDomainSurfaceDecompositionCount: The approximate number of regions to split the surface
         *     of the micro surface into.
         * :param volumeReconstruction::volumeReconstructionBase &reconstructedVolume: The reconstructed volume
         *     ready to have surface integrals computed over.
         */

        //Extract the required micro-scale values
        const floatVector *microDensities = _inputProcessor.getMicroDensities( );
        const floatVector *microDisplacements = _inputProcessor.getMicroDisplacements( );
        const floatVector *microReferencePositions = _inputProcessor.getMicroNodeReferencePositions( );
        const floatVector *microStresses = _inputProcessor.getMicroStresses( );

        //Check the micro-stress size
        if ( microStresses->size( ) != _dim * _dim * microDensities->size( ) ){

            return new errorNode( "computeDomainSurfaceAverages",
                                  "The micro-stress vector is not consistent with the dimension and number of points" );

        }

        /*=====================================================================
        |           Compute the reconstructed domain's surface area           |
        =====================================================================*/

        uIntType dataCountAtPoint = 1; //Total surface area

        floatVector dataAtMicroPoints( microDomainNodeIDs.size( ), 1 );

        floatVector integratedValue;

        errorOut error = reconstructedVolume->performSurfaceIntegration( dataAtMicroPoints, 1, integratedValue );

        if ( error ){

            errorOut result = new errorNode( "computeDomainSurfaceAverages",
                                             "Error in the computation of the domain's surface area" );
            result->addNext( error );
            return result;

        }

        if ( homogenizedSurfaceAreas.find( macroCellID ) == homogenizedSurfaceAreas.end( ) ){

            homogenizedSurfaceAreas[ macroCellID ] = { integratedValue[ 0 ] };

        }
        else{

            homogenizedSurfaceAreas[ macroCellID ].push_back( integratedValue[ 0 ] );

        }

        /*=====================================================================
        |           Compute the properties of the surface subdomains          |
        =====================================================================*/

        uIntVector subdomainNodeCounts;
        uIntVector subdomainNodeIDs;

        floatType minSurfaceSpacing
            = std::sqrt( homogenizedSurfaceAreas[ macroCellID ].back( ) / ( 3.14159 * microDomainSurfaceDecompositionCount ) );

        error = reconstructedVolume->getSurfaceSubdomains( minSurfaceSpacing, subdomainNodeCounts, subdomainNodeIDs );

        if ( error ){

            errorOut result = new errorNode( "computeDomainSurfaceAverages", "Error in extracting of the reconstructed volume's surface subdomains" );
            result->addNext( error );
            return result;

        }

        //Get the centers of mass of the surface regions
        floatVector surfaceRegionCentersOfMass( _dim * subdomainNodeCounts.size( ) );

        dataCountAtPoint = 1     //Surface area of region
                         + 1     //Mass of region
                         + _dim; //Micro point position

        dataAtMicroPoints.clear( );
        dataAtMicroPoints.reserve( dataCountAtPoint * microDomainNodeIDs.size( ) );

        for ( auto node = microDomainNodeIDs.begin( ); node != microDomainNodeIDs.end( ); node++ ){

            dataAtMicroPoints.push_back( 1 );
            dataAtMicroPoints.push_back( ( *microDensities )[ *node ] );

            floatVector microPoint = floatVector( microReferencePositions->begin( ) + _dim * ( *node ),
                                                  microReferencePositions->begin( ) + _dim * ( ( *node ) + 1 ) )
                                   + floatVector( microDisplacements->begin( ) + _dim * ( *node ),
                                                  microDisplacements->begin( ) + _dim * ( ( *node ) + 1 ) );

            for ( unsigned int i = 0; i < microPoint.size( ); i++ ){

                dataAtMicroPoints.push_back( ( *microDensities )[ *node ] * microPoint[ i ] );

            }

        }

        unsigned int startPoint = 0;

        if ( homogenizedSurfaceRegionCentersOfMass.find( macroCellID ) == homogenizedSurfaceRegionCentersOfMass.end( ) ){

            homogenizedSurfaceRegionAreas.emplace( macroCellID, floatVector( 0 ) );
            homogenizedSurfaceRegionDensities.emplace( macroCellID, floatVector( 0 ) );
            homogenizedSurfaceRegionCentersOfMass.emplace( macroCellID, floatVector( 0 ) );
            homogenizedSurfaceRegionTractions.emplace( macroCellID, floatVector( 0 ) );
            homogenizedSurfaceRegionCouples.emplace( macroCellID, floatVector( 0 ) );

        }

        for ( auto sNC = subdomainNodeCounts.begin( ); sNC != subdomainNodeCounts.end( ); sNC++ ){

            //Perform the surface integration
            uIntVector nodesInDomain( subdomainNodeIDs.begin( ) + startPoint,
                                      subdomainNodeIDs.begin( ) + startPoint + *sNC );

            error = reconstructedVolume->performSurfaceIntegration( dataAtMicroPoints, dataCountAtPoint,
                                                                    integratedValue, &nodesInDomain );

            if ( error ){

                errorOut result = new errorNode( "computeDomainSurfaceAverages",
                                                 "Error in the integration of the micro region ( "
                                                 + std::to_string( sNC - subdomainNodeCounts.begin( ) ) + " )" );
                result->addNext( error );
                return result;

            }

            //Extract the region surface areas and the region surface densities
            homogenizedSurfaceRegionAreas[ macroCellID ].push_back( integratedValue[ 0 ] );
            homogenizedSurfaceRegionDensities[ macroCellID ].push_back( integratedValue[ 1 ] / integratedValue[ 0 ] );

            //Extract the centers of mass of the surface regions
            floatVector regionCenterOfMass( integratedValue.begin( ) + 2,
                                            integratedValue.begin( ) + dataCountAtPoint );
            regionCenterOfMass /= integratedValue[ 1 ];

            homogenizedSurfaceRegionCentersOfMass[ macroCellID ]
                = vectorTools::appendVectors( { homogenizedSurfaceRegionCentersOfMass[ macroCellID ],
                                                regionCenterOfMass } );

            startPoint += *sNC;
        }

        /*===================================================================================
        | Compute the surface tractions and couples over the micro domain's surface regions |
        ===================================================================================*/

        dataCountAtPoint = _dim * _dim;

        dataAtMicroPoints.clear( );
        dataAtMicroPoints.reserve( dataCountAtPoint * microDomainNodeIDs.size( ) );

        for ( auto node = microDomainNodeIDs.begin( ); node != microDomainNodeIDs.end( ); node++ ){

            floatVector microStress( microStresses->begin( ) + _dim * _dim * ( *node ),
                                     microStresses->begin( ) + _dim * _dim * ( ( *node ) + 1 ) );

            for ( unsigned int i = 0; i < _dim * _dim; i++ ){

                dataAtMicroPoints.push_back( microStress[ i ] );

            }

        }

        startPoint = 0;

        for ( auto sNC = subdomainNodeCounts.begin( ); sNC != subdomainNodeCounts.end( ); sNC++ ){

            uIntVector nodesInDomain( subdomainNodeIDs.begin( ) + startPoint,
                                      subdomainNodeIDs.begin( ) + startPoint + *sNC );

            //Compute the tractions
            error = reconstructedVolume->performSurfaceFluxIntegration( dataAtMicroPoints, dataCountAtPoint,
                                                                        integratedValue, &nodesInDomain );

            if ( error ){

                errorOut result = new errorNode( "computeDomainSurfaceAverages",
                                                 "Error in the computation of the surface traction of the the micro region ( "
                                                 + std::to_string( sNC - subdomainNodeCounts.begin( ) ) + " )" );
                result->addNext( error );
                return result;

            }

            unsigned int regionID = homogenizedSurfaceRegionTractions[ macroCellID ].size( ) / _dim;

            floatType regionSurfaceArea = homogenizedSurfaceRegionAreas[ macroCellID ][ regionID ];

            integratedValue /= regionSurfaceArea;

            homogenizedSurfaceRegionTractions[ macroCellID ]
                = vectorTools::appendVectors( { homogenizedSurfaceRegionTractions[ macroCellID ],
                                                integratedValue } );

            //Compute the couples
            floatVector regionCenterOfMass( homogenizedSurfaceRegionCentersOfMass[ macroCellID ].begin( ) + _dim * regionID,
                                            homogenizedSurfaceRegionCentersOfMass[ macroCellID ].begin( ) + _dim * ( regionID + 1 ) );

            error = reconstructedVolume->performRelativePositionSurfaceFluxIntegration( dataAtMicroPoints, dataCountAtPoint,
                                                                                        regionCenterOfMass, integratedValue,
                                                                                        &nodesInDomain );

            if ( error ){

                errorOut result = new errorNode( "computeDomainSurfaceAverages",
                                                 "Error in the computation of the surface couple of the micro region ( "
                                                 + std::to_string( sNC - subdomainNodeCounts.begin( ) ) + " )" );
                result->addNext( error );
                return result;

            }

            integratedValue /= regionSurfaceArea;

            homogenizedSurfaceRegionCouples[ macroCellID ]
                = vectorTools::appendVectors( { homogenizedSurfaceRegionCouples[ macroCellID ],
                                                integratedValue } );
            

            startPoint += *sNC;

        }

        return NULL;
    }

    errorOut overlapCoupling::computeHomogenizedStresses( const uIntType &macroCellID ){
        /*!
         * Compute the homogenized stresses for the macro cell
         *
         * :param const uIntType &macroCellID: The ID of the macro cell to have the 
         *     homogenized stresses computed at the quadrature points.
         */

        ( void ) macroCellID;

        return new errorNode( "computeHomogenizedStresses", "Error: Not implemented" );
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
