/*!============================================================================
|                               overlapCoupling                               |
===============================================================================
| The implementation of the Overlap coupling method                           |
=============================================================================*/

#include<overlapCoupling.h>
#include<Eigen/SparseQR>
#include<micromorphic_tools.h>
#include<balance_equations.h>

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

        //Assemble the mass matrix for the free micromorphic domians
        error = assembleFreeMicromorphicMassMatrix( );

        if ( error ){

            errorOut result = new errorNode( "processIncrement", "Error in the assembly of the mass matrix for the free macro domains" );
            result->addNext( error );
            return result;

        }

        //Assemble the coupling mass and damping matrices
        error = assembleCouplingMassAndDampingMatrices( );

        if ( error ){

            errorOut result = new errorNode( "processIncrement", "Error in the construction of the coupling mass and damping matrices" );
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

    errorOut overlapCoupling::buildMacroDomainElement( const unsigned int cellID,
                                                       const floatVector &nodeLocations,
                                                       const uIntVector &connectivity,
                                                       const uIntVector &connectivityCellIndices,
                                                       std::unique_ptr< elib::Element > &element ){
        /*!
         * Construct a finite element representation of the macro domain
         *
         * :param const unsigned int cellID: The macro cell ID number
         * :param const floatVector &nodeLocations: The nodal location vector
         * :param const uIntVector &connectivity: The connectivity vector
         * :param const uIntVector &connectivityCellIndices: The indices of the different cells
         *     in the connectivity vector.
         * :param std::unique_ptr< elib::Element > &element: The element representation of the macro domain
         */

        //Make sure the cellID is allowable
        if ( cellID >= connectivityCellIndices.size( ) ){

            return new errorNode( "buildMacroDomainElement",
                                  "The cellID is too large for the connectivity cell indices vector" );

        }

        //Get the XDMF cell type
        unsigned int index0 = connectivityCellIndices[ cellID ];
        unsigned int cellType = connectivity[ index0 ];

        //Get the element name
        auto it = elib::XDMFTypeToElementName.find( cellType );

        if ( it == elib::XDMFTypeToElementName.end( ) ){

            return new errorNode( "buildMacroDomainElement",
                                  "The cell type " + std::to_string(cellType) + " is not supported" );

        }
        
        //Get the number of nodes
        auto it2 = elib::XDMFTypeToNodeCount.find( cellType );

        if ( it2 == elib::XDMFTypeToNodeCount.end( ) ){

            return new errorNode( "buildMacroDomainElement",
                                  "The cell type " + std::to_string( cellType ) + " is not found in the node count map" );

        }

        //Get the nodes from the file
        elib::vecOfvec nodes( it2->second, elib::vec( _dim, 0 ) );
        uIntVector globalNodeIds( connectivity.begin( ) + index0 + 1,
                                  connectivity.begin( ) + index0 + 1 + it2->second );
        for ( unsigned int n = 0; n < it2->second; n++ ){

            for ( unsigned int i = 0; i < _dim; i++ ){

                nodes[ n ][ i ] = nodeLocations[ _dim * connectivity[ index0 + 1 + n ] + i ];

            }

        }
        
        //Get the element
        auto qrule = elib::default_qrules.find( it->second );
        if ( qrule == elib::default_qrules.end( ) ){

            return new errorNode( "buildMacroDomainElement",
                                  "The element type " + it->second + " is not found in the default quadrature rules map" );

        }

        element = elib::build_element_from_string( it->second, globalNodeIds, nodes, qrule->second );

        return NULL;
    }

    errorOut overlapCoupling::buildMacroDomainElement( const unsigned int cellID,
                                                       const floatVector &nodeReferenceLocations,
                                                       const floatVector &nodeDisplacements,
                                                       const uIntVector &connectivity,
                                                       const uIntVector &connectivityCellIndices,
                                                       std::unique_ptr< elib::Element > &element ){
        /*!
         * Construct a finite element representation of the macro domain
         *
         * :param const unsigned int cellID: The macro cell ID number
         * :param const floatVector &nodeReferenceLocations: The nodal reference location vector
         * :param const floatVector &nodeDisplacements: The nodal displacement vector
         * :param const uIntVector &connectivity: The connectivity vector
         * :param const uIntVector &connectivityCellIndices: The indices of the different cells
         *     in the connectivity vector.
         * :param std::unique_ptr< elib::Element > &element: The element representation of the macro domain
         */

        //Make sure the cellID is allowable
        if ( cellID >= connectivityCellIndices.size( ) ){

            return new errorNode( "buildMacroDomainElement",
                                  "The cellID is too large for the connectivity cell indices vector" );

        }

        //Get the XDMF cell type
        unsigned int index0 = connectivityCellIndices[ cellID ];
        unsigned int cellType = connectivity[ index0 ];

        //Get the element name
        auto it = elib::XDMFTypeToElementName.find( cellType );

        if ( it == elib::XDMFTypeToElementName.end( ) ){

            return new errorNode( "buildMacroDomainElement",
                                  "The cell type " + std::to_string(cellType) + " is not supported" );

        }
        
        //Get the number of nodes
        auto it2 = elib::XDMFTypeToNodeCount.find( cellType );

        if ( it2 == elib::XDMFTypeToNodeCount.end( ) ){

            return new errorNode( "buildMacroDomainElement",
                                  "The cell type " + std::to_string( cellType ) + " is not found in the node count map" );

        }

        //Get the nodes from the file
        elib::vecOfvec referenceNodes( it2->second, elib::vec( _dim, 0 ) );
        elib::vecOfvec displacements( it2->second, elib::vec( _dim, 0 ) );
        uIntVector globalNodeIds( connectivity.begin( ) + index0 + 1,
                                  connectivity.begin( ) + index0 + 1 + it2->second );
        for ( unsigned int n = 0; n < it2->second; n++ ){

            for ( unsigned int i = 0; i < _dim; i++ ){

                referenceNodes[ n ][ i ] = nodeReferenceLocations[ _dim * connectivity[ index0 + 1 + n ] + i ];
                displacements[ n ][ i ] = nodeDisplacements[ _dim * connectivity[ index0 + 1 + n ] + i ];

            }

        }
        
        //Get the element
        auto qrule = elib::default_qrules.find( it->second );
        if ( qrule == elib::default_qrules.end( ) ){

            return new errorNode( "buildMacroDomainElement",
                                  "The element type " + it->second + " is not found in the default quadrature rules map" );

        }

        element = elib::build_element_from_string( it->second, globalNodeIds, referenceNodes, qrule->second );
        element->update_node_positions( displacements );

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

        //Build the element representing the macro-scale domain
        std::unique_ptr< elib::Element > element;
        errorOut error = overlapCoupling::buildMacroDomainElement( cellID, nodeLocations, connectivity,
                                                                   connectivityCellIndices, element );

        //Make sure the number of points and the size of the output are consistent
        unsigned nPoints = points.size( ) / _dim;
        if ( ( points.size( ) % _dim ) > 0 ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "The points vector is inconsistent with the dimension\n"
                                  "    points.size( ): " + std::to_string( points.size( ) ) + "\n" +
                                  "    nPoints: " + std::to_string( nPoints ) );

        }

        shapeFunctions.clear();
        shapeFunctions.reserve( element->reference_nodes.size( ) * nPoints );

        //Compute the shape functions at each point
        floatVector point;
        floatVector localPosition, pointShapeFunctions;

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

        //Build the element representing the macro-scale domain
        std::unique_ptr< elib::Element > element;
        errorOut error = overlapCoupling::buildMacroDomainElement( cellID, nodeReferenceLocations, nodeDisplacements,
                                                                   connectivity, connectivityCellIndices, element );

        //Make sure the number of points and the size of the output are consistent
        unsigned nPoints = points.size( ) / _dim;
        if ( ( points.size( ) % _dim ) > 0 ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "The points vector is inconsistent with the dimension\n"
                                  "    points.size( ): " + std::to_string( points.size( ) ) + "\n" +
                                  "    nPoints: " + std::to_string( nPoints ) );

        }

        shapeFunctions.clear();
        shapeFunctions.reserve( element->reference_nodes.size( ) * nPoints );

        //Compute the shape functions at each point
        floatVector point;
        floatVector localPosition, pointShapeFunctions;

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

    errorOut overlapCoupling::computeShapeFunctionGradientsAtPoints( const unsigned int cellID,
                                                                     const floatVector &nodeReferenceLocations,
                                                                     const floatVector &nodeDisplacements,
                                                                     const uIntVector &connectivity,
                                                                     const uIntVector &connectivityCellIndices,
                                                                     const floatVector &points,
                                                                     floatVector &shapeFunctionGradients ){
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
         * :param floatVector &shapeFunctionGradients: The gradients of the shapefunctions at the points. The 
         *     output vector is organized [ N11_x, N11_y, N11_z, N12_x, ... ] where the first index is the
         *     point number, the second index is the node number, and the _ indicates a gradient w.r.t. the
         *     third index.
         */

        //Build the element representing the macro-scale domain
        std::unique_ptr< elib::Element > element;
        errorOut error = overlapCoupling::buildMacroDomainElement( cellID, nodeReferenceLocations, nodeDisplacements,
                                                                   connectivity, connectivityCellIndices, element );

        //Make sure the number of points and the size of the output are consistent
        unsigned nPoints = points.size( ) / _dim;
        if ( ( points.size( ) % _dim ) > 0 ){

            return new errorNode( "computeShapeFunctionGradientsAtPoints",
                                  "The points vector is inconsistent with the dimension\n"
                                  "    points.size( ): " + std::to_string( points.size( ) ) + "\n" +
                                  "    nPoints: " + std::to_string( nPoints ) );

        }

        shapeFunctionGradients.clear();
        shapeFunctionGradients.reserve( element->reference_nodes.size( ) * nPoints * _dim );

        //Compute the shape functions at each point
        floatVector point;
        floatMatrix dNdx;
        floatVector localPosition, pointShapeFunctionGradientsVec;

        //Loop over the output vector
        for ( unsigned int p = 0; p < nPoints; p++ ){

            point = floatVector( points.begin( ) + _dim * p, points.begin( ) + _dim * ( p + 1 ) );

            error = element->compute_local_coordinates( point, localPosition );

            if ( !element->local_point_inside( localPosition ) ){

                shapeFunctionGradients.push_back( 0. );
                continue;

            }

            if ( error ) {

                return new errorNode( "computeShapeFunctionGradientsAtPoints",
                                      "Error in computing the local coordinates for point " + std::to_string( p ) );

            }

            error = element->get_global_shapefunction_gradients( localPosition, dNdx );

            if ( error ) {

                return new errorNode( "computeShapeFunctionGradientsAtPoints",
                                      "Error in the computation of the shape functions for point " + std::to_string( p ) );

            }

            pointShapeFunctionGradientsVec = vectorTools::appendVectors( dNdx );

            for ( unsigned int i = 0; i < pointShapeFunctionGradientsVec.size( ); i++ ){

                shapeFunctionGradients.push_back( pointShapeFunctionGradientsVec[ i ] );

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

        //Clear all of the previous values at the micro domains
        homogenizedVolumes.clear( );
        homogenizedSurfaceAreas.clear( );
        homogenizedDensities.clear( );
        homogenizedMicroInertias.clear( );
        homogenizedCentersOfMass.clear( );
        homogenizedBodyForces.clear( );
        homogenizedBodyForceCouples.clear( );
        homogenizedAccelerations.clear( );
        homogenizedMicroSpinInertias.clear( );
        homogenizedSymmetricMicroStresses.clear( );
        homogenizedSurfaceRegionAreas.clear( );
        homogenizedSurfaceRegionDensities.clear( );
        homogenizedSurfaceRegionCentersOfMass.clear( );
        homogenizedSurfaceRegionTractions.clear( );
        homogenizedSurfaceRegionCouples.clear( );

        //Clear all of the previous values at the quadrature points
        quadraturePointDensities.clear( );
        quadraturePointBodyForce.clear( );
        quadraturePointAccelerations.clear( );
        quadraturePointMicroInertias.clear( );
        quadraturePointBodyCouples.clear( );
        quadraturePointMicroSpinInertias.clear( );
        quadraturePointSymmetricMicroStress.clear( );
        quadraturePointCauchyStress.clear( );
        quadraturePointHigherOrderStress.clear( );

        //Clear the external forces at the nodes
        externalForcesAtNodes.clear( );
        externalCouplesAtNodes.clear( );

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

        //Loop through the ghost macro-scale cells
        microDomainStartIndex = 0;
        for ( auto macroCell  = _inputProcessor.getGhostMacroCellIds( )->begin( );
                   macroCell != _inputProcessor.getGhostMacroCellIds( )->end( );
                   macroCell++ ){

            //Set the macro index
            unsigned int macroIndex = macroCell - _inputProcessor.getGhostMacroCellIds( )->begin( );

            //Get the number of micro domain in this macro cell
            unsigned int nCellMicroDomains = ( *_inputProcessor.getGhostMacroCellMicroDomainCounts( ) )[ macroIndex ];

            //Domain micro centers of mass
            floatVector microDomainCentersOfMass( _freeMicroDomainCentersOfMass.begin( ) + _dim * microDomainStartIndex,
                                                  _freeMicroDomainCentersOfMass.begin( ) + _dim * ( microDomainStartIndex + nCellMicroDomains ) );

            //Domain surface appproximate number of decompositions
            const uIntVector *microDomainSurfaceDecompositions = _inputProcessor.getFreeMicroSurfaceApproximateSplitCount( );

            unsigned int microIndex = microDomainStartIndex;

            for ( auto microDomain  = _inputProcessor.getFreeMicroDomainNames( )->begin( ) + microDomainStartIndex;
                       microDomain != _inputProcessor.getFreeMicroDomainNames( )->begin( ) + microDomainStartIndex + nCellMicroDomains;
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

        //Compute the homogenized force vectors and mass matrices
        error = assembleHomogenizedMatricesAndVectors( );

        if ( error ){

            errorOut result = new errorNode( "homogenizeMicroScale",
                                             "Error in the computation of the homogenized forces and mass matrix" );
            result->addNext( error );
            return result;

        }

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

        unsigned int initialOffset = 11;

        if ( _inputProcessor.useReconstructedMassCenters( ) ){
            dataCountAtPoint += _dim; //Domain center of mass
        }

        //Add the micro body force if it is defined
        if ( _inputProcessor.microBodyForceDefined( ) ){

            dataCountAtPoint += _dim; //Add the micro body force

        }

        //Add the micro spin inertia term if the acceleration is defined
        if ( _inputProcessor.microAccelerationDefined( ) ){

            dataCountAtPoint += _dim; //Add the micro spin inertia term

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

            localIndex = initialOffset;

            if ( _inputProcessor.useReconstructedMassCenters( ) ){

                //Integrate for the domain's center of mass
                for ( unsigned int i = 0; i < _dim; i++ ){
    
                    dataAtMicroPoints[ dataCountAtPoint * index + localIndex + i ] =
                        ( *microDensities )[ *node ] * ( ( *microReferencePositions )[ _dim * ( *node )  + i ]
                                                       + ( *microDisplacements )[ _dim * ( *node ) + i ] );
    
                }

                //Local index
                localIndex += _dim;

            }

            //Add the micro body forces
            if ( _inputProcessor.microBodyForceDefined( ) ){

                for ( unsigned int i = 0; i < _dim; i++ ){

                    dataAtMicroPoints[ dataCountAtPoint * index + localIndex + i ]
                        = ( *microDensities )[ *node ] * ( *microBodyForces )[ _dim * ( *node ) + i ]; //Integrate the body forces of the domain

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

        //Splice in the micro domain's center of mass if required
        if ( !_inputProcessor.useReconstructedMassCenters( ) ){

            integratedValues
                = vectorTools::appendVectors( { floatVector( integratedValues.begin( ), integratedValues.begin( ) + initialOffset ),
                                                integratedValues[ 1 ] * ( *microDomainCenterOfMass ),
                                                floatVector( integratedValues.begin( ) + 11, integratedValues.end( ) ) } );

        }

        //Splice in the micro-body force if required
        if ( !_inputProcessor.microBodyForceDefined( ) ){

            integratedValues
                = vectorTools::appendVectors( { floatVector( integratedValues.begin( ),
                                                             integratedValues.begin( ) + initialOffset + _dim ),
                                                floatVector( _dim, 0 ),
                                                floatVector( integratedValues.begin( ) + initialOffset + _dim,
                                                             integratedValues.end( ) ) } );

        }

        //Splice in the micro acceleration if required
        if ( !_inputProcessor.microAccelerationDefined( ) ){

            integratedValues
                = vectorTools::appendVectors( { integratedValues, floatVector( _dim, 0 ) } );

        }

        //Save the values
        if ( homogenizedVolumes.find( macroCellID ) == homogenizedVolumes.end( ) ){

            //Save the volume
            floatVector tmp = { integratedValues[ 0 ] };
            homogenizedVolumes.emplace( macroCellID, tmp );

            //Save the density
            tmp = { integratedValues[ 1 ] / integratedValues[ 0 ] };
            homogenizedDensities.emplace( macroCellID, tmp );

            //Save the symmetric micro stress
            tmp = floatVector( integratedValues.begin( ) + 2,
                               integratedValues.begin( ) + 2 + _dim * _dim ) / integratedValues[ 0 ];
            homogenizedSymmetricMicroStresses.emplace( macroCellID, tmp );

            //Save the center of mass
            tmp = floatVector( integratedValues.begin( ) + initialOffset,
                               integratedValues.begin( ) + initialOffset + _dim );
            tmp /= integratedValues[ 1 ];
            homogenizedCentersOfMass.emplace( macroCellID, tmp );

            //Save the body force
            tmp = floatVector( integratedValues.begin( ) + initialOffset + _dim,
                               integratedValues.begin( ) + initialOffset + 2 * _dim );
            tmp /= integratedValues[ 1 ];
            homogenizedBodyForces.emplace( macroCellID, tmp );

            //Save the acceleration
            tmp = floatVector( integratedValues.begin( ) + initialOffset + 2 * _dim,
                               integratedValues.end( ) );
            tmp /= integratedValues[ 1 ];
            homogenizedAccelerations.emplace( macroCellID, tmp );

        }
        else{

            //Save the volume
            homogenizedVolumes[ macroCellID ].push_back( integratedValues[ 0 ] );

            //Save the density
            homogenizedDensities[ macroCellID ].push_back( integratedValues[ 1 ] / integratedValues[ 0 ] );

            //Save the symmetric micro-stress
            homogenizedSymmetricMicroStresses[ macroCellID ]
                = vectorTools::appendVectors( { homogenizedSymmetricMicroStresses[ macroCellID ],
                                                floatVector( integratedValues.begin( ) + 2,
                                                             integratedValues.begin( ) + 2 + _dim * _dim ) / integratedValues[ 0 ] } );

            //Save the centers of mass
            homogenizedCentersOfMass[ macroCellID ]
                = vectorTools::appendVectors( { homogenizedCentersOfMass[ macroCellID ],
                                                floatVector( integratedValues.begin( ) + initialOffset,
                                                             integratedValues.begin( ) + initialOffset + _dim ) / integratedValues[ 1 ] } );

            //Save the body forces
            homogenizedBodyForces[ macroCellID ]
                = vectorTools::appendVectors( { homogenizedBodyForces[ macroCellID ],
                                                floatVector( integratedValues.begin( ) + initialOffset + _dim,
                                                              integratedValues.begin( ) + initialOffset + 2 * _dim ) / integratedValues[ 1 ] } );

            //Save the accelerations
            homogenizedAccelerations[ macroCellID ]
                = vectorTools::appendVectors( { homogenizedAccelerations[ macroCellID ],
                                                floatVector( integratedValues.begin( ) + initialOffset + 2 * _dim,
                                                              integratedValues.end( ) ) / integratedValues[ 1 ] } );

        }

        //Perform the relative position volume integrations
        dataCountAtPoint = _dim * _dim; //The micro inertia

        initialOffset = dataCountAtPoint; //Set the initial offset

        //Add the micro body force couple if it is defined
        if ( _inputProcessor.microBodyForceDefined( ) ){

            dataCountAtPoint += _dim * _dim; //Add the micro body force couple

        }

        //Add the micro spin inertia term if the acceleration is defined
        if ( _inputProcessor.microAccelerationDefined( ) ){

            dataCountAtPoint += _dim * _dim; //Add the micro spin inertia term

        }

        integratedValues.clear( );
        integratedValues.resize( 0 );
        floatVector microRelativePosition( _dim );

        if ( dataCountAtPoint > 0 ){

            index = 0;
            floatVector dataAtMicroPoints( dataCountAtPoint * microDomainNodeIDs.size( ), 0 );

            //Set the center of mass
            floatVector centerOfMass( homogenizedCentersOfMass[ macroCellID ].end( ) - _dim,
                                      homogenizedCentersOfMass[ macroCellID ].end( ) );
    
            for ( auto node = microDomainNodeIDs.begin( ); node != microDomainNodeIDs.end( ); node++, index++ ){

                //Extract the micro relative position
                microRelativePosition = floatVector( microReferencePositions->begin( ) + _dim * ( *node ),
                                                     microReferencePositions->begin( ) + _dim * ( ( *node ) + 1 ) )
                                      + floatVector( microDisplacements->begin( ) + _dim * ( *node ),
                                                     microDisplacements->begin( ) + _dim * ( ( *node ) + 1 ) )
                                      - centerOfMass;

                floatVector integrand
                    = ( *microDensities )[ *node ]
                    * vectorTools::appendVectors( vectorTools::dyadic( microRelativePosition, microRelativePosition ) );

                //Add the contributions to the micro inertia
                for ( unsigned int i = 0; i < _dim * _dim; i++ ){

                    dataAtMicroPoints[ dataCountAtPoint * index + i ]
                        = integrand[ i ];

                }

                localIndex = initialOffset;

                //Add the contributions to the micro body couple
                if ( _inputProcessor.microBodyForceDefined( ) ){

                    floatVector microBodyForce( microBodyForces->begin( ) + _dim * ( *node ),
                                                microBodyForces->begin( ) + _dim * ( *node + 1 ) );

                    floatVector integrand
                        = ( *microDensities )[ *node ]
                        * vectorTools::appendVectors( vectorTools::dyadic( microBodyForce, microRelativePosition ) );
    
                    for ( unsigned int i = 0; i < _dim * _dim; i++ ){
    
                        dataAtMicroPoints[ dataCountAtPoint * index + localIndex + i ]
                            = integrand[ i ]; //Integrate the body force couple over the domain
    
                    }
    
                    localIndex += _dim * _dim;
    
                }
    
                //Add the contributions to the micro spin inertia
                if ( _inputProcessor.microAccelerationDefined( ) ){

                    floatVector microRelativeAcceleration( microAccelerations->begin( ) + _dim * ( *node ),
                                                           microAccelerations->begin( ) + _dim * ( *node + 1 ) );
                    microRelativeAcceleration -= floatVector( homogenizedAccelerations[ macroCellID ].end( ) - _dim,
                                                              homogenizedAccelerations[ macroCellID ].end( ) );

                    floatVector integrand
                        = ( *microDensities )[ *node ]
                        * vectorTools::appendVectors( vectorTools::dyadic( microRelativeAcceleration, microRelativePosition ) );
    
                    for ( unsigned int i = 0; i < _dim * _dim; i++ ){
    
                        dataAtMicroPoints[ dataCountAtPoint * index + localIndex + i ]
                            = integrand[ i ]; //Integrate the accelerations of the domain
    
                    }
    
                    localIndex += _dim * _dim;
    
                }

            }

            error = reconstructedVolume->performVolumeIntegration( dataAtMicroPoints, dataCountAtPoint, integratedValues );

            if ( error ){

                errorOut result = new errorNode( "computeDomainVolumeAverages",
                                                 "Error in the computation of the relative position volume integrals" );
                result->addNext( error );
                return result;

            }

        }

        if ( !_inputProcessor.microBodyForceDefined( ) ){

            integratedValues = vectorTools::appendVectors( { floatVector( integratedValues.begin( ),
                                                                          integratedValues.begin( ) + initialOffset ),
                                                             floatVector( _dim * _dim, 0 ),
                                                             floatVector( integratedValues.begin( ) + initialOffset,
                                                                          integratedValues.end( ) ) } );

        }

        if ( !_inputProcessor.microAccelerationDefined( ) ){

            integratedValues = vectorTools::appendVectors( { integratedValues, floatVector( _dim * _dim, 0 ) } );

        }

        if ( homogenizedBodyForceCouples.find( macroCellID ) == homogenizedBodyForceCouples.end( ) ){

            floatVector tmp( integratedValues.begin( ),
                             integratedValues.begin( ) + initialOffset );
            tmp /= ( homogenizedVolumes[ macroCellID ].back( ) * homogenizedDensities[ macroCellID ].back( ) );
            homogenizedMicroInertias.emplace( macroCellID, tmp );

            tmp = floatVector( integratedValues.begin( ) + initialOffset,
                               integratedValues.begin( ) + initialOffset + _dim * _dim );
            tmp /= ( homogenizedVolumes[ macroCellID ].back( ) * homogenizedDensities[ macroCellID ].back( ) );
            homogenizedBodyForceCouples.emplace( macroCellID, tmp );

            tmp = floatVector( integratedValues.begin( ) + initialOffset + _dim * _dim,
                               integratedValues.begin( ) + initialOffset + 2 * _dim * _dim );
            tmp /= ( homogenizedVolumes[ macroCellID ].back( ) * homogenizedDensities[ macroCellID ].back( ) );
            homogenizedMicroSpinInertias.emplace( macroCellID, tmp );

        }
        else{
            floatVector tmp( integratedValues.begin( ),
                             integratedValues.begin( ) + initialOffset );
            tmp /= ( homogenizedVolumes[ macroCellID ].back( ) * homogenizedDensities[ macroCellID ].back( ) );
            homogenizedMicroInertias[ macroCellID ] = vectorTools::appendVectors( { homogenizedMicroInertias[ macroCellID ],
                                                                                    tmp } );

            tmp = floatVector( integratedValues.begin( ) + initialOffset,
                               integratedValues.begin( ) + initialOffset + _dim * _dim );
            tmp /= ( homogenizedVolumes[ macroCellID ].back( ) * homogenizedDensities[ macroCellID ].back( ) );
            homogenizedBodyForceCouples[ macroCellID ] = vectorTools::appendVectors( { homogenizedBodyForceCouples[ macroCellID ],
                                                                                       tmp } );
            tmp = floatVector( integratedValues.begin( ) + initialOffset + _dim * _dim,
                               integratedValues.begin( ) + initialOffset + 2 * _dim * _dim );
            tmp /= ( homogenizedVolumes[ macroCellID ].back( ) * homogenizedDensities[ macroCellID ].back( ) );
            homogenizedMicroSpinInertias[ macroCellID ] = vectorTools::appendVectors( { homogenizedMicroSpinInertias[ macroCellID ],
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

        //Get the pointers to the values
        const floatVector *macroNodeReferenceLocations = _inputProcessor.getMacroNodeReferencePositions( );
        const floatVector *macroDisplacements = _inputProcessor.getMacroDisplacements( );
        const uIntVector *macroConnectivity = _inputProcessor.getMacroNodeReferenceConnectivity( );
        const uIntVector *macroConnectivityCellIndices = _inputProcessor.getMacroNodeReferenceConnectivityCellIndices( );

        //Form the finite element representation of the macro-scale
        std::unique_ptr< elib::Element > element;
        errorOut error = overlapCoupling::buildMacroDomainElement( macroCellID, *macroNodeReferenceLocations,
                                                                   *macroDisplacements, *macroConnectivity,
                                                                   *macroConnectivityCellIndices, element );

        if ( error ){

            errorOut result = new errorNode( "computeHomogenizedStresses",
                                             "Error in the formation of the finite element representation of the macro-scale" );
            result->addNext( error );
            return result;

        }

        //Get the shape functions at the micro-domain centroids
        floatVector shapefunctionsAtCentersOfMass;
        error = overlapCoupling::computeShapeFunctionsAtPoints( macroCellID, *macroNodeReferenceLocations, *macroDisplacements,
                                                                *macroConnectivity, *macroConnectivityCellIndices,
                                                                homogenizedCentersOfMass[ macroCellID ],
                                                                shapefunctionsAtCentersOfMass );

        if ( error ){

            errorOut result = new errorNode( "computeHomogenizedStresses",
                                             "Error in the computation of the shapefunctions at the micro domain centers of mass for macro cell " + std::to_string( macroCellID ) );
            result->addNext( error );
            return result;

        }

        uIntType nMicroDomains = homogenizedCentersOfMass[ macroCellID ].size( ) / _dim;
        uIntType nMacroCellNodes = element->nodes.size( );

        if ( shapefunctionsAtCentersOfMass.size( ) != nMicroDomains * nMacroCellNodes ){
            
            std::string output;
            output += "The number of shape-function defined is not consistent with the number of micro domains\n";
            output += "and the number of nodes in the macro element for macro-cell " + std::to_string( macroCellID ) + ".\n";
            output += "This is likely because one of the micro-domains center of mass is located outside of the macro cell";

            return new errorNode( "computeHomogenizedStresses", output );

        }

        floatVector linearMomentumRHS( _dim * nMacroCellNodes, 0 );
        floatVector firstMomentRHS( _dim * _dim * nMacroCellNodes, 0 );

        //Project values ot the nodes
        floatVector volumeAtNodes( nMacroCellNodes, 0 );
        floatVector densityAtNodes( nMacroCellNodes, 0 );
        floatMatrix bodyForceAtNodes( nMacroCellNodes, floatVector( _dim, 0 ) );
        floatMatrix accelerationAtNodes( nMacroCellNodes, floatVector( _dim, 0 ) );
        floatMatrix microInertiaAtNodes( nMacroCellNodes, floatVector( _dim * _dim, 0 ) );
        floatMatrix bodyCoupleAtNodes( nMacroCellNodes, floatVector( _dim * _dim, 0 ) );
        floatMatrix microSpinInertiaAtNodes( nMacroCellNodes, floatVector( _dim * _dim, 0 ) );
        floatMatrix symmetricMicroStressAtNodes( nMacroCellNodes, floatVector( _dim * _dim, 0 ) );

        //Add the volume integral components of the right hand side vectors
        //Also project the integral components to the nodes
        for ( unsigned int i = 0; i < nMicroDomains; i++ ){

            floatType density = homogenizedDensities[ macroCellID ][ i ];

            floatType volume = homogenizedVolumes[ macroCellID ][ i ];

            floatVector bodyForce( homogenizedBodyForces[ macroCellID ].begin( ) + _dim * i,
                                   homogenizedBodyForces[ macroCellID ].begin( ) + _dim * ( i + 1 ) );

            floatVector acceleration( homogenizedAccelerations[ macroCellID ].begin( ) + _dim * i,
                                      homogenizedAccelerations[ macroCellID ].begin( ) + _dim * ( i + 1 ) );

            floatVector microInertia( homogenizedMicroInertias[ macroCellID ].begin( ) + _dim * _dim * i,
                                      homogenizedMicroInertias[ macroCellID ].begin( ) + _dim * _dim * ( i + 1 ) );

            floatVector bodyCouple( homogenizedBodyForceCouples[ macroCellID ].begin( ) + _dim * _dim * i,
                                    homogenizedBodyForceCouples[ macroCellID ].begin( ) + _dim * _dim * ( i + 1 ) );

            floatVector microSpinInertia( homogenizedMicroSpinInertias[ macroCellID ].begin( ) + _dim * _dim * i,
                                          homogenizedMicroSpinInertias[ macroCellID ].begin( ) + _dim * _dim * ( i + 1 ) );

            floatVector symmetricMicroStress( homogenizedSymmetricMicroStresses[ macroCellID ].begin( ) + _dim * _dim * i,
                                              homogenizedSymmetricMicroStresses[ macroCellID ].begin( ) + _dim * _dim * ( i + 1 ) );

            floatVector symmetricMicroStress_T( _dim * _dim );

            for ( unsigned int _i = 0; _i < _dim; _i++ ){

                for ( unsigned int _j = 0; _j < _dim; _j++ ){

                    symmetricMicroStress_T[ _dim * _j + _i ] = symmetricMicroStress[ _dim * _i + _j ];

                }

            }

            for ( unsigned int j = 0; j < nMacroCellNodes; j++ ){

                //Get the shapefunction value for the node
                floatType N = shapefunctionsAtCentersOfMass[ nMacroCellNodes * i + j ];

                //Compute the contribution to the node
                floatVector nLinearMomentumRHS = N * density * ( bodyForce - acceleration ) * volume;

                floatVector nFirstMomentRHS = N * ( density * ( bodyCouple - microSpinInertia ) - symmetricMicroStress_T ) * volume;

                //Add the contribution to the overall RHS vectors
                for ( auto it = nLinearMomentumRHS.begin( ); it != nLinearMomentumRHS.end( ); it++ ){

                    linearMomentumRHS[ _dim * j + it - nLinearMomentumRHS.begin( ) ] += *it;

                }

                for ( auto it = nFirstMomentRHS.begin( ); it != nFirstMomentRHS.end( ); it++ ){

                    firstMomentRHS[ _dim * _dim * j + it - nFirstMomentRHS.begin( ) ] += *it;

                }

                //Project values to the nodes
                volumeAtNodes[ j ]               += N * volume;
                densityAtNodes[ j ]              += N * density * volume;
                bodyForceAtNodes[ j ]            += N * density * bodyForce * volume;
                accelerationAtNodes[ j ]         += N * density * acceleration * volume;
                microInertiaAtNodes[ j ]         += N * density * microInertia * volume;
                bodyCoupleAtNodes[ j ]           += N * density * bodyCouple * volume;
                microSpinInertiaAtNodes[ j ]     += N * density * microSpinInertia * volume;
                symmetricMicroStressAtNodes[ j ] += N * symmetricMicroStress * volume;

            }

        }

        //Save the contributions of the body forces and couples to the external force at the nodes
        externalForcesAtNodes.emplace( macroCellID, vectorTools::appendVectors( bodyForceAtNodes ) );
        externalCouplesAtNodes.emplace( macroCellID, vectorTools::appendVectors( bodyCoupleAtNodes ) );

        //De-weight the projected values at the nodes
        for ( unsigned int n = 0; n < nMacroCellNodes; n++ ){

            densityAtNodes[ n ]              /= volumeAtNodes[ n ];
            bodyForceAtNodes[ n ]            /= ( densityAtNodes[ n ] * volumeAtNodes[ n ] );
            accelerationAtNodes[ n ]         /= ( densityAtNodes[ n ] * volumeAtNodes[ n ] );
            microInertiaAtNodes[ n ]         /= ( densityAtNodes[ n ] * volumeAtNodes[ n ] );
            bodyCoupleAtNodes[ n ]           /= ( densityAtNodes[ n ] * volumeAtNodes[ n ] );
            microSpinInertiaAtNodes[ n ]     /= ( densityAtNodes[ n ] * volumeAtNodes[ n ] );
            symmetricMicroStressAtNodes[ n ] /= volumeAtNodes[ n ];

        }

        //Add the surface integral components of the right hand side vectors
        uIntType nMicroSurfaceRegions = homogenizedSurfaceRegionAreas[ macroCellID ].size( );

        //Compute the shape functions at the surface region centers
        floatVector shapefunctionsAtSurfaceRegionCentersOfMass;

        //TODO: The following computation of the shape-function values may cause issues if the re-constructed surface is found
        //      to be outside of the macro Cell's domain. This will likely need to be addressed
        error = overlapCoupling::computeShapeFunctionsAtPoints( macroCellID, *macroNodeReferenceLocations, *macroDisplacements,
                                                                *macroConnectivity, *macroConnectivityCellIndices,
                                                                homogenizedSurfaceRegionCentersOfMass[ macroCellID ],
                                                                shapefunctionsAtSurfaceRegionCentersOfMass );

        if ( error ){

            errorOut result = new errorNode( "computeHomogenizedStresses",
                                             "Error in the computation of the shapefunctions at the micro domain surface region centers of mass for macro cell " + std::to_string( macroCellID ) );
            result->addNext( error );
            return result;

        }

        if ( ( shapefunctionsAtSurfaceRegionCentersOfMass.size( ) / element->nodes.size( ) ) !=
             ( homogenizedSurfaceRegionCentersOfMass[ macroCellID ].size( ) / _dim ) ){

            std::string output = "The number of shape-function defined is not consistent with the number of micro surface regions\n";
            output += "and the number of nodes in the macro element for macro-cell " + std::to_string( macroCellID ) + ".\n";
            output += "This is likely because one of the micro surface regions center of mass is located outside of the macro cell.";
            output += "\nA future workaround would be to project the surfaces back to the nearest point on the macro-element's surface";

            return new errorNode( "computeHomogenizedStresses", output );


        }

        //Add the surface integral components of the right hand side vectors
        for ( unsigned int i = 0; i < nMicroSurfaceRegions; i++ ){

            floatType area = homogenizedSurfaceRegionAreas[ macroCellID ][ i ];

            floatVector traction( homogenizedSurfaceRegionTractions[ macroCellID ].begin( ) + _dim * i,
                                  homogenizedSurfaceRegionTractions[ macroCellID ].begin( ) + _dim * ( i + 1 ) );

            floatVector couple( homogenizedSurfaceRegionCouples[ macroCellID ].begin( ) + _dim * _dim * i,
                                homogenizedSurfaceRegionCouples[ macroCellID ].begin( ) + _dim * _dim * ( i + 1 ) );

            for ( unsigned int j = 0; j < nMacroCellNodes; j++ ){

                //Get the shapefunction value for the surface region
                floatType N = shapefunctionsAtSurfaceRegionCentersOfMass[ nMacroCellNodes * i + j ];

                //Compute the contribution to the node
                floatVector nLinearMomentumRHS = N * traction * area;

                floatVector nFirstMomentRHS = N * couple * area;

                //Add the contribution to the overall RHS vectors
                for ( auto it = nLinearMomentumRHS.begin( ); it != nLinearMomentumRHS.end( ); it++ ){

                    linearMomentumRHS[ _dim * j + it - nLinearMomentumRHS.begin( ) ] += *it;
                    externalForcesAtNodes[ macroCellID ][ _dim * j + it - nLinearMomentumRHS.begin( ) ] += *it;

                }

                for ( auto it = nFirstMomentRHS.begin( ); it != nFirstMomentRHS.end( ); it++ ){

                    firstMomentRHS[ _dim * _dim * j + it - nFirstMomentRHS.begin( ) ] += *it;
                    externalCouplesAtNodes[ macroCellID ][ _dim * _dim * j + it - nFirstMomentRHS.begin( ) ] += *it;

                }

            }

        }

        //Assemble the LHS matrix

        //Loop over the quadrature points adding the contribution of each to the LHS matrix
        //This formulation does not construct the stresses at each of the micro domains but 
        //rather computes the stresses at the macro domain's Gauss points which enables the
        //solution of the macro balance equations.

        floatType Jxw;
        floatVector shapeFunctions;
        floatMatrix dNdx, jacobian;

        std::vector< DOFProjection::T > coefficients;
        coefficients.reserve( ( 2 * _dim * _dim + 3 * _dim * _dim ) * element->nodes.size( ) * element->qrule.size( ) );

        //Quadrature point interpolated values
        floatVector densities( element->qrule.size( ), 0 );
        floatMatrix bodyForces( element->qrule.size( ), floatVector( _dim, 0 ) );
        floatMatrix accelerations( element->qrule.size( ), floatVector( _dim, 0 ) );
        floatMatrix microInertias( element->qrule.size( ), floatVector( _dim * _dim, 0 ) );
        floatMatrix bodyCouples( element->qrule.size( ), floatVector( _dim * _dim, 0 ) );
        floatMatrix microSpinInertias( element->qrule.size( ), floatVector( _dim * _dim, 0 ) );
        floatMatrix symmetricMicroStress( element->qrule.size( ), floatVector( _dim * _dim, 0 ) ); 

        for ( auto qpt = element->qrule.begin( ); qpt != element->qrule.end( ); qpt++ ){

            //Set the index
            uIntType qptIndex = qpt - element->qrule.begin( );

            //Set the column
            uIntType col0 = ( _dim * _dim + _dim * _dim * _dim ) * qptIndex;

            //Get the values of the shape function and the gradients
            error = element->get_shape_functions( qpt->first, shapeFunctions );

            if ( error ){

                errorOut result = new errorNode( "computeHomogenizedStresses",
                                                 "Error in the computation of the shape functions\n" );
                result->addNext( error );
                return result;

            }

            //Get the values of the shape function gradients
            error = element->get_global_shapefunction_gradients( qpt->first, dNdx );

            if ( error ){

                errorOut result = new errorNode( "computeHomogenizedStresses",
                                                 "Error in the computation of the shape function gradients\n" );
                result->addNext( error );
                return result;

            }

            //Get the Jacobian of transformation
            error = element->get_local_gradient( element->nodes, qpt->first, jacobian );

            if ( error ){

                errorOut result = new errorNode( "computeHomogenizedStresses",
                                                 "Error in the computation of the local gradient\n" );
                result->addNext( error );
                return result;

            }

            Jxw = vectorTools::determinant( vectorTools::appendVectors( jacobian ), _dim, _dim ) * qpt->second;

            for ( unsigned int n = 0; n < element->nodes.size( ); n++ ){

                //Set the row
                uIntType row0 = n * ( _dim + _dim * _dim );

                //Add the balance of linear momentum contributions
                for ( unsigned int i = 0; i < _dim; i++ ){

                    for ( unsigned int j = 0; j < _dim; j++ ){

                        coefficients.push_back( DOFProjection::T( row0 + i, col0 + i + _dim * j, dNdx[ n ][ j ] * Jxw ) );

                    }

                }

                //Add the balance of the first moment of momentum contributions
                row0 += _dim;

                //Cauchy stress contribution
                for ( unsigned int i = 0; i < _dim; i++ ){

                    for ( unsigned int j = 0; j < _dim; j++ ){

                        coefficients.push_back( DOFProjection::T( row0 + _dim * j + i, col0 + _dim * i + j, -shapeFunctions[ n ] * Jxw ) );

                    }

                }

                //Higher order stress contribution
                for ( unsigned int i = 0; i < _dim * _dim; i++ ){

                    for ( unsigned int j = 0; j < _dim; j++ ){

                        coefficients.push_back( DOFProjection::T( row0 + i, col0 + _dim * _dim + _dim * _dim * j + i, dNdx[ n ][ j ] * Jxw ) );

                    }

                }

                //Interpolate the nodal values to the quadrature points
                densities[ qptIndex ]            += shapeFunctions[ n ] * densityAtNodes[ n ];
                bodyForces[ qptIndex ]           += shapeFunctions[ n ] * bodyForceAtNodes[ n ];
                accelerations[ qptIndex ]        += shapeFunctions[ n ] * accelerationAtNodes[ n ];
                microInertias[ qptIndex ]        += shapeFunctions[ n ] * microInertiaAtNodes[ n ];
                bodyCouples[ qptIndex ]          += shapeFunctions[ n ] * bodyCoupleAtNodes[ n ];
                microSpinInertias[ qptIndex ]    += shapeFunctions[ n ] * microSpinInertiaAtNodes[ n ];
                symmetricMicroStress[ qptIndex ] += shapeFunctions[ n ] * symmetricMicroStressAtNodes[ n ];

            }

        }

        //Form the left-hand side sparse matrix
        SparseMatrix LHS( ( _dim + _dim * _dim ) * element->nodes.size( ), _dim * _dim * ( 1 + _dim ) * element->qrule.size( ) );
        LHS.setFromTriplets( coefficients.begin( ), coefficients.end( ) );

        //Perform the SVD decomposition
        Eigen::JacobiSVD< Eigen::MatrixXd > svd( LHS.toDense( ), Eigen::ComputeThinU | Eigen::ComputeThinV );
       
        //Compute the threshold for the SVD 
        floatVector logSVec( LHS.rows( ), 0 );

        Eigen::Map< Eigen::MatrixXd > logS( logSVec.data(), logSVec.size(), 1 );

        //Compute the singular values
        logS = svd.singularValues();

        for ( unsigned int i = 0; i < logSVec.size( ); i++ ){

            logSVec[ i ] = std::log10( logSVec[ i ] + _absoluteTolerance );

        }

        //Determine where the "shelf" in the singular values occurs
        std::vector< unsigned int > outliers;

        MADOutlierDetection( logSVec, outliers, 10 );

        if ( outliers.size( ) > 0 ){

            svd.setThreshold( std::max( pow( 10, logSVec[ outliers[ 0 ] ] ), _absoluteTolerance ) );

        }
        else{

            svd.setThreshold( _absoluteTolerance );

        }

        floatVector rhsVec = vectorTools::appendVectors( { linearMomentumRHS, firstMomentRHS } );

        Eigen::Map< Eigen::MatrixXd > RHS( rhsVec.data( ), rhsVec.size( ), 1 ); 

        //Solve for the stresses
        Eigen::MatrixXd x = svd.solve( RHS );

        //Extract the stresses at the evaluation points
        uIntType nCauchy = _dim * _dim;
        uIntType nHigherOrder = _dim * _dim * _dim;

        uIntType nEvaluationPoints = x.size( ) / ( nCauchy + nHigherOrder );

        floatVector cauchyStresses( _dim * _dim * nEvaluationPoints );
        floatVector higherOrderStresses( _dim * _dim * _dim * nEvaluationPoints );

        for ( unsigned int n = 0; n < nEvaluationPoints; n++ ){

            for ( unsigned int i = 0; i < nCauchy; i++ ){

                cauchyStresses[ nCauchy * n + i ] = x( ( nCauchy + nHigherOrder ) * n + i );

            }

            for ( unsigned int i = 0; i < nHigherOrder; i++ ){

                higherOrderStresses[ nHigherOrder * n + i ] = x( ( nCauchy + nHigherOrder ) * n + nCauchy + i );

            }

        }

        quadraturePointCauchyStress.emplace( macroCellID, cauchyStresses );
        quadraturePointHigherOrderStress.emplace( macroCellID, higherOrderStresses );

        //Emplace the values at the quadrature points
        quadraturePointDensities.emplace( macroCellID, densities );
        quadraturePointBodyForce.emplace( macroCellID, vectorTools::appendVectors( bodyForces ) );
        quadraturePointAccelerations.emplace( macroCellID, vectorTools::appendVectors( accelerations ) );
        quadraturePointMicroInertias.emplace( macroCellID, vectorTools::appendVectors( microInertias ) );
        quadraturePointBodyCouples.emplace( macroCellID, vectorTools::appendVectors( bodyCouples ) );
        quadraturePointMicroSpinInertias.emplace( macroCellID, vectorTools::appendVectors( microSpinInertias ) );
        quadraturePointSymmetricMicroStress.emplace( macroCellID, vectorTools::appendVectors( symmetricMicroStress ) );

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

    errorOut MADOutlierDetection( const floatVector &x, uIntVector &outliers, const floatType threshold,
                                  const floatType eps ){
        /*!
         * Detect outliers using median absolute deviation
         * MAD = median ( | X_i - median(X) | )
         *
         * :param const floatVector &x: The x vector to search for outliers
         * :param uIntVector> &outliers: The vector of outliers
         * :param const floatType threshold: The threshold with which to identify an outlier. Defaults to 10.
         * :param const floatType eps: The minimum allowable value for MAD
         */

        floatType median = vectorTools::median(x);

        std::vector< floatType > absDeviations = vectorTools::abs(x - median);

        floatType MAD = vectorTools::median(absDeviations) + eps;

        absDeviations /= MAD;

        outliers.resize(0);
        outliers.reserve(x.size() / 10);

        for ( unsigned int i = 0; i < absDeviations.size( ); i++ ){

            if ( absDeviations[ i ] > threshold ){

                outliers.push_back(i);

            }

        }

        return NULL;
    }

    errorOut formMicromorphicElementMassMatrix( const std::unique_ptr< elib::Element > &element,
                                                const floatVector &degreeOfFreedomValues,
                                                const floatVector &momentOfInertia,
                                                const floatVector &density,
                                                const DOFMap *nodeIDToIndex,
                                                std::vector< DOFProjection::T > &coefficients ){
        /*!
         * Form the micromorphic mass matrix for an element
         *
         * :param const std::unique_ptr< elib::Element > &element: The element to form the mass matrix of
         * :param const floatVector &degreeOfFreedomValues: The degree of freedom values at the element nodes
         * :param const floatVector &momentOfInertia: The moment of inertia in the current configuration
         *     at the quadrature points ordered as [ i1_11, i1_12, i1_13, i1_21, ... , i2_11, i2_12, ... ] 
         *     where the first index is the quadrature point and the second indices are the indices of the
         *     moment of inertia tensor.
         * :param const floatVector &density: The density in the current configuration
         *     at the quadrature points
         * :param const DOFMap &nodeIDToIndex: A map from the node id's to the DOF index
         * :param std::vector< DOFProjection::T > &coefficients: The coefficients of the mass matrix
         */

        //Get the dimension of the element
        uIntType dim = element->nodes[ 0 ].size( );

        const uIntType uSize   = dim;
        const uIntType phiSize = dim * dim;

        //Check that the degree of freedom value vector's length is consistent with the element
        if ( degreeOfFreedomValues.size( ) != ( uSize + phiSize ) * element->nodes.size( ) ){

            return new errorNode( "formMicromorphicElementMassMatrix",
                                  "The degree of freedom vector size is not consistent with the element dimension" );

        }

        if ( momentOfInertia.size( ) != element->qrule.size( ) * phiSize ){

            return new errorNode( "formMicromorphicElementMassMatrix",
                                  "The moment of inertia vector size is not consistent with the quadrature rule and element dimension" );

        }

        if ( density.size( ) != element->qrule.size( ) ){

            return new errorNode( "formMicromorphicElementMassMatrix",
                                  "The density vector size is not consistent with the quadrature rule" );

        }

        if ( element->global_node_ids.size( ) != element->nodes.size( )  ){

            return new errorNode( "formMicromorphicElementMassMatrix",
                                  "The size of the global node id in the element are not the same size as the number of nodes" );

        }

        //Reshape the degree of freedom values to a matrix of values where the rows are the values at the nodes
        floatMatrix reshapedDOFValues = vectorTools::inflate( degreeOfFreedomValues, element->nodes.size( ), uSize + phiSize );

        //Variable initialize
        floatVector shapeFunctions;
        floatVector interpolatedValues, deformationGradient;
        floatVector qptMomentOfInertia;
        floatVector uQpt, XiQpt, invXiQpt, referenceMomentOfInertia, inertiaTerm;
        floatMatrix gradShapeFunctions;

        floatVector eye( dim * dim );
        vectorTools::eye( eye );

        floatType J, Jxw, sFo, sFp;
        uIntType qptIndex, row0, col0;
        errorOut error = NULL;

        //Loop over the quadrature points
        for ( auto qpt = element->qrule.begin( ); qpt != element->qrule.end( ); qpt++ ){

            //Set the quadrature point index
            qptIndex = qpt - element->qrule.begin( );

            //Compute the base micromorphic element terms
            error = computeMicromorphicElementRequiredValues( element, qpt, dim, reshapedDOFValues, true,
                                                              shapeFunctions, gradShapeFunctions,
                                                              deformationGradient, J, Jxw, uQpt, XiQpt );

            if ( error ){

                errorOut result = new errorNode( "formMicromorphicElementInternalForceVector",
                                                 "Error in the computation of the required values for the element" );
                result->addNext( error );
                return result;

            }

            invXiQpt = vectorTools::inverse( XiQpt, dim, dim );

            //Compute the moment of inertia in the reference configuration
            qptMomentOfInertia = floatVector( momentOfInertia.begin( ) + dim * dim * qptIndex,
                                              momentOfInertia.begin( ) + dim * dim * ( qptIndex + 1 ) );

            referenceMomentOfInertia
                = vectorTools::matrixMultiply( vectorTools::matrixMultiply( invXiQpt, qptMomentOfInertia, dim, dim, dim, dim ),
                                               invXiQpt, dim, dim, dim, dim, false, true );

            //Evaluate the integrand term
            inertiaTerm = density[ qptIndex ] * J * referenceMomentOfInertia * Jxw;

            //Add the contrubutions to the mass matrix
            for ( uIntType o = 0; o < shapeFunctions.size( ); o++ ){

                sFo = shapeFunctions[ o ];

                auto gni1 = nodeIDToIndex->find( element->global_node_ids[ o ] );

                if ( gni1 == nodeIDToIndex->end( ) ){

                    return new errorNode( "formMicromorphicElementMassMatrix",
                                          "Node " + std::to_string( element->global_node_ids[ o ] ) + " not found in the ID map" );
                                          

                }

                row0 = ( uSize + phiSize ) * gni1->second;

                for ( uIntType p = 0; p < shapeFunctions.size( ); p++ ){

                    sFp = shapeFunctions[ p ];

                    auto gni2 = nodeIDToIndex->find( element->global_node_ids[ p ] );

                    if ( gni2 == nodeIDToIndex->end( ) ){
    
                        return new errorNode( "formMicromorphicElementMassMatrix",
                                              "Node " + std::to_string( element->global_node_ids[ p ] ) + " not found in the ID map" );
                                              
    
                    }
    
                    col0 = ( uSize + phiSize ) * gni2->second;

                    for ( unsigned int j = 0; j < dim; j++ ){

                        for ( unsigned int k = 0; k < dim; k++ ){

                            coefficients.push_back( DOFProjection::T( row0 + j,
                                                                      col0 + k,
                                                                      eye[ dim * j + k ] * density[ qptIndex ] * J * sFo * sFp * Jxw ) );
    
                            for ( unsigned int K = 0; K < dim; K++ ){
    
                                for ( unsigned int L = 0; L < dim; L++ ){
    
                                    coefficients.push_back( DOFProjection::T( row0 + dim + dim * j + K,
                                                                              col0 + dim + dim * k + L,
                                                                              eye[ dim * j + k ] * sFo * sFp * inertiaTerm[ dim * K + L ] ) );

                                }
    
                            }

                        }

                    }

                }

            }

        }

        return NULL;
    }

    errorOut formMicromorphicElementInternalForceVector( const std::unique_ptr< elib::Element > &element,
                                                         const floatVector &degreeOfFreedomValues,
                                                         const floatVector &cauchyStress,
                                                         const floatVector &symmetricMicroStress,
                                                         const floatVector &higherOrderStress,
                                                         const DOFMap *nodeIDToIndex,
                                                         Eigen::MatrixXd &internalForceVector ){
        /*!
         * Add the contribution of the micromorphic element to the internal force vector
         *
         * :param const std::unique_ptr< elib::Element > &element: The FEA representation of the micromorphic element
         * :param const floatVector &degreeOfFreedomValues: The values of the degrees of freedom at the nodes of the element
         * :param const floatVector &cauchyStress: The cauchy stress at the quadrature points of the element
         *     ( current configuration )
         * :param const floatVector &symmetricMicroStress: The symmetric micro-stress at the quadrature points of the element
         *     ( current configuration )
         * :param const floatVector &higherOrderStress: The higher-order stress at the quadrature points of the element
         *     ( current configuraiton )
         * :param const DOFMap *nodeIDToIndex: The map from the global node ids to the ordering in the internal force vector
         * :param Eigen::MatrixXd &internalForceVector: The internal force vector to be populated.
         */

        //Get the dimension of the element
        uIntType dim = element->nodes[ 0 ].size( );

        const uIntType uSize   = dim;
        const uIntType phiSize = dim * dim;

        if ( dim != 3 ){

            std::string output  = "The dimension of the problem is required to be 3. This only matters ( it is believed )\n";
            output             += "because of dNdX, fint, and cint which are currently consistend with a 3D problem as required\n";
            output             += "by balance_equations.h";
            return new errorNode( "formMicromorphicElementInternalForceVector", output );

        }

        //Check that the degree of freedom value vector's length is consistent with the element
        if ( degreeOfFreedomValues.size( ) != ( uSize + phiSize ) * element->nodes.size( ) ){

            return new errorNode( "formMicromorphicElementInternalForceVector",
                                  "The degree of freedom vector size is not consistent with the element dimension" );

        }

        if ( cauchyStress.size( ) != element->qrule.size( ) * dim * dim ){

            return new errorNode( "formMicromorphicElementInternalForceVector",
                                  "The Cauchy stress vector size is not consistent with the quadrature rule and element dimension" );

        }

        if ( symmetricMicroStress.size( ) != element->qrule.size( ) * dim * dim ){

            return new errorNode( "formMicromorphicElementInternalForceVector",
                                  "The symmetric micro-stress vector size is not consistent with the quadrature rule" );

        }

        if ( higherOrderStress.size( ) != element->qrule.size( ) * dim * dim * dim ){

            return new errorNode( "formMicromoprhicElementInternalForceVector",
                                  "The higher-order stress vector size is not consistent with the quadrature rule" );

        }

        if ( element->global_node_ids.size( ) != element->nodes.size( )  ){

            return new errorNode( "formMicromorphicElementInternalForceVector",
                                  "The size of the global node id in the element are not the same size as the number of nodes" );

        }

        //Reshape the degree of freedom values to a matrix of values where the rows are the values at the nodes
        floatMatrix reshapedDOFValues = vectorTools::inflate( degreeOfFreedomValues, element->nodes.size( ), uSize + phiSize );

        //Initialize variables
        floatType N, J, Jxw;
        floatType dNdX[ 3 ], fint[ 3 ], cint[ 9 ];
        int errorCode;

        floatVector shapeFunctions, deformationGradient, uQpt, XiQpt;
        floatVector cauchyQpt, sQpt, mQpt;
        floatVector pk2Qpt, referenceMicroStressQpt, referenceHigherOrderStressQpt;

        floatMatrix gradShapeFunctions;

        uIntType qptIndex, row0;

        errorOut error = NULL;

        //Loop over the quadrature points
        for ( auto qpt = element->qrule.begin( ); qpt != element->qrule.end( ); qpt++ ){

            //Set the quadrature point index
            qptIndex = qpt - element->qrule.begin( );

            //Evaluate the properties required for the integration
            error = computeMicromorphicElementRequiredValues( element, qpt, dim, reshapedDOFValues, true,
                                                              shapeFunctions, gradShapeFunctions,
                                                              deformationGradient, J, Jxw, uQpt, XiQpt );

            if ( error ){

                errorOut result = new errorNode( "formMicromorphicElementInternalForceVector",
                                                 "Error in the computation of the required values for the element" );
                result->addNext( error );
                return result;

            }

            //Pull back the stresses at the quadrature points to the reference configuration
            cauchyQpt = floatVector( cauchyStress.begin( ) + dim * dim * qptIndex,
                                     cauchyStress.begin( ) + dim * dim * ( qptIndex + 1 ) );

            sQpt = floatVector( symmetricMicroStress.begin( ) + dim * dim * qptIndex,
                                symmetricMicroStress.begin( ) + dim * dim * ( qptIndex + 1 ) );

            mQpt = floatVector( higherOrderStress.begin( ) + dim * dim * dim * qptIndex,
                                higherOrderStress.begin( ) + dim * dim * dim * ( qptIndex + 1 ) );

            //Pull back the Cauchy stress
            error = micromorphicTools::pullBackCauchyStress( cauchyQpt, deformationGradient, pk2Qpt );

            if ( error ){

                errorOut result = new errorNode( "formMicromorphicElementInternalForceVector",
                                                 "Error in the pull-back operation on the Cauchy stress" );
                result->addNext( error );
                return result;

            }

            //Pull back the symmetric micro-stress
            error = micromorphicTools::pullBackMicroStress( sQpt, deformationGradient, referenceMicroStressQpt );

            if ( error ){

                errorOut result = new errorNode( "formMicromorphicElementInternalForceVector",
                                                 "Error in the pull-back operation on the symmetric micro-stress" );
                result->addNext( error );
                return result;

            }

            //Pull back the higher order stress
            error = micromorphicTools::pullBackHigherOrderStress( mQpt, deformationGradient, XiQpt, referenceHigherOrderStressQpt );

            if ( error ){

                errorOut result = new errorNode( "formMicromorphicElementInternalForceVector",
                                                 "Error in the pull-back operation on the higher order stress" );
                result->addNext( error );
                return result;

            }

            //Loop over the nodes
            for ( unsigned int n = 0; n < shapeFunctions.size( ); n++ ){

                //Set the shape function and gradient of the shape function values
                N = shapeFunctions[ n ];

                for ( unsigned int i = 0; i < dim; i++ ){

                    dNdX[ i ] = gradShapeFunctions[ n ][ i ];

                }

                //Compute the terms for the balance of linear momentum
                errorCode = balance_equations::compute_internal_force( dNdX, deformationGradient, pk2Qpt, fint );

                if ( errorCode != 0 ){

                    return new errorNode( "formMicromorphicElementInternalForceVector",
                                          "The internal force term returned an error code: " + std::to_string( errorCode ) );

                }

                //Compute the terms for the balance of first moment of momentum
                errorCode = balance_equations::compute_internal_couple( N, dNdX, deformationGradient, XiQpt,
                                                                        pk2Qpt, referenceMicroStressQpt,
                                                                        referenceHigherOrderStressQpt, cint );

                if ( errorCode != 0 ){

                    return new errorNode( "formMicromorphicElementInternalForceVector",
                                          "The internal couple term returned an error code: " + std::to_string( errorCode ) );

                }

                //Get the initial index
                auto it = nodeIDToIndex->find( element->global_node_ids[ n ] );
                
                if ( it == nodeIDToIndex->end( ) ){

                    return new errorNode( "formMicromorphicElementInternalForceVector",
                                          "The global node id " + std::to_string( element->global_node_ids[ n ] ) +
                                          " is not found in the id to index map" );

                }

                //Set the row index
                row0 = ( uSize + phiSize ) * it->second;

                if ( ( row0 + uSize + phiSize ) > internalForceVector.rows( ) ){

                    return new errorNode( "formMicromorphicElementInternalForceVector",
                                          "The global node id " + std::to_string( element->global_node_ids[ n ] ) +
                                          " has an index ( " + std::to_string( it->second ) + " ) which results in a index larger than" +
                                          " the internal force vector size ( " + std::to_string( internalForceVector.rows( ) ) + ")" );

                }

                for ( unsigned int i = 0; i < dim; i++ ){

                    internalForceVector( row0 + i, 0 ) -= fint[ i ] * Jxw;

                }

                for ( unsigned int i = 0; i < ( dim * dim ); i++ ){

                    internalForceVector( row0 + i + dim, 0 ) -= cint[ i ] * Jxw;

                }

            }

        }

        return NULL;
    }

    errorOut computeMicromorphicElementRequiredValues( const std::unique_ptr< elib::Element > &element,
                                                       const elib::quadrature_rule::iterator &qpt,
                                                       const uIntType dim,
                                                       const floatMatrix &reshapedDOFValues,
                                                       const bool useReference,
                                                       floatVector &shapeFunctions,
                                                       floatMatrix &gradShapeFunctions,
                                                       floatVector &deformationGradient,
                                                       floatType &J, floatType &Jxw,
                                                       floatVector &uQpt, floatVector &XiQpt ){
        /*!
         * Compute the required values for the integration of a micromorphic element
         */

        //Initialize internal variables
        floatVector interpolatedValues, eye( dim * dim );
        vectorTools::eye( eye );

        floatMatrix jacobian;

        //Evaluate the shape function values
        errorOut error = element->get_shape_functions( qpt->first, shapeFunctions );

        if ( error ){

            errorOut result = new errorNode( "computeMicromorphicElementRequiredValues",
                                             "Error in the computation of the shape functions" );
            result->addNext( error );
            return result;

        }

        //Evaluate the gradients of the shape functions
        error = element->get_global_shapefunction_gradients( qpt->first, gradShapeFunctions, useReference );

        if ( error ){

            errorOut result = new errorNode( "computeMicromorphicElementRequiredValues",
                                             "Error in the computation of the shape function gradients" );
            result->addNext( error );
            return result;

        }

        //Get the deformation gradient between the reference and current configurations
        error = element->get_jacobian( qpt->first, element->reference_nodes, jacobian );

        if ( error ){

            errorOut result = new errorNode( "computeMicromorphicElementRequiredValues",
                                             "Error in the computation of the jacobian" );
            result->addNext( error );
            return result;

        }

        deformationGradient = vectorTools::appendVectors( jacobian );
        J = vectorTools::determinant( vectorTools::appendVectors( jacobian ), dim, dim );

        //Get the Jacobian between the local and reference configurations
        if ( useReference ){

            error = element->get_local_gradient( element->reference_nodes, qpt->first, jacobian );

        }
        else{

            error = element->get_local_gradient( element->nodes, qpt->first, jacobian );

        }

        if ( error ){

            errorOut result = new errorNode( "computeMicromorphicElementRequiredValues",
                                             "Error in the computation of the local gradient\n" );
            result->addNext( error );
            return result;

        }

        Jxw = vectorTools::determinant( vectorTools::appendVectors( jacobian ), dim, dim ) * qpt->second;

        //Interpolate the DOF nodes to the node
        error = element->interpolate( reshapedDOFValues, qpt->first, interpolatedValues );

        if ( error ){

            errorOut result = new errorNode( "computeMicromorphicElementRequiredValues",
                                             "Error in the interpolation of the degree of freedom values" );
            result->addNext( error );
            return result;

        }

        if ( interpolatedValues.size( ) < ( dim + dim * dim ) ){

            std::string output = "The interpolated values shape is not consistent with the required dimension for the displacement ";
            output            += "and micro-displacement interpolation";
            return new errorNode( "computeMicromorphicElementRequiredValues", output );

        }

        uQpt = floatVector( interpolatedValues.begin( ), interpolatedValues.begin( ) + dim );
        XiQpt = eye + floatVector( interpolatedValues.begin( ) + dim, interpolatedValues.begin( ) + ( dim + dim * dim ) );

        return NULL;
    }

    errorOut overlapCoupling::assembleHomogenizedExternalForceVector( ){
        /*!
         * Assemble the homogenized external force vector. This vector doesn't have
         * any scaling from the coefficients but is just the raw vector.
         */

        //Loop over the elements in the external force container
        std::unique_ptr< elib::Element > element;
        errorOut error = NULL;
        const DOFMap *nodeIDToIndex = _inputProcessor.getMacroGlobalToLocalDOFMap( );

        //Resize the output vector
        homogenizedFEXT = Eigen::MatrixXd::Zero( ( _dim + _dim * _dim ) * nodeIDToIndex->size( ), 1 );

        //Collect all of the cells in the overlapping domain
        const uIntVector *freeMacroCellIds = _inputProcessor.getFreeMacroCellIds( );
        const uIntVector *ghostMacroCellIds = _inputProcessor.getGhostMacroCellIds( );

        uIntVector macroCellIDVector( freeMacroCellIds->begin( ), freeMacroCellIds->end( ) );
        macroCellIDVector = vectorTools::appendVectors( { macroCellIDVector, *ghostMacroCellIds } );

        for ( auto macroCellID = macroCellIDVector.begin( ); macroCellID != macroCellIDVector.end( ); macroCellID++ ){

            //Make sure that the macroCellID is stored in the external force vector
            if ( externalForcesAtNodes.find( *macroCellID ) == externalForcesAtNodes.end( ) ){

                return new errorNode( "assembleHomogenizedExternalForceVector",
                                      "Macro cell ID " + std::to_string( *macroCellID ) +
                                      " not found in external forces at nodes." );

            }

            //Make sure that the macroCellID is stored in the external couple vector
            if ( externalCouplesAtNodes.find( *macroCellID ) == externalCouplesAtNodes.end( ) ){

                return new errorNode( "assembleHomogenizedExternalForceVector",
                                      "Macro cell ID " + std::to_string( *macroCellID ) +
                                      " not found in external couples at nodes." );

            }

            //Form the macro element
            error = buildMacroDomainElement( *macroCellID,
                                             *_inputProcessor.getMacroNodeReferencePositions( ),
                                             *_inputProcessor.getMacroDisplacements( ),
                                             *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                             *_inputProcessor.getMacroNodeReferenceConnectivityCellIndices( ),
                                             element );

            if ( error ){

                errorOut result = new errorNode( "assembleHomogenizedExternalForceVector",
                                                 "Error in the construction of the macro domain element for macro cell " +
                                                 std::to_string( *macroCellID ) );
                result->addNext( error );
                return result;

            }

            //Loop over the external node ids

            for ( auto globalNodeID = element->global_node_ids.begin( );
                      globalNodeID != element->global_node_ids.end( );
                      globalNodeID++ ){

                //Get the element nodal index
                uIntType elementNodeIndex = globalNodeID - element->global_node_ids.begin( );

                auto index = nodeIDToIndex->find( *globalNodeID );

                if ( index == nodeIDToIndex->end( ) ){

                    return new errorNode( "assembleHomogenizedExternalForceVector",
                                          "Macro global node " + std::to_string( *globalNodeID ) +
                                          " not found in the id to index map" );

                }

                for ( unsigned int i = 0; i < _dim; i++ ){

                    homogenizedFEXT( ( _dim + _dim * _dim ) * index->second + i, 0 )
                        += externalForcesAtNodes[ *macroCellID ][ _dim * elementNodeIndex + i ];

                }

                for ( unsigned int i = 0; i < _dim * _dim; i++ ){

                    homogenizedFEXT( ( _dim + _dim * _dim ) * index->second + i + _dim, 0 )
                        += externalCouplesAtNodes[ *macroCellID ][ _dim * _dim * elementNodeIndex + i ];

                }

            }

        }

        return NULL;

    }

    errorOut overlapCoupling::assembleHomogenizedInternalForceVector( ){
        /*!
         * Assemble the homogenized internal force vector.
         *
         */

        //Loop over the elements in the external force container
        std::unique_ptr< elib::Element > element;
        errorOut error = NULL;
        const DOFMap *nodeIDToIndex = _inputProcessor.getMacroGlobalToLocalDOFMap( );

        //Set the number of displacement degrees of freedom
        unsigned int nMacroDispDOF = _dim + _dim * _dim;

        //Get the free and ghost macro node ids
        const uIntVector *freeMacroNodeIds = _inputProcessor.getFreeMacroNodeIds( );
        const uIntVector *ghostMacroNodeIds = _inputProcessor.getGhostMacroNodeIds( );

        //Assemble the free macro node degree of freedom vector
        floatVector freeMacroDisplacements( nMacroDispDOF * freeMacroNodeIds->size( ) );

        const floatVector *macroDispDOFVector = _inputProcessor.getMacroDispDOFVector( );
        for ( auto it = freeMacroNodeIds->begin( ); it != freeMacroNodeIds->end( ); it++ ){

            auto map = _inputProcessor.getMacroGlobalToLocalDOFMap( )->find( *it );

            if ( map == _inputProcessor.getMacroGlobalToLocalDOFMap( )->end( ) ){

                return new errorNode( "assembleHomogenizedInternalForceVector",
                                      "Global degree of freedom '" + std::to_string( *it ) + "' not found in degree of freedom map" );

            }

            //Set the macro displacements
            for ( unsigned int i = 0; i < nMacroDispDOF; i++ ){

                freeMacroDisplacements[ nMacroDispDOF * ( map->second ) + i ]
                    = ( *macroDispDOFVector )[ nMacroDispDOF * map->first + i ];

            }

        }

        //Resize the output vector
        homogenizedFINT = Eigen::MatrixXd::Zero( nMacroDispDOF * nodeIDToIndex->size( ), 1 );

        //Loop over the macro cells
        floatVector elementDOFVector;

        //Collect all of the cells in the overlapping domain
        const uIntVector *freeMacroCellIds = _inputProcessor.getFreeMacroCellIds( );
        const uIntVector *ghostMacroCellIds = _inputProcessor.getGhostMacroCellIds( );

        uIntVector macroCellIDVector( freeMacroCellIds->begin( ), freeMacroCellIds->end( ) );
        macroCellIDVector = vectorTools::appendVectors( { macroCellIDVector, *ghostMacroCellIds } );

        for ( auto macroCellID = macroCellIDVector.begin( ); macroCellID != macroCellIDVector.end( ); macroCellID++ ){

            //Form the macro element
            error = buildMacroDomainElement( *macroCellID,
                                             *_inputProcessor.getMacroNodeReferencePositions( ),
                                             *_inputProcessor.getMacroDisplacements( ),
                                             *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                             *_inputProcessor.getMacroNodeReferenceConnectivityCellIndices( ),
                                             element );

            if ( error ){

                errorOut result = new errorNode( "assembleHomogenizedInternalForceVector",
                                                 "Error in the construction of the macro domain element for macro cell " +
                                                 std::to_string( *macroCellID ) );
                result->addNext( error );
                return result;

            }
            
            //Extract the homogenized degrees of freedom for the element

            floatVector nodeDOF;
            elementDOFVector.clear( );
            elementDOFVector.reserve( nMacroDispDOF * element->nodes.size( ) );

            for ( auto nodeID = element->global_node_ids.begin( ); nodeID != element->global_node_ids.end( ); nodeID++ ){

                //Find the node in the global to local map
                auto index = nodeIDToIndex->find( *nodeID );

                if ( index == nodeIDToIndex->end( ) ){

                    return new errorNode( "assembleHomogenizedInternalForceVector",
                                          "Macro-scale node with global id " + std::to_string( *nodeID ) +
                                          " is not found in the global ID to the local index" );

                }

                //Check if the macro-node is free
                auto freeNodeID = std::find( freeMacroNodeIds->begin( ), freeMacroNodeIds->end( ), *nodeID );

                if ( freeNodeID != freeMacroNodeIds->end( ) ){

                    //Extract the degrees of freedom from the free node
                    nodeDOF = floatVector( freeMacroDisplacements.begin( ) + ( nMacroDispDOF ) * ( index->second ),
                                           freeMacroDisplacements.begin( ) + ( nMacroDispDOF ) * ( index->second + 1 ) );

                }
                else{

                    auto ghostNodeID = std::find( ghostMacroNodeIds->begin( ), ghostMacroNodeIds->end( ), *nodeID );

                    if ( ghostNodeID != ghostMacroNodeIds->end( ) ){

                        //Extract the degrees of freedom from the ghost node
                        nodeDOF = floatVector( _projected_ghost_macro_displacement.begin( ) + nMacroDispDOF * ( index->second - _inputProcessor.getFreeMacroNodeIds( )->size( ) ),
                                               _projected_ghost_macro_displacement.begin( ) + nMacroDispDOF * ( index->second + 1 - _inputProcessor.getFreeMacroNodeIds( )->size( ) ) );

                    }
                    else{

                        return new errorNode( "assembleHomogenizedInternalForceVector",
                                              "The macro node " + std::to_string( *nodeID ) +
                                              " is not found in either the ghost or free macro node IDs" );

                    }

                }

                //Save the element DOF vector
                for ( unsigned int i = 0; i < nodeDOF.size( ); i++ ){

                    elementDOFVector.push_back( nodeDOF[ i ] );

                }

            }

            error = formMicromorphicElementInternalForceVector( element, elementDOFVector,
                                                                quadraturePointCauchyStress[ *macroCellID ],
                                                                quadraturePointSymmetricMicroStress[ *macroCellID ],
                                                                quadraturePointHigherOrderStress[ *macroCellID ],
                                                                nodeIDToIndex, homogenizedFINT );

            if ( error ){

                errorOut result = new errorNode( "assembleHomogenizedInternalForceVector",
                                                 "Error in the assembly of the terms of the internal force vector for element " +
                                                 std::to_string( *macroCellID ) );
                result->addNext( error );
                return result;

            }

        }

        return NULL;

    }

    errorOut overlapCoupling::assembleHomogenizedMassMatrix( ){
        /*!
         * Assemble the homogenized mass matrix.
         */

        //Loop over the elements in the external force container
        std::unique_ptr< elib::Element > element;
        errorOut error = NULL;
        const DOFMap *nodeIDToIndex = _inputProcessor.getMacroGlobalToLocalDOFMap( );

        //Set the number of displacement degrees of freedom
        unsigned int nMacroDispDOF = _dim + _dim * _dim;

        //Get the free and ghost macro node ids
        const uIntVector *freeMacroNodeIds = _inputProcessor.getFreeMacroNodeIds( );
        const uIntVector *ghostMacroNodeIds = _inputProcessor.getGhostMacroNodeIds( );

        //Assemble the free macro node degree of freedom vector
        floatVector freeMacroDisplacements( nMacroDispDOF * freeMacroNodeIds->size( ) );

        const floatVector *macroDispDOFVector = _inputProcessor.getMacroDispDOFVector( );
        for ( auto it = freeMacroNodeIds->begin( ); it != freeMacroNodeIds->end( ); it++ ){

            auto map = _inputProcessor.getMacroGlobalToLocalDOFMap( )->find( *it );

            if ( map == _inputProcessor.getMacroGlobalToLocalDOFMap( )->end( ) ){

                return new errorNode( "assembleHomogenizedInternalForceVector",
                                      "Global degree of freedom '" + std::to_string( *it ) + "' not found in degree of freedom map" );

            }

            //Set the macro displacements
            for ( unsigned int i = 0; i < nMacroDispDOF; i++ ){

                freeMacroDisplacements[ nMacroDispDOF * ( map->second ) + i ]
                    = ( *macroDispDOFVector )[ nMacroDispDOF * map->first + i ];

            }

        }

        //Loop over the macro cells adding the contributions to the mass matrix
        floatVector elementDOFVector;

        //Collect all of the cells in the overlapping domain
        const uIntVector *freeMacroCellIds = _inputProcessor.getFreeMacroCellIds( );
        const uIntVector *ghostMacroCellIds = _inputProcessor.getGhostMacroCellIds( );

        uIntVector macroCellIDVector( freeMacroCellIds->begin( ), freeMacroCellIds->end( ) );
        macroCellIDVector = vectorTools::appendVectors( { macroCellIDVector, *ghostMacroCellIds } );

        std::vector< DOFProjection::T > coefficients;

        uIntType numCoefficients = 0;
        for ( auto it = externalForcesAtNodes.begin( ); it != externalForcesAtNodes.end( ); it++ ){

            //Get the number of quadrature points in the element
            uIntType elementQuadraturePointCount = quadraturePointDensities[ it->first ].size( );

            //Get the number of nodes in the element
            uIntType elementNodeCount = it->second.size( ) / _dim;
            numCoefficients += elementQuadraturePointCount * elementNodeCount * elementNodeCount * _dim * _dim * ( 1 + _dim * _dim );

        }

        coefficients.reserve( numCoefficients );

        for ( auto macroCellID = macroCellIDVector.begin( ); macroCellID != macroCellIDVector.end( ); macroCellID++ ){

            //Form the macro element
            error = buildMacroDomainElement( *macroCellID,
                                             *_inputProcessor.getMacroNodeReferencePositions( ),
                                             *_inputProcessor.getMacroDisplacements( ),
                                             *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                             *_inputProcessor.getMacroNodeReferenceConnectivityCellIndices( ),
                                             element );

            if ( error ){

                errorOut result = new errorNode( "assembleHomogenizedInternalForceVector",
                                                 "Error in the construction of the macro domain element for macro cell " +
                                                 std::to_string( *macroCellID ) );
                result->addNext( error );
                return result;

            }

            //Extract the homogenized degrees of freedom for the element
            floatVector nodeDOF;
            elementDOFVector.clear( );
            elementDOFVector.reserve( nMacroDispDOF * element->nodes.size( ) );

            for ( auto nodeID = element->global_node_ids.begin( ); nodeID != element->global_node_ids.end( ); nodeID++ ){

                //Find the node in the global to local map
                auto index = nodeIDToIndex->find( *nodeID );

                if ( index == nodeIDToIndex->end( ) ){

                    return new errorNode( "assembleHomogenizedInternalForceVector",
                                          "Macro-scale node with global id " + std::to_string( *nodeID ) +
                                          " is not found in the global ID to the local index" );

                }

                //Check if the macro-node is free
                auto freeNodeID = std::find( freeMacroNodeIds->begin( ), freeMacroNodeIds->end( ), *nodeID );

                if ( freeNodeID != freeMacroNodeIds->end( ) ){

                    //Extract the degrees of freedom from the free node
                    nodeDOF = floatVector( freeMacroDisplacements.begin( ) + ( nMacroDispDOF ) * ( index->second ),
                                           freeMacroDisplacements.begin( ) + ( nMacroDispDOF ) * ( index->second + 1 ) );

                }
                else{

                    auto ghostNodeID = std::find( ghostMacroNodeIds->begin( ), ghostMacroNodeIds->end( ), *nodeID );

                    if ( ghostNodeID != ghostMacroNodeIds->end( ) ){

                        //Extract the degrees of freedom from the ghost node
                        nodeDOF = floatVector( _projected_ghost_macro_displacement.begin( ) + nMacroDispDOF * ( index->second - _inputProcessor.getFreeMacroNodeIds( )->size( ) ),
                                               _projected_ghost_macro_displacement.begin( ) + nMacroDispDOF * ( index->second + 1 - _inputProcessor.getFreeMacroNodeIds( )->size( ) ) );

                    }
                    else{

                        return new errorNode( "assembleHomogenizedInternalForceVector",
                                              "The macro node " + std::to_string( *nodeID ) +
                                              " is not found in either the ghost or free macro node IDs" );

                    }

                }

                //Save the element DOF vector
                for ( unsigned int i = 0; i < nodeDOF.size( ); i++ ){

                    elementDOFVector.push_back( nodeDOF[ i ] );

                }

            }

            //Add the contributions to the mass matrix
            error = formMicromorphicElementMassMatrix( element, elementDOFVector,
                                                       quadraturePointMicroInertias[ *macroCellID ],
                                                       quadraturePointDensities[ *macroCellID ],
                                                       nodeIDToIndex, coefficients );

        }

        homogenizedMassMatrix = SparseMatrix( ( _dim + _dim * _dim ) * nodeIDToIndex->size( ),
                                              ( _dim + _dim * _dim ) * nodeIDToIndex->size( ) );
        homogenizedMassMatrix.setFromTriplets( coefficients.begin( ), coefficients.end( ) );

        return NULL;

    }

    errorOut overlapCoupling::assembleHomogenizedMatricesAndVectors( ){
        /*!
         * Assemble the homogenized mass matrices and force vectors
         */

        errorOut error = assembleHomogenizedExternalForceVector( );

        if ( error ){

            errorOut result = new errorNode( "assembleHomogenizedMatricesAndVectors",
                                             "Error in the construction of the homogenized external force vector" );
            result->addNext( error );
            return result;

        }

        error = assembleHomogenizedInternalForceVector( );

        if ( error ){

            errorOut result = new errorNode( "assembleHomogenizedMatricesAndVectors",
                                             "Error in the construction of the homogenized internal force vector" );
            result->addNext( error );
            return result;

        }

        error = assembleHomogenizedMassMatrix( );

        if ( error ){

            errorOut result = new errorNode( "assembleHomogenizedMatricesAndVectors",
                                             "Error in the construction of the homogenized mass matrix" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut overlapCoupling::assembleFreeMicromorphicMassMatrix( ){
        /*!
         * Assemble the micromorphic mass matrix for the free micromorphic domains.
         * It is noted that while we are only processing the free micromorphic domains,
         * there WILL be ghost nodes because even though a micromorphic domain may be
         * marked as free, some of the nodes can ( and will in the current framework )
         * be ghost.
         */

        //Set the displacement degrees of freedom for the element
        const unsigned int nMacroDOF = _dim + _dim * _dim;

        //Get the micromorphic densities in the reference configuration
        const std::unordered_map< uIntType, std::string > *macroReferenceDensityTypes = _inputProcessor.getMacroReferenceDensityTypes( );
        const std::unordered_map< uIntType, floatVector > *macroReferenceDensities = _inputProcessor.getMacroReferenceDensities( );

        //Get the micromorphic moments of inertia in the reference configuration
        const std::unordered_map< uIntType, std::string > *macroReferenceMomentOfInertiaTypes
            = _inputProcessor.getMacroReferenceMomentOfInertiaTypes( );
        const std::unordered_map< uIntType, floatVector > *macroReferenceMomentsOfInertia
            = _inputProcessor.getMacroReferenceMomentsOfInertia( );

        //Initialize the coefficients vector
        std::vector< DOFProjection::T > coefficients;

        uIntType numCoefficients = 0;
        for ( auto it = externalForcesAtNodes.begin( ); it != externalForcesAtNodes.end( ); it++ ){

            //Get the number of quadrature points in the element
            uIntType elementQuadraturePointCount = quadraturePointDensities[ it->first ].size( );

            //Get the number of nodes in the element
            uIntType elementNodeCount = it->second.size( ) / _dim;
            numCoefficients += elementQuadraturePointCount * elementNodeCount * elementNodeCount * _dim * _dim * ( 1 + _dim * _dim );

        }

        coefficients.reserve( numCoefficients );

        //Loop over the free micromorphic elements
        for ( auto macroCellID  = _inputProcessor.getFreeMacroCellIds( )->begin( );
                   macroCellID != _inputProcessor.getFreeMacroCellIds( )->end( );
                   macroCellID++ ){

            //Construct the macro-domain element
            std::unique_ptr< elib::Element > element;
            errorOut error = buildMacroDomainElement( *macroCellID,
                                                      *_inputProcessor.getMacroNodeReferencePositions( ),
                                                      *_inputProcessor.getMacroDisplacements( ),
                                                      *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                      *_inputProcessor.getMacroNodeReferenceConnectivityCellIndices( ),
                                                      element );

            if ( error ){

                errorOut result = new errorNode( "assembleFreeMicromorphicMassMatrix",
                                                 "Error in the construction of the macro element " + std::to_string( *macroCellID ) );
                result->addNext( error );
                return result;

            }

            //Get the degree of freedom values for the free macro-domain element
            const floatVector* macroDispDOFVector = _inputProcessor.getMacroDispDOFVector( );
            floatVector elementDOFVector( 0 );

            for ( auto nodeID  = element->global_node_ids.begin( );
                       nodeID != element->global_node_ids.end( );
                       nodeID++ ){

                elementDOFVector
                    = vectorTools::appendVectors( { elementDOFVector, 
                                                    floatVector( macroDispDOFVector->begin( ) + nMacroDOF * *nodeID,
                                                                 macroDispDOFVector->begin( ) + nMacroDOF * ( ( *nodeID ) + 1 ) )
                                                  } );

            }

            //Extract the density and moment of inertia in the reference configuration
            auto densityType = macroReferenceDensityTypes->find( *macroCellID );

            if ( densityType == macroReferenceDensityTypes->end( ) ){

                return new errorNode( "assembleFreeMicromorphicMassMatrix",
                                      "The macro cell with ID " + std::to_string( *macroCellID ) +
                                      " was not found in the density type map" );

            }

            auto momentOfInertiaType = macroReferenceDensityTypes->find( *macroCellID );

            if ( momentOfInertiaType == macroReferenceMomentOfInertiaTypes->end( ) ){

                return new errorNode( "assembleFreeMicromorphicMassMatrix",
                                      "The macro cell with ID " + std::to_string( *macroCellID ) +
                                      " was not found in the moment of inertia type map" );

            }

            if ( densityType->second.compare( "constant" ) != 0 ){

                return new errorNode( "assembleFreeMicromorphicMassMatrix",
                                      "Only constant densities for the macro-scale are allowed currently. This is not true for macro cell ID " + std::to_string( *macroCellID ) );

            }

            if ( momentOfInertiaType->second.compare( "constant" ) != 0 ){

                return new errorNode( "assembleFreeMicromorphicMassMatrix",
                                      "Only constant moments of inertia for the macro-scale are allowed currently. This is not true for macro cell ID " + std::to_string( *macroCellID ) );

            }

            auto macroDensities = macroReferenceDensities->find( *macroCellID );

            if ( macroDensities == macroReferenceDensities->end( ) ){

                return new errorNode( "assembleFreeMicromorphicMassMatrix",
                                      "Macro cell ID " + std::to_string( *macroCellID ) +
                                      " is not in the macro reference density map" );

            }

            if ( macroDensities->second.size( ) != 1 ){

                return new errorNode( "assembleFreeMicromorphicMassMatrix",
                                      "The macro densities for macro cell " + std::to_string( *macroCellID ) +
                                      "Define " + std::to_string( macroDensities->second.size( ) ) +
                                      " values when only 1 can be defined" );

            }

            auto macroMomentsOfInertia = macroReferenceMomentsOfInertia->find( *macroCellID );

            if ( macroMomentsOfInertia == macroReferenceMomentsOfInertia->end( ) ){

                return new errorNode( "assembleFreeMicromorphicMassMatrix",
                                      "Macro cell ID " + std::to_string( *macroCellID ) +
                                      " is not in the macro reference moments of inertia map" );

            }

            if ( macroMomentsOfInertia->second.size( ) != _dim * _dim ){

                return new errorNode( "assembleFreeMicromorphicMassMatrix",
                                      "The macro moments of inertia for macro cell " + std::to_string( *macroCellID ) +
                                      "Define " + std::to_string( macroDensities->second.size( ) ) +
                                      " values when only " + std::to_string( _dim * _dim ) + " can be defined" );

            }


            floatVector densities( element->qrule.size( ), macroDensities->second[ 0 ] );

            floatVector momentsOfInertia
                = vectorTools::appendVectors( floatMatrix( element->qrule.size( ), macroMomentsOfInertia->second ) );

            //Construct the kinetic energy partitioning coefficient at each quadrature point
            floatVector res;
            error = constructKineticEnergyPartitioningCoefficient( *macroCellID, element, res );

            if ( error ){

                errorOut result = new errorNode( "assembleFreeMicromorphicMassMatrix",
                                                 "Error in the construction of the kinetic energy partitoning coefficient for macro cell " + std::to_string( *macroCellID ) );
                result->addNext( error );
                return result;

            }

            for ( auto re = res.begin( ); re != res.end( ); re++ ){

                unsigned int re_indx = re - res.begin( );

                densities[ re_indx ] *= *re;

                for ( auto mOI  = momentsOfInertia.begin( ) + _dim * _dim * re_indx;
                           mOI != momentsOfInertia.begin( ) + _dim * _dim * ( re_indx + 1 );
                           mOI++ ){

                    unsigned int mOI_indx = mOI - momentsOfInertia.begin( );
                    momentsOfInertia[ mOI_indx ] *= *re;

                }

            }

            error = formMicromorphicElementMassMatrix( element, elementDOFVector, momentsOfInertia, densities, 
                                                       _inputProcessor.getMacroGlobalToLocalDOFMap( ), coefficients );

            if ( error ){

                std::string outstr  = "Error in the construction of the contributions of the macro element to ";
                            outstr += "the free micromorphic mass matrix";

                errorOut result = new errorNode( "assembleFreeMicromorphicMassMatrix", outstr );
                result->addNext( error );
                return result;

            }


        }

        const DOFMap *nodeIDToIndex = _inputProcessor.getMacroGlobalToLocalDOFMap( );
        freeMicromorphicMassMatrix = SparseMatrix( nMacroDOF * nodeIDToIndex->size( ), nMacroDOF * nodeIDToIndex->size( ) );
        freeMicromorphicMassMatrix.setFromTriplets( coefficients.begin( ), coefficients.end( ) );

        return NULL;
    }

    errorOut overlapCoupling::assembleCouplingMassAndDampingMatrices( ){
        /*!
         * Assemble the mass matrix for the coupling equations
         *
         * It is assumed that ghost nodes cannot also be free or non-overlapped.
         * It is assumed that free nodes cannot also be non-overlapped
         * If a node is defined in multiple locations ( i.e. it is on the boundary ) then:
         * For the micro-scale the order of preference is:
         *     non-overlapped > free > ghost ( we assume the open window DNS is the best representation of the PDE )
         * For the macro-scale the order of preference is:
         *     ghost > free > non-overlapped ( we assume the micro-scale is the best representation of the PDE )
         * It is also assumed that the micro-scale can be represented by a diagonal mass matrix while the 
         * consistent mass matrix is used at the macro-scale.
         */

        //Set the number of displacement degrees of freedom
        uIntType nMacroDispDOF = _dim + _dim * _dim;

        //Get the configuration of the coupling
        const YAML::Node config = _inputProcessor.getCouplingInitialization( );
        floatType rhat = config[ "kinetic_energy_weighting_factor" ].as< floatType >( );
        floatType qhat = config[ "potential_energy_weighting_factor" ].as< floatType >( );
        floatType aQ = config[ "micro_proportionality_coefficient" ].as< floatType >( );
        floatType aD = config[ "macro_proportionality_coefficient" ].as< floatType >( );

        //Get the micro densities and volumes
        const floatVector *microVolumes   = _inputProcessor.getMicroVolumes( );
        const floatVector *microDensities = _inputProcessor.getMicroVolumes( );

        //Get the global to local micro node mapping
        const DOFMap *microGlobalToLocalDOFMap = _inputProcessor.getMicroGlobalToLocalDOFMap( );

        //Get the IDs of the ghost and free micro nodes
        const uIntVector *ghostMicroNodeIDs = _inputProcessor.getGhostMicroNodeIds( );
        const uIntVector *freeMicroNodeIDs = _inputProcessor.getFreeMicroNodeIds( );

        //Get the IDs of the ghost and free macro nodes
        const uIntVector *ghostMacroNodeIDs = _inputProcessor.getGhostMacroNodeIds( );
        const uIntVector *freeMacroNodeIDs = _inputProcessor.getFreeMacroNodeIds( );

        //Determine the offset of the free micro nodes
        const uIntType nFreeMicroNodes = freeMicroNodeIDs->size( );

        //Determine the offset of the free macro nodes
        const uIntType nFreeMacroNodes = freeMacroNodeIDs->size( );
        const uIntType nGhostMacroNodes = ghostMacroNodeIDs->size( );
                
        //Get the micro mass vectors
        floatVector ghostMicroMasses( ghostMicroNodeIDs->size( ), 0 );
        floatVector freeMicroMasses( ghostMicroNodeIDs->size( ), 0 );

        //Assemble the free micro mass vector
        for ( auto microID = freeMicroNodeIDs->begin( ); microID != freeMicroNodeIDs->end( ); microID++ ){

            auto localMicroNodeIDMap = microGlobalToLocalDOFMap->find( *microID );

            if ( localMicroNodeIDMap == microGlobalToLocalDOFMap->end( ) ){

                return new errorNode( "assembleMacroMassAndDampingMatrices",
                                      "Free micro node: " + std::to_string( *microID ) + " not found in global to local map\n" );

            }

            freeMicroMasses[ localMicroNodeIDMap->second ] = ( *microVolumes )[ *microID ] * ( *microDensities )[ *microID ];

        }
        //Assemble the ghost micro mass vector
        for ( auto microID = ghostMicroNodeIDs->begin( ); microID != ghostMicroNodeIDs->end( ); microID++ ){

            auto localMicroNodeIDMap = microGlobalToLocalDOFMap->find( *microID );

            if ( localMicroNodeIDMap == microGlobalToLocalDOFMap->end( ) ){

                return new errorNode( "assembleMacroMassAndDampingMatrices",
                                      "Ghost micro node: " + std::to_string( *microID ) + " not found in global to local map\n" );

            }

            ghostMicroMasses[ localMicroNodeIDMap->second - nFreeMicroNodes ]
                = ( *microVolumes )[ *microID ] * ( *microDensities )[ *microID ];

        }

        //Assemble the mass sub-matrices
        std::vector< DOFProjection::T > c1;
        std::vector< DOFProjection::T > c2;

        c1.reserve( _dim * ghostMicroMasses.size( ) ); 
        c2.reserve( _dim * freeMicroMasses.size( ) ); 

        uIntType mIndex = 0;

        for ( auto m = ghostMicroMasses.begin( ); m != ghostMicroMasses.end( ); m++, mIndex++ ){

            for ( unsigned int i = 0; i < _dim; i++ ){

                c1.push_back( DOFProjection::T( _dim * mIndex + i, _dim * mIndex + i, ( 1 - rhat ) * ( *m ) ) );

            }

        }

        mIndex = 0;

        for ( auto m = freeMicroMasses.begin( ); m != freeMicroMasses.end( ); m++, mIndex++ ){

            for ( unsigned int i = 0; i < _dim; i++ ){

                c2.push_back( DOFProjection::T( _dim * mIndex + i, _dim * mIndex + i, ( 1 - rhat ) * ( *m ) ) );

            }

        }

        SparseMatrix MQ( _dim * freeMicroMasses.size( ), _dim * freeMicroMasses.size( ) );
        MQ.setFromTriplets( c2.begin( ), c2.end( ) );
        std::cout << "MQ rows x cols: " << MQ.rows( ) << " x " << MQ.cols( ) << "\n";

        SparseMatrix MQhat( _dim * ghostMicroMasses.size( ), _dim * ghostMicroMasses.size( ) );
        MQhat.setFromTriplets( c1.begin( ), c1.end( ) );
        std::cout << "MQhat rows x cols: " << MQhat.rows( ) << " x " << MQhat.cols( ) << "\n";

        std::cout << "homogenizedMassMatrix " << homogenizedMassMatrix.rows( ) << " x " << homogenizedMassMatrix.cols( ) << "\n";
        std::cout << "freeMicromorphicMassMatrix " << freeMicromorphicMassMatrix.rows( ) << " x " << freeMicromorphicMassMatrix.cols( ) << "\n";
        SparseMatrix MTildeDBreve = rhat * homogenizedMassMatrix + freeMicromorphicMassMatrix;
        //Note: kinetic partitioning coefficient applied when the matrix was formed
        std::cout << "MTildeDBreve " << MTildeDBreve.rows( ) << " x " << MTildeDBreve.cols( ) << "\n";
    
        std::cout << nMacroDispDOF * nFreeMacroNodes << "\n";
        SparseMatrix MD    = MTildeDBreve.block( 0, 0, nMacroDispDOF * nFreeMacroNodes, nMacroDispDOF * nFreeMacroNodes );
        std::cout << "MD " << MD.rows( ) << " x " << MD.cols( ) << "\n";
        SparseMatrix MDhat = MTildeDBreve.block( nMacroDispDOF * nFreeMacroNodes, nMacroDispDOF * nFreeMacroNodes,
                                                 nMacroDispDOF * nGhostMacroNodes, nMacroDispDOF * nGhostMacroNodes );
        std::cout << "MDhat " << MDhat.rows( ) << " x " << MDhat.cols( ) << "\n";

        //Due to the restrictions in listed in the comment at the beginning of the function, MBar from Regueiro 2012
        //is an empty matrix as we only handle the coupling domain here.
       
    
        //Assemble Mass matrices for the micro projection equation

        if ( config[ "projection_type" ].as< std::string >( ).compare( "l2_projection" ) == 0 ){

            auto MQQ  = MQ + _L2_BQhatQ.transpose( ) * MQhat * _L2_BQhatQ + _L2_BDhatQ.transpose( ) * MDhat * _L2_BDhatQ;

            auto MQD = _L2_BQhatQ.transpose( ) * MQhat * _L2_BQhatD + _L2_BDhatQ.transpose( ) * MDhat * _L2_BDhatD;
    
            //Assemble Mass matrices for the macro projection equation
            
            auto MDQ = _L2_BQhatD.transpose( ) * MQhat * _L2_BQhatD + _L2_BDhatD.transpose( ) * MDhat * _L2_BDhatQ;

            auto MDD = MD + _L2_BQhatD.transpose( ) * MQhat * _L2_BQhatD + _L2_BDhatD.transpose( ) * MDhat * _L2_BDhatD;
    
            //Assemble the damping matrices for the micro projection equation
            auto CQQ = aQ * MQ + aQ * _L2_BQhatQ.transpose( ) * MQhat * _L2_BQhatQ + aD * _L2_BDhatQ.transpose( ) * MDhat * _L2_BDhatQ;

            auto CQD = aQ * _L2_BQhatQ.transpose( ) * MQhat * _L2_BQhatD;
    
            //Assemble the damping matrices for the macro projection equation
            auto CDQ = aQ * _L2_BQhatD.transpose( ) * MQhat * _L2_BQhatQ;
            auto CDD = aD * MD + aQ * _L2_BQhatD.transpose( ) * MQhat * _L2_BQhatD;

        }
        else if ( config[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){

            std::cout << "MQQ:\n";
            SparseMatrix MQQ  =  MQ + _DP_BQhatQ.transpose( ) * MQhat * _DP_BQhatQ + _DP_BDhatQ.transpose( ) * MDhat * _DP_BDhatQ;
            std::cout << "MQQ nonZeros: " << MQQ.nonZeros( ) << "\n";

            std::cout << "MQD:\n";
            SparseMatrix MQD = _DP_BQhatQ.transpose( ) * MQhat * _DP_BQhatD + _DP_BDhatQ.transpose( ) * MDhat * _DP_BDhatD; //TODO: Verify error in Reguiero 2012 for first term second projection matrix
    
            //Assemble Mass matrices for the macro projection equation
            
            std::cout << "MDQ:\n";
            SparseMatrix MDQ = _DP_BQhatD.transpose( ) * MQhat * _DP_BQhatQ + _DP_BDhatD.transpose( ) * MDhat * _DP_BDhatQ; //TODO: Verify error in Reguiero 2012 for first term second projection matrix

            std::cout << "MDD:\n";
            SparseMatrix MDD = MD + _DP_BQhatD.transpose( ) * MQhat * _DP_BQhatD + _DP_BDhatD.transpose( ) * MDhat * _DP_BDhatD;
    
            //Assemble the damping matrices for the micro projection equation
            std::cout << "CQQ\n";
            SparseMatrix CQQ = aQ * MQ + aQ * _DP_BQhatQ.transpose( ) * MQhat * _DP_BQhatQ +
                               aD * _DP_BDhatQ.transpose( ) * MDhat * _DP_BDhatQ;

            std::cout << "CQD\n";
            SparseMatrix CQD = aQ * _DP_BQhatQ.transpose( ) * MQhat * _DP_BQhatD;
    
            //Assemble the damping matrices for the macro projection equation
            std::cout << "CDQ\n";
            SparseMatrix CDQ = aQ * _DP_BQhatD.transpose( ) * MQhat * _DP_BQhatQ;
            std::cout << "CDD\n";
            SparseMatrix CDD = aD * MD + aQ * _DP_BQhatD.transpose( ) * MQhat * _DP_BQhatD;

        }
        else{

            return new errorNode( "assembleMacroMassAndDampingMatrices",
                                  "The projection type " + config[ "projection_type" ].as< std::string >( ) +
                                  " is not recognized" );

        }

        //Assemble the full mass matrix
//        MASS =; 

        //Assemble the full damping matrix
//        DAMPING =;

        return new errorNode( "assembleMacroMassMatrix", "Not implemented" );
    }

    errorOut overlapCoupling::constructKineticEnergyPartitioningCoefficient( const uIntType &macroCellID,
                                                                             const std::unique_ptr< elib::Element > &element,
                                                                             floatVector &res ){
        /*!
         * Construct the kinetic energy partitioning coefficient
         *
         * :param const uIntType &uIntType: The macro cell's ID number
         * :param const std::unique_ptr< elib::Element > &element: The FEA representation of the macro-scale element
         * :param floatVector &res: The collection of re values at each element quadrature point
         */

        const YAML::Node config = _inputProcessor.getCouplingInitialization( );

        std::string strategy = config[ "kinetic_energy_partitioning_coefficient" ][ "type" ].as< std::string >( );

        if ( strategy.compare( "volume_fraction" ) == 0 ){

            //Compute the volume of the element
            
            floatType elementVolume = 0;
            floatMatrix jacobian;
            
            for ( auto qpt = element->qrule.begin( ); qpt != element->qrule.end( ); qpt++ ){

                //Get the Jacobian of transformation
                errorOut error = element->get_local_gradient( element->nodes, qpt->first, jacobian );

                if ( error ){
 
                    errorOut result = new errorNode( "computeHomogenizedStresses",
                                                     "Error in the computation of the local gradient\n" );
                    result->addNext( error );
                    return result;

                }

                elementVolume += vectorTools::determinant( vectorTools::appendVectors( jacobian ), _dim, _dim ) * qpt->second;

            }

            auto microDomainVolumes = homogenizedVolumes.find( macroCellID );
            if ( microDomainVolumes == homogenizedVolumes.end( ) ){

                return new errorNode( "constructKineticEnergyPartitioningCoefficient",
                                      "The macro cell " + std::to_string( macroCellID ) +
                                      " is not found in the homogenized volumes map" );

            }

            floatType microDomainVolume = 0;

            for ( auto microVolume  = microDomainVolumes->second.begin( );
                       microVolume != microDomainVolumes->second.end( );
                       microVolume++ ){

                microDomainVolume += *microVolume;

            }

            if ( elementVolume < _absoluteTolerance ){

                res = floatVector( element->qrule.size( ), 0 );

            }
            else{

                res = floatVector( element->qrule.size( ), std::fmax( ( elementVolume - microDomainVolume ) / elementVolume, 0 ) );

            }

        }
        else{

            return new errorNode( "constructKineticEnergyPartitioningCoefficient", "Configuration strategy " + strategy + " not recognized" );

        }

        return NULL;

    }

}
