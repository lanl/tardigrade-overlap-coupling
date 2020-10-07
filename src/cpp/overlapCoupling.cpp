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

    errorOut overlapCoupling::processLastIncrements( ){
        /*!
         * Process the final increments of the macro and micro scales
         */

        uIntType numMicroIncrements, numMacroIncrements;

        //Get the number of micro increments
        errorOut error = _inputProcessor._microscale->getNumIncrements( numMicroIncrements );

        if ( error ){

            errorOut result = new errorNode( "processLastIncrements",
                                             "Error in getting the number of micro increments" );
            result->addNext( error );
            return result;

        }

        //Get the number of macro increments
        error = _inputProcessor._macroscale->getNumIncrements( numMacroIncrements );

        if ( error ){

            errorOut result = new errorNode( "processLastIncrements",
                                             "Error in getting the number of macro increments" );
            result->addNext( error );
            return result;

        }

        //Process the final increments
        error = processIncrement( numMicroIncrements - 1, numMacroIncrements - 1 );

        if ( error ){

            std::string outstr = "Error in processing the increments\n";
            outstr += "    macro increment: " + std::to_string( numMacroIncrements - 1 ) + "\n";
            outstr += "    micro increment: " + std::to_string( numMicroIncrements - 1 );
            errorOut result = new errorNode( "processLastIncrements", outstr );
            result->addNext( error );
            return result;

        }

        return NULL;

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

#ifdef TESTACCESS

        _test_initial_projected_ghost_micro_displacement = _projected_ghost_micro_displacement;
        _test_initial_projected_ghost_macro_displacement = _projected_ghost_macro_displacement;

#endif

        //Homogenize the material properties at the micro-scale to the macro-scale
        std::cout << "homogenizing material properties\n";
        error = homogenizeMicroScale( microIncrement );

        if ( error ){

            errorOut result = new errorNode( "processIncrement", "Error in the homogenization of the micro-scale to the macro-scale" );
            result->addNext( error );
            return result;

        }

        return NULL; //REMOVE THIS

        YAML::Node couplingConfiguration = _inputProcessor.getCouplingInitialization( );

        if ( !couplingConfiguration[ "update_displacement" ].IsScalar( ) ){

            //Assemble the mass matrix for the free micromorphic domians
            std::cout << "assembling the free micromorphic mass matrix\n";
            error = assembleFreeMicromorphicMassMatrix( );
    
            if ( error ){
    
                errorOut result = new errorNode( "processIncrement", "Error in the assembly of the mass matrix for the free macro domains" );
                result->addNext( error );
                return result;
    
            }

            //Assemble the coupling mass and damping matrices
            std::cout << "assembling the coupling mass and damping matrices\n";
            error = assembleCouplingMassAndDampingMatrices( );
    
            if ( error ){
    
                errorOut result = new errorNode( "processIncrement", "Error in the construction of the coupling mass and damping matrices" );
                result->addNext( error );
                return result;
    
            }
    
            //Assemble the coupling force vector
            std::cout << "assembling the coupling force vector\n";
            error = assembleCouplingForceVector( );
    
            if ( error ){
    
                errorOut result = new errorNode( "processIncrement", "Error in the construction of the coupling force vector" );
                result->addNext( error );
                return result;
    
            }
    
            //Solve for the free displacements
            std::cout << "solving for the free displacements\n";
            error = solveFreeDisplacement( true );
    
            if ( error ){
    
                errorOut result = new errorNode( "processIncrement", "Error when solving for the free displacements" );
                result->addNext( error );
                return result;
    
            }

        }

        if ( !couplingConfiguration [ "output_homogenized_response" ].IsScalar( ) ){

            //Output the homogenized material response to a data file
            std::cout << "outputting the homogenized response\n";
            error = outputHomogenizedResponse( );
            if ( error ){

                errorOut result = new errorNode( "processIncrement", "Error when writing the homogenized response out to file" );
                result->addNext( error );
                return result;

            }

        }

        if ( !couplingConfiguration[ "output_updated_dof" ].IsScalar( ) ){

            //Output the updated dof values to a data file
            std::cout << "writing the updated DOF to file\n";
            error = overlapCoupling::writeUpdatedDOFToFile( );

            if ( error ){

                errorOut result = new errorNode( "processIncrement", "Error when writing the updated dof information to file" );
                result->addNext( error );
                return result;

            }

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

            if ( error ){

                errorOut result = new errorNode( "initializeCoupling", "Error in setting the initial reference state" );
                result->addNext( error );
                return result;

            }

            bool save_reference_positions = false;

            if ( couplingInitialization[ "projection_type" ].as< std::string >( ).compare(  "direct_projection" ) == 0 ){

                save_reference_positions = true;

            }

            if ( save_reference_positions ){

                _macroReferencePositions.reserve( _inputProcessor.getMacroNodeReferencePositions( )->size( ) );

                for ( auto n = _inputProcessor.getMacroNodeReferencePositions( )->begin( );
                           n != _inputProcessor.getMacroNodeReferencePositions( )->end( );
                           n++ ){

                    auto m = _inputProcessor.getMacroDisplacements( )->find( n->first );

                    if ( m == _inputProcessor.getMacroDisplacements( )->end( ) ){

                        return new errorNode( "initializeCoupling", "Macro node " + std::to_string( n->first ) + " not found in the macro displacements map. Fatal error in the input processor" );

                    }

                    _macroReferencePositions.emplace( n->first, n->second + m->second );

                }

                _microReferencePositions.reserve( _inputProcessor.getMicroNodeReferencePositions( )->size( ) );

                for ( auto n = _inputProcessor.getMicroNodeReferencePositions( )->begin( );
                           n != _inputProcessor.getMicroNodeReferencePositions( )->end( );
                           n++ ){

                    auto m = _inputProcessor.getMicroDisplacements( )->find( n->first );

                    if ( m == _inputProcessor.getMicroDisplacements( )->end( ) ){

                        return new errorNode( "initializeCoupling", "Micro node " + std::to_string( n->first ) + " not found in the micro displacements map. Fatal error in the input processor" );

                    }

                    _microReferencePositions.emplace( n->first, n->second + m->second );

                }

            }

            if ( !couplingInitialization[ "output_homogenized_response" ].IsScalar( ) ){

                error = writeReferenceMeshDataToFile( 0 );

            }

        }
        else if ( couplingInitialization[ "type" ].as< std::string >( ).compare( "from_file" ) == 0 ){

            error = extractProjectionMatricesFromFile( );

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

        //Save the reference state if required
        if ( _inputProcessor.outputReferenceInformation( ) ){

            error = outputReferenceInformation( );

            if ( error ){

                errorOut result = new errorNode( "initializeCoupling", "Error in the output of the reference information" );
                result->addNext( error );
                return result;

            }

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
            errorOut result = new errorNode( "setReferenceStateFromIncrement", "Error in initialization of the input processor" );
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

        //Get the domains encompassed by the macro cells
        const std::unordered_map< uIntType, stringVector > *macroCellToMicroDomainMap = _inputProcessor.getMacroCellToDomainMap( );

//        //Set the output vector sizes
//        unsigned int numFreeMicroDomains  = _inputProcessor.getFreeMicroDomainNames( )->size( );
//        unsigned int numGhostMicroDomains = _inputProcessor.getGhostMicroDomainNames( )->size( );

        //Set the reference micro domain mass vector sizes
//        _referenceFreeMicroDomainMasses = floatVector( numFreeMicroDomains, 0 );
//        _referenceGhostMicroDomainMasses = floatVector( numGhostMicroDomains, 0 );

        //Set the reference micro domain center of mass vector sizes
//        _referenceFreeMicroDomainCentersOfMass = floatVector( _dim * numFreeMicroDomains, 0 );
//        _referenceGhostMicroDomainCentersOfMass = floatVector( _dim * numGhostMicroDomains, 0 );

//        if ( _inputProcessor.getCouplingInitialization( )[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){
//
//            _macroNodeProjectedMass
//                = floatVector( _inputProcessor.getMacroGlobalToLocalDOFMap( )->size( ), 0 );
//
//            _macroNodeProjectedMassMomentOfInertia
//                = floatVector( _dim * _dim * _inputProcessor.getMacroGlobalToLocalDOFMap( )->size( ), 0 );
//
//            _macroNodeMassRelativePositionConstant
//                = floatVector( _dim * _inputProcessor.getMacroGlobalToLocalDOFMap( )->size( ), 0 );
//
//        }

        //Loop over the free macro-scale cells
        unsigned int cellIndex;
        unsigned int nMicroDomains;

        uIntVector macroNodes;

        std::unordered_map< uIntType, floatVector > domainReferenceXiVectors;
        floatVector domainCenterOfMassShapeFunctionValues;
        std::unordered_map< uIntType, floatVector > domainMicroPositionShapeFunctionValues;

        _referenceGhostMicroDomainMasses.clear( );
        _referenceFreeMicroDomainMasses.clear( );
        _referenceGhostMicroDomainCentersOfMass.clear( );
        _referenceFreeMicroDomainCentersOfMass.clear( );
        _referenceCellDomainCenterOfMassShapefunctions.clear( );

        _homogenizationMatrix_initialized = false;
        bool centerOfMassN_initialized = false;

        for ( auto cellID  = freeMacroCellIDs->begin( );
                   cellID != freeMacroCellIDs->end( );
                   cellID++ ){

                //Set the index
                cellIndex = cellID - freeMacroCellIDs->begin( );

                //Set the number of micro-domains encompassed by the cell
                auto microDomains = macroCellToMicroDomainMap->find( *cellID );

                if ( microDomains == macroCellToMicroDomainMap->end( ) ){

                    return new errorNode( "setReferenceStateFromIncrement",
                                          "Macro cell " + std::to_string( *cellID ) + " not found in the macro cell to micro domain map" );

                }

                //Get the macro-node set
                error = _inputProcessor._macroscale->getSubDomainNodes( macroIncrement, ( *freeMacroDomainNames )[ cellIndex ], macroNodes );

                if ( error ){

                    errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                                     "Error in extracting the free macro-node set" );
                    result->addNext( error );
                    return result;

                }

                //Loop over the micro-domains
                domainFloatMap       domainMass;
                domainFloatVectorMap domainCentersOfMass;
                domainFloatVectorMap domainMomentsOfInertia;

                if ( _referenceCellDomainCenterOfMassShapefunctions.find( *cellID ) !=
                     _referenceCellDomainCenterOfMassShapefunctions.end( ) ){

                    return new errorNode( "setReferenceStateFromIncrement",
                                          "Macro cell " + std::to_string( *cellID ) + " was found twice in the reference cell domain center of mass shapefunctions map" );

                }
                else{

                    domainFloatVectorMap tmpDomainFloatVectorMap;
                    _referenceCellDomainCenterOfMassShapefunctions.emplace( *cellID, tmpDomainFloatVectorMap );

                }

#ifdef TESTACCESS
                domainFloatMap _test_floatMap;
                domainFloatVectorMap _test_floatVectorMap;
                std::unordered_map< std::string, std::unordered_map < uIntType, floatVector > > _test_XiMap;

                _test_domainMass.emplace( *cellID, _test_floatMap );
                _test_domainCOM.emplace( *cellID, _test_floatVectorMap );
                _test_domainXi.emplace( *cellID, _test_XiMap );

                if ( _inputProcessor.getCouplingInitialization( )[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){
                    std::unordered_map< std::string, std::unordered_map< uIntType, floatVector > > _test_domainIntFloatVectorMap;
                    _test_domainMUP.emplace( *cellID, _test_domainIntFloatVectorMap );

                }
#endif

                for ( auto domain  = microDomains->second.begin( ); domain != microDomains->second.end( ); domain++ ){

                    error = processDomainReference( microIncrement, *domain,
                                                    *cellID, macroNodes,
                                                    domainMass,
                                                    domainCentersOfMass,
                                                    domainMomentsOfInertia,
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

                _referenceGhostMicroDomainMasses.emplace( *cellID, domainMass );
                _referenceGhostMicroDomainCentersOfMass.emplace( *cellID, domainCentersOfMass );
                _referenceGhostMicroDomainMomentsOfInertia.emplace( *cellID, domainMomentsOfInertia );

                SparseMatrix domainCOMN;
                error = DOFProjection::constructCellCenterOfMassInterpolationMatrixContribution( 1, *cellID, macroNodes,
                                                                                                 *macroCellToMicroDomainMap,
                                                                                                 _referenceCellDomainCenterOfMassShapefunctions,
                                                                                                 *_inputProcessor.getMacroGlobalToLocalDOFMap( ),
                                                                                                 *_inputProcessor.getMicroDomainIDMap( ),
                                                                                                 domainCOMN );

                if ( error ){

                    errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                                     "Error in forming the contribution of macro element " + std::to_string( *cellID ) +
                                                     " to the center of mass shapefunction matrix" );
                    result->addNext( error );
                    return result;

                }

                if ( centerOfMassN_initialized ){

                    _centerOfMassN += domainCOMN;

                }
                else{

                    _centerOfMassN = domainCOMN;
                    centerOfMassN_initialized = true;

                } 
        }

        //Loop over the ghost macro-scale cells

        for ( auto cellID  = ghostMacroCellIDs->begin( );
                   cellID != ghostMacroCellIDs->end( );
                   cellID++ ){

                //Set the index
                cellIndex = cellID - ghostMacroCellIDs->begin( );

                //Set the number of micro-domains encompassed by the cell
                auto microDomains = macroCellToMicroDomainMap->find( *cellID );

                if ( microDomains == macroCellToMicroDomainMap->end( ) ){

                    return new errorNode( "setReferenceStateFromIncrement",
                                          "Macro cell " + std::to_string( *cellID ) + " not found in the macro cell to micro domain map" );

                }

                //Get the macro-node set
                error = _inputProcessor._macroscale->getSubDomainNodes( macroIncrement, ( *ghostMacroDomainNames )[ cellIndex ], macroNodes );

                if ( error ){

                    errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                                     "Error in extracting the ghost macro-node set" );
                    result->addNext( error );
                    return result;

                }

                if ( _referenceCellDomainCenterOfMassShapefunctions.find( *cellID ) !=
                     _referenceCellDomainCenterOfMassShapefunctions.end( ) ){

                    return new errorNode( "setReferenceStateFromIncrement",
                                          "Macro cell " + std::to_string( *cellID ) + " was found twice in the reference cell domain center of mass shapefunctions map" );

                }
                else{

                    domainFloatVectorMap tmpDomainFloatVectorMap;
                    _referenceCellDomainCenterOfMassShapefunctions.emplace( *cellID, tmpDomainFloatVectorMap );

                }

                //Loop over the micro-domains
                domainFloatMap       domainMass;
                domainFloatVectorMap domainCentersOfMass;
                domainFloatVectorMap domainMomentsOfInertia;

                for ( auto domain  = microDomains->second.begin( ); domain != microDomains->second.end( ); domain++ ){

                    error = processDomainReference( microIncrement, *domain,
                                                    *cellID, macroNodes,
                                                    domainMass, domainCentersOfMass,
                                                    domainMomentsOfInertia,
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

                _referenceFreeMicroDomainMasses.emplace( *cellID, domainMass );
                _referenceFreeMicroDomainCentersOfMass.emplace( *cellID, domainCentersOfMass );
                _referenceFreeMicroDomainMomentsOfInertia.emplace( *cellID, domainMomentsOfInertia );

                SparseMatrix domainCOMN;
                error = DOFProjection::constructCellCenterOfMassInterpolationMatrixContribution( 1, *cellID, macroNodes,
                                                                                                 *macroCellToMicroDomainMap,
                                                                                                 _referenceCellDomainCenterOfMassShapefunctions,
                                                                                                 *_inputProcessor.getMacroGlobalToLocalDOFMap( ),
                                                                                                 *_inputProcessor.getMicroDomainIDMap( ),
                                                                                                 domainCOMN );

                if ( error ){

                    errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                                     "Error in forming the contribution of macro element " + std::to_string( *cellID ) +
                                                     " to the center of mass shapefunction matrix" );
                    result->addNext( error );
                    return result;

                }

                if ( centerOfMassN_initialized ){

                    _centerOfMassN += domainCOMN;

                }
                else{

                    _centerOfMassN = domainCOMN;
                    centerOfMassN_initialized = true;

                } 
        }

        //Compress the shape-function matrix
        _N.makeCompressed( );

        //Compute the center of mass to domain projector
        error = DOFProjection::formMoorePenrosePseudoInverse( _centerOfMassN.toDense( ), _centerOfMassProjector );

        if ( error ){

            errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                             "Error in the formation of the center of mass to macro node projector" );
            result->addNext( error );
            return result;

        }

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
        else if ( config[ "projection_type" ].as< std::string >( ).compare( "averaged_l2_projection" ) == 0 ){

            errorOut error = formAveragedL2Projectors( );

            if ( error ){

                errorOut result = new errorNode( "formTheProjectors",
                                                 "Error in the formation of the averaged L2 projectors" );
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
        std::cout << "PERFORMING SVD DECOMPOSITION OF NQDhat\n";
        SparseMatrix NQDhat = _N.block( 0, nFreeMacroDOF, nFreeMicroDOF, nGhostMacroDOF );
        NQDhat.makeCompressed( );
        errorOut error = DOFProjection::formMoorePenrosePseudoInverse( NQDhat.toDense( ), _L2_BDhatQ );

        if ( error ){

            errorOut result = new errorNode( "formL2Projectors", "Error in solving for _L2_BDhatQ" );
            result->addNext( error );
            return result;

        }

        _L2_BDhatD = -_L2_BDhatQ * _N.topLeftCorner( nFreeMicroDOF, nFreeMacroDOF );

        _L2_BQhatQ = _N.bottomRightCorner( nGhostMicroDOF, nGhostMacroDOF ) * _L2_BDhatQ;
        _L2_BQhatD = _N.bottomLeftCorner( nGhostMicroDOF, nFreeMacroDOF )
                   + _N.bottomRightCorner( nGhostMicroDOF, nGhostMacroDOF ) * _L2_BDhatD;

        return NULL;
    }

    errorOut overlapCoupling::formAveragedL2Projectors( ){
        /*!
         * Form the projectors using the averaged micro domain values at the centers of mass.
         *
         * This is the currently recommended projection method
         */

        //Form the micro to macro projection matrix
        //Note: This is the FULL matrix not just the part we need for the overlap coupling.
        //      We will select sub-blocks of it to use in the actual projectors. This may or
        //      may not be correct.

        Eigen::MatrixXd microMacroProjector;
        SparseMatrix S, T;

        uIntType nMicroDOF = _dim;
        uIntType nMacroDOF = _dim + _dim * _dim;

        errorOut error;

        std::cerr << "ASSEMBLING MICRO-TO-MACRO PROJECTOR\n";
        for ( uIntType i = 0; i < nMacroDOF; i++ ){

            error = DOFProjection::formDomainSelectionMatrix( i, nMacroDOF, *_inputProcessor.getMicroDomainIDMap( ), S );

            if ( error ){

                errorOut result = new errorNode( "formAveragedL2Projectors", "Error in the formation of the selection matrix" );
                result->addNext( error );
                return result;

            }

            error = DOFProjection::formMacroNodeExpansionMatrix( i, nMacroDOF, *_inputProcessor.getMacroGlobalToLocalDOFMap( ), T );

            if ( error ){

                errorOut result = new errorNode( "formAveragedL2Projectors", "Error in the formation of the expansion matrix" );
                result->addNext( error );
                return result;

            }

            if ( i == 0 ){

                microMacroProjector = T * _centerOfMassProjector * S;

            }
            else{

                microMacroProjector += T * _centerOfMassProjector * S;

            }

        }

        //Add the homogenization matrix
        microMacroProjector *= _homogenizationMatrix;

        //Set the DOF type sizes
        uIntType nFreeMacroDOF  = nMacroDOF * _inputProcessor.getFreeMacroNodeIds( )->size( );        
        uIntType nGhostMacroDOF = nMacroDOF * _inputProcessor.getGhostMacroNodeIds( )->size( );

        uIntType nFreeMicroDOF  = nMicroDOF * _inputProcessor.getFreeMicroNodeIds( )->size( );
        uIntType nGhostMicroDOF = nMicroDOF * _inputProcessor.getGhostMicroNodeIds( )->size( );

        //Compute the projectors
        _L2_BDhatQ = microMacroProjector.bottomLeftCorner( nGhostMacroDOF, nFreeMicroDOF );

        _L2_BDhatD = -_L2_BDhatQ * _N.topLeftCorner( nFreeMicroDOF, nFreeMacroDOF );

        _L2_BQhatQ = _N.bottomRightCorner( nGhostMicroDOF, nGhostMacroDOF ) * _L2_BDhatQ;
        _L2_BQhatD = _N.bottomLeftCorner( nGhostMicroDOF, nFreeMacroDOF )
                   + _N.bottomRightCorner( nGhostMicroDOF, nGhostMacroDOF ) * _L2_BDhatD;

        return NULL;
    }

    errorOut overlapCoupling::formDirectProjectionProjectors( const unsigned int &microIncrement, const unsigned int &macroIncrement ){
        /*!
         * Form the projectors if the direct projection is to be used
         *
         * :param const unsigned int &microIncrement: The micro increment at which to form the projectors
         * :param const unsigned int &macroIncrement: The macro increment at which to form the projectors
         */

        return new errorNode( "formDirectProjectionProjectors",
                              "This subroutine, and the routines it calls, requires extensive re-working to obtain the expected results. The method is not recommended so this has not been done yet." );

        //Get the ghost macro cell IDs
        const uIntVector *ghostMacroCellIDs = _inputProcessor.getGhostMacroCellIds( );
        const std::unordered_map< uIntType, stringVector > *macroCellToMicroDomainMap = _inputProcessor.getMacroCellToDomainMap( );

        const stringVector *ghostMacroDomainNames = _inputProcessor.getGhostMacroDomainNames( );
        const stringVector *freeMicroDomainNames = _inputProcessor.getFreeMicroDomainNames( );

        unsigned int cellIndex;
        unsigned int nMicroDomains;

        errorOut error = NULL;

        uIntVector macroNodes;

        //Form the projector from the free micro-scale to the ghost macro-scale
        for ( auto cellID  = ghostMacroCellIDs->begin( );
                   cellID != ghostMacroCellIDs->end( );
                   cellID++ ){

            //Set the index
            cellIndex = cellID - ghostMacroCellIDs->begin( );

            //Set the number of micro-domains encompassed by the cell
            auto microDomains = macroCellToMicroDomainMap->find( *cellID );

            if ( microDomains == macroCellToMicroDomainMap->end( ) ){

                return new errorNode( "setReferenceStateFromIncrement",
                                      "Macro cell " + std::to_string( *cellID ) + " not found in the macro cell to micro domain map" );

            }

            //Get the macro-node set
            error = _inputProcessor._macroscale->getSubDomainNodes( macroIncrement, ( *ghostMacroDomainNames )[ cellIndex ], macroNodes );

            if ( error ){

                errorOut result = new errorNode( "formDirectProjectionProjectors",
                                                 "Error in extracting the ghost macro-node set" );
                result->addNext( error );
                return result;

            }

            //Loop over the free micro-domains
            for ( auto domain  = microDomains->second.begin( ); domain != microDomains->second.end( ); domain++ ){

                error = addDomainContributionToDirectFreeMicroToGhostMacroProjector( cellIndex, *cellID, microIncrement, 
                                                                                     *domain, macroNodes );

                if ( error ){

                    errorOut result = new errorNode( "formDirectProjectionProjectors",
                                                     "Error in processing free micro-scale domain '" + *domain + "' for a ghost macro domain reference state" );
                    result->addNext( error );
                    return result;

                }

            }

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

        _DP_BDhatD = -_DP_BDhatQ * _N.topLeftCorner( nFreeMicroDOF, nFreeMacroDOF );

        _DP_BQhatQ = _N.bottomRightCorner( nGhostMicroDOF, nGhostMacroDOF ) * _DP_BDhatQ;

        _DP_BQhatD = _N.bottomLeftCorner( nGhostMicroDOF, nFreeMacroDOF )
                   + _N.bottomRightCorner( nGhostMicroDOF, nGhostMacroDOF ) * _DP_BDhatD;

        return NULL;
    }

    errorOut overlapCoupling::processDomainReference( const unsigned int &microIncrement,
                                                      const std::string &domainName,
                                                      const unsigned int cellID, const uIntVector &macroNodes,
                                                      domainFloatMap   &referenceMicroDomainMass,
                                                      domainFloatVectorMap &referenceMicroDomainCentersOfMass,
                                                      domainFloatVectorMap &referenceMicroDomainMomentsOfInertia,
                                                      std::unordered_map< uIntType, floatVector > &domainReferenceXiVectors,
                                                      floatVector &domainCenterOfMassShapeFunctionValues,
                                                      std::unordered_map< uIntType, floatVector > &domainMicroPositionShapeFunctionValues ){
        /*!
         * Process the domain for use with preparing the reference configuration
         *
         * :param const unsigned int &microIncrement: The micro-scale increment at which to compute the reference
         * :param const std::string &domainName: The name of the domain
         * :param const unsigned int cellID: The global cell ID number
         * :param const uIntVector &macroNodes: The nodes of the macro domain
         * :param domainFloatMap   &referenceMicroDomainMass: The reference mass of the micro domain. This
         *     should be a reference to the micro-domain mass vector
         * :param domainFloatVectorMap &referenceMicroDomainCentersOfMass: The reference micro-domain center of mass
         *     vector for all of the micro domains
         * :param domainFloatVectorMap &referenceDomainMomentsOfInertia: The reference moment of inertia for all of the micro domains
         * :param std::unordered_map< uIntType, floatVector > &domainReferenceXiVectors: The reference Xi vectors for the domain.
         * :param floatVector &domainCenterOfMassShapeFunctionValues: The shape function values at the 
         *     center of mass
         * :param std::unordered_map< uIntType, floatVector > &domainMicroPositionShapeFunctionValues: The shape function values at the
         *     micro node positions inside of the domain. This is only computed if indicated by the 
         *     configuration file.
         */

        //Clear the existing Xi and shapefunction maps
        domainReferenceXiVectors.clear( );
        domainMicroPositionShapeFunctionValues.clear( );

        //Process the domain mass data
        floatVector domainCenterOfMass;

        errorOut error = processDomainMassData( microIncrement, domainName, referenceMicroDomainMass,
                                                referenceMicroDomainCentersOfMass, referenceMicroDomainMomentsOfInertia,
                                                domainReferenceXiVectors );
#ifdef TESTACCESS
        _test_domainMass[ cellID ].emplace( domainName, referenceMicroDomainMass[ domainName ] );
        _test_domainCOM[ cellID ].emplace( domainName, referenceMicroDomainCentersOfMass[ domainName ] );
        _test_domainXi[ cellID ].emplace( domainName, domainReferenceXiVectors );
#endif

        if ( error ){
            
            errorOut result = new errorNode( "processDomainReference",
                                             "Error in processing the mass data for the micro domain '" + domainName + "'" );
            result->addNext( error );
            return result;

        }

        //Compute the domain's shape function information
        error = computeDomainShapeFunctionInformation( cellID, domainName, microIncrement,
                                                       referenceMicroDomainCentersOfMass[ domainName ],
                                                       domainCenterOfMassShapeFunctionValues,
                                                       domainMicroPositionShapeFunctionValues );

        _referenceCellDomainCenterOfMassShapefunctions[ cellID ].emplace( domainName, domainCenterOfMassShapeFunctionValues );

#ifdef TESTACCESS
        if ( _inputProcessor.getCouplingInitialization( )[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){
            _test_domainMUP[ cellID ].emplace( domainName, domainMicroPositionShapeFunctionValues );
        }
#endif

        if ( error ){
            errorOut result = new errorNode( "processDomainReference",
                                             "Error in computing the shape function values for the domain center of mass or the domainMicroPositionShapeFunctionValues" );
            result->addNext( error );
            return result;
        }

        //Get the domain node ids
        uIntVector domainNodes;
        error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, domainName, domainNodes );

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
        if ( _inputProcessor.getCouplingInitialization( )[ "projection_type" ].as< std::string >( ).compare( "averaged_l2_projection" ) == 0 ){

            SparseMatrix domainE;
            error = DOFProjection::assembleMicroDomainHomogenizationMatrixContribution( domainName, domainNodes,
                                                                                        *_inputProcessor.getMicroDensities( ),
                                                                                        *_inputProcessor.getMicroVolumes( ),
                                                                                        *_inputProcessor.getMicroWeights( ),
                                                                                        domainReferenceXiVectors,
                                                                                        *_inputProcessor.getMicroGlobalToLocalDOFMap( ),
                                                                                        referenceMicroDomainMass, referenceMicroDomainMomentsOfInertia,
                                                                                        *_inputProcessor.getMicroDomainIDMap( ),
                                                                                        domainE );

            if ( _homogenizationMatrix_initialized ){

                _homogenizationMatrix += domainE;

            }
            else{

                _homogenizationMatrix = domainE;
                _homogenizationMatrix_initialized = true;

            }

        }


        return NULL;

    }

    errorOut overlapCoupling::processDomainMassData( const unsigned int &microIncrement, const std::string &domainName,
                                                     domainFloatMap &domainMass, domainFloatVectorMap &domainCenterOfMass,
                                                     domainFloatVectorMap &domainMomentOfInertia,
                                                     std::unordered_map< uIntType, floatVector > &domainXiVectors ){
        /*!
         * Process a micro-scale domain
         *
         * :param const unsigned int microIncrement: The micro increment to process
         * :param const std::string &domainName: The name of the domain
         * :param domainFloatMap &domainMass: The mass of the domain
         * :param domainFloatVectorMap &domainCenterOfMass: The center of mass of the domain
         * :param domainFloatVectorMap &domainMomentOfInertia: The moment of inertia of the domain
         * :param std::unordered_map< uIntType, floatVector > &domainXiVectors: The Xi vectors of the domain
         */

        //Get the domain's nodes
        uIntVector domainNodes;

        errorOut error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, domainName, domainNodes );

        if ( error ){

            errorOut result = new errorNode( "processDomain",
                                             "Error in getting the nodes from the micro domain '" + domainName + "'" );
            result->addNext( error );
            return result;

        }

        //Compute the center of mass of the domain
        floatType   mass;
        floatVector centerOfMass;
        error = DOFProjection::computeDomainCenterOfMass( _dim, domainNodes, *_inputProcessor.getMicroVolumes( ),
                                                          *_inputProcessor.getMicroDensities( ),
                                                          *_inputProcessor.getMicroNodeReferencePositions( ),
                                                          *_inputProcessor.getMicroDisplacements( ),
                                                          *_inputProcessor.getMicroWeights( ),
                                                          mass, centerOfMass );

        if ( error ){

            errorOut result = new errorNode( "processDomain", "Error in calculation of '" + domainName + "' center of mass" );
            result->addNext( error );
            return result;

        }

        domainMass.emplace( domainName, mass );
        domainCenterOfMass.emplace( domainName, centerOfMass );

        //Compute the relative position vectors
        floatVector momentOfInertia;
        error = DOFProjection::computeDomainXis( _dim, domainNodes,
                                                 *_inputProcessor.getMicroNodeReferencePositions( ),
                                                 *_inputProcessor.getMicroDisplacements( ),
                                                 *_inputProcessor.getMicroVolumes( ),
                                                 *_inputProcessor.getMicroDensities( ),
                                                 *_inputProcessor.getMicroWeights( ),
                                                 domainCenterOfMass[ domainName ], domainXiVectors, momentOfInertia );

        if ( error ){
            
            errorOut result = new errorNode( "processDomain", "Error in calculation of '" + domainName + "' xi vectors" );
            result->addNext( error );
            return result;

        }

        domainMomentOfInertia.emplace( domainName, momentOfInertia );

        return NULL;

    }

    errorOut overlapCoupling::computeDomainShapeFunctionInformation( const unsigned int &cellID,
                                                                     const std::string &domainName,
                                                                     const unsigned int &microIncrement,
                                                                     const floatVector &domainCenterOfMass,
                                                                     floatVector &domainCenterOfMassShapeFunctionValues,
                                                                     std::unordered_map< uIntType, floatVector > &domainMicroPositionShapeFunctionValues ){
        /*!
         * Compute the shape function values at the required locations
         *
         * :param const unsigned int &cellID: The ID number of the cell which contains the domain
         * :param const std::string &domainName: The name of the nodeset which defines the micro-scale domain
         * :param const unsigned int &microIncrement: The micro-scale increment to analyze
         * :param const floatVector &domainCenterOfMass: The center of mass of the domain
         * :param floatVector &domainCenterOfMassShapeFunctionValues: The shapefunction values at the center of mass
         * :param std::unordered_map< uIntType, floatVector > &domainMicroPositionShapeFunctionValues:
         *     The shapefunction values at all of the micro nodes contained within the domain.
         */

        //Compute the shape functions of the domain's center of mass
        errorOut error = computeShapeFunctionsAtPoint( cellID,
                                                       *_inputProcessor.getMacroNodeReferencePositions( ),
                                                       *_inputProcessor.getMacroDisplacements( ),
                                                       *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                       domainCenterOfMass, domainCenterOfMassShapeFunctionValues );

        if ( error ){

            errorOut result = new errorNode( "computeDomainShapeFunctionInformation",
                                             "Error in the computation of the shape function at the center of mass for a micro domain" );
            result->addNext( error );
            return result;

        }

        //Get the domain's nodes
        uIntVector domainNodes;

        error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, domainName, domainNodes );

        if ( error ){

            errorOut result = new errorNode( "computeDomainShapeFunctionInformation",
                                             "Error in the extraction of the nodes in the micro domain" );
            result->addNext( error );
            return result;

        }

        if ( _inputProcessor.computeMicroShapeFunctions( ) ){

            //Get the micro-node positions
            std::unordered_map< uIntType, floatVector > microNodePositions;
    
            const std::unordered_map< uIntType, floatVector > *microReferencePositions
                = _inputProcessor.getMicroNodeReferencePositions( );
            const std::unordered_map< uIntType, floatVector > *microDisplacements
                = _inputProcessor.getMicroDisplacements( );

            microNodePositions.reserve( microReferencePositions->size( ) );
            for ( auto it = domainNodes.begin( ); it != domainNodes.end( ); it++ ){
  
                auto microReferencePosition = microReferencePositions->find( *it );

                if ( microReferencePosition == microReferencePositions->end( ) ){

                    return new errorNode( "computeDomainShapeFunctionInformation",
                                          "Micro node " + std::to_string( *it ) + " was not found in the reference position map" );

                }

                auto microDisplacement = microDisplacements->find( *it );
                if ( microDisplacement == microDisplacements->end( ) ){

                    return new errorNode( "computeDomainShapeFunctionInformation",
                                          "Micro node " + std::to_string( *it ) + " was not found in the displacement map" );

                }

                microNodePositions.emplace( *it, microReferencePosition->second + microDisplacement->second );
    
            }
    
            //Compute the shape function values at the micro positions
            error = computeShapeFunctionsAtPoints( cellID,
                                                   *_inputProcessor.getMacroNodeReferencePositions( ),
                                                   *_inputProcessor.getMacroDisplacements( ),
                                                   *_inputProcessor.getMacroNodeReferenceConnectivity( ),
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
                                                             domainFloatMap &freeDomainMass, domainFloatMap &ghostDomainMass,
                                                             domainFloatVectorMap &freeDomainCM,
                                                             domainFloatVectorMap &ghostDomainCM ){
        /*!
         * Compute the centers of mass for an micro increment. Also computes the mass of the micro-scale domains
         *
         * :param const unsigned int microIncrement: The micro increment at which to compute the centers of mass
         * :param const unsigned int macroIncrement: The macro increment at which to compute the centers of mass
         * :param domainFloatMap &freeDomainMass: The mass of the free domains
         * :param domainFloatMap &ghostDomainMass: The mass of the ghost domains
         * :param domainFloatVectorMap &freeDomainCM: The center of mass of the free domains
         * :param domainFloatVectorMap &ghostDomainCM: The center of mass of the ghost domains
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

        freeDomainMass.clear( );
        freeDomainCM.clear( );

        unsigned int indx = 0;

        for ( auto name = freeDomains->begin( ); name != freeDomains->end( ); name++ ){

            error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, *name, domainNodes );

            if ( error ){

                errorOut result = new errorNode( "computeIncrementCentersOfMass", "Error in extraction of the free domain's nodes" );
                result->addNext( error );
                return result;

            }

            floatType mass;
            floatVector centerOfMass;
            error = DOFProjection::computeDomainCenterOfMass( _dim, domainNodes, *_inputProcessor.getMicroVolumes( ),
                                                              *_inputProcessor.getMicroDensities( ),
                                                              *_inputProcessor.getMicroNodeReferencePositions( ),
                                                              *_inputProcessor.getMicroDisplacements( ),
                                                              *_inputProcessor.getMicroWeights( ),
                                                              mass, centerOfMass );
            if ( error ){

                errorOut result = new errorNode( "computeIncrementCentersOfMass", "Error in calculation of '" + *name + "' center of mass" );
                result->addNext( error );
                return result;

            }
            
            freeDomainMass.emplace( *name, mass );
            freeDomainCM.emplace( *name, centerOfMass );

            indx++;

        }

        //Loop over the ghost micro-domains
        const stringVector* ghostDomains = _inputProcessor.getGhostMicroDomainNames( );

        ghostDomainMass.clear( );
        ghostDomainCM.clear( );

        indx = 0;

        for ( auto name = ghostDomains->begin( ); name != ghostDomains->end( ); name++ ){

            error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, *name, domainNodes );

            if ( error ){

                errorOut result = new errorNode( "computeIncrementCentersOfMass", "Error in extraction of the ghost domain's nodes" );
                result->addNext( error );
                return result;

            }

            floatType mass;
            floatVector centerOfMass;
            error = DOFProjection::computeDomainCenterOfMass( _dim, domainNodes, *_inputProcessor.getMicroVolumes( ),
                                                              *_inputProcessor.getMicroDensities( ),
                                                              *_inputProcessor.getMicroNodeReferencePositions( ),
                                                              *_inputProcessor.getMicroDisplacements( ),
                                                              *_inputProcessor.getMicroWeights( ),
                                                              mass, centerOfMass );

            if ( error ){

                errorOut result = new errorNode( "computeIncrementCentersOfMass", "Error in calculation of '" + *name + "' center of mass" );
                result->addNext( error );
                return result;

            }

            ghostDomainMass.emplace( *name, mass );
            ghostDomainCM.emplace( *name, centerOfMass );

            indx++;

        }

        return NULL;
    }

    errorOut overlapCoupling::buildMacroDomainElement( const unsigned int cellID,
                                                       const std::unordered_map< uIntType, floatVector > &nodeLocations,
                                                       const std::unordered_map< uIntType, uIntVector > &connectivity,
                                                       std::unique_ptr< elib::Element > &element ){
        /*!
         * Construct a finite element representation of the macro domain
         *
         * :param const unsigned int cellID: The macro cell ID number
         * :param const std::unordered_map< uIntType, floatVector > &nodeLocations: The nodal location vector
         * :param const std::unordered_map< uIntType, uIntVector > &connectivity: The connectivity vector
         * :param std::unique_ptr< elib::Element > &element: The element representation of the macro domain
         */

        //Get the XDMF cell type
        auto connectivityCellIndices = connectivity.find( cellID );
        if ( connectivityCellIndices == connectivity.end( ) ){

            return new errorNode( "buildMacroDomainElement",
                                  "Cell ID " + std::to_string( cellID ) + " was not found in the connectivity map" );

        }
        unsigned int cellType = connectivityCellIndices->second[ 0 ];

        //Get the element name
        auto it = elib::XDMFTypeToElementName.find( cellType );

        if ( it == elib::XDMFTypeToElementName.end( ) ){

            return new errorNode( "buildMacroDomainElement",
                                  "The cell type " + std::to_string(cellType) + " is not supported" );

        }
        
        //Get the nodes from the file
        elib::vecOfvec nodes( connectivityCellIndices->second.size( ) - 1, elib::vec( _dim, 0 ) );
        uIntVector globalNodeIds( connectivityCellIndices->second.begin( ) + 1,
                                  connectivityCellIndices->second.end( ) );
        for ( auto nodeId = globalNodeIds.begin( ); nodeId != globalNodeIds.end( ); nodeId++ ){

            auto nodeLocation = nodeLocations.find( *nodeId );

            if ( nodeLocation == nodeLocations.end( ) ){

                return new errorNode( "Node " + std::to_string( *nodeId ) + " was not found in the node locations map" );

            }

            uIntType index = nodeId - globalNodeIds.begin( );

            nodes[ index ] = nodeLocation->second;

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
                                                       const std::unordered_map< uIntType, floatVector > &nodeReferenceLocations,
                                                       const std::unordered_map< uIntType, floatVector > &nodeDisplacements,
                                                       const std::unordered_map< uIntType, uIntVector > &connectivity,
                                                       std::unique_ptr< elib::Element > &element ){
        /*!
         * Construct a finite element representation of the macro domain
         *
         * :param const unsigned int cellID: The macro cell ID number
         * :param const std::unordered_map< uIntType, floatVector > &nodeReferenceLocations: The nodal reference location vector
         * :param const std::unordered_map< uIntType, floatVector > &nodeDisplacements: The nodal displacement vector
         * :param const std::unordered_map< uIntType, uIntVector > &connectivity: The connectivity vector
         * :param std::unique_ptr< elib::Element > &element: The element representation of the macro domain
         */

        //Make sure the cellID is allowable
        auto connectivityCellIndices = connectivity.find( cellID );

        if ( connectivityCellIndices == connectivity.end( ) ){

            return new errorNode( "buildMacroDomainElement", 
                                  "Cell " + std::to_string( cellID ) + " was not found in the connectivity map" );

        }

        //Get the XDMF cell type
        unsigned int cellType = connectivityCellIndices->second[ 0 ];

        //Get the element name
        auto it = elib::XDMFTypeToElementName.find( cellType );

        if ( it == elib::XDMFTypeToElementName.end( ) ){

            return new errorNode( "buildMacroDomainElement",
                                  "The cell type " + std::to_string(cellType) + " is not supported" );

        }

        //Get the nodes from the file
        elib::vecOfvec referenceNodes( connectivityCellIndices->second.size( ) - 1, elib::vec( _dim, 0 ) );
        elib::vecOfvec displacements( connectivityCellIndices->second.size( ) - 1, elib::vec( _dim, 0 ) );
        uIntVector globalNodeIds( connectivityCellIndices->second.begin( ) + 1,
                                  connectivityCellIndices->second.end( ) );
        for ( auto nodeId = globalNodeIds.begin( ); nodeId != globalNodeIds.end( ); nodeId++ ){

            uIntType index = nodeId - globalNodeIds.begin( );

            auto nodeReferenceLocation = nodeReferenceLocations.find( *nodeId );

            if ( nodeReferenceLocation == nodeReferenceLocations.end( ) ){

                return new errorNode( "buildMacroDomainElement",
                                      "The node " + std::to_string( *nodeId ) + " was not found in the node reference location map" );

            }

            auto nodeDisplacement = nodeDisplacements.find( *nodeId );

            if ( nodeDisplacement == nodeDisplacements.end( ) ){

                return new errorNode( "buildMacroDomainElement",
                                      "The node " + std::to_string( *nodeId ) + " was not found in the node displacement map" );

            }

            referenceNodes[ index ] = nodeReferenceLocation->second;
            displacements[ index ]  = nodeDisplacement->second;

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

    errorOut overlapCoupling::computeShapeFunctionsAtPoint( const unsigned int cellID,
                                                            const std::unordered_map< uIntType, floatVector > &nodeLocations,
                                                            const std::unordered_map< uIntType, uIntVector > &connectivity,
                                                            const floatVector &point,
                                                            floatVector &shapeFunctions ){
        /*!
         * Compute the shape functions of a given macro-scale domain at a single point
         *
         * :param const unsigned int &cellID: The cell ID at which to compute the shape functions
         * :param const std::unordered_map< uIntType, floatVector > &nodeLocations: The nodal location map
         * :param const std::unordered_map< uIntType, uIntVector > &connectivity: The connectivity map
         * :param const floatVector &point: The point at which to compute the shape functions
         * :param floatVector &shapeFunctions: The shapefunctions at the point
         */

        if ( point.size( ) != _dim ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "This function only works for a single point of dimension " + std::to_string( _dim ) );

        }

        std::unordered_map< uIntType, floatVector > pointMap;
        pointMap.emplace( 0, point );

        std::unordered_map< uIntType, floatVector > shapefunctionMap;

        errorOut error = computeShapeFunctionsAtPoints( cellID, nodeLocations, connectivity, pointMap, shapefunctionMap );

        if ( error ){

            errorOut result = new errorNode( "computeShapeFunctionsAtPoints",
                                             "Error when computing shape functions" );
            result->addNext( error );
            return result;

        }

        shapeFunctions = shapefunctionMap[ 0 ];

        return NULL;

    }

    errorOut overlapCoupling::computeShapeFunctionsAtPoint( const unsigned int cellID,
                                                            const std::unordered_map< uIntType, floatVector > &nodeReferenceLocations,
                                                            const std::unordered_map< uIntType, floatVector > &nodeDisplacements,
                                                            const std::unordered_map< uIntType, uIntVector > &connectivity,
                                                            const floatVector &point,
                                                            floatVector &shapeFunctions ){
        /*!
         * Compute the shape functions of a given macro-scale domain at a single point
         *
         * :param const unsigned int &cellID: The cell ID at which to compute the shape functions
         * :param const std::unordered_map< uIntType, floatVector > &nodeLocations: The nodal location map
         * :param const std::unordered_map< uIntType, floatVector > &nodeLocations: The nodal displacement map
         * :param const std::unordered_map< uIntType, uIntVector > &connectivity: The connectivity map
         * :param const floatVector &point: The point at which to compute the shape functions
         * :param floatVector &shapeFunctions: The shapefunctions at the point
         */

        if ( point.size( ) != _dim ){

            return new errorNode( "computeShapeFunctionsAtPoints",
                                  "This function only works for a single point of dimension " + std::to_string( _dim ) );

        }

        std::unordered_map< uIntType, floatVector > pointMap;
        pointMap.emplace( 0, point );

        std::unordered_map< uIntType, floatVector > shapefunctionMap;

        errorOut error = computeShapeFunctionsAtPoints( cellID, nodeReferenceLocations, nodeDisplacements, 
                                                        connectivity, pointMap, shapefunctionMap );

        if ( error ){

            errorOut result = new errorNode( "computeShapeFunctionsAtPoints",
                                             "Error when computing shape functions" );
            result->addNext( error );
            return result;

        }

        shapeFunctions = shapefunctionMap[ 0 ];

        return NULL;

    }

    errorOut overlapCoupling::computeShapeFunctionsAtPoints( const unsigned int cellID,
                                                             const std::unordered_map< uIntType, floatVector > &nodeLocations,
                                                             const std::unordered_map< uIntType, uIntVector > &connectivity,
                                                             const std::unordered_map< uIntType, floatVector > &points,
                                                             std::unordered_map< uIntType, floatVector > &shapeFunctions ){
        /*!
         * Compute the shape functions of a given macro-scale domain at the given points
         *
         * :param const unsigned int &cellID: The cell ID at which to compute the shape functions
         * :param const std::unordered_map< uIntType, floatVector > &nodeLocations: The nodal location map
         * :param const std::unordered_map< uIntType, uIntVector > &connectivity: The connectivity map
         * :param const std::unordered_map< uIntType, floatVector > &points: The points at which to compute the shape functions
         * :param std::unordered_map< uIntType, floatVector > &shapeFunctions: The shapefunctions at the points
         */

        //Build the element representing the macro-scale domain
        std::unique_ptr< elib::Element > element;
        errorOut error = overlapCoupling::buildMacroDomainElement( cellID, nodeLocations,
                                                                   connectivity, element );

        //Make sure the number of points and the size of the output are consistent
        unsigned nPoints = points.size( );

        shapeFunctions.clear();
        shapeFunctions.reserve( nPoints );

        //Compute the shape functions at each point
        floatVector point;
        floatVector localPosition, pointShapeFunctions;

        //Loop over the output vector
        for ( auto p = points.begin( ); p != points.end( ); p++ ){

            error = element->compute_local_coordinates( p->second, localPosition );

            if ( !element->local_point_inside( localPosition ) ){

                pointShapeFunctions = { 0 };
                shapeFunctions.emplace( p->first, pointShapeFunctions );
                continue;

            }

            if ( error ) {

                return new errorNode( "computeShapeFunctionsAtPoints",
                                      "Error in computing the local coordinates for point " + std::to_string( p->first ) );

            }

            error = element->get_shape_functions( localPosition, pointShapeFunctions );

            if ( error ) {

                return new errorNode( "computeShapeFunctionsAtPoints",
                                      "Error in the computation of the shape functions for point " + std::to_string( p->first ) );

            }

            shapeFunctions.emplace( p->first, pointShapeFunctions );

        }

        return NULL;
    }

    errorOut overlapCoupling::computeShapeFunctionsAtPoints( const unsigned int cellID,
                                                             const std::unordered_map< uIntType, floatVector > &nodeReferenceLocations,
                                                             const std::unordered_map< uIntType, floatVector > &nodeDisplacements,
                                                             const std::unordered_map< uIntType, uIntVector > &connectivity,
                                                             const std::unordered_map< uIntType, floatVector > &points,
                                                             std::unordered_map< uIntType, floatVector > &shapeFunctions ){
        /*!
         * Compute the shape functions of a given macro-scale domain at the given points
         *
         * :param const unsigned int &cellID: The cell ID at which to compute the shape functions
         * :param const std::unordered_map< uIntType, floatVector > &nodeReferenceLocations: The nodal reference location map
         * :param const std::unordered_map< uIntType, floatVector > &nodeDisplacements: The nodal displacement map
         * :param const std::unordered_map< uIntType, uIntVector > &connectivity: The connectivity map
         * :param const std::unordered_map< uIntType, floatVector > &points: The points at which to compute the shape functions
         * :param std::unordered_map< uIntType, floatVector > &shapeFunctions: The shapefunctions at the points
         */

        //Build the element representing the macro-scale domain
        std::unique_ptr< elib::Element > element;
        errorOut error = overlapCoupling::buildMacroDomainElement( cellID, nodeReferenceLocations, nodeDisplacements,
                                                                   connectivity, element );

        //Make sure the number of points and the size of the output are consistent
        unsigned nPoints = points.size( );

        shapeFunctions.clear();
        shapeFunctions.reserve( nPoints );

        //Compute the shape functions at each point
        floatVector point;
        floatVector localPosition, pointShapeFunctions;

        //Loop over the output vector
        for ( auto p = points.begin( ); p != points.end( ); p++ ){

            error = element->compute_local_coordinates( p->second, localPosition );

            if ( !element->local_point_inside( localPosition ) ){

                pointShapeFunctions = { 0 };
                shapeFunctions.emplace( p->first, pointShapeFunctions );
                continue;

            }

            if ( error ) {

                return new errorNode( "computeShapeFunctionsAtPoints",
                                      "Error in computing the local coordinates for point " + std::to_string( p->first ) );

            }

            error = element->get_shape_functions( localPosition, pointShapeFunctions );

            if ( error ) {

                return new errorNode( "computeShapeFunctionsAtPoints",
                                      "Error in the computation of the shape functions for point " + std::to_string( p->first ) );

            }

            shapeFunctions.emplace( p->first, pointShapeFunctions );

        }

        return NULL;
            
    }

    errorOut overlapCoupling::computeShapeFunctionGradientsAtPoints( const unsigned int cellID,
                                                       const std::unordered_map< uIntType, floatVector > &nodeReferenceLocations,
                                                       const std::unordered_map< uIntType, floatVector > &nodeDisplacements,
                                                       const std::unordered_map< uIntType, uIntVector > &connectivity,
                                                       const std::unordered_map< uIntType, floatVector > points,
                                                       std::unordered_map< uIntType, floatVector > &shapeFunctionGradients ){
        /*!
         * Compute the shape functions of a given macro-scale domain at the given points
         *
         * :param const unsigned int &cellID: The cell ID at which to compute the shape functions
         * :param const std::unordered_map< uIntType, floatVector > &nodeReferenceLocations: The nodal reference location vector
         * :param const std::unordered_map< uIntType, floatVector > &nodeDisplacements: The nodal reference location vector
         * :param const std::unordered_map< uIntType, uIntVector > &connectivity: The connectivity vector
         * :param const std::unordered_map< uIntType, floatVector > &points: The points at which to compute the shape functions
         * :param std::unordered_map< uIntType, floatVector > &shapeFunctionGradients: The gradients of the shapefunctions
         *     at the points. The output vector is organized [ N11_x, N11_y, N11_z, N12_x, ... ] where the first index is the
         *     point number, the second index is the node number, and the _ indicates a gradient w.r.t. the
         *     third index.
         */

        //Build the element representing the macro-scale domain
        std::unique_ptr< elib::Element > element;
        errorOut error = overlapCoupling::buildMacroDomainElement( cellID, nodeReferenceLocations, nodeDisplacements,
                                                                   connectivity, element );

        //Make sure the number of points and the size of the output are consistent
        unsigned nPoints = points.size( ) / _dim;
        if ( ( points.size( ) % _dim ) > 0 ){

            return new errorNode( "computeShapeFunctionGradientsAtPoints",
                                  "The points vector is inconsistent with the dimension\n"
                                  "    points.size( ): " + std::to_string( points.size( ) ) + "\n" +
                                  "    nPoints: " + std::to_string( nPoints ) );

        }

        shapeFunctionGradients.clear();
        shapeFunctionGradients.reserve( element->reference_nodes.size( ) );

        //Compute the shape functions at each point
        floatVector point;
        floatMatrix dNdx;
        floatVector localPosition, pointShapeFunctionGradientsVec;

        //Loop over the output vector
        for ( auto p = points.begin( ); p != points.end( ); p++ ){

            error = element->compute_local_coordinates( p->second, localPosition );

            if ( !element->local_point_inside( localPosition ) ){

                pointShapeFunctionGradientsVec = { 0 };
                shapeFunctionGradients.emplace( p->first, pointShapeFunctionGradientsVec );
                continue;

            }

            if ( error ) {

                return new errorNode( "computeShapeFunctionGradientsAtPoints",
                                      "Error in computing the local coordinates for point " + std::to_string( p->first ) );

            }

            error = element->get_global_shapefunction_gradients( localPosition, dNdx );

            if ( error ) {

                return new errorNode( "computeShapeFunctionGradientsAtPoints",
                                      "Error in the computation of the shape functions for point " + std::to_string( p->first ) );

            }

            pointShapeFunctionGradientsVec = vectorTools::appendVectors( dNdx );
            shapeFunctionGradients.emplace( p->first, pointShapeFunctionGradientsVec );

        }

        return NULL;
            
    }

    errorOut overlapCoupling::computeShapeFunctionsAtReferenceCentersOfMass( ){
        /*!
         * Compute the shape functions at the reference centers of mass
         */

        //Loop over the free domains
        stringVector domainNames;
        std::unordered_map< uIntType, floatVector > domainCOMs;

        std::unordered_map< uIntType, floatVector > macroDomainShapeFunctions;
        errorOut error = NULL;

        const uIntVector *freeMacroCellIds = _inputProcessor.getFreeMacroCellIds( );
        const std::unordered_map< uIntType, stringVector > *macroCellToMicroDomainMap = _inputProcessor.getMacroCellToDomainMap( );

        _referenceGhostMicroDomainCenterOfMassShapeFunctions.clear( );

        for ( auto cellID = freeMacroCellIds->begin( ); cellID != freeMacroCellIds->end( ); cellID++ ){

            //Get the centers of mass of the macro-domain
            domainCOMs.clear( );
            auto cellCentersOfMass = _referenceGhostMicroDomainCentersOfMass.find( *cellID );

            if ( cellCentersOfMass == _referenceGhostMicroDomainCentersOfMass.end( ) ){

                return new errorNode( "computeShapeFunctionsAtReferenceCentersOfMass",
                                      "The macro cell " + std::to_string( *cellID ) +
                                      " was not found in the reference ghost micro domain centers of mass" );

            }

            stringVector domainNames( cellCentersOfMass->second.size( ) );
            domainCOMs.clear( );

            uIntType dindex = 0;
            for ( auto domain = cellCentersOfMass->second.begin( ); domain != cellCentersOfMass->second.end( ); domain++, dindex++ ){

                domainNames[ dindex ] = domain->first;
                domainCOMs.emplace( dindex, domain->second );

            }

            error = computeShapeFunctionsAtPoints( *cellID,
                                                   *_inputProcessor.getMacroNodeReferencePositions( ),
                                                   *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                   domainCOMs, macroDomainShapeFunctions );

            if ( error ){

                errorOut result = new errorNode( "computeShapeFunctionsAtReferenceCentersOfMass",
                                                 "Error in computation of the shape functions at the reference ghost micro centers of mass" );
                result->addNext( error );
                return result;

            }

            //Append the shape functions to the storage vector
            domainFloatVectorMap temp;
            temp.reserve( macroDomainShapeFunctions.size( ) );
            for ( auto domainName = domainNames.begin( ); domainName != domainNames.end( ); domainName++ ){

                uIntType index = domainName - domainNames.begin( );

                temp.emplace( *domainName, macroDomainShapeFunctions[ index ] ); 

            }

            _referenceGhostMicroDomainCenterOfMassShapeFunctions.emplace( *cellID, temp );

        }

        //Loop over the ghost domains

        const uIntVector *ghostMacroCellIds = _inputProcessor.getGhostMacroCellIds( );

        _referenceFreeMicroDomainCenterOfMassShapeFunctions.clear( );

        for ( auto cellID = ghostMacroCellIds->begin( ); cellID != ghostMacroCellIds->end( ); cellID++ ){

            //Get the centers of mass of the macro-domain
            auto cellCentersOfMass = _referenceFreeMicroDomainCentersOfMass.find( *cellID );

            if ( cellCentersOfMass == _referenceFreeMicroDomainCentersOfMass.end( ) ){

                return new errorNode( "computeShapeFunctionsAtReferenceCentersOfMass",
                                      "The macro cell " + std::to_string( *cellID ) +
                                      " was not found in the reference free micro domain centers of mass" );

            }

            stringVector domainNames( cellCentersOfMass->second.size( ) );
            domainCOMs.clear( );

            uIntType dindex = 0;
            for ( auto domain = cellCentersOfMass->second.begin( ); domain != cellCentersOfMass->second.end( ); domain++, dindex++ ){

                domainNames[ dindex ] = domain->first;
                domainCOMs.emplace( dindex, domain->second );

            }

            error = computeShapeFunctionsAtPoints( *cellID,
                                                   *_inputProcessor.getMacroNodeReferencePositions( ),
                                                   *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                   domainCOMs, macroDomainShapeFunctions );

            if ( error ){

                errorOut result = new errorNode( "computeShapeFunctionsAtReferenceCentersOfMass",
                                                 "Error in computation of the shape functions at the reference free micro centers of mass" );
                result->addNext( error );
                return result;

            }

            //Append the shape functions to the storage vector
            domainFloatVectorMap temp;
            temp.reserve( macroDomainShapeFunctions.size( ) );
            for ( auto domainName = domainNames.begin( ); domainName != domainNames.end( ); domainName++ ){

                uIntType index = domainName - domainNames.begin( );

                temp.emplace( *domainName, macroDomainShapeFunctions[ index ] ); 

            }

            _referenceFreeMicroDomainCenterOfMassShapeFunctions.emplace( *cellID, temp );

        }

        return NULL;
        
    }

    errorOut overlapCoupling::addDomainContributionToInterpolationMatrix( const uIntVector  &domainNodes,
                                                                          const uIntVector  &macroNodes,
                                                                          const std::unordered_map< uIntType, floatVector > &domainReferenceXis,
                                                                          const floatVector &domainCenterOfMassShapeFunctionValues ){
        /*!
         * Add the contribution of a domain to the interpolation matrices
         *
         * :param const std::string &domainName: The name of the domain
         * :param const std::unordered_map< uIntType, floatVector > &domainReferenceXis: The micro-position vectors of 
         *     the nodes within the domain
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

    errorOut overlapCoupling::projectDegreesOfFreedom( const bool useUpdatedFreeDOF ){
        /*!
         * Project the degrees of freedom of the ghost nodes
         * for the current increment
         *
         * :const bool useUpdatedFreeDOF: A flag indicating if the free displacements vectors to be
         *     used should be the updated values or the values from the input file
         */

        //Get the displacement vectors
        const std::unordered_map< uIntType, floatVector > *macroDispDOFVector = _inputProcessor.getMacroDispDOFVector( );
        const std::unordered_map< uIntType, floatVector > *microDisplacements = _inputProcessor.getMicroDisplacements( );

        //Get the free and ghost node ids
        const uIntVector *freeMacroNodeIds = _inputProcessor.getFreeMacroNodeIds( );
        const uIntVector *ghostMacroNodeIds = _inputProcessor.getGhostMacroNodeIds( );

        const uIntVector *freeMicroNodeIds = _inputProcessor.getFreeMicroNodeIds( );
        const uIntVector *ghostMicroNodeIds = _inputProcessor.getGhostMicroNodeIds( );

        //Set the number of displacement degrees of freedom
        unsigned int nMacroDispDOF = _dim + _dim * _dim;
        unsigned int nMicroDispDOF = _dim;

        //Assemble the free displacements
        floatVector *freeMacroDisplacements;
        floatVector *freeMicroDisplacements;

        floatVector store1, store2;

        if ( useUpdatedFreeDOF ){

            freeMicroDisplacements = &_updatedFreeMicroDispDOFValues;
            freeMacroDisplacements = &_updatedFreeMacroDispDOFValues;

        }
        else{

            store1.resize( nMacroDispDOF * freeMacroNodeIds->size( ) );
            store2.resize( nMicroDispDOF * freeMicroNodeIds->size( ) );

            freeMacroDisplacements = &store1;
            freeMicroDisplacements = &store2;

            const DOFMap *macroGlobalToLocalDOFMap = _inputProcessor.getMacroGlobalToLocalDOFMap( );
            const DOFMap *microGlobalToLocalDOFMap = _inputProcessor.getMicroGlobalToLocalDOFMap( );
   
            for ( auto it = freeMacroNodeIds->begin( ); it != freeMacroNodeIds->end( ); it++ ){

                auto macroDispDOF = macroDispDOFVector->find( *it );

                if ( macroDispDOF == macroDispDOFVector->end( ) ){

                    return new errorNode( "projectDegreesOfFreedom",
                                          "Macro node " + std::to_string( *it ) +
                                          " was not found in the macro displacement dof vector map" );

                }

                if ( macroDispDOF->second.size( ) != nMacroDispDOF ){

                    return new errorNode( "projectDegreesOfFreedom",
                                          "Macro node " + std::to_string( *it ) +
                                          " does not have a dimensionally consistent number of degrees of freedom" );

                }

                auto map = macroGlobalToLocalDOFMap->find( *it );

                if ( map == macroGlobalToLocalDOFMap->end( ) ){

                    return new errorNode( "projectDegreesOfFreedom",
                                          "Macro node " + std::to_string( *it ) +
                                          " was not found in the macro global-to-local node map" );

                }
    
                //Set the macro displacements
                for ( unsigned int i = 0; i < nMacroDispDOF; i++ ){
    
                    ( *freeMacroDisplacements )[ nMacroDispDOF * ( map->second ) + i ] = macroDispDOF->second[ i ];
    
                }
    
            }
    
            for ( auto it = freeMicroNodeIds->begin( ); it != freeMicroNodeIds->end( ); it++ ){
   
                auto microDispDOF = microDisplacements->find( *it );

                if ( microDispDOF == microDisplacements->end(  ) ){

                    return new errorNode( "projectDegreesOfFreedom",
                                          "Micro node " + std::to_string( *it ) +
                                          " was not found in the micro displacement dof vector map" );

                }

                if ( microDispDOF->second.size( ) != nMicroDispDOF ){

                    return new errorNode( "projectDegreesOfFreedom",
                                          "Micro node " + std::to_string( *it ) +
                                          " does not have a dimensionally consistent number of degrees of freedom" );

                }

                auto map = microGlobalToLocalDOFMap->find( *it );

                if ( map == microGlobalToLocalDOFMap->end( ) ){

                    return new errorNode( "projectDegreesOfFreedom",
                                          "Micro node " + std::to_string( *it ) +
                                          " was not found in the micro global-to-local node map" );

                }
    
                //Set the micro displacements
                for ( unsigned int i = 0; i < nMicroDispDOF; i++ ){
    
                    ( *freeMicroDisplacements )[ nMicroDispDOF * ( map->second ) + i ] = microDispDOF->second[ i ];
    
                }
    
            }

        }

        //Map the macro and micro free displacements to Eigen matrices
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > Q( freeMicroDisplacements->data(), freeMicroDisplacements->size( ), 1 );
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > D( freeMacroDisplacements->data(), freeMacroDisplacements->size( ), 1 );

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

        if ( ( config[ "projection_type" ].as< std::string >( ).compare( "l2_projection" ) == 0 ) ||
             ( config[ "projection_type" ].as< std::string >( ).compare( "averaged_l2_projection" ) == 0 ) ){

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
                                                                          const std::unordered_map< uIntType, floatVector > &domainReferenceXiVectors,
                                                                          const std::unordered_map< uIntType, floatVector > &domainMicroPositionShapeFunctionValues ){
        /*!
         * Add the current domain information to the direct projection reference values
         *
         * :param const uIntVector &domainNodes: The nodes associated with the current domain
         * :param const uIntVector &macroNodes: The macro nodes associated with the current domain
         * :param const std::unordered_map< uIntType, floatVector > &domainReferenceXiVectors: The relative position
         *     vectors of the domain
         * :param const std::unordered_map< uIntType, floatVector > &domainMicroPositionShapeFunctionValues: The shape function values at the micro
         *     node positions.
         */

        const std::unordered_map< uIntType, floatType > *microDensities = _inputProcessor.getMicroDensities( );
        const std::unordered_map< uIntType, floatType > *microVolumes   = _inputProcessor.getMicroVolumes( );
        const std::unordered_map< uIntType, floatType > *microWeights   = _inputProcessor.getMicroWeights( );

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

            auto microDensity = microDensities->find( *microNode );

            if ( microDensity == microDensities->end( ) ){

                return new errorNode( "addDomainToDirectProjectionReferenceValues",
                                      "Micro node " + std::to_string( *microNode ) + " was not found in the micro density map" );

            }

            auto microVolume = microVolumes->find( *microNode );

            if ( microVolume == microVolumes->end( ) ){

                return new errorNode( "addDomainToDirectProjectionReferenceValues",
                                      "Micro node " + std::to_string( *microNode ) + " was not found in the micro volume map" );

            }

            auto microWeight = microWeights->find( *microNode );

            if ( microWeight == microWeights->end( ) ){

                return new errorNode( "addDomainToDirectProjectionReferenceValues",
                                      "Micro node " + std::to_string( *microNode ) + " was not found in the micro weight map" );

            }

            auto referenceXi = domainReferenceXiVectors.find( *microNode );

            if ( referenceXi == domainReferenceXiVectors.end( ) ){

                return new errorNode( "addDomainToDirectProjectionReferenceValues",
                                      "Micro node " + std::to_string( *microNode ) + " was not found in the reference Xi map" );

            }

            auto shapefunctions = domainMicroPositionShapeFunctionValues.find( *microNode );

            if ( shapefunctions == domainMicroPositionShapeFunctionValues.end( ) ){

                return new errorNode( "addDomainToDirectProjectionReferenceValues",
                                      "Micro node " + std::to_string( *microNode ) + " was not found in the shape function values map" );

            }

            //Compute the micro-mass
            microMass = microDensity->second * microVolume->second;

            //Extract the micro Xi vector
            Xi = referenceXi->second;

            //Extract the weighting value
            weight = microWeight->second;

            //Loop through the macro nodes
            for ( auto macroNode  = macroNodes.begin( );
                       macroNode != macroNodes.end( );
                       macroNode++ ){

                //Set the index
                n = macroNode - macroNodes.begin( );

                auto indx = _inputProcessor.getMacroGlobalToLocalDOFMap( )->find( *macroNode );

                if ( indx == _inputProcessor.getMacroGlobalToLocalDOFMap( )->end( ) ){

                    return new errorNode( "addDomainToDirectProjectionReferenceValues",
                                          "Macro node '" + std::to_string( *macroNode ) + "' not found in global to local macro node map" );

                }
                else{

                    p = indx->second;

                }

                //Get the shape function
                sf = shapefunctions->second[ n ];

                //Add the contribution to the nodal mass
                if ( _macroNodeProjectedMass.find( *macroNode ) == _macroNodeProjectedMass.end( ) ){

                    _macroNodeProjectedMass.emplace( *macroNode, 0. );
                    _macroNodeProjectedMassMomentOfInertia.emplace( *macroNode, floatVector( _dim * _dim, 0 ) );
                    _macroNodeMassRelativePositionConstant.emplace( *macroNode, floatVector( _dim, 0 ) );

                }

                _macroNodeProjectedMass[ *macroNode ] += microMass * sf * weight;
                _macroNodeMassRelativePositionConstant[ *macroNode ] += microMass * sf * weight * Xi;

                for ( unsigned int I = 0; I < _dim; I++ ){

                    for ( unsigned int J = 0; J < _dim; J++ ){

                        _macroNodeProjectedMassMomentOfInertia[ *macroNode ][ _dim * I + J ] += microMass * sf * weight * Xi[ I ] * Xi[ J ];

                    }

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
        errorOut error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, domainName, domainNodes );

        if ( error ){

            errorOut result = new errorNode( "addDomainContributionToDirectFreeMicroToGhostMacroProjector",
                                             "Error in extracting the domain ( " + domainName + " ) nodes" );
            result->addNext( error );
            return result;

        }

        //Get the micro-node positions and relative position vectors
        std::unordered_map< uIntType, floatVector > microNodePositions;
 
        unsigned int index;
        const std::unordered_map< uIntType, floatVector > *microReferencePositions = _inputProcessor.getMicroNodeReferencePositions( );
        const std::unordered_map< uIntType, floatVector > *microDisplacements      = _inputProcessor.getMicroDisplacements( );
        std::unordered_map< uIntType, floatVector > domainReferenceXiVectors;

        auto cellDomainCentersOfMass = _referenceFreeMicroDomainCentersOfMass.find( cellID );
        if ( cellDomainCentersOfMass == _referenceFreeMicroDomainCentersOfMass.end( ) ){

            return new errorNode( "addDomainContributionToDirectFreeMicroToGhostMacroProjector",
                                  "Macro cell " + std::to_string( cellID ) + " not found in reference domain centers of mass map" );

        }

        auto domainCenterOfMass = cellDomainCentersOfMass->second.find( domainName );

        if ( domainCenterOfMass == cellDomainCentersOfMass->second.end( ) ){

            std::string outstr = "Micro domain ";
            outstr += domainName;
            outstr += " not found in the micro domain centers of mass";
            return new errorNode( "addDomainContributionToDirectFreeMicroToGhostMacroProjector", outstr );

        }
 
        for ( auto it = domainNodes.begin( ); it != domainNodes.end( ); it++ ){

            auto microReferencePosition = microReferencePositions->find( *it );

            if ( microReferencePosition == microReferencePositions->end( ) ){

                return new errorNode( "addDomainContributionToDirectFreeMicroToGhostMacroProjector",
                                      "Micro node " + std::to_string( *it ) + " was not found in the micro reference position map" );

            }

            auto microDisplacement = microDisplacements->find( *it );

            if ( microDisplacement == microDisplacements->end( ) ){

                return new errorNode( "addDomainContributionToDirectFreeMicroToGhostMacroProjector",
                                      "Micro node " + std::to_string( *it ) + " was not found in the micro displacement map" );

            }

            microNodePositions.emplace( *it, microReferencePosition->second + microDisplacement->second );

            domainReferenceXiVectors.emplace( *it, microNodePositions[ *it ] - domainCenterOfMass->second );

        }
 
        //Compute the shape function values at the micro positions
        std::unordered_map< uIntType, floatVector > domainMicroPositionShapeFunctionValues;
        error = computeShapeFunctionsAtPoints( cellID,
                                               *_inputProcessor.getMacroNodeReferencePositions( ),
                                               *_inputProcessor.getMacroDisplacements( ),
                                               *_inputProcessor.getMacroNodeReferenceConnectivity( ),
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

            if ( indx == macroGlobalToLocalDOFMap->end( ) ){

                return new errorNode( "addDomainContributionToDirectFreeMicroToGhostMacroProjector",
                                      "'" + std::to_string( *md ) + "' not found in the DOF map" );

            }

            projectorMacroGlobalToLocalDOFMap.emplace( indx->first, indx->second - _inputProcessor.getFreeMacroNodeIds( )->size( ) );

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
                                                                       _macroNodeProjectedMass,
                                                                       _macroNodeProjectedMassMomentOfInertia,
                                                                       _macroNodeMassRelativePositionConstant,
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
//        homogenizedSurfaceRegionDensities.clear( );
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

        const std::unordered_map< uIntType, stringVector > *macroCellToMicroDomainMap = _inputProcessor.getMacroCellToDomainMap( );
        const std::unordered_map< std::string, uIntType > *microDomainSurfaceSplitCount = _inputProcessor.getMicroDomainSurfaceApproximateSplitCount( );

        std::cout << "  looping through the free macro cells\n";
        for ( auto macroCell  = _inputProcessor.getFreeMacroCellIds( )->begin( );
                   macroCell != _inputProcessor.getFreeMacroCellIds( )->end( );
                   macroCell++ ){

            //Set the macro index
            unsigned int macroIndex = macroCell - _inputProcessor.getFreeMacroCellIds( )->begin( );

            //Get the micro domain names within this cell
            auto microDomains = macroCellToMicroDomainMap->find( *macroCell );
            if ( microDomains == macroCellToMicroDomainMap->end( ) ){

                return new errorNode( "homogenizedMicroScale",
                                      "Macro cell " + std::to_string( *macroCell ) + " not found in the macro cell to micro domain map" ) ;

            }

            std::cout << "    looping over the micro domains\n";
            for ( auto microDomain  = microDomains->second.begin( ); microDomain != microDomains->second.end( ); microDomain++ ){

                std::cout << "      " + *microDomain << "\n";
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
                auto domainCenterOfMass = _ghostMicroDomainCentersOfMass.find( *microDomain );
                if ( domainCenterOfMass == _ghostMicroDomainCentersOfMass.end( ) ){

                    std::string outstr = "Ghost micro domain ";
                    outstr += *microDomain;
                    outstr += " not found in the center of mass map";
                    return new errorNode( "homogenizedMicroScale", outstr );

                }

                std::cout << "        computing volume averages\n";
                error = computeDomainVolumeAverages( *macroCell, *microDomain, microDomainNodeIds,
                                                     reconstructedVolume, &domainCenterOfMass->second );

                if ( error ){

                    errorOut result = new errorNode( "homogenizeMicroscale",
                                                     "Error in the computation of the volume averages of the microscale domain" );
                    result->addNext( error );
                    return result;

                }

                auto domainSurfaceCount = microDomainSurfaceSplitCount->find( *microDomain );
                if ( domainSurfaceCount == microDomainSurfaceSplitCount->end( ) ){

                    return new errorNode( "homogenizeMicroscale",
                                          "The micro domain " + *microDomain + " was not found in the domain surface split count map" );

                }
                
                //Compute the surface averages
                std::cout << "        computing surface averages\n";
                error = computeDomainSurfaceAverages( *macroCell, *microDomain, microDomainNodeIds,
                                                      domainSurfaceCount->second,
                                                      reconstructedVolume );

                if ( error ){

                    errorOut result = new errorNode( "homogenizeMicroScale",
                                                     "Error in the computation of the surface averages of the microscale domain" );
                    result->addNext( error );
                    return result;

                }
                
            }

            //Compute the approximate stresses
            std::cout << "    computing the homogenized stresses\n";
            return NULL; //Remove this!
            error = computeHomogenizedStresses( *macroCell );

            if ( error ){

                errorOut result = new errorNode( "homogenizeMicroScale",
                                                 "Error in the computation of the homogenized stresses" );
                result->addNext( error );
                return result;

            }

        }

        //Loop through the ghost macro-scale cells
        std::cout << "  looping through the ghost macro cells\n";
        microDomainStartIndex = 0;
        for ( auto macroCell  = _inputProcessor.getGhostMacroCellIds( )->begin( );
                   macroCell != _inputProcessor.getGhostMacroCellIds( )->end( );
                   macroCell++ ){

            //Set the macro index
            unsigned int macroIndex = macroCell - _inputProcessor.getGhostMacroCellIds( )->begin( );

            //Get the micro domain names within this cell
            auto microDomains = macroCellToMicroDomainMap->find( *macroCell );
            if ( microDomains == macroCellToMicroDomainMap->end( ) ){

                return new errorNode( "homogenizedMicroScale",
                                      "Macro cell " + std::to_string( *macroCell ) + " not found in the macro cell to micro domain map" ) ;
            }

            unsigned int microIndex = microDomainStartIndex;

            std::cout << "    looping over the micro domains\n";
            for ( auto microDomain = microDomains->second.begin( ); microDomain != microDomains->second.end( ); microDomain++ ){

                std::cout << "      " + *microDomain << "\n";
                microNodePositions.clear( );
                reconstructedVolume.reset( );

                //Reconstruct the micro-domain's volume
                std::cout << "        reconstructing the domain\n";
                error = reconstructDomain( microIncrement, *microDomain, microDomainNodeIds, microNodePositions, reconstructedVolume );

                if ( error ){

                    errorOut result = new errorNode( "homogenizeMicroScale",
                                                     "Error in the reconstruction of the microscale domain" );
                    result->addNext( error );
                    return result;

                }

                //Compute the volume averages
                auto domainCenterOfMass = _freeMicroDomainCentersOfMass.find( *microDomain );
                if ( domainCenterOfMass == _freeMicroDomainCentersOfMass.end( ) ){

                    std::string outstr = "Free micro domain ";
                    outstr += *microDomain;
                    outstr += " not found in the center of mass map";
                    return new errorNode( "homogenizedMicroScale", outstr );

                }

                std::cout << "        computing volume averages\n";
                error = computeDomainVolumeAverages( *macroCell, *microDomain, microDomainNodeIds,
                                                     reconstructedVolume, &domainCenterOfMass->second );

                if ( error ){

                    errorOut result = new errorNode( "computeDomainVolumeAverages",
                                                     "Error in the computation of the volume averages of the microscale domain" );
                    result->addNext( error );
                    return result;

                }
                
                auto domainSurfaceCount = microDomainSurfaceSplitCount->find( *microDomain );
                if ( domainSurfaceCount == microDomainSurfaceSplitCount->end( ) ){

                    return new errorNode( "homogenizeMicroscale",
                                          "The micro domain " + *microDomain + " was not found in the domain surface split count map" );

                }

                //Compute the surface averages
                std::cout << "        computing surface averages\n";
                error = computeDomainSurfaceAverages( *macroCell, *microDomain, microDomainNodeIds,
                                                      domainSurfaceCount->second,
                                                      reconstructedVolume );

                if ( error ){

                    errorOut result = new errorNode( "homogenizeMicroScale",
                                                     "Error in the computation of the surface averages of the microscale domain" );
                    result->addNext( error );
                    return result;

                }
                
            }

            //Compute the approximate stresses
            std::cout << "    computing homogenized stresses\n";
            error = computeHomogenizedStresses( *macroCell );

            if ( error ){

                errorOut result = new errorNode( "homogenizeMicroScale",
                                                 "Error in the computation of the homogenized stresses" );
                result->addNext( error );
                return result;

            }

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
        errorOut error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, microDomainName, microDomainNodes );

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
        const std::unordered_map< uIntType, floatVector > *microReferencePositions = _inputProcessor.getMicroNodeReferencePositions( );
        const std::unordered_map< uIntType, floatVector > *microDisplacements      = _inputProcessor.getMicroDisplacements( );
 
        for ( auto it = microDomainNodes.begin( ); it != microDomainNodes.end( ); it++ ){
 
            index = it - microDomainNodes.begin( );

            auto microReferencePosition = microReferencePositions->find( *it );

            if ( microReferencePosition == microReferencePositions->end( ) ){

                return new errorNode( "reconstructDomain", "Micro node " + std::to_string( *it ) +
                                      " was not found in the micro reference position map" );

            }

            auto microDisplacement = microDisplacements->find( *it );

            if ( microDisplacement == microDisplacements->end( ) ){

                return new errorNode( "reconstructDomain", "Micro node " + std::to_string( *it ) +
                                      " was not found in the micro displacement map" );

            }


            for ( unsigned int i = 0; i < _dim; i++ ){

                microNodePositions[ _dim * index + i ] = microReferencePosition->second[ i ] + microDisplacement->second[ i ];

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

    errorOut overlapCoupling::computeDomainVolumeAverages( const uIntType &macroCellID, const std::string &microDomainName,
                                                           const uIntVector &microDomainNodeIDs,
                                                           std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume,
                                                           const floatVector *microDomainCenterOfMass ){
        /*!
         * Compute the required volume averages over the micro-domain.
         *
         * :param const uIntType &macroCellID: The ID number of the macro-cell associated with the micro domain.
         * :param const std::string &microDomainName: The name of the micro domain
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

        const std::unordered_map< uIntType, floatType > *microDensities = _inputProcessor.getMicroDensities( );

        const std::unordered_map< uIntType, floatVector > *microBodyForces = _inputProcessor.getMicroBodyForces( );

        const std::unordered_map< uIntType, floatVector > *microAccelerations = _inputProcessor.getMicroAccelerations( );

        const std::unordered_map< uIntType, floatVector > *microReferencePositions = _inputProcessor.getMicroNodeReferencePositions( );

        const std::unordered_map< uIntType, floatVector > *microDisplacements = _inputProcessor.getMicroDisplacements( );

        const std::unordered_map< uIntType, floatVector > *microStresses = _inputProcessor.getMicroStresses( );

        unsigned int index = 0;
        unsigned int localIndex = 0;

        for ( auto node = microDomainNodeIDs.begin( ); node != microDomainNodeIDs.end( ); node++, index++ ){

            dataAtMicroPoints[ dataCountAtPoint * index + 0 ] = 1.;                           //Integrate the volume of the domain

            auto microDensity = microDensities->find( *node );
            if ( microDensity == microDensities->end( ) ){
                return new errorNode( "computeDomainVolumeAverages",
                                      "Micro node " + std::to_string( *node ) + " was not found in the micro density map" );
            }

            dataAtMicroPoints[ dataCountAtPoint * index + 1 ] = microDensity->second; //Integrate the density of the domain

            //Integrate the micro stresses
            auto microStress = microStresses->find( *node );
            if ( microStress == microStresses->end( ) ){
                return new errorNode( "computeDomainVolumeAverages",
                                      "Micro node " + std::to_string( *node ) + " was not found in the micro stress map" );
            }
            for ( unsigned int i = 0; i < _dim * _dim; i++ ){
                dataAtMicroPoints[ dataCountAtPoint * index + 2 + i ] = microStress->second[ i ];
            }

            localIndex = initialOffset;

            if ( _inputProcessor.useReconstructedMassCenters( ) ){

                auto microReferencePosition = microReferencePositions->find( *node );
                if ( microReferencePosition == microReferencePositions->end( ) ){
                    return new errorNode( "computeDomainVolumeAverages",
                                          "Micro node " + std::to_string( *node ) + " was not found in the micro reference position map" );
                }

                auto microDisplacement = microDisplacements->find( *node );
                if ( microDisplacement == microDisplacements->end( ) ){
                    return new errorNode( "computeDomainVolumeAverages",
                                          "Micro node " + std::to_string( *node ) + " was not found in the micro displacement map" );
                }

                //Integrate for the domain's center of mass
                for ( unsigned int i = 0; i < _dim; i++ ){
    
                    dataAtMicroPoints[ dataCountAtPoint * index + localIndex + i ] =
                        microDensity->second * ( microReferencePosition->second[ i ] + microDisplacement->second[ i ] );
    
                }

                //Local index
                localIndex += _dim;

            }

            //Add the micro body forces
            if ( _inputProcessor.microBodyForceDefined( ) ){

                auto microBodyForce = microBodyForces->find( *node );
                if ( microBodyForce == microBodyForces->end( ) ){
                    return new errorNode( "computeDomainVolumeAverages",
                                          "Micro node " + std::to_string( *node ) + " was not found in the micro body force map" );
                }

                for ( unsigned int i = 0; i < _dim; i++ ){

                    dataAtMicroPoints[ dataCountAtPoint * index + localIndex + i ]
                        = microDensity->second * microBodyForce->second[ i ]; //Integrate the body forces of the domain

                }

                localIndex += _dim;

            }

            //Add the micro acclerations
            if ( _inputProcessor.microAccelerationDefined( ) ){

                auto microAcceleration = microAccelerations->find( *node );
                if ( microAcceleration == microAccelerations->end( ) ){
                    return new errorNode( "computeDomainVolumeAverages",
                                          "Micro node " + std::to_string( *node ) + " was not found in the micro acceleration map" );
                }

                for ( unsigned int i = 0; i < _dim; i++ ){

                    dataAtMicroPoints[ dataCountAtPoint * index + localIndex + i ]
                        = microDensity->second * microAcceleration->second[ i ]; //Integrate the accelerations of the domain

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

        //Initialize
        if ( homogenizedVolumes.find( macroCellID ) == homogenizedVolumes.end( ) ){

            domainFloatMap tmpFloat;
            domainFloatVectorMap tmpFloatVector;

            homogenizedVolumes.emplace( macroCellID, tmpFloat );
            homogenizedDensities.emplace( macroCellID, tmpFloat );
            homogenizedSymmetricMicroStresses.emplace( macroCellID, tmpFloatVector );
            homogenizedCentersOfMass.emplace( macroCellID, tmpFloatVector );
            homogenizedBodyForces.emplace( macroCellID, tmpFloatVector );
            homogenizedAccelerations.emplace( macroCellID, tmpFloatVector );
            homogenizedMicroInertias.emplace( macroCellID, tmpFloatVector );
            homogenizedBodyForceCouples.emplace( macroCellID, tmpFloatVector );
            homogenizedMicroSpinInertias.emplace( macroCellID, tmpFloatVector );

        }

        homogenizedVolumes[ macroCellID ].emplace( microDomainName, integratedValues[ 0 ] );
        homogenizedDensities[ macroCellID ].emplace( microDomainName, integratedValues[ 1 ] / integratedValues[ 0 ] );
        homogenizedSymmetricMicroStresses[ macroCellID ].emplace( microDomainName, floatVector( integratedValues.begin( ) + 2,
                                                                                   integratedValues.begin( ) + 2 + _dim * _dim ) / integratedValues[ 0 ] );
        homogenizedCentersOfMass[ macroCellID ].emplace( microDomainName, floatVector( integratedValues.begin( ) + initialOffset,
                                                                                       integratedValues.begin( ) + initialOffset + _dim ) / integratedValues[ 1 ] );
        homogenizedBodyForces[ macroCellID ].emplace( microDomainName, floatVector( integratedValues.begin( ) + initialOffset + _dim,
                                                                                    integratedValues.begin( ) + initialOffset + 2 * _dim ) / integratedValues[ 1 ] );
        homogenizedAccelerations[ macroCellID ].emplace( microDomainName, floatVector( integratedValues.begin( ) + initialOffset + 2 * _dim,
                                                                                       integratedValues.end( ) ) / integratedValues[ 1 ] );


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
            floatVector centerOfMass = homogenizedCentersOfMass[ macroCellID ][ microDomainName ];
    
            for ( auto node = microDomainNodeIDs.begin( ); node != microDomainNodeIDs.end( ); node++, index++ ){

                auto microDensity = microDensities->find( *node );
                if ( microDensity == microDensities->end( ) ){
                    return new errorNode( "computeDomainVolumeAverages",
                                          "Micro node " + std::to_string( *node ) + " was not found in the micro density map" );
                }

                auto microReferencePosition = microReferencePositions->find( *node );
                if ( microReferencePosition == microReferencePositions->end( ) ){
                    return new errorNode( "computeDomainVolumeAverages",
                                          "Micro node " + std::to_string( *node ) + " was not found in the micro reference position map" );
                }

                auto microDisplacement = microDisplacements->find( *node );
                if ( microDisplacement == microDisplacements->end( ) ){
                    return new errorNode( "computeDomainVolumeAverages",
                                          "Micro node " + std::to_string( *node ) + " was not found in the micro displacement map" );
                }

                //Extract the micro relative position
                microRelativePosition = microReferencePosition->second + microDisplacement->second - centerOfMass;

                floatVector integrand
                    = microDensity->second
                    * vectorTools::appendVectors( vectorTools::dyadic( microRelativePosition, microRelativePosition ) );

                //Add the contributions to the micro inertia
                for ( unsigned int i = 0; i < _dim * _dim; i++ ){

                    dataAtMicroPoints[ dataCountAtPoint * index + i ]
                        = integrand[ i ];

                }

                localIndex = initialOffset;

                //Add the contributions to the micro body couple
                if ( _inputProcessor.microBodyForceDefined( ) ){

                    auto microBodyForce = microBodyForces->find( *node );
                    if ( microBodyForce == microBodyForces->end( ) ){
                        return new errorNode( "computeDomainVolumeAverages",
                                              "Micro node " + std::to_string( *node ) + " was not found in the micro body force map" );
                    }

                    floatVector integrand
                        = microDensity->second
                        * vectorTools::appendVectors( vectorTools::dyadic( microBodyForce->second, microRelativePosition ) );
    
                    for ( unsigned int i = 0; i < _dim * _dim; i++ ){
    
                        dataAtMicroPoints[ dataCountAtPoint * index + localIndex + i ]
                            = integrand[ i ]; //Integrate the body force couple over the domain
    
                    }
    
                    localIndex += _dim * _dim;
    
                }
    
                //Add the contributions to the micro spin inertia
                if ( _inputProcessor.microAccelerationDefined( ) ){

                    auto microAcceleration = microAccelerations->find( *node );
                    if ( microAcceleration == microAccelerations->end( ) ){
                        return new errorNode( "computeDomainVolumeAverages",
                                              "Micro node " + std::to_string( *node ) + " was not found in the micro acceleration map" );
                    }

                    floatVector microRelativeAcceleration = microAcceleration->second
                                                          - homogenizedAccelerations[ macroCellID ][ microDomainName ];

                    floatVector integrand
                        = microDensity->second
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

        floatType domainMass = homogenizedDensities[ macroCellID ][ microDomainName ]
                             * homogenizedVolumes[ macroCellID ][ microDomainName ];
        homogenizedMicroInertias[ macroCellID ].emplace( microDomainName, floatVector( integratedValues.begin( ),
                                                                                       integratedValues.begin( ) + initialOffset ) / domainMass );
        homogenizedBodyForceCouples[ macroCellID ].emplace( microDomainName, floatVector( integratedValues.begin( ) + initialOffset,
                                                                                          integratedValues.begin( ) + initialOffset + _dim * _dim ) / domainMass );
        homogenizedMicroSpinInertias[ macroCellID ].emplace( microDomainName, floatVector( integratedValues.begin( ) + initialOffset + _dim * _dim,
                                                                                           integratedValues.begin( ) + initialOffset + 2 * _dim * _dim ) / domainMass );

        return NULL;

    }

    errorOut overlapCoupling::computeDomainSurfaceAverages( const uIntType &macroCellID, const std::string &microDomainName,
                                                            const uIntVector &microDomainNodeIDs,
                                                            const uIntType &microDomainSurfaceDecompositionCount,
                                                            std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume ){
        /*!
         * Compute the required surface averages over the micro-domain.
         *
         * :param const uIntType &macroCellID: The ID of the macro-cell associated with the micro domain.
         * :param const std::string &microDomainName: The name of the micro domain
         * :param const std::string &microDomainNodeIDs: The IDs of the nodes in the micro-domain to have the surface averages
         *     computed over.
         * :param const uIntType &microDomainSurfaceDecompositionCount: The approximate number of regions to split the surface
         *     of the micro surface into.
         * :param volumeReconstruction::volumeReconstructionBase &reconstructedVolume: The reconstructed volume
         *     ready to have surface integrals computed over.
         */

        //Extract the required micro-scale values
        const std::unordered_map< uIntType, floatType > *microDensities = _inputProcessor.getMicroDensities( );
        const std::unordered_map< uIntType, floatVector > *microDisplacements = _inputProcessor.getMicroDisplacements( );
        const std::unordered_map< uIntType, floatVector > *microReferencePositions = _inputProcessor.getMicroNodeReferencePositions( );
        const std::unordered_map< uIntType, floatVector > *microStresses = _inputProcessor.getMicroStresses( );

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

            domainFloatMap tmpFloatMap;
            domainFloatVectorMap tmpFloatVectorMap;
            homogenizedSurfaceAreas.emplace( macroCellID, tmpFloatMap );
            homogenizedSurfaceRegionAreas.emplace( macroCellID, tmpFloatVectorMap );
            homogenizedSurfaceRegionTractions.emplace( macroCellID, tmpFloatVectorMap );
            homogenizedSurfaceRegionCouples.emplace( macroCellID, tmpFloatVectorMap );

        }

        homogenizedSurfaceAreas[ macroCellID ].emplace( microDomainName, integratedValue[ 0 ] );

        /*=====================================================================
        |           Compute the properties of the surface subdomains          |
        =====================================================================*/

        uIntVector subdomainNodeCounts;
        uIntVector subdomainNodeIDs;

        floatType minSurfaceSpacing
            = std::sqrt( homogenizedSurfaceAreas[ macroCellID ][ microDomainName ] / ( 3.14159 * microDomainSurfaceDecompositionCount ) );

        error = reconstructedVolume->getSurfaceSubdomains( minSurfaceSpacing, subdomainNodeCounts, subdomainNodeIDs );

        if ( error ){

            errorOut result = new errorNode( "computeDomainSurfaceAverages", "Error in extracting of the reconstructed volume's surface subdomains" );
            result->addNext( error );
            return result;

        }

        //Get the centers of mass of the surface regions

        dataCountAtPoint = 1     //Surface area of region
                         + 1     //Mass of region
                         + _dim; //Micro point position

        dataAtMicroPoints.clear( );
        dataAtMicroPoints.reserve( dataCountAtPoint * microDomainNodeIDs.size( ) );

        for ( auto node = microDomainNodeIDs.begin( ); node != microDomainNodeIDs.end( ); node++ ){

            auto microDensity = microDensities->find( *node );
            if ( microDensity == microDensities->end( ) ){
                return new errorNode( "computeDomainSurfaceAverages",
                                      "The micro node " + std::to_string( *node ) + " was not found in the micro density map" );
            }

            dataAtMicroPoints.push_back( 1 );
            dataAtMicroPoints.push_back( microDensity->second );

            auto microReferencePosition = microReferencePositions->find( *node );
            if ( microReferencePosition == microReferencePositions->end( ) ){
                return new errorNode( "computeDomainVolumeAverages",
                                      "Micro node " + std::to_string( *node ) + " was not found in the micro reference position map" );
            }

            auto microDisplacement = microDisplacements->find( *node );
            if ( microDisplacement == microDisplacements->end( ) ){
                return new errorNode( "computeDomainVolumeAverages",
                                      "Micro node " + std::to_string( *node ) + " was not found in the micro displacement map" );
            }

            floatVector microPoint = microReferencePosition->second + microDisplacement->second;

            for ( unsigned int i = 0; i < microPoint.size( ); i++ ){

                dataAtMicroPoints.push_back( microDensity->second * microPoint[ i ] );

            }

        }

        unsigned int startPoint = 0;

        //Initialize storage values for homogenization
        homogenizedSurfaceRegionAreas[ macroCellID ].emplace( microDomainName, floatVector( subdomainNodeCounts.size( ), 0 ) );
        floatVector regionDensities( subdomainNodeCounts.size( ) );
        homogenizedSurfaceRegionCentersOfMass[ macroCellID ].emplace( microDomainName, floatVector( _dim * subdomainNodeCounts.size( ), 0 ) );

        //Initialize the homogenized values
        homogenizedSurfaceRegionTractions[ macroCellID ].emplace( microDomainName, floatVector( _dim * subdomainNodeCounts.size( ), 0 ) );
        homogenizedSurfaceRegionCouples[ macroCellID ].emplace( microDomainName, floatVector( _dim * _dim * subdomainNodeCounts.size( ), 0 ) );

        uIntType index = 0;
        for ( auto sNC = subdomainNodeCounts.begin( ); sNC != subdomainNodeCounts.end( ); sNC++, index++ ){

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
            homogenizedSurfaceRegionAreas[  macroCellID ][ microDomainName ][ index ] = integratedValue[ 0 ];
            regionDensities[ index ] = integratedValue[ 1 ] / integratedValue[ 0 ];

            for ( unsigned int i = 0; i < _dim; i++ ){
                homogenizedSurfaceRegionCentersOfMass[ macroCellID ][ microDomainName ][ _dim * index + i ]
                    = integratedValue[ 2 + i ] / integratedValue[ 1 ];
            }

            startPoint += *sNC;
        }


        /*===================================================================================
        | Compute the surface tractions and couples over the micro domain's surface regions |
        ===================================================================================*/

        dataCountAtPoint = _dim * _dim;

        dataAtMicroPoints.clear( );
        dataAtMicroPoints.reserve( dataCountAtPoint * microDomainNodeIDs.size( ) );

        for ( auto node = microDomainNodeIDs.begin( ); node != microDomainNodeIDs.end( ); node++ ){

            auto microStress = microStresses->find( *node );
            if ( microStress == microStresses->end( ) ){
                return new errorNode( "computeDomainVolumeAverages",
                                      "Micro node " + std::to_string( *node ) + " was not found in the micro stress map" );
            }

            for ( unsigned int i = 0; i < _dim * _dim; i++ ){

                dataAtMicroPoints.push_back( microStress->second[ i ] );

            }

        }

        startPoint = 0;

        index = 0;
        for ( auto sNC = subdomainNodeCounts.begin( ); sNC != subdomainNodeCounts.end( ); sNC++, index++ ){

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

//            unsigned int regionID = homogenizedSurfaceRegionTractions[ macroCellID ].size( ) / _dim;
//
//            floatType regionSurfaceArea = homogenizedSurfaceRegionAreas[ macroCellID ][ regionID ];
//
//            integratedValue /= regionSurfaceArea;
//
//            homogenizedSurfaceRegionTractions[ macroCellID ]
//                = vectorTools::appendVectors( { homogenizedSurfaceRegionTractions[ macroCellID ],
//                                                integratedValue } );
            for ( uIntType i = 0; i < _dim; i++ ){
                homogenizedSurfaceRegionTractions[ macroCellID ][ microDomainName ][ _dim * index + i ]
                    = integratedValue[ i ] / homogenizedSurfaceRegionAreas[ macroCellID ][ microDomainName ][ index ];
            }

            //Compute the couples
            floatVector regionCenterOfMass ( homogenizedSurfaceRegionCentersOfMass[ macroCellID ][ microDomainName ].begin( ) + _dim * index,
                                             homogenizedSurfaceRegionCentersOfMass[ macroCellID ][ microDomainName ].begin( ) + _dim * ( index + 1 ) );

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

//            integratedValue /= regionSurfaceArea;
//
//            homogenizedSurfaceRegionCouples[ macroCellID ]
//                = vectorTools::appendVectors( { homogenizedSurfaceRegionCouples[ macroCellID ],
//                                                integratedValue } );

            for ( uIntType i = 0; i < _dim * _dim; i++ ){

                homogenizedSurfaceRegionCouples[ macroCellID ][ microDomainName ][ _dim * _dim * index + i ]
                    = integratedValue[ i ] / homogenizedSurfaceRegionAreas[ macroCellID ][ microDomainName ][ index ];

            }

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

        //Get the pointers to the values
        const std::unordered_map< uIntType, floatVector > *macroNodeReferenceLocations
            = _inputProcessor.getMacroNodeReferencePositions( );
        const std::unordered_map< uIntType, floatVector > *macroDisplacements
            = _inputProcessor.getMacroDisplacements( );
        const std::unordered_map< uIntType, uIntVector > *macroConnectivity
            = _inputProcessor.getMacroNodeReferenceConnectivity( );

        //Form the finite element representation of the macro-scale
        std::unique_ptr< elib::Element > element;
        errorOut error = overlapCoupling::buildMacroDomainElement( macroCellID, *macroNodeReferenceLocations,
                                                                   *macroDisplacements, *macroConnectivity,
                                                                   element );

        if ( error ){

            errorOut result = new errorNode( "computeHomogenizedStresses",
                                             "Error in the formation of the finite element representation of the macro-scale" );
            result->addNext( error );
            return result;

        }

        //Get the shape functions at the micro-domain centroids
        std::unordered_map< uIntType, floatVector > centerOfMassMap;
        stringVector domainNames;
        centerOfMassMap.reserve( homogenizedCentersOfMass[ macroCellID ].size( ) );

        uIntType index = 0;
        for ( auto name  = homogenizedCentersOfMass[ macroCellID ].begin( );
                   name != homogenizedCentersOfMass[ macroCellID ].end( ); name++, index++ ){

            domainNames.push_back( name->first );
            centerOfMassMap.emplace( index, name->second );

        }
        

        std::unordered_map< uIntType, floatVector > shapefunctionsAtCentersOfMassByID;
        error = overlapCoupling::computeShapeFunctionsAtPoints( macroCellID, *macroNodeReferenceLocations, *macroDisplacements,
                                                                *macroConnectivity,
                                                                centerOfMassMap,
                                                                shapefunctionsAtCentersOfMassByID );

        if ( error ){

            errorOut result = new errorNode( "computeHomogenizedStresses",
                                             "Error in the computation of the shapefunctions at the micro domain centers of mass for macro cell " + std::to_string( macroCellID ) );
            result->addNext( error );
            return result;

        }

        domainFloatVectorMap shapefunctionsAtCentersOfMass;
        index = 0;
        for ( auto name = domainNames.begin( ); name != domainNames.end( ); name++, index++ ){

            shapefunctionsAtCentersOfMass.emplace( *name, shapefunctionsAtCentersOfMassByID[ index ] );

            if ( shapefunctionsAtCentersOfMassByID[ index ].size( ) != element->nodes.size( ) ){

                std::string output;
                output += "The number of shape-function defined is not consistent with the number of micro domains\n";
                output += "and the number of nodes in the macro element for macro-cell " + std::to_string( macroCellID ) + ".\n";
                output += "This is likely because one of the micro-domains center of mass is located outside of the macro cell";
    
                return new errorNode( "computeHomogenizedStresses", output );

            }

        }

        uIntType nMacroCellNodes = element->nodes.size( );

        floatVector linearMomentumRHS( _dim * nMacroCellNodes, 0 );
        floatVector firstMomentRHS( _dim * _dim * nMacroCellNodes, 0 );

        //Project values to the nodes
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
        for ( auto microDomainName = domainNames.begin( ); microDomainName != domainNames.end( ); microDomainName++ ){

            floatType density = homogenizedDensities[ macroCellID ][ *microDomainName ];

            floatType volume = homogenizedVolumes[ macroCellID ][ *microDomainName ];

            floatVector bodyForce = homogenizedBodyForces[ macroCellID ][ *microDomainName ];

            floatVector acceleration = homogenizedAccelerations[ macroCellID ][ *microDomainName ];

            floatVector microInertia = homogenizedMicroInertias[ macroCellID ][ *microDomainName ];

            floatVector bodyCouple = homogenizedBodyForceCouples[ macroCellID ][ *microDomainName ];

            floatVector microSpinInertia = homogenizedMicroSpinInertias[ macroCellID ][ *microDomainName ];

            floatVector symmetricMicroStress = homogenizedSymmetricMicroStresses[ macroCellID ][ *microDomainName ];

            floatVector symmetricMicroStress_T( _dim * _dim );

            for ( unsigned int _i = 0; _i < _dim; _i++ ){

                for ( unsigned int _j = 0; _j < _dim; _j++ ){

                    symmetricMicroStress_T[ _dim * _j + _i ] = symmetricMicroStress[ _dim * _i + _j ];

                }

            }

            for ( unsigned int j = 0; j < nMacroCellNodes; j++ ){

                //Get the shapefunction value for the node
                floatType N = shapefunctionsAtCentersOfMass[ *microDomainName ][ j ];

                //Compute the contribution to the node
                floatVector nLinearMomentumRHS = N * density * ( bodyForce - acceleration ) * volume;

                floatVector nFirstMomentRHS = N * ( density * ( bodyCouple - microSpinInertia ) - symmetricMicroStress_T ) * volume;

                //Add the contribution to the overall RHS vectors
                for ( auto it = nLinearMomentumRHS.begin( ); it != nLinearMomentumRHS.end( ); it++ ){

                    uIntType index = _dim * j + ( it - nLinearMomentumRHS.begin( ) );

                    linearMomentumRHS[ index ] += *it;

                }

                for ( auto it = nFirstMomentRHS.begin( ); it != nFirstMomentRHS.end( ); it++ ){

                    uIntType index = _dim * _dim * j + ( it - nFirstMomentRHS.begin( ) );

                    firstMomentRHS[ index ] += *it;

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
        for ( auto domain = domainNames.begin( ); domain != domainNames.end( ); domain++ ){

            //Get the number of micro surface regions in the domain
            uIntType nMicroSurfaceRegions = homogenizedSurfaceRegionAreas[ macroCellID ][ *domain ].size( );

            //Compute the shape functions at the micro surface region area centers of mass
            std::unordered_map< uIntType, floatVector > domainCentersOfMass;
            for ( unsigned int i = 0; i < nMicroSurfaceRegions; i++ ){

                domainCentersOfMass.emplace( i, floatVector( homogenizedSurfaceRegionCentersOfMass[ macroCellID ][ *domain ].begin( ) + _dim * i,
                                                             homogenizedSurfaceRegionCentersOfMass[ macroCellID ][ *domain ].begin( ) + _dim * ( i + 1 ) ) );

            }

            std::unordered_map< uIntType, floatVector > shapefunctionsAtSurfaceRegionCentersOfMass;

            //TODO: The following computation of the shape-function values may cause issues if the re-constructed surface is
            //      found to be outside of the macro Cell's domain. This will likely need to be addressed
            error = overlapCoupling::computeShapeFunctionsAtPoints( macroCellID, *macroNodeReferenceLocations, *macroDisplacements,
                                                                *macroConnectivity, domainCentersOfMass,
                                                                shapefunctionsAtSurfaceRegionCentersOfMass );
            
            if ( error ){
    
                errorOut result = new errorNode( "computeHomogenizedStresses",
                                                 "Error in the computation of the shapefunctions at the micro domain surface region centers of mass for macro cell " + std::to_string( macroCellID ) );
                result->addNext( error );
                return result;
    
            }
    
            for ( unsigned int i = 0; i < nMicroSurfaceRegions; i++ ){
    
                if ( shapefunctionsAtSurfaceRegionCentersOfMass[ index ].size( ) != element->nodes.size( ) ){
    
                    std::string output;
                    output += "The number of shape-function defined is not consistent with the number of micro domains\n";
                    output += "and the number of nodes in the macro element for macro-cell " + std::to_string( macroCellID ) + ".\n";
                    output += "This is likely because one of the surface region's center of mass is located outside of the macro cell";
        
                    return new errorNode( "computeHomogenizedStresses", output );
    
                }
    
            }

            //Add the surface integral components to the reight hand side vectors
            for ( unsigned int i = 0; i < nMicroSurfaceRegions; i++ ){

                floatType area = homogenizedSurfaceRegionAreas[ macroCellID ][ *domain ][ i ];

                floatVector traction( homogenizedSurfaceRegionTractions[ macroCellID ][ *domain ].begin( ) + _dim * i,
                                      homogenizedSurfaceRegionTractions[ macroCellID ][ *domain ].begin( ) + _dim * ( i + 1 ) );

                floatVector couple( homogenizedSurfaceRegionTractions[ macroCellID ][ *domain ].begin( ) + _dim * _dim * i,
                                    homogenizedSurfaceRegionTractions[ macroCellID ][ *domain ].begin( ) + _dim * _dim * ( i + 1 ) );

                floatVector shapefunctions = shapefunctionsAtSurfaceRegionCentersOfMass[ i ];

                for ( unsigned int j = 0; j < nMacroCellNodes; j++ ){

                    //Get the shapefunction value for the surface region
                    floatType N = shapefunctions[ j ];

                    //Compute the contribution to the node
                    floatVector nLinearMomentumRHS = N * traction * area;

                    floatVector nFirstMomentRHS = N * couple * area;

                    //Add the contribution to the overall RHS vectors
                    for ( auto it = nLinearMomentumRHS.begin( ); it != nLinearMomentumRHS.end( ); it++ ){
    
                        uIntType index = _dim * j + ( it - nLinearMomentumRHS.begin( ) );
    
                        linearMomentumRHS[ index ] += *it;
                        externalForcesAtNodes[ macroCellID ][ index ] += *it;
    
                    }
    
                    for ( auto it = nFirstMomentRHS.begin( ); it != nFirstMomentRHS.end( ); it++ ){
    
                        uIntType index = _dim * _dim * j + ( it - nFirstMomentRHS.begin( ) );
    
                        firstMomentRHS[ index ] += *it;
                        externalCouplesAtNodes[ macroCellID ][ index ] += *it;
    
                    }

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

        tripletVector coefficients;
        coefficients.reserve( ( 2 * _dim * _dim + 3 * _dim * _dim ) * element->nodes.size( ) * element->qrule.size( ) );

        //Quadrature point interpolated values
        floatVector densities( element->qrule.size( ), 0 );
        floatMatrix bodyForces( element->qrule.size( ), floatVector( _dim, 0 ) );
        floatMatrix accelerations( element->qrule.size( ), floatVector( _dim, 0 ) );
        floatMatrix microInertias( element->qrule.size( ), floatVector( _dim * _dim, 0 ) );
        floatMatrix bodyCouples( element->qrule.size( ), floatVector( _dim * _dim, 0 ) );
        floatMatrix microSpinInertias( element->qrule.size( ), floatVector( _dim * _dim, 0 ) );
        floatMatrix symmetricMicroStress( element->qrule.size( ), floatVector( _dim * _dim, 0 ) );
        floatMatrix nodalDOFValues( element->nodes.size( ), floatVector( _dim + _dim * _dim, 0 ) );

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
        uIntVector outliers;

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

    const cellDomainFloatMap* overlapCoupling::getReferenceFreeMicroDomainMasses( ){
        /*!
         * Get access to the reference free micro-domain mass
         */

        return &_referenceFreeMicroDomainMasses;
    }

    const cellDomainFloatMap* overlapCoupling::getReferenceGhostMicroDomainMasses( ){
        /*!
         * Get access to the reference ghost micro-domain masses
         */

        return &_referenceGhostMicroDomainMasses;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getReferenceFreeMicroDomainCentersOfMass( ){
        /*!
         * Get access to the reference free micro-domain centers of mass
         */

        return &_referenceFreeMicroDomainCentersOfMass;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getReferenceGhostMicroDomainCentersOfMass( ){
        /*!
         * Get access to the reference ghost micro-domain centers of mass
         */

        return &_referenceGhostMicroDomainCentersOfMass;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getReferenceFreeMicroDomainMomentsOfInertia( ){
        /*!
         * Get access to the reference free micro-domain moments of inertia
         */

        return &_referenceFreeMicroDomainMomentsOfInertia;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getReferenceGhostMicroDomainMomentsOfInertia( ){
        /*!
         * Get access to the reference ghost micro-domain moments of inertia
         */

        return &_referenceGhostMicroDomainMomentsOfInertia;
    }

    const domainFloatMap* overlapCoupling::getFreeMicroDomainMasses( ){
        /*!
         * Get access to the free micro-domain mass
         */

        return &_freeMicroDomainMasses;
    }

    const domainFloatMap* overlapCoupling::getGhostMicroDomainMasses( ){
        /*!
         * Get access to the ghost micro-domain masses
         */

        return &_ghostMicroDomainMasses;
    }

    const domainFloatVectorMap* overlapCoupling::getFreeMicroDomainCentersOfMass( ){
        /*!
         * Get access to the free micro-domain centers of mass
         */

        return &_freeMicroDomainCentersOfMass;
    }

    const domainFloatVectorMap* overlapCoupling::getGhostMicroDomainCentersOfMass( ){
        /*!
         * Get access to the ghost micro-domain centers of mass
         */

        return &_ghostMicroDomainCentersOfMass;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getReferenceFreeMicroDomainCenterOfMassShapeFunctions( ){
        /*!
         * Get access to the shapefunction values of the reference free micro domain centers of mass
         */

        return &_referenceFreeMicroDomainCenterOfMassShapeFunctions;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getReferenceGhostMicroDomainCenterOfMassShapeFunctions( ){
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

        floatVector absDeviations = vectorTools::abs(x - median);

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
                                                tripletVector &coefficients ){
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
         * :param tripletVector &coefficients: The coefficients of the mass matrix
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

        const std::unordered_map< uIntType, floatVector > *macroDispDOFVector = _inputProcessor.getMacroDispDOFVector( );
        for ( auto it = freeMacroNodeIds->begin( ); it != freeMacroNodeIds->end( ); it++ ){

            auto map = _inputProcessor.getMacroGlobalToLocalDOFMap( )->find( *it );

            if ( map == _inputProcessor.getMacroGlobalToLocalDOFMap( )->end( ) ){

                return new errorNode( "assembleHomogenizedInternalForceVector",
                                      "Global degree of freedom '" + std::to_string( *it ) + "' not found in degree of freedom map" );

            }

            auto macroDisplacement = macroDispDOFVector->find( *it );

            if ( macroDisplacement == macroDispDOFVector->end( ) ){

                return new errorNode( "assembleHomogenizedInternalForceVector",
                                      "Global macro degree of freedom '" + std::to_string( *it ) + "' not found in the macro displacement dof map" );

            }

            //Set the macro displacements
            for ( unsigned int i = 0; i < nMacroDispDOF; i++ ){

                freeMacroDisplacements[ nMacroDispDOF * ( map->second ) + i ] = macroDisplacement->second[ i ];

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

        const std::unordered_map< uIntType, floatVector > *macroDispDOFVector = _inputProcessor.getMacroDispDOFVector( );
        for ( auto it = freeMacroNodeIds->begin( ); it != freeMacroNodeIds->end( ); it++ ){

            auto map = _inputProcessor.getMacroGlobalToLocalDOFMap( )->find( *it );

            if ( map == _inputProcessor.getMacroGlobalToLocalDOFMap( )->end( ) ){

                return new errorNode( "assembleHomogenizedInternalForceVector",
                                      "Global degree of freedom '" + std::to_string( *it ) + "' not found in degree of freedom map" );

            }

            auto macroDisplacement = macroDispDOFVector->find( *it );

            if ( macroDisplacement == macroDispDOFVector->end( ) ){

                return new errorNode( "assembleHomogenizedInternalForceVector",
                                      "Global macro degree of freedom '" + std::to_string( *it ) + "' not found in the macro displacement dof map" );

            }

            //Set the macro displacements
            for ( unsigned int i = 0; i < nMacroDispDOF; i++ ){

                freeMacroDisplacements[ nMacroDispDOF * ( map->second ) + i ] = macroDisplacement->second[ i ];

            }

        }

        //Loop over the macro cells adding the contributions to the mass matrix
        floatVector elementDOFVector;

        //Collect all of the cells in the overlapping domain
        const uIntVector *freeMacroCellIds = _inputProcessor.getFreeMacroCellIds( );
        const uIntVector *ghostMacroCellIds = _inputProcessor.getGhostMacroCellIds( );

        uIntVector macroCellIDVector( freeMacroCellIds->begin( ), freeMacroCellIds->end( ) );
        macroCellIDVector = vectorTools::appendVectors( { macroCellIDVector, *ghostMacroCellIds } );

        tripletVector coefficients;

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
        tripletVector coefficients;

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
                                                      element );

            if ( error ){

                errorOut result = new errorNode( "assembleFreeMicromorphicMassMatrix",
                                                 "Error in the construction of the macro element " + std::to_string( *macroCellID ) );
                result->addNext( error );
                return result;

            }

            //Get the degree of freedom values for the free macro-domain element
            const std::unordered_map< uIntType, floatVector >* macroDispDOFVector = _inputProcessor.getMacroDispDOFVector( );
            floatVector elementDOFVector( 0 );

            for ( auto nodeID  = element->global_node_ids.begin( );
                       nodeID != element->global_node_ids.end( );
                       nodeID++ ){

                auto macroDisplacement = macroDispDOFVector->find( *nodeID );

                if ( macroDisplacement == macroDispDOFVector->end( ) ){

                    return new errorNode( "assembleFreeMicromorphicMassMatrix",
                                          "Macro node " + std::to_string( *nodeID ) + " was not found in the macro displacement DOF vector map" );

                }

                elementDOFVector
                    = vectorTools::appendVectors( { elementDOFVector, macroDisplacement->second } );

            }

            //Extract the density and moment of inertia in the reference configuration
            auto densityType = macroReferenceDensityTypes->find( *macroCellID );

            if ( densityType == macroReferenceDensityTypes->end( ) ){

                return new errorNode( "assembleFreeMicromorphicMassMatrix",
                                      "The macro cell with ID " + std::to_string( *macroCellID ) +
                                      " was not found in the density type map" );

            }

            auto momentOfInertiaType = macroReferenceMomentOfInertiaTypes->find( *macroCellID );

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
        floatType aQ = config[ "micro_proportionality_coefficient" ].as< floatType >( );
        floatType aD = config[ "macro_proportionality_coefficient" ].as< floatType >( );

        //Get the micro densities and volumes
        const std::unordered_map< uIntType, floatType > *microVolumes   = _inputProcessor.getMicroVolumes( );
        const std::unordered_map< uIntType, floatType > *microDensities = _inputProcessor.getMicroVolumes( );

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

            auto microDensity = microDensities->find( *microID );

            if ( microDensity == microDensities->end( ) ){

                return new errorNode( "assembleCouplingMassAndDampingMatrices",
                                      "Free micro node " + std::to_string( *microID ) + " was not found in the micro density map" );

            }

            auto microVolume = microVolumes->find( *microID );

            if ( microVolume == microVolumes->end( ) ){

                return new errorNode( "assembleCouplingMassAndDampingMatrices",
                                      "Free micro node " + std::to_string( *microID ) + " was not found in the micro volume map" );

            }


            freeMicroMasses[ localMicroNodeIDMap->second ] = microVolume->second * microDensity->second;

        }
        //Assemble the ghost micro mass vector
        for ( auto microID = ghostMicroNodeIDs->begin( ); microID != ghostMicroNodeIDs->end( ); microID++ ){

            auto localMicroNodeIDMap = microGlobalToLocalDOFMap->find( *microID );

            if ( localMicroNodeIDMap == microGlobalToLocalDOFMap->end( ) ){

                return new errorNode( "assembleMacroMassAndDampingMatrices",
                                      "Ghost micro node: " + std::to_string( *microID ) + " not found in global to local map\n" );

            }

            auto microDensity = microDensities->find( *microID );

            if ( microDensity == microDensities->end( ) ){

                return new errorNode( "assembleCouplingMassAndDampingMatrices",
                                      "Free micro node " + std::to_string( *microID ) + " was not found in the micro density map" );

            }

            auto microVolume = microVolumes->find( *microID );

            if ( microVolume == microVolumes->end( ) ){

                return new errorNode( "assembleCouplingMassAndDampingMatrices",
                                      "Free micro node " + std::to_string( *microID ) + " was not found in the micro volume map" );

            }

            ghostMicroMasses[ localMicroNodeIDMap->second - nFreeMicroNodes ] = microVolume->second * microDensity->second;

        }

        //Assemble the mass sub-matrices
        tripletVector c1;
        tripletVector c2;

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

        SparseMatrix MQhat( _dim * ghostMicroMasses.size( ), _dim * ghostMicroMasses.size( ) );
        MQhat.setFromTriplets( c1.begin( ), c1.end( ) );
        SparseMatrix MTildeDBreve = rhat * homogenizedMassMatrix + freeMicromorphicMassMatrix;
        //Note: kinetic partitioning coefficient applied when the matrix was formed
    
        SparseMatrix MD    = MTildeDBreve.block( 0, 0, nMacroDispDOF * nFreeMacroNodes, nMacroDispDOF * nFreeMacroNodes );
        SparseMatrix MDhat = MTildeDBreve.block( nMacroDispDOF * nFreeMacroNodes, nMacroDispDOF * nFreeMacroNodes,
                                                 nMacroDispDOF * nGhostMacroNodes, nMacroDispDOF * nGhostMacroNodes );

        //Due to the restrictions in listed in the comment at the beginning of the function, MBar from Regueiro 2012
        //is an empty matrix as we only handle the coupling domain here.
       
    
        //Assemble Mass matrices for the micro projection equation

        if ( config[ "projection_type" ].as< std::string >( ).compare( "l2_projection" ) == 0 ){

            //TODO: Improve efficiency

            auto MQQ  = MQ + _L2_BQhatQ.transpose( ) * MQhat * _L2_BQhatQ + _L2_BDhatQ.transpose( ) * MDhat * _L2_BDhatQ;

            auto MQD = _L2_BQhatQ.transpose( ) * MQhat * _L2_BQhatD + _L2_BDhatQ.transpose( ) * MDhat * _L2_BDhatD;
    
            //Assemble Mass matrices for the macro projection equation
            
            auto MDQ = _L2_BQhatD.transpose( ) * MQhat * _L2_BQhatQ + _L2_BDhatD.transpose( ) * MDhat * _L2_BDhatQ; //TODO: Verify error in Regueiro 2012 for first term second projection matrix

            auto MDD = MD + _L2_BQhatD.transpose( ) * MQhat * _L2_BQhatD + _L2_BDhatD.transpose( ) * MDhat * _L2_BDhatD;
    
            //Assemble the damping matrices for the micro projection equation
            auto CQQ = aQ * MQ + aQ * _L2_BQhatQ.transpose( ) * MQhat * _L2_BQhatQ + aD * _L2_BDhatQ.transpose( ) * MDhat * _L2_BDhatQ;

            auto CQD = aQ * _L2_BQhatQ.transpose( ) * MQhat * _L2_BQhatD;
    
            //Assemble the damping matrices for the macro projection equation
            auto CDQ = aQ * _L2_BQhatD.transpose( ) * MQhat * _L2_BQhatQ;
            auto CDD = aD * MD + aQ * _L2_BQhatD.transpose( ) * MQhat * _L2_BQhatD;

            //Assemble the full mass matrix
            _L2_MASS = Eigen::MatrixXd( MQQ.rows( ) + MDQ.rows( ), MQQ.cols( ) + MQD.cols( ) );
            _L2_MASS.topLeftCorner(     MQQ.rows( ), MQQ.cols( ) ) = MQQ;
            _L2_MASS.topRightCorner(    MQD.rows( ), MQD.cols( ) ) = MQD;
            _L2_MASS.bottomLeftCorner(  MDQ.rows( ), MDQ.cols( ) ) = MDQ;
            _L2_MASS.bottomRightCorner( MDD.rows( ), MDD.cols( ) ) = MDD;

            //Assemble the full damping matrix
            _L2_DAMPING = Eigen::MatrixXd( CQQ.rows( ) + CDQ.rows( ), CQQ.cols( ) + CQD.cols( ) );
            _L2_DAMPING.topLeftCorner(     CQQ.rows( ), CQQ.cols( ) ) = CQQ;
            _L2_DAMPING.topRightCorner(    CQD.rows( ), CQD.cols( ) ) = CQD;
            _L2_DAMPING.bottomLeftCorner(  CDQ.rows( ), CDQ.cols( ) ) = CDQ;
            _L2_DAMPING.bottomRightCorner( CDD.rows( ), CDD.cols( ) ) = CDD;

        }
        else if ( config[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){

            std::cout << "ASSEMBLING MASS BLOCK MATRICES\n";
            SparseMatrix MQQ  =  MQ;
            MQQ += _DP_BQhatQ.transpose( ) * MQhat * _DP_BQhatQ;
            MQQ += _DP_BDhatQ.transpose( ) * MDhat * _DP_BDhatQ;

            SparseMatrix MQD = _DP_BQhatQ.transpose( ) * MQhat * _DP_BQhatD;
            MQD += _DP_BDhatQ.transpose( ) * MDhat * _DP_BDhatD;
    
            //Assemble Mass matrices for the macro projection equation
            
            SparseMatrix MDQ = _DP_BQhatD.transpose( ) * MQhat * _DP_BQhatQ;
            MDQ += _DP_BDhatD.transpose( ) * MDhat * _DP_BDhatQ; //TODO: Verify error in Reguiero 2012 for first term second projection matrix

            SparseMatrix MDD = MD;
            MDD += _DP_BQhatD.transpose( ) * MQhat * _DP_BQhatD;
            MDD += _DP_BDhatD.transpose( ) * MDhat * _DP_BDhatD;

            std::cout << "ASSEMBLING DAMPING BLOCK MATRICES\n";
    
            //Assemble the damping matrices for the micro projection equation
            SparseMatrix CQQ = aQ * MQ;
            CQQ += aQ * _DP_BQhatQ.transpose( ) * MQhat * _DP_BQhatQ;
            CQQ += aD * _DP_BDhatQ.transpose( ) * MDhat * _DP_BDhatQ;

            SparseMatrix CQD = aQ * _DP_BQhatQ.transpose( ) * MQhat * _DP_BQhatD;
    
            //Assemble the damping matrices for the macro projection equation
            SparseMatrix CDQ = aQ * _DP_BQhatD.transpose( ) * MQhat * _DP_BQhatQ;
            SparseMatrix CDD = aD * MD;
            CDD += aQ * _DP_BQhatD.transpose( ) * MQhat * _DP_BQhatD;

            //Assemble the full mass and damping matrices
            std::cout << "ASSEMBLING FULL MASS AND DAMPING MATRICES\n";
            _DP_MASS    = SparseMatrix( MQQ.rows( ) + MDQ.rows( ), MQQ.cols( ) + MQD.cols( ) );
            _DP_DAMPING = SparseMatrix( CQQ.rows( ) + CDQ.rows( ), CQQ.cols( ) + CQD.cols( ) );

            //The sparse matrices are in column major format so we loop over the columns to assemble the matrix
            _DP_MASS.reserve( MQQ.nonZeros( ) + MQD.nonZeros( ) + MDQ.nonZeros( ) + MDD.nonZeros( ) );
            _DP_DAMPING.reserve( CQQ.nonZeros( ) + CQD.nonZeros( ) + CDQ.nonZeros( ) + CDD.nonZeros( ) );
            for ( uIntType c = 0; c < MQQ.cols( ); ++c ){

                //Add terms to the mass matrix
                _DP_MASS.startVec( c );
                for ( SparseMatrix::InnerIterator itMQQ( MQQ, c ); itMQQ; ++itMQQ )
                    _DP_MASS.insertBack( itMQQ.row( ), c ) = itMQQ.value( );
                for ( SparseMatrix::InnerIterator itMDQ( MDQ, c ); itMDQ; ++itMDQ )
                    _DP_MASS.insertBack( itMDQ.row( ) + MQQ.rows( ), c ) = itMDQ.value( );

                //Add terms to the damping matrix
                _DP_DAMPING.startVec( c );
                for ( SparseMatrix::InnerIterator itCQQ( CQQ, c ); itCQQ; ++itCQQ )
                    _DP_DAMPING.insertBack( itCQQ.row( ), c ) = itCQQ.value( );
                for ( SparseMatrix::InnerIterator itCDQ( CDQ, c ); itCDQ; ++itCDQ )
                    _DP_DAMPING.insertBack( itCDQ.row( ) + CQQ.rows( ), c ) = itCDQ.value( );

            }

            for ( uIntType c = 0; c < MDD.cols( ); ++c ){

                //Add terms to the mass matrix
                _DP_MASS.startVec( c + MQQ.cols( ) );
                for ( SparseMatrix::InnerIterator itMQD( MQD, c ); itMQD; ++itMQD )
                    _DP_MASS.insertBack( itMQD.row( ), c + MQQ.cols( ) ) = itMQD.value( );
                for ( SparseMatrix::InnerIterator itMDD( MDD, c ); itMDD; ++itMDD )
                    _DP_MASS.insertBack( itMDD.row( ) + MQD.rows( ), c + MDQ.cols( ) ) = itMDD.value( );

                //Add terms to the damping matrix
                _DP_DAMPING.startVec( c + CQQ.cols( ) );
                for ( SparseMatrix::InnerIterator itCQD( CQD, c ); itCQD; ++itCQD )
                    _DP_DAMPING.insertBack( itCQD.row( ), c + CQQ.cols( ) ) = itCQD.value( );
                for ( SparseMatrix::InnerIterator itCDD( CDD, c ); itCDD; ++itCDD )
                    _DP_DAMPING.insertBack( itCDD.row( ) + CQD.rows( ), c + CDQ.cols( ) ) = itCDD.value( );

            }

        }
        else{

            return new errorNode( "assembleMacroMassAndDampingMatrices",
                                  "The projection type " + config[ "projection_type" ].as< std::string >( ) +
                                  " is not recognized" );

        }

        return NULL;

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

                microDomainVolume += microVolume->second;

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

    errorOut overlapCoupling::assembleCouplingForceVector( ){
        /*!
         * Assemble the force vector. In Regueiro 2012, the external forces are on the right hand side
         * and the internal forces are on the left. I will call the, "force vector," the external forces
         * minus the internal forces. This is distinct from the LHS vector which will include the inertial
         * and damping forces from the previous increment.
         *
         * I also note that, unless we re-compute the micromorphic internal forces, it is not clear how
         * to get the element-wise contributions for the assembly of the vector. Consequently, we will do
         * an approximate scheme to handle the potential energy weighting factor.
         *
         * Note: If the signs in either the macro ( or micro ) scales are flipped, this can be accounted
         * for by using the keywords in the coupling initialization section of the configuration file. The
         * value in "( )" is the default value.
         *
         *     macro_internal_force_sign: ( -1 )
         *     macro_external_force_sign: ( 1 )
         *     micro_internal_force_sign: ( 1 )
         *     micro_external_force_sign: ( 1 )
         */

        const YAML::Node config = _inputProcessor.getCouplingInitialization( );
        floatType qhat = config[ "potential_energy_weighting_factor" ].as< floatType >( );
        std::string projection_type = config[ "projection_type" ].as< std::string >( );

        const uIntType nMacroNodeForces = _dim + _dim * _dim;

        //Get the maps from the global to the local degrees of freedom
        const DOFMap *microGlobalToLocalDOFMap = _inputProcessor.getMicroGlobalToLocalDOFMap( );
        const DOFMap *macroGlobalToLocalDOFMap = _inputProcessor.getMacroGlobalToLocalDOFMap( );

        //Get the free and ghost node ids
        const uIntVector *freeMacroNodeIds = _inputProcessor.getFreeMacroNodeIds( );
        const uIntVector *ghostMacroNodeIds = _inputProcessor.getGhostMacroNodeIds( );

        const uIntVector *freeMicroNodeIds = _inputProcessor.getFreeMicroNodeIds( );
        const uIntVector *ghostMicroNodeIds = _inputProcessor.getGhostMicroNodeIds( );

        const uIntType nFreeMicroNodes = freeMicroNodeIds->size( );
        const uIntType nFreeMacroNodes = freeMacroNodeIds->size( );

        //Extract the internal and external forces
        const std::unordered_map< uIntType, floatVector > *microInternalForces = _inputProcessor.getMicroInternalForces( );
        const std::unordered_map< uIntType, floatVector > *macroInternalForces = _inputProcessor.getMacroInternalForces( );

        const std::unordered_map< uIntType, floatVector > *microExternalForces = _inputProcessor.getMicroExternalForces( );
        const std::unordered_map< uIntType, floatVector > *macroExternalForces = _inputProcessor.getMacroExternalForces( );

        //Assemble the micro internal and external forces
        floatVector FintQhat( _dim * ghostMicroNodeIds->size( ), 0 );
        floatVector FextQhat( _dim * ghostMicroNodeIds->size( ), 0 );
        for ( auto microID  = ghostMicroNodeIds->begin( );
                   microID != ghostMicroNodeIds->end( );
                   microID++ ){

            //Get the local ID
            auto idMap = microGlobalToLocalDOFMap->find( *microID );

            if ( idMap == microGlobalToLocalDOFMap->end( ) ){

                return new errorNode( "assembleCouplingForceVector",
                                      "ghost micro node " + std::to_string( *microID ) +
                                      " is not found in the local to global DOF map" );

            }

            if ( idMap->second < nFreeMicroNodes ){

                return new errorNode( "assembleCouplingForceVector", "The local index is smaller than the number of free micro nodes" );

            }

            if ( ( _dim * ( idMap->second - nFreeMicroNodes ) + _dim ) > FintQhat.size( ) ){

                return new errorNode( "assembleCouplingForceVector", "Local index is larger than the force vector" );

            }

            auto internalForce = microInternalForces->find( *microID );

            if ( internalForce == microInternalForces->end( ) ){

                return new errorNode( "assembleCouplingForceVector", "Micro node " + std::to_string( *microID ) + " not found in internal force vector" );

            }

            auto externalForce = microExternalForces->find( *microID );

            if ( externalForce == microInternalForces->end( ) ){

                return new errorNode( "assembleCouplingForceVector", "Micro node " + std::to_string( *microID ) + " not found in external force vector" );

            }


            for ( unsigned int i = 0; i < _dim; i++ ){

                if ( _inputProcessor.microInternalForceDefined( ) ){

                    FintQhat[ _dim * ( idMap->second - nFreeMicroNodes ) + i ] = ( 1 - qhat ) * internalForce->second[ i ];

                }

                if ( _inputProcessor.microExternalForceDefined( ) ){

                    FextQhat[ _dim * ( idMap->second - nFreeMicroNodes ) + i ] = ( 1 - qhat ) * externalForce->second[ i ];

                }

            }

        }

        floatVector FintQ( _dim * freeMicroNodeIds->size( ), 0 );
        floatVector FextQ( _dim * freeMicroNodeIds->size( ), 0 );
        for ( auto microID  = freeMicroNodeIds->begin( );
                   microID != freeMicroNodeIds->end( );
                   microID++ ){

            //Get the local ID
            auto idMap = microGlobalToLocalDOFMap->find( *microID );

            if ( idMap == microGlobalToLocalDOFMap->end( ) ){

                return new errorNode( "assembleCouplingForceVector",
                                      "free micro node " + std::to_string( *microID ) +
                                      " is not found in the local to global DOF map" );

            }

            auto internalForce = microInternalForces->find( *microID );

            if ( internalForce == microInternalForces->end( ) ){

                return new errorNode( "assembleCouplingForceVector", "Micro node " + std::to_string( *microID ) + " not found in internal force vector" );

            }

            auto externalForce = microExternalForces->find( *microID );

            if ( externalForce == microInternalForces->end( ) ){

                return new errorNode( "assembleCouplingForceVector", "Micro node " + std::to_string( *microID ) + " not found in external force vector" );

            }

            for ( unsigned int i = 0; i < _dim; i++ ){

                if ( _inputProcessor.microInternalForceDefined( ) ){

                    FintQ[ _dim * idMap->second + i ] = ( 1 - qhat ) * internalForce->second[ i ];

                }

                if ( _inputProcessor.microExternalForceDefined( ) ){

                    FextQ[ _dim * idMap->second + i ] = ( 1 - qhat ) * externalForce->second[ i ];

                }

            }

        }

        //Compute the potential energy partitioning coefficients
        std::unordered_map< uIntType, floatType > qes;
        errorOut error = constructPotentialEnergyPartitioningCoefficient( qes );

        if ( error ){

            errorOut result = new errorNode( "assembleCouplingForceVector",
                                             "Error in the construction of the potential energy partitioning coefficients" );
            result->addNext( error );
            return result;

        }

        //Assemble the macro internal and external forces
        floatVector FintDhat( nMacroNodeForces * ghostMacroNodeIds->size( ), 0 );
        floatVector FextDhat( nMacroNodeForces * ghostMacroNodeIds->size( ), 0 );
        for ( auto nodeID = ghostMacroNodeIds->begin( ); nodeID != ghostMacroNodeIds->end( ); nodeID++ ){

            //Get the local ID
            auto idMap = macroGlobalToLocalDOFMap->find( *nodeID );

            if ( idMap == macroGlobalToLocalDOFMap->end( ) ){

                return new errorNode( "assembleCouplingForceVector",
                                      "ghost macro node " + std::to_string( *nodeID ) +
                                      " is not found in the local to global DOF map" );

            }

            auto nodeQes = qes.find( *nodeID );

            if ( idMap->second < nFreeMacroNodes ){

                return new errorNode( "assembleCouplingForceVector",
                                      "ghost macro node " + std::to_string( *nodeID ) +
                                      " has a local position not consistent with the number of free macro nodes" );

            }

            if ( ( nMacroNodeForces * ( idMap->second - nFreeMacroNodes ) + nMacroNodeForces ) > FintDhat.size( ) ){

                return new errorNode( "assembleCouplingForceVector",
                                      "ghost macro node " + std::to_string( *nodeID ) +
                                      " has a local position larger than allocated in the coupling force vector" );

            }

            auto internalForce = macroInternalForces->find( *nodeID );

            if ( internalForce == macroInternalForces->end( ) ){

                return new errorNode( "assembleCouplingForceVector", "Macro node " + std::to_string( *nodeID ) +
                                      " not found in internal force vector" );

            }

            auto externalForce = macroExternalForces->find( *nodeID );

            if ( externalForce == macroInternalForces->end( ) ){

                return new errorNode( "assembleCouplingForceVector", "Macro node " + std::to_string( *nodeID ) +
                                      " not found in external force vector" );

            }

            for ( unsigned int i = 0; i < nMacroNodeForces; i++ ){

                FintDhat[ nMacroNodeForces * ( idMap->second - nFreeMacroNodes ) + i ]
                    = qhat * homogenizedFINT( nMacroNodeForces * idMap->second + i );
                FextDhat[ nMacroNodeForces * ( idMap->second - nFreeMacroNodes ) + i ]
                    = qhat * homogenizedFEXT( nMacroNodeForces * idMap->second + i );

                if ( nodeQes != qes.end( ) ){

                    if ( _inputProcessor.macroInternalForceDefined( ) ){

                        FintDhat[ nMacroNodeForces * ( idMap->second - nFreeMacroNodes ) + i ]
                            += nodeQes->second * internalForce->second[ i ];

                    }

                    if ( _inputProcessor.macroExternalForceDefined( ) ){

                        FintDhat[ nMacroNodeForces * ( idMap->second - nFreeMacroNodes ) + i ]
                            += nodeQes->second * externalForce->second[ i ];

                    }

                }

            }

        }

        floatVector FintD( nMacroNodeForces * freeMacroNodeIds->size( ), 0 );
        floatVector FextD( nMacroNodeForces * freeMacroNodeIds->size( ), 0 );
        for ( auto nodeID = freeMacroNodeIds->begin( ); nodeID != freeMacroNodeIds->end( ); nodeID++ ){

            //Get the local ID
            auto idMap = macroGlobalToLocalDOFMap->find( *nodeID );

            if ( idMap == macroGlobalToLocalDOFMap->end( ) ){

                return new errorNode( "assembleCouplingForceVector",
                                      "free macro node " + std::to_string( *nodeID ) +
                                      " is not found in the local to global DOF map" );

            }

            auto nodeQes = qes.find( *nodeID );

            auto internalForce = macroInternalForces->find( *nodeID );

            if ( internalForce == macroInternalForces->end( ) ){

                return new errorNode( "assembleCouplingForceVector", "Macro node " + std::to_string( *nodeID ) +
                                      " not found in internal force vector" );

            }

            auto externalForce = macroExternalForces->find( *nodeID );

            if ( externalForce == macroInternalForces->end( ) ){

                return new errorNode( "assembleCouplingForceVector", "Macro node " + std::to_string( *nodeID ) +
                                      " not found in external force vector" );

            }

            for ( unsigned int i = 0; i < nMacroNodeForces; i++ ){

                FintD[ nMacroNodeForces * idMap->second + i ]
                    = qhat * homogenizedFINT( nMacroNodeForces * idMap->second + i );
                FextD[ nMacroNodeForces * idMap->second + i ]
                    = qhat * homogenizedFEXT( nMacroNodeForces * idMap->second + i );

                if ( nodeQes != qes.end( ) ){

                    if ( _inputProcessor.macroInternalForceDefined( ) ){

                        FintD[ nMacroNodeForces * idMap->second + i ]
                            += nodeQes->second * internalForce->second[ i ];

                    }

                    if ( _inputProcessor.macroExternalForceDefined( ) ){

                        FextD[ nMacroNodeForces * idMap->second + i ]
                            += nodeQes->second * externalForce->second[ i ];

                    }

                }

            }

        }

        //Map the forces to Eigen Vectors to allow for easy computation
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > _FintQ( FintQ.data(), FintQ.size( ), 1 );
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > _FextQ( FextQ.data(), FextQ.size( ), 1 );
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > _FintD( FintD.data(), FintD.size( ), 1 );
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > _FextD( FextD.data(), FextD.size( ), 1 );

        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > _FintQhat( FintQhat.data(), FintQhat.size( ), 1 );
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > _FextQhat( FextQhat.data(), FextQhat.size( ), 1 );
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > _FintDhat( FintDhat.data(), FintDhat.size( ), 1 );
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > _FextDhat( FextDhat.data(), FextDhat.size( ), 1 );

        //Initialize the force vectors
        Eigen::MatrixXd _FQ, _FD;

        

        if ( projection_type.compare( "l2_projection" ) == 0 ){

            //Assemble the micro force vector
            _FQ  = _FextQ;
            _FQ += _L2_BQhatQ.transpose( ) * _FextQhat;
            _FQ -= _FintQ;
            _FQ -= _L2_BQhatQ.transpose( ) * _FintQhat;
            _FQ -= _L2_BDhatQ.transpose( ) * _FintDhat;
    
            //Assemble the macro force vector
            _FD  = _FextD;
            _FD += _L2_BDhatD.transpose( ) * _FextDhat;
            _FD -= _FintD;
            _FD -= _L2_BQhatD.transpose( ) * _FintQhat;
            _FD -= _L2_BDhatD.transpose( ) * _FintDhat;

        }
        else if ( projection_type.compare( "direct_projection" ) == 0 ){

            //Assemble the micro force vector
            _FQ  = _FextQ;
            _FQ += _DP_BQhatQ.transpose( ) * _FextQhat;
            _FQ -= _FintQ;
            _FQ -= _DP_BQhatQ.transpose( ) * _FintQhat;
            _FQ -= _DP_BDhatQ.transpose( ) * _FintDhat;
    
            //Assemble the macro force vector
            _FD  = _FextD;
            _FD += _DP_BDhatD.transpose( ) * _FextDhat;
            _FD -= _FintD;
            _FD -= _DP_BQhatD.transpose( ) * _FintQhat;
            _FD -= _DP_BDhatD.transpose( ) * _FintDhat;

        }
        else{

            return new errorNode( "assembleCouplingForceVector",
                                  "The projection type: " + projection_type + " is not recognized" );

        }

        //Assemble the full force vector
        _FORCE = Eigen::MatrixXd( _FQ.rows( ) + _FD.rows( ), 1 );
        _FORCE << _FQ, _FD;

//        std::cout << "FintQ:\n";
//        for ( unsigned int i = 0; i < FintQ.size( ); i++ ){
//            std::cout << _FintQ( i, 0 ) << " ";
//            if ( ( ( i + 1 ) % 12 ) == 0 ){
//                std::cout << "\n";
//            }
//        }
//        std::cout << "FintD:\n";
//        for ( unsigned int i = 0; i < FintD.size( ); i++ ){
//            std::cout << _FintD( i, 0 ) << " ";
//            if ( ( ( i + 1 ) % 12 ) == 0 ){
//                std::cout << "\n";
//            }
//        }

//        std::cout << "FextQ:\n";
//        for ( unsigned int i = 0; i < FextQ.size( ); i++ ){
//            std::cout << _FextQ( i, 0 ) << " ";
//            if ( ( ( i + 1 ) % 12 ) == 0 ){
//                std::cout << "\n";
//            }
//        }
//        std::cout << "FextD:\n";
//        for ( unsigned int i = 0; i < FextD.size( ); i++ ){
//            std::cout << _FextD( i, 0 ) << " ";
//            if ( ( ( i + 1 ) % 12 ) == 0 ){
//                std::cout << "\n";
//            }
//        }

//        std::cout << "FextQhat:\n";
//        for ( unsigned int i = 0; i < FextQhat.size( ); i++ ){
//            std::cout << _FextQhat( i, 0 ) << " ";
//            if ( ( ( i + 1 ) % 12 ) == 0 ){
//                std::cout << "\n";
//            }
//        }
//        std::cout << "FextDhat:\n";
//        for ( unsigned int i = 0; i < FextDhat.size( ); i++ ){
//            std::cout << _FextDhat( i, 0 ) << " ";
//            if ( ( ( i + 1 ) % 12 ) == 0 ){
//                std::cout << "\n";
//            }
//        }

//        std::cout << "FintQhat:\n";
//        for ( unsigned int i = 0; i < FintQhat.size( ); i++ ){
//            std::cout << _FintQhat( i, 0 ) << " ";
//            if ( ( ( i + 1 ) % 12 ) == 0 ){
//                std::cout << "\n";
//            }
//        }
//        std::cout << "FintDhat:\n";
//        for ( unsigned int i = 0; i < FintDhat.size( ); i++ ){
//            std::cout << _FintDhat( i, 0 ) << " ";
//            if ( ( ( i + 1 ) % 12 ) == 0 ){
//                std::cout << "\n";
//            }
//        }

        std::cout << "FQ:\n";
//        for ( unsigned int i = 0; i < _FQ.size( ); i++ ){
//            std::cout << _FQ( i, 0 ) << " ";
//            if ( ( ( i + 1 ) % 12 ) == 0 ){
//                std::cout << "\n";
//            }
//        }
//        std::cout << "FD:\n";
//        for ( unsigned int i = 0; i < _FD.size( ); i++ ){
//            std::cout << _FD( i, 0 ) << " ";
//            if ( ( ( i + 1 ) % 12 ) == 0 ){
//                std::cout << "\n";
//            }
//        }
//
//        floatVector tmpVec( 3, 0 );
//        for ( unsigned int i = 0; i < _FQ.size( ); i++ ){
//            tmpVec[ ( i + 1 ) % 3 ] += fabs( _FQ( i, 0 ) );
//        }
//        std::cout << "Sum of FQ: "; vectorTools::print( tmpVec );
//
//        tmpVec = floatVector( 12, 0 );
//        for ( unsigned int i = 0; i < _FD.size( ); i++ ){
//            tmpVec[ ( i + 1 ) % 12 ] += fabs( _FD( i, 0 ) );
//        }
//        std::cout << "Sum of FD: "; vectorTools::print( tmpVec );
//
//        return new errorNode( "assembleCouplingForceVector", "derp" );
//
//        std::cout << "FextQ NORM:    " << _FextQ.norm( ) << "\n";
//        std::cout << "FextD NORM:    " << _FextD.norm( ) << "\n";
//        std::cout << "FextQhat NORM: " << _FextQhat.norm( ) << "\n";
//        std::cout << "FextDhat NORM: " << _FextDhat.norm( ) << "\n";
//        std::cout << "FintQ NORM:    " << _FintQ.norm( ) << "\n";
//        std::cout << "FintD NORM:    " << _FintD.norm( ) << "\n";
//        std::cout << "FintQhat NORM: " << _FintQhat.norm( ) << "\n";
//        std::cout << "FintDhat NORM: " << _FintDhat.norm( ) << "\n";
//        std::cout << "FQ NORM:       " << _FQ.norm( ) << "\n";
//        std::cout << "FD NORM:       " << _FD.norm( ) << "\n";
//        std::cout << "FORCE NORM:    " << _FORCE.norm( ) << "\n";

        return NULL;
    }

    errorOut overlapCoupling::constructPotentialEnergyPartitioningCoefficient( std::unordered_map< uIntType, floatType > &qes ){
        /*!
         * Construct the potential energy partitioning coefficient ( qe in Regueiro 2012 )
         *
         * Currently only volume_fraction is implemented
         *
         * :param std::unordered_map< uIntType, floatType > &qes: The potential energy partitioning coefficents
         *     at the free nodes in the free macro cells.
         */

        const YAML::Node config = _inputProcessor.getCouplingInitialization( );
        std::string strategy = config[ "potential_energy_partitioning_coefficient" ][ "type" ].as< std::string >( );

        if ( strategy.compare( "volume_fraction" ) == 0 ){

            const uIntVector *freeMacroNodeIds = _inputProcessor.getFreeMacroNodeIds( );
            const uIntVector *ghostMacroNodeIds = _inputProcessor.getGhostMacroNodeIds( );
            const uIntVector *freeMacroCellIds = _inputProcessor.getFreeMacroCellIds( );
    
            //Compute the volumes of the overlaped macro domains
            std::unordered_map< uIntType, floatType > macroCellNodeTotalVolumes;
            macroCellNodeTotalVolumes.reserve( freeMacroNodeIds->size( ) + ghostMacroNodeIds->size( ) );
    
            qes.reserve( freeMacroNodeIds->size( ) );
    
            for ( auto macroCellID  = freeMacroCellIds->begin( );
                       macroCellID != freeMacroCellIds->end( );
                       macroCellID++ ){
    
                //Construct the macro-domain element
                std::unique_ptr< elib::Element > element;
                errorOut error = buildMacroDomainElement( *macroCellID,
                                                          *_inputProcessor.getMacroNodeReferencePositions( ),
                                                          *_inputProcessor.getMacroDisplacements( ),
                                                          *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                          element );
    
                if ( error ){
    
                    errorOut result = new errorNode( "assembleCouplingForceVector",
                                                     "Error in the construction of the macro element " + std::to_string( *macroCellID ) );
                    result->addNext( error );
                    return result;
    
                }
    
                //Compute the volume of the macro element
                floatType elementVolume = 0;
                floatMatrix jacobian;
                
                for ( auto qpt = element->qrule.begin( ); qpt != element->qrule.end( ); qpt++ ){
    
                    //Get the Jacobian of transformation
                    errorOut error = element->get_local_gradient( element->nodes, qpt->first, jacobian );
    
                    if ( error ){
     
                        errorOut result = new errorNode( "assembleCouplingForceVector",
                                                         "Error in the computation of the local gradient\n" );
                        result->addNext( error );
                        return result;
    
                    }
    
                    elementVolume += vectorTools::determinant( vectorTools::appendVectors( jacobian ), _dim, _dim ) * qpt->second;
    
                }
    
                //Compute the volumes of the overlapped micro domains for the macro cell
                auto microDomainVolumes = homogenizedVolumes.find( *macroCellID );
                if ( microDomainVolumes == homogenizedVolumes.end( ) ){
    
                    return new errorNode( "assembleCouplingForceVector",
                                          "The macro cell " + std::to_string( *macroCellID ) +
                                          " is not found in the homogenized volumes map" );
    
                }
    
                floatType microDomainVolume = 0;
    
                for ( auto microVolume  = microDomainVolumes->second.begin( );
                           microVolume != microDomainVolumes->second.end( );
                           microVolume++ ){
    
                    microDomainVolume += microVolume->second;
    
                }
    
                //Determine the amount of the macro element which is non-overlapped
                floatType openMacroVolume;
                if ( elementVolume < _absoluteTolerance ){
                    openMacroVolume = 0;
                }
                else{
                    openMacroVolume = std::fmax( elementVolume - microDomainVolume, 0. );
                }
    
                //Store the contributions of the open macro volume and element volume to the element's free nodes
                for ( auto nodeID  = element->global_node_ids.begin( );
                           nodeID != element->global_node_ids.end( );
                           nodeID++ ){
    
                    //Store the contribution to the nodes
                    if ( qes.find( *nodeID ) != qes.end( ) ){
    
                        qes[ *nodeID ] += openMacroVolume;
                        macroCellNodeTotalVolumes[ *nodeID ] += elementVolume;
    
                    }
                    else{

                        qes.emplace( *nodeID, openMacroVolume );
                        macroCellNodeTotalVolumes[ *nodeID ] += elementVolume;

                    }
    
                }
    
            }
    
            for ( auto nodeID  = macroCellNodeTotalVolumes.begin( );
                       nodeID != macroCellNodeTotalVolumes.end( );
                       nodeID++ ){
    
                qes[ nodeID->first ] /= nodeID->second;
    
            }

        }
        else{

            return new errorNode( "assembleCouplingForceVector",
                                  "The potential energy partitioning strategy: " + strategy + " is not recognized." );

        }

        return NULL;
    }

    errorOut overlapCoupling::solveFreeDisplacement( const bool updateGhostDOF ){
        /*!
         * Solve for the free displacement. This is done using the Newmark-Beta method ( from Wikipedia ):
         *
         * u_{ n + 1 } = u_n + \Delta t * \dot{ u }_n + \frac{ \Delta t^2 }{ 2 }\left( \left( 1 - 2 \beta \right) \ddot{ u }_n + 2 \beta \ddot{ u }_{n+1 } \right)
         * \dot{ u }_{ n + 1 } = \dot{ u }_n + ( 1 - \gamma ) \Delta t \ddot{ u }_n + \gamma \Delta t \ddot{ u }_{n + 1}
         * M \ddot{ u }_{n + 1} + C \dot{ u }_{n + 1} + f^{int}\left( u_{n + 1} \right) = f^{ext}_{n+1}
         *
         * Where the Explicit central difference scheme is obtained by letting $\gamma = 0.5$ and $\beta = 0$
         *
         * and
         *
         * Average constant acceleration is obtained by setting $\gamma = 0.5$ and $\beta = 0.25$
         *
         * NOTE: The current implementation cannot do implicit solves here but it does solve for what the
         * displacement values *should* be given the current values of the force vectors and those can be
         * used to construct residual equations. It is hoped that Jacobian Free Newton Krylov will be 
         * capable of solving the coupled PDEs without having to explicitly form the Jacobian.
         *
         * :param const bool updateGhostDOF: A flag for whether the ghost degrees of freedom should be updated as well
         */

        //Set the configuration
        const YAML::Node config = _inputProcessor.getCouplingInitialization( );

        std::string projection_type = config[ "projection_type" ].as< std::string >( );
        const floatType gamma = *_inputProcessor.getNewmarkGamma( );
        const floatType beta = *_inputProcessor.getNewmarkBeta( );

        //Get the timestep
        const floatType *dt = _inputProcessor.getDt( );

        //Set the number of displacement degrees of freedom for each scale
        const uIntType nMicroDispDOF = _dim;
        const uIntType nMacroDispDOF = _dim + _dim * _dim;

        //Get the maps from the global to the local degrees of freedom
        const DOFMap *microGlobalToLocalDOFMap = _inputProcessor.getMicroGlobalToLocalDOFMap( );
        const DOFMap *macroGlobalToLocalDOFMap = _inputProcessor.getMacroGlobalToLocalDOFMap( );

        //Get the free degree of freedom values at both scale
        const uIntVector *freeMicroNodeIds = _inputProcessor.getFreeMicroNodeIds( );
        const uIntVector *freeMacroNodeIds = _inputProcessor.getFreeMacroNodeIds( );

        const uIntType microOffset = nMicroDispDOF * freeMicroNodeIds->size( );

        //Get the previous values of the degrees of freedom and their velocities and accelerations
        const std::unordered_map< uIntType, floatVector > *previousMicroDispDOFVector = _inputProcessor.getPreviousMicroDisplacements( );
        const std::unordered_map< uIntType, floatVector > *previousMicroVelocities    = _inputProcessor.getPreviousMicroVelocities( );
        const std::unordered_map< uIntType, floatVector > *previousMicroAccelerations = _inputProcessor.getPreviousMicroAccelerations( );

        const std::unordered_map< uIntType, floatVector > *previousMacroDispDOFVector = _inputProcessor.getPreviousMacroDispDOFVector( );
        const std::unordered_map< uIntType, floatVector > *previousMacroVelocities    = _inputProcessor.getPreviousMacroVelocities( );
        const std::unordered_map< uIntType, floatVector > *previousMacroAccelerations = _inputProcessor.getPreviousMacroAccelerations( );

        //Set the degree of freedom vectors
        //To try and save some memory, we will overwrite DOF and DotDOF as the integrated values become available
        floatVector       FreeDOF( nMicroDispDOF * freeMicroNodeIds->size( ) + nMacroDispDOF * freeMacroNodeIds->size( ), 0 );
        floatVector        DotDOF( nMicroDispDOF * freeMicroNodeIds->size( ) + nMacroDispDOF * freeMacroNodeIds->size( ), 0 );
        floatVector   DotDotDOF_t( nMicroDispDOF * freeMicroNodeIds->size( ) + nMacroDispDOF * freeMacroNodeIds->size( ), 0 );
        floatVector DotDotDOF_tp1( nMicroDispDOF * freeMicroNodeIds->size( ) + nMacroDispDOF * freeMacroNodeIds->size( ), 0 );

        //Extract the micro points
        for ( auto nodeId = freeMicroNodeIds->begin( ); nodeId != freeMicroNodeIds->end( ); nodeId++ ){

            //Find the map from global to local node ids
            auto indexMap = microGlobalToLocalDOFMap->find( *nodeId );

            if ( indexMap == microGlobalToLocalDOFMap->end( ) ){

                return new errorNode( "solveFreeDisplacement",
                                      "Micro node " + std::to_string( *nodeId ) + " not found in global to local map" );

            }

            auto previousMicroDisp = previousMicroDispDOFVector->find( *nodeId );

            if ( previousMicroDisp == previousMicroDispDOFVector->end( ) ){

                return new errorNode( "solveFreeDisplacement",
                                      "The micro node " + std::to_string( *nodeId ) +
                                      " is not found in the previous micro displacement dof vector" );

            }

            auto previousMicroVel = previousMicroVelocities->find( *nodeId );

            if ( previousMicroVel == previousMicroVelocities->end( ) ){

                return new errorNode( "solveFreeDisplacement",
                                      "The micro node " + std::to_string( *nodeId ) +
                                      " is not found in the previous micro velocities vector" );

            }

            auto previousMicroAccel = previousMicroAccelerations->find( *nodeId );

            if ( previousMicroAccel == previousMicroAccelerations->end( ) ){

                return new errorNode( "solveFreeDisplacement",
                                      "The micro node " + std::to_string( *nodeId ) +
                                      " is not found in the previous micro accelerations vector" );

            }

            //Extract the micro points
            for ( unsigned int i = 0; i < nMicroDispDOF; i++ ){

                FreeDOF[ nMicroDispDOF * indexMap->second + i ] = previousMicroDisp->second[ i ];

                if ( _inputProcessor.microVelocitiesDefined( ) ){

                    DotDOF[ nMicroDispDOF * indexMap->second + i ] = previousMicroVel->second[ i ];

                }

                if ( _inputProcessor.microAccelerationDefined( ) ){

                    DotDotDOF_t[ nMicroDispDOF * indexMap->second + i ] = previousMicroAccel->second[ i ];

                }

            }

        }

        //Extract macro points
        for ( auto nodeId = freeMacroNodeIds->begin( ); nodeId != freeMacroNodeIds->end( ); nodeId++ ){

            //Find the map from global to local node ids
            auto indexMap = macroGlobalToLocalDOFMap->find( *nodeId );

            if ( indexMap == macroGlobalToLocalDOFMap->end( ) ){

                return new errorNode( "solveFreeDisplacement",
                                      "Macro node " + std::to_string( *nodeId ) + " not found in global to local map" );

            }

            auto previousMacroDisp = previousMacroDispDOFVector->find( *nodeId );

            if ( previousMacroDisp == previousMacroDispDOFVector->end( ) ){

                return new errorNode( "solveFreeDisplacement",
                                      "The macro node " + std::to_string( *nodeId ) +
                                      " is not found in the previous macro displacement dof vector" );

            }

            auto previousMacroVel = previousMacroVelocities->find( *nodeId );

            if ( previousMacroVel == previousMacroVelocities->end( ) ){

                return new errorNode( "solveFreeDisplacement",
                                      "The macro node " + std::to_string( *nodeId ) +
                                      " is not found in the previous macro velocities vector" );

            }

            auto previousMacroAccel = previousMacroAccelerations->find( *nodeId );

            if ( previousMacroAccel == previousMacroAccelerations->end( ) ){

                return new errorNode( "solveFreeDisplacement",
                                      "The macro node " + std::to_string( *nodeId ) +
                                      " is not found in the previous macro accelerations vector" );

            }

            //Extract the macro points
            for ( unsigned int i = 0; i < nMacroDispDOF; i++ ){

                FreeDOF[ nMacroDispDOF * indexMap->second + i + microOffset ] = previousMacroDisp->second[ i ];

                if ( _inputProcessor.macroVelocitiesDefined( ) ){

                    DotDOF[ nMacroDispDOF * indexMap->second + i + microOffset ] = previousMacroVel->second[ i ];

                }

                if ( _inputProcessor.macroAccelerationDefined( ) ){

                    DotDotDOF_t[ nMicroDispDOF * indexMap->second + i + microOffset ] = previousMacroAccel->second[ i ];

                }

            }

        }

        //Map the vectors to Eigen matrices
        Eigen::Map< Eigen::Matrix< floatType, -1,  1 > > _DOF( FreeDOF.data(), FreeDOF.size( ), 1 );
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > _DotDOF( DotDOF.data(), DotDOF.size( ), 1 );
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > _DotDotDOF_t( DotDotDOF_t.data(), DotDotDOF_t.size( ), 1 );
        Eigen::Map< Eigen::Matrix< floatType, -1,  1 > > _DotDotDOF_tp1( DotDotDOF_tp1.data(), DotDotDOF_tp1.size( ), 1 );

        //Instantiate the QR solver
        std::cout << "Performing QR decomposition of the Free DOF LHS matrix\n";

        Eigen::MatrixXd RHS;
        if ( projection_type.compare( "l2_projection" ) == 0 ){

            Eigen::MatrixXd LHS;
            LHS = _L2_MASS;
            LHS += gamma * ( *dt ) * _L2_DAMPING;

            RHS = _FORCE;
            RHS -= _L2_DAMPING * ( _DotDOF + ( 1 - gamma ) * ( *dt ) * _DotDotDOF_t );

            _DotDotDOF_tp1 = LHS.colPivHouseholderQr( ).solve( RHS );
        }
        else if ( projection_type.compare( "direct_projection" ) == 0 ){

            SparseMatrix LHS( _DP_MASS.rows( ), _DP_MASS.cols( ) );
            LHS = _DP_MASS;
            LHS += gamma * ( *dt ) * _DP_DAMPING;
            LHS.makeCompressed( );

            RHS = _FORCE;
            RHS -= _DP_DAMPING * ( _DotDOF + ( 1 - gamma ) * ( *dt ) * _DotDotDOF_t );

            Eigen::SparseQR< SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
            solver.compute( LHS );
            _DotDotDOF_tp1 = solver.solve( RHS );

//            return new errorNode( "solveFreeDisplacement", "derp3" );
//
//            std::cout << "_DotDof:\n";
//            std::cout << "microscale\n";
//            for ( unsigned int i = 0; i < microOffset; i++ ){
//                std::cout << _DotDOF( i, 0 ) << " ";
//                if ( ( ( i + 1 ) % 12 ) == 0 ){
//                    std::cout << "\n";
//                }
//            }
//            std::cout << "macroscale\n";
//            for ( unsigned int i = microOffset; i < _DotDOF.size( ); i++ ){
//                std::cout << _DotDOF( i, 0 ) << " ";
//                if ( ( ( i + 1 ) % 12 ) == 0 ){
//                    std::cout << "\n";
//                }
//            }
//
//            return new errorNode( "solveFreeDisplacement", "derp" );
//
//            std::cout << "_DotDotDof_t:\n";
//            std::cout << "microscale\n";
//            for ( unsigned int i = 0; i < microOffset; i++ ){
//                std::cout << _DotDotDOF_t( i, 0 ) << " ";
//                if ( ( ( i + 1 ) % 12 ) == 0 ){
//                    std::cout << "\n";
//                }
//            }
//            std::cout << "macroscale\n";
//            for ( unsigned int i = microOffset; i < _DotDotDOF_t.size( ); i++ ){
//                std::cout << _DotDotDOF_t( i, 0 ) << " ";
//                if ( ( ( i + 1 ) % 12 ) == 0 ){
//                    std::cout << "\n";
//                }
//            }
//
//            std::cout << "_DotDotDof_tp1:\n";
//            std::cout << "microscale\n";
//            for ( unsigned int i = 0; i < microOffset; i++ ){
//                std::cout << _DotDotDOF_tp1( i, 0 ) << " ";
//                if ( ( ( i + 1 ) % 12 ) == 0 ){
//                    std::cout << "\n";
//                }
//            }
//            std::cout << "macroscale\n";
//            for ( unsigned int i = microOffset; i < _DotDotDOF_tp1.size( ); i++ ){
//                std::cout << _DotDotDOF_tp1( i, 0 ) << " ";
//                if ( ( ( i + 1 ) % 12 ) == 0 ){
//                    std::cout << "\n";
//                }
//            }
//
//            return new errorNode( "solveFreeDisplacement", "derp2" );

        }
        else{

            return new errorNode( "solveFreeDisplacement", "Projection type " + projection_type + " not recognized" );

        }

        //Update the free degrees of freedom
        _DOF += ( *dt ) * _DotDOF + 0.5 * ( ( *dt ) * ( *dt ) ) * ( ( 1 - 2 * beta ) * _DotDotDOF_t + 2 * beta * _DotDotDOF_tp1 );

        std::cout << "_DOF:\n";
        std::cout << "microscale\n";
        for ( unsigned int i = 0; i < microOffset; i++ ){
            std::cout << _DOF( i, 0 ) << " ";
            if ( ( ( i + 1 ) % 12 ) == 0 ){
                std::cout << "\n";
            }
        }
        std::cout << "macroscale\n";
        for ( unsigned int i = microOffset; i < _DOF.size( ); i++ ){
            std::cout << _DOF( i, 0 ) << " ";
            if ( ( ( i + 1 ) % 12 ) == 0 ){
                std::cout << "\n";
            }
        }

        return new errorNode( "solveFreeDisplacement", "derp4" );

        //Store the free degrees of freedom
        _updatedFreeMicroDispDOFValues = floatVector( FreeDOF.begin( ), FreeDOF.begin( ) + microOffset );
        _updatedFreeMacroDispDOFValues = floatVector( FreeDOF.begin( ) + microOffset, FreeDOF.end( ) );
        _freeDOFValuesUpdated = true;

        if ( updateGhostDOF ){

            errorOut error = projectDegreesOfFreedom( updateGhostDOF );

            if ( error ){
    
                errorOut result = new errorNode( "solveFreeDisplacement",
                                                 "Error in the projection of the ghost degrees of freedom" );
                result->addNext( error );
                return result;
    
            }            

        }

        return NULL;

    }

    errorOut overlapCoupling::outputReferenceInformation( ){
        /*!
         * Output the reference information to file
         *
         * For now, I'm using a XDMF file which may, or may not, be the most 
         * appropriate file type to use but it does allow me to store the 
         * data in a binary file ( rather than ASCII ). I want to store the
         * data in hdf5 ( which XDMF also does ) so the end result of doing
         * this versus a custom approach is probably not a big deal.
         */

        //Get the coupling initialization
        YAML::Node couplingInitialization = _inputProcessor.getCouplingInitialization( );

        //Initialize the XDMF domain
        shared_ptr< XdmfDomain > domain = XdmfDomain::New( );
        std::string domainInfoDescription = "This is not a mesh-based XDMF file and should only be used ";
        domainInfoDescription += "/ read by overlapCoupling and associated file readers";

        shared_ptr< XdmfInformation > domainInformation = XdmfInformation::New( "REFERENCE_INFORMATION_DOMAIN", domainInfoDescription );
        domain->insert( domainInformation );
        shared_ptr< XdmfUnstructuredGrid > grid = XdmfUnstructuredGrid::New( );
        domain->insert( grid );

        //Get the output filename
        std::string reference_filename
            = couplingInitialization[ "output_reference_information" ][ "filename" ].as< std::string >( );

        //Try to remove the files
        std::string xdmfFilename = reference_filename + ".xdmf";
        std::string h5Filename = reference_filename + ".h5";
        remove( xdmfFilename.c_str( ) );
        remove( h5Filename.c_str( ) );

        //Save the interpolation matrix if required
        if ( couplingInitialization[ "output_reference_information" ][ "save_interpolation_matrix" ] ){

            errorOut error = writeSparseMatrixToXDMF( _N, "N", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( "outputReferenceInformation",
                                                 "Error when writing out the interpolation matrix N" );
                result->addNext( error );
                return result;

            }


        }

        errorOut error;

        //Save the projection matrices
        error = writeSparseMatrixToXDMF( _centerOfMassN, "centerOfMassInterpolator", reference_filename, domain, grid );

        if ( error ){

            errorOut result = new errorNode( "outputReferenceInformation",
                                             "Error when writing out the center of mass projection matrix" );
            result->addNext( error );
            return result;

        }

        error = writeDenseMatrixToXDMF( _centerOfMassProjector, "centerOfMassProjector", reference_filename, domain, grid );

        if ( error ){

            errorOut result = new errorNode( "outputReferenceInformation",
                                             "Error when writing out the center of mass projection matrix" );
            result->addNext( error );
            return result;

        }

        std::string projectionType = couplingInitialization[ "projection_type" ].as< std::string >( );
        if ( ( projectionType.compare( "l2_projection" ) == 0 ) || ( projectionType.compare( "averaged_l2_projection" ) == 0 ) ){

            //Initialize the projection matrix attributes
            shared_ptr< XdmfAttribute > BQhatQ = XdmfAttribute::New( );
            shared_ptr< XdmfAttribute > BQhatD = XdmfAttribute::New( );
            shared_ptr< XdmfAttribute > BDhatQ = XdmfAttribute::New( );
            shared_ptr< XdmfAttribute > BDhatD = XdmfAttribute::New( );

            //Write the matrix type to the output file
            shared_ptr< XdmfInformation > projectionType = XdmfInformation::New( "EIGEN_MATRIX_TYPE", "DENSE" );
            grid->insert( projectionType );

            //Write BQhatQ
            error = writeDenseMatrixToXDMF( _L2_BQhatQ, "BQhatQ", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( "outputReferenceInformation",
                                                 "Error when writing out BQhatQ" );
                result->addNext( error );
                return result;

            }

            //Write BQhatD
            error = writeDenseMatrixToXDMF( _L2_BQhatD, "BQhatD", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( "outputReferenceInformation",
                                                 "Error when writing out BQhatD" );
                result->addNext( error );
                return result;

            }

            //Write BDhatQ
            error = writeDenseMatrixToXDMF( _L2_BDhatQ, "BDhatQ", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( "outputReferenceInformation",
                                                 "Error when writing out BDhatQ" );
                result->addNext( error );
                return result;

            }

            //Write BDhatD
            error = writeDenseMatrixToXDMF( _L2_BDhatD, "BDhatD", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( "outputReferenceInformation",
                                                 "Error when writing out BDhatD" );
                result->addNext( error );
                return result;

            }

        }
        else if ( projectionType.compare( "direct_projection" ) == 0 ){

            //Write the matrix type to the output file
            shared_ptr< XdmfInformation > projectionType = XdmfInformation::New( "EIGEN_MATRIX_TYPE", "SPARSE" );
            domain->insert( projectionType );

            //Write BQhatQ
            error = writeSparseMatrixToXDMF( _DP_BQhatQ, "BQhatQ", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( "outputReferenceInformation",
                                                 "Error when writing out BQhatQ" );
                result->addNext( error );
                return result;

            }

            //Write BQhatD
            error = writeSparseMatrixToXDMF( _DP_BQhatD, "BQhatD", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( "outputReferenceInformation",
                                                 "Error when writing out BQhatD" );
                result->addNext( error );
                return result;

            }

            //Write BDhatQ
            error = writeSparseMatrixToXDMF( _DP_BDhatQ, "BDhatQ", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( "outputReferenceInformation",
                                                 "Error when writing out BDhatQ" );
                result->addNext( error );
                return result;

            }

            //Write BDhatD
            error = writeSparseMatrixToXDMF( _DP_BDhatD, "BDhatD", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( "outputReferenceInformation",
                                                 "Error when writing out BDhatD" );
                result->addNext( error );
                return result;

            }

        }
        else{

            return new errorNode( "outputReferenceInformation",
                                  "The projection type " + projectionType + " is not recognized" );

        }

        return NULL;

    }

    errorOut writeDenseMatrixToXDMF( const Eigen::MatrixXd &A, const std::string matrixName,
                                     const std::string &filename, shared_ptr< XdmfDomain > &domain,
                                     shared_ptr< XdmfUnstructuredGrid > &grid ){
        /*!
         * Write a dense matrix to a XDMF file ( not a standard XDMF format )
         *
         * :param const Eigen::MatrixXd &A: The matrix to write
         * :param const std::string matrixName: The name of the matrix MUST BE UNIQUE IN THE OUTPUT FILE
         * :param const std::string &filename: The name of the output file
         * :param shared_ptr< XdmfDomain > &domain: The XDMF domain
         * :param shared_ptr< XdmfDomain > &grid: The XDMF grid
         */

        //Initialize the writer
        shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( filename + ".h5", false );
        heavyWriter->setReleaseData( true );
        shared_ptr< XdmfWriter > writer = XdmfWriter::New( filename + ".xdmf", heavyWriter );

        shared_ptr< XdmfAttribute > _A = XdmfAttribute::New( );
        _A->setName( matrixName );

        shared_ptr< XdmfInformation > A_info
            = XdmfInformation::New( matrixName + "_shape",
                                    std::to_string( A.rows( ) ) + "," + std::to_string( A.cols( ) ) );
        _A->insert( A_info );
        _A->insert( 0, A.data( ), A.size( ), 1, 1 );
        grid->insert( _A );

        domain->accept( writer );

        return NULL;

    }

    errorOut writeSparseMatrixToXDMF( const SparseMatrix &A, const std::string matrixName,
                                      const std::string &filename, shared_ptr< XdmfDomain > &domain,
                                      shared_ptr< XdmfUnstructuredGrid > &grid ){
        /*!
         * Write a sparse matrix to a XDMF file ( not a standard XDMF format )
         *
         * :param const SparseMatrix &A: The matrix to write
         * :param const std::string matrixName: The name of the matrix MUST BE UNIQUE IN THE OUTPUT FILE
         * :param const std::string &filename: The name of the output file
         * :param shared_ptr< XdmfDomain > &domain: The XDMF domain
         * :param shared_ptr< XdmfDomain > &grid: The XDMF grid
         */

        //Initialize the writer
        shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( filename + ".h5", false );
        heavyWriter->setReleaseData( true );
        shared_ptr< XdmfWriter > writer = XdmfWriter::New( filename + ".xdmf", heavyWriter );

        shared_ptr< XdmfAttribute > _A = XdmfAttribute::New( );
        _A->setName( matrixName );
        shared_ptr< XdmfInformation > _A_info
            = XdmfInformation::New( matrixName + "_shape",
                                    std::to_string( A.rows( ) ) + "," + std::to_string( A.cols( ) ) );
        _A->insert( _A_info );
        grid->insert( _A );

        uIntVector rowIndices, colIndices;
        uIntType nonZeros = A.nonZeros( );
        rowIndices.reserve( nonZeros );
        colIndices.reserve( nonZeros );

        floatVector values;
        values.reserve( nonZeros );

        for (int k=0; k < A.outerSize(); ++k) {

            for (SparseMatrix::InnerIterator it(A,k); it; ++it){

                rowIndices.push_back( it.row( ) );
                colIndices.push_back( it.col( ) );
                values.push_back( it.value( ) );

            }

        }

        shared_ptr< XdmfAttribute > _A_rows = XdmfAttribute::New( );
        _A_rows->setName( matrixName + "_rows" );
        _A_rows->insert( 0, rowIndices.data( ), rowIndices.size( ), 1, 1 );
        grid->insert( _A_rows );

        shared_ptr< XdmfAttribute > _A_cols = XdmfAttribute::New( );
        _A_cols->setName( matrixName + "_cols" );
        _A_cols->insert( 0, colIndices.data( ), colIndices.size( ), 1, 1 );
        grid->insert( _A_cols );

        shared_ptr< XdmfAttribute > _A_values = XdmfAttribute::New( );
        _A_values->setName( matrixName + "_values" );
        _A_values->insert( 0, values.data( ), values.size( ), 1, 1 );
        grid->insert( _A_values );

        domain->accept( writer );
        return NULL;
    }

    errorOut readDenseMatrixFromXDMF( const shared_ptr< XdmfUnstructuredGrid > &grid, const std::string &matrixName, Eigen::MatrixXd &A ){
        /*!
         * Read a dense matrix from a XDMF file attribute. This is not a standard use of a XDMF file.
         *
         * NOTE: Assumes default Eigen storage ( column major );
         *
         * :param const shared_ptr< XdmfUnstructuredGrid > &grid: The pointer to the grid object that contains the matrix
         * :param const std::string &matrixName: The name of the matrix
         * :param Eigen::MatrixXdx &A: The re-constructed matrix
         */

        shared_ptr< XdmfAttribute > _A = grid->getAttribute( matrixName );

        if ( !_A ){

            std::string outstr = matrixName;
            outstr += " does not appear as an attribute in the provided grid\n";
            return new errorNode( "readDenseMatrixFromXDMF", outstr );

        }

        shared_ptr< XdmfInformation > shapeInfo = _A->getInformation( 0 );

        if ( !shapeInfo ){

            std::string outstr = "There is no information defined for the matrix " + matrixName;
            return new errorNode( "readDenseMatrixFromXDMF", outstr );

        }

        if ( shapeInfo->getKey( ).compare( matrixName + "_shape" ) != 0 ){

            std::string outstr = matrixName;
            outstr += "_shape is not in the information key";
            return new errorNode( "readDenseMatrixFromXDMF", outstr );

        }

        std::string matrixDimensionString = shapeInfo->getValue( );

        //Extract the rows and columns
        uIntType rows = std::stoi( matrixDimensionString.substr( 0, matrixDimensionString.find( "," ) ) );
        uIntType cols = std::stoi( matrixDimensionString.substr( matrixDimensionString.find( "," ) + 1, matrixDimensionString.size( ) ) );

        _A->read( );

        //Initialize the matrix
        A = Eigen::MatrixXd( rows, cols );

        for ( unsigned int i = 0; i < cols; i++ ){
            for ( unsigned int j = 0; j < rows; j++ ){
                A( j, i ) = _A->getValue< floatType >( rows * i + j );
            }
        }

        return NULL;
    }

    errorOut readSparseMatrixFromXDMF( const shared_ptr< XdmfUnstructuredGrid > &grid, const std::string &matrixName, SparseMatrix &A ){
        /*!
         * Read a sparse matrix from a XDMF file attribute. This is not a standard use of a XDMF file.
         *
         * :param const shared_ptr< XdmfUnstructuredGrid > &grid: The pointer to the grid object that contains the matrix
         * :param const std::string &matrixName: The name of the matrix
         * :param SparseMatrix &A: The re-constructed sparse matrix
         */

        shared_ptr< XdmfAttribute > _A = grid->getAttribute( matrixName );

        if ( !_A ){

            std::string outstr = matrixName;
            outstr += " does not appear as an attribute in the provided grid\n";
            return new errorNode( "readSparseMatrixFromXDMF", outstr );

        }

        shared_ptr< XdmfInformation > shapeInfo = _A->getInformation( 0 );

        if ( !shapeInfo ){

            std::string outstr = "There is no information defined for the SparseMatrix " + matrixName;
            return new errorNode( "readSparseMatrixFromXDMF", outstr );

        }

        if ( shapeInfo->getKey( ).compare( matrixName + "_shape" ) != 0 ){

            std::string outstr = matrixName;
            outstr += "_shape is not in the information key";
            return new errorNode( "readSparseMatrixFromXDMF", outstr );

        }

        std::string matrixDimensionString = shapeInfo->getValue( );

        //Extract the rows and cols
        uIntType rows = std::stoi( matrixDimensionString.substr( 0, matrixDimensionString.find( "," ) ) );
        uIntType cols = std::stoi( matrixDimensionString.substr( matrixDimensionString.find( "," ) + 1, matrixDimensionString.size( ) ) );

        //Construct the triplet vector
        shared_ptr< XdmfAttribute > _rows = grid->getAttribute( matrixName + "_rows" );

        if ( !_rows ){

            std::string outstr = matrixName;
            outstr += "_rows attribute is not found";
            return new errorNode( "readSparseMatrixFromXDMF", outstr );

        }

        shared_ptr< XdmfAttribute > _cols = grid->getAttribute( matrixName + "_cols" );

        if ( !_cols ){

            std::string outstr = matrixName;
            outstr += "_cols attribute is not found";
            return new errorNode( "readSparseMatrixFromXDMF", outstr );

        }

        shared_ptr< XdmfAttribute > _vals = grid->getAttribute( matrixName + "_values" );

        if ( !_cols ){

            std::string outstr = matrixName;
            outstr += "_values attribute is not found";
            return new errorNode( "readSparseMatrixFromXDMF", outstr );

        }

        if ( ( _rows->getSize( ) != _cols->getSize( ) ) && ( _rows->getSize( ) != _vals->getSize( ) ) ){

            std::string outstr = matrixName;
            outstr += " attributes rows, cols, and values don't have consistent sizes";
            return new errorNode( "readSparseMatrixFromXDMF", outstr );

        }

        _rows->read( );
        _cols->read( );
        _vals->read( );

        tripletVector triplets;
        triplets.reserve( _rows->getSize( ) );

        for ( uIntType i = 0; i != _rows->getSize( ); i++ ){

            uIntType _r = _rows->getValue< uIntType >( i );
            uIntType _c = _cols->getValue< uIntType >( i );
            floatType _v = _vals->getValue< floatType >( i );

            triplets.push_back( DOFProjection::T( _r, _c, _v ) );

        }

        A = SparseMatrix( rows, cols );
        A.setFromTriplets( triplets.begin( ), triplets.end( ) );

        return NULL;

    }

    errorOut overlapCoupling::extractProjectionMatricesFromFile( ){
        /*!
         * Extract the projection matrices from the storage file
         */

        //Set the filename
        YAML::Node config = _inputProcessor.getCouplingInitialization( );
        std::string filename = config[ "reference_filename" ].as< std::string >( );

        //Initialize the XDMF reader
        shared_ptr< XdmfReader > reader = XdmfReader::New( );
        shared_ptr< XdmfDomain > _readDomain = shared_dynamic_cast< XdmfDomain >( reader->read( filename ) );
        shared_ptr< XdmfUnstructuredGrid > _readGrid = _readDomain->getUnstructuredGrid( 0 );

        std::string projectionType = config[ "type" ].as< std::string >( );

        errorOut error = readSparseMatrixFromXDMF( _readGrid, "centerOfMassInterpolator", _centerOfMassN );

        if ( error ){

            errorOut result = new errorNode( "extractProjectionMatricesFromFile",
                                             "Error in extracting the center of mass interpolation" );
            result->addNext( error );
            return result;

        }

        error = readDenseMatrixFromXDMF( _readGrid, "centerOfMassProjector", _centerOfMassProjector );

        if ( error ){

            errorOut result = new errorNode( "extractProjectionMatricesFromFile",
                                             "Error in extracting the center of mass projector" );
            result->addNext( error );
            return result;

        }

        if ( ( projectionType.compare( "l2_projection" ) ) || ( projectionType.compare( "averaged_l2_projection" ) ) ){

            error = readDenseMatrixFromXDMF( _readGrid, "BQhatQ", _L2_BQhatQ );

            if ( error ){

                errorOut result
                    = new errorNode( "extractProjectionMatricesFromFile", "Error when extracting BQhatQ from the XDMF file" );
                result->addNext( error );
                return result;

            }

            error = readDenseMatrixFromXDMF( _readGrid, "BQhatD", _L2_BQhatD );

            if ( error ){

                errorOut result
                    = new errorNode( "extractProjectionMatricesFromFile", "Error when extracting BQhatD from the XDMF file" );
                result->addNext( error );
                return result;

            }

            error = readDenseMatrixFromXDMF( _readGrid, "BDhatQ", _L2_BDhatQ );

            if ( error ){

                errorOut result
                    = new errorNode( "extractProjectionMatricesFromFile", "Error when extracting BDhatQ from the XDMF file" );
                result->addNext( error );
                return result;

            }

            error = readDenseMatrixFromXDMF( _readGrid, "BDhatD", _L2_BDhatD );

            if ( error ){

                errorOut result
                    = new errorNode( "extractProjectionMatricesFromFile", "Error when extracting BDhatD from the XDMF file" );
                result->addNext( error );
                return result;

            }

        }
        else if ( projectionType.compare( "direct_projection" ) ){

            error = readSparseMatrixFromXDMF( _readGrid, "BQhatQ", _DP_BQhatQ );

            if ( error ){

                errorOut result
                    = new errorNode( "extractProjectionMatricesFromFile", "Error when extracting BQhatQ from the XDMF file" );
                result->addNext( error );
                return result;

            }

            error = readSparseMatrixFromXDMF( _readGrid, "BQhatD", _DP_BQhatD );

            if ( error ){

                errorOut result
                    = new errorNode( "extractProjectionMatricesFromFile", "Error when extracting BQhatD from the XDMF file" );
                result->addNext( error );
                return result;

            }

            error = readSparseMatrixFromXDMF( _readGrid, "BDhatQ", _DP_BDhatQ );

            if ( error ){

                errorOut result
                    = new errorNode( "extractProjectionMatricesFromFile", "Error when extracting BDhatQ from the XDMF file" );
                result->addNext( error );
                return result;

            }

            error = readSparseMatrixFromXDMF( _readGrid, "BDhatD", _DP_BDhatD );

            if ( error ){

                errorOut result
                    = new errorNode( "extractProjectionMatricesFromFile", "Error when extracting BDhatD from the XDMF file" );
                result->addNext( error );
                return result;

            }

        }
        else{

            return new errorNode( "extractProjectionMatricesFromFile", "Not implemented" );

        }

        return NULL;

    }

    errorOut overlapCoupling::outputHomogenizedResponse( const uIntType collectionNumber ){
        /*!
         * Output the homogenized response to the data file
         *
         * :param const uIntType collectionNumber: The collection to place the reference state in ( defaults to 0 )
         */

        //Get the configuration
        const YAML::Node config = _inputProcessor.getCouplingInitialization( );

        //Form the writer
        std::shared_ptr< dataFileInterface::dataFileBase > writer
            = dataFileInterface::dataFileBase( config[ "output_homogenized_response" ] ).create( );

        if ( writer->_error ){

            errorOut result = new errorNode( "outputHomogenizedResponse",
                                             "Error when initializing the writer" );
            result->addNext( writer->_error );
            return result;

        }

        const floatType *time = _inputProcessor.getMacroTime( );

        uIntType increment;
        errorOut error = writer->initializeIncrement( *time, _currentReferenceOutputIncrement, collectionNumber, increment );

        if ( error ){

            errorOut result = new errorNode( "outputHomogenizedResponse",
                                             "Error in the initialization of the increment for the homogenized output" );
            result->addNext( error );
            return result;

        }

        //Write the mesh data to file. This increment references a previous increment for the mesh data so we don't need to do it again
        error = writer->writeIncrementMeshData( increment, collectionNumber, { }, { { } }, { }, { }, { }, { { } }, { }, { } );

        if ( error ){

            errorOut result = new errorNode( "outputHomogenizedResponse",
                                             "Error in the initialization of the mesh data for the homogenized output" );
            result->addNext( error );
            return result;

        }

        //Assemble the macro-scale cell ids
        const uIntVector cellIds = vectorTools::appendVectors( { *_inputProcessor.getFreeMacroCellIds( ),
                                                                 *_inputProcessor.getGhostMacroCellIds( ) } );

        //Determine the maximum number of quadrature points
        uIntType maxQP = 0;
        for ( auto cellId = cellIds.begin( ); cellId != cellIds.end( ); cellId++ ){

            maxQP = std::max( ( uIntType )quadraturePointDensities[ *cellId ].size( ), maxQP );

        }

        //Assemble the projected displacement vectors
        floatVector projectedMacroDisplacement;
        if ( !_currentReferenceOutputIncrement ){

            projectedMacroDisplacement
                = vectorTools::appendVectors(
                    {
                        floatVector( ( _dim + _dim * _dim ) * _inputProcessor.getFreeMacroNodeIds( )->size( ), 0 ),
                        _projected_ghost_macro_displacement
                    } );

        }
        else{

            projectedMacroDisplacement =
                vectorTools::appendVectors( { _updatedFreeMacroDispDOFValues, _projected_ghost_macro_displacement } );

        }

        //Loop over the quadrature points
        for ( unsigned int qp = 0; qp < maxQP; qp++ ){

            floatVector densityOut(                                      cellIds.size( ), 0 );
            floatVector bodyForceOut(                             _dim * cellIds.size( ), 0 );
            floatVector accelerationsOut(                         _dim * cellIds.size( ), 0 );
            floatVector microInertiasOut(                  _dim * _dim * cellIds.size( ), 0 );
            floatVector bodyCouplesOut(                    _dim * _dim * cellIds.size( ), 0 );
            floatVector microSpinInertiasOut(              _dim * _dim * cellIds.size( ), 0 );
            floatVector symmetricMicroStressOut(           _dim * _dim * cellIds.size( ), 0 );
            floatVector cauchyStressOut(                   _dim * _dim * cellIds.size( ), 0 );
            floatVector higherOrderStressOut(       _dim * _dim * _dim * cellIds.size( ), 0 );
            floatVector dofValuesOut(           ( _dim + _dim * _dim ) * cellIds.size( ), 0 );
            floatVector dofGradientsOut( _dim * ( _dim + _dim * _dim ) * cellIds.size( ), 0 );

            //Loop over the cells
            for ( auto cellId = cellIds.begin( ); cellId != cellIds.end( ); cellId++ ){

                //Set the cell index
                uIntType index = cellId - cellIds.begin( );

                //Determine the number of quadrature points for the cell
                uIntType cellQP = quadraturePointDensities[ *cellId ].size( );

                if ( qp >= cellQP ){

                    //If the current quadrature point is too big then skip the element
                    continue;

                }

                //Set the outputs for the current quadrature point
                densityOut[ index ] = quadraturePointDensities[ *cellId ][ qp ];

                for ( unsigned int i = 0; i < _dim; i++ ){

                    bodyForceOut[ _dim * index + i ] 
                        = quadraturePointBodyForce[     *cellId ][        _dim * qp + i ];
                    accelerationsOut[ _dim * index + i ] 
                        = quadraturePointAccelerations[ *cellId ][        _dim * qp + i ];

                }

                for ( unsigned int i = 0; i < _dim * _dim; i++ ){

                    microInertiasOut[ _dim * _dim * index + i ]
                        = quadraturePointMicroInertias[ *cellId ][ _dim * _dim * qp + i ];
                    bodyCouplesOut[ _dim * _dim * index + i ]
                        = quadraturePointBodyCouples[        *cellId ][ _dim * _dim * qp + i ];
                    microSpinInertiasOut[ _dim * _dim * index + i ]
                        = quadraturePointMicroSpinInertias[ *cellId ][ _dim * _dim * qp + i ];
                    symmetricMicroStressOut[ _dim * _dim * index + i ]
                        = quadraturePointSymmetricMicroStress[ *cellId ][ _dim * _dim * qp + i ];
                    cauchyStressOut[ _dim * _dim * index + i ]
                        = quadraturePointCauchyStress[ *cellId ][ _dim * _dim * qp + i ];

                }

                for ( unsigned int i = 0; i < _dim * _dim * _dim; i++ ){

                    higherOrderStressOut[ _dim * _dim * _dim * index + i ]
                        = quadraturePointHigherOrderStress[ *cellId ][ _dim * _dim * _dim * qp + i ];

                }

                //Form the element
                std::unique_ptr< elib::Element > element;
                error = buildMacroDomainElement( *cellId,
                                                 *_inputProcessor.getMacroNodeReferencePositions( ),
                                                 *_inputProcessor.getMacroDisplacements( ),
                                                 *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                 element );

                //Build the DOF vector
                floatMatrix dofMatrix( element->qrule.size( ), floatVector( _dim + _dim * _dim, 0 ) );

                for ( auto node = element->global_node_ids.begin( ); node != element->global_node_ids.end( ); node++ ){

                    auto localNode = _inputProcessor.getMacroGlobalToLocalDOFMap( )->find( *node );

                    if ( localNode == _inputProcessor.getMacroGlobalToLocalDOFMap( )->end( ) ){

                        return new errorNode( "outputHomogenizedResponse",
                                              "Error in finding the global node " + std::to_string( *node ) +
                                              " in the macro global to local DOF map" );

                    }

                    dofMatrix[ node - element->global_node_ids.begin( ) ] =
                        floatVector( projectedMacroDisplacement.begin( ) + ( _dim + _dim * _dim ) * localNode->second,
                                     projectedMacroDisplacement.begin( ) + ( _dim + _dim * _dim ) * ( localNode->second + 1 ) );

                }

                error = element->interpolate( dofMatrix, element->qrule[ qp ].first, dofValuesOut );

                if ( error ){

                    errorOut result = new errorNode( "outputHomogenizedResponse",
                                                     "Error in the interpolation of the DOF values" );
                    result->addNext( error );
                    return result;

                }

                floatMatrix qptDOFGradient;

                error = element->get_global_gradient( dofMatrix, element->qrule[ qp ].first, element->reference_nodes, qptDOFGradient );

                if ( error ){

                    errorOut result = new errorNode( "outputHomogenizedResponse",
                                                     "Error in the interpolation of the DOF values" );
                    result->addNext( error );
                    return result;

                }

                dofGradientsOut = vectorTools::appendVectors( qptDOFGradient );

            }

            //Write quadrature point information to file
            stringVector outputNames = { "density_" + std::to_string( qp ) };
            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", densityOut );

            if ( error ){

                errorOut result = new errorNode( "outputHomogenizedResponse", "Error in outputting the density" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim );

            for ( unsigned int i = 0; i < _dim; i++ ){

                outputNames[ i ] = "acceleration_" + std::to_string( i + 1 ) + "_" + std::to_string( qp );

            }

            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", accelerationsOut );

            if ( error ){

                errorOut result = new errorNode( "outputHomogenizedResponse", "Error in outputting the acceleration" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim );
            for ( unsigned int i = 0; i < _dim; i++ ){

                outputNames[ i ] = "body_force_" + std::to_string( i + 1 ) + "_" + std::to_string( qp );

            }
            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", bodyForceOut );

            if ( error ){

                errorOut result = new errorNode( "outputHomogenizedResponse", "Error in outputting the body force" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim * _dim );
            for ( unsigned int i = 0; i < _dim; i++ ){

                for ( unsigned int j = 0; j < _dim; j++ ){

                    outputNames[ _dim * i + j ] = "micro_inertia_" + std::to_string( i + 1 )  + std::to_string( j + 1 ) + "_" + std::to_string( qp );

                }

            }
            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", microInertiasOut );

            if ( error ){

                errorOut result = new errorNode( "outputHomogenizedResponse", "Error in outputting the micro inertias" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim * _dim );
            for ( unsigned int i = 0; i < _dim; i++ ){

                for ( unsigned int j = 0; j < _dim; j++ ){

                    outputNames[ _dim * i + j ] = "body_couple_" + std::to_string( i + 1 )  + std::to_string( j + 1 ) + "_" + std::to_string( qp );

                }

            }
            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", bodyCouplesOut );

            if ( error ){

                errorOut result = new errorNode( "outputHomogenizedResponse", "Error in outputting the body couples" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim * _dim );
            for ( unsigned int i = 0; i < _dim; i++ ){

                for ( unsigned int j = 0; j < _dim; j++ ){

                    outputNames[ _dim * i + j ] = "micro_spin_inertia_" + std::to_string( i + 1 )  + std::to_string( j + 1 ) + "_" + std::to_string( qp );

                }

            }
            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", microSpinInertiasOut );

            if ( error ){

                errorOut result = new errorNode( "outputHomogenizedResponse", "Error in outputting the body couples" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim * _dim );
            for ( unsigned int i = 0; i < _dim; i++ ){

                for ( unsigned int j = 0; j < _dim; j++ ){

                    outputNames[ _dim * i + j ] = "symmetric_micro_stress_" + std::to_string( i + 1 )  + std::to_string( j + 1 ) + "_" + std::to_string( qp );

                }

            }
            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", symmetricMicroStressOut );

            if ( error ){

                errorOut result = new errorNode( "outputHomogenizedResponse", "Error in outputting the symmetric micro stress" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim * _dim );
            for ( unsigned int i = 0; i < _dim; i++ ){

                for ( unsigned int j = 0; j < _dim; j++ ){

                    outputNames[ _dim * i + j ] = "cauchy_stress_" + std::to_string( i + 1 )  + std::to_string( j + 1 ) + "_" + std::to_string( qp );

                }

            }
            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", cauchyStressOut );

            if ( error ){

                errorOut result = new errorNode( "outputHomogenizedResponse", "Error in outputting the cauchy stress" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim * _dim * _dim );
            for ( unsigned int i = 0; i < _dim; i++ ){

                for ( unsigned int j = 0; j < _dim; j++ ){

                    for ( unsigned int k = 0; k < _dim; k++ ){

                        outputNames[ _dim * _dim * i + _dim * j + k ] = "higher_order_stress_" + std::to_string( i + 1 )  + std::to_string( j + 1 ) + std::to_string( k + 1 ) + "_" + std::to_string( qp );

                    }

                }

            }
            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", higherOrderStressOut );

            if ( error ){

                errorOut result = new errorNode( "outputHomogenizedResponse", "Error in outputting the higher order stress" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim + _dim * _dim );
            for ( unsigned int i = 0; i < ( _dim + _dim * _dim ); i++ ){

                    outputNames[ i ] = "DOF_" + std::to_string( i + 1 ) + "_" + std::to_string( qp );

            }
            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", dofValuesOut );

            if ( error ){

                errorOut result = new errorNode( "outputHomogenizedResponse", "Error in outputting the degree of freedom values" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim * ( _dim + _dim * _dim ) );
            for ( unsigned int i = 0; i < ( _dim + _dim * _dim ); i++ ){

                for ( unsigned int j = 0; j < _dim; j++ ){

                    outputNames[ _dim * i + j ]
                        = "DOF_" + std::to_string( i + 1 ) + "," + std::to_string( j + 1 ) + "_" + std::to_string( qp );

                }

            }

            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", dofGradientsOut );

            if ( error ){

                errorOut result = new errorNode( "outputHomogenizedResponse", "Error in outputting the degree of freedom values" );
                result->addNext( error );
                return result;

            }

        }
        
        return NULL;
    }

    errorOut overlapCoupling::writeReferenceMeshDataToFile( const uIntType collectionNumber ){
        /*!
         * Write the reference mesh data to the output file
         *
         * :param const uIntType collectionNumber: The collection to place the reference state in ( defaults to 0 )
         */

        YAML::Node config = _inputProcessor.getCouplingInitialization( )[ "output_homogenized_response" ]; 

        //Form a writer object
        std::shared_ptr< dataFileInterface::dataFileBase > writer
            = dataFileInterface::dataFileBase( config ).create( );

        if ( writer->_error ){

            errorOut result = new errorNode( "writeReferenceMeshDataToFile", "Error in construction of writer" );
            result->addNext( writer->_error );
            return result;

        }

        //Get the free and ghost micro nodes global to local map
        const DOFMap *macroGlobalToLocalDOFMap = _inputProcessor.getMacroGlobalToLocalDOFMap( );

        //Get the macro node ids
        const uIntVector *freeMacroNodeIds = _inputProcessor.getFreeMacroNodeIds( );
        const uIntVector *ghostMacroNodeIds = _inputProcessor.getGhostMacroNodeIds( );

        //Get the macro cell ids
        const uIntVector elementIds = vectorTools::appendVectors( { *_inputProcessor.getFreeMacroCellIds( ),
                                                                    *_inputProcessor.getGhostMacroCellIds( ) } );

        //Get the mesh information
        const std::unordered_map< uIntType, uIntVector > *macroNodeReferenceConnectivity
            = _inputProcessor.getMacroNodeReferenceConnectivity( );

        //Assemble the mesh data
        const std::unordered_map< uIntType, floatVector > *macroNodeReferencePositions
            = _inputProcessor.getMacroNodeReferencePositions( );
        const std::unordered_map< uIntType, floatVector > *macroDisplacements = _inputProcessor.getMacroDisplacements( );

        //Assemble the node ids
        uIntVector nodeIds( macroGlobalToLocalDOFMap->size( ) );

        //Assemble the node position vectors
        floatVector nodePositions( _dim * macroGlobalToLocalDOFMap->size( ) );

        for ( auto node  = macroGlobalToLocalDOFMap->begin( );
                   node != macroGlobalToLocalDOFMap->end( );
                   node++ ){

            nodeIds[ node->second ] = node->first;

            auto referencePosition = macroNodeReferencePositions->find( node->first );

            if ( referencePosition == macroNodeReferencePositions->end( ) ){

                return new errorNode( "writeReferenceMeshDataToFile", "The macro node " + std::to_string( node->first ) +
                                      " was not found in the reference positions vector" );

            }

            auto displacement = macroDisplacements->find( node->first );

            if ( displacement == macroDisplacements->end( ) ){

                return new errorNode( "writeReferenceMeshDataToFile", "The macro node " + std::to_string( node->first ) +
                                      " was not found in the displacements vector" );

            }

            for ( unsigned int i = 0; i < _dim; i++ ){

                nodePositions[ _dim * node->second + i ] = referencePosition->second[ i ] + displacement->second[ i ];

            }

        }

        //Assemble the node sets
        stringVector nodeSetNames = { "free_macro_nodes", "ghost_macro_nodes" };
        uIntMatrix nodeSets( 2 );
        nodeSets[ 0 ].resize( freeMacroNodeIds->size( ) );
        nodeSets[ 1 ].resize( ghostMacroNodeIds->size( ) );

        for ( auto node = freeMacroNodeIds->begin( ); node != freeMacroNodeIds->end( ); node++ ){

            auto localNode = macroGlobalToLocalDOFMap->find( *node );

            if ( localNode == macroGlobalToLocalDOFMap->end( ) ){

                return new errorNode( "writeReferenceMeshDataToFile",
                                      "The free macro node " + std::to_string( *node ) + " was not found in the macro global to local DOF map" );

            }

            nodeSets[ 0 ][ node - freeMacroNodeIds->begin( ) ] = localNode->second;

        }

        for ( auto node = ghostMacroNodeIds->begin( ); node != ghostMacroNodeIds->end( ); node++ ){

            auto localNode = macroGlobalToLocalDOFMap->find( *node );

            if ( localNode == macroGlobalToLocalDOFMap->end( ) ){

                return new errorNode( "writeReferenceMeshDataToFile",
                                      "The ghost macro node " + std::to_string( *node ) + " was not found in the macro global to local DOF map" );

            }

            nodeSets[ 1 ][ node - ghostMacroNodeIds->begin( ) ] = localNode->second;

        }

        //Assemble the element sets
        stringVector elementSetNames = { "free_macro_elements", "ghost_macro_elements" };
        uIntMatrix elementSets( 2 );
        elementSets[ 0 ].resize( _inputProcessor.getFreeMacroCellIds( )->size( ) );
        elementSets[ 1 ].resize( _inputProcessor.getGhostMacroCellIds( )->size( ) );

        for ( uIntType i = 0; i < elementSets[ 0 ].size( ); i++ ){
            elementSets[ 0 ][ i ] = i;
        }

        for ( uIntType i = 0; i < elementSets[ 1 ].size( ); i++ ){
            elementSets[ 1 ][ i ] = elementSets[ 0 ].size( ) + i;
        }

        //Assemble the connectivity vector
        uIntVector connectivity( 0 );
        for ( auto cell = elementIds.begin( ); cell != elementIds.end( ); cell++ ){

            std::unique_ptr< elib::Element > element;

            errorOut error = buildMacroDomainElement( *cell, *macroNodeReferencePositions, *macroDisplacements,
                                                      *macroNodeReferenceConnectivity, element );

            if ( error ){

                errorOut result = new errorNode( "writeReferenceMeshDataToFile",
                                                 "Error in construction of the micromorphic element" );
                result->addNext( error );
                return result;

            }

            //Get the XDMF cell type
            auto elementConnectivity = macroNodeReferenceConnectivity->find( *cell );

            if ( elementConnectivity == macroNodeReferenceConnectivity->end( ) ){

                return new errorNode( "writeReferenceMeshDataToFile",
                                      "Macro cell " + std::to_string( *cell ) + " was not found in the macro mesh connectivity" );

            }

            uIntType cellType = elementConnectivity->second[ 0 ];

            uIntVector localNodeIds( element->global_node_ids.size( ) + 1 );
            localNodeIds[ 0 ] = cellType;

            for ( auto gN = element->global_node_ids.begin( ); gN != element->global_node_ids.end( ); gN++ ){

                auto node = macroGlobalToLocalDOFMap->find( *gN );

                if ( node == macroGlobalToLocalDOFMap->end( ) ){

                    return new errorNode( "writeReferenceMeshDataToFile",
                                          "The global macro node " + std::to_string( *gN ) +
                                          " can't be found in the global to local map" );

                }

                localNodeIds[ gN - element->global_node_ids.begin( ) + 1 ] = node->second;

            }

            connectivity = vectorTools::appendVectors( { connectivity, localNodeIds } );

        }

        const floatType *time = _inputProcessor.getMacroTime( );

        uIntType numIncrements;
        writer->getNumIncrements( numIncrements );

        errorOut error = writer->initializeIncrement( *time, numIncrements, collectionNumber, _currentReferenceOutputIncrement ); 

        if ( error ){

            errorOut result = new errorNode( "writeReferenceMeshDataToFile",
                                             "Error in initialization of the increment in the writer" );
            result->addNext( error );
            return result;

        }

        error = writer->writeIncrementMeshData( _currentReferenceOutputIncrement, collectionNumber,
                                                nodeIds, nodeSets, nodeSetNames, nodePositions,
                                                elementIds, elementSets, elementSetNames, connectivity );

        if ( error ){

            errorOut result = new errorNode( "writeReferenceMeshDataToFile",
                                             "Error in writing the mesh data to a file" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut overlapCoupling::writeUpdatedDOFToFile( const uIntType collectionNumber ){
        /*!
         * Write the updated degree of freedom information to a output file
         *
         * :param const uIntType collectionNumber: The collection to place the reference state in ( defaults to 0 )
         */

        //Get the configuration for the writer object
        YAML::Node config = _inputProcessor.getCouplingInitialization( )[ "output_updated_dof" ];

        if ( !config ){

            return new errorNode( "writeUpdatedDOFToFile",
                                  "'output_updated_dof' is not defined in the configuration file" );

        }

        std::string macro_config_string = "filename: ";
        macro_config_string += config[ "macroscale_filename" ].as< std::string >( ) + "\n";
        macro_config_string += "mode: write\n";
        macro_config_string += "filetype: ";
        macro_config_string += config[ "macroscale_filetype" ].as< std::string >( ) + "\n";
        macro_config_string += "append_to_existing_file: false\n";

        std::string micro_config_string = "filename: ";
        micro_config_string += config[ "microscale_filename" ].as< std::string >( ) + "\n";
        micro_config_string += "mode: write\n";
        micro_config_string += "filetype: ";
        micro_config_string += config[ "microscale_filetype" ].as< std::string >( ) + "\n";
        micro_config_string += "append_to_existing_file: false\n";

        YAML::Node macro_config = YAML::Load( macro_config_string.c_str( ) ); 
        YAML::Node micro_config = YAML::Load( micro_config_string.c_str( ) ); 

        //Form the macro writer object
        std::shared_ptr< dataFileInterface::dataFileBase > macro_writer
            = dataFileInterface::dataFileBase( macro_config ).create( );

        if ( macro_writer->_error ){

            errorOut result = new errorNode( "writeReferenceMeshDataToFile", "Error in construction of writer" );
            result->addNext( macro_writer->_error );
            return result;

        }

        //Form the micro writer object
        std::shared_ptr< dataFileInterface::dataFileBase > micro_writer
            = dataFileInterface::dataFileBase( micro_config ).create( );

        if ( micro_writer->_error ){

            errorOut result = new errorNode( "writeReferenceMeshDataToFile", "Error in construction of writer" );
            result->addNext( micro_writer->_error );
            return result;

        }

        //Write out the macro DOF information
        floatVector outputDOF = vectorTools::appendVectors( { _updatedFreeMacroDispDOFValues, _projected_ghost_macro_displacement } );

//        std::cout << "Macroscale\n";
//        for ( unsigned int i = 0; i < outputDOF.size( ); i++ ){
//            std::cout << outputDOF[ i ] << " ";
//            if ( ( ( i + 1 ) % 12 ) == 0 ){
//                std::cout << "\n";
//            }
//        }

        uIntType increment;
        errorOut error = macro_writer->initializeIncrement( 1.0, 0, collectionNumber, increment );

        if ( error ){

            errorOut result = new errorNode( "writeUpdatedDOFToFile", "Error in initializating the increment of the macro output file" );
            result->addNext( error );
            return result;

        }

        error = macro_writer->writeScalarSolutionData( increment, collectionNumber, "updated_DOF", "Node", outputDOF ); 
        
        if ( error ){

            errorOut result = new errorNode( "writeUpdatedDOFToFile", "Error in outputting the updated macro DOF to the output file" );
            result->addNext( error );
            return result;

        }

        floatVector nodeIds( _inputProcessor.getMacroGlobalToLocalDOFMap( )->size( ) );
        for ( auto node  = _inputProcessor.getMacroGlobalToLocalDOFMap( )->begin( );
                   node != _inputProcessor.getMacroGlobalToLocalDOFMap( )->end( );
                   node++ ){

            nodeIds[ node->second ] = node->first;

        }

        error = macro_writer->writeScalarSolutionData( increment, collectionNumber, "node_ids", "Node", nodeIds );

        if ( error ){

            errorOut result = new errorNode( "writeUpdatedDOFToFile", "Error in outputting the updated node ids to the micro output file" );
            result->addNext( error );
            return result;

        }

        //Write out the micro DOF information
        outputDOF = vectorTools::appendVectors( { _updatedFreeMicroDispDOFValues, _projected_ghost_micro_displacement } );

//        std::cout << "\nMicroscale\n";
//        for ( unsigned int i = 0; i < outputDOF.size( ); i++ ){
//            std::cout << outputDOF[ i ] << " ";
//            if ( ( ( i + 1 ) % 3 ) == 0 ){
//                std::cout << "\n";
//            }
//        }

        error = micro_writer->initializeIncrement( 1.0, 0, collectionNumber, increment );

        if ( error ){

            errorOut result = new errorNode( "writeUpdatedDOFToFile", "Error in initializating the increment of the micro output file" );
            result->addNext( error );
            return result;

        }

        error = micro_writer->writeScalarSolutionData( increment, collectionNumber, "updated_DOF", "Node", outputDOF ); 
        
        if ( error ){

            errorOut result = new errorNode( "writeUpdatedDOFToFile", "Error in outputting the updated micro DOF to the output file" );
            result->addNext( error );
            return result;

        }

        nodeIds = floatVector( _inputProcessor.getMicroGlobalToLocalDOFMap( )->size( ) );
        for ( auto node  = _inputProcessor.getMicroGlobalToLocalDOFMap( )->begin( );
                   node != _inputProcessor.getMicroGlobalToLocalDOFMap( )->end( );
                   node++ ){

            nodeIds[ node->second ] = node->first;

        }

        error = micro_writer->writeScalarSolutionData( increment, collectionNumber, "node_ids", "Node", nodeIds );

        if ( error ){

            errorOut result = new errorNode( "writeUpdatedDOFToFile", "Error in outputting the updated node ids to the micro output file" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    DOFMap overlapCoupling::getMicroGlobalLocalNodeMap( ){
        /*!
         * Return a copy of the micro global to local node map
         */

        return *_inputProcessor.getMicroGlobalToLocalDOFMap( );
    }

    DOFMap overlapCoupling::getMacroGlobalLocalNodeMap( ){
        /*!
         * Return a copy of the macro global to local node map
         */

        return *_inputProcessor.getMacroGlobalToLocalDOFMap( );
    }

    floatVector overlapCoupling::getUpdatedMicroDisplacementDOF( ){
        /*!
         * Return a copy of the micro displacement dof vector
         */

        return vectorTools::appendVectors( { _updatedFreeMicroDispDOFValues, _projected_ghost_micro_displacement } );
    }

    floatVector overlapCoupling::getUpdatedMacroDisplacementDOF( ){
        /*!
         * Return a copy of the macro displacement dof vector
         */

        return vectorTools::appendVectors( { _updatedFreeMacroDispDOFValues, _projected_ghost_macro_displacement } );
    }

    errorOut runOverlapCoupling( const std::string &filename,
                                 DOFMap &microGlobalLocalNodeMap, floatVector &updatedMicroDisplacementDOF,
                                 DOFMap &macroGlobalLocalNodeMap, floatVector &updatedMacroDisplacementDOF 
                               ){
        /*!
         * Run the overlap coupling method
         *
         * :param const std::string &filename: The name of the input YAML file
         * :param DOFMap &microGlobalLocalNodeMap: The map from global to local node numbers for the micro nodes
         * :param floatVector &updatedMicroDisplacementDOF: The updated micro displacement degrees of freedom
         * :param DOFMap &macroGlobalLocalNodeMap: The map from global to local node numbers for the macro nodes
         * :param floatVector &updatedMacroDisplacementDOF: The updated macro displacement degrees of freedom
         */

        //Construct the overlap coupling object
        overlapCoupling oc( filename );

        if ( oc.getConstructorError( ) ){

            errorOut result
                = new errorNode( "runOverlapCoupling", "Error in construction of overlapCoupling object" );

            result->addNext( oc.getConstructorError( ) );

            return result;

        }

        //Initialize the overlap coupling object
        errorOut error = oc.initializeCoupling( );
    
        if ( error ){
    
            errorOut result
                = new errorNode( "runOverlapCoupling", "Error in the initialization of the overlapCoupling object" );
    
            result->addNext( error );

            return result;
    
        }

        //Process the final increments of both the macro and micro-scales
        error = oc.processLastIncrements( );
    
        if ( error ){
    
            errorOut result
                = new errorNode( "runOverlapCoupling", "Error in processing the data" );
    
            result->addNext( error );

            return result;
    
        }

        //Return the updated DOF values
        microGlobalLocalNodeMap = oc.getMicroGlobalLocalNodeMap( );
        updatedMicroDisplacementDOF = oc.getUpdatedMicroDisplacementDOF( );

        macroGlobalLocalNodeMap = oc.getMacroGlobalLocalNodeMap( );
        updatedMacroDisplacementDOF = oc.getUpdatedMacroDisplacementDOF( );

        return NULL;

    }

    const cellDomainFloatVectorMap* overlapCoupling::getReferenceCellDomainCenterOfMassShapeFunctions( ){
        /*!
         * Get a constant reference to the shapefunctions of the centers of mass at each macro cell
         */

        return &_referenceCellDomainCenterOfMassShapefunctions;

    }

#ifdef TESTACCESS
    const std::unordered_map< uIntType, floatType >* overlapCoupling::getMacroNodeProjectedMass( ){
        /*!
         * Test access to macro node projected mass
         */

        return &_macroNodeProjectedMass;

    }
    const std::unordered_map< uIntType, floatVector >* overlapCoupling::getMacroNodeProjectedMassMomentOfInertia( ){
        /*!
         * Test access to macro node projected mass moment of inertia
         */

        return &_macroNodeProjectedMassMomentOfInertia;

    }
    const std::unordered_map< uIntType, floatVector >* overlapCoupling::getMacroNodeMassRelativePositionConstant( ){
        /*!
         * Test access to macro node mass relative position constant
         */

        return &_macroNodeMassRelativePositionConstant;

    }
    const std::unordered_map< uIntType, floatVector >* overlapCoupling::getMacroReferencePositions( ){
        /*!
         * Get the macro reference positions
         */

        return &_macroReferencePositions;

    }
    const std::unordered_map< uIntType, floatVector >* overlapCoupling::getMicroReferencePositions( ){
        /*!
         * Get the micro reference positions
         */

        return &_microReferencePositions;

    }
    const SparseMatrix *overlapCoupling::getCenterOfMassNMatrix( ){
        /*!
         * Get a constant reference to the center of mass interpolation matrix
         */

        return &_centerOfMassN;

    }
    const Eigen::MatrixXd *overlapCoupling::getCenterOfMassProjector( ){
        /*!
         * Get a constant reference to the center of mass projector
         */

        return &_centerOfMassProjector;

    }
    const SparseMatrix *overlapCoupling::getHomogenizationMatrix( ){
        /*!
         * Get a constant reference to the homogenization matrix
         */

        return &_homogenizationMatrix;

    }

    const cellDomainFloatMap* overlapCoupling::getHomogenizedVolumes( ){
        /*!
         * Get the homogenized volumes at the micro domains
         */

        return &homogenizedVolumes;
    }

    const cellDomainFloatMap* overlapCoupling::getHomogenizedDensities( ){
        /*!
         * Get the homogenized densities at the micro domains
         */

        return &homogenizedDensities;  
    }

    const cellDomainFloatVectorMap* overlapCoupling::getHomogenizedSymmetricMicroStresses( ){
        /*!
         * Get the homogenized symmetric micro stresses
         */

        return &homogenizedSymmetricMicroStresses;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getHomogenizedCentersOfMass( ){
        /*!
         * Get the homogenized centers of mass
         */

        return &homogenizedCentersOfMass;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getHomogenizedBodyForces( ){
        /*!
         * Get the homogenized body forces
         */

        return &homogenizedBodyForces;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getHomogenizedAccelerations( ){
        /*!
         * Get the homogenized accelerations
         */

        return &homogenizedAccelerations;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getHomogenizedMicroInertias( ){
        /*!
         * Get the homogenized micro inertias
         */

        return &homogenizedMicroInertias;
    }
    
    const cellDomainFloatVectorMap* overlapCoupling::getHomogenizedBodyForceCouples( ){
        /*!
         * Get the homogenized body force couples
         */

        return &homogenizedBodyForceCouples;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getHomogenizedMicroSpinInertias( ){
        /*!
         * Get the homogenized micro spin inertias
         */

        return &homogenizedMicroSpinInertias;
    }

#endif

}
