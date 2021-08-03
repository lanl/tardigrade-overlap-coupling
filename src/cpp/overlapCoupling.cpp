/*!============================================================================
|                               overlapCoupling                               |
===============================================================================
| The implementation of the Overlap coupling method                           |
=============================================================================*/

#include<overlapCoupling.h>
#include<Eigen/SparseQR>
#include<micromorphic_tools.h>
#include<balance_equations.h>

#include<boost/format.hpp>

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
            _error = new errorNode( __func__, "Error when setting the configuration filename" );
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
            errorOut result = new errorNode( __func__,
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

            errorOut result = new errorNode( __func__,
                                             "Error in getting the number of micro increments" );
            result->addNext( error );
            return result;

        }

        //Get the number of macro increments
        error = _inputProcessor._macroscale->getNumIncrements( numMacroIncrements );

        if ( error ){

            errorOut result = new errorNode( __func__,
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
            errorOut result = new errorNode( __func__, outstr );
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
        std::cerr << "INITIALIZE INCREMENT\n";
        errorOut error = _inputProcessor.initializeIncrement( microIncrement, macroIncrement );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in initialization of the input processor" );
            result->addNext( error );
            return result;

        }

        //Compute the centers of mass of the free and ghost domains
        std::cerr << "COMPUTE CENTERS OF MASS\n";
        error = computeIncrementCentersOfMass( microIncrement, macroIncrement,
                                               _freeMicroDomainMasses, _ghostMicroDomainMasses,
                                               _freeMicroDomainCentersOfMass, _ghostMicroDomainCentersOfMass );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in computation of the domain centers of mass" );
            result->addNext( error );
            return result;

        }

        //Project the degrees of freedom
        std::cerr << "PROJECT THE DEGREES OF FREEDOM\n";
        YAML::Node couplingConfiguration = _inputProcessor.getCouplingInitialization( );

        domainFloatVectorMap centerOfMassDisplacements;
        domainFloatVectorMap centerOfMassPhis;

        if ( _inputProcessor.isFiltering( ) ){

            // Solve for displacements directly
            floatVector domainAMatrix( _dim * _dim, 0);
            
            // Loop through the Ghost Cells
            const uIntVector *ghostMacroCellIds = _inputProcessor.getGhostMacroCellIds( );
    
            for ( auto cellID = ghostMacroCellIds->begin( ); cellID != ghostMacroCellIds->end( ); cellID++ ){
    
                // Get the centers of mass of the macro-domain
                auto referenceCellCentersOfMass = _referenceFreeMicroDomainCentersOfMass.find( *cellID );
    
                if ( referenceCellCentersOfMass == _referenceFreeMicroDomainCentersOfMass.end( ) ){
    
                    return new errorNode( __func__,
                                          "The macro cell " + std::to_string( *cellID ) +
                                          " was not found in the reference free micro domain centers of mass" );
    
                }

                // Get the moments of inertia of the macro-domain
                auto referenceMomentsOfInertia = _referenceFreeMicroDomainMomentsOfInertia.find( *cellID );
    
                if ( referenceMomentsOfInertia == _referenceFreeMicroDomainMomentsOfInertia.end( ) ){
    
                    return new errorNode( __func__,
                                          "The macro cell " + std::to_string( *cellID ) +
                                          " was not found in the reference free micro domain moments of inertia" );
    
                }
    
                // Loop through the micro domains
                for ( auto domain = referenceCellCentersOfMass->second.begin( ); domain != referenceCellCentersOfMass->second.end( ); domain++ ){

                    auto currentCenterOfMass = _freeMicroDomainCentersOfMass.find( domain->first );

                    if ( currentCenterOfMass == _freeMicroDomainCentersOfMass.end( ) ){

                        return new errorNode( __func__, "Micro domain " + domain->first +
                                                                  " not found in the current free micro centers of mass" );

                    }

                    auto domainMass = _freeMicroDomainMasses.find( domain->first );

                    if ( domainMass == _freeMicroDomainMasses.end( ) ){

                        return new errorNode( __func__, "Micro domain " + domain->first +
                                                                  " not found in the current free micro domain mass" );

                    }

                    auto domainMomentOfInertia = referenceMomentsOfInertia->second.find( domain->first );

                    if ( domainMomentOfInertia == referenceMomentsOfInertia->second.end( ) ){

                        return new errorNode( __func__, "Micro domain " + domain->first +
                                                                  " not found in the reference free reference moments of inertia" );

                    }

                    floatVector invI = vectorTools::inverse( domainMomentOfInertia->second, _dim, _dim );

                    // Compute the center of mass displacement
                    centerOfMassDisplacements.emplace( domain->first, currentCenterOfMass->second - domain->second );

                    // Loop through the micro nodes

                    //Get the domain node ids
                    uIntVector microDomainNodes;
                    errorOut error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, domain->first, microDomainNodes );
    
                    if ( error ){
    
                        errorOut result = new errorNode( __func__,
                                                         "Error in getting the node ids for the domain ( " + domain->first + " )" );
                        result->addNext( error );
                        return result;
    
                    }
    
                    const std::unordered_map< uIntType, floatVector > *microDisplacements = _inputProcessor.getMicroDisplacements( );
                    const std::unordered_map< uIntType, floatType >   *microVolumes       = _inputProcessor.getMicroVolumes( );
                    const std::unordered_map< uIntType, floatType >   *microDensities     = _inputProcessor.getMicroDensities( );
                    const std::unordered_map< uIntType, floatType >   *microWeights       = _inputProcessor.getMicroWeights( );
                      
                    domainAMatrix = floatVector( _dim * _dim, 0 );
                    for ( auto it = microDomainNodes.begin( ); it != microDomainNodes.end( ); it++ ){
    
                        auto microDisplacement = microDisplacements->find( *it );
    
                        if ( microDisplacement == microDisplacements->end( ) ){
    
                            return new errorNode( __func__, "Micro node " + std::to_string( *it ) +
                                                  " was not found in the micro displacement map" );
    
                        }

                        auto microVolume = microVolumes->find( *it );

                        if ( microVolume == microVolumes->end( ) ){

                            return new errorNode( __func__, "Micro node " + std::to_string( *it ) +
                                                                      " was not found in the micro volume map" );

                        }

                        auto microDensity = microDensities->find( *it );

                        if ( microDensity == microDensities->end( ) ){

                            return new errorNode( __func__, "Micro node " + std::to_string( *it ) +
                                                                      " was not found in the micro density map" );

                        }

                        auto microWeight = microWeights->find( *it );

                        if ( microWeight == microWeights->end( ) ){

                            return new errorNode( __func__, "Micro node " + std::to_string( *it ) +
                                                                      " was not found in the micro weight map" );

                        }

                        auto microReferencePosition = _inputProcessor.getMicroNodeReferencePositions( )->find( *it );

                        if ( microReferencePosition == _inputProcessor.getMicroNodeReferencePositions( )->end( ) ){

                            return new errorNode( __func__, "Micro node " + std::to_string( *it ) +
                                                                      " was not found in the micro reference position map" );

                        }

                        domainAMatrix += vectorTools::appendVectors( vectorTools::dyadic( microDisplacement->second - centerOfMassDisplacements[domain->first],
                                                                                          microReferencePosition->second - domain->second ) )
                                       * microVolume->second * microDensity->second * microWeight->second / domainMass->second;

                    }

                    centerOfMassPhis.emplace( domain->first, vectorTools::matrixMultiply( domainAMatrix, invI, _dim, _dim, _dim, _dim ) );

                }

            }

            // Project the center of mass DOF values to the nodes
            auto microDomainIDMap = _inputProcessor.getMicroDomainIDMap( );
            unsigned int nMacroDOF = _dim + ( _dim * _dim );

            Eigen::MatrixXd comValues( microDomainIDMap->size( ), 1 );
            Eigen::MatrixXd nodeValues;

            // Project the displacements
            _projected_ghost_macro_displacement.clear( );
            _projected_ghost_macro_displacement.resize( _inputProcessor.getGhostMacroNodeIds( )->size( ) * nMacroDOF );

            for ( unsigned int index = 0; index < _dim; index++ ){

                for ( auto domain = microDomainIDMap->begin( ); domain != microDomainIDMap->end( ); domain++ ){
    
                    auto us = centerOfMassDisplacements.find( domain->first );
    
                    if ( us == centerOfMassDisplacements.end( ) ){
    
                        return new errorNode( __func__, "Micro domain " + domain->first + " was not found in the center of mass displacement map" );
    
                    }
    
                    comValues( domain->second, 0 ) = us->second[ index ];
    
                }
   
                nodeValues = _centerOfMassProjector * comValues;

                // Assign the values
                for ( unsigned int i = 0; i < nodeValues.size( ); i++ ){

                    _projected_ghost_macro_displacement[ nMacroDOF * i + index ] = nodeValues( i );

                }

            }

            for ( unsigned int index = 0; index < _dim * _dim; index++ ){

                for ( auto domain = microDomainIDMap->begin( ); domain != microDomainIDMap->end( ); domain++ ){

                    auto phis = centerOfMassPhis.find( domain->first );

                    if ( phis == centerOfMassPhis.end( ) ){

                        return new errorNode( __func__, "Micro domain " + domain->first + " was not found in the phi map" );

                    }

                    comValues( domain->second, 0 ) = phis->second[ index ];

                }

                nodeValues = _centerOfMassProjector * comValues;

                // Assign the values
                for ( unsigned int i = 0; i < nodeValues.size( ); i++ ){

                    _projected_ghost_macro_displacement[ nMacroDOF * i + index + _dim ] = nodeValues( i );

                }

            }

            // Homogenize the response

            error = homogenizeMicroScale( microIncrement );

            if ( error ){
    
                errorOut result = new errorNode( __func__, "Error in the homogenization of the micro-scale to the macro-scale" );
                result->addNext( error );
                return result;
    
            }

        }
        else{

            if ( !couplingConfiguration[ "projection_type" ].as< std::string >( ).compare( "arlequin" ) == 0 ){
    
                error = applyKZProjection( macroIncrement, microIncrement );
    
                if ( error ){
        
                    errorOut result = new errorNode( __func__, "Error in applying the Klein-Zimmerman projection" );
                    result->addNext( error );
                    return result;
        
                }
    
            }
            else{
    
                error = applyArlequinProjection( macroIncrement, microIncrement );
    
                if ( error ){
    
                    errorOut result = new errorNode( __func__, "Error in applying the Arlequin projection" );
                    result->addNext( error );
                    return result;
    
                }
    
            }

        }

        if ( !couplingConfiguration [ "output_homogenized_response" ].IsScalar( ) ){

            //Output the homogenized material response to a data file
            error = outputHomogenizedResponse( );
            if ( error ){

                errorOut result = new errorNode( __func__, "Error when writing the homogenized response out to file" );
                result->addNext( error );
                return result;

            }

        }

        if ( !couplingConfiguration[ "output_updated_dof" ].IsScalar( ) ){

            //Output the updated dof values to a data file
            error = overlapCoupling::writeUpdatedDOFToFile( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error when writing the updated dof information to file" );
                result->addNext( error );
                return result;

            }

        }

        return NULL;
    }

    errorOut overlapCoupling::applyArlequinProjection( const uIntType &macroIncrement, const uIntType &microIncrement ){
        /*!
         * Apply the projection of Arlequin
         * to the current increment
         *
         * :param const uIntType &macroIncrement: The increment number for the macro-scale
         * :param const uIntType &microIncrement: The increment number for the micro-scale
         */

        //Compute the weighting factors for the micro nodes
        errorOut error = computeArlequinMicroWeightingFactors( microIncrement );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in computing the micro weighting factors" );
            result->addNext( error );
            return result;

        }

        // Compute the Arlequin coupling force

        error = computeArlequinCouplingForce( );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in computing the coupling force" );
            result->addNext( error );
            return result;

        }

//        //Compute the micromorphic mass matrix
//        error = computeArlequinMicromorphicMassMatrix( );
//
//        if ( error ){
//
//            errorOut result = new errorNode( __func__, "Error in the computation of the micromorphic mass matrix" );
//            result->addNext( error );
//            return result;
//
//        }
//
//        //Compute the remaining terms required for Arlequin
//        error = computeArlequinForceAndErrorVectors( );
//
//        if ( error ){
//
//            errorOut result = new errorNode( __func__, "Error in the computation of the force and error vectors" );
//            result->addNext( error );
//            return result;
//
//        }
//
//        //Update the degrees of freedom
//        error = computeArlequinDeformationUpdate( );
//
//        if ( error ){
//
//            errorOut result = new errorNode( __func__, "Error in the computation of the updated deformation values" );
//            result->addNext( error );
//            return result;
//
//        }

        //Homogenize the material properties at the micro-scale to the macro-scale
        std::cerr << "HOMOGENIZE THE MICROSCALE\n";
        error = homogenizeMicroScale( microIncrement );
    
        if ( error ){
    
            errorOut result = new errorNode( __func__, "Error in the homogenization of the micro-scale to the macro-scale" );
            result->addNext( error );
            return result;
    
        }

        return NULL;

    }

    errorOut overlapCoupling::computeArlequinMicroWeightingFactors( const uIntType &microIncrement ){
        /*!
         * Compute the weighting factor at each of the micro nodes
         *
         * Note that we will be using the basic structure from the Klein-Zimmerman formulation
         * ( i.e. ghost and free nodes ) but Arlequin doesn't make this distinction.
         *
         * :param const uIntType microIncrement: The current micro increment
         */

        //Get the macro cell ids
        const uIntVector *freeMacroCellIds = _inputProcessor.getFreeMacroCellIds( );
        const uIntVector *ghostMacroCellIds = _inputProcessor.getGhostMacroCellIds( );

        const uIntVector macroCellIds = vectorTools::appendVectors( { *freeMacroCellIds, *ghostMacroCellIds } );

        const std::unordered_map< uIntType, floatType > *macroArlequinWeights = _inputProcessor.getMacroArlequinWeights( );

        const std::unordered_map< uIntType, uIntType > *microGlobalToLocalDOFMap = _inputProcessor.getMicroGlobalToLocalDOFMap( );

        YAML::Node volumeReconstructionConfig = _inputProcessor.getVolumeReconstructionConfig( );

        errorOut error;

        arlequinMicroWeightingFactors.clear( );
        arlequinMicroWeightingFactors.reserve( microGlobalToLocalDOFMap->size( ) );

        //Loop over the macro cells
        for ( auto cellID = macroCellIds.begin( ); cellID != macroCellIds.end( ); cellID++ ){

            //Construct the element
            std::unique_ptr< elib::Element > element;
            error = buildMacroDomainElement( *cellID, *_inputProcessor.getMacroNodeReferencePositions( ),
                                             *_inputProcessor.getMacroDisplacements( ),
                                             *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                             element );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in building the macro-scale element" );
                result->addNext( error );
                return result;

            }

            //Reset the position to the element reference position
            element->nodes = element->reference_nodes;

            //Get the micro weights for the cell
            floatVector elementNodalWeights( element->global_node_ids.size( ) );

            for ( auto node = element->global_node_ids.begin( ); node != element->global_node_ids.end( ); node++ ){

                auto weight = macroArlequinWeights->find( *node );

                if ( weight == macroArlequinWeights->end( ) ){

                    return new errorNode( __func__,
                                          "Macro node " + std::to_string( *node ) +
                                          " was not found in the macro node to weighting factor map" );

                }

                elementNodalWeights[ node - element->global_node_ids.begin( ) ] = weight->second;

            }

            //Get the micro domains in the macro-cell
            auto domains = _inputProcessor.getMacroCellToDomainMap( )->find( *cellID );
            if ( domains == _inputProcessor.getMacroCellToDomainMap( )->end( ) ){

                return new errorNode( __func__,
                                      "The macro cell " + std::to_string( *cellID ) +
                                      " was not found in the cell ID to micro-domain map" );

            }

            //Loop over the domains contained in the macro-cell
            for ( auto domain = domains->second.begin( ); domain != domains->second.end( ); domain++ ){

                //Get the domain node ids
                uIntVector microDomainNodes;
                errorOut error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, *domain, microDomainNodes );

                if ( error ){

                    errorOut result = new errorNode( __func__,
                                                     "Error in getting the node ids for the domain ( " + *domain + " )" );
                    result->addNext( error );
                    return result;

                }

                unsigned int index = 0;
                const std::unordered_map< uIntType, floatVector > *microReferencePositions = _inputProcessor.getMicroNodeReferencePositions( );
                const std::unordered_map< uIntType, floatVector > *microDisplacements      = _inputProcessor.getMicroDisplacements( );
                uIntVector interiorNodes;

                for ( auto it = microDomainNodes.begin( ); it != microDomainNodes.end( ); it++, index++ ){

                    auto microReferencePosition = microReferencePositions->find( *it );

                    if ( microReferencePosition == microReferencePositions->end( ) ){

                        return new errorNode( __func__, "Micro node " + std::to_string( *it ) +
                                              " was not found in the micro reference position map" );

                    }

                    auto microDisplacement = microDisplacements->find( *it );

                    if ( microDisplacement == microDisplacements->end( ) ){

                        return new errorNode( __func__, "Micro node " + std::to_string( *it ) +
                                              " was not found in the micro displacement map" );

                    }

                    //compute the position
//                    floatVector globalCoordinates = microDisplacement->second + microReferencePosition->second;

                    //Get the local coordinates of the micro-node
                    floatVector localCoordinates;
                    std::unique_ptr< errorNode > tmp;
                    tmp.reset( element->compute_local_coordinates( microReferencePosition->second, localCoordinates ) );

                    if ( error ){
                           
                        if ( arlequinMicroWeightingFactors.find( *it ) == arlequinMicroWeightingFactors.end( ) ){
                        
                            std::cout << "error in local coordinates\n"; 
                            arlequinMicroWeightingFactors.emplace( *it, 0.001 ); //The micro-node is free

                        }

                        continue;

                    }

                    //See if the micro point is inside of the element
                    if ( !element->local_point_inside( localCoordinates, volumeReconstructionConfig[ "element_contain_tolerance" ].as< floatType >( ) ) ){
                        if ( arlequinMicroWeightingFactors.find( *it ) == arlequinMicroWeightingFactors.end( ) ){

                            arlequinMicroWeightingFactors.emplace( *it, 0.001 ); //The micro-node is free

                        }

                        continue;

                    }

                    //The point is contained in the element. We need to interpolate the weighting factors
                    floatType value;
                    element->interpolate( elementNodalWeights, localCoordinates, value );
                    value = std::fmax( value, 0.001 );
                    value = std::fmin( value, 0.999 );
                    arlequinMicroWeightingFactors.emplace( *it, value );
 
                }

            }

        }

        return NULL;

    }

    errorOut overlapCoupling::computeArlequinMicromorphicMassMatrix( ){
        /*!
         * Compute the Arlequin micromorphic mass matrix
         *
         * Note that we will be using the basic structure from the Klein-Zimmerman formulation
         * ( i.e. ghost and free nodes ) but Arlequin doesn't make this distinction.
         */

        const uIntVector *freeMacroCellIds = _inputProcessor.getFreeMacroCellIds( );
        const uIntVector *ghostMacroCellIds = _inputProcessor.getGhostMacroCellIds( );

        const uIntVector macroCellIds = vectorTools::appendVectors( { *freeMacroCellIds, *ghostMacroCellIds } );

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

        //Get the macro nodal Arlequin weights
        const std::unordered_map< uIntType, floatType > *macroArlequinWeights = _inputProcessor.getMacroArlequinWeights( );

        errorOut error;
//        floatVector arlequinWeights;
        std::vector< DOFProjection::T > coefficients;
        uIntType index = 0;
        for ( auto cellID = macroCellIds.begin(); cellID != macroCellIds.end(); cellID++, index++){

            //Construct the element
            std::unique_ptr< elib::Element > element;
            error = buildMacroDomainElement( *cellID, *_inputProcessor.getMacroNodeReferencePositions( ),
                                             *_inputProcessor.getMacroDisplacements( ),
                                             *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                             element );
    
            if ( error ){
    
                errorOut result = new errorNode( __func__,
                                                 "Error in building the macro-scale element" );
                result->addNext( error );
                return result;
    
            }

//            //Get the degree of freedom values for the free macro-domain element
//            arlequinWeights = floatVector( element->global_node_ids.size( ), 0 );
            const std::unordered_map< uIntType, floatVector >* macroDispDOFVector = _inputProcessor.getMacroDispDOFVector( );
            floatVector elementDOFVector( 0 );

            for ( auto nodeID  = element->global_node_ids.begin( );
                       nodeID != element->global_node_ids.end( );
                       nodeID++ ){

                auto macroDisplacement = macroDispDOFVector->find( *nodeID );

                if ( macroDisplacement == macroDispDOFVector->end( ) ){

                    return new errorNode( __func__,
                                          "Macro node " + std::to_string( *nodeID ) + " was not found in the macro displacement DOF vector map" );

                }

                elementDOFVector
                    = vectorTools::appendVectors( { elementDOFVector, macroDisplacement->second } );

                auto weight = macroArlequinWeights->find( *nodeID );

                if ( weight == macroArlequinWeights->end( ) ){

                    return new errorNode( __func__,
                                          "Macro node " + std::to_string( *nodeID ) +
                                          " was not found in the macro node to weighting factor map" );

                }
//
//                //Set the Arlequin weights at the nodes. We do force the weights to be between 0.001 and 0.999 so that
//                //The mass matrix won't be degenerate
//                arlequinWeights[ nodeID - element->global_node_ids.begin( ) ] = std::fmin( std::fmax( weight->second, 0.001 ), 0.999 );

            }

            //Extract the density and moment of inertia in the reference configuration
            auto densityType = macroReferenceDensityTypes->find( *cellID );

            if ( densityType == macroReferenceDensityTypes->end( ) ){

                return new errorNode( __func__,
                                      "The macro cell with ID " + std::to_string( *cellID ) +
                                      " was not found in the density type map" );

            }

            auto momentOfInertiaType = macroReferenceMomentOfInertiaTypes->find( *cellID );

            if ( momentOfInertiaType == macroReferenceMomentOfInertiaTypes->end( ) ){

                return new errorNode( __func__,
                                      "The macro cell with ID " + std::to_string( *cellID ) +
                                      " was not found in the moment of inertia type map" );

            }

            if ( densityType->second.compare( "constant" ) != 0 ){

                return new errorNode( __func__,
                                      "Only constant densities for the macro-scale are allowed currently. This is not true for macro cell ID " + std::to_string( *cellID ) );

            }

            if ( momentOfInertiaType->second.compare( "constant" ) != 0 ){

                return new errorNode( __func__,
                                      "Only constant moments of inertia for the macro-scale are allowed currently. This is not true for macro cell ID " + std::to_string( *cellID ) );

            }

            auto macroDensities = macroReferenceDensities->find( *cellID );

            if ( macroDensities == macroReferenceDensities->end( ) ){

                return new errorNode( __func__,
                                      "Macro cell ID " + std::to_string( *cellID ) +
                                      " is not in the macro reference density map" );

            }

            if ( macroDensities->second.size( ) != 1 ){

                return new errorNode( __func__,
                                      "The macro densities for macro cell " + std::to_string( *cellID ) +
                                      "Define " + std::to_string( macroDensities->second.size( ) ) +
                                      " values when only 1 can be defined" );

            }

            auto macroMomentsOfInertia = macroReferenceMomentsOfInertia->find( *cellID );

            if ( macroMomentsOfInertia == macroReferenceMomentsOfInertia->end( ) ){

                return new errorNode( __func__,
                                      "Macro cell ID " + std::to_string( *cellID ) +
                                      " is not in the macro reference moments of inertia map" );

            }

            if ( macroMomentsOfInertia->second.size( ) != _dim * _dim ){

                return new errorNode( __func__,
                                      "The macro moments of inertia for macro cell " + std::to_string( *cellID ) +
                                      "Define " + std::to_string( macroDensities->second.size( ) ) +
                                      " values when only " + std::to_string( _dim * _dim ) + " can be defined" );

            }

            bool quantitiesInReference = true; //TODO: This is hard coded for now but, in the future when functionally-graded micromorphic is implemented, we will want to change this.

            floatVector densities( element->qrule.size( ), macroDensities->second[ 0 ] );

            floatVector momentsOfInertia
                = vectorTools::appendVectors( floatMatrix( element->qrule.size( ), macroMomentsOfInertia->second ) );

            //Initialize the coefficients vector
            coefficients.clear( );

            uIntType numCoefficients = 0;
            for ( auto it = element->global_node_ids.begin( ); it != element->global_node_ids.end( ); it++ ){
    
                //Get the number of nodes in the element
                uIntType elementNodeCount = element->global_node_ids.size( );
                numCoefficients += element->qrule.size( ) * elementNodeCount * elementNodeCount * _dim * _dim * ( 1 + _dim * _dim );
    
            }
    
            coefficients.reserve( numCoefficients );

            error = formMicromorphicElementMassMatrix( element, elementDOFVector, momentsOfInertia, densities, 
                                                       _inputProcessor.getMacroGlobalToLocalDOFMap( ), coefficients,
//                                                       &arlequinWeights, quantitiesInReference );
                                                       NULL, quantitiesInReference );

            if ( error ){

                std::string outstr  = "Error in the construction of the contributions of the macro element to ";
                            outstr += "the micromorphic mass matrix";

                errorOut result = new errorNode( __func__, outstr );
                result->addNext( error );
                return result;

            }

            SparseMatrix _MDe( nMacroDOF * _inputProcessor.getMacroGlobalToLocalDOFMap( )->size( ),
                               nMacroDOF * _inputProcessor.getMacroGlobalToLocalDOFMap( )->size( ) );
            _MDe.setFromTriplets( coefficients.begin( ), coefficients.end( ) );
            if( index==0 ){
                _MD = _MDe;
            }
            else{
                _MD += _MDe;
            }

        }

        return NULL;

    }

    errorOut overlapCoupling::computeArlequinCouplingForce( ){
        /*!
         * Compute the Arlequin coupling force vector
         */

        //Set the configuration
        const YAML::Node config = _inputProcessor.getCouplingInitialization( );

        std::string projection_type = config[ "projection_type" ].as< std::string >( );

        //Set the number of displacement degrees of freedom for each scale
        const uIntType nMicroDispDOF = _dim;
        const uIntType nMacroDispDOF = _dim + _dim * _dim;

        //Get the maps from the global to the local degrees of freedom
        const DOFMap *microGlobalToLocalDOFMap = _inputProcessor.getMicroGlobalToLocalDOFMap( );
        const DOFMap *macroGlobalToLocalDOFMap = _inputProcessor.getMacroGlobalToLocalDOFMap( );

        //Get the macro properties
        const std::unordered_map< uIntType, floatType > *macroArlequinWeights = _inputProcessor.getMacroArlequinWeights( );

        //Get the forces
        const std::unordered_map< uIntType, floatVector > *microInternalForces = _inputProcessor.getMicroInternalForces( );
        const std::unordered_map< uIntType, floatVector > *microInertialForces = _inputProcessor.getMicroInertialForces( );
        const std::unordered_map< uIntType, floatVector > *microExternalForces = _inputProcessor.getMicroExternalForces( );

        //Get the degree of freedom vectors
        const std::unordered_map< uIntType, floatVector > *microDisplacements = _inputProcessor.getMicroDisplacements( );

        const std::unordered_map< uIntType, floatVector > *macroInternalForces = _inputProcessor.getMacroInternalForces( );
        const std::unordered_map< uIntType, floatVector > *macroInertialForces = _inputProcessor.getMacroInertialForces( );
        const std::unordered_map< uIntType, floatVector > *macroExternalForces = _inputProcessor.getMacroExternalForces( );
        const std::unordered_map< uIntType, floatVector > *macroDispDOFVector = _inputProcessor.getMacroDispDOFVector( );

        const floatType mu = *_inputProcessor.getArlequinPenaltyParameter( );

        //Size the Lagrangian
        floatVector RHS( nMicroDispDOF * microGlobalToLocalDOFMap->size( ), 0 );

        floatVector G( nMicroDispDOF * microGlobalToLocalDOFMap->size( ), 0 );

        //Assemble the micro-force vector
//        std::cout << "Assembling the micro-force vector:\n";
        for( auto node = microGlobalToLocalDOFMap->begin( ); node != microGlobalToLocalDOFMap->end( ); node++ ){

//            std::cerr << "  node: " << node->first << "\n";// << ": " << node->second << "\n";
//            std::cerr << "    reference position: "; vectorTools::print( _inputProcessor.getMicroNodeReferencePositions( )->find( node->first )->second );
            floatVector nodeForce( nMicroDispDOF, 0 );

            //Get the external force contribution
            if( _inputProcessor.microExternalForceDefined( ) ){

                auto externalForce = microExternalForces->find( node->first );
//                std::cerr << "    external force: "; vectorTools::print( externalForce->second );
    
                if ( externalForce == microExternalForces->end( ) ){
    
                    return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                          " not found in external force vector" );
    
                }

                nodeForce -= externalForce->second;

            }

            //Get the internal force contribution
            if( _inputProcessor.microInternalForceDefined( ) ){

                auto internalForce = microInternalForces->find( node->first );
//                std::cerr << "    internal force: "; vectorTools::print( internalForce->second );
    
                if ( internalForce == microInternalForces->end( ) ){
    
                    return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                          " not found in internal force vector" );
    
                }

                nodeForce += internalForce->second;

            }

            if( _inputProcessor.microInertialForceDefined( ) ){

                auto inertialForce = microInertialForces->find( node->first );

                if ( inertialForce == microInertialForces->end( ) ){

                    return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                          " not found in inertial force vector" );

                }

                nodeForce += inertialForce->second;

            }

            auto microDisplacement = microDisplacements->find( node->first );

            if ( microDisplacement == microDisplacements->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                      " not found in micro displacement vector" );

            }

//            std::cerr << "    displacement: "; vectorTools::print( microDisplacement->second );

            //Get the arlequin weight
            auto arlequinWeight = arlequinMicroWeightingFactors.find( node->first );
            if ( arlequinWeight == arlequinMicroWeightingFactors.end( ) ){

                return new errorNode( __func__,
                                      "Micro node " + std::to_string( node->first ) + " not found in arlequin weight vector" );

            }

            for ( uIntType i = 0; i < nMicroDispDOF; i++ ){

                RHS[ nMicroDispDOF * node->second + i ] += ( 1 - arlequinWeight->second ) * nodeForce[ i ] + mu * microDisplacement->second[i];
                G[ nMicroDispDOF * node->second + i ] -= microDisplacement->second[ i ];

            }

        }

        //Assemble the macro-force and error vector
        floatVector RHS_D( nMacroDispDOF * macroGlobalToLocalDOFMap->size( ), 0 );
        floatVector D( nMacroDispDOF * macroGlobalToLocalDOFMap->size( ), 0 );

        // Assemble the macro-displacement vector
        for( auto node = macroGlobalToLocalDOFMap->begin( ); node != macroGlobalToLocalDOFMap->end( ); node++ ){

            floatVector nodeForce( nMacroDispDOF, 0 );

//            std::cout << "node: " << node->first << "\n";
            //Add the external force contribution
            if ( _inputProcessor.macroExternalForceDefined( ) ){

                auto externalForce = macroExternalForces->find( node->first );

                if ( externalForce == macroExternalForces->end( ) ){

                    return new errorNode( __func__,
                                          "Macro node " + std::to_string( node->first ) + " not found in external force map" );

                }

                nodeForce -= externalForce->second;

            }

            //Add the internal force contribution
            if ( _inputProcessor.macroInternalForceDefined( ) ){

                auto internalForce = macroInternalForces->find( node->first );
//                std::cerr << "    internal force: "; vectorTools::print( internalForce->second );

                if ( internalForce == macroInternalForces->end( ) ){

                    return new errorNode( __func__,
                                          "Macro node " + std::to_string( node->first ) + " not found in internal force map" );

                }

                nodeForce += internalForce->second;

            }

            // Add the inertial force contribution
            if ( _inputProcessor.macroInertialForceDefined( ) ){

                auto inertialForce = macroInertialForces->find( node->first );

                if ( inertialForce == macroInertialForces->end( ) ){

                    return new errorNode( __func__,
                                          "Macro node " + std::to_string( node->first ) + " not found in inertial force map" );

                }

                nodeForce += inertialForce->second;

            }

            // Add the deformation force
            auto macroDispDOF = macroDispDOFVector->find( node->first );

            if ( macroDispDOF == macroDispDOFVector->end( ) ){

                return new errorNode( __func__, "Macro node " + std::to_string( node->first ) +
                                      " not found in macro displacement DOF map" );

            }

//            std::cerr << "    deformation: "; vectorTools::print(macroDispDOF->second);

            //Get the Macro Arlequin weight
            auto arlequinWeight =  macroArlequinWeights->find( node->first );
            if ( arlequinWeight == macroArlequinWeights->end( ) ){

                return new errorNode( __func__,
                                      "Macro node " + std::to_string( node->first ) + " not found in arlequin weight vector" );

            }

            for ( uIntType i = 0; i < nMacroDispDOF; i++ ){

                RHS_D[ nMacroDispDOF * node->second + i ] += arlequinWeight->second * nodeForce[ i ];
                D[ nMacroDispDOF * node->second + i ] = macroDispDOF->second[ i ];

            }

        }

//        std::cerr << "out of loop\n";

        Eigen::Map< Eigen::Matrix< floatType, -1,  1 > > _RHS( RHS.data(), RHS.size( ), 1 );
        Eigen::Map< Eigen::Matrix< floatType, -1,  1 > > _RHS_D( RHS_D.data(), RHS_D.size( ), 1 );

        // Assemble the coupling force vectors
        FALQ.resize( RHS.size( ) );
        FALD.resize( RHS_D.size( ) );

        Eigen::Map< Eigen::Matrix< floatType, -1,  1 > > microCouplingForce( FALQ.data(), FALQ.size( ), 1 );
        Eigen::Map< Eigen::Matrix< floatType, -1,  1 > > macroCouplingForce( FALD.data(), FALD.size( ), 1 );
        Eigen::Map< Eigen::Matrix< floatType, -1, 1 > > _G( G.data( ), G.size( ), 1 );
        Eigen::Map< Eigen::Matrix< floatType, -1, 1 > > _D( D.data( ), D.size( ), 1 );
        std::cerr << "forming G\n";
        _G += _N * _D;

        std::cerr << "forming RHS\n";
        _RHS -= ( _N * _RHS_D + mu * ( _N * _N.transpose( ) * _G + _G ) );

        std::cerr << "forming LHS\n";
        SparseMatrix _LHS = _N * _N.transpose( );
        _LHS.diagonal( ).array( ) += 1;

        std::cerr << "diagonalizing LHS\n";
        Eigen::Matrix< floatType, -1, 1 > _LHS_vec = _LHS * Eigen::VectorXd::Ones( _LHS.cols( ) );

        std::cerr << "building micro coupling force\n";
        microCouplingForce = _RHS.array( ) / _LHS_vec.array( );
        std::cerr << "building macro coupling force\n";
        macroCouplingForce = _N.transpose( ) * microCouplingForce;

//        std::cerr << "outputting G\n";
//        for ( unsigned int i = 0; i < microGlobalToLocalDOFMap->size( ); i++ ){
//
//            for ( unsigned int j = 0; j < nMicroDispDOF; j++ ){
//                std::cerr << G[ nMicroDispDOF * i + j ] << ", ";
//            }
//            std::cerr << "\n";
//
//        }
//
//        std::cerr << "microCouplingForce:\n";
//        for ( unsigned int i = 0; i < microGlobalToLocalDOFMap->size( ); i++ ){
//
//            for ( unsigned int j = 0; j < nMicroDispDOF; j++ ){
//                std::cerr << microCouplingForce[ nMicroDispDOF * i + j ] << ", ";
//            }
//            std::cerr << "\n";
//
//        }
//
//        std::cerr << "macroCouplingForce:\n";
//        for ( unsigned int i = 0; i < macroGlobalToLocalDOFMap->size( ); i++ ){
//
//            for ( unsigned int j = 0; j < nMacroDispDOF; j++ ){
//                std::cerr << macroCouplingForce[ nMacroDispDOF * i + j ] << ", ";
//            }
//            std::cerr << "\n";
//
//        }
//
//        std::cerr << "G norm: " << vectorTools::l2norm( G ) << "\n";
//
//        std::cerr << "micro-domains:\n";
//        for ( auto domain = _inputProcessor.getMicroDomainIDMap( )->begin( ); domain != _inputProcessor.getMicroDomainIDMap( )->end( ); domain++ ){
//
//            uIntVector domainNodes;
//            errorOut error = _inputProcessor._microscale->getSubDomainNodes( 0, domain->first, domainNodes );
//            if ( error ){
//                return error;
//            }
//            std::cerr << domain->first << ": "; vectorTools::print( domainNodes );
//            std::cerr << "    reference positions\n";
//            auto referencePositions = _inputProcessor.getMicroNodeReferencePositions( );
//            for ( auto d = domainNodes.begin( ); d != domainNodes.end( ); d++ ){std::cerr << "        "; auto p = referencePositions->find( *d ); vectorTools::print( p->second ); };
//
//            std::cerr << "    center of mass: ";
//            for ( auto cell = _referenceFreeMicroDomainCentersOfMass.begin( ); cell != _referenceFreeMicroDomainCentersOfMass.end( ); cell++ ){
//
//                std::cout << "cell: " << cell->first << "\n";
//
//                auto centerOfMass = cell->second.find( domain->first );
//
//                if ( centerOfMass != cell->second.end( ) ){
//
//                    vectorTools::print( centerOfMass->second );
//
//                }
//
//            }
//
//            for ( auto cell = _referenceGhostMicroDomainCentersOfMass.begin( ); cell != _referenceGhostMicroDomainCentersOfMass.end( ); cell++ ){
//
//                std::cout << "cell: " << cell->first << "\n";
//
//                auto centerOfMass = cell->second.find( domain->first );
//
//                if ( centerOfMass != cell->second.end( ) ){
//
//                    vectorTools::print( centerOfMass->second );
//
//                }
//
//            }
//
//            std::cout << "DONE\n";
//
//        }
        
//
//        std::cerr << "microCouplingForce:\n" << microCouplingForce << "\n";
//        std::cerr << "macroCouplingForce:\n" << macroCouplingForce << "\n";
//        assert( 1 == 0 );

        return NULL;

    }

    errorOut overlapCoupling::computeArlequinForceAndErrorVectors( ){
        /*!
         * Compute the Arlequin force vectors along with the error vector
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

        //Get the micro properties
        const std::unordered_map< uIntType, floatType > *microDensities = _inputProcessor.getMicroDensities( );
        const std::unordered_map< uIntType, floatType > *microVolumes = _inputProcessor.getMicroVolumes( );

        //Get the macro properties
        const std::unordered_map< uIntType, floatType > *macroArlequinWeights = _inputProcessor.getMacroArlequinWeights( );

        //Get the forces
        const std::unordered_map< uIntType, floatVector > *microInternalForces = _inputProcessor.getMicroInternalForces( );
        const std::unordered_map< uIntType, floatVector > *microExternalForces = _inputProcessor.getMicroExternalForces( );

        const std::unordered_map< uIntType, floatVector > *macroInternalForces = _inputProcessor.getMacroInternalForces( );
        const std::unordered_map< uIntType, floatVector > *macroExternalForces = _inputProcessor.getMacroExternalForces( );


        //Get the degree of freedom vectors
        const std::unordered_map< uIntType, floatVector > *previousMicroDisplacements = _inputProcessor.getPreviousMicroDisplacements( );
        const std::unordered_map< uIntType, floatVector > *previousMicroVelocities    = _inputProcessor.getPreviousMicroVelocities( );
        const std::unordered_map< uIntType, floatVector > *previousMicroAccelerations = _inputProcessor.getPreviousMicroAccelerations( );

        const std::unordered_map< uIntType, floatVector > *previousMacroDispDOFVector = _inputProcessor.getPreviousMacroDispDOFVector( );
        const std::unordered_map< uIntType, floatVector > *previousMacroVelocities    = _inputProcessor.getPreviousMacroVelocities( );
        const std::unordered_map< uIntType, floatVector > *previousMacroAccelerations = _inputProcessor.getPreviousMacroAccelerations( );

        floatType aQ = config[ "micro_proportionality_coefficient" ].as< floatType >( );
        floatType aD = config[ "macro_proportionality_coefficient" ].as< floatType >( );

        Qe = floatVector( nMicroDispDOF * microGlobalToLocalDOFMap->size( ), 0 );
        De = floatVector( nMacroDispDOF * macroGlobalToLocalDOFMap->size( ), 0 );

        FQ = floatVector( nMicroDispDOF * microGlobalToLocalDOFMap->size( ), 0 );
        FD = floatVector( nMacroDispDOF * macroGlobalToLocalDOFMap->size( ), 0 );

        //Assemble the micro-force vector
//        std::cout << "Assembling the micro-force vector:\n";
        for( auto node = microGlobalToLocalDOFMap->begin( ); node != microGlobalToLocalDOFMap->end( ); node++ ){

//            std::cerr << "  node: " << node->first << ": " << node->second << "\n";
//            std::cerr << "    reference position: "; vectorTools::print( _inputProcessor.getMicroNodeReferencePositions( )->find( node->first )->second );
            floatVector nodeForce( nMicroDispDOF, 0 );

            //Get the external force contribution
            if( _inputProcessor.microExternalForceDefined( ) ){

                auto externalForce = microExternalForces->find( node->first );
//                std::cerr << "    external force: "; vectorTools::print( externalForce->second );
    
                if ( externalForce == microExternalForces->end( ) ){
    
                    return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                          " not found in external force vector" );
    
                }

                nodeForce += externalForce->second;

            }

            //Get the internal force contribution
            if( _inputProcessor.microInternalForceDefined( ) ){

                auto internalForce = microInternalForces->find( node->first );
//                std::cerr << "    internal force: "; vectorTools::print( internalForce->second );
    
                if ( internalForce == microInternalForces->end( ) ){
    
                    return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                          " not found in internal force vector" );
    
                }

                nodeForce -= internalForce->second;

            }

            auto microDensity = microDensities->find( node->first );

            if ( microDensity == microDensities->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                      " not found in micro densities vector" );

            }

            auto microVolume = microVolumes->find( node->first );

            if ( microVolume == microVolumes->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                      " not found in micro volumes vector" );

            }

            floatVector qe( 3, 0 );
            auto microDisplacement = previousMicroDisplacements->find( node->first );

            if ( microDisplacement == previousMicroDisplacements->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                      " not found in micro displacement vector" );

            }

            qe = microDisplacement->second;

            //Get the velocity contribution to the damping force and error
            if ( _inputProcessor.microVelocitiesDefined( ) ){

                auto microVelocity = previousMicroVelocities->find( node->first );

                if ( microVelocity == previousMicroVelocities->end( ) ){

                    return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                          " not found in previous micro velocities vector" );

                }

                nodeForce -= aQ * microDensity->second * microVolume->second * microVelocity->second;

                qe += ( *dt ) * microVelocity->second;

            }

            //Get the acceleration contribution to the damping force
            if ( _inputProcessor.microAccelerationDefined( ) ){

                auto microAcceleration = previousMicroAccelerations->find( node->first );

                if ( microAcceleration == previousMicroAccelerations->end( ) ){

                    return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                          " not found in previous micro accelerations vector" );

                }

                nodeForce -= aQ * microDensity->second * microVolume->second * ( 1 - gamma ) * ( *dt ) * microAcceleration->second;

                qe += 0.5 * ( *dt ) * ( *dt ) * ( 1 - 2 * beta ) * microAcceleration->second;

            }

            //Get the arlequin weight
            auto arlequinWeight = arlequinMicroWeightingFactors.find( node->first );
//            std::cout << "    arlequin weight: " << arlequinWeight->second << "\n";
            if ( arlequinWeight == arlequinMicroWeightingFactors.end( ) ){

                return new errorNode( __func__,
                                      "Micro node " + std::to_string( node->first ) + " not found in arlequin weight vector" );

            }

            for ( uIntType i = 0; i < nMicroDispDOF; i++ ){

                FQ[ nMicroDispDOF * node->second + i ] += ( 1 - arlequinWeight->second ) * nodeForce[ i ];
                Qe[ nMicroDispDOF * node->second + i ] += qe[ i ];

            }

        }

        //Assemble the macro-force and error vector
        floatVector explicitVelocity;

        if ( ( _inputProcessor.macroVelocitiesDefined( ) ) || ( _inputProcessor.macroAccelerationDefined( ) ) ){
            explicitVelocity = floatVector( nMacroDispDOF * macroGlobalToLocalDOFMap->size( ), 0 );
        }

//        std::cout << "Assembling the macro-force vector:\n";
        for( auto node = macroGlobalToLocalDOFMap->begin( ); node != macroGlobalToLocalDOFMap->end( ); node++ ){

            floatVector nodeForce( nMacroDispDOF, 0 );
//            std::cerr << "  node: " << node->first << ": " << node->second << "\n";

//            std::cerr << "    reference position: "; vectorTools::print( _inputProcessor.getMacroNodeReferencePositions( )->find( node->first )->second );

            //Add the external force contribution
            if ( _inputProcessor.macroExternalForceDefined( ) ){

                auto externalForce = macroExternalForces->find( node->first );
                std::cerr << "    external force: "; vectorTools::print( externalForce->second );

                if ( externalForce == macroExternalForces->end( ) ){

                    return new errorNode( __func__,
                                          "Macro node " + std::to_string( node->first ) + " not found in external force map" );

                }

                nodeForce += externalForce->second;

            }

            //Add the internal force contribution
            if ( _inputProcessor.macroInternalForceDefined( ) ){

                auto internalForce = macroInternalForces->find( node->first );
//                std::cerr << "    internal force: "; vectorTools::print( internalForce->second );

                if ( internalForce == macroInternalForces->end( ) ){

                    return new errorNode( __func__,
                                          "Macro node " + std::to_string( node->first ) + " not found in internal force map" );

                }

                nodeForce -= internalForce->second;

            }

            auto macroDisplacement = previousMacroDispDOFVector->find( node->first );

            if ( macroDisplacement == previousMacroDispDOFVector->end( ) ){

                return new errorNode( __func__, "Macro node " + std::to_string( node->first ) +
                                      " not found in macro displacement vector" );

            }

            for ( uIntType i = 0; i < macroDisplacement->second.size( ); i++ ){

                De[ nMacroDispDOF * node->second + i ] += macroDisplacement->second[ i ];

            }

            //Add the contribution from the previous velocity term
            if ( _inputProcessor.macroVelocitiesDefined( ) ){

                auto previousVelocity = previousMacroVelocities->find( node->first );

                if ( previousVelocity == previousMacroVelocities->end( ) ){

                    return new errorNode( __func__,
                                          "Macro node " + std::to_string( node->first ) + " not found in previous macro velocity map" );

                }

                for ( uIntType i = 0; i < previousVelocity->second.size( ); i++ ){

                    explicitVelocity[ nMacroDispDOF * node->second + i ] += previousVelocity->second[ i ];
                    De[ nMacroDispDOF * node->second + i ] += ( *dt ) * previousVelocity->second[ i ];

                }

            }

            //Add the contribution from the previous acceleration term
            if ( _inputProcessor.macroAccelerationDefined( ) ){

                auto previousAcceleration = previousMacroAccelerations->find( node->first );

                if ( previousAcceleration == previousMacroAccelerations->end( ) ){

                    return new errorNode( __func__,
                                          "Macro node " + std::to_string( node->first ) + " not found in previous macro acceleration map" );

                }

                for ( uIntType i = 0; i < previousAcceleration->second.size( ); i++ ){

                    explicitVelocity[ nMacroDispDOF * node->second + i ]
                        += ( 1 - gamma ) * ( *dt ) * previousAcceleration->second[ i ];
                    De[ nMacroDispDOF * node->second + i ] += 0.5 * ( *dt ) * ( *dt ) * ( 1 - 2 * beta )
                                                            * previousAcceleration->second[ i ];

                }

            }

            //Add the lumped inertia contribution
            for ( uIntType i = 0; i < nMacroDispDOF; i++ ){
//                nodeForce[ i ] -= ( 1 - aD * gamma * ( *dt ) ) * explicitVelocity[ nMacroDispDOF * node->second + i ]
//                                * _MD.row( nMacroDispDOF * node->second + i ).sum( );
                nodeForce[ i ] -= aD * explicitVelocity[ nMacroDispDOF * node->second + i ]
                                * _MD.row( nMacroDispDOF * node->second + i ).sum( );
            }

            //Get the Arlequin weight
            auto arlequinWeight = macroArlequinWeights->find( node->first );
//            std::cerr << "    arlequin weight: " << arlequinWeight->second << "\n";

            if ( arlequinWeight == macroArlequinWeights->end( ) ){

                return new errorNode( __func__,
                                      "Macro node " + std::to_string( node->first ) + " not found in arlequin weights map" );

            }

            for ( uIntType i = 0; i < nMacroDispDOF; i++ ){

                FD[ nMacroDispDOF * node->second + i ]
                    += std::fmin( std::fmax( arlequinWeight->second, 0.001 ), 0.999 ) * nodeForce[ i ];

            }


        }

        //I don't think this should be here
//        if ( ( _inputProcessor.macroVelocitiesDefined( ) ) || ( _inputProcessor.macroAccelerationDefined( ) ) ){
//
//            
//            Eigen::Map< Eigen::Matrix< floatType, -1,  1 > > _explicitVelocity( explicitVelocity.data(), explicitVelocity.size( ), 1 );
//            Eigen::Map< Eigen::Matrix< floatType, -1,  1 > > _FD( FD.data( ), FD.size( ), 1 );
//
//            _FD -= aD * _MD * _explicitVelocity;
//
//        }

        return NULL;
    }

    errorOut overlapCoupling::computeArlequinDeformationUpdate( ){
        /*!
         * Compute the deformation update using the Arlequin method
         */

        //Set the configuration
        const YAML::Node config = _inputProcessor.getCouplingInitialization( );

        std::string projection_type = config[ "projection_type" ].as< std::string >( );
        const floatType gamma = *_inputProcessor.getNewmarkGamma( );
        const floatType beta = *_inputProcessor.getNewmarkBeta( );
        const floatType mu = *_inputProcessor.getArlequinPenaltyParameter( );
        const floatType mu_update = *_inputProcessor.getArlequinUpdatePenaltyParameter( );

        //Get the timestep
        const floatType *dt = _inputProcessor.getDt( );

        //Get the proportionality constants
        floatType aQ = config[ "micro_proportionality_coefficient" ].as< floatType >( );
        floatType aD = config[ "macro_proportionality_coefficient" ].as< floatType >( );

        //Set the micro mass vector
        uIntType nMicroDispDOF = _dim;
        uIntType nMacroDispDOF = _dim + _dim * _dim;

        //Get the DOF maps
        const std::unordered_map< uIntType, uIntType > *microGlobalToLocalDOFMap = _inputProcessor.getMicroGlobalToLocalDOFMap( );
        const std::unordered_map< uIntType, uIntType > *macroGlobalToLocalDOFMap = _inputProcessor.getMacroGlobalToLocalDOFMap( );

        //Required properties to compute the micro mass vector
        const std::unordered_map< uIntType, floatType > *microDensities = _inputProcessor.getMicroDensities( );
        const std::unordered_map< uIntType, floatType > *microVolumes = _inputProcessor.getMicroVolumes( );

        //Get the required DOF values
        const std::unordered_map< uIntType, floatVector > *microDisplacements         = _inputProcessor.getMicroDisplacements( );
        const std::unordered_map< uIntType, floatVector > *microVelocities            = _inputProcessor.getMicroVelocities( );
        const std::unordered_map< uIntType, floatVector > *microAccelerations         = _inputProcessor.getMicroAccelerations( );
        const std::unordered_map< uIntType, floatVector > *previousMicroDisplacements = _inputProcessor.getPreviousMicroDisplacements( );
        const std::unordered_map< uIntType, floatVector > *previousMicroVelocities    = _inputProcessor.getPreviousMicroVelocities( );
        const std::unordered_map< uIntType, floatVector > *previousMicroAccelerations = _inputProcessor.getPreviousMicroAccelerations( );

        const std::unordered_map< uIntType, floatVector > *macroDispDOFVector         = _inputProcessor.getMacroDispDOFVector( );
        const std::unordered_map< uIntType, floatVector > *macroVelocities            = _inputProcessor.getMacroVelocities( );
        const std::unordered_map< uIntType, floatVector > *macroAccelerations         = _inputProcessor.getMacroAccelerations( );
        const std::unordered_map< uIntType, floatVector > *previousMacroDispDOFVector = _inputProcessor.getPreviousMacroDispDOFVector( );
        const std::unordered_map< uIntType, floatVector > *previousMacroVelocities    = _inputProcessor.getPreviousMacroVelocities( );
        const std::unordered_map< uIntType, floatVector > *previousMacroAccelerations = _inputProcessor.getPreviousMacroAccelerations( );
        const std::unordered_map< uIntType, floatType > *macroArlequinWeights = _inputProcessor.getMacroArlequinWeights( );

        Eigen::VectorXd Q( nMicroDispDOF * microGlobalToLocalDOFMap->size( ) );
        Eigen::VectorXd QDot( nMicroDispDOF * microGlobalToLocalDOFMap->size( ) );
        Eigen::VectorXd QDotDot( nMicroDispDOF * microGlobalToLocalDOFMap->size( ) );

        Eigen::VectorXd MQ( nMicroDispDOF * microGlobalToLocalDOFMap->size( ) );
        Eigen::VectorXd WQ( nMicroDispDOF * microGlobalToLocalDOFMap->size( ) );
        std::cerr << "Building Q and MQ\n";
        for( auto node = microGlobalToLocalDOFMap->begin( ); node != microGlobalToLocalDOFMap->end( ); node++ ){

            auto microDensity = microDensities->find( node->first );

            if ( microDensity == microDensities->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                      " not found in micro densities vector" );

            }

            auto microVolume = microVolumes->find( node->first );

            if ( microVolume == microVolumes->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                      " not found in micro volumes vector" );

            }

            auto microWeight = arlequinMicroWeightingFactors.find( node->first );

            if ( microWeight == arlequinMicroWeightingFactors.end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                      " not found in micro weighting factor vector" );

            }

            for ( uIntType i = 0; i < nMicroDispDOF; i++ ){

                MQ( nMicroDispDOF * node->second + i ) = ( 1 - microWeight->second ) * microDensity->second * microVolume->second;

            }

            auto microDisplacement = microDisplacements->find( node->first );

            if ( microDisplacement == microDisplacements->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                      " not found in micro displacements vector" );

            }

            for ( uIntType i = 0; i < nMicroDispDOF; i++ ){

                Q( nMicroDispDOF * node->second + i ) = microDisplacement->second[ i ];

            }

            auto microVelocity = microVelocities->find( node->first );

            if ( microVelocity == microVelocities->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                      " not found in micro velocity vector" );

            }

            for ( uIntType i = 0; i < nMicroDispDOF; i++ ){

                QDot( nMicroDispDOF * node->second + i ) = microVelocity->second[ i ];

            }

            auto microAcceleration = microAccelerations->find( node->first );

            if ( microAcceleration == microAccelerations->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( node->first ) +
                                      " not found in micro acceleration vector" );

            }

            for ( uIntType i = 0; i < nMicroDispDOF; i++ ){

                QDotDot( nMicroDispDOF * node->second + i ) = microAcceleration->second[ i ];

            }

            for ( uIntType i = 0; i < nMicroDispDOF; i++ ){

                WQ( nMicroDispDOF * node->second + i ) = ( 1 - microWeight->second );

            }

        }

        Eigen::VectorXd D( nMacroDispDOF * macroGlobalToLocalDOFMap->size( ) );
        Eigen::VectorXd DDot( nMacroDispDOF * macroGlobalToLocalDOFMap->size( ) );
        Eigen::VectorXd DDotDot( nMacroDispDOF * macroGlobalToLocalDOFMap->size( ) );

        Eigen::VectorXd WD( nMacroDispDOF * macroGlobalToLocalDOFMap->size( ) );
        for ( auto node = macroGlobalToLocalDOFMap->begin( ); node != macroGlobalToLocalDOFMap->end( ); node++ ){

            auto macroDisplacement = macroDispDOFVector->find( node->first );

            if ( macroDisplacement == macroDispDOFVector->end( ) ){

                return new errorNode( __func__, "Macro node " + std::to_string( node->first ) +
                                      " not found in macro displacement DOF vector" );

            }

            for ( uIntType i = 0; i < nMacroDispDOF; i++ ){

                D( nMacroDispDOF * node->second + i ) = macroDisplacement->second[ i ];

            }

            auto macroVelocity = macroVelocities->find( node->first );

            if ( macroVelocity == macroVelocities->end( ) ){

                return new errorNode( __func__, "Macro node " + std::to_string( node->first ) +
                                      " not found in macro velocity DOF vector" );

            }

            for ( uIntType i = 0; i < nMacroDispDOF; i++ ){

                DDot( nMacroDispDOF * node->second + i ) = macroVelocity->second[ i ];

            }

            auto macroAcceleration = macroAccelerations->find( node->first );

            if ( macroAcceleration == macroAccelerations->end( ) ){

                return new errorNode( __func__, "Macro node " + std::to_string( node->first ) +
                                      " not found in macro acceleration DOF vector" );

            }

            for ( uIntType i = 0; i < nMacroDispDOF; i++ ){

                DDotDot( nMacroDispDOF * node->second + i ) = macroAcceleration->second[ i ];

            }

            auto macroWeight = macroArlequinWeights->find( node->first );

            if ( macroWeight == macroArlequinWeights->end( ) ){

                return new errorNode( __func__, "Macro node " + std::to_string( node->first ) +
                                      " not found in macro Arlequin weight map" );

            }

            for ( uIntType i = 0; i < nMacroDispDOF; i++ ){

                WD( nMacroDispDOF * node->second + i ) = std::fmin(std::fmax(macroWeight->second, 0.001), 0.999);

            }

        }

        // Print statements for debugging
        std::cerr << "Initial micro deformations\n";
        for ( auto node = microGlobalToLocalDOFMap->begin( ); node != microGlobalToLocalDOFMap->end( ); node++ ){

            auto mnd = _inputProcessor.getMicroDisplacements( )->find( node->first );
            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
            for ( unsigned int i = 0; i < 3; i++ ){

                std::cerr << boost::format( "%+1.7f ") % mnd->second[ i ];

            }
            std::cout << "\n";

        }

        std::cout << "\nInitial macro deformations\n";
        for ( auto node = macroGlobalToLocalDOFMap->begin( ); node != macroGlobalToLocalDOFMap->end( ); node++ ){

            auto mnd = _inputProcessor.getMacroDispDOFVector( )->find( node->first );
            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
            for ( unsigned int i = 0; i < 12; i++ ){

                std::cerr << boost::format( "%+1.7f " ) % mnd->second[ i ];

            }
            std::cerr << "\n";

        }

//        // Compute the penalty force
//        Eigen::VectorXd FPenalty = mu * ( _N * D - Q );
//        std::cout << "FPenalty:\n";
//        for ( auto node = microGlobalToLocalDOFMap->begin( ); node != microGlobalToLocalDOFMap->end( ); node++ ){
//            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
//            std::cerr << boost::format( "%1.7f " ) % MQ[ _dim * node->second ];
//            std::cerr << boost::format( "%1.7f " ) % (1 - arlequinMicroWeightingFactors.find( node->first )->second );
//
//            for ( unsigned int j = 0; j < _dim; j++ ){
//
//                    std::cerr << boost::format( "%+1.7f " ) % FPenalty[ nMicroDispDOF * node->second + j ];
//
//            }
//
//            std::cerr << "\n";
//
//        }

        //Compute the diagonal Micromorphic mass matrix ( HRZ homogenization )
        floatType massWeightingFactor = 0;
        floatType totalMass = 0;
        Eigen::VectorXd _MD_Diag = _MD.diagonal( );

//        //HRZ homogenization
//        for ( uIntType n = 0; n < macroGlobalToLocalDOFMap->size( ); n++ ){
//
//            for ( uIntType i = 0; i < _dim; i++ ){
//                massWeightingFactor += _MD_Diag( nMacroDispDOF * n + i);
//                totalMass += _MD.row( nMacroDispDOF * n + i ).sum( );
//            }
//
//        }
//
//        _MD_Diag *= ( totalMass / massWeightingFactor );

        //Row sum homogenization
        for ( auto node = macroGlobalToLocalDOFMap->begin( ); node != macroGlobalToLocalDOFMap->end( ); node++ ){

            auto macroArlequinWeight = macroArlequinWeights->find( node->first );

            if ( macroArlequinWeight == macroArlequinWeights->end( ) ){

                return new errorNode( __func__, "Macro node " + std::to_string( node->first ) +
                                      " not found in macro Arlequin weight vector" );

            }

            for ( uIntType i = 0; i < nMacroDispDOF; i++ ){

                _MD_Diag( nMacroDispDOF * node->second + i )
                    = std::fmin( 0.999, std::fmax( 0.001, macroArlequinWeight->second ) ) *_MD.row( nMacroDispDOF * node->second + i ).sum( );

            }

        }

        // Form the right hand vector for the computation of the Lagrange multiplier values
        Eigen::Map< Eigen::VectorXd > _FD( FD.data( ), FD.size( ) ); 
        Eigen::Map< Eigen::VectorXd > _FQ( FQ.data( ), FQ.size( ) );
        Eigen::Map< Eigen::VectorXd > _De( De.data( ), De.size( ) ); 
        Eigen::Map< Eigen::VectorXd > _Qe( Qe.data( ), Qe.size( ) );

//        // Compute the value of the Lagrangian from the micro-scale equation
//        Eigen::VectorXd _Lagrange = MQ.asDiagonal( ) * QDotDot + aQ * MQ.asDiagonal( ) * QDot - _FQ - mu * ( _N * D - Q );

        // Compute AD and AQ
//        std::cout << "WD:\n" << WD << "\n";
//        std::cout << "WQ:\n" << WQ << "\n";
        Eigen::VectorXd AD = ( ( 1 + aD * gamma * ( *dt ) ) * _MD_Diag + mu_update * beta * ( *dt ) * ( *dt ) * WD ).cwiseInverse( );
        Eigen::VectorXd AQ = ( ( 1 + aQ * gamma * ( *dt ) ) * MQ + mu_update * beta * ( *dt ) * ( *dt ) * WQ ).cwiseInverse( );

//        std::cout << "AD:\n" << AD << "\n";
//        std::cout << "AQ:\n" << AQ << "\n";
//
//
        Eigen::VectorXd RHS = _N * _De - _Qe
                            + beta * ( *dt ) * ( *dt ) * ( _N * AD.asDiagonal( ) * ( _FD + mu_update * WD.asDiagonal( ) * ( D - _De ) ) - AQ.asDiagonal( ) * ( _FQ + mu_update * WQ.asDiagonal( ) * ( Q - _Qe ) ) );

//        std::cout << "RHS:\n" << RHS << "\n";

        // Form the left hand matrix for the computation of the Lagrange multiplier values
//        SparseMatrix LHS = beta * ( *dt ) * ( *dt ) * _N * AD.asDiagonal( ) * _N.transpose();
//        LHS += beta * ( *dt ) * ( *dt ) * AQ.asDiagonal( );
//
//        // Solve for the Lagrange multipliers
//        Eigen::SparseQR< SparseMatrix, Eigen::COLAMDOrdering<int> > lagrangeSolver;
//        lagrangeSolver.compute( LHS );
//        Eigen::VectorXd _Lagrange = lagrangeSolver.solve( RHS );

        // Compute the diagonal LHS matrix
        Eigen::VectorXd LHS = ( beta * ( *dt ) * ( *dt ) * _N * AD.asDiagonal( ) * _N.transpose( ) ) * Eigen::VectorXd::Ones( _N.rows() );
        LHS += beta * ( *dt ) * ( *dt ) * AQ;

        // Solve for the Lagrange multipliers
        Eigen::VectorXd _Lagrange = RHS.cwiseProduct( LHS.cwiseInverse( ) );

        std::cout << "LAGRANGE:\n";
        uIntType _indx = 0;
        std::cerr << "gid, lid: mass      weight    R_1        R_2        R_3        L_1        L_2        L_3\n";
        for ( auto node = microGlobalToLocalDOFMap->begin( ); node != microGlobalToLocalDOFMap->end( ); node++, _indx++ ){

            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
            std::cerr << boost::format( "%1.7f " ) % MQ[ _dim * node->second ];
            std::cerr << boost::format( "%1.7f " ) % (1 - arlequinMicroWeightingFactors.find( node->first )->second );

            auto mnp = _inputProcessor.getMicroNodeReferencePositions( )->find( node->first );
            for ( unsigned int j = 0; j < _dim; j++ ){

                    std::cerr << boost::format( "%+1.7f " ) % mnp->second[ j ];

            }

            for ( unsigned int j = 0; j < _dim; j++ ){

                    std::cerr << boost::format( "%+1.7f " ) % _Lagrange( _dim * node->second + j );

            }
            std::cerr << "\n";

        }

        // Compute the Lagrangian forces on the macro and micro scales
        FALD = floatVector( _N.cols( ), 0 );
        FALQ = floatVector( _Lagrange.size( ), 0 );

        Eigen::Map< Eigen::VectorXd > _FALD( FALD.data( ), FALD.size( ) ); 
        Eigen::Map< Eigen::VectorXd > _FALQ( FALQ.data( ), FALQ.size( ) ); 

        _FALD = -( _N.transpose( ) * _Lagrange ).cwiseProduct( WD.cwiseInverse( ) );
        _FALQ = _Lagrange.cwiseProduct( WQ.cwiseInverse( ) );

        std::cout << "FALQ:\n";
        _indx = 0;
        std::cerr << "gid, lid: mass      weight    R_1        R_2        R_3        L_1        L_2        L_3\n";
        for ( auto node = microGlobalToLocalDOFMap->begin( ); node != microGlobalToLocalDOFMap->end( ); node++, _indx++ ){

            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
            std::cerr << boost::format( "%1.7f " ) % MQ[ _dim * node->second ];
            std::cerr << boost::format( "%1.7f " ) % (1 - arlequinMicroWeightingFactors.find( node->first )->second );

            auto mnp = _inputProcessor.getMicroNodeReferencePositions( )->find( node->first );
            for ( unsigned int j = 0; j < _dim; j++ ){

                    std::cerr << boost::format( "%+1.7f " ) % mnp->second[ j ];

            }

            for ( unsigned int j = 0; j < _dim; j++ ){

                    std::cerr << boost::format( "%+1.7f " ) % _FALQ( _dim * node->second + j );

            }
            std::cerr << "\n";

        }

        std::cout << "FALD:\n";
        _indx = 0;
        std::cerr << "gid, lid: mass      weight    R_1        R_2        R_3        L_1        L_2        L_3\n";
        for ( auto node = macroGlobalToLocalDOFMap->begin( ); node != macroGlobalToLocalDOFMap->end( ); node++, _indx++ ){

            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
            std::cerr << boost::format( "%1.7f " ) % _MD_Diag[ ( _dim + _dim * _dim ) * node->second ];
            std::cerr << boost::format( "%1.7f " ) % WD[ ( _dim + _dim * _dim ) * node->second ];

            auto mnp = _inputProcessor.getMacroNodeReferencePositions( )->find( node->first );
            for ( unsigned int j = 0; j < _dim; j++ ){

                    std::cerr << boost::format( "%+1.7f " ) % mnp->second[ j ];

            }

            for ( unsigned int j = 0; j < ( _dim + _dim * _dim ); j++ ){

                    std::cerr << boost::format( "%+1.7f " ) % _FALD( ( _dim + _dim * _dim ) * node->second + j );

            }
            std::cerr << "\n";

        }

//        std::cerr << "Micro Lagrange\n";
//        std::cerr << _FALQ << "\n";
//
//        std::cerr << "Macro Lagrange\n";
//        std::cerr << _FALD << "\n";

        //TODO: We probably don't need anything left in this function if this works

        // Compute the new accelerations
//        std::cerr << "DDotDot_tp1\n";
        Eigen::VectorXd DDotDot_tp1 = AD.asDiagonal( ) * ( _FD + mu_update * WD.asDiagonal( ) * ( D - _De ) - _N.transpose( ) * _Lagrange );
//        std::cerr << DDotDot_tp1 << "\n";
//        std::cerr << "QDotDot_tp1\n";
        Eigen::VectorXd QDotDot_tp1 = AQ.asDiagonal( ) * ( _FQ + mu_update * WQ.asDiagonal( ) * ( Q - _Qe ) + _Lagrange );
//        std::cerr << QDotDot_tp1 << "\n";

        std::cerr << "\nUpdated micro accelerations\n";
        std::cerr << "gid, lid: R_1        R_2        R_3        U_1        U_2        U_3\n";
        for ( auto node = microGlobalToLocalDOFMap->begin( ); node != microGlobalToLocalDOFMap->end( ); node++ ){

            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
            auto mnp = _inputProcessor.getMicroNodeReferencePositions( )->find( node->first );
            for ( unsigned int j = 0; j < _dim; j++ ){

                    std::cerr << boost::format( "%+1.7f " ) % mnp->second[ j ];

            }

            for ( unsigned int i = 0; i < 3; i++ ){

                std::cerr << boost::format( "%+1.7f ") % QDotDot_tp1[ 3 * node->second + i ];

            }
            std::cout << "\n";

        }

        std::cerr << "\nUpdated macro accelerations\n";
        std::cerr << "gid, gid: R_1        R_2        R_3        U_1        U_2        U_3        PHI_11     PHI_12     PHI_13     PHI_21     PHI_22     PHI_23     PHI_31     PHI_32     PHI_33\n";
        for ( auto node = macroGlobalToLocalDOFMap->begin( ); node != macroGlobalToLocalDOFMap->end( ); node++ ){

            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
            auto mnp = _inputProcessor.getMacroNodeReferencePositions( )->find( node->first );
            for ( unsigned int j = 0; j < _dim; j++ ){

                    std::cerr << boost::format( "%+1.7f " ) % mnp->second[ j ];

            }

            for ( unsigned int i = 0; i < 12; i++ ){

                std::cerr << boost::format( "%+1.7f " ) % DDotDot_tp1[ 12 * node->second + i ];

            }
            std::cerr << "\n";

        }

        // Compute the new deformations
        Eigen::VectorXd _D_tp1 = _De + beta * ( *dt ) * ( *dt ) * DDotDot_tp1;
        Eigen::VectorXd _Q_tp1 = _Qe + beta * ( *dt ) * ( *dt ) * QDotDot_tp1;

//        return new errorNode( "derp", "derp" );
//
//        std::cout << "_MD_Diag:\n" << _MD_Diag << "\n";
//        std::cout << "MQ:\n" << MQ << "\n";
//
//        // Compute the trial acceleration
//        Eigen::Map< Eigen::VectorXd > _FD( FD.data( ), FD.size( ) ); 
//        Eigen::Map< Eigen::VectorXd > _FQ( FQ.data( ), FQ.size( ) );
//        std::cout << "FQ:\n";
//        for ( auto node = microGlobalToLocalDOFMap->begin( ); node != microGlobalToLocalDOFMap->end( ); node++ ){
//
//            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
//            std::cerr << boost::format( "%1.7f " ) % MQ[ _dim * node->second ];
//            std::cerr << boost::format( "%1.7f " ) % (1 - arlequinMicroWeightingFactors.find( node->first )->second );
//
//            for ( unsigned int j = 0; j < _dim; j++ ){
//
//                    std::cerr << boost::format( "%+1.7f " ) % FQ[ nMicroDispDOF * node->second + j ];
//
//            }
//
//            std::cerr << "\n";
//        }
//
//        std::cout << "FD:\n";
//        for ( auto node = macroGlobalToLocalDOFMap->begin( ); node != macroGlobalToLocalDOFMap->end( ); node++ ){
//            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
//
//            for ( unsigned int j = 0; j < nMacroDispDOF; j++ ){
//
//                    std::cerr << boost::format( "%+1.7f " ) % FD[ nMacroDispDOF * node->second + j ];
//
//            }
//
//            std::cerr << "\n";
//        }
//
//        std::cout << "Computing trial accelerations\n";
//        std::cout << "QDotDotTrial:\n";
//        Eigen::VectorXd QDotDotTrial = MQ.cwiseInverse( ).cwiseProduct( _FQ + FPenalty ) / ( 1 + aQ * gamma * ( *dt ) );
//        Eigen::VectorXd DDotDotTrial
//            = _MD_Diag.cwiseInverse( ).cwiseProduct( _FD - _N.transpose( ) * FPenalty ) / ( 1 + aD * gamma * ( *dt ) );
//
//        for ( auto node = microGlobalToLocalDOFMap->begin( ); node != microGlobalToLocalDOFMap->end( ); node++ ){
//            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
//            std::cerr << boost::format( "%1.7f " ) % MQ[ _dim * node->second ];
//            std::cerr << boost::format( "%1.7f " ) % (1 - arlequinMicroWeightingFactors.find( node->first )->second );
//
//            for ( unsigned int j = 0; j < _dim; j++ ){
//
//                    std::cerr << boost::format( "%+1.7f " ) % QDotDotTrial[ nMicroDispDOF * node->second + j ];
//
//            }
//
//            std::cerr << "\n";
//
//        }
//
//        std::cout << "DDotDotTrial:\n";
//        uIntType _indx = 0;
//        for ( auto node = macroGlobalToLocalDOFMap->begin( ); node != macroGlobalToLocalDOFMap->end( ); node++, _indx++ ){
//            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
//
//            for ( unsigned int j = 0; j < nMacroDispDOF; j++ ){
//
//                    std::cerr << boost::format( "%+1.7f " ) % DDotDotTrial[ nMacroDispDOF * node->second + j ];
//
//            }
//
//            std::cerr << "\n";
//
//        }
//
//        // Compute the trial values of D and Q
//        Eigen::Map< Eigen::VectorXd > _De( De.data( ), De.size( ) ); 
//        Eigen::Map< Eigen::VectorXd > _Qe( Qe.data( ), Qe.size( ) ); 
//
//        Eigen::VectorXd Dtrial = _De + ( ( *dt ) * ( *dt ) * beta ) * DDotDotTrial;
//        Eigen::VectorXd Qtrial = _Qe + ( ( *dt ) * ( *dt ) * beta ) * QDotDotTrial;
//
//        std::cout << "Dtrial:\n";
//        for ( auto node = macroGlobalToLocalDOFMap->begin( ); node != macroGlobalToLocalDOFMap->end( ); node++, _indx++ ){
//            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
//
//            for ( unsigned int j = 0; j < nMacroDispDOF; j++ ){
//
//                    std::cerr << boost::format( "%+1.7f " ) % Dtrial[ nMacroDispDOF * node->second + j ];
//
//            }
//
//            std::cerr << "\n";
//
//        }
//        std::cout << "Qtrial:\n";
//        for ( auto node = microGlobalToLocalDOFMap->begin( ); node != microGlobalToLocalDOFMap->end( ); node++ ){
//            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
//            std::cerr << boost::format( "%1.7f " ) % MQ[ _dim * node->second ];
//            std::cerr << boost::format( "%1.7f " ) % (1 - arlequinMicroWeightingFactors.find( node->first )->second );
//
//            for ( unsigned int j = 0; j < _dim; j++ ){
//
//                    std::cerr << boost::format( "%+1.7f " ) % Qtrial[ nMicroDispDOF * node->second + j ];
//
//            }
//
//            std::cerr << "\n";
//
//        }
//
//        //Construct the A matrix
//        floatType AFactor1 = ( ( *dt ) * ( *dt ) * beta ) / ( 1 + aD * gamma * ( *dt ) );
//        floatType AFactor2 = ( ( *dt ) * ( *dt ) * beta ) / ( 1 + aQ * gamma * ( *dt ) );
//        SparseMatrix _A = _N * ( AFactor1 * _MD_Diag.cwiseInverse( ) ).asDiagonal( ) * _N.transpose( );
//        _A += ( AFactor2 * MQ.cwiseInverse( ) ).asDiagonal( );
//
//        //Solve for the Lagrange multipliers
//        Eigen::SparseQR< SparseMatrix, Eigen::COLAMDOrdering<int> > lagrangeSolver;
//        lagrangeSolver.compute( _A );
//        Eigen::VectorXd _Lagrange = lagrangeSolver.solve( mu * ( _N * Dtrial - Qtrial ) );
//
//        std::cout << "LAGRANGE:\n";
//        _indx = 0;
//        std::cerr << "gid, lid: mass      weight    R_1        R_2        R_3        L_1        L_2        L_3\n";
//        for ( auto node = microGlobalToLocalDOFMap->begin( ); node != microGlobalToLocalDOFMap->end( ); node++, _indx++ ){
//
//            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
//            std::cerr << boost::format( "%1.7f " ) % MQ[ _dim * node->second ];
//            std::cerr << boost::format( "%1.7f " ) % (1 - arlequinMicroWeightingFactors.find( node->first )->second );
//
//            auto mnp = _inputProcessor.getMicroNodeReferencePositions( )->find( node->first );
//            for ( unsigned int j = 0; j < _dim; j++ ){
//
//                    std::cerr << boost::format( "%+1.7f " ) % mnp->second[ j ];
//
//            }
//
//            for ( unsigned int j = 0; j < _dim; j++ ){
//
//                    std::cerr << boost::format( "%+1.7f " ) % _Lagrange( _dim * node->second + j );
//
//            }
//            std::cerr << "\n";
//
//        }
//
//        //Compute the new values of the degrees of freedom
//        Eigen::VectorXd _D_tp1 = Dtrial - ( AFactor1 * _MD_Diag.cwiseInverse( ).asDiagonal( ) ) * ( _N.transpose( ) * _Lagrange );
//        Eigen::VectorXd _Q_tp1 = Qtrial + ( AFactor2 * MQ.cwiseInverse( ).asDiagonal( ) ) * _Lagrange;

//        //Solve for the updated DOF vectors
//        _updatedFreeMicroDispDOFValues
//            = floatVector( _Q_tp1.begin( ), _Q_tp1.begin( ) + nMicroDispDOF * _inputProcessor.getFreeMicroNodeIds( )->size( ) );
//        _updatedFreeMacroDispDOFValues
//            = floatVector( _D_tp1.begin( ), _D_tp1.begin( ) + nMacroDispDOF * _inputProcessor.getFreeMacroNodeIds( )->size( ) );
//
//        _projected_ghost_micro_displacement
//            = floatVector( _Q_tp1.begin( ) + nMicroDispDOF * _inputProcessor.getFreeMicroNodeIds( )->size( ), _Q_tp1.end( ) );
//        _projected_ghost_macro_displacement
//            = floatVector( _D_tp1.begin( ) + nMacroDispDOF * _inputProcessor.getFreeMacroNodeIds( )->size( ), _D_tp1.end( ) );
//
//        // Print statements for debugging
//        std::cerr << "\nUpdated micro deformations\n";
//        floatVector microDef = vectorTools::appendVectors( {_updatedFreeMicroDispDOFValues, _projected_ghost_micro_displacement} );
//        std::cerr << "gid, lid: R_1        R_2        R_3        U_1        U_2        U_3\n";
//        for ( auto node = microGlobalToLocalDOFMap->begin( ); node != microGlobalToLocalDOFMap->end( ); node++ ){
//
//            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
//            auto mnp = _inputProcessor.getMicroNodeReferencePositions( )->find( node->first );
//            for ( unsigned int j = 0; j < _dim; j++ ){
//
//                    std::cerr << boost::format( "%+1.7f " ) % mnp->second[ j ];
//
//            }
//
//            for ( unsigned int i = 0; i < 3; i++ ){
//
//                std::cerr << boost::format( "%+1.7f ") % microDef[ 3 * node->second + i ];
//
//            }
//            std::cout << "\n";
//
//        }
//
//        std::cerr << "\nUpdated macro deformations\n";
//        floatVector macroDef = vectorTools::appendVectors( {_updatedFreeMacroDispDOFValues, _projected_ghost_macro_displacement} );
//        std::cerr << "gid, gid: R_1        R_2        R_3        U_1        U_2        U_3        PHI_11     PHI_12     PHI_13     PHI_21     PHI_22     PHI_23     PHI_31     PHI_32     PHI_33\n";
//        for ( auto node = macroGlobalToLocalDOFMap->begin( ); node != macroGlobalToLocalDOFMap->end( ); node++ ){
//
//            std::cerr << boost::format( "%3i, %3i: ") % node->first % node->second;
//            auto mnp = _inputProcessor.getMacroNodeReferencePositions( )->find( node->first );
//            for ( unsigned int j = 0; j < _dim; j++ ){
//
//                    std::cerr << boost::format( "%+1.7f " ) % mnp->second[ j ];
//
//            }
//
//            for ( unsigned int i = 0; i < 12; i++ ){
//
//                std::cerr << boost::format( "%+1.7f " ) % macroDef[ 12 * node->second + i ];
//
//            }
//            std::cerr << "\n";
//
//        }

        return NULL;
    }

    errorOut overlapCoupling::applyKZProjection( const uIntType &macroIncrement, const uIntType &microIncrement ){
        /*!
         * Apply the projection of Klein and Zimmerman
         * to the current increment
         *
         * :param const uIntType &macroIncrement: The increment number for the macro-scale
         * :param const uIntType &microIncrement: The increment number for the micro-scale
         */

        errorOut error = projectDegreesOfFreedom( );
    
        if ( error ){
    
            errorOut result = new errorNode( __func__, "Error in the projection of the ghost degrees of freedom" );
            result->addNext( error );
            return result;
    
        }

#ifdef TESTACCESS

        _test_initial_projected_ghost_micro_displacement = _projected_ghost_micro_displacement;
        _test_initial_projected_ghost_macro_displacement = _projected_ghost_macro_displacement;

#endif

        //Homogenize the material properties at the micro-scale to the macro-scale
        std::cerr << "HOMOGENIZE THE MICROSCALE\n";
        error = homogenizeMicroScale( microIncrement );
    
        if ( error ){
    
            errorOut result = new errorNode( __func__, "Error in the homogenization of the micro-scale to the macro-scale" );
            result->addNext( error );
            return result;
    
        }

        YAML::Node couplingConfiguration = _inputProcessor.getCouplingInitialization( );
        if ( !couplingConfiguration[ "update_displacement" ].IsScalar( ) ){

            //Assemble the mass matrix for the free micromorphic domians
            std::cerr << "ASSEMBLING THE FREE MICROMORPHIC MASS MATRIX\n";
            error = assembleFreeMicromorphicMassMatrix( );
    
            if ( error ){
    
                errorOut result = new errorNode( __func__, "Error in the assembly of the mass matrix for the free macro domains" );
                result->addNext( error );
                return result;
    
            }

            //Assemble the coupling mass and damping matrices
            std::cerr << "ASSEMBLING THE COUPLING MASS AND DAMPING MATRICES\n";
            error = assembleCouplingMassAndDampingMatrices( );
    
            if ( error ){
    
                errorOut result = new errorNode( __func__, "Error in the construction of the coupling mass and damping matrices" );
                result->addNext( error );
                return result;
    
            }

            //Assemble the coupling force vector
            std::cerr << "ASSEMBLING THE COUPLING FORCE VECTOR\n";
            error = assembleCouplingForceVector( );
    
            if ( error ){
    
                errorOut result = new errorNode( __func__, "Error in the construction of the coupling force vector" );
                result->addNext( error );
                return result;
    
            }
    
            //Solve for the free displacements
            std::cout << "SOLVE FOR THE FREE DISPLACEMENT\n";
            error = solveFreeDisplacement( true );
    
            if ( error ){
    
                errorOut result = new errorNode( __func__, "Error when solving for the free displacements" );
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
            return new errorNode ( __func__, "The coupling initialization configuration is not defined" );
        }

        errorOut error = NULL;
        if ( couplingInitialization[ "type" ].as< std::string >( ).compare( "use_first_increment" ) == 0 ){

            error = setReferenceStateFromIncrement( 0, 0 );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error in setting the initial reference state" );
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

                        return new errorNode( __func__, "Macro node " + std::to_string( n->first ) + " not found in the macro displacements map. Fatal error in the input processor" );

                    }

                    _macroReferencePositions.emplace( n->first, n->second + m->second );

                }

                _microReferencePositions.reserve( _inputProcessor.getMicroNodeReferencePositions( )->size( ) );

                for ( auto n = _inputProcessor.getMicroNodeReferencePositions( )->begin( );
                           n != _inputProcessor.getMicroNodeReferencePositions( )->end( );
                           n++ ){

                    auto m = _inputProcessor.getMicroDisplacements( )->find( n->first );

                    if ( m == _inputProcessor.getMicroDisplacements( )->end( ) ){

                        return new errorNode( __func__, "Micro node " + std::to_string( n->first ) + " not found in the micro displacements map. Fatal error in the input processor" );

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
            return new errorNode( __func__,
                                  "The coupling initialization type '" + couplingInitialization[ "type" ].as< std::string >( )
                                + "' is not recognized" );
        }
        
        if ( error ){
            errorOut result = new errorNode( __func__, "Error in initialization of the coupling" );
            result->addNext( error );
            return result;
        }

        //Save the reference state if required
        if ( _inputProcessor.outputReferenceInformation( ) ){

            std::cerr << "OUTPUTTING REFERENCE INFORMATION\n";
            error = outputReferenceInformation( );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error in the output of the reference information" );
                result->addNext( error );
                return result;

            }

        }

        std::cerr << "COUPLING INITIALIZATION COMPLETE\n";

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
            errorOut result = new errorNode( __func__, "Error in initialization of the input processor" );
            result->addNext( error );
            return result;
        }

        //Get the macro cell ids
        const uIntVector *freeMacroCellIDs = _inputProcessor.getFreeMacroCellIds( );
        const uIntVector *ghostMacroCellIDs = _inputProcessor.getGhostMacroCellIds( );

        //Get the macro domain names
        const stringVector *freeMacroDomainNames = _inputProcessor.getFreeMacroDomainNames( );
        const stringVector *ghostMacroDomainNames = _inputProcessor.getGhostMacroDomainNames( );

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

                    return new errorNode( __func__,
                                          "Macro cell " + std::to_string( *cellID ) + " not found in the macro cell to micro domain map" );

                }

                //Get the macro-node set
                error = _inputProcessor._macroscale->getSubDomainNodes( macroIncrement, ( *freeMacroDomainNames )[ cellIndex ], macroNodes );

                if ( error ){

                    errorOut result = new errorNode( __func__,
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

                    return new errorNode( __func__,
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

                        errorOut result = new errorNode( __func__,
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

                    errorOut result = new errorNode( __func__,
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

                    return new errorNode( __func__,
                                          "Macro cell " + std::to_string( *cellID ) + " not found in the macro cell to micro domain map" );

                }

                //Get the macro-node set
                error = _inputProcessor._macroscale->getSubDomainNodes( macroIncrement, ( *ghostMacroDomainNames )[ cellIndex ], macroNodes );

                if ( error ){

                    errorOut result = new errorNode( __func__,
                                                     "Error in extracting the ghost macro-node set" );
                    result->addNext( error );
                    return result;

                }

                if ( _referenceCellDomainCenterOfMassShapefunctions.find( *cellID ) !=
                     _referenceCellDomainCenterOfMassShapefunctions.end( ) ){

                    return new errorNode( __func__,
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

                        errorOut result = new errorNode( __func__,
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

                    errorOut result = new errorNode( __func__,
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

            errorOut result = new errorNode( __func__,
                                             "Error in the formation of the center of mass to macro node projector" );
            result->addNext( error );
            return result;

        }

        //Form the projectors
        error = formTheProjectors( microIncrement, macroIncrement );

        if ( error ){

            errorOut result = new errorNode( __func__,
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

                errorOut result = new errorNode( __func__,
                                                 "Error in the formation of the L2 projectors" );
                result->addNext( error );
                return result;

            }

        }
        else if ( config[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){

            errorOut error = formDirectProjectionProjectors( microIncrement, macroIncrement );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in the formation of the direct projection projectors" );
                result->addNext( error );
                return result;

            }

        }
        else if ( config[ "projection_type" ].as< std::string >( ).compare( "averaged_l2_projection" ) == 0 ){

            errorOut error = formAveragedL2Projectors( );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in the formation of the averaged L2 projectors" );
                result->addNext( error );
                return result;

            }

        }
        else if ( config[ "projection_type" ].as< std::string >( ).compare( "arlequin" ) == 0 ){

            return NULL;

        }
        else{

            return new errorNode( __func__,
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
        errorOut error = DOFProjection::formMoorePenrosePseudoInverse( NQDhat.toDense( ), _dense_BDhatQ );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in solving for _dense_BDhatQ" );
            result->addNext( error );
            return result;

        }

        _dense_BDhatD = -_dense_BDhatQ * _N.topLeftCorner( nFreeMicroDOF, nFreeMacroDOF );

        _dense_BQhatQ = _N.bottomRightCorner( nGhostMicroDOF, nGhostMacroDOF ) * _dense_BDhatQ;
        _dense_BQhatD = _N.bottomLeftCorner( nGhostMicroDOF, nFreeMacroDOF )
                   + _N.bottomRightCorner( nGhostMicroDOF, nGhostMacroDOF ) * _dense_BDhatD;

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

                errorOut result = new errorNode( __func__, "Error in the formation of the selection matrix" );
                result->addNext( error );
                return result;

            }

            error = DOFProjection::formMacroNodeExpansionMatrix( i, nMacroDOF, *_inputProcessor.getMacroGlobalToLocalDOFMap( ), T );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error in the formation of the expansion matrix" );
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

        //Set the DOF type sizes
        uIntType nFreeMacroDOF  = nMacroDOF * _inputProcessor.getFreeMacroNodeIds( )->size( );        
        uIntType nGhostMacroDOF = nMacroDOF * _inputProcessor.getGhostMacroNodeIds( )->size( );

        uIntType nFreeMicroDOF  = nMicroDOF * _inputProcessor.getFreeMicroNodeIds( )->size( );
        uIntType nGhostMicroDOF = nMicroDOF * _inputProcessor.getGhostMicroNodeIds( )->size( );

        std::cerr << "ASSEMBLING THE PROJECTORS\n";
        floatType sparseFactor
        = 1e-4 * 0.5 * ( std::fabs( microMacroProjector.maxCoeff( ) ) + std::fabs( microMacroProjector.minCoeff( ) ) ); //TODO: Make the 1e-4 settable by the user
    
        //Add the homogenization matrix
        microMacroProjector *= _homogenizationMatrix;

        //Compute the projectors
        std::cerr << "  BDhatQ\n";
    
        _sparse_BDhatQ = microMacroProjector.bottomLeftCorner( nGhostMacroDOF, nFreeMicroDOF ).sparseView( 1, sparseFactor );
    
        std::cerr << "  BQhatQ\n";
        _sparse_BQhatQ = _N.bottomRightCorner( nGhostMicroDOF, nGhostMacroDOF ) * _sparse_BDhatQ;
    
        std::cerr << "  BDhatD\n";
        _sparse_BDhatD = -_sparse_BDhatQ * _N.topLeftCorner( nFreeMicroDOF, nFreeMacroDOF );

        std::cerr << "  BQhatD\n";
        _sparse_BQhatD = _N.bottomLeftCorner( nGhostMicroDOF, nFreeMacroDOF );
        _sparse_BQhatD += _N.bottomRightCorner( nGhostMicroDOF, nGhostMacroDOF ) * _sparse_BDhatD;

        std::cerr << "PROJECTORS ASSEMBLED\n";
        return NULL;
    }

    errorOut overlapCoupling::formDirectProjectionProjectors( const unsigned int &microIncrement, const unsigned int &macroIncrement ){
        /*!
         * Form the projectors if the direct projection is to be used
         *
         * :param const unsigned int &microIncrement: The micro increment at which to form the projectors
         * :param const unsigned int &macroIncrement: The macro increment at which to form the projectors
         */

        return new errorNode( __func__,
                              "This subroutine, and the routines it calls, requires extensive re-working to obtain the expected results. The method is not recommended so this has not been done yet." );

        //Get the ghost macro cell IDs
        const uIntVector *ghostMacroCellIDs = _inputProcessor.getGhostMacroCellIds( );
        const std::unordered_map< uIntType, stringVector > *macroCellToMicroDomainMap = _inputProcessor.getMacroCellToDomainMap( );

        const stringVector *ghostMacroDomainNames = _inputProcessor.getGhostMacroDomainNames( );

        unsigned int cellIndex;

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

                return new errorNode( __func__,
                                      "Macro cell " + std::to_string( *cellID ) + " not found in the macro cell to micro domain map" );

            }

            //Get the macro-node set
            error = _inputProcessor._macroscale->getSubDomainNodes( macroIncrement, ( *ghostMacroDomainNames )[ cellIndex ], macroNodes );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in extracting the ghost macro-node set" );
                result->addNext( error );
                return result;

            }

            //Loop over the free micro-domains
            for ( auto domain  = microDomains->second.begin( ); domain != microDomains->second.end( ); domain++ ){

                error = addDomainContributionToDirectFreeMicroToGhostMacroProjector( cellIndex, *cellID, microIncrement, 
                                                                                     *domain, macroNodes );

                if ( error ){

                    errorOut result = new errorNode( __func__,
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

        _sparse_BDhatD = -_sparse_BDhatQ * _N.topLeftCorner( nFreeMicroDOF, nFreeMacroDOF );

        _sparse_BQhatQ = _N.bottomRightCorner( nGhostMicroDOF, nGhostMacroDOF ) * _sparse_BDhatQ;

        _sparse_BQhatD = _N.bottomLeftCorner( nGhostMicroDOF, nFreeMacroDOF )
                   + _N.bottomRightCorner( nGhostMicroDOF, nGhostMacroDOF ) * _sparse_BDhatD;

        return NULL;
    }

    errorOut overlapCoupling::processDomainReference( const unsigned int &microIncrement,
                                                      const std::string &domainName,
                                                      const unsigned int cellID, uIntVector &macroNodes,
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
         * :param uIntVector &macroNodes: The nodes of the macro domain
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

        //Form the macro element if required
        std::unique_ptr< elib::Element > element;
        element.reset(NULL);

        if ( _inputProcessor.isFiltering( ) && _inputProcessor.useReconstructedVolumeForMassMatrix( ) ){

            std::cerr << "#########################################################\n"
                      << "###                      Warning                      ###\n"
                      << "#########################################################\n"
                      << "# Using the reconstructed volume for the center of mass #\n"
                      << "# can result in drastically inconsistent deformation    #\n"
                      << "# when filtering. Overriding to use direct computation  #\n"
                      << "#########################################################\n";

        }
        else if ( _inputProcessor.useReconstructedVolumeForMassMatrix( ) ){

            errorOut error = buildMacroDomainElement( cellID,
                                                      *_inputProcessor.getMacroNodeReferencePositions( ),
                                                      *_inputProcessor.getMacroDisplacements( ),
                                                      *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                      element );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in the construction of the macro-scale element " + std::to_string( cellID ) );
                result->addNext( error );
                return result;

            }

        }

        errorOut error = processDomainMassData( microIncrement, domainName, referenceMicroDomainMass,
                                                referenceMicroDomainCentersOfMass, referenceMicroDomainMomentsOfInertia,
                                                domainReferenceXiVectors, element );
#ifdef TESTACCESS
        _test_domainMass[ cellID ].emplace( domainName, referenceMicroDomainMass[ domainName ] );
        _test_domainCOM[ cellID ].emplace( domainName, referenceMicroDomainCentersOfMass[ domainName ] );
        _test_domainXi[ cellID ].emplace( domainName, domainReferenceXiVectors );
#endif

        if ( error ){
            
            errorOut result = new errorNode( __func__,
                                             "Error in processing the mass data for the micro domain '" + domainName + "'" );
            result->addNext( error );
            return result;

        }

        //Compute the domain's shape function information
        error = computeDomainShapeFunctionInformation( cellID, domainName, microIncrement,
                                                       referenceMicroDomainCentersOfMass[ domainName ],
                                                       macroNodes,
                                                       domainCenterOfMassShapeFunctionValues,
                                                       domainMicroPositionShapeFunctionValues );

        _referenceCellDomainCenterOfMassShapefunctions[ cellID ].emplace( domainName, domainCenterOfMassShapeFunctionValues );

#ifdef TESTACCESS
        if ( _inputProcessor.getCouplingInitialization( )[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ){
            _test_domainMUP[ cellID ].emplace( domainName, domainMicroPositionShapeFunctionValues );
        }
#endif

        if ( error ){
            errorOut result = new errorNode( __func__,
                                             "Error in computing the shape function values for the domain center of mass or the domainMicroPositionShapeFunctionValues" );
            result->addNext( error );
            return result;
        }

        //Get the domain node ids
        uIntVector domainNodes;
        error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, domainName, domainNodes );

        if ( error ){
            errorOut result = new errorNode( __func__,
                                             "Error in extracting the micro-node set" );
            result->addNext( error );
            return result;
        }

        //Add the domain's contribution to the shape function matrix
        error = addDomainContributionToInterpolationMatrix( domainNodes, macroNodes, domainReferenceXiVectors,
                                                            domainCenterOfMassShapeFunctionValues );

        if ( error ){
            errorOut result = new errorNode( __func__,
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

                errorOut result = new errorNode( __func__,
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
                                                     std::unordered_map< uIntType, floatVector > &domainXiVectors,
                                                     std::unique_ptr< elib::Element > &element ){
        /*!
         * Process a micro-scale domain
         *
         * :param const unsigned int microIncrement: The micro increment to process
         * :param const std::string &domainName: The name of the domain
         * :param domainFloatMap &domainMass: The mass of the domain
         * :param domainFloatVectorMap &domainCenterOfMass: The center of mass of the domain
         * :param domainFloatVectorMap &domainMomentOfInertia: The moment of inertia of the domain
         * :param std::unordered_map< uIntType, floatVector > &domainXiVectors: The Xi vectors of the domain
         * :param std::unique_ptr< elib::Element > &element: The macroscale element. Defualts to NULL;
         */

        //Get the domain's nodes
        uIntVector domainNodes;

        errorOut error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, domainName, domainNodes );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error in getting the nodes from the micro domain '" + domainName + "'" );
            result->addNext( error );
            return result;

        }

        floatType   mass;
        floatVector centerOfMass;

        if ( element ){

            std::shared_ptr< volumeReconstruction::volumeReconstructionBase > reconstructedVolume;

            //Reconstruct the domain
            floatVector microNodePositions;
            error = reconstructDomain( microIncrement, domainName, domainNodes, microNodePositions,
                                       element, reconstructedVolume );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in the reconstruction of the microscale domain" );
                result->addNext( error );
                return result;

            }

            const std::unordered_map< uIntType, floatType > *microDensities = _inputProcessor.getMicroDensities( );
            const std::unordered_map< uIntType, floatType > *microVolumes = _inputProcessor.getMicroVolumes( );
            const std::unordered_map< uIntType, floatVector > *microReferencePositions = _inputProcessor.getMicroNodeReferencePositions( );
            const std::unordered_map< uIntType, floatVector > *microDisplacements = _inputProcessor.getMicroDisplacements( );

            //Compute the domain mass properties
            uIntType dataCountAtPoint = 2 + _dim;
            uIntType index = 0;
            floatVector dataAtMicroPoints( dataCountAtPoint * domainNodes.size( ), 0 );
            for ( auto node = domainNodes.begin( ); node != domainNodes.end( ); node++, index++ ){

                dataAtMicroPoints[ dataCountAtPoint * index + 0 ] = 1.; //Integrate the volume of the domain
    
                auto microDensity = microDensities->find( *node );
                if ( microDensity == microDensities->end( ) ){
                    return new errorNode( __func__,
                                          "Micro node " + std::to_string( *node ) + " was not found in the micro density map" );
                }
    
                dataAtMicroPoints[ dataCountAtPoint * index + 1 ] = microDensity->second; //Integrate the density of the domain
    
                auto microReferencePosition = microReferencePositions->find( *node );
                if ( microReferencePosition == microReferencePositions->end( ) ){
                    return new errorNode( __func__,
                                          "Micro node " + std::to_string( *node ) +
                                          " was not found in the micro reference position map" );
                }
    
                auto microDisplacement = microDisplacements->find( *node );
                if ( microDisplacement == microDisplacements->end( ) ){
                    return new errorNode( __func__,
                                          "Micro node " + std::to_string( *node ) +
                                          " was not found in the micro displacement map" );
                }
    
                //Integrate for the domain's center of mass
                for ( unsigned int i = 0; i < _dim; i++ ){
        
                    dataAtMicroPoints[ dataCountAtPoint * index + 2 + i ] =
                        microDensity->second * ( microReferencePosition->second[ i ] + microDisplacement->second[ i ] );
        
                }

            }

            floatVector integratedValues;
            error = reconstructedVolume->performVolumeIntegration( dataAtMicroPoints, dataCountAtPoint, integratedValues );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error in performing volume integration" );
                result->addNext( error );
                return result;
            }

            mass = integratedValues[ 1 ];
            centerOfMass = floatVector( integratedValues.begin( ) + 2, integratedValues.end( ) ) / mass;

        }
        else{

            //Compute the center of mass of the domain
            error = DOFProjection::computeDomainCenterOfMass( _dim, domainNodes, *_inputProcessor.getMicroVolumes( ),
                                                              *_inputProcessor.getMicroDensities( ),
                                                              *_inputProcessor.getMicroNodeReferencePositions( ),
                                                              *_inputProcessor.getMicroDisplacements( ),
                                                              *_inputProcessor.getMicroWeights( ),
                                                              mass, centerOfMass );
    
            if ( error ){
    
                errorOut result = new errorNode( __func__, "Error in calculation of '" + domainName + "' center of mass" );
                result->addNext( error );
                return result;
    
            }

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
            
            errorOut result = new errorNode( __func__, "Error in calculation of '" + domainName + "' xi vectors" );
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
                                                                     uIntVector &macroNodes,
                                                                     floatVector &domainCenterOfMassShapeFunctionValues,
                                                                     std::unordered_map< uIntType, floatVector > &domainMicroPositionShapeFunctionValues ){
        /*!
         * Compute the shape function values at the required locations
         *
         * :param const unsigned int &cellID: The ID number of the cell which contains the domain
         * :param const std::string &domainName: The name of the nodeset which defines the micro-scale domain
         * :param const unsigned int &microIncrement: The micro-scale increment to analyze
         * :param const floatVector &domainCenterOfMass: The center of mass of the domain
         * :param uIntVector &macroNodes: The IDs of the macro nodes corresponding with the shape functions
         * :param floatVector &domainCenterOfMassShapeFunctionValues: The shapefunction values at the center of mass
         * :param std::unordered_map< uIntType, floatVector > &domainMicroPositionShapeFunctionValues:
         *     The shapefunction values at all of the micro nodes contained within the domain.
         */

        //Compute the shape functions of the domain's center of mass
        errorOut error = computeShapeFunctionsAtPoint( cellID,
                                                       *_inputProcessor.getMacroNodeReferencePositions( ),
                                                       *_inputProcessor.getMacroDisplacements( ),
                                                       *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                       domainCenterOfMass, macroNodes, domainCenterOfMassShapeFunctionValues );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error in the computation of the shape function at the center of mass for a micro domain" );
            result->addNext( error );
            return result;

        }

        //Get the domain's nodes
        uIntVector domainNodes;

        error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, domainName, domainNodes );

        if ( error ){

            errorOut result = new errorNode( __func__,
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

                    return new errorNode( __func__,
                                          "Micro node " + std::to_string( *it ) + " was not found in the reference position map" );

                }

                auto microDisplacement = microDisplacements->find( *it );
                if ( microDisplacement == microDisplacements->end( ) ){

                    return new errorNode( __func__,
                                          "Micro node " + std::to_string( *it ) + " was not found in the displacement map" );

                }

                microNodePositions.emplace( *it, microReferencePosition->second + microDisplacement->second );
    
            }
    
            //Compute the shape function values at the micro positions
            error = computeShapeFunctionsAtPoints( cellID,
                                                   *_inputProcessor.getMacroNodeReferencePositions( ),
                                                   *_inputProcessor.getMacroDisplacements( ),
                                                   *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                   microNodePositions, macroNodes,
                                                   domainMicroPositionShapeFunctionValues );

            if ( error ){
    
                errorOut result = new errorNode( __func__,
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
         * Compute the centers of mass for a micro increment. Also computes the mass of the micro-scale domains
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
            errorOut result = new errorNode( __func__, "Error in initialization of the initial increment" );
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

                errorOut result = new errorNode( __func__, "Error in extraction of the free domain's nodes" );
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

                errorOut result = new errorNode( __func__, "Error in calculation of '" + *name + "' center of mass" );
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

                errorOut result = new errorNode( __func__, "Error in extraction of the ghost domain's nodes" );
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

                errorOut result = new errorNode( __func__, "Error in calculation of '" + *name + "' center of mass" );
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

            return new errorNode( __func__,
                                  "Cell ID " + std::to_string( cellID ) + " was not found in the connectivity map" );

        }
        unsigned int cellType = connectivityCellIndices->second[ 0 ];

        //Get the element name
        auto it = elib::XDMFTypeToElementName.find( cellType );

        if ( it == elib::XDMFTypeToElementName.end( ) ){

            return new errorNode( __func__,
                                  "The cell type " + std::to_string(cellType) + " is not supported" );

        }
        
        //Get the nodes from the file
        elib::vecOfvec nodes( connectivityCellIndices->second.size( ) - 1, elib::vec( _dim, 0 ) );
        uIntVector globalNodeIds( connectivityCellIndices->second.begin( ) + 1,
                                  connectivityCellIndices->second.end( ) );
        for ( auto nodeId = globalNodeIds.begin( ); nodeId != globalNodeIds.end( ); nodeId++ ){

            auto nodeLocation = nodeLocations.find( *nodeId );

            if ( nodeLocation == nodeLocations.end( ) ){

                return new errorNode( __func__, "Node " + std::to_string( *nodeId ) + " was not found in the node locations map" );

            }

            uIntType index = nodeId - globalNodeIds.begin( );

            nodes[ index ] = nodeLocation->second;

        }
        
        //Get the element
        auto qrule = elib::default_qrules.find( it->second );
        if ( qrule == elib::default_qrules.end( ) ){

            return new errorNode( __func__,
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

            return new errorNode( __func__, 
                                  "Cell " + std::to_string( cellID ) + " was not found in the connectivity map" );

        }

        //Get the XDMF cell type
        unsigned int cellType = connectivityCellIndices->second[ 0 ];

        //Get the element name
        auto it = elib::XDMFTypeToElementName.find( cellType );

        if ( it == elib::XDMFTypeToElementName.end( ) ){

            return new errorNode( __func__,
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

                return new errorNode( __func__,
                                      "The node " + std::to_string( *nodeId ) + " was not found in the node reference location map" );

            }

            auto nodeDisplacement = nodeDisplacements.find( *nodeId );

            if ( nodeDisplacement == nodeDisplacements.end( ) ){

                return new errorNode( __func__,
                                      "The node " + std::to_string( *nodeId ) + " was not found in the node displacement map" );

            }

            referenceNodes[ index ] = nodeReferenceLocation->second;
            displacements[ index ]  = nodeDisplacement->second;

        }
        
        //Get the element
        auto qrule = elib::default_qrules.find( it->second );
        if ( qrule == elib::default_qrules.end( ) ){

            return new errorNode( __func__,
                                  "The element type " + it->second + " is not found in the default quadrature rules map" );

        }

        element = elib::build_element_from_string( it->second, globalNodeIds, referenceNodes, qrule->second );
        element->update_node_positions( displacements );

        return NULL;

    }

    errorOut overlapCoupling::computeShapeFunctionsAtPoint( const unsigned int cellID,
                                                            const std::unordered_map< uIntType, floatVector > &nodeLocations,
                                                            const std::unordered_map< uIntType, uIntVector > &connectivity,
                                                            const floatVector &point, uIntVector &macroNodes,
                                                            floatVector &shapeFunctions ){
        /*!
         * Compute the shape functions of a given macro-scale domain at a single point
         *
         * :param const unsigned int &cellID: The cell ID at which to compute the shape functions
         * :param const std::unordered_map< uIntType, floatVector > &nodeLocations: The nodal location map
         * :param const std::unordered_map< uIntType, uIntVector > &connectivity: The connectivity map
         * :param const floatVector &point: The point at which to compute the shape functions
         * :param uIntVector &macroNodes: The macro node ids
         * :param floatVector &shapeFunctions: The shapefunctions at the point
         */

        if ( point.size( ) != _dim ){

            return new errorNode( __func__,
                                  "This function only works for a single point of dimension " + std::to_string( _dim ) );

        }

        std::unordered_map< uIntType, floatVector > pointMap;
        pointMap.emplace( 0, point );

        std::unordered_map< uIntType, floatVector > shapefunctionMap;

        errorOut error = computeShapeFunctionsAtPoints( cellID, nodeLocations, connectivity, pointMap, macroNodes, shapefunctionMap );

        if ( error ){

            errorOut result = new errorNode( __func__,
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
                                                            const floatVector &point, uIntVector &macroNodes,
                                                            floatVector &shapeFunctions ){
        /*!
         * Compute the shape functions of a given macro-scale domain at a single point
         *
         * :param const unsigned int &cellID: The cell ID at which to compute the shape functions
         * :param const std::unordered_map< uIntType, floatVector > &nodeLocations: The nodal location map
         * :param const std::unordered_map< uIntType, floatVector > &nodeLocations: The nodal displacement map
         * :param const std::unordered_map< uIntType, uIntVector > &connectivity: The connectivity map
         * :param const floatVector &point: The point at which to compute the shape functions
         * :param uIntVector &macroNodes: The macro node ids
         * :param floatVector &shapeFunctions: The shapefunctions at the point
         */

        if ( point.size( ) != _dim ){

            return new errorNode( __func__,
                                  "This function only works for a single point of dimension " + std::to_string( _dim ) );

        }

        std::unordered_map< uIntType, floatVector > pointMap;
        pointMap.emplace( 0, point );

        std::unordered_map< uIntType, floatVector > shapefunctionMap;

        errorOut error = computeShapeFunctionsAtPoints( cellID, nodeReferenceLocations, nodeDisplacements, 
                                                        connectivity, pointMap, macroNodes, shapefunctionMap );

        if ( error ){

            errorOut result = new errorNode( __func__,
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
                                                             uIntVector &macroNodes,
                                                             std::unordered_map< uIntType, floatVector > &shapeFunctions ){
        /*!
         * Compute the shape functions of a given macro-scale domain at the given points
         *
         * :param const unsigned int &cellID: The cell ID at which to compute the shape functions
         * :param const std::unordered_map< uIntType, floatVector > &nodeLocations: The nodal location map
         * :param const std::unordered_map< uIntType, uIntVector > &connectivity: The connectivity map
         * :param const std::unordered_map< uIntType, floatVector > &points: The points at which to compute the shape functions
         * :param uIntVector &macroNodes: The macro node ids
         * :param std::unordered_map< uIntType, floatVector > &shapeFunctions: The shapefunctions at the points
         */

        //Build the element representing the macro-scale domain
        std::unique_ptr< elib::Element > element;
        errorOut error = overlapCoupling::buildMacroDomainElement( cellID, nodeLocations,
                                                                   connectivity, element );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in element formulation" );
            result->addNext( error );
            return result;

        }

        macroNodes = element->global_node_ids;

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

                return new errorNode( __func__,
                                      "Error in computing the local coordinates for point " + std::to_string( p->first ) );

            }

            error = element->get_shape_functions( localPosition, pointShapeFunctions );

            if ( error ) {

                return new errorNode( __func__,
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
                                                             uIntVector &macroNodes,
                                                             std::unordered_map< uIntType, floatVector > &shapeFunctions ){
        /*!
         * Compute the shape functions of a given macro-scale domain at the given points
         *
         * :param const unsigned int &cellID: The cell ID at which to compute the shape functions
         * :param const std::unordered_map< uIntType, floatVector > &nodeReferenceLocations: The nodal reference location map
         * :param const std::unordered_map< uIntType, floatVector > &nodeDisplacements: The nodal displacement map
         * :param const std::unordered_map< uIntType, uIntVector > &connectivity: The connectivity map
         * :param const std::unordered_map< uIntType, floatVector > &points: The points at which to compute the shape functions
         * :param uIntVector &macroNodes: The macro node ids
         * :param std::unordered_map< uIntType, floatVector > &shapeFunctions: The shapefunctions at the points
         */

        //Build the element representing the macro-scale domain
        std::unique_ptr< elib::Element > element;
        errorOut error = overlapCoupling::buildMacroDomainElement( cellID, nodeReferenceLocations, nodeDisplacements,
                                                                   connectivity, element );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in element formulation" );
            result->addNext( error );
            return result;

        }

        macroNodes = element->global_node_ids;

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

                return new errorNode( __func__,
                                      "Error in computing the local coordinates for point " + std::to_string( p->first ) );

            }

            error = element->get_shape_functions( localPosition, pointShapeFunctions );

            if ( error ) {

                return new errorNode( __func__,
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
                                                       uIntVector &macroNodes,
                                                       std::unordered_map< uIntType, floatVector > &shapeFunctionGradients ){
        /*!
         * Compute the shape functions of a given macro-scale domain at the given points
         *
         * :param const unsigned int &cellID: The cell ID at which to compute the shape functions
         * :param const std::unordered_map< uIntType, floatVector > &nodeReferenceLocations: The nodal reference location vector
         * :param const std::unordered_map< uIntType, floatVector > &nodeDisplacements: The nodal reference location vector
         * :param const std::unordered_map< uIntType, uIntVector > &connectivity: The connectivity vector
         * :param const std::unordered_map< uIntType, floatVector > &points: The points at which to compute the shape functions
         * :param uIntVector &macroNodes: The macro node ids
         * :param std::unordered_map< uIntType, floatVector > &shapeFunctionGradients: The gradients of the shapefunctions
         *     at the points. The output vector is organized [ N11_x, N11_y, N11_z, N12_x, ... ] where the first index is the
         *     point number, the second index is the node number, and the _ indicates a gradient w.r.t. the
         *     third index.
         */

        //Build the element representing the macro-scale domain
        std::unique_ptr< elib::Element > element;
        errorOut error = overlapCoupling::buildMacroDomainElement( cellID, nodeReferenceLocations, nodeDisplacements,
                                                                   connectivity, element );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in element formulation" );
            result->addNext( error );
            return result;

        }

        macroNodes = element->global_node_ids;

        //Make sure the number of points and the size of the output are consistent
        unsigned nPoints = points.size( ) / _dim;
        if ( ( points.size( ) % _dim ) > 0 ){

            return new errorNode( __func__,
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

                return new errorNode( __func__,
                                      "Error in computing the local coordinates for point " + std::to_string( p->first ) );

            }

            error = element->get_global_shapefunction_gradients( localPosition, dNdx );

            if ( error ) {

                return new errorNode( __func__,
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

        _referenceGhostMicroDomainCenterOfMassShapeFunctions.clear( );

        for ( auto cellID = freeMacroCellIds->begin( ); cellID != freeMacroCellIds->end( ); cellID++ ){

            //Get the centers of mass of the macro-domain
            domainCOMs.clear( );
            auto cellCentersOfMass = _referenceGhostMicroDomainCentersOfMass.find( *cellID );

            if ( cellCentersOfMass == _referenceGhostMicroDomainCentersOfMass.end( ) ){

                return new errorNode( __func__,
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

            uIntVector macroNodes;
            error = computeShapeFunctionsAtPoints( *cellID,
                                                   *_inputProcessor.getMacroNodeReferencePositions( ),
                                                   *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                   domainCOMs, macroNodes, macroDomainShapeFunctions );

            if ( error ){

                errorOut result = new errorNode( __func__,
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

                return new errorNode( __func__,
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

            uIntVector macroNodes;
            error = computeShapeFunctionsAtPoints( *cellID,
                                                   *_inputProcessor.getMacroNodeReferencePositions( ),
                                                   *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                   domainCOMs, macroNodes, macroDomainShapeFunctions );

            if ( error ){

                errorOut result = new errorNode( __func__,
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

            errorOut result = new errorNode( __func__,
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

                    return new errorNode( __func__,
                                          "Macro node " + std::to_string( *it ) +
                                          " was not found in the macro displacement dof vector map" );

                }

                if ( macroDispDOF->second.size( ) != nMacroDispDOF ){

                    return new errorNode( __func__,
                                          "Macro node " + std::to_string( *it ) +
                                          " does not have a dimensionally consistent number of degrees of freedom" );

                }

                auto map = macroGlobalToLocalDOFMap->find( *it );

                if ( map == macroGlobalToLocalDOFMap->end( ) ){

                    return new errorNode( __func__,
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

                    return new errorNode( __func__,
                                          "Micro node " + std::to_string( *it ) +
                                          " was not found in the micro displacement dof vector map" );

                }

                if ( microDispDOF->second.size( ) != nMicroDispDOF ){

                    return new errorNode( __func__,
                                          "Micro node " + std::to_string( *it ) +
                                          " does not have a dimensionally consistent number of degrees of freedom" );

                }

                auto map = microGlobalToLocalDOFMap->find( *it );

                if ( map == microGlobalToLocalDOFMap->end( ) ){

                    return new errorNode( __func__,
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

        if ( ( config[ "projection_type" ].as< std::string >( ).compare( "l2_projection" ) == 0 ) ){

            Dhat = _dense_BDhatQ * Q + _dense_BDhatD * D;
            Qhat = _dense_BQhatQ * Q + _dense_BQhatD * D;

        }
        else if ( ( config[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ) ||
                  ( config[ "projection_type" ].as< std::string >( ).compare( "averaged_l2_projection" ) == 0 ) ){

            Dhat = _sparse_BDhatQ * Q + _sparse_BDhatD * D;
            Qhat = _sparse_BQhatQ * Q + _sparse_BQhatD * D;

        }
        else{

            return new errorNode( __func__,
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
        unsigned int n;

        floatType microMass, weight, sf;
        floatVector Xi;

        //Loop through the micro nodes
        for ( auto microNode  = domainNodes.begin( );
                   microNode != domainNodes.end( );
                   microNode++ ){

            auto microDensity = microDensities->find( *microNode );

            if ( microDensity == microDensities->end( ) ){

                return new errorNode( __func__,
                                      "Micro node " + std::to_string( *microNode ) + " was not found in the micro density map" );

            }

            auto microVolume = microVolumes->find( *microNode );

            if ( microVolume == microVolumes->end( ) ){

                return new errorNode( __func__,
                                      "Micro node " + std::to_string( *microNode ) + " was not found in the micro volume map" );

            }

            auto microWeight = microWeights->find( *microNode );

            if ( microWeight == microWeights->end( ) ){

                return new errorNode( __func__,
                                      "Micro node " + std::to_string( *microNode ) + " was not found in the micro weight map" );

            }

            auto referenceXi = domainReferenceXiVectors.find( *microNode );

            if ( referenceXi == domainReferenceXiVectors.end( ) ){

                return new errorNode( __func__,
                                      "Micro node " + std::to_string( *microNode ) + " was not found in the reference Xi map" );

            }

            auto shapefunctions = domainMicroPositionShapeFunctionValues.find( *microNode );

            if ( shapefunctions == domainMicroPositionShapeFunctionValues.end( ) ){

                return new errorNode( __func__,
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

                    return new errorNode( __func__,
                                          "Macro node '" + std::to_string( *macroNode ) + "' not found in global to local macro node map" );

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
                                                                                           uIntVector &macroNodes ){
        /*!
         * Compute the current domain's contribution to the direct free micro to ghost macro
         * projection matrix
         */

        //Compute the shape functions at the micro-nodes for the domain

        //Get the domain node ids
        uIntVector domainNodes;
        errorOut error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, domainName, domainNodes );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error in extracting the domain ( " + domainName + " ) nodes" );
            result->addNext( error );
            return result;

        }

        //Get the micro-node positions and relative position vectors
        std::unordered_map< uIntType, floatVector > microNodePositions;
 
        const std::unordered_map< uIntType, floatVector > *microReferencePositions = _inputProcessor.getMicroNodeReferencePositions( );
        const std::unordered_map< uIntType, floatVector > *microDisplacements      = _inputProcessor.getMicroDisplacements( );
        std::unordered_map< uIntType, floatVector > domainReferenceXiVectors;

        auto cellDomainCentersOfMass = _referenceFreeMicroDomainCentersOfMass.find( cellID );
        if ( cellDomainCentersOfMass == _referenceFreeMicroDomainCentersOfMass.end( ) ){

            return new errorNode( __func__,
                                  "Macro cell " + std::to_string( cellID ) + " not found in reference domain centers of mass map" );

        }

        auto domainCenterOfMass = cellDomainCentersOfMass->second.find( domainName );

        if ( domainCenterOfMass == cellDomainCentersOfMass->second.end( ) ){

            std::string outstr = "Micro domain ";
            outstr += domainName;
            outstr += " not found in the micro domain centers of mass";
            return new errorNode( __func__, outstr );

        }
 
        for ( auto it = domainNodes.begin( ); it != domainNodes.end( ); it++ ){

            auto microReferencePosition = microReferencePositions->find( *it );

            if ( microReferencePosition == microReferencePositions->end( ) ){

                return new errorNode( __func__,
                                      "Micro node " + std::to_string( *it ) + " was not found in the micro reference position map" );

            }

            auto microDisplacement = microDisplacements->find( *it );

            if ( microDisplacement == microDisplacements->end( ) ){

                return new errorNode( __func__,
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
                                               microNodePositions, macroNodes,
                                               domainMicroPositionShapeFunctionValues );

        if ( error ){
 
            errorOut result = new errorNode( __func__,
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

                return new errorNode( __func__,
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
            errorOut result = new errorNode( __func__,
                                             "Error in the computation of the domain's contribution to the micro to macro projection matrix" );
            result->addNext( error );
            return result;
        }

        if ( _sparse_BQhatQ.nonZeros( ) == 0 ){

            _sparse_BDhatQ = domainProjector;

        }
        else{

            _sparse_BDhatQ += domainProjector;

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
        cellDomainMacroSurfaces.clear( );
        homogenizedSurfaceRegionAreas.clear( );
//        homogenizedSurfaceRegionDensities.clear( );
        homogenizedSurfaceRegionCentersOfMass.clear( );
        homogenizedSurfaceRegionProjectedLocalCentersOfMass.clear( );
        homogenizedSurfaceRegionProjectedCentersOfMass.clear( );
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
        errorOut error = NULL;

        uIntVector microDomainNodeIds;
        floatVector microNodePositions;
        std::shared_ptr< volumeReconstruction::volumeReconstructionBase > reconstructedVolume;

        const std::unordered_map< uIntType, stringVector > *macroCellToMicroDomainMap = _inputProcessor.getMacroCellToDomainMap( );
        const std::unordered_map< std::string, uIntType > *microDomainSurfaceSplitCount = _inputProcessor.getMicroDomainSurfaceApproximateSplitCount( );

        const std::unordered_map< uIntType, floatVector > *macroNodeReferenceLocations
            = _inputProcessor.getMacroNodeReferencePositions( );
        const std::unordered_map< uIntType, floatVector > *macroDisplacements
            = _inputProcessor.getMacroDisplacements( );
        const std::unordered_map< uIntType, uIntVector > *macroConnectivity
            = _inputProcessor.getMacroNodeReferenceConnectivity( );

        std::unique_ptr< elib::Element > element;

        std::cerr << "HOMOGENIZING THE FREE MACRO CELLS\n";
        for ( auto macroCell  = _inputProcessor.getFreeMacroCellIds( )->begin( );
                   macroCell != _inputProcessor.getFreeMacroCellIds( )->end( );
                   macroCell++ ){

            //Build the macro element
            error = overlapCoupling::buildMacroDomainElement( *macroCell, *macroNodeReferenceLocations,
                                                              *macroDisplacements, *macroConnectivity,
                                                              element );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error constructing the macro element" );
                result->addNext( error );
                return result;

            }

            //Get the micro domain names within this cell
            auto microDomains = macroCellToMicroDomainMap->find( *macroCell );
            if ( microDomains == macroCellToMicroDomainMap->end( ) ){

                return new errorNode( __func__,
                                      "Macro cell " + std::to_string( *macroCell ) + " not found in the macro cell to micro domain map" ) ;

            }

            for ( auto microDomain  = microDomains->second.begin( ); microDomain != microDomains->second.end( ); microDomain++ ){

                microNodePositions.clear( );
                reconstructedVolume.reset( );

                //Reconstruct the micro-domain's volume
                std::cerr << "  RECONSTRUCTING THE DOMAIN\n";
                error = reconstructDomain( microIncrement, *microDomain, microDomainNodeIds, microNodePositions,
                                           element, reconstructedVolume );

                if ( error ){

                    errorOut result = new errorNode( __func__,
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
                    return new errorNode( __func__, outstr );

                }

                std::cerr << "  COMPUTING THE VOLUME AVERAGES\n";
                error = computeDomainVolumeAverages( *macroCell, *microDomain, microDomainNodeIds,
                                                     reconstructedVolume, &domainCenterOfMass->second );

                if ( error ){

                    errorOut result = new errorNode( __func__,
                                                     "Error in the computation of the volume averages of the microscale domain" );
                    result->addNext( error );
                    return result;

                }

                auto domainSurfaceCount = microDomainSurfaceSplitCount->find( *microDomain );
                if ( domainSurfaceCount == microDomainSurfaceSplitCount->end( ) ){

                    return new errorNode( __func__,
                                          "The micro domain " + *microDomain + " was not found in the domain surface split count map" );

                }
                
                //Compute the surface averages
                std::cerr << "  COMPUTING THE SURFACE AVERAGES\n";
                error = computeDomainSurfaceAverages( *macroCell, *microDomain, microDomainNodeIds,
                                                      domainSurfaceCount->second,
                                                      reconstructedVolume, element );

                if ( error ){

                    errorOut result = new errorNode( __func__,
                                                     "Error in the computation of the surface averages of the microscale domain" );
                    result->addNext( error );
                    return result;

                }
                
            }

            //Compute the approximate stresses
            error = computeHomogenizedStresses( *macroCell );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in the computation of the homogenized stresses" );
                result->addNext( error );
                return result;

            }

        }

        //Loop through the ghost macro-scale cells
        std::cerr << "HOMOGENIZING THE GHOST MACRO CELLS\n";
        for ( auto macroCell  = _inputProcessor.getGhostMacroCellIds( )->begin( );
                   macroCell != _inputProcessor.getGhostMacroCellIds( )->end( );
                   macroCell++ ){

    
            if ( _inputProcessor.isFiltering( ) ){
    
                std::unordered_map< uIntType, floatVector > projectedMacroDisplacements;
    
                auto cellConnectivity = macroConnectivity->find( *macroCell );
    
                if ( cellConnectivity == macroConnectivity->end( ) ){
    
                    return new errorNode( __func__, "Macro cell " + std::to_string( *macroCell ) + " not found in connectivity map" );
    
                }

                auto macroGlobalToLocalDOFMap = _inputProcessor.getMacroGlobalToLocalDOFMap( );

                unsigned int nMacroFreeNodes = _inputProcessor.getFreeMacroNodeIds( )->size( ); // I anticipate this should always be zero

                unsigned int nMacroDOF = _dim + _dim * _dim;

                for ( auto node = cellConnectivity->second.begin( ) + 1; node != cellConnectivity->second.end( ); node++ ){

                    auto local_node = macroGlobalToLocalDOFMap->find( *node );

                    if ( local_node == macroGlobalToLocalDOFMap->end( ) ){

                        return new errorNode( __func__, "Micro node " + std::to_string( *node) + " not found in the global to local map" );

                    }

                    projectedMacroDisplacements.emplace( *node, floatVector( _projected_ghost_macro_displacement.begin( ) + nMacroDOF * local_node->second,
                                                                             _projected_ghost_macro_displacement.begin( ) + nMacroDOF * local_node->second + 3 ) );

                }

                //Build the macro element
                error = overlapCoupling::buildMacroDomainElement( *macroCell, *macroNodeReferenceLocations,
                                                                  projectedMacroDisplacements, *macroConnectivity,
                                                                  element );
 
            }
            else{
    
                //Build the macro element
                error = overlapCoupling::buildMacroDomainElement( *macroCell, *macroNodeReferenceLocations,
                                                                  *macroDisplacements, *macroConnectivity,
                                                                  element );

            }

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error constructing the macro element" );
                result->addNext( error );
                return result;

            }

            //Get the micro domain names within this cell
            auto microDomains = macroCellToMicroDomainMap->find( *macroCell );
            if ( microDomains == macroCellToMicroDomainMap->end( ) ){

                return new errorNode( __func__,
                                      "Macro cell " + std::to_string( *macroCell ) + " not found in the macro cell to micro domain map" ) ;
            }

            for ( auto microDomain = microDomains->second.begin( ); microDomain != microDomains->second.end( ); microDomain++ ){

                microNodePositions.clear( );
                reconstructedVolume.reset( );

                //Reconstruct the micro-domain's volume
                std::cerr << "  RECONSTRUCTING THE DOMAIN\n";
                error = reconstructDomain( microIncrement, *microDomain, microDomainNodeIds, microNodePositions,
                                           element, reconstructedVolume );

                if ( error ){

                    errorOut result = new errorNode( __func__,
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
                    return new errorNode( __func__, outstr );

                }

                std::cerr << "  COMPUTING THE VOLUME AVERAGES\n";
                error = computeDomainVolumeAverages( *macroCell, *microDomain, microDomainNodeIds,
                                                     reconstructedVolume, &domainCenterOfMass->second );

                if ( error ){

                    errorOut result = new errorNode( __func__,
                                                     "Error in the computation of the volume averages of the microscale domain" );
                    result->addNext( error );
                    return result;

                }
                
                auto domainSurfaceCount = microDomainSurfaceSplitCount->find( *microDomain );
                if ( domainSurfaceCount == microDomainSurfaceSplitCount->end( ) ){

                    return new errorNode( __func__,
                                          "The micro domain " + *microDomain + " was not found in the domain surface split count map" );

                }

                //Compute the surface averages
                std::cerr << "  COMPUTING THE SURFACE AVERAGES\n";
                error = computeDomainSurfaceAverages( *macroCell, *microDomain, microDomainNodeIds,
                                                      domainSurfaceCount->second,
                                                      reconstructedVolume, element );

                if ( error ){

                    errorOut result = new errorNode( __func__,
                                                     "Error in the computation of the surface averages of the microscale domain" );
                    result->addNext( error );
                    return result;

                }
                
            }

            //Compute the approximate stresses
            error = computeHomogenizedStresses( *macroCell );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in the computation of the homogenized stresses" );
                result->addNext( error );
                return result;

            }

        }

        //Compute the homogenized force vectors and mass matrices
        const YAML::Node config = _inputProcessor.getCouplingInitialization( );

        if ( config[ "projection_type" ].as< std::string >( ).compare( "arlequin" ) != 0 ){

            std::cerr << "ASSEMBLE HOMOGENIZED MATRICES AND VECTORS\n";
            error = assembleHomogenizedMatricesAndVectors( );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in the computation of the homogenized forces and mass matrix" );
                result->addNext( error );
                return result;

            }

        }

        std::cerr << "MICROSCALE HOMOGENIZATION COMPLETED\n";

        return NULL;
    }

    errorOut overlapCoupling::reconstructDomain( const unsigned int &microIncrement, const std::string &microDomainName,
                                                 uIntVector &microDomainNodes, floatVector &microNodePositions,
                                                 const std::unique_ptr< elib::Element > &element,
                                                 std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume ){
        /*!
         * Reconstruct the micro-domain's volume to perform volume and surface integrals over that
         * domain.
         *
         * :param const unsigned int &microIncrement: The increment at which to extract the micro-positions
         * :param const std::string &microDomainName: The name of the micro-domain to be re-constructed.
         * :param uIntVector &microDomainNodes: The nodes associated with the micro domain
         * :param floatVector &microNodePositions: The positions of the micro nodes for the current domain
         * :param const std::unique_ptr< elib::element > &element: The macro-scale element that contains the domain.
         * :param std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume: The reconstructed
         *     volume ready for additional processing.
         */

        //Get the volume reconstruction configuration
        YAML::Node volumeReconstructionConfig = _inputProcessor.getVolumeReconstructionConfig( );

        //Get the domain node ids
        errorOut error = _inputProcessor._microscale->getSubDomainNodes( microIncrement, microDomainName, microDomainNodes );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error in getting the node ids for the domain ( " + microDomainName + " )" );
            result->addNext( error );
            return result;

        }

        //Get the micro-node positions
        microNodePositions.clear( );
        microNodePositions.reserve( _dim * microDomainNodes.size( ) );
 
        unsigned int index = 0;
        const std::unordered_map< uIntType, floatVector > *microReferencePositions = _inputProcessor.getMicroNodeReferencePositions( );
        const std::unordered_map< uIntType, floatVector > *microDisplacements      = _inputProcessor.getMicroDisplacements( );
        uIntVector interiorNodes;
 
        for ( auto it = microDomainNodes.begin( ); it != microDomainNodes.end( ); it++, index++ ){

            auto microReferencePosition = microReferencePositions->find( *it );

            if ( microReferencePosition == microReferencePositions->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( *it ) +
                                      " was not found in the micro reference position map" );

            }

            auto microDisplacement = microDisplacements->find( *it );

            if ( microDisplacement == microDisplacements->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( *it ) +
                                      " was not found in the micro displacement map" );

            }

            //Check that the given micro-node is located inside of the macro-scale element.
            //If not, we will ignore the micro node
            floatVector cp = microDisplacement->second + microReferencePosition->second;
            if ( element->contains_point( cp, volumeReconstructionConfig[ "element_contain_tolerance" ].as< floatType >( ) ) ){
            
                interiorNodes.push_back( *it );

                for ( unsigned int i = 0; i < _dim; i++ ){
    
                    microNodePositions.push_back( microReferencePosition->second[ i ] + microDisplacement->second[ i ] );
    
                }

            }
 
        }

        //Compute the macro-domain's surface planes in the current configuration
        floatMatrix surfacePoints( element->local_surface_points.size( ), floatVector( _dim, 0 ) );

        floatMatrix surfaceNormals( element->local_surface_normals.size( ), floatVector( _dim, 0 ) );

        for ( unsigned int i = 0; i < element->local_surface_normals.size( ), floatVector( _dim, 0 ) ){

            floatVector lN = element->local_surface_normals[ i ];

            floatVector lP = element->local_surface_points[ i ];

            element->transform_local_vector( lP, lN, surfaceNormals[ i ], true );

            surfaceNormals[ i ] /= vectorTools::l2norm( surfaceNormals[ i ] );

            element->interpolate( element->nodes, lP, surfacePoints[ i ] );

        }

        //Reset micro domain nodes
        microDomainNodes = interiorNodes;

        //Pass the base name of the output file to the volume reconstruction configuration to be used if output has been requested
        volumeReconstructionConfig[ "baseOutputFilename" ] = microDomainName + "_" + std::to_string( microIncrement );

        //Get the volume reconstruction object
        reconstructedVolume
            = volumeReconstruction::volumeReconstructionBase( volumeReconstructionConfig ).create( );

        if ( reconstructedVolume->getError( ) ){

            errorOut result = new errorNode( __func__,
                                             "Error in creating the volume reconstruction object for " + microDomainName );

            result->addNext( reconstructedVolume->getError( ) );

            return result;

        }

        //Load the micro points
        error = reconstructedVolume->loadPoints( &microNodePositions );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error in loading the micro-scale points for " + microDomainName );
            result->addNext( error );

            return result;

        }

        //Add the element's bounding planes
        error = reconstructedVolume->addBoundingPlanes( surfacePoints, surfaceNormals );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error in loading the boundary planes for " + microDomainName );

            result->addNext( error );

            return result;

        }

        //Reconstruct the volume
        error = reconstructedVolume->evaluate( );

        if ( error ){

            errorOut result = new errorNode( __func__,
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
                return new errorNode( __func__,
                                      "Micro node " + std::to_string( *node ) + " was not found in the micro density map" );
            }

            dataAtMicroPoints[ dataCountAtPoint * index + 1 ] = microDensity->second; //Integrate the density of the domain

            //Integrate the micro stresses
            auto microStress = microStresses->find( *node );
            if ( microStress == microStresses->end( ) ){
                return new errorNode( __func__,
                                      "Micro node " + std::to_string( *node ) + " was not found in the micro stress map" );
            }
            for ( unsigned int i = 0; i < _dim * _dim; i++ ){
                dataAtMicroPoints[ dataCountAtPoint * index + 2 + i ] = microStress->second[ i ];
            }

            localIndex = initialOffset;

            if ( _inputProcessor.useReconstructedMassCenters( ) ){

                auto microReferencePosition = microReferencePositions->find( *node );
                if ( microReferencePosition == microReferencePositions->end( ) ){
                    return new errorNode( __func__,
                                          "Micro node " + std::to_string( *node ) + " was not found in the micro reference position map" );
                }

                auto microDisplacement = microDisplacements->find( *node );
                if ( microDisplacement == microDisplacements->end( ) ){
                    return new errorNode( __func__,
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
                    return new errorNode( __func__,
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
                    return new errorNode( __func__,
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

            errorOut result = new errorNode( __func__,
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
                    return new errorNode( __func__,
                                          "Micro node " + std::to_string( *node ) + " was not found in the micro density map" );
                }

                auto microReferencePosition = microReferencePositions->find( *node );
                if ( microReferencePosition == microReferencePositions->end( ) ){
                    return new errorNode( __func__,
                                          "Micro node " + std::to_string( *node ) + " was not found in the micro reference position map" );
                }

                auto microDisplacement = microDisplacements->find( *node );
                if ( microDisplacement == microDisplacements->end( ) ){
                    return new errorNode( __func__,
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
                        return new errorNode( __func__,
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
                        return new errorNode( __func__,
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

                errorOut result = new errorNode( __func__,
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
                                                            std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume,
                                                            std::unique_ptr< elib::Element > &element ){
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
         * :param std::unique_ptr< elib::Element > &element: The enclosing macro element
         */

        //Get the volume reconstruction configuration
        YAML::Node volumeReconstructionConfig = _inputProcessor.getVolumeReconstructionConfig( );

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

            errorOut result = new errorNode( __func__,
                                             "Error in the computation of the domain's surface area" );
            result->addNext( error );
            return result;

        }

        if ( homogenizedSurfaceAreas.find( macroCellID ) == homogenizedSurfaceAreas.end( ) ){

            domainFloatMap tmpFloatMap;
            domainFloatVectorMap tmpFloatVectorMap;
            homogenizedSurfaceAreas.emplace( macroCellID, tmpFloatMap );
            homogenizedSurfaceRegionAreas.emplace( macroCellID, tmpFloatVectorMap );
            homogenizedSurfaceRegionCentersOfMass.emplace( macroCellID, tmpFloatVectorMap );
            homogenizedSurfaceRegionProjectedLocalCentersOfMass.emplace( macroCellID, tmpFloatVectorMap );
            homogenizedSurfaceRegionProjectedCentersOfMass.emplace( macroCellID, tmpFloatVectorMap );
            homogenizedSurfaceRegionTractions.emplace( macroCellID, tmpFloatVectorMap );
            homogenizedSurfaceRegionCouples.emplace( macroCellID, tmpFloatVectorMap );

        }

        homogenizedSurfaceAreas[ macroCellID ].emplace( microDomainName, integratedValue[ 0 ] );

        /*=====================================================================
        |           Compute the properties of the surface subdomains          |
        =====================================================================*/

        //Get the boundary points from the reconstructed domain
        const uIntVector *boundaryIDs = reconstructedVolume->getBoundaryIDs( );
        const floatVector *boundaryPoints = reconstructedVolume->getBoundaryPoints( );

        //Determine which element surface each point belongs to ( if any )
        uIntType bpIndex = 0;
        floatType tol = volumeReconstructionConfig[ "element_contain_tolerance" ].as< floatType >( );
        bool useMacroNormals = volumeReconstructionConfig[ "use_macro_normals" ].as< bool >( );

        uIntVector surf;
        uIntVector tmp( 0 );
        floatVector ftmp( 0 );
        uIntMatrix subdomainNodeIDs( element->local_surface_points.size( ), tmp );
        floatMatrix subdomainNodeNormals( element->local_surface_points.size( ), ftmp );

        auto sIndex = 0;
        for ( auto s = subdomainNodeIDs.begin( ); s != subdomainNodeIDs.end( ); s++, sIndex++ ){

            ( *s ).reserve( boundaryIDs->size( ) / subdomainNodeIDs.size( ) );
            subdomainNodeNormals[ sIndex ].reserve( _dim * boundaryIDs->size( ) / subdomainNodeIDs.size( ) );

        }

        for ( auto id = boundaryIDs->begin( ); id != boundaryIDs->end( ); id++, bpIndex++ ){

            floatVector p( boundaryPoints->begin( ) + _dim * bpIndex,
                           boundaryPoints->begin( ) + _dim * ( bpIndex + 1 ) );

            floatVector xi;

            //Compute the local coordinates of the point
            std::unique_ptr< errorNode > tmpError;
            tmpError.reset( element->compute_local_coordinates( p, xi ) );

            if ( tmpError ){

                continue;

            }

            //Determine if the point is on any of the element's surfaces
            if ( element->local_point_on_surface( xi, surf, tol ) ){

                for ( auto s = surf.begin( ); s != surf.end( ); s++ ){

                    subdomainNodeIDs[ *s ].push_back( *id );

                    floatVector n;
                    element->transform_local_vector( xi, element->local_surface_normals[ *s ], n );

                    n /= vectorTools::l2norm( n );

                    for ( uIntType i = 0; i < _dim; i++ ){

                        subdomainNodeNormals[ *s ].push_back( n[ i ] );

                    }

                }

            }

        }

        uIntVector macroSurfaces;
        for ( uIntType i = 0; i < subdomainNodeIDs.size( ); i++ ){

            if ( subdomainNodeIDs[ i ].size( ) > 0 ){

                macroSurfaces.push_back( i );

            }

        }

        cellDomainMacroSurfaces[ macroCellID ].emplace( microDomainName, macroSurfaces );

        //Get the centers of mass of the surface regions

        dataCountAtPoint = 1     //Surface area of region
                         + 1;    //Mass of region
//                         + _dim; //Micro point position

        dataAtMicroPoints.clear( );
        dataAtMicroPoints.reserve( dataCountAtPoint * microDomainNodeIDs.size( ) );

        for ( auto node = microDomainNodeIDs.begin( ); node != microDomainNodeIDs.end( ); node++ ){

            auto microDensity = microDensities->find( *node );
            if ( microDensity == microDensities->end( ) ){
                return new errorNode( __func__,
                                      "The micro node " + std::to_string( *node ) + " was not found in the micro density map" );
            }

            dataAtMicroPoints.push_back( 1 );
            dataAtMicroPoints.push_back( microDensity->second );

            auto microReferencePosition = microReferencePositions->find( *node );
            if ( microReferencePosition == microReferencePositions->end( ) ){
                return new errorNode( __func__,
                                      "Micro node " + std::to_string( *node ) + " was not found in the micro reference position map" );
            }

//            auto microDisplacement = microDisplacements->find( *node );
//            if ( microDisplacement == microDisplacements->end( ) ){
//                return new errorNode( __func__,
//                                      "Micro node " + std::to_string( *node ) + " was not found in the micro displacement map" );
//            }
//
//            floatVector microPoint = microReferencePosition->second + microDisplacement->second;
//
//            for ( unsigned int i = 0; i < microPoint.size( ); i++ ){
//
//                dataAtMicroPoints.push_back( microDensity->second * microPoint[ i ] );
//
//            }

        }

        //Initialize storage values for homogenization
        homogenizedSurfaceRegionAreas[ macroCellID ].emplace( microDomainName, floatVector( subdomainNodeIDs.size( ), 0 ) );
        floatVector regionDensities( subdomainNodeIDs.size( ) );
        homogenizedSurfaceRegionCentersOfMass[ macroCellID ].emplace( microDomainName, floatVector( _dim * subdomainNodeIDs.size( ), 0 ) );
        homogenizedSurfaceRegionProjectedLocalCentersOfMass[ macroCellID ].emplace( microDomainName, floatVector( _dim * subdomainNodeIDs.size( ), 0 ) );
        homogenizedSurfaceRegionProjectedCentersOfMass[ macroCellID ].emplace( microDomainName, floatVector( _dim * subdomainNodeIDs.size( ), 0 ) );

        //Initialize the homogenized values
        homogenizedSurfaceRegionTractions[ macroCellID ].emplace( microDomainName, floatVector( _dim * subdomainNodeIDs.size( ), 0 ) );
        homogenizedSurfaceRegionCouples[ macroCellID ].emplace( microDomainName, floatVector( _dim * _dim * subdomainNodeIDs.size( ), 0 ) );

        uIntType index = 0;
        for ( auto sN = subdomainNodeIDs.begin( ); sN != subdomainNodeIDs.end( ); sN++, index++ ){

            if ( ( *sN ).size( ) > 0 ){

                // Perform the surface integration
                uIntVector *nodesOnSurface = &( *sN );
                error = reconstructedVolume->performSurfaceIntegration( dataAtMicroPoints, dataCountAtPoint,
                                                                        integratedValue, nodesOnSurface,
                                                                        NULL,
                                                                        &subdomainNodeNormals[ index ] );
    
                if ( error ){
    
                    errorOut result = new errorNode( __func__,
                                                     "Error in the integration of the macro surface ( "
                                                     + std::to_string( sN - subdomainNodeIDs.begin( ) ) + " )" );
                    result->addNext( error );
                    return result;
    
                }

                // Extract the region surface areas and the region surface densities
                homogenizedSurfaceRegionAreas[ macroCellID ][ microDomainName ][ index ] = integratedValue[ 0 ];
                regionDensities[ index ] = integratedValue[ 1 ] / integratedValue[ 0 ];

                // Perform the position weighted surface integration
                error = reconstructedVolume->performPositionWeightedSurfaceIntegration( dataAtMicroPoints, dataCountAtPoint,
                                                                                        integratedValue, nodesOnSurface,
                                                                                        NULL,
                                                                                        &subdomainNodeNormals[ index ] );

                if ( error ){
    
                    errorOut result = new errorNode( __func__,
                                                     "Error in the integration of the position weighted macro surface ( "
                                                     + std::to_string( sN - subdomainNodeIDs.begin( ) ) + " )" );
                    result->addNext( error );
                    return result;
    
                }

                floatVector centerOfMass( integratedValue.begin( ) + _dim,
                                          integratedValue.begin( ) + 2 * _dim );
                centerOfMass /= ( homogenizedSurfaceRegionAreas[ macroCellID ][ microDomainName ][ index ] * regionDensities[ index ] );

                for ( uIntType i = 0; i < _dim; i++ ){

                    homogenizedSurfaceRegionCentersOfMass[ macroCellID ][ microDomainName ][ index * _dim + i ] = centerOfMass[ i ];

                }

                //Compute the local coordinates of the center of mass
                floatVector localCenterOfMass; 
                error = element->compute_local_coordinates( centerOfMass, localCenterOfMass );

                if ( error ){

                    errorOut result = new errorNode( __func__,
                                                     "Error in the computation of the local surface center of mass" );
                    result->addNext( error );
                    return result;

                }

                //Find the closest point on the surface of the element in local coordinates
                floatType distance = vectorTools::dot( element->local_surface_normals[ index ],
                                                       localCenterOfMass - element->local_surface_points[ index ] );

                floatVector projectedLocalCenterOfMass = localCenterOfMass - distance * element->local_surface_normals[ index ];

                for ( unsigned int i = 0; i < _dim; i++ ){
                    homogenizedSurfaceRegionProjectedLocalCentersOfMass[ macroCellID ][ microDomainName ][ index * _dim + i ]
                        = projectedLocalCenterOfMass[ i ];
                }

                //Interpolate this point back to global coordinates
                floatVector projectedCenterOfMass;
                error = element->interpolate( element->nodes, projectedLocalCenterOfMass, projectedCenterOfMass );

                if ( error ){

                    errorOut result = new errorNode( __func__,
                                                     "Error in the interpolation of the local projected surface center of mass to global coordinates" );
                    result->addNext( error );
                    return result;

                }

                for ( unsigned int i = 0; i < _dim; i++ ){
                    homogenizedSurfaceRegionProjectedCentersOfMass[ macroCellID ][ microDomainName ][ index * _dim + i ]
                        = projectedCenterOfMass[ i ];
                }

            }

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
                return new errorNode( __func__,
                                      "Micro node " + std::to_string( *node ) + " was not found in the micro stress map" );
            }

            for ( unsigned int i = 0; i < _dim * _dim; i++ ){

                dataAtMicroPoints.push_back( microStress->second[ i ] );

            }

        }

        index = 0;
        for ( auto sN = subdomainNodeIDs.begin( ); sN != subdomainNodeIDs.end( ); sN++, index++ ){

            if ( ( *sN ).size( ) > 0 ){

                //Compute the tractions
                uIntVector *nodesOnSurface = &( *sN );
                error = reconstructedVolume->performSurfaceFluxIntegration( dataAtMicroPoints, dataCountAtPoint,
                                                                            integratedValue, nodesOnSurface,
                                                                            NULL,
                                                                            &subdomainNodeNormals[ index ], useMacroNormals );

                if ( error ){
    
                    errorOut result = new errorNode( __func__,
                                                     "Error in the computation of the surface traction of the the macro surface ( "
                                                     + std::to_string( sN - subdomainNodeIDs.begin( ) ) + " )" );
                    result->addNext( error );
                    return result;
    
                }

                for ( uIntType i = 0; i < _dim; i++ ){
                    homogenizedSurfaceRegionTractions[ macroCellID ][ microDomainName ][ _dim * index + i ]
                        = integratedValue[ i ] / homogenizedSurfaceRegionAreas[ macroCellID ][ microDomainName ][ index ];
                }

                //Compute the couples
                floatVector regionCenterOfMass( homogenizedSurfaceRegionCentersOfMass[ macroCellID ][ microDomainName ].begin( ) + _dim * index,
                                                homogenizedSurfaceRegionCentersOfMass[ macroCellID ][ microDomainName ].begin( ) + _dim * ( index + 1 ) );

                error = reconstructedVolume->performRelativePositionSurfaceFluxIntegration( dataAtMicroPoints, dataCountAtPoint,
                                                                                            regionCenterOfMass, integratedValue,
                                                                                            nodesOnSurface,
                                                                                            NULL,
                                                                                            &subdomainNodeNormals[ index ], useMacroNormals );
    
                if ( error ){
    
                    errorOut result = new errorNode( __func__,
                                                     "Error in the computation of the surface couple of the micro region ( "
                                                     + std::to_string( sN - subdomainNodeIDs.begin( ) ) + " )" );
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

            }

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
        errorOut error;

        std::unordered_map< uIntType, floatVector > projectedMacroDisplacements;
    
        if ( _inputProcessor.isFiltering( ) ){
    
            auto cellConnectivity = macroConnectivity->find( macroCellID );
    
            if ( cellConnectivity == macroConnectivity->end( ) ){
    
                return new errorNode( __func__, "Macro cell " + std::to_string( macroCellID ) + " not found in connectivity map" );
    
            }

            auto macroGlobalToLocalDOFMap = _inputProcessor.getMacroGlobalToLocalDOFMap( );

            unsigned int nMacroFreeNodes = _inputProcessor.getFreeMacroNodeIds( )->size( ); // I anticipate this should always be zero

            unsigned int nMacroDOF = _dim + _dim * _dim;

            for ( auto node = cellConnectivity->second.begin( ) + 1; node != cellConnectivity->second.end( ); node++ ){

                auto local_node = macroGlobalToLocalDOFMap->find( *node );

                if ( local_node == macroGlobalToLocalDOFMap->end( ) ){

                    return new errorNode( __func__, "Micro node " + std::to_string( *node) + " not found in the global to local map" );

                }

                projectedMacroDisplacements.emplace( *node, floatVector( _projected_ghost_macro_displacement.begin( ) + nMacroDOF * local_node->second,
                                                                         _projected_ghost_macro_displacement.begin( ) + nMacroDOF * local_node->second + 3 ) );

            }

            //Build the macro element
            error = overlapCoupling::buildMacroDomainElement( macroCellID, *macroNodeReferenceLocations,
                                                              projectedMacroDisplacements, *macroConnectivity,
                                                              element );
 
        }
        else{
    
            //Build the macro element
            error = overlapCoupling::buildMacroDomainElement(  macroCellID, *macroNodeReferenceLocations,
                                                              *macroDisplacements, *macroConnectivity,
                                                              element );

        }

        if ( error ){

            errorOut result = new errorNode( __func__,
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

        uIntVector macroNodes;

        if ( _inputProcessor.isFiltering( ) ){

            error = overlapCoupling::computeShapeFunctionsAtPoints( macroCellID, *macroNodeReferenceLocations,
                                                                    projectedMacroDisplacements, *macroConnectivity,
                                                                    centerOfMassMap, macroNodes, shapefunctionsAtCentersOfMassByID );

        }
        else{

            error = overlapCoupling::computeShapeFunctionsAtPoints( macroCellID, *macroNodeReferenceLocations, *macroDisplacements,
                                                                    *macroConnectivity,
                                                                    centerOfMassMap, macroNodes,
                                                                    shapefunctionsAtCentersOfMassByID );

        }

        if ( error ){

            errorOut result = new errorNode( __func__,
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
    
                return new errorNode( __func__, output );

            }

        }

        //Compute element nodal volumes
        floatVector elementNodalVolumes( element->nodes.size( ), 0 );
        for ( auto qpt = element->qrule.begin( ); qpt != element->qrule.end( ); qpt++ ){

            //Get the values of the shape function and the gradients
            floatVector shapeFunctions;
            error = element->get_shape_functions( qpt->first, shapeFunctions );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in the computation of the shape functions\n" );
                result->addNext( error );
                return result;

            }

            //Get the Jacobian of transformation
            elib::vecOfvec jacobian;
            error = element->get_local_gradient( element->nodes, qpt->first, jacobian );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in the computation of the local gradient\n" );
                result->addNext( error );
                return result;

            }

            floatType Jxw = vectorTools::determinant( vectorTools::appendVectors( jacobian ), _dim, _dim ) * qpt->second;

            for ( uIntType n = 0; n < shapeFunctions.size( ); n++ ){

                elementNodalVolumes[ n ] += shapeFunctions[ n ] * Jxw;

            }

        }

        //Compute element nodal surface areas
        floatMatrix elementNodalSurfaceAreas( element->nodes.size( ), floatVector( element->local_surface_points.size( ), 0 ) );
        uIntType surfaceNum = 0;
        for ( auto element_surface = element->surface_quadrature_rules.begin( ); element_surface != element->surface_quadrature_rules.end( ); element_surface++, surfaceNum++ ){

            for ( auto s_qpt = element_surface->begin( ); s_qpt != element_surface->end( ); s_qpt++ ){

                floatMatrix surface_jacobian;
                error = element->get_local_gradient( element->nodes, s_qpt->first, surface_jacobian );

                if ( error ){

                    std::string message = "Error in computing the local gradient of the shape functions for quadrature point " + std::to_string( ( unsigned int )( s_qpt - element_surface->begin( ) ) )
                                        + " on surface " + std::to_string( surfaceNum ) + " for element " + std::to_string( macroCellID );
                    return new errorNode( __func__, message );

                }

                floatVector inv_jacobian = vectorTools::inverse( vectorTools::appendVectors( surface_jacobian ), _dim, _dim );
                floatType Jxw = vectorTools::determinant( vectorTools::appendVectors( surface_jacobian ), _dim, _dim ) * s_qpt->second;
                floatType surfaceAreaContribution = vectorTools::l2norm( vectorTools::Tdot( vectorTools::inflate( inv_jacobian, _dim, _dim ), element->local_surface_normals[ surfaceNum ] ) * Jxw );

                floatVector shapeFunctions;
                element->get_shape_functions( s_qpt->first, shapeFunctions );

                for ( uIntType n = 0; n < shapeFunctions.size( ); n++ ){
                
                    elementNodalSurfaceAreas[ n ][ surfaceNum ] += shapeFunctions[ n ] * surfaceAreaContribution;
                
                }
                
            }

        }

#ifdef TESTACCESS

        _test_elementNodalVolumes.emplace( macroCellID, elementNodalVolumes );

#endif

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

        std::vector< floatVector > areaAtNodes( nMacroCellNodes, floatVector( element->local_surface_normals.size( ), 0 ) );
        std::vector< floatMatrix > tractionAtNodes( nMacroCellNodes, floatMatrix( element->local_surface_normals.size( ), floatVector( _dim, 0 ) ) );
        std::vector< floatMatrix > coupleAtNodes( nMacroCellNodes, floatMatrix( element->local_surface_normals.size( ), floatVector( _dim * _dim, 0 ) ) );

        //TODO: Consider doing a least-squares projection rather than a nodal averaging

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

            //Don't need to do this because it is symmetric. We do it for testing purposes
            //since to ring out the code we sometimes pass in non-symmetric tensors for the
            //micro-stresses just to check.
            floatVector symmetricMicroStress_T( _dim * _dim );

            for ( unsigned int _i = 0; _i < _dim; _i++ ){

                for ( unsigned int _j = 0; _j < _dim; _j++ ){

                    symmetricMicroStress_T[ _dim * _j + _i ] = symmetricMicroStress[ _dim * _i + _j ];

                }

            }

            for ( unsigned int j = 0; j < nMacroCellNodes; j++ ){

                //Get the shapefunction value for the node
                floatType N = shapefunctionsAtCentersOfMass[ *microDomainName ][ j ];

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

            // Project the traction and couples
            uIntVector domainMacroSurfaces = cellDomainMacroSurfaces[ macroCellID ][ *microDomainName ];
            for ( auto face = domainMacroSurfaces.begin( ); face != domainMacroSurfaces.end( ); face++ ){

                floatVector lcom( homogenizedSurfaceRegionProjectedLocalCentersOfMass[ macroCellID ][ *microDomainName ].begin( ) +
                                                                                                                _dim * ( *face ),
                                  homogenizedSurfaceRegionProjectedLocalCentersOfMass[ macroCellID ][ *microDomainName ].begin( ) +
                                                                                                                _dim * ( ( *face ) + 1 ) );
                floatVector sfs;

                error = element->get_shape_functions( lcom, sfs );

                if ( error ){

                    errorOut result = new errorNode( __func__,
                                                     "Error in the computation of the shapefunctions at the micro-domain surface center of mass for macro cell " + std::to_string( macroCellID ) + " on face " + std::to_string( *face ) );
                    result->addNext( error );
                    return result;

                }

                floatType area = homogenizedSurfaceRegionAreas[ macroCellID ][ *microDomainName ][ ( *face ) ];

                floatVector traction( homogenizedSurfaceRegionTractions[ macroCellID ][ *microDomainName ].begin( ) + _dim * ( *face ),
                                      homogenizedSurfaceRegionTractions[ macroCellID ][ *microDomainName ].begin( ) + _dim * ( ( *face ) + 1 ) );

                floatVector couple( homogenizedSurfaceRegionCouples[ macroCellID ][ *microDomainName ].begin( ) + _dim * _dim * ( *face ),
                                    homogenizedSurfaceRegionCouples[ macroCellID ][ *microDomainName ].begin( ) + _dim * _dim * ( ( *face ) + 1 ) );

                for ( unsigned int j = 0; j < nMacroCellNodes; j++ ){

                    areaAtNodes[ j ][ *face ] += sfs[ j ] * area;
                    tractionAtNodes[ j ][ *face ] += sfs[ j ] * traction * area;
                    coupleAtNodes[ j ][ *face ] += sfs[ j ] * couple * area;

                }

            }

        }

        //Save the contributions of the body forces and couples to the external force at the nodes
        externalForcesAtNodes.emplace( macroCellID, vectorTools::appendVectors( bodyForceAtNodes ) );
        externalCouplesAtNodes.emplace( macroCellID, vectorTools::appendVectors( bodyCoupleAtNodes ) );

        //De-weight the projected values at the nodes
        for ( unsigned int n = 0; n < nMacroCellNodes; n++ ){

            floatType volume;
            if ( _inputProcessor.assumeVoidlessBody( ) ){
                volume = std::fmax( _absoluteTolerance, volumeAtNodes[ n ] );
            }
            else{
                volume = std::fmax( elementNodalVolumes[ n ], volumeAtNodes[ n ] );
            }

            densityAtNodes[ n ]              /= volume;
            bodyForceAtNodes[ n ]            /= ( densityAtNodes[ n ] * volume );
            accelerationAtNodes[ n ]         /= ( densityAtNodes[ n ] * volume );
            microInertiaAtNodes[ n ]         /= ( densityAtNodes[ n ] * volume );
            bodyCoupleAtNodes[ n ]           /= ( densityAtNodes[ n ] * volume );
            microSpinInertiaAtNodes[ n ]     /= ( densityAtNodes[ n ] * volume );
            symmetricMicroStressAtNodes[ n ] /= volumeAtNodes[ n ];

            for ( unsigned int face = 0; face < element->local_surface_normals.size( ); face++ ){

                floatType area;
                if ( _inputProcessor.assumeVoidlessBody( ) ){
                    area = std::fmax( _absoluteTolerance, areaAtNodes[ n ][ face ] );
                }
                else{
                    area = std::fmax( elementNodalSurfaceAreas[ n ][ face ], areaAtNodes[ n ][ face ] );
                }

                tractionAtNodes[ n ][ face ] /= area;
                coupleAtNodes[ n ][ face ]   /= area;

            }

        }

#ifdef TESTACCESS

        _test_volumeAtNodes.emplace( macroCellID, volumeAtNodes );
        _test_densityAtNodes.emplace( macroCellID, densityAtNodes );
        _test_bodyForceAtNodes.emplace( macroCellID, bodyForceAtNodes );
        _test_accelerationAtNodes.emplace( macroCellID, accelerationAtNodes );
        _test_microInertiaAtNodes.emplace( macroCellID, microInertiaAtNodes );
        _test_bodyCoupleAtNodes.emplace( macroCellID, bodyCoupleAtNodes );
        _test_microSpinInertiaAtNodes.emplace( macroCellID, microSpinInertiaAtNodes );
        _test_symmetricMicroStressAtNodes.emplace( macroCellID, symmetricMicroStressAtNodes );

#endif

        //Add the volume integral components of the right hand side vectors
        for ( auto qpt = element->qrule.begin( ); qpt != element->qrule.end( ); qpt++ ){

            // Compute the shapefunction values
            floatVector shapefunctions;
            element->get_shape_functions( qpt->first, shapefunctions );

            // Compute the weighting factor
            floatMatrix jacobian;
            error = element->get_local_gradient( element->nodes, qpt->first, jacobian );

            if ( error ){

                std::string message = "Error in the computation of the local gradient of the shape functions for volume quadrature point " + std::to_string( qpt - element->qrule.begin( ) );
                return new errorNode( __func__, message );

            }

            floatVector inv_jacobian = vectorTools::inverse( vectorTools::appendVectors( jacobian ), _dim, _dim );
            floatType Jxw = vectorTools::determinant( vectorTools::appendVectors( jacobian ), _dim, _dim ) * qpt->second;

            // Interpolate the nodal values
            floatType   density;
            floatVector bodyForce;
            floatVector acceleration;

            floatVector bodyCouple;
            floatVector microSpinInertia;
            floatVector symmetricMicroStress;

            element->interpolate( densityAtNodes,      qpt->first, density );
            element->interpolate( bodyForceAtNodes,    qpt->first, bodyForce );
            element->interpolate( accelerationAtNodes, qpt->first, acceleration );

            element->interpolate( bodyCoupleAtNodes,           qpt->first, bodyCouple );
            element->interpolate( microSpinInertiaAtNodes,     qpt->first, microSpinInertia );
            element->interpolate( symmetricMicroStressAtNodes, qpt->first, symmetricMicroStress );

            for ( unsigned int j = 0; j < nMacroCellNodes; j++ ){

                //Get the shapefunction value for the node
                floatType N = shapefunctions[ j ];

                //Compute the contribution to the node
                floatVector nLinearMomentumRHS = N * density * ( bodyForce - acceleration ) * Jxw;

                floatVector nFirstMomentRHS = N * ( density * ( bodyCouple - microSpinInertia ) - symmetricMicroStress ) * Jxw;

                //Add the contribution to the overall RHS vectors
                for ( auto it = nLinearMomentumRHS.begin( ); it != nLinearMomentumRHS.end( ); it++ ){

                    uIntType index = _dim * j + ( it - nLinearMomentumRHS.begin( ) );

                    linearMomentumRHS[ index ] += *it;

                }

                for ( auto it = nFirstMomentRHS.begin( ); it != nFirstMomentRHS.end( ); it++ ){

                    uIntType index = _dim * _dim * j + ( it - nFirstMomentRHS.begin( ) );

                    firstMomentRHS[ index ] += *it;

                }

            }

        }

        //Add the surface integral components of the right hand side vectors
        uIntType surface_index = 0;
        for ( auto face_nodes = element->local_surface_node_ids.begin( ); face_nodes != element->local_surface_node_ids.end( ); face_nodes++, surface_index++ ){
            
            floatMatrix tractions( element->nodes.size( ), floatVector( _dim, 0 ) );
            floatMatrix couples( element->nodes.size( ), floatVector( _dim * _dim, 0 ) );

            // Initialize the nodes
            for ( auto fn = face_nodes->begin( ); fn != face_nodes->end( ); fn++ ){

                tractions[ *fn ] = tractionAtNodes[ *fn ][ surface_index ];
                couples[ *fn ] = coupleAtNodes[ *fn ][ surface_index ];

            }

            std::cerr << "surface: " << surface_index << "\n";
            // Begin the surface integration
            for ( auto s_qpt = element->surface_quadrature_rules[ surface_index ].begin( ); s_qpt != element->surface_quadrature_rules[ surface_index ].end( ); s_qpt++ ){

                // Interpolate the traction and couple
                floatVector traction;
                floatVector couple;

                element->interpolate( tractions, s_qpt->first, traction );
                element->interpolate( couples, s_qpt->first, couple );

                std::cerr << "    traction: "; vectorTools::print( traction );
                std::cerr << "    couple:   "; vectorTools::print( couple );

                // Compute the shapefunction values
                floatVector shapefunctions;
                element->get_shape_functions( s_qpt->first, shapefunctions );

                // Compute the weighting factor
                floatMatrix surface_jacobian;
                error = element->get_local_gradient( element->nodes, s_qpt->first, surface_jacobian );
                std::cerr << "    surface_jacobian: "; vectorTools::print( vectorTools::appendVectors( surface_jacobian ) );

                if ( error ){

                    std::string message = "Error in computing the local gradient of the shape functions for quadrature point "
                                        + std::to_string( ( unsigned int )( s_qpt - element->surface_quadrature_rules[ surface_index ].begin( ) ) )
                                        + " on surface " + std::to_string( surfaceNum ) + " for element " + std::to_string( macroCellID );
                    return new errorNode( __func__, message );

                }

                floatVector inv_jacobian = vectorTools::inverse( vectorTools::appendVectors( surface_jacobian ), _dim, _dim );
                std::cerr << "    inv_jacobian: "; vectorTools::print( inv_jacobian ); 
                floatType Jxw = vectorTools::determinant( vectorTools::appendVectors( surface_jacobian ), _dim, _dim ) * s_qpt->second;
                std::cerr << "    Jxw: " << Jxw << "\n";
                floatType da = vectorTools::l2norm( vectorTools::Tdot( vectorTools::inflate( inv_jacobian, _dim, _dim ), element->local_surface_normals[ surface_index ] ) * Jxw );
                std::cerr << "    da:  " << da << "\n";

                for ( unsigned int j = 0; j < element->nodes.size( ); j++ ){

                    // Compute the contribution to the node
                    floatVector nLinearMomentumRHS = shapefunctions[ j ] * traction * da;

                    floatVector nFirstMomentRHS = shapefunctions[ j ] * couple * da;

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

#ifdef TESTACCESS

        _test_cellLinearMomentumRHS.emplace( macroCellID, linearMomentumRHS );
        _test_cellFirstMomentRHS.emplace( macroCellID, firstMomentRHS );

#endif

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
            uIntType col0lm = _dim * _dim * qptIndex;
            uIntType col0fm = _dim * _dim * element->nodes.size( ) + _dim * _dim * _dim * qptIndex;

            //Get the values of the shape function and the gradients
            error = element->get_shape_functions( qpt->first, shapeFunctions );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in the computation of the shape functions\n" );
                result->addNext( error );
                return result;

            }

            //Get the values of the shape function gradients
            error = element->get_global_shapefunction_gradients( qpt->first, dNdx );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in the computation of the shape function gradients\n" );
                result->addNext( error );
                return result;

            }

            //Get the Jacobian of transformation
            error = element->get_local_gradient( element->nodes, qpt->first, jacobian );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in the computation of the local gradient\n" );
                result->addNext( error );
                return result;

            }

            Jxw = vectorTools::determinant( vectorTools::appendVectors( jacobian ), _dim, _dim ) * qpt->second;

            for ( unsigned int n = 0; n < element->nodes.size( ); n++ ){

                //Set the row
                uIntType row0 = _dim * n;

                //Add the balance of linear momentum contributions
                for ( unsigned int i = 0; i < _dim; i++ ){

                    for ( unsigned int j = 0; j < _dim; j++ ){

                        coefficients.push_back( DOFProjection::T( row0 + i, col0lm + i + _dim * j, dNdx[ n ][ j ] * Jxw ) );

                    }

                }

                //Add the balance of the first moment of momentum contributions
                row0 = _dim * element->nodes.size( ) + _dim * _dim * n;

                //Cauchy stress contribution
                for ( unsigned int i = 0; i < _dim; i++ ){

                    for ( unsigned int j = 0; j < _dim; j++ ){

                        coefficients.push_back( DOFProjection::T( row0 + _dim * i + j, col0lm + _dim * j + i, -shapeFunctions[ n ] * Jxw ) );

                    }

                }

                //Higher order stress contribution
                for ( unsigned int i = 0; i < _dim * _dim; i++ ){

                    for ( unsigned int j = 0; j < _dim; j++ ){

                        coefficients.push_back( DOFProjection::T( row0 + i, col0fm + _dim * _dim * j + i, dNdx[ n ][ j ] * Jxw ) );

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

#ifdef TESTACCESS

        _test_stressProjectionLHS.emplace( macroCellID, LHS.toDense( ) );

#endif

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

                cauchyStresses[ nCauchy * n + i ] = x( nCauchy * n + i );

            }

            for ( unsigned int i = 0; i < nHigherOrder; i++ ){

                higherOrderStresses[ nHigherOrder * n + i ] = x( nEvaluationPoints * nCauchy + nHigherOrder * n + i );

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
                                                tripletVector &coefficients,
                                                const floatVector *arlequinNodalWeights,
                                                const bool quantitiesInReference ){
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
         * :param const floatVector *arlequinNodalWeights: The weights of the nodes for the Arlequin method.
         *     Defaults to NULL.
         * :param const bool floatVector: Flag for if the density and moments of inerta are actually defined
         *     in the reference configuration.
         */

        //Get the dimension of the element
        uIntType dim = element->nodes[ 0 ].size( );

        const uIntType uSize   = dim;
        const uIntType phiSize = dim * dim;

        //Check that the degree of freedom value vector's length is consistent with the element
        if ( degreeOfFreedomValues.size( ) != ( uSize + phiSize ) * element->nodes.size( ) ){

            return new errorNode( __func__,
                                  "The degree of freedom vector size is not consistent with the element dimension" );

        }

        if ( momentOfInertia.size( ) != element->qrule.size( ) * phiSize ){

            return new errorNode( __func__,
                                  "The moment of inertia vector size is not consistent with the quadrature rule and element dimension" );

        }

        if ( density.size( ) != element->qrule.size( ) ){

            return new errorNode( __func__,
                                  "The density vector size is not consistent with the quadrature rule" );

        }

        if ( element->global_node_ids.size( ) != element->nodes.size( )  ){

            return new errorNode( __func__,
                                  "The size of the global node id in the element are not the same size as the number of nodes" );

        }

        //Reshape the degree of freedom values to a matrix of values where the rows are the values at the nodes
        floatMatrix reshapedDOFValues = vectorTools::inflate( degreeOfFreedomValues, element->nodes.size( ), uSize + phiSize );

        //Variable initialize
        floatVector shapeFunctions;
        floatVector interpolatedValues, deformationGradient;
        floatType qptDensity, referenceDensity;
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

                errorOut result = new errorNode( __func__,
                                                 "Error in the computation of the required values for the element" );
                result->addNext( error );
                return result;

            }

            invXiQpt = vectorTools::inverse( XiQpt, dim, dim );

            //Get the current density
            qptDensity = density[ qptIndex ];

            //Compute the moment of inertia in the reference configuration
            qptMomentOfInertia = floatVector( momentOfInertia.begin( ) + dim * dim * qptIndex,
                                              momentOfInertia.begin( ) + dim * dim * ( qptIndex + 1 ) );

            if ( quantitiesInReference ){

                referenceDensity = qptDensity;
                referenceMomentOfInertia = qptMomentOfInertia;

            }
            else{

                referenceDensity = J * qptDensity;
                referenceMomentOfInertia
                    = vectorTools::matrixMultiply( vectorTools::matrixMultiply( invXiQpt, qptMomentOfInertia, dim, dim, dim, dim ),
                                                   invXiQpt, dim, dim, dim, dim, false, true );

            }

            //Evaluate the integrand term
            inertiaTerm = referenceDensity * referenceMomentOfInertia * Jxw;

            floatType arlequinWeight = 1;
            if ( arlequinNodalWeights ){

                element->interpolate( *arlequinNodalWeights, qpt->first, arlequinWeight );

            }

            //Add the contrubutions to the mass matrix
            for ( uIntType o = 0; o < shapeFunctions.size( ); o++ ){

                sFo = shapeFunctions[ o ];

                auto gni1 = nodeIDToIndex->find( element->global_node_ids[ o ] );

                if ( gni1 == nodeIDToIndex->end( ) ){

                    return new errorNode( __func__,
                                          "Node " + std::to_string( element->global_node_ids[ o ] ) + " not found in the ID map" );
                                          

                }

                row0 = ( uSize + phiSize ) * gni1->second;

                for ( uIntType p = 0; p < shapeFunctions.size( ); p++ ){

                    sFp = shapeFunctions[ p ];

                    auto gni2 = nodeIDToIndex->find( element->global_node_ids[ p ] );

                    if ( gni2 == nodeIDToIndex->end( ) ){
    
                        return new errorNode( __func__,
                                              "Node " + std::to_string( element->global_node_ids[ p ] ) + " not found in the ID map" );
                                              
    
                    }
    
                    col0 = ( uSize + phiSize ) * gni2->second;

                    for ( unsigned int j = 0; j < dim; j++ ){

                        for ( unsigned int k = 0; k < dim; k++ ){

                            coefficients.push_back( DOFProjection::T( row0 + j,
                                                                      col0 + k,
                                                                      arlequinWeight * eye[ dim * j + k ] * referenceDensity * sFo * sFp * Jxw ) );
    
                            for ( unsigned int K = 0; K < dim; K++ ){
    
                                for ( unsigned int L = 0; L < dim; L++ ){
    
                                    coefficients.push_back( DOFProjection::T( row0 + dim + dim * j + K,
                                                                              col0 + dim + dim * k + L,
                                                                              arlequinWeight * eye[ dim * j + k ] * sFo * sFp * inertiaTerm[ dim * K + L ] ) );

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
         * TODO: This function is considerably more complicated than it needs to be so that it
         *       uses the functions in the micromorphic balance equations library. We don't have to
         *       do this and it would make the function shorter and probably faster.
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
            output             += "because of dNdX, fint, and cint which are currently consistent with a 3D problem as required\n";
            output             += "by balance_equations.h";
            return new errorNode( __func__, output );

        }

        //Check that the degree of freedom value vector's length is consistent with the element
        if ( degreeOfFreedomValues.size( ) != ( uSize + phiSize ) * element->nodes.size( ) ){

            return new errorNode( __func__,
                                  "The degree of freedom vector size is not consistent with the element dimension" );

        }

        if ( cauchyStress.size( ) != element->qrule.size( ) * dim * dim ){

            return new errorNode( __func__,
                                  "The Cauchy stress vector size is not consistent with the quadrature rule and element dimension" );

        }

        if ( symmetricMicroStress.size( ) != element->qrule.size( ) * dim * dim ){

            return new errorNode( __func__,
                                  "The symmetric micro-stress vector size is not consistent with the quadrature rule" );

        }

        if ( higherOrderStress.size( ) != element->qrule.size( ) * dim * dim * dim ){

            return new errorNode( __func__,
                                  "The higher-order stress vector size is not consistent with the quadrature rule" );

        }

        if ( element->global_node_ids.size( ) != element->nodes.size( )  ){

            return new errorNode( __func__,
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

                errorOut result = new errorNode( __func__,
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

                errorOut result = new errorNode( __func__,
                                                 "Error in the pull-back operation on the Cauchy stress" );
                result->addNext( error );
                return result;

            }

            //Pull back the symmetric micro-stress
            error = micromorphicTools::pullBackMicroStress( sQpt, deformationGradient, referenceMicroStressQpt );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error in the pull-back operation on the symmetric micro-stress" );
                result->addNext( error );
                return result;

            }

            //Pull back the higher order stress
            error = micromorphicTools::pullBackHigherOrderStress( mQpt, deformationGradient, XiQpt, referenceHigherOrderStressQpt );

            if ( error ){

                errorOut result = new errorNode( __func__,
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

                    return new errorNode( __func__,
                                          "The internal force term returned an error code: " + std::to_string( errorCode ) );

                }

                //Compute the terms for the balance of first moment of momentum
                errorCode = balance_equations::compute_internal_couple( N, dNdX, deformationGradient, XiQpt,
                                                                        pk2Qpt, referenceMicroStressQpt,
                                                                        referenceHigherOrderStressQpt, cint );

                if ( errorCode != 0 ){

                    return new errorNode( __func__,
                                          "The internal couple term returned an error code: " + std::to_string( errorCode ) );

                }

                //Get the initial index
                auto it = nodeIDToIndex->find( element->global_node_ids[ n ] );
                
                if ( it == nodeIDToIndex->end( ) ){

                    return new errorNode( __func__,
                                          "The global node id " + std::to_string( element->global_node_ids[ n ] ) +
                                          " is not found in the id to index map" );

                }

                //Set the row index
                row0 = ( uSize + phiSize ) * it->second;

                if ( ( row0 + uSize + phiSize ) > internalForceVector.rows( ) ){

                    return new errorNode( __func__,
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

            errorOut result = new errorNode( __func__,
                                             "Error in the computation of the shape functions" );
            result->addNext( error );
            return result;

        }

        //Evaluate the gradients of the shape functions
        error = element->get_global_shapefunction_gradients( qpt->first, gradShapeFunctions, useReference );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error in the computation of the shape function gradients" );
            result->addNext( error );
            return result;

        }

        //Get the deformation gradient between the reference and current configurations
        error = element->get_jacobian( qpt->first, element->reference_nodes, jacobian );

        if ( error ){

            errorOut result = new errorNode( __func__,
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

            errorOut result = new errorNode( __func__,
                                             "Error in the computation of the local gradient\n" );
            result->addNext( error );
            return result;

        }

        Jxw = vectorTools::determinant( vectorTools::appendVectors( jacobian ), dim, dim ) * qpt->second;

        //Interpolate the DOF nodes to the node
        error = element->interpolate( reshapedDOFValues, qpt->first, interpolatedValues );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error in the interpolation of the degree of freedom values" );
            result->addNext( error );
            return result;

        }

        if ( interpolatedValues.size( ) < ( dim + dim * dim ) ){

            std::string output = "The interpolated values shape is not consistent with the required dimension for the displacement ";
            output            += "and micro-displacement interpolation";
            return new errorNode( __func__, output );

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

                return new errorNode( __func__,
                                      "Macro cell ID " + std::to_string( *macroCellID ) +
                                      " not found in external forces at nodes." );

            }

            //Make sure that the macroCellID is stored in the external couple vector
            if ( externalCouplesAtNodes.find( *macroCellID ) == externalCouplesAtNodes.end( ) ){

                return new errorNode( __func__,
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

                errorOut result = new errorNode( __func__,
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

                    return new errorNode( __func__,
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

                return new errorNode( __func__,
                                      "Global degree of freedom '" + std::to_string( *it ) + "' not found in degree of freedom map" );

            }

            auto macroDisplacement = macroDispDOFVector->find( *it );

            if ( macroDisplacement == macroDispDOFVector->end( ) ){

                return new errorNode( __func__,
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

                errorOut result = new errorNode( __func__,
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

                    return new errorNode( __func__,
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

                        return new errorNode( __func__,
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

                errorOut result = new errorNode( __func__,
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

                return new errorNode( __func__,
                                      "Global degree of freedom '" + std::to_string( *it ) + "' not found in degree of freedom map" );

            }

            auto macroDisplacement = macroDispDOFVector->find( *it );

            if ( macroDisplacement == macroDispDOFVector->end( ) ){

                return new errorNode( __func__,
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

                errorOut result = new errorNode( __func__,
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

                    return new errorNode( __func__,
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

                        return new errorNode( __func__,
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

            errorOut result = new errorNode( __func__,
                                             "Error in the construction of the homogenized external force vector" );
            result->addNext( error );
            return result;

        }

        error = assembleHomogenizedInternalForceVector( );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error in the construction of the homogenized internal force vector" );
            result->addNext( error );
            return result;

        }

        error = assembleHomogenizedMassMatrix( );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error in the construction of the homogenized mass matrix" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut overlapCoupling::assembleFreeMicromorphicMassMatrix( ){
        /*!
         * Assemble the micromorphic mass matrix for the free micromorphic domains.
         *
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

                errorOut result = new errorNode( __func__,
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

                    return new errorNode( __func__,
                                          "Macro node " + std::to_string( *nodeID ) + " was not found in the macro displacement DOF vector map" );

                }

                elementDOFVector
                    = vectorTools::appendVectors( { elementDOFVector, macroDisplacement->second } );

            }

            //Extract the density and moment of inertia in the reference configuration
            auto densityType = macroReferenceDensityTypes->find( *macroCellID );

            if ( densityType == macroReferenceDensityTypes->end( ) ){

                return new errorNode( __func__,
                                      "The macro cell with ID " + std::to_string( *macroCellID ) +
                                      " was not found in the density type map" );

            }

            auto momentOfInertiaType = macroReferenceMomentOfInertiaTypes->find( *macroCellID );

            if ( momentOfInertiaType == macroReferenceMomentOfInertiaTypes->end( ) ){

                return new errorNode( __func__,
                                      "The macro cell with ID " + std::to_string( *macroCellID ) +
                                      " was not found in the moment of inertia type map" );

            }

            if ( densityType->second.compare( "constant" ) != 0 ){

                return new errorNode( __func__,
                                      "Only constant densities for the macro-scale are allowed currently. This is not true for macro cell ID " + std::to_string( *macroCellID ) );

            }

            if ( momentOfInertiaType->second.compare( "constant" ) != 0 ){

                return new errorNode( __func__,
                                      "Only constant moments of inertia for the macro-scale are allowed currently. This is not true for macro cell ID " + std::to_string( *macroCellID ) );

            }

            bool quantitiesInReference = true; //TODO: We currently assume all of the quantities are defined in the reference configuration

            auto macroDensities = macroReferenceDensities->find( *macroCellID );

            if ( macroDensities == macroReferenceDensities->end( ) ){

                return new errorNode( __func__,
                                      "Macro cell ID " + std::to_string( *macroCellID ) +
                                      " is not in the macro reference density map" );

            }

            if ( macroDensities->second.size( ) != 1 ){

                return new errorNode( __func__,
                                      "The macro densities for macro cell " + std::to_string( *macroCellID ) +
                                      "Define " + std::to_string( macroDensities->second.size( ) ) +
                                      " values when only 1 can be defined" );

            }

            auto macroMomentsOfInertia = macroReferenceMomentsOfInertia->find( *macroCellID );

            if ( macroMomentsOfInertia == macroReferenceMomentsOfInertia->end( ) ){

                return new errorNode( __func__,
                                      "Macro cell ID " + std::to_string( *macroCellID ) +
                                      " is not in the macro reference moments of inertia map" );

            }

            if ( macroMomentsOfInertia->second.size( ) != _dim * _dim ){

                return new errorNode( __func__,
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

                errorOut result = new errorNode( __func__,
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
                                                       _inputProcessor.getMacroGlobalToLocalDOFMap( ), coefficients,
                                                       NULL, quantitiesInReference );

            if ( error ){

                std::string outstr  = "Error in the construction of the contributions of the macro element to ";
                            outstr += "the free micromorphic mass matrix";

                errorOut result = new errorNode( __func__, outstr );
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
        const std::unordered_map< uIntType, floatType > *microDensities = _inputProcessor.getMicroDensities( );

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
        floatVector freeMicroMasses( freeMicroNodeIDs->size( ), 0 );

        //Assemble the free micro mass vector
        for ( auto microID = freeMicroNodeIDs->begin( ); microID != freeMicroNodeIDs->end( ); microID++ ){

            auto localMicroNodeIDMap = microGlobalToLocalDOFMap->find( *microID );

            if ( localMicroNodeIDMap == microGlobalToLocalDOFMap->end( ) ){

                return new errorNode( __func__,
                                      "Free micro node: " + std::to_string( *microID ) + " not found in global to local map\n" );

            }

            auto microDensity = microDensities->find( *microID );

            if ( microDensity == microDensities->end( ) ){

                return new errorNode( __func__,
                                      "Free micro node " + std::to_string( *microID ) + " was not found in the micro density map" );

            }

            auto microVolume = microVolumes->find( *microID );

            if ( microVolume == microVolumes->end( ) ){

                return new errorNode( __func__,
                                      "Free micro node " + std::to_string( *microID ) + " was not found in the micro volume map" );

            }

            freeMicroMasses[ localMicroNodeIDMap->second ] = microVolume->second * microDensity->second;

        }

#ifdef TESTACCESS

        _test_freeMicroMasses = freeMicroMasses;

#endif

        //Assemble the ghost micro mass vector
        for ( auto microID = ghostMicroNodeIDs->begin( ); microID != ghostMicroNodeIDs->end( ); microID++ ){

            auto localMicroNodeIDMap = microGlobalToLocalDOFMap->find( *microID );

            if ( localMicroNodeIDMap == microGlobalToLocalDOFMap->end( ) ){

                return new errorNode( __func__,
                                      "Ghost micro node: " + std::to_string( *microID ) + " not found in global to local map\n" );

            }

            auto microDensity = microDensities->find( *microID );

            if ( microDensity == microDensities->end( ) ){

                return new errorNode( __func__,
                                      "Free micro node " + std::to_string( *microID ) + " was not found in the micro density map" );

            }

            auto microVolume = microVolumes->find( *microID );

            if ( microVolume == microVolumes->end( ) ){

                return new errorNode( __func__,
                                      "Free micro node " + std::to_string( *microID ) + " was not found in the micro volume map" );

            }

            ghostMicroMasses[ localMicroNodeIDMap->second - nFreeMicroNodes ] = microVolume->second * microDensity->second;

        }

#ifdef TESTACCESS

        _test_ghostMicroMasses = ghostMicroMasses;

#endif

        //Assemble the mass sub-matrices
//        Eigen::VectorXd mq( 3 * freeMicroMasses.size( ) );
//        Eigen::VectorXd mqhat( 3 * ghostMicroMasses.size( ) );
        tripletVector c1;
        tripletVector c2;

        c1.reserve( _dim * ghostMicroMasses.size( ) ); 
        c2.reserve( _dim * freeMicroMasses.size( ) ); 

        uIntType mIndex = 0;

        for ( auto m = ghostMicroMasses.begin( ); m != ghostMicroMasses.end( ); m++, mIndex++ ){

            for ( unsigned int i = 0; i < _dim; i++ ){

                c1.push_back( DOFProjection::T( _dim * mIndex + i, _dim * mIndex + i, ( 1 - rhat ) * ( *m ) ) );
//                mqhat( _dim * mIndex + i ) = ( 1 - rhat ) * ( *m );

            }

        }

        mIndex = 0;

        for ( auto m = freeMicroMasses.begin( ); m != freeMicroMasses.end( ); m++, mIndex++ ){

            for ( unsigned int i = 0; i < _dim; i++ ){

                c2.push_back( DOFProjection::T( _dim * mIndex + i, _dim * mIndex + i, ( 1 - rhat ) * ( *m ) ) );
//                mq( _dim * mIndex + i ) = ( 1 - rhat ) * ( *m );

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

#ifdef TESTACCESS
       
        _test_MQ = MQ;
        _test_MQhat = MQhat;

#endif
        //Due to the restrictions in listed in the comment at the beginning of the function, MBar from Regueiro 2012
        //is an empty matrix as we only handle the coupling domain here.
       
    
        //Assemble Mass matrices for the micro projection equation

        if ( ( config[ "projection_type" ].as< std::string >( ).compare( "l2_projection" ) == 0 ) ){

            //Determine the coefficients where we cut off the projectors TODO: Expose factor to user
            floatType BQhatQ_coeff = 1e-4 * ( std::fabs( _dense_BQhatQ.minCoeff( ) ) + std::fabs( _dense_BQhatQ.maxCoeff( ) ) );
            floatType BQhatD_coeff = 1e-4 * ( std::fabs( _dense_BQhatD.minCoeff( ) ) + std::fabs( _dense_BQhatD.maxCoeff( ) ) );
            floatType BDhatQ_coeff = 1e-4 * ( std::fabs( _dense_BDhatQ.minCoeff( ) ) + std::fabs( _dense_BDhatQ.maxCoeff( ) ) );
            floatType BDhatD_coeff = 1e-4 * ( std::fabs( _dense_BDhatD.minCoeff( ) ) + std::fabs( _dense_BDhatD.maxCoeff( ) ) );

            SparseMatrix _sparse_BQhatQ = _dense_BQhatQ.sparseView( 1, BQhatQ_coeff );
            SparseMatrix _sparse_BQhatD = _dense_BQhatD.sparseView( 1, BQhatD_coeff );
            SparseMatrix _sparse_BDhatQ = _dense_BDhatQ.sparseView( 1, BDhatQ_coeff );
            SparseMatrix _sparse_BDhatD = _dense_BDhatD.sparseView( 1, BDhatD_coeff );

            //TODO: Improve efficiency

            std::cerr << "  ASSEMBLING MQQ\n";
//            std::cout << "coeffs:\n" << BQhatQ_coeff << " " << BQhatD_coeff << " " << BDhatQ_coeff << " " << BDhatD_coeff << "\n";
//
//            std::cout << "Zero ratios\n";
//            std::cerr << ( ( floatType )_sparse_BQhatQ.nonZeros( ) ) / ( ( floatType )_sparse_BQhatQ.size( ) ) << "\n";
//            std::cerr << ( ( floatType )_sparse_BQhatD.nonZeros( ) ) / ( ( floatType )_sparse_BQhatD.size( ) ) << "\n";
//            std::cerr << ( ( floatType )_sparse_BDhatQ.nonZeros( ) ) / ( ( floatType )_sparse_BDhatQ.size( ) ) << "\n";
//            std::cerr << ( ( floatType )_sparse_BDhatD.nonZeros( ) ) / ( ( floatType )_sparse_BDhatD.size( ) ) << "\n";

            SparseMatrix MQQ = MQ;
//            std::cerr << "mul1\n";
            MQQ += _sparse_BQhatQ.transpose( ) * MQhat * _sparse_BQhatQ;
//            std::cerr << "mul2\n";
            MQQ += _sparse_BDhatQ.transpose( ) * MDhat * _sparse_BDhatQ;

            std::cerr << "  ASSEMBLING MQD\n";
            SparseMatrix MQD = _sparse_BQhatQ.transpose( ) * MQhat * _sparse_BQhatD;
            MQD += _sparse_BDhatQ.transpose( ) * MDhat * _sparse_BDhatD;
    
            //Assemble Mass matrices for the macro projection equation
            
            std::cerr << "  ASSEMBLING MDQ\n";
            SparseMatrix MDQ = _sparse_BQhatD.transpose( ) * MQhat * _sparse_BQhatQ;//MQhat * _dense_BQhatQ;
            MDQ += _sparse_BDhatD.transpose( ) * MDhat * _sparse_BDhatQ;

            std::cerr << "  ASSEMBLING MDD\n";
            SparseMatrix MDD = MD;
            MDD += _sparse_BQhatD.transpose( ) * MQhat * _sparse_BQhatD;//MQhat * _dense_BQhatD;
            MDD += _sparse_BDhatD.transpose( ) * MDhat * _sparse_BDhatD;
    
            //Assemble the damping matrices for the micro projection equation
            std::cerr << "  ASSEMBLING CQQ\n";
            SparseMatrix CQQ = aQ * MQ;
            CQQ += aQ * _sparse_BQhatQ.transpose( ) * MQhat * _sparse_BQhatQ;//MQhat * _dense_BQhatQ;
            CQQ += aD * _sparse_BDhatQ.transpose( ) * MDhat * _sparse_BDhatQ;

            std::cerr << "  ASSEMBLING CQD\n";
            SparseMatrix CQD = aQ * _sparse_BQhatQ.transpose( ) * MQhat * _sparse_BQhatD;//MQhat * _dense_BQhatD;
    
            //Assemble the damping matrices for the macro projection equation
            std::cerr << "  ASSEMBLING CDQ\n";
            SparseMatrix CDQ = aQ * _sparse_BQhatD.transpose( ) * MQhat * _sparse_BQhatQ;//MQhat * _dense_BQhatQ;

            std::cerr << "  ASSEMBLING CDD\n";
            SparseMatrix CDD = aD * MD;
            CDD += aQ * _sparse_BQhatD.transpose( ) * MQhat * _sparse_BQhatD;//MQhat * _dense_BQhatD;

//            //Assemble the full mass matrix
//            std::cerr << "  ASSEMBLING THE FULL MASS MATRIX\n";
//            _dense_MASS = Eigen::MatrixXd( MQQ.rows( ) + MDQ.rows( ), MQQ.cols( ) + MQD.cols( ) );
//            _dense_MASS.topLeftCorner(     MQQ.rows( ), MQQ.cols( ) ) = MQQ;
//            _dense_MASS.topRightCorner(    MQD.rows( ), MQD.cols( ) ) = MQD;
//            _dense_MASS.bottomLeftCorner(  MDQ.rows( ), MDQ.cols( ) ) = MDQ;
//            _dense_MASS.bottomRightCorner( MDD.rows( ), MDD.cols( ) ) = MDD;
//
//            //Assemble the full damping matrix
//            std::cerr << "  ASSEMBLING THE FULL DAMPING MATRIX\n";
//            _dense_DAMPING = Eigen::MatrixXd( CQQ.rows( ) + CDQ.rows( ), CQQ.cols( ) + CQD.cols( ) );
//            _dense_DAMPING.topLeftCorner(     CQQ.rows( ), CQQ.cols( ) ) = CQQ;
//            _dense_DAMPING.topRightCorner(    CQD.rows( ), CQD.cols( ) ) = CQD;
//            _dense_DAMPING.bottomLeftCorner(  CDQ.rows( ), CDQ.cols( ) ) = CDQ;
//            _dense_DAMPING.bottomRightCorner( CDD.rows( ), CDD.cols( ) ) = CDD;

#ifdef TESTACCESS

            _test_dense_MQQ = MQQ;
            _test_dense_MQD = MQD;
            _test_dense_MDQ = MDQ;
            _test_dense_MDD = MDD;

            _test_dense_CQQ = CQQ;
            _test_dense_CQD = CQD;
            _test_dense_CDQ = CDQ;
            _test_dense_CDD = CDD;

#endif
            //The sparse matrices are in column major format so we loop over the columns to assemble the matrix
            std::cerr << "  ASSEMBLING FULL MASS AND DAMPING MATRICES\n";
            _sparse_MASS    = SparseMatrix( MQQ.rows( ) + MDQ.rows( ), MQQ.cols( ) + MQD.cols( ) );
            _sparse_DAMPING = SparseMatrix( CQQ.rows( ) + CDQ.rows( ), CQQ.cols( ) + CQD.cols( ) );

            _sparse_MASS.reserve( MQQ.nonZeros( ) + MQD.nonZeros( ) + MDQ.nonZeros( ) + MDD.nonZeros( ) );
            _sparse_DAMPING.reserve( CQQ.nonZeros( ) + CQD.nonZeros( ) + CDQ.nonZeros( ) + CDD.nonZeros( ) );
            for ( uIntType c = 0; c < MQQ.cols( ); ++c ){

                //Add terms to the mass matrix
                _sparse_MASS.startVec( c );
                for ( SparseMatrix::InnerIterator itMQQ( MQQ, c ); itMQQ; ++itMQQ )
                    _sparse_MASS.insertBack( itMQQ.row( ), c ) = itMQQ.value( );
                for ( SparseMatrix::InnerIterator itMDQ( MDQ, c ); itMDQ; ++itMDQ )
                    _sparse_MASS.insertBack( itMDQ.row( ) + MQQ.rows( ), c ) = itMDQ.value( );

                //Add terms to the damping matrix
                _sparse_DAMPING.startVec( c );
                for ( SparseMatrix::InnerIterator itCQQ( CQQ, c ); itCQQ; ++itCQQ )
                    _sparse_DAMPING.insertBack( itCQQ.row( ), c ) = itCQQ.value( );
                for ( SparseMatrix::InnerIterator itCDQ( CDQ, c ); itCDQ; ++itCDQ )
                    _sparse_DAMPING.insertBack( itCDQ.row( ) + CQQ.rows( ), c ) = itCDQ.value( );

            }

            for ( uIntType c = 0; c < MDD.cols( ); ++c ){

                //Add terms to the mass matrix
                _sparse_MASS.startVec( c + MQQ.cols( ) );
                for ( SparseMatrix::InnerIterator itMQD( MQD, c ); itMQD; ++itMQD )
                    _sparse_MASS.insertBack( itMQD.row( ), c + MQQ.cols( ) ) = itMQD.value( );
                for ( SparseMatrix::InnerIterator itMDD( MDD, c ); itMDD; ++itMDD )
                    _sparse_MASS.insertBack( itMDD.row( ) + MQD.rows( ), c + MDQ.cols( ) ) = itMDD.value( );

                //Add terms to the damping matrix
                _sparse_DAMPING.startVec( c + CQQ.cols( ) );
                for ( SparseMatrix::InnerIterator itCQD( CQD, c ); itCQD; ++itCQD )
                    _sparse_DAMPING.insertBack( itCQD.row( ), c + CQQ.cols( ) ) = itCQD.value( );
                for ( SparseMatrix::InnerIterator itCDD( CDD, c ); itCDD; ++itCDD )
                    _sparse_DAMPING.insertBack( itCDD.row( ) + CQD.rows( ), c + CDQ.cols( ) ) = itCDD.value( );

            }

        }
        else if ( ( config[ "projection_type" ].as< std::string >( ).compare( "direct_projection" ) == 0 ) ||
                  ( config[ "projection_type" ].as< std::string >( ).compare( "averaged_l2_projection" ) == 0 )
                ){

            std::cerr << "ASSEMBLING MASS BLOCK MATRICES\n";
            SparseMatrix MQQ  =  MQ;
            std::cerr << "  MQQ\n";
            MQQ += _sparse_BQhatQ.transpose( ) * MQhat * _sparse_BQhatQ;
            MQQ += _sparse_BDhatQ.transpose( ) * MDhat * _sparse_BDhatQ;

            std::cerr << "  MQD\n";
            SparseMatrix MQD = _sparse_BQhatQ.transpose( ) * MQhat * _sparse_BQhatD;
            MQD += _sparse_BDhatQ.transpose( ) * MDhat * _sparse_BDhatD;
    
            //Assemble Mass matrices for the macro projection equation
            
            std::cerr << "  MDQ\n";
            SparseMatrix MDQ = _sparse_BQhatD.transpose( ) * MQhat * _sparse_BQhatQ;
            MDQ += _sparse_BDhatD.transpose( ) * MDhat * _sparse_BDhatQ;

            std::cerr << "  MDD\n";
            SparseMatrix MDD = MD;
            MDD += _sparse_BQhatD.transpose( ) * MQhat * _sparse_BQhatD;
            MDD += _sparse_BDhatD.transpose( ) * MDhat * _sparse_BDhatD;

            std::cout << "ASSEMBLING DAMPING BLOCK MATRICES\n";
    
            //Assemble the damping matrices for the micro projection equation
            std::cerr << "  CQQ\n";
            SparseMatrix CQQ = aQ * MQ;
            CQQ += aQ * _sparse_BQhatQ.transpose( ) * MQhat * _sparse_BQhatQ;
            CQQ += aD * _sparse_BDhatQ.transpose( ) * MDhat * _sparse_BDhatQ;

            std::cerr << "  CQD\n";
            SparseMatrix CQD = aQ * _sparse_BQhatQ.transpose( ) * MQhat * _sparse_BQhatD;
    
            //Assemble the damping matrices for the macro projection equation
            std::cerr << "  CDQ\n";
            SparseMatrix CDQ = aQ * _sparse_BQhatD.transpose( ) * MQhat * _sparse_BQhatQ;
            std::cerr << "  CDD\n";
            SparseMatrix CDD = aD * MD;
            CDD += aQ * _sparse_BQhatD.transpose( ) * MQhat * _sparse_BQhatD;

#ifdef TESTACCESS

            _test_dense_MQQ = MQQ;
            _test_dense_MQD = MQD;
            _test_dense_MDQ = MDQ;
            _test_dense_MDD = MDD;

            _test_dense_CQQ = CQQ;
            _test_dense_CQD = CQD;
            _test_dense_CDQ = CDQ;
            _test_dense_CDD = CDD;

#endif

            //Assemble the full mass and damping matrices
            std::cout << "ASSEMBLING FULL MASS AND DAMPING MATRICES\n";
            _sparse_MASS    = SparseMatrix( MQQ.rows( ) + MDQ.rows( ), MQQ.cols( ) + MQD.cols( ) );
            _sparse_DAMPING = SparseMatrix( CQQ.rows( ) + CDQ.rows( ), CQQ.cols( ) + CQD.cols( ) );

            //The sparse matrices are in column major format so we loop over the columns to assemble the matrix
            _sparse_MASS.reserve( MQQ.nonZeros( ) + MQD.nonZeros( ) + MDQ.nonZeros( ) + MDD.nonZeros( ) );
            _sparse_DAMPING.reserve( CQQ.nonZeros( ) + CQD.nonZeros( ) + CDQ.nonZeros( ) + CDD.nonZeros( ) );
            for ( uIntType c = 0; c < MQQ.cols( ); ++c ){

                //Add terms to the mass matrix
                _sparse_MASS.startVec( c );
                for ( SparseMatrix::InnerIterator itMQQ( MQQ, c ); itMQQ; ++itMQQ )
                    _sparse_MASS.insertBack( itMQQ.row( ), c ) = itMQQ.value( );
                for ( SparseMatrix::InnerIterator itMDQ( MDQ, c ); itMDQ; ++itMDQ )
                    _sparse_MASS.insertBack( itMDQ.row( ) + MQQ.rows( ), c ) = itMDQ.value( );

                //Add terms to the damping matrix
                _sparse_DAMPING.startVec( c );
                for ( SparseMatrix::InnerIterator itCQQ( CQQ, c ); itCQQ; ++itCQQ )
                    _sparse_DAMPING.insertBack( itCQQ.row( ), c ) = itCQQ.value( );
                for ( SparseMatrix::InnerIterator itCDQ( CDQ, c ); itCDQ; ++itCDQ )
                    _sparse_DAMPING.insertBack( itCDQ.row( ) + CQQ.rows( ), c ) = itCDQ.value( );

            }

            for ( uIntType c = 0; c < MDD.cols( ); ++c ){

                //Add terms to the mass matrix
                _sparse_MASS.startVec( c + MQQ.cols( ) );
                for ( SparseMatrix::InnerIterator itMQD( MQD, c ); itMQD; ++itMQD )
                    _sparse_MASS.insertBack( itMQD.row( ), c + MQQ.cols( ) ) = itMQD.value( );
                for ( SparseMatrix::InnerIterator itMDD( MDD, c ); itMDD; ++itMDD )
                    _sparse_MASS.insertBack( itMDD.row( ) + MQD.rows( ), c + MDQ.cols( ) ) = itMDD.value( );

                //Add terms to the damping matrix
                _sparse_DAMPING.startVec( c + CQQ.cols( ) );
                for ( SparseMatrix::InnerIterator itCQD( CQD, c ); itCQD; ++itCQD )
                    _sparse_DAMPING.insertBack( itCQD.row( ), c + CQQ.cols( ) ) = itCQD.value( );
                for ( SparseMatrix::InnerIterator itCDD( CDD, c ); itCDD; ++itCDD )
                    _sparse_DAMPING.insertBack( itCDD.row( ) + CQD.rows( ), c + CDQ.cols( ) ) = itCDD.value( );

            }

        }
        else{

            return new errorNode( __func__,
                                  "The projection type " + config[ "projection_type" ].as< std::string >( ) +
                                  " is not recognized" );

        }

        std::cerr << "MASS AND DAMPING MATRICES ASSEMBLED\n";

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
 
                    errorOut result = new errorNode( __func__,
                                                     "Error in the computation of the local gradient\n" );
                    result->addNext( error );
                    return result;

                }

                elementVolume += vectorTools::determinant( vectorTools::appendVectors( jacobian ), _dim, _dim ) * qpt->second;

            }

            auto microDomainVolumes = homogenizedVolumes.find( macroCellID );
            if ( microDomainVolumes == homogenizedVolumes.end( ) ){

                return new errorNode( __func__,
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

            return new errorNode( __func__, "Configuration strategy " + strategy + " not recognized" );

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

                return new errorNode( __func__,
                                      "ghost micro node " + std::to_string( *microID ) +
                                      " is not found in the local to global DOF map" );

            }

            if ( idMap->second < nFreeMicroNodes ){

                return new errorNode( __func__, "The local index is smaller than the number of free micro nodes" );

            }

            if ( ( _dim * ( idMap->second - nFreeMicroNodes ) + _dim ) > FintQhat.size( ) ){

                return new errorNode( __func__, "Local index is larger than the force vector" );

            }

            auto internalForce = microInternalForces->find( *microID );

            if ( internalForce == microInternalForces->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( *microID ) + " not found in internal force vector" );

            }

            auto externalForce = microExternalForces->find( *microID );

            if ( externalForce == microInternalForces->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( *microID ) + " not found in external force vector" );

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

                return new errorNode( __func__,
                                      "free micro node " + std::to_string( *microID ) +
                                      " is not found in the local to global DOF map" );

            }

            auto internalForce = microInternalForces->find( *microID );

            if ( internalForce == microInternalForces->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( *microID ) + " not found in internal force vector" );

            }

            auto externalForce = microExternalForces->find( *microID );

            if ( externalForce == microInternalForces->end( ) ){

                return new errorNode( __func__, "Micro node " + std::to_string( *microID ) + " not found in external force vector" );

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

            errorOut result = new errorNode( __func__,
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

                return new errorNode( __func__,
                                      "ghost macro node " + std::to_string( *nodeID ) +
                                      " is not found in the local to global DOF map" );

            }

            auto nodeQes = qes.find( *nodeID );

            if ( idMap->second < nFreeMacroNodes ){

                return new errorNode( __func__,
                                      "ghost macro node " + std::to_string( *nodeID ) +
                                      " has a local position not consistent with the number of free macro nodes" );

            }

            if ( ( nMacroNodeForces * ( idMap->second - nFreeMacroNodes ) + nMacroNodeForces ) > FintDhat.size( ) ){

                return new errorNode( __func__,
                                      "ghost macro node " + std::to_string( *nodeID ) +
                                      " has a local position larger than allocated in the coupling force vector" );

            }

            auto internalForce = macroInternalForces->find( *nodeID );

            if ( internalForce == macroInternalForces->end( ) ){

                return new errorNode( __func__, "Macro node " + std::to_string( *nodeID ) +
                                      " not found in internal force vector" );

            }

            auto externalForce = macroExternalForces->find( *nodeID );

            if ( externalForce == macroInternalForces->end( ) ){

                return new errorNode( __func__, "Macro node " + std::to_string( *nodeID ) +
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

                return new errorNode( __func__,
                                      "free macro node " + std::to_string( *nodeID ) +
                                      " is not found in the local to global DOF map" );

            }

            auto nodeQes = qes.find( *nodeID );

            auto internalForce = macroInternalForces->find( *nodeID );

            if ( internalForce == macroInternalForces->end( ) ){

                return new errorNode( __func__, "Macro node " + std::to_string( *nodeID ) +
                                      " not found in internal force vector" );

            }

            auto externalForce = macroExternalForces->find( *nodeID );

            if ( externalForce == macroInternalForces->end( ) ){

                return new errorNode( __func__, "Macro node " + std::to_string( *nodeID ) +
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

#ifdef TESTACCESS

        _test_FintQ = _FintQ;
        _test_FextQ = _FextQ;

        _test_FintD = _FintD;
        _test_FextD = _FextD;

        _test_FintQhat = _FintQhat;
        _test_FextQhat = _FextQhat;

        _test_FintDhat = _FintDhat;
        _test_FextDhat = _FextDhat;

#endif

        if ( ( projection_type.compare( "l2_projection" ) == 0 ) ){
            //Assemble the micro force vector
            _FQ  = _FextQ;
            _FQ += _dense_BQhatQ.transpose( ) * _FextQhat;
            _FQ -= _FintQ;
            _FQ -= _dense_BQhatQ.transpose( ) * _FintQhat;
            _FQ -= _dense_BDhatQ.transpose( ) * _FintDhat;
    
            //Assemble the macro force vector
            _FD  = _FextD;
            _FD += _dense_BDhatD.transpose( ) * _FextDhat;
            _FD -= _FintD;
            _FD -= _dense_BQhatD.transpose( ) * _FintQhat;
            _FD -= _dense_BDhatD.transpose( ) * _FintDhat;

        }
        else if ( ( projection_type.compare( "direct_projection" ) == 0 ) ||
                  ( projection_type.compare( "averaged_l2_projection" ) == 0 ) ){

            //Assemble the micro force vector
            _FQ  = _FextQ;
            _FQ += _sparse_BQhatQ.transpose( ) * _FextQhat;
            _FQ -= _FintQ;
            _FQ -= _sparse_BQhatQ.transpose( ) * _FintQhat;
            _FQ -= _sparse_BDhatQ.transpose( ) * _FintDhat;
    
            //Assemble the macro force vector
            _FD  = _FextD;
            _FD += _sparse_BDhatD.transpose( ) * _FextDhat;
            _FD -= _FintD;
            _FD -= _sparse_BQhatD.transpose( ) * _FintQhat;
            _FD -= _sparse_BDhatD.transpose( ) * _FintDhat;

        }
        else{

            return new errorNode( __func__,
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

//        std::cout << "FQ:\n";
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
//        return new errorNode( __func__, "derp" );
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
    
                    errorOut result = new errorNode( __func__,
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
     
                        errorOut result = new errorNode( __func__,
                                                         "Error in the computation of the local gradient\n" );
                        result->addNext( error );
                        return result;
    
                    }
    
                    elementVolume += vectorTools::determinant( vectorTools::appendVectors( jacobian ), _dim, _dim ) * qpt->second;
    
                }
    
                //Compute the volumes of the overlapped micro domains for the macro cell
                auto microDomainVolumes = homogenizedVolumes.find( *macroCellID );
                if ( microDomainVolumes == homogenizedVolumes.end( ) ){
    
                    return new errorNode( __func__,
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

            return new errorNode( __func__,
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

                return new errorNode( __func__,
                                      "Micro node " + std::to_string( *nodeId ) + " not found in global to local map" );

            }

            auto previousMicroDisp = previousMicroDispDOFVector->find( *nodeId );

            if ( previousMicroDisp == previousMicroDispDOFVector->end( ) ){

                return new errorNode( __func__,
                                      "The micro node " + std::to_string( *nodeId ) +
                                      " is not found in the previous micro displacement dof vector" );

            }

            auto previousMicroVel = previousMicroVelocities->find( *nodeId );

            if ( previousMicroVel == previousMicroVelocities->end( ) ){

                return new errorNode( __func__,
                                      "The micro node " + std::to_string( *nodeId ) +
                                      " is not found in the previous micro velocities vector" );

            }

            auto previousMicroAccel = previousMicroAccelerations->find( *nodeId );

            if ( previousMicroAccel == previousMicroAccelerations->end( ) ){

                return new errorNode( __func__,
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

                return new errorNode( __func__,
                                      "Macro node " + std::to_string( *nodeId ) + " not found in global to local map" );

            }

            auto previousMacroDisp = previousMacroDispDOFVector->find( *nodeId );

            if ( previousMacroDisp == previousMacroDispDOFVector->end( ) ){

                return new errorNode( __func__,
                                      "The macro node " + std::to_string( *nodeId ) +
                                      " is not found in the previous macro displacement dof vector" );

            }

            auto previousMacroVel = previousMacroVelocities->find( *nodeId );

            if ( previousMacroVel == previousMacroVelocities->end( ) ){

                return new errorNode( __func__,
                                      "The macro node " + std::to_string( *nodeId ) +
                                      " is not found in the previous macro velocities vector" );

            }

            auto previousMacroAccel = previousMacroAccelerations->find( *nodeId );

            if ( previousMacroAccel == previousMacroAccelerations->end( ) ){

                return new errorNode( __func__,
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

                    DotDotDOF_t[ nMacroDispDOF * indexMap->second + i + microOffset ] = previousMacroAccel->second[ i ];

                }

            }

        }

        //Map the vectors to Eigen matrices
        Eigen::Map< Eigen::Matrix< floatType, -1,  1 > > _DOF( FreeDOF.data(), FreeDOF.size( ), 1 );
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > _DotDOF( DotDOF.data(), DotDOF.size( ), 1 );
        Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > _DotDotDOF_t( DotDotDOF_t.data(), DotDotDOF_t.size( ), 1 );
        Eigen::Map< Eigen::Matrix< floatType, -1,  1 > > _DotDotDOF_tp1( DotDotDOF_tp1.data(), DotDotDOF_tp1.size( ), 1 );

#ifdef TESTACCESS

        _test_DOF_t = _DOF;
        _test_DotDOF_t = _DotDOF;
        _test_DotDotDOF_t = _DotDotDOF_t;

#endif

        //Instantiate the QR solver
        std::cout << "Performing QR decomposition of the Free DOF LHS matrix\n";

        Eigen::MatrixXd RHS;
        if ( ( projection_type.compare( "l2_projection" ) == 0 ) ){

            SparseMatrix LHS( _sparse_MASS.rows( ), _sparse_MASS.cols( ) );
            LHS = _sparse_MASS;
            LHS += gamma * ( *dt ) * _sparse_DAMPING;

            RHS = _FORCE;
            RHS -= _sparse_DAMPING * ( _DotDOF + ( 1 - gamma ) * ( *dt ) * _DotDotDOF_t );

//            _DotDotDOF_tp1 = LHS.colPivHouseholderQr( ).solve( RHS );
            Eigen::SparseQR< SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
            solver.compute( LHS );
            _DotDotDOF_tp1 = solver.solve( RHS );
        }
        else if ( ( projection_type.compare( "direct_projection" ) == 0 ) ||
                  ( projection_type.compare( "averaged_l2_projection" ) == 0 ) ){

            SparseMatrix LHS( _sparse_MASS.rows( ), _sparse_MASS.cols( ) );
            LHS = _sparse_MASS;
            LHS += gamma * ( *dt ) * _sparse_DAMPING;
            LHS.makeCompressed( );

            RHS = _FORCE;
            RHS -= _sparse_DAMPING * ( _DotDOF + ( 1 - gamma ) * ( *dt ) * _DotDotDOF_t );

            Eigen::SparseQR< SparseMatrix, Eigen::COLAMDOrdering<int> > solver;
            solver.compute( LHS );
            _DotDotDOF_tp1 = solver.solve( RHS );

//            return new errorNode( __func__, "derp3" );
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
//            return new errorNode( __func__, "derp" );
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
//            return new errorNode( __func__, "derp2" );

        }
        else{

            return new errorNode( __func__, "Projection type " + projection_type + " not recognized" );

        }

        //Update the free degrees of freedom
        _DOF += ( *dt ) * _DotDOF + 0.5 * ( ( *dt ) * ( *dt ) ) * ( ( 1 - 2 * beta ) * _DotDotDOF_t + 2 * beta * _DotDotDOF_tp1 );

#ifdef TESTACCESS

        _test_DOF_tp1 = _DOF;
        _test_DotDotDOF_tp1 = _DotDotDOF_tp1;

#endif

        //Store the free degrees of freedom
        _updatedFreeMicroDispDOFValues = floatVector( FreeDOF.begin( ), FreeDOF.begin( ) + microOffset );
        _updatedFreeMacroDispDOFValues = floatVector( FreeDOF.begin( ) + microOffset, FreeDOF.end( ) );
        _freeDOFValuesUpdated = true;

        if ( updateGhostDOF ){

            errorOut error = projectDegreesOfFreedom( updateGhostDOF );

            if ( error ){
    
                errorOut result = new errorNode( __func__,
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

                errorOut result = new errorNode( __func__,
                                                 "Error when writing out the interpolation matrix N" );
                result->addNext( error );
                return result;

            }


        }

        errorOut error;

        //Save the projection matrices
        error = writeSparseMatrixToXDMF( _centerOfMassN, "centerOfMassInterpolator", reference_filename, domain, grid );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error when writing out the center of mass projection matrix" );
            result->addNext( error );
            return result;

        }

        error = writeDenseMatrixToXDMF( _centerOfMassProjector, "centerOfMassProjector", reference_filename, domain, grid );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error when writing out the center of mass projection matrix" );
            result->addNext( error );
            return result;

        }

        std::string projectionType = couplingInitialization[ "projection_type" ].as< std::string >( );
        if ( ( projectionType.compare( "l2_projection" ) == 0 ) ){

            //Initialize the projection matrix attributes
            shared_ptr< XdmfAttribute > BQhatQ = XdmfAttribute::New( );
            shared_ptr< XdmfAttribute > BQhatD = XdmfAttribute::New( );
            shared_ptr< XdmfAttribute > BDhatQ = XdmfAttribute::New( );
            shared_ptr< XdmfAttribute > BDhatD = XdmfAttribute::New( );

            //Write the matrix type to the output file
            shared_ptr< XdmfInformation > projectionType = XdmfInformation::New( "EIGEN_MATRIX_TYPE", "DENSE" );
            grid->insert( projectionType );

            //Write BQhatQ
            error = writeDenseMatrixToXDMF( _dense_BQhatQ, "BQhatQ", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error when writing out BQhatQ" );
                result->addNext( error );
                return result;

            }

            //Write BQhatD
            error = writeDenseMatrixToXDMF( _dense_BQhatD, "BQhatD", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error when writing out BQhatD" );
                result->addNext( error );
                return result;

            }

            //Write BDhatQ
            error = writeDenseMatrixToXDMF( _dense_BDhatQ, "BDhatQ", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error when writing out BDhatQ" );
                result->addNext( error );
                return result;

            }

            //Write BDhatD
            error = writeDenseMatrixToXDMF( _dense_BDhatD, "BDhatD", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error when writing out BDhatD" );
                result->addNext( error );
                return result;

            }

        }
        else if ( ( projectionType.compare( "direct_projection" ) == 0 ) || ( projectionType.compare( "averaged_l2_projection" ) == 0 ) ){

            //Write the matrix type to the output file
            shared_ptr< XdmfInformation > projectionType = XdmfInformation::New( "EIGEN_MATRIX_TYPE", "SPARSE" );
            domain->insert( projectionType );

            //Write BQhatD
            error = writeSparseMatrixToXDMF( _sparse_BQhatD, "BQhatD", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error when writing out BQhatD" );
                result->addNext( error );
                return result;

            }

            //Write BDhatD
            error = writeSparseMatrixToXDMF( _sparse_BDhatD, "BDhatD", reference_filename, domain, grid );
    
            if ( error ){
    
                errorOut result = new errorNode( __func__,
                                                 "Error when writing out BDhatD" );
                result->addNext( error );
                return result;
    
            }

            //Write BQhatQ
            error = writeSparseMatrixToXDMF( _sparse_BQhatQ, "BQhatQ", reference_filename, domain, grid );

            if ( error ){

                errorOut result = new errorNode( __func__,
                                                 "Error when writing out BQhatQ" );
                result->addNext( error );
                return result;

            }

            //Write BDhatQ
            error = writeSparseMatrixToXDMF( _sparse_BDhatQ, "BDhatQ", reference_filename, domain, grid );
    
            if ( error ){
    
                errorOut result = new errorNode( __func__,
                                                 "Error when writing out BDhatQ" );
                result->addNext( error );
                return result;
    
            }

        }
        else if ( projectionType.compare( "arlequin" ) == 0 ){

            return NULL;

        }
        else{

            return new errorNode( __func__,
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
            return new errorNode( __func__, outstr );

        }

        shared_ptr< XdmfInformation > shapeInfo = _A->getInformation( 0 );

        if ( !shapeInfo ){

            std::string outstr = "There is no information defined for the matrix " + matrixName;
            return new errorNode( __func__, outstr );

        }

        if ( shapeInfo->getKey( ).compare( matrixName + "_shape" ) != 0 ){

            std::string outstr = matrixName;
            outstr += "_shape is not in the information key";
            return new errorNode( __func__, outstr );

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
            return new errorNode( __func__, outstr );

        }

        shared_ptr< XdmfInformation > shapeInfo = _A->getInformation( 0 );

        if ( !shapeInfo ){

            std::string outstr = "There is no information defined for the SparseMatrix " + matrixName;
            return new errorNode( __func__, outstr );

        }

        if ( shapeInfo->getKey( ).compare( matrixName + "_shape" ) != 0 ){

            std::string outstr = matrixName;
            outstr += "_shape is not in the information key";
            return new errorNode( __func__, outstr );

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
            return new errorNode( __func__, outstr );

        }

        shared_ptr< XdmfAttribute > _cols = grid->getAttribute( matrixName + "_cols" );

        if ( !_cols ){

            std::string outstr = matrixName;
            outstr += "_cols attribute is not found";
            return new errorNode( __func__, outstr );

        }

        shared_ptr< XdmfAttribute > _vals = grid->getAttribute( matrixName + "_values" );

        if ( !_cols ){

            std::string outstr = matrixName;
            outstr += "_values attribute is not found";
            return new errorNode( __func__, outstr );

        }

        if ( ( _rows->getSize( ) != _cols->getSize( ) ) && ( _rows->getSize( ) != _vals->getSize( ) ) ){

            std::string outstr = matrixName;
            outstr += " attributes rows, cols, and values don't have consistent sizes";
            return new errorNode( __func__, outstr );

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

            errorOut result = new errorNode( __func__,
                                             "Error in extracting the center of mass interpolation" );
            result->addNext( error );
            return result;

        }

        error = readDenseMatrixFromXDMF( _readGrid, "centerOfMassProjector", _centerOfMassProjector );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error in extracting the center of mass projector" );
            result->addNext( error );
            return result;

        }

        if ( ( projectionType.compare( "l2_projection" ) ) ) {

            error = readDenseMatrixFromXDMF( _readGrid, "BQhatQ", _dense_BQhatQ );

            if ( error ){

                errorOut result
                    = new errorNode( __func__, "Error when extracting BQhatQ from the XDMF file" );
                result->addNext( error );
                return result;

            }

            error = readDenseMatrixFromXDMF( _readGrid, "BQhatD", _dense_BQhatD );

            if ( error ){

                errorOut result
                    = new errorNode( __func__, "Error when extracting BQhatD from the XDMF file" );
                result->addNext( error );
                return result;

            }

            error = readDenseMatrixFromXDMF( _readGrid, "BDhatQ", _dense_BDhatQ );

            if ( error ){

                errorOut result
                    = new errorNode( __func__, "Error when extracting BDhatQ from the XDMF file" );
                result->addNext( error );
                return result;

            }

            error = readDenseMatrixFromXDMF( _readGrid, "BDhatD", _dense_BDhatD );

            if ( error ){

                errorOut result
                    = new errorNode( __func__, "Error when extracting BDhatD from the XDMF file" );
                result->addNext( error );
                return result;

            }

        }
        else if ( ( projectionType.compare( "direct_projection" ) ) || ( projectionType.compare( "averaged_l2_projection" ) ) ){

            error = readSparseMatrixFromXDMF( _readGrid, "BQhatQ", _sparse_BQhatQ );

            if ( error ){

                errorOut result
                    = new errorNode( __func__, "Error when extracting BQhatQ from the XDMF file" );
                result->addNext( error );
                return result;

            }

            error = readSparseMatrixFromXDMF( _readGrid, "BQhatD", _sparse_BQhatD );

            if ( error ){

                errorOut result
                    = new errorNode( __func__, "Error when extracting BQhatD from the XDMF file" );
                result->addNext( error );
                return result;

            }

            error = readSparseMatrixFromXDMF( _readGrid, "BDhatQ", _sparse_BDhatQ );

            if ( error ){

                errorOut result
                    = new errorNode( __func__, "Error when extracting BDhatQ from the XDMF file" );
                result->addNext( error );
                return result;

            }

            error = readSparseMatrixFromXDMF( _readGrid, "BDhatD", _sparse_BDhatD );

            if ( error ){

                errorOut result
                    = new errorNode( __func__, "Error when extracting BDhatD from the XDMF file" );
                result->addNext( error );
                return result;

            }

        }
        else{

            return new errorNode( __func__, "Not implemented" );

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

            errorOut result = new errorNode( __func__,
                                             "Error when initializing the writer" );
            result->addNext( writer->_error );
            return result;

        }

        const floatType *time = _inputProcessor.getMacroTime( );

        uIntType increment;
        errorOut error = writer->initializeIncrement( *time, _currentReferenceOutputIncrement, collectionNumber, increment );

        if ( error ){

            errorOut result = new errorNode( __func__,
                                             "Error in the initialization of the increment for the homogenized output" );
            result->addNext( error );
            return result;

        }

        //Write the mesh data to file. This increment references a previous increment for the mesh data so we don't need to do it again
        error = writer->writeIncrementMeshData( increment, collectionNumber, { }, { { } }, { }, { }, { }, { { } }, { }, { } );

        if ( error ){

            errorOut result = new errorNode( __func__,
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

        //Loop over the macro displacement degrees of freedom
        stringVector node_dof_names( _dim + _dim * _dim );
        for ( unsigned int dof_i = 0; dof_i < ( _dim + _dim * _dim ); dof_i++ ){

            node_dof_names[ dof_i ] = "NODE_DOF_" + std::to_string( dof_i + 1 );

        }

        
        error = writer->writeSolutionData( increment, collectionNumber, node_dof_names, "Node", projectedMacroDisplacement );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in writing the nodal DOF solution data" );
            result->addNext( error );
            return result;

        }

        //Loop over the quadrature points
        for ( unsigned int qp = 0; qp < maxQP; qp++ ){

            floatVector densityOut(                                      cellIds.size( ), 0 );
            floatVector bodyForceOut(                             _dim * cellIds.size( ), 0 );
            floatVector positionOut(                              _dim * cellIds.size( ), 0 );
            floatVector accelerationsOut(                         _dim * cellIds.size( ), 0 );
            floatVector microInertiasOut(                  _dim * _dim * cellIds.size( ), 0 );
            floatVector bodyCouplesOut(                    _dim * _dim * cellIds.size( ), 0 );
            floatVector microSpinInertiasOut(              _dim * _dim * cellIds.size( ), 0 );
            floatVector symmetricMicroStressOut(           _dim * _dim * cellIds.size( ), 0 );
            floatVector cauchyStressOut(                   _dim * _dim * cellIds.size( ), 0 );
            floatVector higherOrderStressOut(       _dim * _dim * _dim * cellIds.size( ), 0 );
            floatVector dofValuesOut(           ( _dim + _dim * _dim ) * cellIds.size( ), 0 );
            floatVector dofGradientsOut( _dim * ( _dim + _dim * _dim ) * cellIds.size( ), 0 );

            floatVector elementDofValuesOut(                      ( _dim + _dim * _dim ), 0 );
            floatVector elementDofGradientsOut(            _dim * ( _dim + _dim * _dim ), 0 );

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

                if ( _inputProcessor.isFiltering( ) ){
        
                    std::unordered_map< uIntType, floatVector > projectedMacroDisplacements;
        
                    auto cellConnectivity = _inputProcessor.getMacroNodeReferenceConnectivity( )->find( *cellId );
        
                    if ( cellConnectivity == _inputProcessor.getMacroNodeReferenceConnectivity( )->end( ) ){
        
                        return new errorNode( __func__, "Macro cell " + std::to_string( *cellId ) + " not found in connectivity map" );
        
                    }
    
                    auto macroGlobalToLocalDOFMap = _inputProcessor.getMacroGlobalToLocalDOFMap( );
    
                    unsigned int nMacroDOF = _dim + _dim * _dim;
    
                    for ( auto node = cellConnectivity->second.begin( ) + 1; node != cellConnectivity->second.end( ); node++ ){
    
                        auto local_node = macroGlobalToLocalDOFMap->find( *node );
    
                        if ( local_node == macroGlobalToLocalDOFMap->end( ) ){
    
                            return new errorNode( __func__, "Micro node " + std::to_string( *node) + " not found in the global to local map" );
    
                        }
    
                        projectedMacroDisplacements.emplace( *node, floatVector( _projected_ghost_macro_displacement.begin( ) + nMacroDOF * local_node->second,
                                                                                 _projected_ghost_macro_displacement.begin( ) + nMacroDOF * local_node->second + 3 ) );
    
                    }

                    error = buildMacroDomainElement( *cellId,
                                                     *_inputProcessor.getMacroNodeReferencePositions( ),
                                                     projectedMacroDisplacements,
                                                     *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                     element );

                }
                else{

                    error = buildMacroDomainElement( *cellId,
                                                     *_inputProcessor.getMacroNodeReferencePositions( ),
                                                     *_inputProcessor.getMacroDisplacements( ),
                                                     *_inputProcessor.getMacroNodeReferenceConnectivity( ),
                                                     element );

                }

                //Compute the quadrature point position
                floatVector qpPosition;
                error = element->interpolate( element->nodes, element->qrule[ qp ].first, qpPosition );
                
                if ( error ){
                
                    errorOut result = new errorNode( __func__, "Error in the computation of the gauss point location" );
                    result->addNext( error );
                    return result;
                
                }
                
                for ( unsigned int i = 0; i < _dim; i++ ){
                
                    positionOut[ _dim * index + i ] = qpPosition[ i ];
                
                }

                //Build the DOF vector
                floatMatrix dofMatrix( element->qrule.size( ), floatVector( _dim + _dim * _dim, 0 ) );

                for ( auto node = element->global_node_ids.begin( ); node != element->global_node_ids.end( ); node++ ){

                    auto localNode = _inputProcessor.getMacroGlobalToLocalDOFMap( )->find( *node );

                    if ( localNode == _inputProcessor.getMacroGlobalToLocalDOFMap( )->end( ) ){

                        return new errorNode( __func__,
                                              "Error in finding the global node " + std::to_string( *node ) +
                                              " in the macro global to local DOF map" );

                    }

                    dofMatrix[ node - element->global_node_ids.begin( ) ] =
                        floatVector( projectedMacroDisplacement.begin( ) + ( _dim + _dim * _dim ) * localNode->second,
                                     projectedMacroDisplacement.begin( ) + ( _dim + _dim * _dim ) * ( localNode->second + 1 ) );

                }

                error = element->interpolate( dofMatrix, element->qrule[ qp ].first, elementDofValuesOut );

                if ( error ){

                    errorOut result = new errorNode( __func__,
                                                     "Error in the interpolation of the DOF values" );
                    result->addNext( error );
                    return result;

                }

                // Load the Dof values
                for ( unsigned int i = 0; i < elementDofValuesOut.size( ); i++ ){

                    dofValuesOut[ elementDofValuesOut.size( ) * index + i ] = elementDofValuesOut[ i ];

                }

                floatMatrix qptDOFGradient;

                error = element->get_global_gradient( dofMatrix, element->qrule[ qp ].first, element->reference_nodes, qptDOFGradient );

                if ( error ){

                    errorOut result = new errorNode( __func__,
                                                     "Error in the interpolation of the DOF values" );
                    result->addNext( error );
                    return result;

                }

                elementDofGradientsOut = vectorTools::appendVectors( qptDOFGradient );

                // Load the gradients of the DOF values
                for ( unsigned int i = 0; i < elementDofGradientsOut.size( ); i++ ){

                    dofGradientsOut[ elementDofGradientsOut.size( ) * index + i ] = elementDofGradientsOut[ i ];

                }

            }

            //Write quadrature point information to file
            stringVector outputNames = { "density_" + std::to_string( qp ) };
            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", densityOut );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error in outputting the density" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim );

            for ( unsigned int i = 0; i < _dim; i++ ){

                outputNames[ i ] = "acceleration_" + std::to_string( i + 1 ) + "_" + std::to_string( qp );

            }

            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", accelerationsOut );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error in outputting the acceleration" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim );
            for ( unsigned int i = 0; i < _dim; i++ ){

                outputNames[ i ] = "body_force_" + std::to_string( i + 1 ) + "_" + std::to_string( qp );

            }
            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", bodyForceOut );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error in outputting the body force" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim );
            for ( unsigned int i = 0; i < _dim; i++ ){
            
                outputNames[ i ] = "position_" + std::to_string( i + 1 ) + "_" + std::to_string( qp );
            
            }
            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", positionOut );
            
            if ( error ){
            
                errorOut result = new errorNode( __func__, "Error in outputting the gauss point position" );
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

                errorOut result = new errorNode( __func__, "Error in outputting the micro inertias" );
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

                errorOut result = new errorNode( __func__, "Error in outputting the body couples" );
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

                errorOut result = new errorNode( __func__, "Error in outputting the body couples" );
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

                errorOut result = new errorNode( __func__, "Error in outputting the symmetric micro stress" );
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

                errorOut result = new errorNode( __func__, "Error in outputting the cauchy stress" );
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

                errorOut result = new errorNode( __func__, "Error in outputting the higher order stress" );
                result->addNext( error );
                return result;

            }

            outputNames = stringVector( _dim + _dim * _dim );
            for ( unsigned int i = 0; i < ( _dim + _dim * _dim ); i++ ){

                    outputNames[ i ] = "DOF_" + std::to_string( i + 1 ) + "_" + std::to_string( qp );

            }
            error = writer->writeSolutionData( increment, collectionNumber, outputNames, "Cell", dofValuesOut );

            if ( error ){

                errorOut result = new errorNode( __func__, "Error in outputting the degree of freedom values" );
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

                errorOut result = new errorNode( __func__, "Error in outputting the degree of freedom values" );
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

            errorOut result = new errorNode( __func__, "Error in construction of writer" );
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

                return new errorNode( __func__, "The macro node " + std::to_string( node->first ) +
                                      " was not found in the reference positions vector" );

            }

            auto displacement = macroDisplacements->find( node->first );

            if ( displacement == macroDisplacements->end( ) ){

                return new errorNode( __func__, "The macro node " + std::to_string( node->first ) +
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

                return new errorNode( __func__,
                                      "The free macro node " + std::to_string( *node ) + " was not found in the macro global to local DOF map" );

            }

            nodeSets[ 0 ][ node - freeMacroNodeIds->begin( ) ] = localNode->second;

        }

        for ( auto node = ghostMacroNodeIds->begin( ); node != ghostMacroNodeIds->end( ); node++ ){

            auto localNode = macroGlobalToLocalDOFMap->find( *node );

            if ( localNode == macroGlobalToLocalDOFMap->end( ) ){

                return new errorNode( __func__,
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

                errorOut result = new errorNode( __func__,
                                                 "Error in construction of the micromorphic element" );
                result->addNext( error );
                return result;

            }

            //Get the XDMF cell type
            auto elementConnectivity = macroNodeReferenceConnectivity->find( *cell );

            if ( elementConnectivity == macroNodeReferenceConnectivity->end( ) ){

                return new errorNode( __func__,
                                      "Macro cell " + std::to_string( *cell ) + " was not found in the macro mesh connectivity" );

            }

            uIntType cellType = elementConnectivity->second[ 0 ];

            uIntVector localNodeIds( element->global_node_ids.size( ) + 1 );
            localNodeIds[ 0 ] = cellType;

            for ( auto gN = element->global_node_ids.begin( ); gN != element->global_node_ids.end( ); gN++ ){

                auto node = macroGlobalToLocalDOFMap->find( *gN );

                if ( node == macroGlobalToLocalDOFMap->end( ) ){

                    return new errorNode( __func__,
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

            errorOut result = new errorNode( __func__,
                                             "Error in initialization of the increment in the writer" );
            result->addNext( error );
            return result;

        }

        error = writer->writeIncrementMeshData( _currentReferenceOutputIncrement, collectionNumber,
                                                nodeIds, nodeSets, nodeSetNames, nodePositions,
                                                elementIds, elementSets, elementSetNames, connectivity );

        if ( error ){

            errorOut result = new errorNode( __func__,
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

            return new errorNode( __func__,
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

            errorOut result = new errorNode( __func__, "Error in construction of writer" );
            result->addNext( macro_writer->_error );
            return result;

        }

        //Form the micro writer object
        std::shared_ptr< dataFileInterface::dataFileBase > micro_writer
            = dataFileInterface::dataFileBase( micro_config ).create( );

        if ( micro_writer->_error ){

            errorOut result = new errorNode( __func__, "Error in construction of writer" );
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

            errorOut result = new errorNode( __func__, "Error in initializating the increment of the macro output file" );
            result->addNext( error );
            return result;

        }

        error = macro_writer->writeScalarSolutionData( increment, collectionNumber, "updated_DOF", "Node", outputDOF ); 
        
        if ( error ){

            errorOut result = new errorNode( __func__, "Error in outputting the updated macro DOF to the output file" );
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

            errorOut result = new errorNode( __func__, "Error in outputting the updated node ids to the micro output file" );
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

            errorOut result = new errorNode( __func__, "Error in initializating the increment of the micro output file" );
            result->addNext( error );
            return result;

        }

        error = micro_writer->writeScalarSolutionData( increment, collectionNumber, "updated_DOF", "Node", outputDOF ); 
        
        if ( error ){

            errorOut result = new errorNode( __func__, "Error in outputting the updated micro DOF to the output file" );
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

            errorOut result = new errorNode( __func__, "Error in outputting the updated node ids to the micro output file" );
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

    floatVector overlapCoupling::getMicroAugmentedLagrangianForce( ){
        /*!
         * Get the micro augmented lagrangian force vector
         */

        return FALQ;
    }

    floatVector overlapCoupling::getMacroAugmentedLagrangianForce( ){
        /*!
         * Get the macro augmented lagrangian force vector
         */

        return FALD;
    }

    floatVector overlapCoupling::getUpdatedMacroDisplacementDOF( ){
        /*!
         * Return a copy of the macro displacement dof vector
         */

        return vectorTools::appendVectors( { _updatedFreeMacroDispDOFValues, _projected_ghost_macro_displacement } );
    }

    errorOut runOverlapCoupling( const std::string &filename,
                                 DOFMap &microGlobalLocalNodeMap, floatVector &updatedMicroDisplacementDOF,
                                 floatVector &Lagrangian_FALQ,
                                 DOFMap &macroGlobalLocalNodeMap, floatVector &updatedMacroDisplacementDOF,
                                 floatVector &Lagrangian_FALD
                               ){
        /*!
         * Run the overlap coupling method
         *
         * :param const std::string &filename: The name of the input YAML file
         * :param DOFMap &microGlobalLocalNodeMap: The map from global to local node numbers for the micro nodes
         * :param floatVector &updatedMicroDisplacementDOF: The updated micro displacement degrees of freedom
         * :param floatVector &Lagrangian_FALQ: The augmented lagrangian force for the micro domain
         * :param DOFMap &macroGlobalLocalNodeMap: The map from global to local node numbers for the macro nodes
         * :param floatVector &updatedMacroDisplacementDOF: The updated macro displacement degrees of freedom
         * :param floatVector &Lagrangian_FALD: The augmented lagrangian force for the macro domain
         */

        //Construct the overlap coupling object
        overlapCoupling oc( filename );

        if ( oc.getConstructorError( ) ){

            errorOut result
                = new errorNode( __func__, "Error in construction of overlapCoupling object" );

            result->addNext( oc.getConstructorError( ) );

            return result;

        }

        //Initialize the overlap coupling object
        std::cerr << "INITIALIZE COUPLING\n";
        errorOut error = oc.initializeCoupling( );

        if ( error ){
    
            errorOut result
                = new errorNode( __func__, "Error in the initialization of the overlapCoupling object" );
    
            result->addNext( error );

            return result;
    
        }

        //Process the final increments of both the macro and micro-scales
        std::cerr << "PROCESS LAST INCREMENTS\n";
        error = oc.processLastIncrements( );
    
        if ( error ){
    
            errorOut result
                = new errorNode( __func__, "Error in processing the data" );
    
            result->addNext( error );

            return result;
    
        }

        //Return the updated DOF values
        std::cerr << "RETURN DOF VALUES\n";
        microGlobalLocalNodeMap = oc.getMicroGlobalLocalNodeMap( );
        updatedMicroDisplacementDOF = oc.getUpdatedMicroDisplacementDOF( );
        Lagrangian_FALQ = oc.getMicroAugmentedLagrangianForce( );

        macroGlobalLocalNodeMap = oc.getMacroGlobalLocalNodeMap( );
        updatedMacroDisplacementDOF = oc.getUpdatedMacroDisplacementDOF( );
        Lagrangian_FALD = oc.getMacroAugmentedLagrangianForce( );

        return NULL;

    }

    const cellDomainFloatVectorMap* overlapCoupling::getReferenceCellDomainCenterOfMassShapeFunctions( ){
        /*!
         * Get a constant reference to the shapefunctions of the centers of mass at each macro cell
         */

        return &_referenceCellDomainCenterOfMassShapefunctions;

    }

    const floatVector* overlapCoupling::getUpdatedFreeMacroDispDOFValues( ){
        /*!
         * Get a constant reference to the macro displacement DOF values
         */

        return &_updatedFreeMacroDispDOFValues;
    }

    const floatVector* overlapCoupling::getUpdatedFreeMicroDispDOFValues( ){
        /*!
         * Get a constant reference to the micro displacement DOF values
         */

        return &_updatedFreeMicroDispDOFValues;
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

    const cellDomainFloatMap* overlapCoupling::getHomogenizedSurfaceAreas( ){
        /*!
         * Get the homogenized surface areas
         */

        return &homogenizedSurfaceAreas;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getHomogenizedSurfaceRegionAreas( ){
        /*!
         * Get the homogenized surface region areas
         */

        return &homogenizedSurfaceRegionAreas;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getHomogenizedSurfaceRegionCentersOfMass( ){
        /*!
         * Get the homogenized surface region centers of mass
         */

        return &homogenizedSurfaceRegionCentersOfMass;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getHomogenizedSurfaceRegionProjectedLocalCentersOfMass( ){
        /*!
         * Get the homogenized projected local centers of mass
         */

        return &homogenizedSurfaceRegionProjectedLocalCentersOfMass;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getHomogenizedSurfaceRegionProjectedCentersOfMass( ){
        /*!
         * Get the homogenized projected centers of mass
         */

        return &homogenizedSurfaceRegionProjectedCentersOfMass;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getHomogenizedSurfaceRegionTractions( ){
        /*!
         * Get the homogenized surface region tractions
         */

        return &homogenizedSurfaceRegionTractions;
    }

    const cellDomainFloatVectorMap* overlapCoupling::getHomogenizedSurfaceRegionCouples( ){
        /*!
         * Get the homogenized surface region couples
         */

        return &homogenizedSurfaceRegionCouples;
    }

    const cellDomainUIntVectorMap* overlapCoupling::getCellDomainMacroSurfaces( ){
        /*!
         * Get the cell domain macro surfaces
         */

        return &cellDomainMacroSurfaces;
    }

    const cellFloatVectorMap* overlapCoupling::getExternalForcesAtNodes( ){
        /*!
         * Get the external forces acting on the nodes
         */

        return &externalForcesAtNodes;
    }

    const cellFloatVectorMap* overlapCoupling::getExternalCouplesAtNodes( ){
        /*!
         * Get the external couples acting on the cells' nodes
         */

        return &externalCouplesAtNodes;
    }

    const std::unordered_map< uIntType, floatVector > *overlapCoupling::getQuadraturePointCauchyStress( ){
        /*!
         * Get a constant reference to the cauchy stress at the quadrature points
         */

        return &quadraturePointCauchyStress;
    }

    const std::unordered_map< uIntType, floatVector > *overlapCoupling::getQuadraturePointHigherOrderStress( ){
        /*!
         * Get a constant reference to the higher order stress at the quadrature points
         */

        return &quadraturePointHigherOrderStress;
    }

    const std::unordered_map< uIntType, floatVector > *overlapCoupling::getQuadraturePointDensities( ){
        /*!
         * Get a constant reference to the densities at the quadrature points
         */

        return &quadraturePointDensities;
    }

    const std::unordered_map< uIntType, floatVector > *overlapCoupling::getQuadraturePointBodyForce( ){
        /*!
         * Get a constant reference to the body force at the quadrature points
         */

        return &quadraturePointBodyForce;
    }

    const std::unordered_map< uIntType, floatVector > *overlapCoupling::getQuadraturePointAccelerations( ){
        /*!
         * Get a constant reference to the accelerations at the quadrature points
         */

        return &quadraturePointAccelerations;
    }

    const std::unordered_map< uIntType, floatVector > *overlapCoupling::getQuadraturePointMicroInertias( ){
        /*!
         * Get a constant reference to the micro inertias at the quadrature points
         */

        return &quadraturePointMicroInertias;
    }

    const std::unordered_map< uIntType, floatVector > *overlapCoupling::getQuadraturePointBodyCouples( ){
        /*!
         * Get a constant reference to the body couples at the quadrature points
         */

        return &quadraturePointBodyCouples;
    }

    const std::unordered_map< uIntType, floatVector > *overlapCoupling::getQuadraturePointMicroSpinInertias( ){
        /*!
         * Get a constant reference to the micro spin inertias at the quadrature points
         */

        return &quadraturePointMicroSpinInertias;
    }

    const std::unordered_map< uIntType, floatVector > *overlapCoupling::getQuadraturePointSymmetricMicroStress( ){
        /*!
         * Get a constant reference to the symmetric micro stress at the quadrature points
         */

        return &quadraturePointSymmetricMicroStress;
    }

    const Eigen::MatrixXd *overlapCoupling::getHomogenizedFEXT( ){
        /*!
         * Get a constant reference to the homogenized external force vector
         */

        return &homogenizedFEXT; 
    }

    const Eigen::MatrixXd *overlapCoupling::getHomogenizedFINT( ){
        /*!
         * Get a constant reference to the homogenized external force vector
         */

        return &homogenizedFINT; 
    }

    const SparseMatrix *overlapCoupling::getHomogenizedMassMatrix( ){
        /*!
         * Get a constant reference to the homogenized mass matrix
         */

        return &homogenizedMassMatrix;
    }

    const SparseMatrix *overlapCoupling::getFreeMicromorphicMassMatrix( ){
        /*!
         * Get a constant reference to the free micromorphic mass matrix
         */

        return &freeMicromorphicMassMatrix;
    }

    const SparseMatrix *overlapCoupling::getMass( ){
        /*!
         * Get the coupled mass matrix
         */

        return &_sparse_MASS;
    }

    const SparseMatrix *overlapCoupling::getDamping( ){
        /*!
         * Get the coupled damping matrix
         */

        return &_sparse_DAMPING;
    }

    const Eigen::MatrixXd *overlapCoupling::getFORCE( ){
        /*!
         * Get a constant reference to the total force vector
         */

        return &_FORCE;
    }

    const std::unordered_map< uIntType, floatType > *overlapCoupling::getArlequinMicroWeightingFactors( ){
        /*!
         * Get a constant reference to the Arlequin micro weighting factors
         */

        return &arlequinMicroWeightingFactors;

    }

    const SparseMatrix *overlapCoupling::getMD( ){
        /*!
         * Get a constant reference to the micromoprhic mass matrix
         */

        return &_MD;
    }

    const floatVector *overlapCoupling::getFQ( ){
        /*!
         * Get a constant reference to the micro-force vector
         */

        return &FQ;
    }

    const floatVector *overlapCoupling::getFD( ){
        /*!
         * Get a constant reference to the macro-force vector
         */

        return &FD;
    }

#endif

}
