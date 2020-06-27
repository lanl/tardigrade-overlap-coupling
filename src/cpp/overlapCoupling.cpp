/*!============================================================================
|                               overlapCoupling                               |
===============================================================================
| The implementation of the Overlap coupling method                           |
=============================================================================*/

#include<overlapCoupling.h>

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

    errorOut overlapCoupling::processIncrement( const unsigned int &increment ){
        /*!
         * Process the indicated increment
         *
         * :param const unsigned int &increment: The increment to process
         */

        //Initialize the input processor
        errorOut error = _inputProcessor.initializeIncrement( increment );
        if ( error ){
            errorOut result = new errorNode( "processIncrement", "Error in initialization of the input processor" );
            result->addNext( error );
            return result;
        }

        //Compute the centers of mass of the free and ghost domains
        error = computeIncrementCentersOfMass( increment,
                                               _freeMicroDomainMasses, _ghostMicroDomainMasses,
                                               _freeMicroDomainCentersOfMass, _ghostMicroDomainCentersOfMass );
        if ( error ){
            errorOut result = new errorNode( "processIncrement", "Error in computation of the domain centers of mass" );
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
            error = setReferenceStateFromIncrement( 0 );
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

    errorOut overlapCoupling::setReferenceStateFromIncrement( const unsigned int &increment ){
        /*!
         * Set the reference state from the indicated increment
         *
         * :param const unsigned int &increment: The increment at which to set the reference state
         */

        //Compute the centers of mass and the total mass
        errorOut error = computeIncrementCentersOfMass( increment, _referenceFreeMicroDomainMasses, _referenceGhostMicroDomainMasses,
                                                        _referenceFreeMicroDomainCentersOfMass, _referenceGhostMicroDomainCentersOfMass );

        if ( error ){
            errorOut result = new errorNode( "setReferenceFromIncrement",
                                             "Error in computing the center of mass of increment " + std::to_string( increment ) );
            result->addNext( error );
            return result;
        }

        //Compute the shape functions at the centers of mass
        error = computeShapeFunctionsAtReferenceCentersOfMass( );

        if ( error ){

            errorOut result = new errorNode( "setReferenceStateFromIncrement",
                                             "Error in computing the shape functions at the centers of mass of increment "
                                             + std::to_string( increment ) );
            result->addNext( error );
            return result;

        }

        return NULL;
    }

    errorOut overlapCoupling::computeIncrementCentersOfMass( const unsigned int increment,
                                                             floatVector &freeDomainMass, floatVector &ghostDomainMass,
                                                             floatVector &freeDomainCM, floatVector &ghostDomainCM ){
        /*!
         * Compute the centers of mass for an increment. Also computes the mass of the micro-scale domains
         *
         * :param const unsigned int increment: The increment at which to compute the centers of mass
         * :param floatVector &freeDomainMass: The mass of the free domains
         * :param floatVector &ghostDomainMass: The mass of the ghost domains
         * :param floatVector &freeDomainCM: The center of mass of the free domains
         * :param floatVector &ghostDomainCM: The center of mass of the ghost domains
         */

        //Compute the centers of mass of each domain
        errorOut error = _inputProcessor.initializeIncrement( increment );
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

            error = _inputProcessor._microscale->getDomainNodes( increment, *name, domainNodes );

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

            error = _inputProcessor._microscale->getDomainNodes( increment, *name, domainNodes );

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
                                  "The points vector is inconsistent with the dimension" );

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

}
