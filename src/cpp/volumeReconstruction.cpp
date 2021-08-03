/*=============================================================================
|                             volumeReconstruction                            |
===============================================================================
| Tools to reconstruct volume information from pointsets                      |
=============================================================================*/

#include<volumeReconstruction.h>
#include<solver_tools.h>

#include "Xdmf.hpp"

#include "XdmfDomain.hpp"
#include "XdmfInformation.hpp"
#include "XdmfReader.hpp"
#include "XdmfWriter.hpp"
#include "XdmfHDF5Writer.hpp"
#include "XdmfGeometry.hpp"
#include "XdmfUnstructuredGrid.hpp"
#include "XdmfGridCollection.hpp"
#include "XdmfGridCollectionType.hpp"

namespace volumeReconstruction{

    KDNode::KDNode( ){
        /*!
         * The default constructor for the KD tree
         */

        return;
    }

    KDNode::KDNode( const floatVector *points, const uIntVector &ownedIndices,
                    const uIntType &depth, const uIntType &dim ){

        /*!
         * The constructor for the KD tree
         *
         * :param const floatVector *points: The points to be sorted into the tree
         * :param const uIntVector &ownedIndices: The indices of the points which this
         *     node is dividing. These are indices to the first values of the 
         *     points in the 'points' vector.
         * :param const uIntType &depth: The depth of the current node
         * :param const uIntType &dim: The dimension of the points
         */

        //Set the points
        _points = points;

        //Set the depth of the node
        _depth = depth;

        //Determine if this is the last node in the chain
        if ( ownedIndices.size( ) == 1 ){
            _index = ownedIndices[ 0 ];
            _axis = 0;
            return;
        }

        //Determine the bounding box of the points
        floatVector lowerBound( _points->begin( ) + ownedIndices[ 0 ],
                                _points->begin( ) + ownedIndices[ 0 ] + dim );

        floatVector upperBound( _points->begin( ) + ownedIndices[ 0 ],
                                _points->begin( ) + ownedIndices[ 0 ] + dim );

        for ( auto index = ownedIndices.begin( ) + 1; index != ownedIndices.end( ); index++ ){

            //Compare the current point against the bounds
            for ( uIntType i = 0; i < dim; i++ ){
                if ( ( *_points )[ *index + i ] > upperBound[ i ] ){
                    upperBound[ i ] = ( *_points )[ *index + i ];
                }
                else if ( ( *_points )[ *index + i ] < lowerBound[ i ] ){
                    lowerBound[ i ] = ( *_points )[ *index + i ];
                }
            }

        }

        //Determine which dimension has the largest variation. That is the axis
        //we want to split our domain by
        floatVector delta = upperBound - lowerBound;
        _axis = 0;
        floatType deltaMax = delta[ _axis ];

        uIntType a = 0;

        for ( auto v = delta.begin( ) + 1; v != delta.end( ); v++ ){

            if ( deltaMax < *v ){

                _axis = a;
                deltaMax = *v;

            }

            a++;

        }

        //Accumulate the dimension indicated by axis
        typedef std::pair< const uIntType*, const floatType* > valueMap;
        
        std::vector< valueMap > values;
        values.reserve( _points->size( ) / dim );

        for ( auto index = ownedIndices.begin( ); index != ownedIndices.end( ); index++ ){

            values.push_back( valueMap( &( *index), &( *_points )[ *index + _axis ] ) );

        }

        //Sort the value vector
        sort( values.begin( ), values.end( ),
              []( const valueMap &a,
                  const valueMap &b ){
                  return *a.second < *b.second;
              }
            );

        //Get the index of the current node
        _index = *values[ values.size( ) / 2 ].first;

        //Split the indices into left and right indices
        uIntVector left_indices;
        left_indices.reserve( values.size( ) / 2 );

        uIntVector right_indices;
        right_indices.reserve( values.size( ) - ( values.size( ) / 2 + 1 ) );

        for ( auto v = values.begin( ); v < values.begin( ) + values.size( ) / 2; v++ ){

            left_indices.push_back( *v->first );

        }

        for ( auto v = values.begin( ) + values.size( ) /2 + 1; v != values.end( ); v++ ){

            right_indices.push_back( *v->first );

        }

        //Form the left and right nodes
        if ( left_indices.size( ) > 0 ){

            left_child  = std::unique_ptr< KDNode >( new KDNode( _points, left_indices, depth + 1, dim ) );

        }

        if ( right_indices.size( ) > 0 ){

            right_child = std::unique_ptr< KDNode >( new KDNode( _points, right_indices, depth + 1, dim ) );

        }

        return;
    }

    const uIntType* KDNode::getIndex( ){
        /*!
         * Get the index associated with this KD node
         */

        return &_index;
    }

    void KDNode::getPointsInRange( const floatVector &upperBounds, const floatVector &lowerBounds,
                                   uIntVector &indices,
                                   floatVector *domainUpperBounds,
                                   floatVector *domainLowerBounds ){
        /*!
         * Get all of the points within the range specified by the rectangular boundary defined by the bounds
         *
         * :param const floatVector &upperBounds: The upper bounds of the boundary
         * :param const floatVector &lowerBounds: The lower bounds of the boundary
         * :param uIntVector &indices: The indices of the points ( or rather their first index in the point vector )
         *     of the points contained within the boundary
         * :param floatVector *domainUpperBounds: The upper boundary of the current domain. If NULL it will be
         *     calculated
         * :param floatVector &domainLowerBounds: The lower boundary of the current domain. If NULL it will be
         *     calculated
         */

        //Get the dimension
        uIntType dim = upperBounds.size( );

        floatVector _domainUpperBounds;
        floatVector _domainLowerBounds;

        //Assemble the current point
        floatVector median( _points->begin( ) + _index,
                            _points->begin( ) + _index + dim );

        floatVector upperDelta, lowerDelta;

        if( !domainUpperBounds ){
            
            _domainUpperBounds = floatVector( dim );
            _domainLowerBounds = floatVector( dim );

            for ( uIntType i = 0; i < dim; i++ ){

                _domainUpperBounds[ i ] = getMaximumValueDimension( i );
                _domainLowerBounds[ i ] = getMinimumValueDimension( i );

            }

            //Reassign the pointers
            domainUpperBounds = &_domainUpperBounds;
            domainLowerBounds = &_domainLowerBounds;

        }

        //Check if the current point is within the range
        upperDelta = upperBounds - median;
        lowerDelta = median - lowerBounds;

        if ( ( std::all_of( upperDelta.begin( ),
                            upperDelta.end( ),
                            [&]( floatType v ){ return v >= 0; } ) ) &&
               std::all_of( lowerDelta.begin( ),
                            lowerDelta.end( ),
                            [&]( floatType v ){ return v >= 0; } ) ){

            indices.push_back( _index );

        }

        //Check the left child exists and is within range for the current axis
        if ( ( left_child ) &&
             lowerDelta[ _axis ] >= 0 ){

             floatVector newDomainUpperBounds = *domainUpperBounds;
             newDomainUpperBounds[ _axis ] = median[ _axis ];

             left_child->getPointsInRange( upperBounds, lowerBounds,
                                           indices,
                                           &newDomainUpperBounds, domainLowerBounds );

        }

        //Check the right child exists and is within range
        if ( ( right_child ) &&
             upperDelta[ _axis ] >= 0 ){

             floatVector newDomainLowerBounds = *domainLowerBounds;
             newDomainLowerBounds[ _axis ] = median[ _axis ];

             right_child->getPointsInRange( upperBounds, lowerBounds,
                                            indices,
                                            domainUpperBounds, &newDomainLowerBounds );

        }

    }

    void KDNode::getPointsWithinRadiusOfOrigin( const floatVector &origin, const floatType &radius,
                                                uIntVector &indices,
                                                floatVector *domainUpperBounds,
                                                floatVector *domainLowerBounds ){
        /*!
         * Find all of the points within a given radius of the origin
         *
         * :param const floatVector &origin: The origin to use as the reference point for the origin
         * :param const floatType &radius: The radius to use as the bound
         * :param uIntVector &indices: The indices contained within the radius
         * :param floatVector *domainUpperBounds: The upper bounds of the domain associated with the node
         * :param floatVector *domainLowerBounds: The lower bounds of the domain associated with the node
         */

        //Get the dimension
        uIntType dim = origin.size( );

        floatVector _domainUpperBounds;
        floatVector _domainLowerBounds;

        //Assemble the current point
        floatVector median( _points->begin( ) + _index,
                            _points->begin( ) + _index + dim );

        if( !domainUpperBounds ){
            
            _domainUpperBounds = floatVector( dim );
            _domainLowerBounds = floatVector( dim );

            for ( uIntType i = 0; i < dim; i++ ){

                _domainUpperBounds[ i ] = getMaximumValueDimension( i );
                _domainLowerBounds[ i ] = getMinimumValueDimension( i );

            }

            //Reassign the pointers
            domainUpperBounds = &_domainUpperBounds;
            domainLowerBounds = &_domainLowerBounds;

        }

        //Check if the current point is within the range
        floatVector deltaVec = median - origin;
        floatType deltaRadiusSquared = vectorTools::dot( deltaVec, deltaVec );

        bool medianInside = deltaRadiusSquared <= ( radius * radius );

        if ( medianInside ){

            indices.push_back( _index );

        }

        //Check the left child exists and is within range for the current axis
        if ( ( left_child ) &&
             ( ( std::fabs( median[ _axis ] - origin[ _axis ] ) <= radius ) ||
               ( std::fabs( ( *domainLowerBounds )[ _axis ] - origin[ _axis ] ) <= radius ) ||
               ( ( median[ _axis ] >= origin[ _axis ] ) && ( origin[ _axis ] >= ( *domainLowerBounds )[ _axis ] ) )
             )
           ){

             floatVector newDomainUpperBounds = *domainUpperBounds;
             newDomainUpperBounds[ _axis ] = median[ _axis ];

             left_child->getPointsWithinRadiusOfOrigin( origin, radius, indices,
                                                        &newDomainUpperBounds, domainLowerBounds );

        }

        //Check the right child exists and is within range
        if ( ( right_child ) &&
             ( ( std::fabs( median[ _axis ] - origin[ _axis ] ) <= radius ) ||
               ( std::fabs( ( *domainUpperBounds )[ _axis ] - origin[ _axis ] ) <= radius ) ||
               ( ( ( *domainUpperBounds )[ _axis ] >= origin[ _axis ] ) && ( origin[ _axis ] >= median[ _axis ] ) ) 
             )
           ){

             floatVector newDomainLowerBounds = *domainLowerBounds;
             newDomainLowerBounds[ _axis ] = median[ _axis ];

             right_child->getPointsWithinRadiusOfOrigin( origin, radius, indices,
                                                         domainUpperBounds, &newDomainLowerBounds );

        }

    }

    floatType KDNode::getMinimumValueDimension( const uIntType &d ){
        /*!
         * Get the minimum value of a given dimension in the tree
         *
         * :param const uIntType &d: The dimension ( starting at zero ) of each point to search
         */

        floatType currentValue = ( *_points )[ _index + d ];

        if ( _axis == d ){

            if ( left_child ){

                return std::fmin( left_child->getMinimumValueDimension( d ),
                                  currentValue );

            }
            else{

                return currentValue;

            }

        }
        else{

            if ( ( left_child ) && ( !right_child ) ){

                return std::fmin( left_child->getMinimumValueDimension( d ),
                                  currentValue );

            }
            else if ( ( !left_child ) && ( right_child ) ){

                return std::fmin( right_child->getMinimumValueDimension( d ),
                                  currentValue );

            }
            else if ( ( left_child ) && ( right_child ) ){

                return std::fmin( std::fmin( left_child->getMinimumValueDimension( d ),
                                             right_child->getMinimumValueDimension( d ) ),
                                  currentValue );


            }
            else{

                return currentValue;

            }

        }

    }

    floatType KDNode::getMaximumValueDimension( const uIntType &d ){
        /*!
         * Get the maximum value of a given dimension in the tree
         *
         * :param const uIntType &d: The dimension ( starting at zero ) of each point to search
         */


        floatType currentValue = ( *_points )[ _index + d ];

        if ( _axis == d ){

            if ( right_child ){

                return std::fmax( right_child->getMaximumValueDimension( d ),
                                  currentValue );

            }
            else{

                return currentValue;

            }

        }
        else{

            if ( ( left_child ) && ( !right_child ) ){

                return std::fmax( left_child->getMaximumValueDimension( d ),
                                  currentValue );

            }
            else if ( ( !left_child ) && ( right_child ) ){

                return std::fmax( right_child->getMaximumValueDimension( d ),
                                  currentValue );

            }
            else if ( ( left_child ) && ( right_child ) ){

                return std::fmax( std::fmax( left_child->getMaximumValueDimension( d ),
                                             right_child->getMaximumValueDimension( d ) ),
                                  currentValue );

            }
            else{

                return currentValue;

            }

        }

    }

    void KDNode::printData( const uIntType &dim ){
        /*!
         * Print the data associated with this node to the terminal
         */

        std::cout << "NODE: " << _index << "\n";
        std::cout << "  depth: " << _depth << "\n";
        std::cout << "  value: "; vectorTools::print( floatVector( _points->begin( ) + _index,
                                                                   _points->begin( ) + _index + dim ) );
        std::cout << "  left: ";
        if ( left_child ){
            std::cout << *left_child->getIndex( ) << "\n";
        }
        else{
            std::cout << "NULL\n";
        }
        std::cout << "  right: ";
        if ( right_child ){
            std::cout << *right_child->getIndex( ) << "\n";
        }
        else{
            std::cout << "NULL\n";
        }
        std::cout << "\n";

        if ( left_child ){
            left_child->printData( dim );
        }

        if ( right_child ){
            right_child->printData( dim );
        }
    }


    volumeReconstructionBase::volumeReconstructionBase( ){
        /*!
         * The base volumeReconstruction constructor
         */

        _error = NULL;

        return;
    }

    volumeReconstructionBase::volumeReconstructionBase( const YAML::Node &config ){
        /*!
         * The base volumeReconstruction constructor
         *
         * type must be defined in the YAML file and will be the type of reconstruction to be used
         *
         * :param const YAML::Node &config: The YAML configuration file associated with the datafile
         */

        //Set the configuration file
        _config = config;

        if ( _config["type"] ){
            _type = _config["type"].as<std::string>();
        }
        else{
            _error = new errorNode( "volumeReconstructionBase", "The type must be defined in the YAML configuration file" );
            return;
        }

        _error = NULL;

        return;
    }

    volumeReconstructionBase::volumeReconstructionBase( const YAML::Node &config, errorOut error ) : volumeReconstructionBase( config ){
        /*!
         * The base volumeReconstruction constructor
         *
         * type must be defined in the YAML file and will be the type of reconstruction to be used
         *
         * :param const YAML::Node &config: The YAML configuration file associated with the datafile
         * :param errorOut error: The current value of the error
         */

        //Copy the error
        _error = error;

        return;
    }

    volumeReconstructionBase::~volumeReconstructionBase( ){
        /*!
         * The base volumeReconstruction destructor
         */

        if ( _config[ "write_config" ] ){
            if ( !_config[ "baseOutputFilename" ].IsScalar( ) ){
                _config[ "write_config" ] = "defaultOutput.yaml";
            }
            else{
                _config[ "write_config" ] = _config[ "baseOutputFilename" ].as< std::string >( ) + ".yaml";
            }

            std::ofstream yamlOut( _config[ "write_config" ].as< std::string >( ) + ".as_evaluated" );
            yamlOut << _config;
        }
    }

    std::shared_ptr<volumeReconstructionBase> volumeReconstructionBase::create( ){
        /*!
         * Create a new volumeReconstruction object using the configuration file
         */

        if ( _config[ "type" ] ){
            return volumeReconstructionBase::create( _config[ "type" ].as< std::string >( ) );
        }

        _error = new errorNode( "create", "The type is not defined" );
        return std::make_shared< volumeReconstructionBase >( _config, _error );
    }

    std::shared_ptr<volumeReconstructionBase> volumeReconstructionBase::create( const std::string &type ){
        /*!
         * Create a new volumeReconstruction object
         *
         * :param const std::string &type: The name of the object to be created.
         */

        auto T = registryMap.find( type );

        if ( T == registryMap.end() ){
            _error = new errorNode( "create", "The filetype ( " + type + " ) is not recognized" );
            return std::make_shared< volumeReconstructionBase >( _config, _error ); //The requested class is not defined
        }
        else{ //Register new volumeReconstruction objects below
            if ( T->second == DUAL_CONTOURING ){
                return std::make_shared< dualContouring >( _config );
            }
            else{
                _error = new errorNode( "create", "The volume reconstruction type ( " + type + " ) is not recognized" );
                return std::make_shared< volumeReconstructionBase >( _config, _error ); //The requested class is not defined
            }
        }

        _error = new errorNode( "create", "You should never get here..." );
        return std::make_shared< volumeReconstructionBase >( _config, _error ); //The requested class is not defined
    }

    errorOut volumeReconstructionBase::getError( ){
        /*!
         * Get the error
         */

        return _error;
    }

    errorOut volumeReconstructionBase::loadPoints( const floatVector *points ){
        /*!
         * Load the datapoints to the file
         *
         * :param const floatVector *points: A pointer to the data points. These points are stored as [ x1, y1, z1, x2, y2, z2, ... ]
         */

        if ( ( points->size( ) % _dim ) != 0 ){
            _error = new errorNode( "loadPoints", "The points vector's size is not consistent with the dimension" );
            return _error;
        }

        _points = points;
        _nPoints = _points->size( ) / _dim;

        uIntVector ownedIndices;
        ownedIndices.reserve( _nPoints );
        for ( uIntType i = 0; i < _dim * _nPoints; i+=_dim ){
            ownedIndices.push_back( i );
        }

        //Form the KD tree
        _pointTree = KDNode( _points, ownedIndices, 0, _dim );

        return NULL;
    }

    errorOut volumeReconstructionBase::computeMedianNeighborhoodDistance( ){
        /*!
         * Compute the median neighborhood distance i.e. the median distance
         * of the closest n points for each incoming point
         */

        floatVector distances( 0, 0 );

        for ( uIntType i = 0; i < _dim * _nPoints; i += _dim ){

            // Compute the n closest distances from the point
            floatVector x0( _points->begin( ) + i, _points->begin( ) + i + _dim );

            floatVector closestDistances( _nNeighborhoodPoints + 1, 0 );

            floatVector xi( _dim, 0 );

            for ( uIntType j = 0; j < _dim * ( _nNeighborhoodPoints + 1 ); j += _dim ){

                xi = floatVector( _points->begin( ) + j, _points->begin( ) + j + _dim );

                closestDistances[ j / _dim ] = vectorTools::l2norm( xi - x0 );

            }

            std::sort( closestDistances.begin( ), closestDistances.end( ) );

            for ( uIntType j = _dim * ( _nNeighborhoodPoints + 1 ); j < _dim * _nPoints; j += _dim ){

                xi = floatVector( _points->begin( ) + j, _points->begin( ) + j + _dim );

                floatType d = vectorTools::l2norm( xi - x0 );

                if ( d < closestDistances[ _nNeighborhoodPoints ] ){

                    closestDistances[ _nNeighborhoodPoints ] = d;
                    std::sort( closestDistances.begin( ), closestDistances.end( ) );

                }

            }

            distances.insert( distances.end(), closestDistances.begin() + 1, closestDistances.end( ) );

        }

        _medianNeighborhoodDistance = vectorTools::median( distances );

        return NULL;

    }

    errorOut volumeReconstructionBase::loadFunction( const floatVector *function ){
        /*!
         * Load the function values at the points to the file
         *
         * :param const floatVector *function: A pointer to the function values at the data points.
         *     These points are stored as [ x1, y1, z1, x2, y2, z2, ... ]
         */

        if ( function->size( ) != _nPoints ){
            _error = new errorNode( "loadPoints", "The function vector and the points vector are not consistent in size" );
            return _error;
        }

        _functionValues = function;

        return NULL;
    }

    errorOut volumeReconstructionBase::initialize( ){
        /*!
         * Base initialization
         */

        errorOut error = setInterpolationConfiguration( );

        if ( error ){
            errorOut result = new errorNode( "initialize", "Error in setting the interpolation configuration" );
            result->addNext( error );
            return result;
        }

        error = computeGeometryInformation( );

        if ( error ){
            errorOut result = new errorNode( "initialize", "Error in computation of the base geometry information" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut volumeReconstructionBase::evaluate( ){
        /*!
         * Base function which evaluates the volume reconstruction and prepares the domain
         * to use the other virtual functions
         *
         * A child class should call the base class' evaluate function to ensure that
         * the general information is computed.
         */

        errorOut error = initialize( );

        if ( error ){
            errorOut result = new errorNode( "evaluate", "Error in the base class initialize function" );
            result->addNext( error );
            return result;
        }

        setEvaluated( true );

        return NULL;
    }

    errorOut volumeReconstructionBase::performVolumeIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                                 floatVector &integratedValue ){
        /*!
         * Integrate a quantity known at the points over the volume returning the value for the domain.
         *
         * TODO: Consider replacing this error out with a default routine that assumes that the volume
         *       reconstruction generates a list of volumes associated with each point and uses those
         *       volumes to compute the overall volume. Leaving this as a virtual function enables the
         *       user to overload it if required.
         *
         * :param const floatVector &valuesAtPoints: A vector of the values at the data points. Stored as
         *     [ v_11, v_12, ..., v_21, v22, ... ] where the first index is the point index in order as 
         *     provided to the volume reconstruction object and the second index is the value of the 
         *     function to be integrated.
         * :param const uIntType valueSize: The size of the subvector associated with each of the datapoints.
         * :param floatVector &integratedValue: The final value of the integral
         */

        ( void ) valuesAtPoints;
        ( void ) valueSize;
        ( void ) integratedValue;

        return new errorNode( "performVolumeIntegration", "Volume integration not implemented" );
    }

    errorOut volumeReconstructionBase::performRelativePositionVolumeIntegration( const floatVector &valuesAtPoints,
                                                                                 const uIntType valueSize,
                                                                                 const floatVector &origin,
                                                                                 floatVector &integratedValue ){
        /*!
         * Integrate a quantity known at the points which has been dyaded with the relative position returning the
         * value for the domain.
         *
         * $V_ij = \int_{\mathcal{D}} v_i \left( x_j' - o_j \right) dv
         *
         * where $o_j$ is the origin and x_j' is the position of the point in the volume
         *
         * :param const floatVector &valuesAtPoints: A vector of the values at the data points. Stored as
         *     [ v_11, v_12, ..., v_21, v22, ... ] where the first index is the point index in order as 
         *     provided to the volume reconstruction object and the second index is the value of the 
         *     function to be integrated.
         * :param const uIntType valueSize: The size of the subvector associated with each of the datapoints.
         * :param const floatVector &origin: The coordinates of the origin.
         * :param floatVector &integratedValue: The final value of the integral
         */

        //Error handling
        if ( ( valuesAtPoints.size( ) / valueSize ) != ( getPoints( )->size( ) / _dim ) ){

            return new errorNode( "performRelativePositionVolumeIntegration",
                                  "The values at points vector is not consistent with the points vector in terms of size" );

        }

        //Compute the dyad at the points
        floatMatrix dyad;
        floatVector dyadVector;
        floatVector integrand;
        integrand.reserve( _dim * valuesAtPoints.size( ) );
        floatVector pointValue;
        floatVector pointPosition;

        uIntType index = 0;

        for ( auto point = getPoints( )->begin( ); point != getPoints( )->end( ); point+=_dim ){

            //Extract the function value at the point
            pointValue = floatVector( valuesAtPoints.begin( ) + ( index / _dim ) * valueSize,
                                      valuesAtPoints.begin( ) + ( index / _dim + 1 ) * valueSize );

            //Extract the point position
            pointPosition = floatVector( point, point + _dim );

            //Compute the dyadic product
            dyad = vectorTools::dyadic( pointValue, pointPosition - origin );

            //Flatten the dyadic product
            dyadVector = vectorTools::appendVectors( dyad );

            //Insert the values of the dyadic vector
            for ( auto v = dyadVector.begin( ); v != dyadVector.end( ); v++ ){

                integrand.push_back( *v );

            }

            index++;

        }

        //Perform the integration
        errorOut error = performVolumeIntegration( integrand, _dim * valueSize, integratedValue );

        if ( error ){

            errorOut result = new errorNode( "performRelativePositionVolumeIntegration",
                                             "Error in performing the volume integration" );
            result->addNext( error );
            return result;

        }

        return NULL;
    }

    errorOut volumeReconstructionBase::performSurfaceIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                                  floatVector &integratedValue, const uIntVector *subdomainIDs,
                                                                  const floatVector *subdomainWeights,
                                                                  const floatVector *macroNormal, const bool useMacroNormal ){
        /*!
         * Integrate a quantity known at the point over the surface and return the value for the domain.
         *
         * :param const floatVector &valuesAtPoints: A vector of the values at the data points. Stored as
         *     [ v_11, v_12, ..., v_21, v22, ... ] where the first index is the point index in order as 
         *     provided to the volume reconstruction object and the second index is the value of the 
         *     function to be integrated.
         * :param const uIntType valueSize: The size of the subvector associated with each of the datapoints.
         * :param floatVector &integratedValue: The final value of the integral
         * :param const uIntVector *subdomainIDs: The pointer to the subdomain of the surface to integrate over
         *     defaults to NULL and so the full domain is integrated over
         * :param const floatVector *subdomainWeights: The weights for the subdomains. Useful if points can be
         *     in multiple subdomains and they aren't small w.r.t. the domain size
         * :param const floatVector *macroNormal: A macro-scale normal vector to use to generate the micro
         *     weight. This can be helpful in cases where some points start to, ``wrap,'' around an edge 
         *     which should be flat. Can either be a single vector of dimension _dim or a collection of vectors
         *     at each boundary point.
         * :param const bool useMacroNormal: Use the macro-scale normal instead of the micro normals. Can help
         *     drive the integral to be more what is expected in some cases.
         */

        ( void ) valuesAtPoints;
        ( void ) valueSize;
        ( void ) integratedValue;
        ( void ) subdomainIDs;
        ( void ) subdomainWeights;
        ( void ) macroNormal;
        ( void ) useMacroNormal;

        return new errorNode( "performSurfaceIntegration", "Surface integration not implemented" );
    }

    errorOut volumeReconstructionBase::performPositionWeightedSurfaceIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                                                  floatVector &integratedValue, const uIntVector *subdomainIDs,
                                                                                  const floatVector *subdomainWeights,
                                                                                  const floatVector *macroNormal, const bool useMacroNormal ){
        /*!
         * Integrate a quantity known at the point over the surface times the position and return the value for the domain.
         *
         * :param const floatVector &valuesAtPoints: A vector of the values at the data points. Stored as
         *     [ v_11, v_12, ..., v_21, v22, ... ] where the first index is the point index in order as 
         *     provided to the volume reconstruction object and the second index is the value of the 
         *     function to be integrated.
         * :param const uIntType valueSize: The size of the subvector associated with each of the datapoints.
         * :param floatVector &integratedValue: The final value of the integral organized as
         *     [ v_11, v_12, v_13, ..., v_21, v_22, ... ] where the first index is the value and the second is the
         *     dimension
         * :param const uIntVector *subdomainIDs: The pointer to the subdomain of the surface to integrate over
         *     defaults to NULL and so the full domain is integrated over
         * :param const floatVector *subdomainWeights: The weights for the subdomains. Useful if points can be
         *     in multiple subdomains and they aren't small w.r.t. the domain size
         * :param const floatVector *macroNormal: A macro-scale normal vector to use to generate the micro
         *     weight. This can be helpful in cases where some points start to, ``wrap,'' around an edge 
         *     which should be flat. Can either be a single vector of dimension _dim or a collection of vectors
         *     at each boundary point.
         * :param const bool useMacroNormal: Use the macro-scale normal instead of the micro normals. Can help
         *     drive the integral to be more what is expected in some cases.
         */

        ( void ) valuesAtPoints;
        ( void ) valueSize;
        ( void ) integratedValue;
        ( void ) subdomainIDs;
        ( void ) subdomainWeights;
        ( void ) macroNormal;
        ( void ) useMacroNormal;

        return new errorNode( "performPositionWeightedSurfaceIntegration", "Surface integration not implemented" );
    }

    errorOut volumeReconstructionBase::performSurfaceFluxIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                                      floatVector &integratedValue, const uIntVector *subdomainIDs,
                                                                      const floatVector *subdomainWeights,
                                                                      const floatVector *macroNormal, const bool useMacroNormal ){
        /*!
         * Integrate the flux of a quantity known at the data points over the surface and return the value for the domain.
         *
         * \int_{\partial\mathcal{B}} n_i v_ij da \approx \sum_{p = 1}^N n_i^p v_ij^p da^p
         *
         * :param const floatVector &valuesAtPoints: A vector of the values at the data points. Stored as
         *     [ v_111, v_112, v_113, ..., v_121, v_122, v_123, ... ] where the first index is the point index in order as 
         *     provided to the volume reconstruction object, the second index is the first index of the v matirx and the
         *     third index is the second index of the value matrix
         *     function to be integrated.
         * :param const uIntType valueSize: The size of the subvector associated with each of the datapoints.
         * :param floatVector &integratedValue: The final value of the integral
         * :param const uIntVector *subdomainIDs: The pointer to the subdomain of the surface to integrate over
         *     defaults to NULL and so the full domain is integrated over
         * :param const floatVector *subdomainWeights: The weights for the subdomains. Useful if points can be
         *     in multiple subdomains and they aren't small w.r.t. the domain size
         * :param const floatVector *macroNormal: A macro-scale normal vector to use to generate the micro
         *     weight. This can be helpful in cases where some points start to, ``wrap,'' around an edge 
         *     which should be flat. Can either be a single vector of dimension _dim or a collection of vectors
         *     at each boundary point.
         * :param const bool useMacroNormal: Use the macro-scale normal instead of the micro normals. Can help
         *     drive the integral to be more what is expected in some cases.
         */

        ( void ) valuesAtPoints;
        ( void ) valueSize;
        ( void ) integratedValue;
        ( void ) subdomainIDs;
        ( void ) subdomainWeights;
        ( void ) macroNormal;
        ( void ) useMacroNormal;

        return new errorNode( "performSurfaceFluxIntegration", "Surface flux integration not implemented" );
    }

    errorOut volumeReconstructionBase::performRelativePositionSurfaceFluxIntegration( const floatVector &valuesAtPoints,
                                                                                      const uIntType valueSize,
                                                                                      const floatVector &origin,
                                                                                      floatVector &integratedValue,
                                                                                      const uIntVector *subdomainIDs,
                                                                                      const floatVector *subdomainWeights,
                                                                                      const floatVector *macroNormal,
                                                                                      const bool useMacroNormal ){
        /*!
         * Integrate the flux of a quantity known at the data points over the surface and return the value for the domain.
         *
         * \int_{\partial\mathcal{B}} n_i v_ij da \approx \sum_{p = 1}^N n_i^p v_ij^p da^p
         *
         * :param const floatVector &valuesAtPoints: A vector of the values at the data points. Stored as
         *     [ v_111, v_112, v_113, ..., v_121, v_122, v_123, ... ] where the first index is the point index in order as 
         *     provided to the volume reconstruction object, the second index is the first index of the v matirx and the
         *     third index is the second index of the value matrix
         *     function to be integrated.
         * :param const uIntType valueSize: The size of the subvector associated with each of the datapoints.
         * :param const floatVector &origin: The origin to compute the surface integral relative to 
         * :param floatVector &integratedValue: The final value of the integral
         * :param uIntVector *subdomainIDs: The IDs of points in the subdomain of the surface to integrate over
         * :param const floatVector *subdomainWeights: The weights for the subdomains. Useful if points can be
         *     in multiple subdomains and they aren't small w.r.t. the domain size
         * :param const floatVector *macroNormal: A macro-scale normal vector to use to generate the micro
         *     weight. This can be helpful in cases where some points start to, ``wrap,'' around an edge 
         *     which should be flat. Can either be a single vector of dimension _dim or a collection of vectors
         *     at each boundary point.
         * :param const bool useMacroNormal: Use the macro-scale normal instead of the micro normals. Can help
         *     drive the integral to be more what is expected in some cases.
         */

        ( void ) valuesAtPoints;
        ( void ) valueSize;
        ( void ) origin;
        ( void ) integratedValue;
        ( void ) subdomainIDs;
        ( void ) subdomainWeights;
        ( void ) macroNormal;
        ( void ) useMacroNormal;

        return new errorNode( "performSurfaceFluxIntegration", "Surface flux integration not implemented" );
    }

    errorOut volumeReconstructionBase::getSurfaceSubdomains( const floatType &minDistance, uIntVector &subdomainNodeCounts,
                                                             uIntVector &subdomainNodes ){
        /*!
         * Break the surface into subdomains which are separated by ( approximately ) minDistance.
         *
         * :param const floatType &minDistance: The minimum distance between the subdomain centers
         * :param uIntVector &subdomainNodeCounts: The number of nodes in each subdomain
         * :param uIntVector &subdomainNodes: The IDs of the nodes in the subdomains in the order of subdomainNodeCounts
         */

        ( void ) minDistance;
        ( void ) subdomainNodeCounts;
        ( void ) subdomainNodes;

        return new errorNode( "getSurfaceSubdomains", "Surface decomposition into subdomains not implemented" );
    }

    errorOut volumeReconstructionBase::writeToXDMF( ){
        /*!
         * Write the volume-reconstruction data to an XDMF file for review
         */

        return new errorNode( "writeToXDMF", "Not implemented" );
    }

    errorOut volumeReconstructionBase::setInterpolationConfiguration( ){
        /*!
         * Set the interpolation configuration
         */

        //Set the default interpolation configuration
        if ( !_config[ "interpolation" ] ){
            _config[ "interpolation" ][ "type" ] = "constant";
            _config[ "interpolation" ][ "constant_value" ] = 1;
            _config[ "interpolation" ][ "nNeighborhoodPoints" ] = 5;
            _nNeighborhoodPoints = 5;
            _functionValue = 1;
        }

        if ( !_config[ "interpolation" ][ "type" ] ){
            _config[ "interpolation" ][ "type" ] = "constant";
            _config[ "interpolation" ][ "constant_value" ] = 1;
            _config[ "interpolation" ][ "nNeighborhoodPoints" ] = 5;
            _nNeighborhoodPoints = 5;
            _functionValue = 1;
        }

        if ( ( _config[ "interpolation" ][ "type" ].as< std::string >( ).compare( "constant" ) == 0 ) &&
             ( !_config[ "interpolation" ][ "constant_value" ] ) ){

            _config[ "interpolation" ][ "constant_value" ] = 1;
            _config[ "interpolation" ][ "nNeighborhoodPoints" ] = 5;
            _nNeighborhoodPoints = 5;
            _functionValue = 1;
            
        }
        
        if ( ( _config[ "interpolation" ][ "type" ].as< std::string >( ).compare( "constant" ) == 0 ) &&
             ( _config[ "interpolation" ][ "constant_value" ] ) ){

            _functionValue = _config[ "interpolation" ][ "constant_value" ].as< floatType >( );
            _nNeighborhoodPoints = _config[ "interpolation" ][ "nNeighborhoodPoints" ].as< uIntType >( );
            
        }

        if ( _config[ "interpolation" ][ "type" ].as< std::string >( ).compare( "from_vector" ) == 0 ){
            
            //Check if _functionValues is null
            if ( !_functionValues ){
                return new errorNode(  "setInterpolationConfiguration",
                                       "'from_vector' is specified in the configuration but the function values have not been set\nThe use order is constructor -> loadPoints -> loadFunction -> evaluate" );
            }
        }

        return NULL;
    }

    errorOut volumeReconstructionBase::computeGeometryInformation( ){
        /*!
         * Compute basic geometry information of the domain
         */

        //Get the domain bounding box
        _upperBounds.resize( _dim );
        _lowerBounds.resize( _dim );

        for ( uIntType i = 0; i < _dim; i++ ){

            _upperBounds[ i ] = _pointTree.getMaximumValueDimension( i );
            _lowerBounds[ i ] = _pointTree.getMinimumValueDimension( i );

        }

        errorOut error = computeMedianNeighborhoodDistance( );

        if ( error ){

            errorOut result = new errorNode( __func__, "Error in computing the median neighborhood distance" );

            result->addNext( error );

            return result;

        }
        
        return NULL;
    }

    const floatVector *volumeReconstructionBase::getPoints( ){
        /*!
         * Get the points used in this class
         */

        return _points;
    }

    const floatVector *volumeReconstructionBase::getFunction( ){
        /*!
         * Get the values of the function at the points
         */

        return _functionValues;
    }

    errorOut volumeReconstructionBase::getFunctionValue( const uIntType i, floatType &value ){
        /*!
         * Get the function value at the provided index. This is the index of the function-value vector
         * corresponding to a given point not the index in the points vector.
         *
         * :param const uIntType i: The index of the point
         * :param floatType &value: The value of the function
         */

        if ( i > _nPoints ){
            
            return new errorNode( "getFunctionValue",
                                  "The index " + std::to_string( i ) + " is outside of the number of points" );

        }

        if ( !_functionValues ){

            value = _functionValue;

        }
        else{

            value = ( *_functionValues )[ i ];

        }

        return NULL;

    }

    const floatVector *volumeReconstructionBase::getLowerBounds( ){
        /*!
         * Get the lower bounds of the domain
         */

        return &_lowerBounds;
    }

    const floatVector *volumeReconstructionBase::getUpperBounds( ){
        /*!
         * Get the upper bounds of the domain
         */

        return &_upperBounds;
    }

    const floatType *volumeReconstructionBase::getMedianNeighborhoodDistance( ){
        /*!
         * Get the median distance between points and their neighborhood
         */

        return &_medianNeighborhoodDistance;
    }

    bool volumeReconstructionBase::getEvaluated( ){
        /*!
         * Get whether the volume reconstruction has been evaluated
         */

        return _isEvaluated;
    }

    void volumeReconstructionBase::setEvaluated( const bool isEvaluated ){
        /*!
         * Set whether the volume construction has been evaluated
         * This should be used in child classes' error handling
         * since the parent class will set this boolean to true 
         * upon successful evaluation but an error may occur in
         * the child class' evaluate function
         */

        _isEvaluated = isEvaluated;

        return;
    }

    errorOut volumeReconstructionBase::addBoundingPlanes( const floatMatrix &boundingPoints, const floatMatrix &boundingNormals ){
        /*!
         * Add planes which bound the domain. These planes MUST NOT form a convex surface!
         * 
         * :param const floatMatrix &boundingPoints: The points on the surfaces of the bounding planes
         * :param const floatMatrix &boundingNormals: The normals which define the bounding planes
         */

        if ( boundingPoints.size( ) != boundingNormals.size( ) ){

            return new errorNode( __func__, "The bounding points and bounding normals have different sizes" );

        }

        _boundingPlanes.clear( );
        _boundingPlanes.reserve( boundingPoints.size( ) );

        for ( uIntType i = 0; i < boundingPoints.size( ); i++ ){

            if ( boundingPoints[ i ].size( ) != _dim ){

                std::string message = "The point on bounding plane " + std::to_string( i ) + " has a dimension of "
                                    + std::to_string( boundingPoints[ i ].size( ) ) + " which is not equal to the dimension ( "
                                    + std::to_string( _dim ) + ")";

                return new errorNode( __func__, message );

            }

            if ( boundingNormals[ i ].size( ) != _dim ){

                std::string message = "The normal on bounding plane " + std::to_string( i ) + " has a dimension of "
                                    + std::to_string( boundingNormals[ i ].size( ) ) + " which is not equal to the dimension ( "
                                    + std::to_string( _dim ) + ")";

                return new errorNode( __func__, message );

            }

            _boundingPlanes.push_back( std::pair< floatVector, floatVector >( boundingPoints[ i ], boundingNormals[ i ] / vectorTools::l2norm( boundingNormals[ i ] ) ) );

        }

        _boundingSurfaces = true;

        return NULL;

    }

    const uIntVector *volumeReconstructionBase::getBoundaryIDs( ){
        /*!
         * Get a constant reference to the collection of boundary point ids
         */

        return NULL;
    }

    const floatVector *volumeReconstructionBase::getBoundaryPoints( ){
        /*!
         * Get a constant reference to the collection of boundary points
         */

        return NULL;
    }

    /*=========================================================================
    |                             dualContouring                              |
    =========================================================================*/

    dualContouring::dualContouring( ) : volumeReconstructionBase( ){
        /*!
         * The dualContouring default constructor
         */
        
        return;
    }

    dualContouring::dualContouring( const YAML::Node &configuration ) : volumeReconstructionBase( configuration ){
        /*!
         * The dualContouring constructor with a YAML node
         *
         * :param const YAMLL::Node &configuration: The YAML configuration file for the dualContouring object
         */

        if ( _error ){
            return;
        }
    }

    dualContouring::~dualContouring( ){
        /*!
         * The destructor of the dual contouring model
         */

        if ( _writeOutput ){
            writeToXDMF( );
        }
    }

    errorOut dualContouring::initialize( ){
        /*!
         * Initialization for the dualContouring method
         */

        //Preserve the base initialization
        errorOut error = volumeReconstructionBase::initialize( );

        if ( error ){

            errorOut result = new errorNode( "initialize", "Error in base initialization" );
            result->addNext( error );
            return result;

        }

        error = processConfigurationFile( );

        if ( error ){

            errorOut result = new errorNode( "initialize", "Error in processing the configuraiton file" );
            result->addNext( error );
            return result;

        }

        error = setGridSpacing( );

        if ( error ){

            errorOut result = new errorNode( "initialize", "Error in setting the grid spacing" );
            result->addNext( error );
            return result;

        }

        error = projectImplicitFunctionToBackgroundGrid( );

        if ( error ){

            errorOut result = new errorNode( "initialize", "Error in the projection of the implicit function to the background grid" );
            result->addNext( error );
            return result;

        }

        error = initializeInternalAndBoundaryCells( );

        if ( error ){

            errorOut result = new errorNode( "initialize", "Error when initializing the interior and boundary cells of the background grid" );
            result->addNext( error );
            return result;

        }

        error = computeBoundaryPointNormalsAndAreas( );

        if ( error ){

            errorOut result = new errorNode( "initialize", "Error when computing the boundary point normals and areas" );
            result->addNext( error );
            return result;

        }

        return NULL;
    }

    errorOut dualContouring::processConfigurationFile( ){
        /*!
         * Process the configuration file reading options and setting defaults.
         */

        if ( !_config[ "interpolation" ][ "discretization_count" ] ){

            if ( !_config[ "interpolation" ][ "grid_factor" ] ){

                _config[ "interpolation" ][ "grid_factor" ] = 1;

            }

            floatType grid_factor;
            if ( !_config[ "interpolation" ][ "grid_factor" ].IsScalar( ) ){

                return new errorNode( __func__, "The interpolation's 'grid_factor' must be a scalar" );

            }

            grid_factor = _config[ "interpolation" ][ "grid_factor" ].as< floatType >( );

            if ( grid_factor < 0 ){

                return new errorNode( __func__, "interpolation's 'grid_factor' must be positive!" );

            }

            floatVector delta = *getUpperBounds( ) - *getLowerBounds( );

            floatVector _discretization_count = ( grid_factor * delta / *getMedianNeighborhoodDistance( ) );

            uIntVector discretization_count( _dim, 0 );

            for ( unsigned int i = 0; i < _dim; i++ ){

                discretization_count[ i ] = ( uIntType )( _discretization_count[ i ] );

            }

            _config[ "interpolation" ][ "discretization_count" ] = discretization_count;

        }

        if ( _config[ "interpolation" ][ "discretization_count" ].IsScalar( ) ){

            uIntType v = _config[ "interpolation" ][ "discretization_count" ].as< uIntType >( );
            _domainDiscretization = { v, v, v };

        }
        else if ( _config[ "interpolation" ][ "discretization_count" ].IsSequence( ) ){

            _domainDiscretization.resize( _config[ "interpolation" ][ "discretization_count" ].size( ) );

            if ( _domainDiscretization.size( ) != _dim ){

                return new errorNode( "processConfigFile",
                                      "The number of discretization indices ( " + std::to_string( _domainDiscretization.size( ) ) 
                                      + " ) is not the same as the dimension ( " + std::to_string( _dim ) + " )" );

            }

            uIntType i = 0;
            for ( auto it  = _config[ "interpolation" ][ "discretization_count" ].begin( );
                       it != _config[ "interpolation" ][ "discretization_count" ].end( );
                       it++ ){

                _domainDiscretization[ i ] = it->as< uIntType>( );

                i++;
            }

        }
        else{
            
            return new errorNode( "processConfigFile", "The type of 'discretization_count' must be a scalar or sequence" );

        }

        if ( _config[ "interpolation" ][ "exterior_relative_delta" ] ){

            if ( _config[ "interpolation" ][ "exterior_relative_delta" ].IsScalar( ) ){

                _exteriorRelativeDelta = _config[ "interpolation" ][ "exterior_relative_delta" ].as< floatType >( );

            }
            else{

                return new errorNode( "processConfigFile", "Exterior relative delta must be a floating point number" );

            }

        }
        else{
            
            _config[ "interpolation" ][ "exterior_relative_delta" ] = _exteriorRelativeDelta;

        }

        if ( _config[ "interpolation" ][ "isosurface_cutoff" ] ){

            if ( _config[ "interpolation" ][ "isosurface_cutoff" ].IsScalar( ) ){

                _isosurfaceCutoff = _config[ "interpolation" ][ "isosurface_cutoff" ].as< floatType >( );

            }
            else{

                return new errorNode( "processConfigFile", "'isosurface_cutoff' must be a floating point number" );

            }

        }
        else{

            _config[ "interpolation" ][ "isosurface_cutoff" ] = _isosurfaceCutoff;

        }

        if ( _config[ "interpolation" ][ "absolute_tolerance" ] ){

            if ( _config[ "interpolation" ][ "absolute_tolerance" ].IsScalar( ) ){

                _absoluteTolerance = _config[ "interpolation" ][ "absolute_tolerance" ].as< floatType >( );

            }
            else{

                return new errorNode( "processConfigFile", "'absolute_tolerance' must be a floating point number" );

            }

        }
        else{

            _config[ "interpolation" ][ "absolute_tolerance" ] = _absoluteTolerance;

        }

        if ( _config[ "write_xdmf_output" ] ){
            
            _writeOutput = true;

            if ( _config[ "baseOutputFilename" ].IsScalar( ) ){

                _config[ "write_xdmf_output" ] = _config[ "baseOutputFilename" ].as< std::string >( );
                _XDMFOutputFilename = _config[ "baseOutputFilename" ].as< std::string >( );

            }
            else{

                _config[ "write_xdmf_output" ] = _XDMFOutputFilename;

            }

        }

        return NULL;
    }

    errorOut dualContouring::evaluate( ){
        /*!
         * Evaluate the dual contouring volume reconstruction
         */

        //Preserve the base class evaluate
        errorOut error = volumeReconstructionBase::evaluate( );

        if ( error ){
            errorOut result = new errorNode( "evaluate", "Error in base class evaluate" );
            result->addNext( error );
            return result;
        }

        setEvaluated( true );
        return NULL;
    }

    errorOut dualContouring::setGridSpacing( ){
        /*!
         * Set the spacing for the background grid
         */

        //The change in spacing
        floatType delta;

        //Set the grid locations
        _gridLocations.resize( _dim ); //Resize the grid to the dimension

        //Get the bounds
        const floatVector *upperBounds = getUpperBounds( );
        const floatVector *lowerBounds = getLowerBounds( );

        for ( uIntType i = 0; i < _dim; i++ ){

            //With the domain being split into n discretizations it requires n + 1 nodes
            //We also add two additional nodes on the outside of the grid to provide
            //a hard boundary. These nodes are located very close to the outside
            _gridLocations[ i ].resize( _domainDiscretization[ i ] + 3 );

            //Set the displacement of the nodes
            delta = ( ( *upperBounds )[ i ] - ( *lowerBounds )[ i ] ) / ( ( floatType )_domainDiscretization[ i ] );

            for ( uIntType j = 0; j < _domainDiscretization[ i ] + 1; j++ ){

                _gridLocations[ i ][ j + 1 ] = ( *lowerBounds )[ i ] + ( ( floatType )j ) * delta;

            }

            _gridLocations[ i ][ 0 ] = ( *lowerBounds )[ i ]
                                     - ( _exteriorRelativeDelta * delta + _absoluteTolerance );
            _gridLocations[ i ][ _domainDiscretization[ i ] + 2 ] = ( *upperBounds )[ i ]
                                                                  + ( _exteriorRelativeDelta * delta + _absoluteTolerance );

        }

        return NULL;

    }

    errorOut dualContouring::projectImplicitFunctionToBackgroundGrid( ){
        /*!
         * Interpolate the implicit function to the background grid.
         *
         * The way this is done, is not by computing the nodal averages but by computing
         * the shape-function weighted value of the implicit function at the nodes. Nodal
         * averages will come into play for the volume and surface integrals.
         */

        if ( _dim != 3 ){ //This must be 3d

            return new errorNode( "projectImplicitFunctionToBackgroundGrid",
                                  "A dimension of 3 is required for this routine" );

        }

        //Resize the implicit function values
        _implicitFunctionValues = floatVector( _gridLocations[ 0 ].size( ) *
                                               _gridLocations[ 1 ].size( ) *
                                               _gridLocations[ 2 ].size( ), 0 );

        uIntVector gridPointCounts( _gridLocations[ 0 ].size( ) *
                                    _gridLocations[ 1 ].size( ) *
                                    _gridLocations[ 2 ].size( ) );

        _length_scale = *getMedianNeighborhoodDistance( ) / ( 2 * std::sqrt( -std::log( 0.5 ) ) );
        _critical_radius = std::sqrt( -std::log( 1e-3 ) ) * 2 * _length_scale;

        //Loop over the elements

        errorOut error;
        uIntVector elementIndices;
        floatVector elementNodalContributions;
        uIntVector globalNodeIds;
        uIntVector pointCounts;

        uIntType ngx = _gridLocations[ 0 ].size( );
        uIntType ngy = _gridLocations[ 1 ].size( );
        uIntType ngz = _gridLocations[ 2 ].size( );
        uIntType nodeID;
        uIntVector pointIndices;

        floatVector node_x( _dim, 0 );

        for ( uIntType i = 1; i < ngx - 1; i++ ){

            for ( uIntType j = 1; j < ngy - 1; j++ ){

                for ( uIntType k = 1; k < ngz - 1; k++ ){

                    // Get the node ID
                    nodeID = ngy * ngz * i + ngz * j + k;

                    node_x = { _gridLocations[ 0 ][ i ], _gridLocations[ 1 ][ j ], _gridLocations[ 2 ][ k ] };

                    // Find the points within the critical radius
                    pointIndices.clear( );
                    _pointTree.getPointsWithinRadiusOfOrigin( node_x, _critical_radius, pointIndices );

                    for ( auto pI = pointIndices.begin( ); pI != pointIndices.end( ); pI++ ){

                        floatType value;

                        floatVector xi( getPoints( )->begin( ) + *pI, getPoints( )->begin( ) + *pI + _dim );

                        error = rbf( node_x, xi, _length_scale, value );

                        if ( error ){

                            errorOut result = new errorNode( __func__, "Error in the computation of the radial basis function" );

                            result->addNext( error );

                            return result;

                        }

//                        std::cerr << value << " ";

                        _implicitFunctionValues[ nodeID ] += value;

//                        std::cout << "value: " << value << "\n";

                    }

//                    std::cerr << "if: " << _implicitFunctionValues[ nodeID ] << "\n";
//
//                    return new errorNode( __func__, "derp" );
// 
//                    //Get the element contribution to the nodal values of the implicit function
//                    elementIndices = { i, j, k };
//                    error = processBackgroundGridElementImplicitFunction( elementIndices, elementNodalContributions,
//                                                                          globalNodeIds, pointCounts );
//
//                    if ( error ){
//
//                        errorOut result = new errorNode( "projectImplicitFunctionToBackgroundGrid",
//                                                         "Error in processing the projection of the implicit function to the nodes for the element with the lower cornrer of indices i, j, k: "
//                                                       + std::to_string( i ) + ", "
//                                                       + std::to_string( j ) + ", "
//                                                       + std::to_string( k ) );
//
//                        result->addNext( error );
//                        return result;
//
//                    }
//
//                    for ( uIntType i = 0; i < globalNodeIds.size( ); i++ ){
//
//                        //Add those values to the grid
//                        _implicitFunctionValues[ globalNodeIds[ i ] ] += elementNodalContributions[ i ];
//    
//                        //Add the point counts
//                        gridPointCounts[ globalNodeIds[ i ] ] += pointCounts[ i ];
//
//                    }

                }

            }

        }

//        for ( uIntType i = 1; i < ngx - 1; i++ ){
//
//            for ( uIntType j = 1; j < ngy - 1; j++ ){
//
//                for ( uIntType k = 1; k < ngz - 1; k++ ){
//
//                    if ( gridPointCounts[ ngy * ngz * i + ngz * j + k ] > 0 ){
//
//                        _implicitFunctionValues[ ngy * ngz * i + ngz * j + k ] /=
//                            ( floatType )gridPointCounts[ ngy * ngz * i + ngz * j + k ];
//
//                    }
//
//                }
//
//            }
//
//        }

        _implicitFunctionValues -= _isosurfaceCutoff;

        return NULL;
    }

    errorOut dualContouring::getGridElement( const uIntVector &indices, std::unique_ptr< elib::Element > &element ){
        /*!
         * Get an element from the grid at the given lower left hand corner nodal indices
         *
         * :param const uIntVector &indices: The indices of the lower-left hand corner of the element
         * :param std::unique_ptr< elib::Element > &element: The pointer to the finite element representation
         */

        if ( _dim != 3 ){
            return new errorNode( "getGridElement",
                                  "A dimension of 3 is required for this routine" );
        }

        if ( indices.size( ) != _dim ){
            return new errorNode( "getGridElement",
                                  "The indices must have the same number of values as the dimension" );
        }

        uIntType ngy = _gridLocations[ 1 ].size( );
        uIntType ngz = _gridLocations[ 2 ].size( );

        //Compute the location of the nodes
        floatVector lbCoordinates( _dim, 0 );
        floatVector ubCoordinates( _dim, 0 );
        for ( auto index = indices.begin( ); index != indices.end( ); index++ ){

            if ( _gridLocations[ index - indices.begin( ) ].size( ) <= *index + 1 ){

                return new errorNode( "getGridElement",
                                      "An index of " + std::to_string( *index )
                                      + " and / or that index plus one is outside the bounds of the defined grid locations" );
                    
            }
            
            lbCoordinates[ index - indices.begin( ) ] = _gridLocations[ index - indices.begin( ) ][ *index ];
            ubCoordinates[ index - indices.begin( ) ] = _gridLocations[ index - indices.begin( ) ][ *index + 1 ];

        }

        //Determine the element's nodes
        floatMatrix nodes( 8, floatVector( _dim, 0 ) );
        nodes = { { lbCoordinates[ 0 ], lbCoordinates[ 1 ], lbCoordinates[ 2 ] },
                  { ubCoordinates[ 0 ], lbCoordinates[ 1 ], lbCoordinates[ 2 ] },
                  { ubCoordinates[ 0 ], ubCoordinates[ 1 ], lbCoordinates[ 2 ] },
                  { lbCoordinates[ 0 ], ubCoordinates[ 1 ], lbCoordinates[ 2 ] },
                  { lbCoordinates[ 0 ], lbCoordinates[ 1 ], ubCoordinates[ 2 ] },
                  { ubCoordinates[ 0 ], lbCoordinates[ 1 ], ubCoordinates[ 2 ] },
                  { ubCoordinates[ 0 ], ubCoordinates[ 1 ], ubCoordinates[ 2 ] },
                  { lbCoordinates[ 0 ], ubCoordinates[ 1 ], ubCoordinates[ 2 ] } };

        //Set the element's global node numbers
        uIntVector global_node_ids( 8, 0 );
        global_node_ids=
            {
                ngy * ngz * ( indices[ 0 ] + 0 ) + ngz * ( indices[ 1 ] + 0 ) + ( indices[ 2 ] + 0 ),
                ngy * ngz * ( indices[ 0 ] + 1 ) + ngz * ( indices[ 1 ] + 0 ) + ( indices[ 2 ] + 0 ),
                ngy * ngz * ( indices[ 0 ] + 1 ) + ngz * ( indices[ 1 ] + 1 ) + ( indices[ 2 ] + 0 ),
                ngy * ngz * ( indices[ 0 ] + 0 ) + ngz * ( indices[ 1 ] + 1 ) + ( indices[ 2 ] + 0 ),
                ngy * ngz * ( indices[ 0 ] + 0 ) + ngz * ( indices[ 1 ] + 0 ) + ( indices[ 2 ] + 1 ),
                ngy * ngz * ( indices[ 0 ] + 1 ) + ngz * ( indices[ 1 ] + 0 ) + ( indices[ 2 ] + 1 ),
                ngy * ngz * ( indices[ 0 ] + 1 ) + ngz * ( indices[ 1 ] + 1 ) + ( indices[ 2 ] + 1 ),
                ngy * ngz * ( indices[ 0 ] + 0 ) + ngz * ( indices[ 1 ] + 1 ) + ( indices[ 2 ] + 1 )
            };

        //Get the element
        auto qrule = elib::default_qrules.find( _elementType );
        if ( qrule == elib::default_qrules.end( ) ){

            return new errorNode( "getGridElement",
                                  "The default quadruature rule for the background grid element ( " + _elementType
                                 +" ) was not found" );

        }

        element = elib::build_element_from_string( _elementType, global_node_ids, nodes, qrule->second );

        return NULL;

    }

    errorOut dualContouring::rbf( const floatVector &x, const floatVector &x0, const floatType &ls, floatType &val ){
        /*!
         * Compute the radial basis function around the provided point. The form of the function is:
         * 
         * v = exp(-(r / (2 * ls) )**2 )
         * 
         * where
         * 
         * r = ||x - x0||
         * 
         * :param const floatVector &x: The point at which to compute the radial basis function
         * :param const floatVector &x0: The reference point for the radial basis function
         * :param const floatType &ls: The lengthscale
         * :param floatType &val: The value of the radial basis function
         */

        if ( x.size( ) != x0.size( ) ){

            return new errorNode( __func__, "The size of x (" + std::to_string( x.size( ) ) + ") and x0 ( " + std::to_string( x0.size( ) ) + ") are not the same" );

        }

        floatType r = vectorTools::l2norm( x - x0 );
        val = std::exp( -std::pow( r / ( 2 * ls ), 2 ) );

        if ( _boundingSurfaces ){

            for ( auto plane = _boundingPlanes.begin( ); plane != _boundingPlanes.end( ); plane++ ){

                floatType d = vectorTools::dot( plane->second, x - plane->first );

                if ( d >= 0 ){

                    val = 0;

                    return NULL;

                }

            }

        }

        return NULL;

    }

    errorOut dualContouring::processBackgroundGridElementImplicitFunction( const uIntVector &indices,
                                                                           floatVector &implicitFunctionNodalValues,
                                                                           uIntVector  &globalNodeIDs,
                                                                           uIntVector  &pointCounts ){
        /*!
         * Project the values of the implicit function of the points contained within an element of
         * the background grid to the background grid's nodes. This projection is done by adding the
         * contribution of the points to the nearest grid node. The resulting implicit function at the 
         * nodes is the average of all of the points contributing to a individual grid
         *
         * :param const uIntVector &indices: The indices of the lower-left corner of the hex element.
         *     The nodes of the element are the indices provided plus the nodes defined at  i + 1, 
         *     j + 1, k + 1 and the other combinations.
         * :param floatVector &implicitFunctionNodalValues: The nodal values of the implicit function.
         * :param uIntVector &globalNodeIds: The global node ids corresponding to the nodal values
         * :param uIntVector &pointCounts: The number of points contributing to each node.
         */

        std::unique_ptr< elib::Element > element;
        errorOut error = getGridElement( indices, element );

        if ( error ){

            errorOut result = new errorNode( "processBackgroundGridElementImplicitFunction",
                                             "Error in getting the element of the current grid indices" );
            result->addNext( error );
            return result;

        }

        //Determine the points which are contained within this element
        uIntVector pointIndices;
        floatVector domainUpperBounds = *getUpperBounds( );
        floatVector domainLowerBounds = *getLowerBounds( );

        _pointTree.getPointsInRange( element->bounding_box[ 1 ], element->bounding_box[ 0 ], pointIndices,
                                     &domainUpperBounds, &domainLowerBounds );

        //If there are no points contained within this element, return
        if ( indices.size( ) == 0 ){

            return NULL;

        }

        //Compute the local coordinates of the nodes and their shapefunctions and project the 
        floatVector globalCoordinates, localCoordinates, nodesSupported( element->nodes.size( ), 0 );
        pointCounts = uIntVector( element->nodes.size( ), 0 );

        //Initialize the implicit function's values at the nodes
        implicitFunctionNodalValues = floatVector( element->nodes.size( ), 0 );
        floatType fxn;

        floatVector distances( element->nodes.size( ) );
        floatVector p;

        uIntType ngy = _gridLocations[ 1 ].size( );
        uIntType ngz = _gridLocations[ 2 ].size( );
        floatType minDistance;

        for ( auto pI = pointIndices.begin( ); pI != pointIndices.end( ); pI++ ){

            //Get the current point's location
            p = floatVector( getPoints( )->begin( ) + *pI, getPoints( )->begin( ) + *pI + _dim );

            //Compute the distances of the point to the nodes
            for ( auto node = element->nodes.begin( ); node != element->nodes.end( ); node++ ){

                distances[ node - element->nodes.begin( ) ] = vectorTools::l2norm( p - *node );

            }

            //Get the minimum distance
            minDistance = *std::min_element( distances.begin( ), distances.end( ) );

            //Determine which values are equal to the smallest value
            uIntType dIndex = 0;
            for ( auto d = distances.begin( ); d != distances.end( ); d++ ){

                nodesSupported[ dIndex ] = ( floatType )( vectorTools::fuzzyEquals( *d, minDistance ) );
                pointCounts[ dIndex ] += ( uIntType )nodesSupported[ dIndex ];

                dIndex++;

            }

            //Get the implicit function value
            errorOut error = getFunctionValue( *pI / _dim, fxn );

            if ( error ){

                errorOut result = new errorNode( "processBackgroundGridElementImplicitFunction",
                                                 "Error in getting the function value" );
                result->addNext( error );
                return result;

            }

            //Project the implicit function to the nodes
            implicitFunctionNodalValues += fxn * nodesSupported;

        }

        //Set the global node ids
        globalNodeIDs = element->global_node_ids;

        return NULL;

    }

    errorOut dualContouring::initializeInternalAndBoundaryCells( ){
        /*!
         * Initialize the cells of the background grid which are internal and on the boundary
         */

        errorOut error = findInternalAndBoundaryCells( );

        if ( error ){

            errorOut result = new errorNode( "initializeInternalAndBoundaryCells",
                                             "Error when finding the internal and boundary cells" );
            result->addNext( error );
            return result;

        }

        error = computeBoundaryPoints( );

        if ( error ){

            errorOut result = new errorNode( "computeBoundaryPoints",
                                             "Error in the computation of the boundary points" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut dualContouring::findInternalAndBoundaryCells( ){
        /*!
         * Find the cells which are internal and on the boundary
         */

        if ( _dim != 3 ){
            
            return new errorNode( "findInternalAndBoundaryCells", "This function requires that the dimension is 3D" );

        }

        uIntType ngx = _gridLocations[ 0 ].size( );
        uIntType ngy = _gridLocations[ 1 ].size( );
        uIntType ngz = _gridLocations[ 2 ].size( );

        //Resize the internal cells vector
        _internalCells.clear( );
        _boundaryCells.clear( );

        _internalCells.reserve( ( ngx - 1 ) * ( ngy - 1 ) * ( ngz - 1 ) );
        _boundaryCells.reserve( ( ngx - 1 ) * ( ngy - 1 ) * ( ngz - 1 ) );

        floatVector cellValues;

        for ( uIntType i = 0; i < ( ngx - 1 ); i++ ){

            for ( uIntType j = 0; j < ( ngy - 1 ); j++ ){

                for ( uIntType k = 0; k < ( ngz - 1 ); k++ ){

                    //Get the values of the implicit function
                    cellValues = { _implicitFunctionValues[ ngy * ngz * ( i + 0 ) + ngz * ( j  + 0 ) + ( k + 0 ) ],
                                   _implicitFunctionValues[ ngy * ngz * ( i + 0 ) + ngz * ( j  + 0 ) + ( k + 1 ) ],
                                   _implicitFunctionValues[ ngy * ngz * ( i + 0 ) + ngz * ( j  + 1 ) + ( k + 0 ) ],
                                   _implicitFunctionValues[ ngy * ngz * ( i + 0 ) + ngz * ( j  + 1 ) + ( k + 1 ) ],
                                   _implicitFunctionValues[ ngy * ngz * ( i + 1 ) + ngz * ( j  + 0 ) + ( k + 0 ) ],
                                   _implicitFunctionValues[ ngy * ngz * ( i + 1 ) + ngz * ( j  + 0 ) + ( k + 1 ) ],
                                   _implicitFunctionValues[ ngy * ngz * ( i + 1 ) + ngz * ( j  + 1 ) + ( k + 0 ) ],
                                   _implicitFunctionValues[ ngy * ngz * ( i + 1 ) + ngz * ( j  + 1 ) + ( k + 1 ) ] };

                    if ( std::any_of( cellValues.begin( ), cellValues.end( ),
                                      []( floatType v ){ return v > 0; } ) ){

                        //The cell contributes to the overall volume of the domain
                        _internalCells.push_back( ngy * ngz * i + ngz * j + k );

                        if ( std::any_of( cellValues.begin( ), cellValues.end( ),
                                          []( floatType v ){ return v <= 0; } ) ){

                            //The cell is on a surface of the body
                            _boundaryCells.push_back( ngy * ngz * i + ngz * j + k );

                        }

                    }

                }

            }

        }

        return NULL;
    }

    errorOut dualContouring::computeBoundaryPoints( ){
        /*!
         * Compute the points which define the boundary
         */

        if ( _dim != 3 ){

            return new errorNode( "computeBoundaryPoints", "This function requires that the dimension is 3D" );

        }

        uIntType ngy = _gridLocations[ 1 ].size( );
        uIntType ngz = _gridLocations[ 2 ].size( );

        //Resize the boundary point vector
        _boundaryPoints.clear( );

        _boundaryPoints.reserve( _dim * _boundaryCells.size( ) );
        _boundaryPointIDToIndex.reserve( _boundaryCells.size( ) );

        //Loop over the boundary cells
        uIntType i, j, k;
        uIntType ri1, rj1, rk1;
        uIntType ri2, rj2, rk2;
        uIntType ri, rj, rk;
        floatVector cellValues;
        errorOut error;
        std::unique_ptr< elib::Element > element;
        std::vector< bool > edgeTransition;

        uIntVector globalNodeIds;

        uIntVector edgeNodes =
            {
                0, 1, 3, 2, 4, 5, 7, 6, //x local node numbers
                1, 2, 0, 3, 5, 6, 4, 7, //y local node numbers
                0, 4, 1, 5, 2, 6, 3, 7  //z local node numbers
            };

        //Reserve the number of potential intersected edges
        _boundaryEdges_x.clear( );
        _boundaryEdges_x.reserve( 8 * _boundaryCells.size( ) ); //This is a worst case scenario

        _boundaryEdges_y.clear( );        
        _boundaryEdges_y.reserve( 8 * _boundaryCells.size( ) ); //This is a worst case scenario

        _boundaryEdges_z.clear( );        
        _boundaryEdges_z.reserve( 8 * _boundaryCells.size( ) ); //This is a worst case scenario

        floatMatrix points, normals;
        floatType m, b;

        floatVector rootNode;

        uIntVector supportingPoints;
        floatVector intersectionPoint( _dim, 0 );
        floatVector localIntersectionPoint;
        floatVector gradient;

        intMatrix intArgs;
        intMatrix intOuts;

        floatMatrix floatArgs;
        floatMatrix floatOuts;

        //Initialize the residual equation
        solverTools::stdFncNLFJ func;
        func = static_cast<solverTools::NonLinearFunctionWithJacobian>(dualContouringInternalPointResidual);

        floatVector X0, X;

        floatVector localBoundaryPoint;
        floatVector boundaryPoint;

        uIntType edgeID;
        uIntVector edgeCells;

        bool flipDirection;

        uIntVector ownedIndices( _boundaryCells.size( ) );

        for ( auto bc = _boundaryCells.begin( ); bc != _boundaryCells.end( ); bc++ ){

            //Determine the lower-left hand corner index
            i = *bc / ( ngy * ngz );
            j = ( *bc - ( ngy * ngz * i ) ) / ngz;
            k = *bc - ngy * ngz * i - ngz * j;

            //Determine the grid element
            element = NULL;
            error = getGridElement( { i, j, k }, element );

            if ( error ){

                errorOut result = new errorNode( "computeBoundaryPoints",
                                                 "Error in construction of the grid element" );
                result->addNext( error );
                return result;

            }

            //Extract the implicit function values
            cellValues =
                {
                    _implicitFunctionValues[ element->global_node_ids[ 0 ] ],
                    _implicitFunctionValues[ element->global_node_ids[ 1 ] ],
                    _implicitFunctionValues[ element->global_node_ids[ 2 ] ],
                    _implicitFunctionValues[ element->global_node_ids[ 3 ] ],
                    _implicitFunctionValues[ element->global_node_ids[ 4 ] ],
                    _implicitFunctionValues[ element->global_node_ids[ 5 ] ],
                    _implicitFunctionValues[ element->global_node_ids[ 6 ] ],
                    _implicitFunctionValues[ element->global_node_ids[ 7 ] ]
                };

            //Determine the edges where transitions occur
            edgeTransition =
                {
                    //Edges oriented along the x axis
                    std::signbit( cellValues[ 0 ] ) != std::signbit( cellValues[ 1 ] ),
                    std::signbit( cellValues[ 2 ] ) != std::signbit( cellValues[ 3 ] ),
                    std::signbit( cellValues[ 4 ] ) != std::signbit( cellValues[ 5 ] ),
                    std::signbit( cellValues[ 6 ] ) != std::signbit( cellValues[ 7 ] ),
                    //Edges oriented along the y axis
                    std::signbit( cellValues[ 1 ] ) != std::signbit( cellValues[ 2 ] ),
                    std::signbit( cellValues[ 3 ] ) != std::signbit( cellValues[ 0 ] ),
                    std::signbit( cellValues[ 5 ] ) != std::signbit( cellValues[ 6 ] ),
                    std::signbit( cellValues[ 7 ] ) != std::signbit( cellValues[ 4 ] ),
                    //Edges oriented along the z axis
                    std::signbit( cellValues[ 0 ] ) != std::signbit( cellValues[ 4 ] ),
                    std::signbit( cellValues[ 1 ] ) != std::signbit( cellValues[ 5 ] ),
                    std::signbit( cellValues[ 2 ] ) != std::signbit( cellValues[ 6 ] ),
                    std::signbit( cellValues[ 3 ] ) != std::signbit( cellValues[ 7 ] )
                };

            points.clear( );
            points.reserve( _dim * edgeTransition.size( ) );
            normals.clear( );
            normals.reserve( _dim * edgeTransition.size( ) );

            for ( auto eT = edgeTransition.begin( ); eT != edgeTransition.end( ); eT++ ){

                //Check if this edge was intersected
                if ( !( *eT ) ){
                    continue;
                }

                //Get the intersection points of the transition edges

                uIntType i2 = edgeNodes[ 2 * ( eT - edgeTransition.begin( ) ) + 1 ];
                uIntType i1 = edgeNodes[ 2 * ( eT - edgeTransition.begin( ) ) + 0 ];

                for ( uIntType i = 0; i < _dim; i++ ){

                    if ( ( element->reference_nodes[ i2 ][ i ] - element->reference_nodes[ i1 ][ i ] ) < _absoluteTolerance ){

                        intersectionPoint[ i ] = element->reference_nodes[ i2 ][ i ];

                    }
                    else{

                        m = ( cellValues[ i2 ] - cellValues[ i1 ] ) / ( element->reference_nodes[ i2 ][ i ] - element->reference_nodes[ i1 ][ i ] );
                        b = cellValues[ i1 ] - m * element->reference_nodes[ i1 ][ i ];

                        intersectionPoint[ i ] = -b / m;

                    }

                }

                //Compute the local coordinates of the intersection point
                error = element->compute_local_coordinates( intersectionPoint, localIntersectionPoint );

                if ( error ){

                    errorOut result = new errorNode( "computeBoundaryPoints", "Error in computation of the local coordinates of the intersection point" );
                    result->addNext( error );
                    return result;

                }

                //Put the local coordinates of the point into the points vector
                points.push_back( localIntersectionPoint );

                //Get the global node ids of the nodes on the edges
                ri1 = element->global_node_ids[ i1 ] / ( ngy * ngz );
                rj1 = ( element->global_node_ids[ i1 ] - ( ngy * ngz * ri1 ) ) / ngz;
                rk1 = element->global_node_ids[ i1 ] - ngy * ngz * ri1 - ngz * rj1;

                ri2 = element->global_node_ids[ i2 ] / ( ngy * ngz );
                rj2 = ( element->global_node_ids[ i2 ] - ( ngy * ngz * ri2 ) ) / ngz;
                rk2 = element->global_node_ids[ i2 ] - ngy * ngz * ri2 - ngz * rj2;

                //Approximate the normal at the transition point
                if ( cellValues[ i2 ] > cellValues[ i1 ] ){

                    rootNode = element->reference_nodes[ i2 ];
                    ri = ri2;
                    rj = rj2;
                    rk = rk2;

                    flipDirection = false;

                }
                else{

                    rootNode = element->reference_nodes[ i1 ];
                    ri = ri1;
                    rj = rj1;
                    rk = rk1;

                    flipDirection = true;

                }

                //Find the points which contribute to the root node
                floatVector upperBounds =
                    {
                        0.5 * ( _gridLocations[ 0 ][ ri + 1 ] - rootNode[ 0 ] ) + rootNode[ 0 ],
                        0.5 * ( _gridLocations[ 1 ][ rj + 1 ] - rootNode[ 1 ] ) + rootNode[ 1 ],
                        0.5 * ( _gridLocations[ 2 ][ rk + 1 ] - rootNode[ 2 ] ) + rootNode[ 2 ]
                    };

                floatVector lowerBounds =
                    {
                        0.5 * ( _gridLocations[ 0 ][ ri - 1 ] - rootNode[ 0 ] ) + rootNode[ 0 ],
                        0.5 * ( _gridLocations[ 1 ][ rj - 1 ] - rootNode[ 1 ] ) + rootNode[ 1 ],
                        0.5 * ( _gridLocations[ 2 ][ rk - 1 ] - rootNode[ 2 ] ) + rootNode[ 2 ]
                    };

                //Find the points that contribute to this root node
                floatVector domainUpperBounds = *getUpperBounds( );
                floatVector domainLowerBounds = *getLowerBounds( );
                supportingPoints.clear( );
                _pointTree.getPointsInRange( upperBounds, lowerBounds, supportingPoints, &domainUpperBounds, &domainLowerBounds );

                //Determine the normal at the intersection point
                if ( ( supportingPoints.size( ) >= _minNormalApproximationCount ) && ( _useMaterialPointsForNormals ) ){

                    return new errorNode( "computeBoundaryPoints", "Using the material points for normals has not been implemented yet" );

                }
                else{

                    //Use the gradient of the shape function to compute the normal
                   
                    //Compute the global gradient
                    error = element->get_global_gradient( cellValues, localIntersectionPoint, gradient );

                    //Compute the normal vector
                    gradient /= vectorTools::l2norm( gradient );

                    normals.push_back( gradient );

                }

                //Store the transition edge
                uIntType edgeIndex = eT - edgeTransition.begin( );
                edgeID = ngy * ngz * ri1 + ngz * rj1 + rk1;

                if ( edgeIndex < 4 ){ //x edge

                    edgeCells =
                        {
                            ngy * ngz * ri1 + ngz * ( rj1 - 0 ) + ( rk1 - 1 ),
                            ngy * ngz * ri1 + ngz * ( rj1 - 1 ) + ( rk1 - 1 ),
                            ngy * ngz * ri1 + ngz * ( rj1 - 1 ) + ( rk1 - 0 ),
                            ngy * ngz * ri1 + ngz * ( rj1 - 0 ) + ( rk1 - 0 )
                        };

                    //Check the direction of the normal and determine if the ordering
                    //needs to be flipped
                    if ( flipDirection ){
                        edgeCells = { edgeCells[ 3 ], edgeCells[ 2 ], edgeCells[ 1 ], edgeCells[ 0 ] };
                    }

                    if ( _boundaryEdges_x.find( edgeID ) == _boundaryEdges_x.end( ) ){

                            _boundaryEdges_x.emplace( edgeID, edgeCells );

                    }

                }
                else if ( edgeIndex < 8 ){ //y edge

                    edgeCells =
                        {
                            ngy * ngz * ( ri1 - 0 ) + ngz * rj1 + ( rk1 - 0 ),
                            ngy * ngz * ( ri1 - 1 ) + ngz * rj1 + ( rk1 - 0 ),
                            ngy * ngz * ( ri1 - 1 ) + ngz * rj1 + ( rk1 - 1 ),
                            ngy * ngz * ( ri1 - 0 ) + ngz * rj1 + ( rk1 - 1 )
                        };

                    //Check the direction of the normal and determine if the ordering
                    //needs to be flipped
                    if ( flipDirection ){
                        edgeCells = { edgeCells[ 3 ], edgeCells[ 2 ], edgeCells[ 1 ], edgeCells[ 0 ] };
                    }

                    if ( _boundaryEdges_y.find( edgeID ) == _boundaryEdges_y.end( ) ){

                            _boundaryEdges_y.emplace( edgeID, edgeCells );

                    }

                }
                else{ //z edge

                    edgeCells =
                        {
                            ngy * ngz * ( ri1 - 0 ) + ngz * ( rj1 - 1 ) + rk1,
                            ngy * ngz * ( ri1 - 1 ) + ngz * ( rj1 - 1 ) + rk1,
                            ngy * ngz * ( ri1 - 1 ) + ngz * ( rj1 - 0 ) + rk1,
                            ngy * ngz * ( ri1 - 0 ) + ngz * ( rj1 - 0 ) + rk1
                        };
                    //Check the direction of the normal and determine if the ordering
                    //needs to be flipped
                    if ( flipDirection ){
                        edgeCells = { edgeCells[ 3 ], edgeCells[ 2 ], edgeCells[ 1 ], edgeCells[ 0 ] };
                    }

                    if ( _boundaryEdges_z.find( edgeID ) == _boundaryEdges_z.end( ) ){

                            _boundaryEdges_z.emplace( edgeID, edgeCells );

                    }

                }

            }

            //Solve for the boundary point
            bool converged;
            bool fatalError;

            //Assemble the floating point argments matrix
            floatArgs.resize( 2 + 2 * points.size( ), floatVector( _dim, 0 ) );
            floatArgs[ 0 ] = {  1,  1,  1 };//element->bounding_box[ 1 ];
            floatArgs[ 1 ] = { -1, -1, -1 };//element->bounding_box[ 0 ];

            for ( uIntType i = 0; i < points.size( ); i++ ){

                floatArgs[ 2 + i ] = points[ i ];
                floatArgs[ 2 + points.size( ) + i ] = normals[ i ];

            }

            //Assemble the integer arguments matrix
            intArgs = { { ( int )_dim, ( int )points.size( ) } }; 

            //Assemble the initial iterate
            X0.resize( 5 * _dim );
            X.resize( 5 * _dim );
            for ( uIntType i = 0; i < _dim; i++ ){

//                X0[            i ] = 0.5 * ( element->bounding_box[ 0 ][ i ] + element->bounding_box[ 1 ][ i ] ); //x
                X0[            i ] = 0.0;
                X0[     _dim + i ] = floatArgs[ 0 ][ i ] - X0[ i ]; //s
                X0[ 2 * _dim + i ] = X0[ i ] - floatArgs[ 1 ][ i ]; //t
                X0[ 3 * _dim + i ] = 0.5;
                X0[ 4 * _dim + i ] = 0.5;

            }

            error = solverTools::newtonRaphson(func, X0, X, converged, fatalError, floatOuts, intOuts, floatArgs, intArgs );

            if ( fatalError ){

                errorOut result = new errorNode( "computeBoundaryPoints",
                                                 "Fatal error in Newton-Raphson solve" );
                result->addNext( error );
                return result;

            
            }
            if ( !converged ){

                localBoundaryPoint = floatVector( X0.begin( ), X0.begin( ) + _dim );
                element->interpolate( element->reference_nodes, localBoundaryPoint, boundaryPoint );

            }
            else{
                
                localBoundaryPoint = floatVector( X.begin( ), X.begin( ) + _dim );
                element->interpolate( element->reference_nodes, localBoundaryPoint, boundaryPoint );

            }

            for ( uIntType i = 0; i < _dim; i++ ){

                _boundaryPoints.push_back( boundaryPoint[ i ] );

            }

            _boundaryPointIDToIndex.emplace( *bc, bc - _boundaryCells.begin( ) );
            ownedIndices[ bc - _boundaryCells.begin( ) ] = _dim * ( bc - _boundaryCells.begin( ) );

        }

        //Form the KD tree of the boundary points
        _boundaryPointTree = KDNode( &_boundaryPoints, ownedIndices, 0, _dim );

        return NULL;
    }

    errorOut dualContouringInternalPointResidual( const floatVector &X, const floatMatrix &floatArgs,
                                                  const intMatrix &intArgs,
                                                  floatVector &residual, floatMatrix &jacobian,
                                                  floatMatrix &floatOuts, intMatrix &intOuts ){
        /*!
         * The residual equation for the computation of the internal point for a boundary cell
         * in the dual contouring method.
         *
         * :param floatVector &X: The solution vector. Ordered as [ x, s, t, lambda_ub, lambda_lb ]
         * :param floatVector &floatArgs: The floating point arguments. Ordered as
         *     [ [ x_ub ], [ x_lb] , [ p1 ], [ p2 ], ... , [ n1 ], [ n2 ], ... ]
         * :param intVector &intArgs: The integer arguments
         *     [ [ dim, nPoints ] ]
         * :param floatVector &residual: The residual vector
         * :param floatMatrix &jacobian: The jacobian matrix
         * :param floatMatrix &floatOuts: Not used
         * :param intMatrix &intOuts: Not used 
         */

        ( void ) floatOuts;
        ( void ) intOuts;

        if ( intArgs.size( ) != 1 ){

            return new errorNode( "internalPointResidual", "The intArgs matrix must have one element" );

        }

        if ( intArgs[ 0 ].size( ) != 2 ){
                
            return new errorNode( "internalPointResidual", "The first value of intArgs must have a length of 2" );

        }


        uIntType dim = intArgs[ 0 ][ 0 ];
        uIntType nPoints = intArgs[ 0 ][ 1 ];

        if ( X.size ( ) != 5 * dim ){

            return new errorNode( "internalPointResidual", "The 'X' vector must have a length of 5 times the dimension" );

        }

        if ( floatArgs.size( ) != ( 2 + 2 * nPoints ) ){

            return new errorNode( "internalPointResidual",
                                  "The floatArgs matrix must have " + std::to_string( 2 + 2 * nPoints ) + " elements" );

        }
        
        //Extract the values from X
        floatVector x( X.begin( ), X.begin( ) + dim );
        floatVector s( X.begin( ) + dim, X.begin( ) + 2 * dim );
        floatVector t( X.begin( ) + 2 * dim, X.begin( ) + 3 * dim );
        floatVector lub( X.begin( ) + 3 * dim, X.begin( ) + 4 * dim );
        floatVector llb( X.begin( ) + 4 * dim, X.begin( ) + 5 * dim );

        //Extract the values from floatArgs
        floatVector xub = floatArgs[ 0 ];
        floatVector xlb = floatArgs[ 1 ];

        floatMatrix points( floatArgs.begin( ) + 2, floatArgs.begin( ) + 2 + nPoints );
        floatMatrix normals( floatArgs.begin( ) + 2 + nPoints, floatArgs.begin( ) + 2 + 2 * nPoints );

        //Form the residual vector and the Jacobian
        residual = floatVector( 5 * dim, 0 );
        jacobian = floatMatrix( 5 * dim, floatVector( 5 * dim, 0 ) );

        for ( uIntType i = 0; i < nPoints; i++ ){

            //Add the contribution to the first residual
            floatType nxmp = vectorTools::dot( normals[ i ], x - points[ i ] );

            for ( uIntType _i = 0; _i < dim; _i++ ){

                residual[ _i ] += nxmp * normals[ i ][ _i ] + x[ _i ];
                jacobian[ _i ][ _i ] += 1;

                for ( uIntType _j = 0; _j < dim; _j++ ){

                    jacobian[ _i ][ _j ] += normals[ i ][ _i ] * normals[ i ][ _j ];

                }

            }

        }

        for ( uIntType i = 0; i < dim; i++ ){

            //Add the terms to the residual
            residual[            i ] +=  lub[ i ] - llb[ i ];
            residual[     dim + i ]  =  2 * lub[ i ] * s[ i ];
            residual[ 2 * dim + i ]  = -2 * llb[ i ] * t[ i ];
            residual[ 3 * dim + i ]  = xub[ i ] - x[ i ] - s[ i ] * s[ i ];
            residual[ 4 * dim + i ]  = x[ i ] - xlb[ i ] - t[ i ] * t[ i ];

            //Assemble the Jacobian
            
            //Remaining jacobians of the residual of the first term
            jacobian[ i ][ 3 * dim + i ] =  1.;
            jacobian[ i ][ 4 * dim + i ] = -1.;

            //Remaining jacobians of the residual of the second term
            jacobian[ dim + i ][     dim + i ] = 2 * lub[ i ];
            jacobian[ dim + i ][ 3 * dim + i ] = 2 * s[ i ];

            //Remaining jacobians of the residual of the third term
            jacobian[ 2 * dim + i ][ 2 * dim + i ] = -2 * llb[ i ];
            jacobian[ 2 * dim + i ][ 4 * dim + i ] = -2 * t[ i ];

            //Remaining jacobians of the residual of the fourth term
            jacobian[ 3 * dim + i ][        i ] = -1;
            jacobian[ 3 * dim + i ][ dim + i ] = -2 * s[ i ];

            //Remaining jacobians of the residual of the fifth term
            jacobian[ 4 * dim + i ][        i ]     =  1;
            jacobian[ 4 * dim + i ][ 2 * dim + i ] = -2 * t[ i ];

        }

        return NULL;
    }

    errorOut dualContouring::writeToXDMF( ){
        /*!
         * write the data to an XDMF file
         */

        //Build the XDMF file
        shared_ptr< XdmfDomain > _domain = XdmfDomain::New( );
        shared_ptr< XdmfInformation > domainInfo
            = XdmfInformation::New( "Domain", "Primary data structure from a volume reconstruction object" );

        _domain->insert( domainInfo );

        shared_ptr< XdmfHDF5Writer > heavyWriter = XdmfHDF5Writer::New( _XDMFOutputFilename + ".h5", true );
        heavyWriter->setReleaseData( true );
        shared_ptr< XdmfWriter > writer = XdmfWriter::New( _XDMFOutputFilename + ".xdmf", heavyWriter );

        //Initialize the grid collection
        shared_ptr< XdmfGridCollection > _gridCollection = XdmfGridCollection::New( );
        _gridCollection->setType( XdmfGridCollectionType::Spatial( ) );
        shared_ptr< XdmfInformation > _gridCollectionInfo = XdmfInformation::New( "Grid Collection", "The collection of grids used in the formation of the reconstructed domain" );
        _gridCollection->insert( _gridCollectionInfo );
        _domain->insert( _gridCollection );

        //Write the source points out
        
        shared_ptr< XdmfUnstructuredGrid > _sourceNodeGrid = XdmfUnstructuredGrid::New( );
        _sourceNodeGrid->setName( "Source Node Grid" );

        //Set the source node geometry
        shared_ptr< XdmfGeometry > _sourceNodeGeometry = XdmfGeometry::New( );
        _sourceNodeGeometry->setType( XdmfGeometryType::XYZ( ) );
        _sourceNodeGeometry->setName( "Source Node Coordinates" );
        _sourceNodeGeometry->insert( 0, getPoints( )->data( ), ( uIntType )3 * _nPoints, 1, 1 );
        shared_ptr< XdmfInformation > _sourceNodeGeometryInfo = XdmfInformation::New( "Source Node Coordinates", "The coordinates of the source nodes ( i.e. the points to be reconstructed ) in x1, y1, z1, x2, ... format" );
        _sourceNodeGrid->setGeometry( _sourceNodeGeometry );

        //Set the source node topology
        shared_ptr< XdmfTopology > _sourceNodeTopology = XdmfTopology::New( );
        _sourceNodeTopology->setType( XdmfTopologyType::Polyvertex( ) );
        _sourceNodeTopology->setName( "Source Node Topology" );
        uIntVector sourceNodeIds( getPoints( )->size( ) );
        for ( uIntType i = 0; i < _nPoints; i++ ){
            sourceNodeIds[ i ] = i;
        } 
        _sourceNodeTopology->insert( 0, sourceNodeIds.data( ), _nPoints, 1, 1 );
        _sourceNodeGrid->setTopology( _sourceNodeTopology );

        //Group the IDs into a set
        shared_ptr< XdmfSet > _sourceNodeSet = XdmfSet::New( );
        _sourceNodeSet->setType( XdmfSetType::Node( ) );
        _sourceNodeSet->setName( "Source Nodes" );
        _sourceNodeSet->insert( 0, sourceNodeIds.data( ), _nPoints, 1, 1 );
        _sourceNodeGrid->insert( _sourceNodeSet );

        _gridCollection->insert( _sourceNodeGrid );

        //Write the boundary points out

        shared_ptr< XdmfUnstructuredGrid > _boundaryPointGrid = XdmfUnstructuredGrid::New ( );
        _boundaryPointGrid->setName( "Boundary Point Grid" );

        //Set the boundary surface geometry
        shared_ptr< XdmfGeometry > _boundaryPointGeometry = XdmfGeometry::New( );
        _boundaryPointGeometry->setType( XdmfGeometryType::XYZ( ) );
        _boundaryPointGeometry->setName( "Boundary Surface Coordinates" );
        _boundaryPointGeometry->insert( 0, _boundaryPoints.data( ), _boundaryPoints.size( ), 1, 1 );
        shared_ptr< XdmfInformation > _boundaryPointsInfo = XdmfInformation::New( "Boundary Surface Coordinates", "The coordinates of the boundary points ( i.e. the points which are joined together to form the mesh ) in x1, y1, z1, x2, ... format" );
        _boundaryPointGeometry->insert( _boundaryPointsInfo );
        _boundaryPointGrid->setGeometry( _boundaryPointGeometry );

        //Set the map from the boundary ID to the XDMF ID

        //Set the boundary surface topology
        shared_ptr< XdmfTopology > _boundaryPointTopology = XdmfTopology::New( );
        _boundaryPointTopology->setType( XdmfTopologyType::Quadrilateral( ) );

        uIntVector _boundaryPointConnectivity;
        _boundaryPointConnectivity.reserve( 4 * ( _boundaryEdges_x.size( ) + _boundaryEdges_y.size( ) + _boundaryEdges_z.size( ) ) );

        for ( auto it = _boundaryEdges_x.begin( ); it != _boundaryEdges_x.end( ); it++ ){

            for ( uIntType i = 0; i < it->second.size( ); i++ ){
                _boundaryPointConnectivity.push_back( _boundaryPointIDToIndex[ it->second[ i ] ] );
            }

        }

        for ( auto it = _boundaryEdges_y.begin( ); it != _boundaryEdges_y.end( ); it++ ){

            for ( uIntType i = 0; i < it->second.size( ); i++ ){
                _boundaryPointConnectivity.push_back( _boundaryPointIDToIndex[ it->second[ i ] ] );
            }

        }

        for ( auto it = _boundaryEdges_z.begin( ); it != _boundaryEdges_z.end( ); it++ ){

            for ( uIntType i = 0; i < it->second.size( ); i++ ){
                _boundaryPointConnectivity.push_back( _boundaryPointIDToIndex[ it->second[ i ] ] );
            }

        }

        _boundaryPointTopology->insert( 0, _boundaryPointConnectivity.data( ), _boundaryPointConnectivity.size( ), 1, 1 );
        shared_ptr< XdmfInformation > _boundaryPointTopologyInfo = XdmfInformation::New( "Boundary Surface Connectivity", "The connectivity of the boundary points" );
        _boundaryPointTopology->insert( _boundaryPointTopologyInfo );
        _boundaryPointGrid->setTopology( _boundaryPointTopology );

        //Export the normals of the boundary cells
        shared_ptr< XdmfAttribute > _boundaryPointNormalsPtr = XdmfAttribute::New( );
        _boundaryPointNormalsPtr->setType( XdmfAttributeType::Vector( ) );
        _boundaryPointNormalsPtr->setCenter( XdmfAttributeCenter::Node( ) );
        _boundaryPointNormalsPtr->setName( "Boundary Point Normal" );

        floatVector _boundaryPointNormalVector;
        _boundaryPointNormalVector.reserve( _dim * _boundaryPointNormals.size( ) );

        for ( auto it = _boundaryPointNormals.begin( ); it != _boundaryPointNormals.end( ); it++ ){

            for ( uIntType i = 0; i < it->second.size( ); i++ ){

                _boundaryPointNormalVector.push_back( it->second[ i ] );

            }

        }

        _boundaryPointNormalsPtr->insert( 0, _boundaryPointNormalVector.data( ), _boundaryPointNormalVector.size( ), 1, 1 );

        shared_ptr< XdmfInformation > _boundaryPointNormalInformation
            = XdmfInformation::New( "Boundary Point Normal", "The average normals at the boundary points" );

        _boundaryPointNormalsPtr->insert( _boundaryPointNormalInformation );

        _boundaryPointGrid->insert( _boundaryPointNormalsPtr );

        _gridCollection->insert( _boundaryPointGrid );

        //Write the output file
        _domain->accept( writer );

        return NULL;
    }

    errorOut dualContouring::computeBoundaryPointNormalsAndAreas( ){
        /*!
         * Compute the normals and surface areas at the boundary points
         */

        if ( _dim != 3 ){
            return new errorNode( "computeBoundaryPointNormals", "This function requires the dimension is 3" );
        }

        _boundaryPointAreas.reserve( _boundaryPoints.size( ) / 3 );
        _boundaryPointNormals.reserve( _boundaryPoints.size( ) );

        //Loop through the edges

        errorOut error = processBoundaryEdges( _boundaryEdges_x );

        if ( error ){

            errorOut result = new errorNode( "computeBoundaryPointNormalsAndAreas",
                                             "Error in processing the x boundary edges" );
            result->addNext( error );
            return result;

        }

        error = processBoundaryEdges( _boundaryEdges_y );

        if ( error ){

            errorOut result = new errorNode( "computeBoundaryPointNormalsAndAreas",
                                             "Error in processing the y boundary edges" );
            result->addNext( error );
            return result;

        }

        error = processBoundaryEdges( _boundaryEdges_z );

        if ( error ){

            errorOut result = new errorNode( "computeBoundaryPointNormalsAndAreas",
                                             "Error in processing the z boundary edges" );
            result->addNext( error );
            return result;

        }

        //Normalize the vectors
        for ( auto it = _boundaryPointAreas.begin( ); it != _boundaryPointAreas.end( ); it++ ){

            _boundaryPointNormals[ it->first ] /= it->second;

        }

        return NULL;
    }


    errorOut dualContouring::processBoundaryEdges( const std::unordered_map< uIntType, uIntVector > &boundaryEdges ){
        /*!
         * Process the boundary edges contributions to the normal and surface area vectors
         *
         * :param const std::unordered_map< uIntType, uIntVector > &boundaryEdges: The boundary edges
         */

        //Loop through the edges
        floatVector n, p1, p2, p3, p4;
        floatType a;

        for ( auto edge = boundaryEdges.begin( ); edge != boundaryEdges.end( ); edge++ ){

            //Get the points
            auto index = _boundaryPointIDToIndex.find( edge->second[ 0 ] );
            if ( index == _boundaryPointIDToIndex.end( ) ){

                return new errorNode( "processBoundaryEdges", "Edge boundary point ID " + std::to_string( edge->second[ 0 ] )
                                      + " not found in boundary point ID to index map." );

            }
            p1 = floatVector( _boundaryPoints.begin( ) + _dim * ( index->second ),
                              _boundaryPoints.begin( ) + _dim * ( index->second + 1 ) );

            index = _boundaryPointIDToIndex.find( edge->second[ 1 ] );
            if ( index == _boundaryPointIDToIndex.end( ) ){

                return new errorNode( "processBoundaryEdges", "Edge boundary point ID " + std::to_string( edge->second[ 1 ] )
                                      + " not found in boundary point ID to index map." );

            }
            p2 = floatVector( _boundaryPoints.begin( ) + _dim * ( index->second ),
                              _boundaryPoints.begin( ) + _dim * ( index->second + 1 ) );

            index = _boundaryPointIDToIndex.find( edge->second[ 2 ] );
            if ( index == _boundaryPointIDToIndex.end( ) ){

                return new errorNode( "processBoundaryEdges", "Edge boundary point ID " + std::to_string( edge->second[ 2 ] )
                                      + " not found in boundary point ID to index map." );

            }
            p3 = floatVector( _boundaryPoints.begin( ) + _dim * ( index->second ),
                              _boundaryPoints.begin( ) + _dim * ( index->second + 1 ) );

            index = _boundaryPointIDToIndex.find( edge->second[ 3 ] );
            if ( index == _boundaryPointIDToIndex.end( ) ){

                return new errorNode( "processBoundaryEdges", "Edge boundary point ID " + std::to_string( edge->second[ 3 ] )
                                      + " not found in boundary point ID to index map." );

            }
            p4 = floatVector( _boundaryPoints.begin( ) + _dim * ( index->second ),
                              _boundaryPoints.begin( ) + _dim * ( index->second + 1 ) );

            //Add the points to the area and normal vectors if required
            for ( uIntType i = 0; i < edge->second.size( ); i++ ){

                if ( _boundaryPointAreas.find( edge->second[ i ] ) == _boundaryPointAreas.end( ) ){
    
                    _boundaryPointAreas.emplace( edge->second[ i ], 0. );
                    _boundaryPointNormals.emplace( edge->second[ i ], floatVector( _dim, 0. ) );
    
                }

            }

            //Compute the first triangle's normal and area
            n = vectorTools::cross( p2 - p1, p4 - p1 );
            a = 0.5 * vectorTools::l2norm( n );
            n /= ( 2 * a );

            for ( uIntType i = 0; i < edge->second.size( ); i++ ){

                _boundaryPointAreas[ edge->second[ i ] ] += 0.25 * a;
                _boundaryPointNormals[ edge->second[ i ] ] += 0.25 * a * n;

            }

            //Compute the second triangle's normal and area
            n = vectorTools::cross( p4 - p3, p2 - p3 );
            a = 0.5 * vectorTools::l2norm( n );
            n /= ( 2 * a );

            for ( uIntType i = 0; i < edge->second.size( ); i++ ){

                _boundaryPointAreas[ edge->second[ i ] ] += 0.25 * a;
                _boundaryPointNormals[ edge->second[ i ] ] += 0.25 * a * n;

            }

        }

        return NULL;
    }

    errorOut dualContouring::interpolateFunctionToBackgroundGrid( const floatVector &functionValuesAtPoints,
                                                                  const uIntType &functionDim,
                                                                  std::unordered_map< uIntType, floatVector > &functionAtGrid
                                                                ){ 
        /*!
         * Interpolate a function defined at the data points to the given element in the background grid.
         * 
         * :param const floatVector &functionValuesAtPoints: The value of the function at the data points
         * :param const uIntType &functionDim: The dimensionality of the function e.g. a function defined by four scalar values
         *     would have a dimensionality of four.
         * :param std::unordered_map< uIntType, floatVector > &functionAtGrid: The value of the function
         *     projected to the nodes in the background grid
         */

        //Error handling
        if ( _points->size( ) / _dim != functionValuesAtPoints.size( ) / functionDim ){
            return new errorNode( "interpolateFunctionToBackgroundGrid",
                                  "The points vector and the function values at points vector are not of compatible sizes" );
        }

        bool usePointwiseProjection = false;

        //Reserve memory for the function at the grid. This will be a worst case.
        functionAtGrid.clear( );
        functionAtGrid.reserve( 8 * functionDim * _internalCells.size( ) );

        std::unordered_map< uIntType, floatType > weights;
        weights.reserve( 8 * _internalCells.size( ) );

        //Get the number of values in the different directions
        uIntType ngy = _gridLocations[ 1 ].size( );
        uIntType ngz = _gridLocations[ 2 ].size( );

        uIntType i, j, k;
        uIntVector indices;
        uIntVector internalPoints;

        std::unique_ptr< elib::Element > element;

        floatVector pointPosition;
        floatVector localCoordinates;
        floatVector shapeFunctions;
        floatVector functionValue;

        errorOut error;

        floatVector upperBounds;
        floatVector lowerBounds;

        for ( auto cell = _internalCells.begin( ); cell != _internalCells.end( ); cell++ ){

            //Get the bottom corner node IDs from the cell id
            i = *cell / ( ngy * ngz );
            j = ( *cell - ngy * ngz * i ) / ngz;
            k = (*cell - ngy * ngz * i - ngz * j );

            indices = { i, j, k };

            //Get the element associated with this cell
            error = getGridElement( indices, element );

            if ( error ){

                errorOut result = new errorNode( "interpolateFunctionToBackgroundGrid",
                                                 "Error in getting the grid element" );
                result->addNext( error );
                return result;

            }

            //Add the global node ids to the map if required
            for ( auto nID = element->global_node_ids.begin( ); nID  != element->global_node_ids.end( ); nID++ ){

                if ( weights.find( *nID ) == weights.end( ) ){

                    functionAtGrid.emplace( *nID, floatVector( functionDim, 0 ) );
                    weights.emplace( *nID, 0 );

                }

            }

            for ( auto node = element->nodes.begin( ); node != element->nodes.end( ); node++ ){

                uIntType   globalNodeID = element->global_node_ids[ node - element->nodes.begin( ) ];
                uIntVector internalNodes;
                _pointTree.getPointsWithinRadiusOfOrigin( *node, _critical_radius, internalNodes );

                for ( auto iN = internalNodes.begin( ); iN != internalNodes.end( ); iN++ ){

                    functionValue = floatVector( functionValuesAtPoints.begin( ) + ( *iN / _dim + 0 ) * functionDim,
                                                 functionValuesAtPoints.begin( ) + ( *iN / _dim + 1 ) * functionDim );

                    pointPosition = floatVector( _points->begin( ) + *iN,
                                                 _points->begin( ) + *iN + _dim );

                    floatType value;
                    rbf( *node, pointPosition, _length_scale, value );

                    functionAtGrid[ globalNodeID ] += value * functionValue;

                    weights[ globalNodeID ] += value;

                }

            }

//            getPointsWithinRadiusOfOrigin( currentSeedPoint, minDistance, internalNodes )
//
//
//            upperBounds = *getUpperBounds( );
//            lowerBounds = *getLowerBounds( );
//
//            //Find the points inside of this element
//            internalPoints.clear( );
//            _pointTree.getPointsInRange( element->bounding_box[ 1 ], element->bounding_box[ 0 ], internalPoints,
//                                         &upperBounds, &lowerBounds );
//
//            //Loop over the points adding their contributions to the function at the current grid points
//
//            std::vector< bool > isDuplicate( internalPoints.size( ), false );
//            floatVector meanFunctionValue( functionDim, 0 );
//            uIntType numPoints = 0;
//            for ( auto p = internalPoints.begin( ); p != internalPoints.end( ); p++ ){
//
//                unsigned int pindex = p - internalPoints.begin( );
//
//                if ( isDuplicate[ pindex ] ){
//                    continue;
//                }
//
//                unsigned int duplicate_count = 0;
//
//                //Get the function value at the current point
//                functionValue = floatVector( functionValuesAtPoints.begin( ) + ( *p / _dim + 0 ) * functionDim,
//                                             functionValuesAtPoints.begin( ) + ( *p / _dim + 1 ) * functionDim );
//
//                pointPosition = floatVector( _points->begin( ) + *p,
//                                             _points->begin( ) + *p + _dim );
//
//                // Check if any of the remaining points are duplicates
//                for ( auto q = internalPoints.begin( ) + pindex + 1; q != internalPoints.end( ); q++ ){
//
//                    unsigned int qindex = q - internalPoints.begin( );
//
//                    floatVector q_pointPosition = floatVector( _points->begin( ) + *q,
//                                                               _points->begin( ) + *q + _dim );
//
//                    if ( vectorTools::l2norm( pointPosition - q_pointPosition ) < _absoluteTolerance ){
//
//                        isDuplicate[ qindex ] = true;
//
//                        functionValue += floatVector( functionValuesAtPoints.begin( ) + ( *q / _dim + 0 ) * functionDim,
//                                                      functionValuesAtPoints.begin( ) + ( *q / _dim + 1 ) * functionDim );
//
//                        duplicate_count++;
//
//                    }
//
//                }
//
//                functionValue /= ( duplicate_count + 1 );
//
//                //Get the local coordinates of the point
//                if ( usePointwiseProjection ){
//
//                    error = element->compute_local_coordinates( pointPosition, localCoordinates );
//
//                    if ( error ){
//
//                        errorOut result = new errorNode( __func__,
//                                                         "Error in the computation of the local coordinates" );
//                        result->addNext( error );
//                        return result;
//
//                    }
//
//                    //Get the shape function values at the point
//                    error = element->get_shape_functions( localCoordinates, shapeFunctions );
//
//                    if ( error ){
//                    
//                        errorOut result = new errorNode( __func__, "Error in the computation of the shape functions" );
//                        result->addNext( error );
//                        return result;
//                    }
//
//                }
//
//                //Compute the contribution of the nodes to the background grid
//
//                if ( usePointwiseProjection ){
//                    for ( auto nID = element->global_node_ids.begin( ); nID != element->global_node_ids.end( ); nID++ ){
//                        functionAtGrid[ *nID ] += shapeFunctions[ nID - element->global_node_ids.begin( ) ] * functionValue;
//                        weights[ *nID ] += shapeFunctions[ nID - element->global_node_ids.begin( ) ];
//                    }
//                }
//                else{
//                    meanFunctionValue += functionValue;
//                    numPoints++;
//                }
//
//            }
//            if ( ( !usePointwiseProjection ) && ( internalPoints.size( ) > 0 ) ){
//                for ( auto nID = element->global_node_ids.begin( ); nID != element->global_node_ids.end( ); nID++ ){
//                    functionAtGrid[ *nID ] += meanFunctionValue / std::fmax( ( floatType )numPoints, 1 );
//                    weights[ *nID ] += 1;
//                }
//            }
        }

        //Normalize the weights by the shape function values
        for ( auto node = functionAtGrid.begin( ); node != functionAtGrid.end( ); node++ ){

            if ( weights[ node->first ] > _absoluteTolerance ){

                functionAtGrid[ node->first ] = node->second / weights[ node->first ];

            }

        }

        return NULL;

    }

    errorOut dualContouring::performVolumeIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                       floatVector &integratedValue ){
        /*!
         * Integrate a quantity known at the points over the volume returning the value for the domain.
         *
         * :param const floatVector &valuesAtPoints: A vector of the values at the data points. Stored as
         *     [ v_11, v_12, ..., v_21, v22, ... ] where the first index is the point index in order as 
         *     provided to the volume reconstruction object and the second index is the value of the 
         *     function to be integrated.
         * :param const uIntType valueSize: The size of the subvector associated with each of the datapoints.
         * :param floatVector &integratedValue: The final value of the integral
         */

        errorOut error;

        //Check if the domain has been constructed yet
        if ( !getEvaluated( ) ){

            error = evaluate( );

            if ( error ){

                errorOut result = new errorNode( "performVolumeIntegration",
                                                 "Error encountered during the reconstruction of the volume" );
                result->addNext( error );
                return result;

            }

        }

        //Interpolate the function to the background grid
        std::unordered_map< uIntType, floatVector > functionAtGrid;
        error = interpolateFunctionToBackgroundGrid( valuesAtPoints, valueSize, functionAtGrid );

        if ( error ){

            errorOut result = new errorNode( "performVolumeIntegration",
                                             "Error encountered during the interpolation of the function to the background grid" );
            result->addNext( error );
            return result;

        }

        //Perform the volume integration
        integratedValue = floatVector( valueSize, 0. );
        floatMatrix jacobian;
        floatType Jxw;
        floatVector qptValue;

        //Get the number of values in the different directions
        uIntType ngy = _gridLocations[ 1 ].size( );
        uIntType ngz = _gridLocations[ 2 ].size( );
        uIntType i, j, k;
        uIntVector indices;

        std::unique_ptr< elib::Element > element;

        floatMatrix nodalFunctionValues;
        floatType fVal;

        for ( auto cell = _internalCells.begin( ); cell != _internalCells.end( ); cell++ ){

            //Get the bottom corner node IDs from the cell id
            i = *cell / ( ngy * ngz );
            j = ( *cell - ngy * ngz * i ) / ngz;
            k = (*cell - ngy * ngz * i - ngz * j );

            indices = { i, j, k };

            //Get the element associated with this cell
            error = getGridElement( indices, element );

            if ( error ){

                errorOut result = new errorNode( "performVolumeIntegration",
                                                 "Error in getting the grid element" );
                result->addNext( error );
                return result;

            }
            
            //Extract the function at the nodes
            nodalFunctionValues = floatMatrix( element->global_node_ids.size( ), floatVector( valueSize ) );
            for ( auto nID = element->global_node_ids.begin( ); nID != element->global_node_ids.end( ); nID++ ){

                if ( functionAtGrid.find( *nID ) == functionAtGrid.end( ) ){

                    return new errorNode( "performVolumeIntegration",
                                          "Node with global ID " + std::to_string( *nID )
                                        + " not found in the grid node to function map" );

                }

                //Get the value of the function at the nodal values multiplied by whether the node is
                //within the domain or not

                if ( *nID > _implicitFunctionValues.size( ) ){

                    return new errorNode( "performVolumeIntegration",
                                          "The nodal ID is too large for the implicit function values vector\n nID: "
                                          + std::to_string( *nID ) );

                }

                
                fVal = _implicitFunctionValues[ *nID ];

                nodalFunctionValues[ nID - element->global_node_ids.begin( ) ]
                    = ( floatType )( fVal > 0 ) * functionAtGrid[ *nID ];

            }

            //Perform the integration
            for ( auto qpt = element->qrule.begin( ); qpt != element->qrule.end( ); qpt++ ){

                element->interpolate( nodalFunctionValues, qpt->first, qptValue );
                element->get_local_gradient( element->reference_nodes, qpt->first, jacobian );

                Jxw = vectorTools::determinant( vectorTools::appendVectors( jacobian ), _dim, _dim ) * qpt->second;
                integratedValue += qptValue * Jxw;

            }

        }

        return NULL;

    }

    errorOut dualContouring::performSurfaceIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                        floatVector &integratedValue, const uIntVector *subdomainIDs,
                                                        const floatVector *subdomainWeights,
                                                        const floatVector *macroNormal, const bool useMacroNormal ){
        /*!
         * Integrate a quantity known at the point over the surface return the value for the domain.
         *
         * :param const floatVector &valuesAtPoints: A vector of the values at the data points. Stored as
         *     [ v_11, v_12, ..., v_21, v22, ... ] where the first index is the point index in order as 
         *     provided to the volume reconstruction object and the second index is the value of the 
         *     function to be integrated.
         * :param const uIntType valueSize: The size of the subvector associated with each of the datapoints.
         * :param floatVector &integratedValue: The final value of the integral
         * :param const uIntVector *subdomainIDs: The pointer to the subdomain of the surface to integrate over
         * :param const floatVector *subdomainWeights: The weights for the subdomains. Useful if points can be
         *     in multiple subdomains and they aren't small w.r.t. the domain size
         * :param const floatVector *macroNormal: A macro-scale normal vector to use to generate the micro
         *     weight. This can be helpful in cases where some points start to, ``wrap,'' around an edge 
         *     which should be flat. Can either be a single vector of dimension _dim or a collection of vectors
         *     at each boundary point.
         * :param const bool useMacroNormal: Use the macro-scale normal instead of the micro normals. Can help
         *     drive the integral to be more what is expected in some cases.
         */

        floatVector origin;
        errorOut error = performSurfaceIntegralMethods( valuesAtPoints, valueSize, origin, integratedValue, false, false, false,
                                                        subdomainIDs, subdomainWeights, macroNormal, useMacroNormal );

        if ( error ){

            errorOut result = new errorNode( "performSurfaceIntegration",
                                             "Error in the computation of the surface integral" );
            result->addNext( error );
            return result;

        }

        return NULL;
    }

    errorOut dualContouring::performPositionWeightedSurfaceIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                                        floatVector &integratedValue, const uIntVector *subdomainIDs,
                                                                        const floatVector *subdomainWeights,
                                                                        const floatVector *macroNormal, const bool useMacroNormal ){
        /*!
         * Integrate a quantity known at the point over the surface return the value for the domain.
         *
         * :param const floatVector &valuesAtPoints: A vector of the values at the data points. Stored as
         *     [ v_11, v_12, ..., v_21, v22, ... ] where the first index is the point index in order as 
         *     provided to the volume reconstruction object and the second index is the value of the 
         *     function to be integrated.
         * :param const uIntType valueSize: The size of the subvector associated with each of the datapoints.
         * :param floatVector &integratedValue: The final value of the integral organized as
         *     [ v_11, v_12, v_13, ..., v_21, v_22, ... ] where the first index is the value and the second is the
         *     dimension
         * :param const uIntVector *subdomainIDs: The pointer to the subdomain of the surface to integrate over
         * :param const floatVector *subdomainWeights: The weights for the subdomains. Useful if points can be
         *     in multiple subdomains and they aren't small w.r.t. the domain size
         * :param const floatVector *macroNormal: A macro-scale normal vector to use to generate the micro
         *     weight. This can be helpful in cases where some points start to, ``wrap,'' around an edge 
         *     which should be flat. Can either be a single vector of dimension _dim or a collection of vectors
         *     at each boundary point.
         * :param const bool useMacroNormal: Use the macro-scale normal instead of the micro normals. Can help
         *     drive the integral to be more what is expected in some cases.
         */

        floatVector origin;
        errorOut error = performSurfaceIntegralMethods( valuesAtPoints, valueSize, origin, integratedValue, false, true, false,
                                                        subdomainIDs, subdomainWeights, macroNormal, useMacroNormal );

        if ( error ){

            errorOut result = new errorNode( "performPositionWeightedSurfaceIntegration",
                                             "Error in the computation of the surface integral" );
            result->addNext( error );
            return result;

        }

        return NULL;
    }


    errorOut dualContouring::performSurfaceFluxIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                            floatVector &integratedValue, const uIntVector *subdomainIDs,
                                                            const floatVector *subdomainWeights,
                                                            const floatVector *macroNormal, const bool useMacroNormal ){
        /*!
         * Integrate the flux of a quantity known at the point over the surface return the value for the domain.
         *
         * :param const floatVector &valuesAtPoints: A vector of the values at the data points. Stored as
         *     [ v_11, v_12, ..., v_21, v22, ... ] where the first index is the point index in order as 
         *     provided to the volume reconstruction object and the second index is the value of the 
         *     function to be integrated.
         * :param const uIntType valueSize: The size of the subvector associated with each of the datapoints.
         * :param floatVector &integratedValue: The final value of the integral
         * :param const uIntVector *subdomainIDs: The pointer to the subdomain of the surface to integrate over
         * :param const floatVector *subdomainWeights: The weights for the subdomains. Useful if points can be
         *     in multiple subdomains and they aren't small w.r.t. the domain size
         * :param const floatVector *macroNormal: A macro-scale normal vector to use to generate the micro
         *     weight. This can be helpful in cases where some points start to, ``wrap,'' around an edge 
         *     which should be flat. Can either be a single vector of dimension _dim or a collection of vectors
         *     at each boundary point.
         * :param const bool useMacroNormal: Use the macro-scale normal instead of the micro normals. Can help
         *     drive the integral to be more what is expected in some cases.
         */

        floatVector origin;
        errorOut error = performSurfaceIntegralMethods( valuesAtPoints, valueSize, origin, integratedValue, true, false, false, subdomainIDs,
                                                        subdomainWeights, macroNormal, useMacroNormal );

        if ( error ){

            errorOut result = new errorNode( "performSurfaceIntegration",
                                             "Error in the computation of the surface integral" );
            result->addNext( error );
            return result;

        }

        return NULL;
    }

    errorOut dualContouring::performRelativePositionSurfaceFluxIntegration( const floatVector &valuesAtPoints,
                                                                            const uIntType valueSize,
                                                                            const floatVector &origin,
                                                                            floatVector &integratedValue,
                                                                            const uIntVector *subdomainIDs,
                                                                            const floatVector *subdomainWeights,
                                                                            const floatVector *macroNormal,
                                                                            const bool useMacroNormal ){
        /*!
         * $\int_{\partial\mathcal{B}} n_i v_ij ( x_k' - o_k ) da \approx \sum_{p = 1}^N n_i^p v_ij^p ( x_k' - o_k ) da^p$
         *
         * where $x_k'$ is the position of the point on the surface, $o_k$ is the origin, and $n_i$ is the normal at the point.
         *
         * :param const floatVector &valuesAtPoints: A vector of the values at the data points. Stored as
         *     [ v_111, v_112, v_113, ..., v_121, v_122, v_123, ... ] where the first index is the point index in order as 
         *     provided to the volume reconstruction object, the second index is the first index of the v matirx and the
         *     third index is the second index of the value matrix
         *     function to be integrated.
         * :param const uIntType valueSize: The size of the subvector associated with each of the datapoints.
         * :param const floatVector &origin: The origin that the relative position vector is computed in relation to.
         * :param floatVector &integratedValue: The final value of the integral
         * :param const uIntVector *subdomainIDs: The IDs of points in the subdomain to integrate over
         * :param const floatVector *subdomainWeights: The weights for the subdomains. Useful if points can be
         *     in multiple subdomains and they aren't small w.r.t. the domain size
         * :param const floatVector *macroNormal: A macro-scale normal vector to use to generate the micro
         *     weight. This can be helpful in cases where some points start to, ``wrap,'' around an edge 
         *     which should be flat. Can either be a single vector of dimension _dim or a collection of vectors
         *     at each boundary point.
         * :param const bool useMacroNormal: Use the macro-scale normal instead of the micro normals. Can help
         *     drive the integral to be more what is expected in some cases.
         */

        errorOut error = performSurfaceIntegralMethods( valuesAtPoints, valueSize, origin, integratedValue, true, false, true, subdomainIDs,
                                                        subdomainWeights, macroNormal, useMacroNormal );

        if ( error ){

            errorOut result = new errorNode( "performRelativePositionSurfaceFluxIntegration",
                                             "Error in computation of the integral of the dyadic product between a flux and the relative position vector" );
            result->addNext( error );
            return result;

        }

        return NULL;
    }
                    

    errorOut dualContouring::performSurfaceIntegralMethods( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                            const floatVector &origin, floatVector &integratedValue,
                                                            bool computeFlux, bool positionWeightedIntegral, bool dyadWithOrigin,
                                                            const uIntVector *subdomainIDs,
                                                            const floatVector *subdomainWeights,
                                                            const floatVector *macroNormal, const bool useMacroNormal ){
        /*!
         * Integrate a quantity known at the points over the surface return the value for the domain.
         *
         * :param const floatVector &valuesAtPoints: A vector of the values at the data points. Stored as
         *     [ v_11, v_12, ..., v_21, v22, ... ] where the first index is the point index in order as 
         *     provided to the volume reconstruction object and the second index is the value of the 
         *     function to be integrated.
         * :param const uIntType valueSize: The size of the subvector associated with each of the datapoints.
         * :param const floatVector &origin: The origin.
         * :param floatVector &integratedValue: The final value of the integral
         * :param bool computeFlux: The flag indicating if the flux needs to be computed of the values at the surface.
         * :param bool positionWeightedIntegral: Flag indicated if the integral is weighted by the position
         * :param bool dyadWithOrigin: The flag indicating if the dyadic product between the values on the surface
         *     ( post flux calculation ) needs to be computed.
         * :param uIntVector *subdomainIDs: The IDs of points in the subdomain to integrate over
         * :param const floatVector *subdomainWeights: The weights for the subdomains. Useful if points can be
         *     in multiple subdomains and they aren't small w.r.t. the domain size
         * :param const floatVector *macroNormal: A macro-scale normal vector to use to generate the micro
         *     weight. This can be helpful in cases where some points start to, ``wrap,'' around an edge 
         *     which should be flat. Can either be a single vector of dimension _dim or a collection of vectors
         *     at each boundary point.
         * :param const bool useMacroNormal: Use the macro-scale normal instead of the micro normals. Can help
         *     drive the integral to be more what is expected in some cases.
         */

        errorOut error;

        if ( ( subdomainIDs ) && ( subdomainWeights ) ){

            if ( subdomainIDs->size( ) != subdomainWeights->size( ) ){

                return new errorNode( "performSurfaceIntegration",
                                      "The size of the subdomain ids and subdomain weights are not consistent" );

            }

        }

        if ( ( !subdomainIDs ) && ( subdomainWeights ) ){

            return new errorNode( "performSurfaceIntegration",
                                  "The subdomain weights are defined but not the subdomain" );

        }

        if ( ( macroNormal ) && ( subdomainWeights ) ){
            return new errorNode( "performSurfaceIntegration",
                                  "Both the macro normal and subdomain weights can't be provided." );
        }

        if ( ( macroNormal ) && ( subdomainIDs ) ){

            if ( ( macroNormal->size( ) != ( subdomainIDs->size( ) * _dim ) ) &&
                 ( macroNormal->size( ) != _dim ) ){
    
                return new errorNode( "performSurfaceIntegration",
                                      "The macro normal and subdomainIDs vector are not of consistent sizes. It must\n either be of length " +
                                      std::to_string( _dim ) + " or " + std::to_string( _dim ) +
                                      " times the number of subdomain IDs" );
    
            }

        }

        if ( ( macroNormal ) && ( !subdomainIDs ) ){
            
            return new errorNode( "performSurfaceIntegration",
                                  "The macro normal and subdomainIDs vector must both be defined together" );

        }

        if ( ( !macroNormal ) && ( useMacroNormal ) ){

            return new errorNode( "performSurfaceIntegration",
                                  "The macro normal is requested to be used for flux calculations but it is not defined" );

        }

        //Check if the domain has been constructed yet
        if ( !getEvaluated( ) ){
            error = evaluate( );

            if ( error ){

                errorOut result = new errorNode( "performSurfaceIntegration",
                                                 "Error encountered during the reconstruction of the volume" );
                result->addNext( error );
                return result;

            }
        }

        //Interpolate the function to the background grid
        std::unordered_map< uIntType, floatVector > functionAtGrid;
        error = interpolateFunctionToBackgroundGrid( valuesAtPoints, valueSize, functionAtGrid );

        if ( error ){

            errorOut result = new errorNode( "performSurfaceIntegration",
                                             "Error encountered during the interpolation of the function to the background grid" );
            result->addNext( error );
            return result;

        }

        //Perform the surface integration
        floatVector boundaryPoint;
        floatVector localBoundaryPoint;
        floatVector valueAtBoundaryPoint;
        floatType weight;

        uIntType ngy = _gridLocations[ 1 ].size( );
        uIntType ngz = _gridLocations[ 2 ].size( );
        uIntType i, j, k;
        uIntVector indices;

        std::unique_ptr< elib::Element > element;

        if ( computeFlux ){

            integratedValue = floatVector( ( _dim * positionWeightedIntegral + !positionWeightedIntegral ) * valueSize / _dim, 0 );

        }
        else{

            integratedValue = floatVector( ( _dim * positionWeightedIntegral + !positionWeightedIntegral ) * valueSize, 0 );

        }

        if ( dyadWithOrigin ){

            if ( origin.size( ) != _dim ){

                return new errorNode( "performSurfaceIntegration",
                                      "The origin must be of dimension: " + std::to_string( _dim ) );

            }

            integratedValue = floatVector( integratedValue.size( ) * _dim, 0 );

        }

        floatMatrix nodalFunctionValues;

        const uIntVector *surfaceCells = &_boundaryCells;

        if ( subdomainIDs ){

            surfaceCells = subdomainIDs;

        }

        for ( auto cell = surfaceCells->begin( ); cell != surfaceCells->end( ); cell++ ){

            //Get the bottom corner node IDs from the cell id
            i = *cell / ( ngy * ngz );
            j = ( *cell - ngy * ngz * i ) / ngz;
            k = (*cell - ngy * ngz * i - ngz * j );

            indices = { i, j, k };

            //Get the element associated with this cell
            error = getGridElement( indices, element );

            if ( error ){

                errorOut result = new errorNode( "performSurfaceIntegration",
                                                 "Error in getting the grid element" );
                result->addNext( error );
                return result;

            }

            //Get the boundary point associated with this cell
            auto index = _boundaryPointIDToIndex.find( *cell );
            if ( index == _boundaryPointIDToIndex.end( ) ){

                return new errorNode( "performSurfaceIntegration",
                                      "The boundary cell is not found in the boundary point ID to index map" );

            }

            boundaryPoint = floatVector( _boundaryPoints.begin( ) + _dim * ( index->second ),
                                         _boundaryPoints.begin( ) + _dim * ( index->second + 1 ) );

            //Compute the local coordinates of the point
            error = element->compute_local_coordinates( boundaryPoint, localBoundaryPoint );

            if ( error ){

                errorOut result = new errorNode( "performSurfaceIntegration",
                                                 "Error in the computation of the local coordinates" );
                result->addNext( error );
                return result;

            }

            //Extract the function at the nodes
            nodalFunctionValues = floatMatrix( element->global_node_ids.size( ), floatVector( valueSize ) );
            for ( auto nID = element->global_node_ids.begin( ); nID != element->global_node_ids.end( ); nID++ ){

                if ( functionAtGrid.find( *nID ) == functionAtGrid.end( ) ){

                    return new errorNode( "performSurfaceIntegration",
                                          "Node with global ID " + std::to_string( *nID )
                                        + " not found in the grid node to function map" );

                }

                nodalFunctionValues[ nID - element->global_node_ids.begin( ) ]
                    = ( floatType )( _implicitFunctionValues[ *nID ] > 0 ) * functionAtGrid[ *nID ];

            }

            //Process the nodes which are outside of the volume

            //Get whether the element's nodes are inside of the reconstructed volume
            floatVector nodalImplicitFunction( element->reference_nodes.size( ), 0 );
            for ( auto nID = element->global_node_ids.begin( ); nID != element->global_node_ids.end( ); nID++ ){

                 nodalImplicitFunction[ nID - element->global_node_ids.begin( ) ]
                     = ( floatType )( _implicitFunctionValues[ *nID ] > 0 );

            }

            //Interpolate the function to the boundary point
            floatVector boundaryPointShapeFunctions;
            element->get_shape_functions( localBoundaryPoint, boundaryPointShapeFunctions );
            floatType normalizationFactor = 0.;
            valueAtBoundaryPoint = floatVector( nodalFunctionValues[ 0 ].size( ), 0 );
            for ( auto N = boundaryPointShapeFunctions.begin( ); N != boundaryPointShapeFunctions.end( ); N++ ){
                unsigned int eIndex = N - boundaryPointShapeFunctions.begin( );
                valueAtBoundaryPoint += nodalFunctionValues[ eIndex ] * ( *N ) * nodalImplicitFunction[ eIndex ];
                normalizationFactor += ( *N ) * nodalImplicitFunction[ eIndex ];
            }
            valueAtBoundaryPoint /= normalizationFactor;

            if ( error ){

                errorOut result = new errorNode( "performSurfaceIntegration",
                                                 "Error in interpolation of the nodal function to the boundary point" );
                result->addNext( error );
                return result;

            }

            //Compute the flux of the value if required
            if ( computeFlux ){

                if ( useMacroNormal ){

                    uIntType tmpIndex = cell - surfaceCells->begin( );
                    if ( macroNormal->size( ) == _dim ){

                        valueAtBoundaryPoint = vectorTools::matrixMultiply( *macroNormal, valueAtBoundaryPoint,
                                                                            1, _dim, _dim, valueSize / _dim );

                    }
                    else{

                        valueAtBoundaryPoint = vectorTools::matrixMultiply( floatVector( macroNormal->begin( ) + _dim * tmpIndex,
                                                                                         macroNormal->begin( ) + _dim * ( tmpIndex + 1 ) ),
                                                                            valueAtBoundaryPoint, 1, _dim, _dim, valueSize / _dim );

                    }

                }
                else{

                    valueAtBoundaryPoint = vectorTools::matrixMultiply( _boundaryPointNormals[ *cell ], valueAtBoundaryPoint,
                                                                        1, _dim, _dim, valueSize / _dim );

                }

            }

            if ( dyadWithOrigin ){

                valueAtBoundaryPoint = vectorTools::appendVectors( vectorTools::dyadic( valueAtBoundaryPoint, boundaryPoint - origin ) );

            }

            if ( _boundaryPointAreas.find( *cell ) == _boundaryPointAreas.end( ) ){

                return new errorNode( "performSurfaceIntegration",
                                      "The current boundary point is not found in the boundary point areas map" );

            }

            weight = 1.;
            if ( subdomainWeights ){

                weight = ( *subdomainWeights )[ cell - surfaceCells->begin( ) ];

            }
            else if ( macroNormal ){

                uIntType tmpIndex = cell - surfaceCells->begin( );
                if ( macroNormal->size( ) == _dim ){

                    weight = std::fabs( vectorTools::dot( *macroNormal, _boundaryPointNormals[ *cell ] ) );

                }
                else{

                    weight = std::fabs( vectorTools::dot( floatVector( macroNormal->begin( ) + _dim * tmpIndex,
                                                                       macroNormal->begin( ) + _dim * ( tmpIndex + 1 ) ),
                                                          _boundaryPointNormals[ *cell ] ) );

                }

            }

            if ( positionWeightedIntegral ){

                floatVector positionWeightedValueAtBoundaryPoint( valueAtBoundaryPoint.size( ) * boundaryPoint.size( ), 0 );

                uIntType bPi_index = 0;

                for ( auto bPi = boundaryPoint.begin( ); bPi != boundaryPoint.end( ); bPi++, bPi_index++ ){

                    uIntType vBI_index = 0;

                    for ( auto vBI = valueAtBoundaryPoint.begin( ); vBI != valueAtBoundaryPoint.end( ); vBI++, vBI_index++ ){

                        positionWeightedValueAtBoundaryPoint[ vBI_index * boundaryPoint.size( ) + bPi_index ] = ( *vBI ) * ( *bPi );

                    }

                }

                integratedValue += weight * positionWeightedValueAtBoundaryPoint * _boundaryPointAreas[ *cell ];

            }
            else{

                integratedValue += weight * valueAtBoundaryPoint * _boundaryPointAreas[ *cell ];

            }

        }

        return NULL;
    }

    errorOut dualContouring::getSurfaceSubdomains( const floatType &minDistance, uIntVector &subdomainNodeCounts,
                                                   uIntVector &subdomainIDs ){
        /*!
         * Break the surface into subdomains which are separated by ( approximately ) minDistance.
         *
         * :param const floatType &minDistance: The minimum distance between the subdomain centers
         * :param uIntVector &subdomainNodeCounts: The number of nodes in each subdomain
         * :param uIntVector &subdomainIDs: The IDs of the nodes in the subdomains in the order of subdomainNodeCounts
         */

        errorOut error;

        //Check if the domain has been constructed yet
        if ( !getEvaluated( ) ){
            error = evaluate( );

            if ( error ){

                errorOut result = new errorNode( "performSurfaceIntegration",
                                                 "Error encountered during the reconstruction of the volume" );
                result->addNext( error );
                return result;

            }
        }

        //Initialize the subdomain nodes
        subdomainIDs.reserve( _boundaryCells.size( ) );

        if ( _boundaryCells.size( ) < 1 ){

            return new errorNode( "getSubsurfaceDomains",
                                  "Boundary points must contain at least one node" );

        }

        /*=====================================================================
        |                       Identify the seed nodes                       |
        =====================================================================*/

        uIntVector remainingNodes = _boundaryCells; //Copy over the boundary cells

        //Begin seed node loop
        uIntVector seedNodeIDs;
        while ( remainingNodes.size( ) > 0 ){ 

            //Store a new subdomain seed node
            seedNodeIDs.push_back( remainingNodes[ 0 ] );
            floatVector currentSeedPoint( _boundaryPoints.begin( ) + _dim * _boundaryPointIDToIndex[ seedNodeIDs.back( ) ],
                                          _boundaryPoints.begin( ) + _dim * ( _boundaryPointIDToIndex[ seedNodeIDs.back( ) ] + 1 ) );

            //Form a KD tree of the remaining nodes
            floatVector remainingNodeCoords;
            remainingNodeCoords.reserve( _dim * remainingNodes.size( ) );
            uIntVector ownedIndices;
            ownedIndices.reserve( remainingNodes.size( ) );

            for ( auto rN = remainingNodes.begin( ); rN != remainingNodes.end( ); rN++ ){

                for ( unsigned int i = 0; i < _dim; i++ ){

                    remainingNodeCoords.push_back( _boundaryPoints[ _dim * _boundaryPointIDToIndex[ *rN ] + i ] );

                }

                ownedIndices.push_back( _dim * ( rN - remainingNodes.begin( ) ) );

            }

            //TODO: Change this to delete the nodes from the tree rather than rebuilding the tree
            KDNode remainingTree( &remainingNodeCoords, ownedIndices, 0, _dim );

            //Find the nodes that are within the min-distance radius of the most recently added seed point
            uIntVector internalNodes;
            remainingTree.getPointsWithinRadiusOfOrigin( currentSeedPoint, minDistance, internalNodes );

            //TODO: This is just the Euclidean distance. If there are multiple surfaces this approximation will fail
            //      This does help to winnow the cells down for a true Geometric distance to be computed using Dijkstra's
            //      algorithm or similar to eliminate nodes that cannot be reached from the seed node

            //Remove those nodes from the remaining node coords
            internalNodes /= _dim;
            std::sort( internalNodes.begin( ), internalNodes.end( ) );

            for ( auto iN = internalNodes.rbegin( ); iN != internalNodes.rend( ); iN++ ){

                if ( *iN < ( remainingNodes.size( ) - 1 ) ){
                    //Swap the element to the end
                    std::swap( remainingNodes[ *iN ], remainingNodes.back( ) );
                }
                remainingNodes.pop_back( );

            }

        }

        /*=====================================================================
        |   Associate the boundary points in the domain with the seed nodes   |
        =====================================================================*/
        
        uIntMatrix seedNodePoints( seedNodeIDs.size( ) );

        for ( auto sNP = seedNodePoints.begin( ); sNP != seedNodePoints.end( ); sNP++ ){

            sNP->reserve( _boundaryCells.size( ) / seedNodePoints.size( ) ); //Assume the points will be distributed evenly

        }

        //Loop through the boundary cells;
        for ( auto bC = _boundaryCells.begin( ); bC!= _boundaryCells.end( ); bC++ ){

            //Get the coordinates of the current boundary cell
            floatVector currentBoundaryPoint( _boundaryPoints.begin( ) + _dim * _boundaryPointIDToIndex[ *bC ],
                                              _boundaryPoints.begin( ) + _dim * ( _boundaryPointIDToIndex[ *bC ] + 1 ) );

            //Get the current seed point
            floatVector currentSeedPoint( _boundaryPoints.begin( ) + _dim * _boundaryPointIDToIndex[ seedNodeIDs[ 0 ] ],
                                          _boundaryPoints.begin( ) + _dim * ( _boundaryPointIDToIndex[ seedNodeIDs[ 0 ] ] + 1 ) );

            //Compute the distance
            floatType bestDistance = vectorTools::l2norm( currentBoundaryPoint - currentSeedPoint );
            uIntType seedNum = 0;

            for ( auto sNP = seedNodeIDs.begin( ) + 1; sNP != seedNodeIDs.end( ); sNP++ ){

                //Get the coordinates of the current seed node
                floatVector currentSeedPoint( _boundaryPoints.begin( ) + _dim * _boundaryPointIDToIndex[ *sNP ],
                                              _boundaryPoints.begin( ) + _dim * ( _boundaryPointIDToIndex[ *sNP ] + 1 ) );

                //Compute the distance
                floatType currentDistance = vectorTools::l2norm( currentBoundaryPoint - currentSeedPoint );

                //Check if the current seed point is the closest
                if ( bestDistance > currentDistance ){

                    //Set the current point as the closest
                    bestDistance = currentDistance;
                    seedNum = sNP - seedNodeIDs.begin( );

                }

            }

            //Add the boundary point to the nearest seed node point
            seedNodePoints[ seedNum ].push_back( *bC );

        }

        //Update the counts for the number of nodes in each subdomain
        subdomainNodeCounts.clear( );
        subdomainNodeCounts.resize( seedNodeIDs.size( ) );

        unsigned int indx = 0;
        for ( auto sNP = seedNodePoints.begin( ); sNP != seedNodePoints.end( ); sNP++, indx++ ){

            subdomainNodeCounts[ indx ] = sNP->size( );

        }

        subdomainIDs = vectorTools::appendVectors( seedNodePoints );

        return NULL;
    }

    const YAML::Node volumeReconstructionBase::exportConfiguration( ){
        /*!
         * Export the configuration
         */

        return YAML::Clone( _config );

    }

    const uIntVector *dualContouring::getBoundaryIDs( ){
        /*!
         * Get a constant reference to the boundary cells
         */

        return &_boundaryCells;

    }

    const floatVector *dualContouring::getBoundaryPoints( ){
        /*!
         * Get a constant reference to the boundary points
         */

        return &_boundaryPoints;

    }

}
