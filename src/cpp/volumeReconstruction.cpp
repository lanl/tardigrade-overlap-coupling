/*=============================================================================
|                             volumeReconstruction                            |
===============================================================================
| Tools to reconstruct volume information from pointsets                      |
=============================================================================*/

#include<volumeReconstruction.h>
#include<solver_tools.h>

namespace volumeReconstruction{

    KDNode::KDNode( ){
        /*!
         * The default constructor for the KD tree
         */

        return;
    }

    KDNode::KDNode( const floatVector *points, const uIntVector &ownedIndices,
                    const unsigned int &depth, const unsigned int &dim ){

        /*!
         * The constructor for the KD tree
         *
         * :param const floatVector *points: The points to be sorted into the tree
         * :param const uIntVector &ownedIndices: The indices of the points which this
         *     node is dividing. These are indices to the first values of the 
         *     points in the 'points' vector.
         * :param const unsigned int &depth: The depth of the current node
         * :param const unsigned int &dim: The dimension of the points
         */

        //Set the points
        _points = points;

        //Set the depth of the node
        _depth = depth;

        //Determine if this is the last node in the chain
        if ( ownedIndices.size( ) == 1 ){
            _index = ownedIndices[ 0 ];
            return;
        }

        //Determine the bounding box of the points
        floatVector lowerBound( _points->begin( ) + ownedIndices[ 0 ],
                                _points->begin( ) + ownedIndices[ 0 ] + dim );

        floatVector upperBound( _points->begin( ) + ownedIndices[ 0 ],
                                _points->begin( ) + ownedIndices[ 0 ] + dim );

        for ( auto index = ownedIndices.begin( ) + 1; index != ownedIndices.end( ); index++ ){

            //Compare the current point against the bounds
            for ( unsigned int i = 0; i < dim; i++ ){
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

        for ( auto v = delta.begin( ) + 1; v != delta.end( ); v++ ){

            if ( deltaMax < *v ){

                _axis = v - delta.begin( );
                deltaMax = *v;

            }            

        }

        //Accumulate the dimension indicated by axis
        typedef std::pair< const unsigned int*, const floatType* > valueMap;
        
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

    const unsigned int* KDNode::getIndex( ){
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
        unsigned int dim = upperBounds.size( );

        floatVector _domainUpperBounds;
        floatVector _domainLowerBounds;

        //Assemble the current point
        floatVector median( _points->begin( ) + _index,
                            _points->begin( ) + _index + dim );

        floatVector upperDelta, lowerDelta;

        if( !domainUpperBounds ){
            
            _domainUpperBounds = floatVector( dim );
            _domainLowerBounds = floatVector( dim );

            for ( unsigned int i = 0; i < dim; i++ ){

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

    floatType KDNode::getMinimumValueDimension( const unsigned int &d ){
        /*!
         * Get the minimum value of a given dimension in the tree
         *
         * :param const unsigned int &d: The dimension ( starting at zero ) of each point to search
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

    floatType KDNode::getMaximumValueDimension( const unsigned int &d ){
        /*!
         * Get the maximum value of a given dimension in the tree
         *
         * :param const unsigned int &d: The dimension ( starting at zero ) of each point to search
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

    void KDNode::printData( const unsigned int &dim ){
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
            if ( !_config[ "write_config" ].IsScalar( ) ){
                _config[ "write_config" ] = "defaultOutput.yaml";
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
        for ( unsigned int i = 0; i < _dim * _nPoints; i+=3 ){
            ownedIndices.push_back( i );
        }

        //Form the KD tree
        _tree = KDNode( _points, ownedIndices, 0, _dim );

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
         * Base function which evaluates the volume reconstruction.
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

        return NULL;
    }

    errorOut volumeReconstructionBase::setInterpolationConfiguration( ){
        /*!
         * Set the interpolation configuration
         */

        //Set the default interpolation configuration
        if ( !_config[ "interpolation" ] ){
            _config[ "interpolation" ][ "type" ] = "constant";
            _config[ "interpolation" ][ "constant_value" ] = 1;
            _functionValue = 1;
        }

        if ( !_config[ "interpolation" ][ "type" ] ){
            _config[ "interpolation" ][ "type" ] = "constant";
            _config[ "interpolation" ][ "constant_value" ] = 1;
            _functionValue = 1;
        }

        if ( ( _config[ "interpolation" ][ "type" ].as< std::string >( ).compare( "constant" ) == 0 ) &&
             ( !_config[ "interpolation" ][ "constant_value" ] ) ){

            _config[ "interpolation" ][ "constant_value" ] = 1;
            _functionValue = 1;
            
        }
        
        if ( ( _config[ "interpolation" ][ "type" ].as< std::string >( ).compare( "constant" ) == 0 ) &&
             ( _config[ "interpolation" ][ "constant_value" ] ) ){

            _functionValue = _config[ "interpolation" ][ "constant_value" ].as< floatType >( );
            
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

        for ( unsigned int i = 0; i < _dim; i++ ){

            _upperBounds[ i ] = _tree.getMaximumValueDimension( i );
            _lowerBounds[ i ] = _tree.getMinimumValueDimension( i );

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

    errorOut volumeReconstructionBase::getFunctionValue( const unsigned int i, floatType &value ){
        /*!
         * Get the function value at the provided index. This is the index of the function-value vector
         * corresponding to a given point not the index in the points vector.
         *
         * :param const unsigned int i: The index of the point
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

        return NULL;
    }

    errorOut dualContouring::processConfigurationFile( ){
        /*!
         * Process the configuration file reading options and setting defaults.
         */

        if ( !_config[ "interpolation" ][ "discretization_count" ] ){

            _config[ "interpolation" ][ "discretization_count" ] = std::max( ( unsigned int )( std::pow( ( floatType )_nPoints, 1. / 3. ) ),
                                                                                                         ( unsigned int )1 );

        }

        if ( _config[ "interpolation" ][ "discretization_count" ].IsScalar( ) ){

            unsigned int v = _config[ "interpolation" ][ "discretization_count" ].as< unsigned int >( );
            _domainDiscretization = { v, v, v };

        }
        else if ( _config[ "interpolation" ][ "discretization_count" ].IsSequence( ) ){

            _domainDiscretization.resize( _config[ "interpolation" ][ "discretization_count" ].size( ) );

            if ( _domainDiscretization.size( ) != _dim ){

                return new errorNode( "processConfigFile",
                                      "The number of discretization indices ( " + std::to_string( _domainDiscretization.size( ) ) 
                                      + " ) is not the same as the dimension ( " + std::to_string( _dim ) + " )" );

            }

            unsigned int i = 0;
            for ( auto it  = _config[ "interpolation" ][ "discretization_count" ].begin( );
                       it != _config[ "interpolation" ][ "discretization_count" ].end( );
                       it++ ){

                _domainDiscretization[ i ] = it->as< unsigned int>( );

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

        for ( unsigned int i = 0; i < _dim; i++ ){

            //With the domain being split into n discretizations it requires n + 1 nodes
            //We also add two additional nodes on the outside of the grid to provide
            //a hard boundary. These nodes are located very close to the outside
            _gridLocations[ i ].resize( _domainDiscretization[ i ] + 3 );

            //Set the displacement of the nodes
            delta = ( ( *upperBounds )[ i ] - ( *lowerBounds )[ i ] ) / _domainDiscretization[ i ];

            for ( unsigned int j = 0; j < _domainDiscretization[ i ] + 1; j++ ){

                _gridLocations[ i ][ j + 1 ] = ( *lowerBounds )[ i ] + j * delta;

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
        _implicitFunctionValues.resize( _gridLocations[ 0 ].size( ) *
                                        _gridLocations[ 1 ].size( ) *
                                        _gridLocations[ 2 ].size( ) );

        uIntVector gridPointCounts( _gridLocations[ 0 ].size( ) *
                                    _gridLocations[ 1 ].size( ) *
                                    _gridLocations[ 2 ].size( ) );

        //Loop over the elements

        errorOut error;
        uIntVector elementIndices;
        floatVector elementNodalContributions;
        uIntVector globalNodeIds;
        uIntVector pointCounts;

        unsigned int ngx = _gridLocations[ 0 ].size( );
        unsigned int ngy = _gridLocations[ 1 ].size( );
        unsigned int ngz = _gridLocations[ 2 ].size( );

        for ( unsigned int i = 1; i < ngx - 2; i++ ){

            for ( unsigned int j = 1; j < ngy - 2; j++ ){

                for ( unsigned int k = 1; k < ngz - 2; k++ ){

                    //Get the element contribution to the nodal values of the implicit function
                    elementIndices = { i, j, k };
                    error = processBackgroundGridElementImplicitFunction( elementIndices, elementNodalContributions,
                                                                          globalNodeIds, pointCounts );

                    if ( error ){

                        errorOut result = new errorNode( "projectImplicitFunctionToBackgroundGrid",
                                                         "Error in processing the projection of the implicit function to the nodes for the element with the lower cornrer of indices i, j, k: "
                                                       + std::to_string( i ) + ", "
                                                       + std::to_string( j ) + ", "
                                                       + std::to_string( k ) );

                        result->addNext( error );
                        return result;

                    }

                    for ( unsigned int i = 0; i < globalNodeIds.size( ); i++ ){

                        //Add those values to the grid
                        _implicitFunctionValues[ globalNodeIds[ i ] ] += elementNodalContributions[ i ];
    
                        //Add the point counts
                        gridPointCounts[ globalNodeIds[ i ] ] += pointCounts[ i ];

                    }

                }

            }

        }

        for ( unsigned int i = 1; i < ngx - 1; i++ ){

            for ( unsigned int j = 1; j < ngy - 1; j++ ){

                for ( unsigned int k = 1; k < ngz - 1; k++ ){

                    if ( gridPointCounts[ ngy * ngz * i + ngz * j + k ] > 0 ){

                        _implicitFunctionValues[ ngy * ngz * i + ngz * j + k ] /=
                            gridPointCounts[ ngy * ngz * i + ngz * j + k ];

                    }

                }

            }

        }

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

        unsigned int ngy = _gridLocations[ 1 ].size( );
        unsigned int ngz = _gridLocations[ 2 ].size( );

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
        floatMatrix nodes = { { lbCoordinates[ 0 ], lbCoordinates[ 1 ], lbCoordinates[ 2 ] },
                              { ubCoordinates[ 0 ], lbCoordinates[ 1 ], lbCoordinates[ 2 ] },
                              { ubCoordinates[ 0 ], ubCoordinates[ 1 ], lbCoordinates[ 2 ] },
                              { lbCoordinates[ 0 ], ubCoordinates[ 1 ], lbCoordinates[ 2 ] },
                              { lbCoordinates[ 0 ], lbCoordinates[ 1 ], ubCoordinates[ 2 ] },
                              { ubCoordinates[ 0 ], lbCoordinates[ 1 ], ubCoordinates[ 2 ] },
                              { ubCoordinates[ 0 ], ubCoordinates[ 1 ], ubCoordinates[ 2 ] },
                              { lbCoordinates[ 0 ], ubCoordinates[ 1 ], ubCoordinates[ 2 ] } };

        //Set the element's global node numbers
        uIntVector global_node_ids=
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

        _tree.getPointsInRange( element->bounding_box[ 1 ], element->bounding_box[ 0 ], pointIndices,
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
            for ( auto d = distances.begin( ); d != distances.end( ); d++ ){

                nodesSupported[ d - distances.begin( ) ] = ( floatType )( vectorTools::fuzzyEquals( *d, minDistance ) );
                pointCounts[ d - distances.begin( ) ] += nodesSupported[ d - distances.begin( ) ];

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

        unsigned int ngx = _gridLocations[ 0 ].size( );
        unsigned int ngy = _gridLocations[ 1 ].size( );
        unsigned int ngz = _gridLocations[ 2 ].size( );

        //Resize the internal cells vector
        _internalCells.clear( );
        _boundaryCells.clear( );

        _internalCells.reserve( ( ngx - 1 ) * ( ngy - 1 ) * ( ngz - 1 ) );
        _boundaryCells.reserve( ( ngx - 1 ) * ( ngy - 1 ) * ( ngz - 1 ) );

        floatVector cellValues;

        for ( unsigned int i = 0; i < ( ngx - 1 ); i++ ){

            for ( unsigned int j = 0; j < ( ngy - 1 ); j++ ){

                for ( unsigned int k = 0; k < ( ngz - 1 ); k++ ){

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

                            //The cell is on the surface of the body
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

        unsigned int ngx = _gridLocations[ 0 ].size( );
        unsigned int ngy = _gridLocations[ 1 ].size( );
        unsigned int ngz = _gridLocations[ 2 ].size( );

        //Resize the boundary point vector
        _boundaryPoints.clear( );

        _boundaryPoints.reserve( _dim * _boundaryCells.size( ) );
        std::cout << "_boundaryPoints.size( ): " << _dim * _boundaryCells.size( ) << "\n";

        //Loop over the boundary cells
        unsigned int i, j, k;
        unsigned int ri1, rj1, rk1;
        unsigned int ri2, rj2, rk2;
        unsigned int ri, rj, rk;
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

        floatVector boundaryPoint;

        unsigned int edgeID;
        uIntVector edgeCells;

        for ( auto bc = _boundaryCells.begin( ); bc != _boundaryCells.end( ); bc++ ){

            //Determine the lower-left hand corner index
            i = *bc / ( ngy * ngz );
            j = ( *bc - ( ngy * ngz * i ) ) / ngz;
            k = *bc - ngy * ngz * i - ngz * j;

            //Determine the grid element
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

                unsigned int i2 = edgeNodes[ 2 * ( eT - edgeTransition.begin( ) ) + 1 ];
                unsigned int i1 = edgeNodes[ 2 * ( eT - edgeTransition.begin( ) ) + 0 ];

                for ( unsigned int i = 0; i < _dim; i++ ){

                    if ( ( element->reference_nodes[ i2 ][ i ] - element->reference_nodes[ i1 ][ i ] ) < _absoluteTolerance ){

                        intersectionPoint[ i ] = element->reference_nodes[ i2 ][ i ];

                    }
                    else{

                        m = ( cellValues[ i2 ] - cellValues[ i1 ] ) / ( element->reference_nodes[ i2 ][ i ] - element->reference_nodes[ i1 ][ i ] );
                        b = cellValues[ i1 ] - m * element->reference_nodes[ i1 ][ i ];

                        intersectionPoint[ i ] = -b / m;

                    }

                }

                points.push_back( intersectionPoint );

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

                }
                else{

                    rootNode = element->reference_nodes[ i1 ];
                    ri = ri1;
                    rj = rj1;
                    rk = rk1;

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
                _tree.getPointsInRange( upperBounds, lowerBounds, supportingPoints, &domainUpperBounds, &domainLowerBounds );

                //Determine the normal at the intersection point
                if ( ( supportingPoints.size( ) >= _minApproximationCount ) && ( _useMaterialPointsForNormals ) ){

                    return new errorNode( "computeBoundaryPoints", "Using the material points for normals has not been implemented yet" );

                }
                else{

                    //Use the gradient of the shape function to compute the normal
                   
                    //Compute the local coordinates
                    error = element->compute_local_coordinates( intersectionPoint, localIntersectionPoint );

                    if ( error ){

                        errorOut result = new errorNode( "computeBoundaryPoints", "Error in computation of the local coordinates of the intersection point" );
                        result->addNext( error );
                        return result;

                    }

                    //Compute the global gradient
                    error = element->get_global_gradient( cellValues, localIntersectionPoint, gradient );

                    //Compute the normal vector
                    gradient /= vectorTools::l2norm( gradient );

                    normals.push_back( gradient );

                }

                //Store the transition edge
                unsigned int edgeIndex = eT - edgeTransition.begin( );
                edgeID = ngy * ngz * ri1 + ngz * rj1 + rk1;

                if ( edgeIndex < 4 ){ //x edge

                    edgeCells =
                        {
                            ngy * ngz * ri1 + ngz * ( rj1 - 0 ) + ( rk1 - 0 ),
                            ngy * ngz * ri1 + ngz * ( rj1 - 1 ) + ( rk1 - 0 ),
                            ngy * ngz * ri1 + ngz * ( rj1 - 1 ) + ( rk1 - 1 ),
                            ngy * ngz * ri1 + ngz * ( rj1 - 0 ) + ( rk1 - 1 )
                        };

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

                    if ( _boundaryEdges_y.find( edgeID ) == _boundaryEdges_y.end( ) ){

                            _boundaryEdges_y.emplace( edgeID, edgeCells );

                    }

                }
                else{ //z edge

                    edgeCells =
                        {
                            ngy * ngz * ( ri1 - 0 ) + ngz * ( rj1 - 0 ) + rk1,
                            ngy * ngz * ( ri1 - 1 ) + ngz * ( rj1 - 0 ) + rk1,
                            ngy * ngz * ( ri1 - 1 ) + ngz * ( rj1 - 1 ) + rk1,
                            ngy * ngz * ( ri1 - 0 ) + ngz * ( rj1 - 1 ) + rk1
                        };

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
            floatArgs[ 0 ] = element->bounding_box[ 1 ];
            floatArgs[ 1 ] = element->bounding_box[ 0 ];

            for ( unsigned int i = 0; i < points.size( ); i++ ){

                floatArgs[ 2 + i ] = points[ i ];
                floatArgs[ 2 + points.size( ) + i ] = normals[ i ];

            }

            //Assemble the integer arguments matrix
            intArgs = { { ( int )_dim, ( int )points.size( ) } }; 

            //Assemble the initial iterate
            X0.resize( 5 * _dim );
            X.resize( 5 * _dim );
            for ( unsigned int i = 0; i < _dim; i++ ){

                X0[            i ] = 0.5 * ( element->bounding_box[ 0 ][ i ] + element->bounding_box[ 1 ][ i ] ); //x
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

                boundaryPoint = floatVector( X0.begin( ), X0.begin( ) + _dim );

            }
            else{
                
                boundaryPoint = floatVector( X.begin( ), X.begin( ) + _dim );

            }

            for ( unsigned int i = 0; i < _dim; i++ ){
                _boundaryPoints.push_back( boundaryPoint[ i ] );
            }

            //Store the transition edge

        }

        std::cout << "boundary points\n";
        std::cout << "{ " << _boundaryCells[ 0 ] << ": [";
        for ( unsigned int i = 0;  i < _boundaryPoints.size( ); i++ ){

            std::cout << _boundaryPoints[ i ] << ", ";

            if ( ( ( i + 1 ) % 3 ) == 0 ){
                std::cout << "]\n";

                if ( ( i + 1 ) / _dim < _boundaryCells.size( ) ){
                    std::cout << _boundaryCells[ ( i + 1 ) / _dim ] << ": [";
                }
            }

        }
        std::cout << "}\n";

        std::cout << "boundaryEdges_x\n";
        std::cout << "{\n";
        for ( auto it = _boundaryEdges_x.begin( ); it != _boundaryEdges_x.end( ); it++ ){

            std::cout << it->first << ": " << "[ ";
            for ( auto it2 = it->second.begin( ); it2 != it->second.end( ); it2++ ){

                std::cout << *it2 << ", ";

            }
            std::cout << "]\n";

        }
        std::cout << "}\n";

        std::cout << "boundaryEdges_y\n";
        std::cout << "{\n";
        for ( auto it = _boundaryEdges_y.begin( ); it != _boundaryEdges_y.end( ); it++ ){

            std::cout << it->first << ": " << "[ ";
            for ( auto it2 = it->second.begin( ); it2 != it->second.end( ); it2++ ){

                std::cout << *it2 << ", ";

            }
            std::cout << "]\n";

        }
        std::cout << "}\n";

        std::cout << "boundaryEdges_z\n";
        std::cout << "{\n";
        for ( auto it = _boundaryEdges_z.begin( ); it != _boundaryEdges_z.end( ); it++ ){

            std::cout << it->first << ": " << "[ ";
            for ( auto it2 = it->second.begin( ); it2 != it->second.end( ); it2++ ){

                std::cout << *it2 << ", ";

            }
            std::cout << "]\n";

        }
        std::cout << "}\n";

        assert( 1 == 0 );

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

        if ( intArgs.size( ) != 1 ){

            return new errorNode( "internalPointResidual", "The intArgs matrix must have one element" );

        }

        if ( intArgs[ 0 ].size( ) != 2 ){
                
            return new errorNode( "internalPointResidual", "The first value of intArgs must have a length of 2" );

        }


        unsigned int dim = intArgs[ 0 ][ 0 ];
        unsigned int nPoints = intArgs[ 0 ][ 1 ];

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

        for ( unsigned int i = 0; i < nPoints; i++ ){

            //Add the contribution to the first residual
            floatType nxmp = vectorTools::dot( normals[ i ], x - points[ i ] );

            for ( unsigned int _i = 0; _i < dim; _i++ ){

                residual[ _i ] += nxmp * normals[ i ][ _i ];

                for ( unsigned int _j = 0; _j < dim; _j++ ){

                    jacobian[ _i ][ _j ] += normals[ i ][ _i ] * normals[ i ][ _j ];

                }

            }

        }

        for ( unsigned int i = 0; i < dim; i++ ){

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

    errorOut interpolateFunctionToBackgroundGrid( );
}