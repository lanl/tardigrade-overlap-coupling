/*=============================================================================
|                             volumeReconstruction                            |
===============================================================================
| Tools to reconstruct volume information from pointsets                      |
=============================================================================*/

#include<volumeReconstruction.h>

namespace volumeReconstruction{

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
        _nPoints = _points->size( );

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

        setInterpolationConfiguration( );

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

    const floatVector *volumeReconstructionBase::getPoints( ){
        /*!
         * Get the points used in this class
         */

        return _points;
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

        return NULL;
    }
}
