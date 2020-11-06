/*=============================================================================
|                             generateXDMFData.cpp                            |
===============================================================================
| Functions used to generate XDMF data which can either be used to test the   |
| overlapCoupling class or to create filtering domains.                       |
=============================================================================*/

#include<generateXDMFData.h>

namespace fileGenerator{

    fileGenerator::fileGenerator( ){
        /*!
         * The default constructor
         */

        return;
    }

    fileGenerator::fileGenerator( const std::string &yamlFilename ){
        /*!
         * The constructor which will read in the YAML configuration file
         *
         * :param const std::string &yamlFilename: The YAML configuration file
         *
         */


        try{

            _config = YAML::LoadFile( yamlFilename );

        }
        catch( std::exception &e ){

            _error.reset( new errorNode( "fileGenerator", e.what( ) ) );

        }

        if ( !_config ){
            _error.reset( new errorNode( "fileGenerator", "File not found" ) );
        }

        //Initialize the output writer
        if ( !_config[ "output_configuration" ] ){
            std::cerr << "'output_configuration' not specified in the YAML file. Using defaults\n";
            _config[ "output_configuration" ][ "filename" ] = "xdmf_out";

        }

        //Remove any existing files
        std::string filename = _config[ "output_configuration" ][ "filename" ].as< std::string >( );
        std::string xmlfilename = filename;
        xmlfilename += ".xdmf";

        std::string h5filename = filename;
        h5filename += ".h5";

        remove( xmlfilename.c_str( ) );
        remove( h5filename.c_str( ) );

        //Force overwrite of 'mode' to be 'write'        
        _config[ "output_configuration" ][ "mode" ] = "write";

        //Force the filetype to be XDMF
        _config[ "output_configuration" ][ "filetype" ] = "XDMF";

        _writer = dataFileInterface::dataFileBase( _config[ "output_configuration" ] ).create( );

        if ( _writer->_error ){

            _error.reset( new errorNode( "fileGenerator", "Error when forming dataFileInterface writer" ) );
            _error->addNext( _writer->_error );
            return;

        }

        return;

    }

    const std::unique_ptr< errorNode > &fileGenerator::getError( ){
        /*!
         * Get a reference to the error object
         */

        return _error;
    }

    int fileGenerator::build( ){
        /*!
         * Build the XDMF file
         */

        if ( _config[ "increments" ] ){

            if ( _config[ "increments" ].IsSequence( ) ){

                uIntType increment_number = 0; 
                for ( auto increment = _config[ "increments" ].begin( ); increment != _config[ "increments" ].end( ); increment++, increment_number++ ){

                    //Initialize the increment
                    errorOut error = _initializeIncrement( *increment );
    
                    if ( error ){
    
                        _error.reset( new errorNode( "build", "Error in initialization of increment " + std::to_string( increment_number ) ) );
                        _error->addNext( error );
                        return 1;
    
                    }

                    //Write the mesh data
                    error = _writeMeshInformation( *increment );

                    if ( error ){

                        _error.reset( new errorNode( "build", "Error in writing the mesh information of increment " + std::to_string( increment_number ) ) );
                        _error->addNext( error );
                        return 1;

                    }

                    //Write the solution data
                    error = _writeSolutionInformation( *increment );

                    if ( error ){

                        _error.reset( new errorNode( "build", "Error in writing the solution information of increment " + std::to_string( increment_number ) ) );

                    }
    
                }
    
            }
            else{
    
                _error.reset( new errorNode( "build", "The keyword 'increments' must be a sequence of values at different increments" ) );
                return 1;
    
            }

        }
        else{

            _error.reset( new errorNode( "build", "The keyword 'increments' was not found in the configuration file" ) );
            return 1;

        }

        return 0;

    }

    errorOut fileGenerator::_initializeIncrement( const YAML::Node &increment ){
        /*!
         * Initialize an increment for output to the XDMF file
         *
         * :param const YAML::Node &increment: The increment information
         */

        //Get the reference increment
        uIntType reference_increment;
        errorOut error = _getPropertyFromYAML( increment, "reference_increment", reference_increment );
        if ( error ){

            errorOut result = new errorNode( "_initializeIncrement", "Error when extracting the reference increment" );
            result->addNext( error );
            return result;

        }

        //Get the time
        floatType time;
        error = _getPropertyFromYAML( increment, "time", time );
        if ( error ){

            errorOut result = new errorNode( "_initializeIncrement", "Error when extracting the time" );
            result->addNext( error );
            return result;

        }

        //Initialize the increment
        error = _writer->initializeIncrement( time, reference_increment, _collectionNumber, _currentIncrement );
        if ( error ){

            errorOut result = new errorNode( "_initializeIncrement", "Error in initialization of increment" );
            result->addNext( error );
            return result;

        }

        return NULL;
    }

    errorOut fileGenerator::_writeMeshInformation( const YAML::Node &increment ){
        /*!
         * Write the mesh information for the current increment to the XDMF file
         *
         * :param const YAML::Node &increment: The increment information
         */

        //Check the if the reference and current increments are the same
        uIntType reference_increment;
        errorOut error = _getPropertyFromYAML( increment, "reference_increment", reference_increment );
        if ( error ){

            errorOut result = new errorNode( "_writeMeshInformation", "Error when extracting the reference increment" );
            result->addNext( error );
            return result;

        }

        //Initialize variables
        uIntVector   nodeIds;
        uIntVector   elementIds;
        floatVector  nodePositions;
        stringVector nodeSetNames;
        stringVector elementSetNames;
        uIntMatrix   nodeSets;
        uIntMatrix   elementSets;
        uIntVector   connectivity;

        if ( reference_increment == _currentIncrement ){

            //Get the node id numbers
            error = _getPropertyFromYAML( increment, "node_ids", nodeIds );
            if ( error ){
    
                errorOut result = new errorNode( "_writeMeshInformation", "Error when extracting the node ids" );
                result->addNext( error );
                return result;
    
            }

            //Get the node positions
            error = _getPropertyFromYAML( increment, "node_positions", nodePositions );
            if ( error ){
    
                errorOut result = new errorNode( "_writeMeshInformation", "Error when extracting the node positions" );
                result->addNext( error );
                return result;
    
            }

            //Get the element id numbers
            error = _getPropertyFromYAML( increment, "element_ids", elementIds );
            if ( error ){
    
                errorOut result = new errorNode( "_writeMeshInformation", "Error when extracting the element ids" );
                result->addNext( error );
                return result;
    
            }

            //Get the connectivity vector
            error = _getPropertyFromYAML( increment, "connectivity", connectivity );
            if ( error ){

                errorOut result = new errorNode( "_writeMeshInformation", "Error when extracting the connectivity" );
                result->addNext( error );
                return result;

            }

            //Get the nodeset names and ids
            if ( increment[ "node_sets" ] ){
                error = _getKeyValuePairsFromYAML( increment, "node_sets", nodeSetNames, nodeSets ); 

                if ( error ){

                    errorOut result = new errorNode( "_writeMeshInformation", "Error when extracting the node sets" );
                    result->addNext( error );
                    return result;

                }

            }

            //Get the element set names and ids
            if ( increment[ "element_sets" ] ){
                error = _getKeyValuePairsFromYAML( increment, "element_sets", elementSetNames, elementSets );

                if ( error ){

                    errorOut result = new errorNode( "_writeMeshInformation", "Error when extracting the element sets" );
                    result->addNext( error );
                    return result;

                }
            }

        }

        error = _writer->writeIncrementMeshData( _currentIncrement, _collectionNumber, nodeIds, nodeSets, nodeSetNames, nodePositions,
                                                 elementIds, elementSets, elementSetNames, connectivity );

        if( error ){

            errorOut result = new errorNode( "_writeMeshInformation", "Error when writing the mesh information" );
            result->addNext( error );
            return result;

        }

        return NULL;

    }

    errorOut fileGenerator::_writeSolutionInformation( const YAML::Node &increment ){
        /*!
         * Write the solution information to the datafile
         *
         * :param const YAML::Node &increment: The increment information
         */

        //Initialize
        stringVector keys;
        floatMatrix values;

        if ( increment[ "node_values" ] ){

            //Extract the nodal solution vector
            errorOut error = _getKeyValuePairsFromYAML( increment, "node_values", keys, values );

            if ( error ){

                errorOut result = new errorNode( "_writeSolutionInformation", "Error in extracting the nodal solution vectors" );
                result->addNext( error );
                return result;

            }

            //Write out the solution vectors
            for ( unsigned int i = 0; i < keys.size( ); i++ ){

                error = _writer->writeScalarSolutionData( _currentIncrement, _collectionNumber, keys[ i ], "Node", values[ i ] );
                if ( error ){

                    errorOut result = new errorNode( "_writeSolutionInformation", "Error in writing out the nodal solution information with key " + keys[ i ] );
                    result->addNext( error );
                    return result;

                }

            }

        }

        if ( increment[ "cell_values" ] ){

            //Extract the cell solution vector
            errorOut error = _getKeyValuePairsFromYAML( increment, "cell_values", keys, values );

            if ( error ){

                errorOut result = new errorNode( "_writeSolutionInformation", "Error in extracting the cell solution vectors" );
                result->addNext( error );
                return result;

            }

            //Write out the solution vectors
            for ( unsigned int i = 0; i < keys.size( ); i++ ){

                error = _writer->writeScalarSolutionData( _currentIncrement, _collectionNumber, keys[ i ], "Node", values[ i ] );
                if ( error ){

                    errorOut result = new errorNode( "_writeSolutionInformation", "Error in writing out the cell solution information with key " + keys[ i ] );
                    result->addNext( error );
                    return result;

                }

            }


        }

        return NULL;
    }

    template< class T >
    errorOut fileGenerator::_getPropertyFromYAML( const YAML::Node &node, const std::string &property_name, T &property ){
        /*!
         * Extract a property by name from the YAML node
         *
         * :param const YAML::Node &node: The YAML node to query
         * :param const std::string &property_name: The name of the property
         * :param T &property: The property value
         */

        if ( !node[ property_name.c_str( ) ] ){

            std::string outstr = "property with name: ";
            outstr += property_name;
            outstr += " not found in YAML node";
            return new errorNode( "_getPropertyFromYAML", outstr );

        }

        try{

            property = node[ property_name ].as< T >( );

        }
        catch( std::exception &e ){

            return new errorNode( "_getPropertyFromYAML", e.what( ) );

        }

        return NULL;

    }

    template< class T >
    errorOut fileGenerator::_getKeyValuePairsFromYAML( const YAML::Node &node, const std::string &property_name,
                                                       stringVector &keys, std::vector< T > &values ){
        /*!
         * Get the key value pairs from the YAML node
         *
         * Note: each of the values must be of the same type
         *
         * :param const YAML::Node &node: The YAML node to query
         * :param const std::string &property_name: The name of the property to query
         * :param const std::string &keys: The keywords
         * :param std::vector< T > &values: The values
         */

        if ( !node[ property_name.c_str( ) ] ){

            std::string outstr = "property with name: ";
            outstr += property_name;
            outstr += " not found in YAML node";
            return new errorNode( "_getKeyValuePairsFromYAML", outstr );

        }

        keys.clear( );
        keys.reserve( node[ property_name.c_str( ) ].size( ) );
        values.reserve( node[ property_name.c_str( ) ].size( ) );

        for ( auto it = node[ property_name.c_str( ) ].begin( ); it != node[ property_name.c_str( ) ].end( ); it++ ){

            try{

                keys.push_back( it->first.as< std::string >( ) );
                values.push_back( it->second.as< T >( ) );

            }
            catch( std::exception &e ){

                return new errorNode( "_getKeyValuePairsFromYAML", e.what( ) );

            }

        }

        return NULL;

    }

    const uIntType *fileGenerator::getCurrentIncrement( ){
        /*!
         * Return a constant reference to the current increment
         */

        return &_currentIncrement;
    }

}
