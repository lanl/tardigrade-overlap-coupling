/*=============================================================================
|                             volumeReconstruction                            |
===============================================================================
| Tools to reconstruct volume information from pointsets                      |
=============================================================================*/

#ifndef VOLUMERECONSTRUCTION_H
#define VOLUMERECONSTRUCTION_H

#define USE_EIGEN
#include<vector_tools.h>
#include<error_tools.h>
#include<yaml-cpp/yaml.h>
#include<unordered_map>

namespace volumeReconstruction{

    typedef errorTools::Node errorNode; //!Redefinition for the error node
    typedef errorNode* errorOut; //!Redefinition for a pointer to the error node
    typedef double floatType; //!Define the float values type.
    typedef std::vector< floatType > floatVector; //! Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!Define a matrix of floats
    typedef std::vector< unsigned int > uIntVector; //!Define a vector of unsigned ints
    typedef std::vector< std::string > stringVector; //!Define a vector of strings
    typedef std::unordered_map< unsigned int, unsigned int > DOFMap;

    //All new volume reconstructors must be defined both in the registry enum 
    //and in the map to allow them to be accessed by strings. They must also 
    //be registerd in volumeReconstructionBase::create

    enum registry { };
    static std::map< std::string, registry > registryMap = {};

    class volumeReconstructionBase {
        /*!
         * The base class for volume reconstruction from pointsets. This allows
         * for many different approaches to be tried without having to 
         * re-implement an interface
         */

        public:
            //Constructors
            volumeReconstructionBase( );
            volumeReconstructionBase( const YAML::Node &config );
            volumeReconstructionBase( const YAML::Node &config, errorOut error );
    
            //The factory for child classes
            std::shared_ptr< volumeReconstructionBase > create( );
            std::shared_ptr< volumeReconstructionBase > create( const std::string &type );

            errorOut getError( );

        private:

            YAML::Node _config;
            errorOut _error;
            std::string _type;
    };
}

#endif
