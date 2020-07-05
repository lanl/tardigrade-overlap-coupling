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

    enum registry { DUAL_CONTOURING };
    static std::map< std::string, registry > registryMap = { { "dual_contouring", registry::DUAL_CONTOURING } };

    //KD tree definitions
    class KDNode{
        /*!
         * Node which forms a simple KD tree
         */

        public:

            //Constructor
            KDNode( const floatVector *points, const uIntVector &ownedIndices,
                    const unsigned int &depth, const unsigned int &dim );

            const unsigned int* getIndex( );

            floatType getMinimumValueDimension( const unsigned int &d );

            floatType getMaximumValueDimension( const unsigned int &d );

            void getPointsInRange( const floatVector &upperBounds, const floatVector &lowerBounds,
                                   uIntVector &indices,
                                   floatVector *domainUpperBounds = NULL,
                                   floatVector *domainLowerBounds = NULL);

            void printData( const unsigned int &dim );

        private:
            const floatVector *_points;
            unsigned int _index;
            unsigned int _depth;
            unsigned int _axis;

            std::unique_ptr< KDNode > left_child = NULL;
            std::unique_ptr< KDNode > right_child = NULL;

    };

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

            ~volumeReconstructionBase( );
    
            //The factory for child classes
            std::shared_ptr< volumeReconstructionBase > create( );

            //Error output
            errorOut getError( );

            //Main functions
            errorOut loadPoints( const floatVector *points );
            errorOut loadFunction( const floatVector *function );

            //Interface functions
            const floatVector *getPoints( );
            const floatVector *getFunction( );

        protected:

            unsigned int _dim = 3; //The dimension is hard coded to 3

            virtual errorOut initialize( );

            //Protected attributes
            YAML::Node _config;
            errorOut _error;
            std::string _type;

            const floatVector *_points;
            floatType _functionValue;
            const floatVector *_functionValues = NULL;
            unsigned int _nPoints;

            errorOut setInterpolationConfiguration( );

        private:

            std::shared_ptr< volumeReconstructionBase > create( const std::string &type );
    };

    //Dual contouring class
    class dualContouring : public volumeReconstructionBase {
        /*!
         * The Dual contouring volume reconstruction method
         *
         * TODO: Consider moving this to its own header at some
         *       point. At this time, there is only one reconstruction
         *       type so it isn't a big deal.
         */

        public:

            //Constructors
            dualContouring( );
            dualContouring( const YAML::Node &configuration );

        protected:

            errorOut initialize( );

    };
}

#endif
