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
#include<element.h>

namespace volumeReconstruction{

    typedef errorTools::Node errorNode; //!Redefinition for the error node
    typedef errorNode* errorOut; //!Redefinition for a pointer to the error node
    typedef double floatType; //!Define the float values type.
    typedef std::vector< floatType > floatVector; //! Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!Define a matrix of floats
    typedef int intType; //!Define the int values type.
    typedef std::vector< intType > intVector; //! Define a vector of ints
    typedef std::vector< std::vector< intType > > intMatrix; //!Define a matrix of ints
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
            KDNode( );
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
            virtual errorOut evaluate( );

            //Interface functions
            const floatVector *getPoints( );
            const floatVector *getFunction( );
            errorOut getFunctionValue( const unsigned int i, floatType &value );

            const floatVector *getLowerBounds( );
            const floatVector *getUpperBounds( );

        protected:

            unsigned int _dim = 3; //The dimension is hard coded to 3
                                   //It will be attempted to make the
                                   //code as general as possible but
                                   //only 3D will be explicitly 
                                   //implemented

            virtual errorOut initialize( );

            //Protected attributes
            YAML::Node _config;
            errorOut _error;
            std::string _type;

            const floatVector *_points;
            KDNode _tree;
            floatType _functionValue = 1;
            const floatVector *_functionValues = NULL;
            unsigned int _nPoints;

        private:

            std::shared_ptr< volumeReconstructionBase > create( const std::string &type );
            floatVector _upperBounds;
            floatVector _lowerBounds;

            errorOut setInterpolationConfiguration( );
            errorOut computeGeometryInformation( );
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
            ~dualContouring( );

            errorOut evaluate( );

        protected:

            errorOut initialize( );

            errorOut setGridSpacing( );

            errorOut interpolateFunctionToBackgroundGrid( );

        private:

            uIntVector _domainDiscretization;
            floatMatrix _gridLocations;
            floatType _exteriorRelativeDelta = 1e-3;

            floatType _absoluteTolerance = 1e-9;

            unsigned int _minApproximationCount = 5;
            bool _useMaterialPointsForNormals = false;

            bool _writeOutput = false;
            std::string _outputFilename = "output";

            errorOut processConfigurationFile( );

            errorOut processBackgroundGridElementImplicitFunction( const uIntVector &indices,
                                                                   floatVector &implicitFunctionNodalValues,
                                                                   uIntVector &globalNodeIds,
                                                                   uIntVector &pointCounts );
            
            errorOut getGridElement( const uIntVector &indices, std::unique_ptr< elib::Element > &element );

            errorOut projectImplicitFunctionToBackgroundGrid( );

            errorOut initializeInternalAndBoundaryCells( );

            errorOut findInternalAndBoundaryCells( );

            errorOut computeBoundaryPoints( );

            errorOut solveBoundLeastSquares( );

            std::string _elementType = "Hex8";
            floatType _isosurfaceCutoff = 0.5;
            floatVector _implicitFunctionValues;

            uIntVector _internalCells;
            uIntVector _boundaryCells;

            std::unordered_map< unsigned int, uIntVector > _boundaryEdges_x;
            std::unordered_map< unsigned int, uIntVector > _boundaryEdges_y;
            std::unordered_map< unsigned int, uIntVector > _boundaryEdges_z;

            floatVector _boundaryPoints;

            errorOut writeToXDMF( );

    };

    errorOut dualContouringInternalPointResidual( const floatVector &X, const floatMatrix &floatArgs, const intMatrix &intArgs,
                                                  floatVector &residual, floatMatrix &jacobian,
                                                  floatMatrix &floatOuts, intMatrix &intOuts );

}

#endif
