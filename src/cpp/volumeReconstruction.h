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
    typedef elib::uitype uIntType;
    typedef std::vector< uIntType > uIntVector; //!Define a vector of unsigned ints
    typedef std::vector< uIntVector > uIntMatrix; //!Define a matrix of unsigned ints
    typedef std::vector< std::string > stringVector; //!Define a vector of strings
    typedef std::unordered_map< uIntType, uIntType > DOFMap;

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
                    const uIntType &depth, const uIntType &dim );

            const uIntType* getIndex( );

            floatType getMinimumValueDimension( const uIntType &d );

            floatType getMaximumValueDimension( const uIntType &d );

            void getPointsInRange( const floatVector &upperBounds, const floatVector &lowerBounds,
                                   uIntVector &indices,
                                   floatVector *domainUpperBounds = NULL,
                                   floatVector *domainLowerBounds = NULL);

            void getPointsWithinRadiusOfOrigin( const floatVector &origin, const floatType &radius,
                                                uIntVector &indices,
                                                floatVector *domainUpperBounds = NULL,
                                                floatVector *domainLowerBounds = NULL );

            void printData( const uIntType &dim );

        private:
            const floatVector *_points;
            uIntType _index;
            uIntType _depth;
            uIntType _axis;

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

            errorOut addBoundingPlanes( const floatMatrix &planePoints, const floatMatrix &planeNormals );

            //Interface functions
            const floatVector *getPoints( );
            const floatVector *getFunction( );
            errorOut getFunctionValue( const uIntType i, floatType &value );

            const floatVector *getLowerBounds( );
            const floatVector *getUpperBounds( );
            const floatType   *getMedianNeighborhoodDistance( );

            //Required overloads
            virtual errorOut evaluate( );
            virtual errorOut performVolumeIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                       floatVector &integratedValue );

            virtual errorOut performRelativePositionVolumeIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                                       const floatVector &origin, floatVector &integratedValue );

            virtual errorOut performSurfaceIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                        floatVector &integratedValue, const uIntVector *subdomainIDs = NULL,
                                                        const floatVector *subdomainWeights = NULL,
                                                        const floatVector *macroNormal = NULL, const bool useMacroNormal = false );

            virtual errorOut performPositionWeightedSurfaceIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                                        floatVector &integratedValue, const uIntVector *subdomainIDs = NULL,
                                                                        const floatVector *subdomainWeights = NULL,
                                                                        const floatVector *macroNormal = NULL, const bool useMacroNormal = false );

            virtual errorOut performSurfaceFluxIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                            floatVector &integratedValue, const uIntVector *subdomainIDs = NULL,
                                                            const floatVector *subdomainWeights = NULL,
                                                            const floatVector *macroNormal = NULL, const bool useMacroNormal = false );

            virtual errorOut performRelativePositionSurfaceFluxIntegration( const floatVector &valuesAtPoints,
                                                                            const uIntType valueSize,
                                                                            const floatVector &origin,
                                                                            floatVector &integratedValue,
                                                                            const uIntVector *subdomainIDs = NULL,
                                                                            const floatVector *subdomainWeights = NULL,
                                                                            const floatVector *macroNormal = NULL,
                                                                            const bool useMacroNormal = false );

            virtual errorOut getSurfaceSubdomains( const floatType &minDistance, uIntVector &subdomainNodeCounts,
                                                   uIntVector &subdomainIDs );

            virtual const uIntVector *getBoundaryIDs( );

            virtual const floatVector *getBoundaryPoints( );

            const YAML::Node exportConfiguration( );

            virtual errorOut writeToXDMF( );

        protected:

            uIntType _dim = 3; //The dimension is hard coded to 3
                               //It will be attempted to make the
                               //code as general as possible but
                               //only 3D will be explicitly 
                               //implemented

            virtual errorOut initialize( );

            void setEvaluated( const bool isEvaluated );
            bool getEvaluated( );

            //Protected attributes
            YAML::Node _config;
            errorOut _error;
            std::string _type;

            const floatVector *_points;
            KDNode _pointTree;
            floatType _functionValue = 1;
            const floatVector *_functionValues = NULL;
            uIntType _nPoints;
            uIntType _nNeighborhoodPoints;
            floatType _length_scale = 1;
            floatType _critical_radius = 1;

            std::vector< std::pair< floatVector, floatVector > > _boundingPlanes;
            bool _boundingSurfaces = false;
            

        private:

            bool _isEvaluated = false;

            std::shared_ptr< volumeReconstructionBase > create( const std::string &type );
            floatVector _upperBounds;
            floatVector _lowerBounds;
            floatType _medianNeighborhoodDistance;

            errorOut setInterpolationConfiguration( );
            errorOut computeGeometryInformation( );
            errorOut computeMedianNeighborhoodDistance( );
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

            errorOut performVolumeIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                               floatVector &integratedValue );

            errorOut performSurfaceIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                floatVector &integratedValue, const uIntVector *subdomainIDs = NULL,
                                                const floatVector *subdomainWeights = NULL,
                                                const floatVector *macroNormal = NULL, const bool useMacroNormal = false );

            errorOut performPositionWeightedSurfaceIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                                floatVector &integratedValue, const uIntVector *subdomainIDs = NULL,
                                                                const floatVector *subdomainWeights = NULL,
                                                                const floatVector *macroNormal = NULL, const bool useMacroNormal = false );

            errorOut performSurfaceFluxIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                    floatVector &integratedValue, const uIntVector *subdomainIDs = NULL,
                                                    const floatVector *subdomainWeights = NULL,
                                                    const floatVector *macroNormal = NULL, const bool useMacroNormal = false );

            errorOut performRelativePositionSurfaceFluxIntegration( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                                    const floatVector &origin, floatVector &integratedValue,
                                                                    const uIntVector *subdomainIDs = NULL,
                                                                    const floatVector *subdomainWeights = NULL,
                                                                    const floatVector *macroNormal = NULL,
                                                                    const bool useMacroNormal = false );

            errorOut getSurfaceSubdomains( const floatType &minDistance, uIntVector &subdomainNodeCounts,
                                           uIntVector &subdomainIDs );

            const uIntVector *getBoundaryIDs( );

            const floatVector *getBoundaryPoints( );

        protected:

            errorOut initialize( );

            errorOut setGridSpacing( );

            errorOut interpolateFunctionToBackgroundGrid( const floatVector &functionValuesAtPoints, const uIntType &functionDim,
                                                          std::unordered_map< uIntType, floatVector > &functionAtGrid );

            errorOut performSurfaceIntegralMethods( const floatVector &valuesAtPoints, const uIntType valueSize,
                                                    const floatVector &origin, floatVector &integratedValue,
                                                    bool computeFlux, bool positionWeightedIntegral, bool dyadWithOrigin, const uIntVector *subdomainIDs = NULL,
                                                    const floatVector *subdomainWeights = NULL, const floatVector *macroNormal = NULL,
                                                    const bool useMacroNormal = false );

            errorOut rbf( const floatVector &x, const floatVector &x0, const floatType &ls, floatType &val );

            errorOut grad_rbf( const floatVector &x, const floatVector &x0, const floatType &ls, floatVector &grad );

        private:

            uIntVector _domainDiscretization;
            floatMatrix _gridLocations;
            floatType _exteriorRelativeDelta = 1e-3;

            floatType _absoluteTolerance = 1e-9;

            uIntType _minPointsPerCell = 2;
            uIntType _minNormalApproximationCount = 5;
            bool _useMaterialPointsForNormals = false;

            bool _writeOutput = false;
            std::string _XDMFOutputFilename = "output";

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

            std::unordered_map< uIntType, uIntVector > _boundaryEdges_x;
            std::unordered_map< uIntType, uIntVector > _boundaryEdges_y;
            std::unordered_map< uIntType, uIntVector > _boundaryEdges_z;

            floatVector _boundaryPoints;
            DOFMap _boundaryPointIDToIndex;

            std::unordered_map< uIntType, floatType > _boundaryPointAreas; 
            std::unordered_map< uIntType, floatVector > _boundaryPointNormals;
            errorOut computeBoundaryPointNormalsAndAreas( );
            errorOut processBoundaryEdges( const std::unordered_map< uIntType, uIntVector > &boundaryEdges );

            KDNode _boundaryPointTree;

            errorOut writeToXDMF( );

    };

    errorOut dualContouringInternalPointResidual( const floatVector &X, const floatMatrix &floatArgs, const intMatrix &intArgs,
                                                  floatVector &residual, floatMatrix &jacobian,
                                                  floatMatrix &floatOuts, intMatrix &intOuts );
}

#endif
