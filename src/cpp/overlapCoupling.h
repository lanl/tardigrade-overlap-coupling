/*!============================================================================
|                               overlapCoupling                               |
===============================================================================
| The implementation of the Overlap coupling method                           |
=============================================================================*/

#ifndef OVERLAPCOUPLING_H
#define OVERLAPCOUPLING_H

#include<DOFProjection.h>
#include<inputFileProcessor.h>
#include<element.h>
#include<volumeReconstruction.h>

namespace overlapCoupling{

    /*=========================================================================
    |                                Typedefs                                 |
    =========================================================================*/

    typedef inputFileProcessor::errorNode errorNode;
    typedef inputFileProcessor::errorOut errorOut;
    typedef inputFileProcessor::floatType floatType;
    typedef inputFileProcessor::floatVector floatVector;
    typedef inputFileProcessor::floatMatrix floatMatrix;
    typedef inputFileProcessor::uIntType uIntType;
    typedef inputFileProcessor::uIntVector uIntVector;
    typedef inputFileProcessor::uIntMatrix uIntMatrix;
    typedef inputFileProcessor::stringVector stringVector;
    typedef inputFileProcessor::DOFMap DOFMap;
    typedef DOFProjection::SparseMatrix SparseMatrix;
    typedef std::vector< DOFProjection::T > tripletVector;

    typedef std::unordered_map< std::string, floatType > domainFloatMap;
    typedef std::unordered_map< std::string, floatVector > domainFloatVectorMap;
    typedef std::unordered_map< uIntType, domainFloatMap > cellDomainFloatMap;
    typedef std::unordered_map< uIntType, domainFloatVectorMap > cellDomainFloatVectorMap;

    //!The different strategies for the partitioning coefficient
    enum partitioningCoefficient { VOLUME_FRACTION };
    const std::map< std::string, partitioningCoefficient > partitioningCoefficientStrategies =
        {
            { "volume_fraction", VOLUME_FRACTION }
        };

    class overlapCoupling{
        /*!
         * The implementation of the overlap coupling
         */

        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            //Constructors
            overlapCoupling( );

            overlapCoupling( const std::string &configurationFileName );

            errorOut setConfigurationFilename( const std::string &configurationFilename );

            errorOut getConstructorError( );

            errorOut initializeCoupling( );

            errorOut processIncrement( const unsigned int &microIncrement,
                                       const unsigned int &macroIncrement );

            errorOut processLastIncrements( );

            //Access functions
            const cellDomainFloatMap* getReferenceFreeMicroDomainMasses( );
            const cellDomainFloatMap* getReferenceGhostMicroDomainMasses( );
            const cellDomainFloatVectorMap* getReferenceFreeMicroDomainCentersOfMass( );
            const cellDomainFloatVectorMap* getReferenceGhostMicroDomainCentersOfMass( );

            const domainFloatMap* getFreeMicroDomainMasses( );
            const domainFloatMap* getGhostMicroDomainMasses( );
            const domainFloatVectorMap* getFreeMicroDomainCentersOfMass( );
            const domainFloatVectorMap* getGhostMicroDomainCentersOfMass( );

            const cellDomainFloatVectorMap* getReferenceFreeMicroDomainCenterOfMassShapeFunctions( );
            const cellDomainFloatVectorMap* getReferenceGhostMicroDomainCenterOfMassShapeFunctions( );

            const floatVector* getProjectedGhostMacroDisplacement( );
            const floatVector* getProjectedGhostMicroDisplacement( );

            DOFMap getMicroGlobalLocalNodeMap( );

            DOFMap getMacroGlobalLocalNodeMap( );

            floatVector getUpdatedMicroDisplacementDOF( );

            floatVector getUpdatedMacroDisplacementDOF( );

#ifdef TESTACCESS
            //Functions and attributes that will ONLY be accessible when compiled with TESTACCESS defined
            cellDomainFloatMap _test_domainMass;
            cellDomainFloatVectorMap _test_domainCOM;
            std::unordered_map< uIntType, std::unordered_map< std::string, std::unordered_map< uIntType, floatVector > > > _test_domainXi;
            cellDomainFloatVectorMap _test_domainCOMSF;
            std::unordered_map< uIntType, std::unordered_map< std::string, std::unordered_map< uIntType, floatVector > > > _test_domainMUP;

            const std::unordered_map< uIntType, floatType >* getMacroNodeProjectedMass( );
            const std::unordered_map< uIntType, floatVector >* getMacroNodeProjectedMassMomentOfInertia( );
            const std::unordered_map< uIntType, floatVector >* getMacroNodeMassRelativePositionConstant( );

            const std::unordered_map< uIntType, floatVector >* getMacroReferencePositions( );
            const std::unordered_map< uIntType, floatVector >* getMicroReferencePositions( );
#endif

        private:

            //Private attributes
            const unsigned int _dim = 3; //The dimension of the problem is hard-coded to be 3D
            errorOut _error;
            inputFileProcessor::inputFileProcessor _inputProcessor;
            floatType _absoluteTolerance = 1e-9;

            //Domain mass properties
//            floatVector _referenceFreeMicroDomainMasses;
//            floatVector _referenceGhostMicroDomainMasses;
//            floatVector _referenceFreeMicroDomainCentersOfMass;
//            floatVector _referenceGhostMicroDomainCentersOfMass;
            cellDomainFloatMap _referenceFreeMicroDomainMasses;
            cellDomainFloatMap _referenceGhostMicroDomainMasses;
            cellDomainFloatVectorMap _referenceFreeMicroDomainCentersOfMass;
            cellDomainFloatVectorMap _referenceGhostMicroDomainCentersOfMass;
            cellDomainFloatVectorMap _referenceFreeMicroDomainMomentsOfInertia;
            cellDomainFloatVectorMap _referenceGhostMicroDomainMomentsOfInertia;

//            floatVector _freeMicroDomainMasses;
//            floatVector _ghostMicroDomainMasses;
//            floatVector _freeMicroDomainCentersOfMass;
//            floatVector _ghostMicroDomainCentersOfMass;
            domainFloatMap       _freeMicroDomainMasses;
            domainFloatMap       _ghostMicroDomainMasses;
            domainFloatVectorMap _freeMicroDomainCentersOfMass;
            domainFloatVectorMap _ghostMicroDomainCentersOfMass;

//            floatVector _referenceFreeMicroDomainCenterOfMassShapeFunctions;
//            floatVector _referenceGhostMicroDomainCenterOfMassShapeFunctions;
            cellDomainFloatVectorMap _referenceFreeMicroDomainCenterOfMassShapeFunctions;
            cellDomainFloatVectorMap _referenceGhostMicroDomainCenterOfMassShapeFunctions;

            floatVector _projected_ghost_macro_displacement;
            floatVector _projected_ghost_micro_displacement;

//            floatVector _macroNodeProjectedMass;
//            floatVector _macroNodeProjectedMassMomentOfInertia;
//            floatVector _macroNodeMassRelativePositionConstant;
            std::unordered_map< uIntType, floatType >   _macroNodeProjectedMass;
            std::unordered_map< uIntType, floatVector > _macroNodeProjectedMassMomentOfInertia;
            std::unordered_map< uIntType, floatVector > _macroNodeMassRelativePositionConstant;


            std::unordered_map< uIntType, floatVector > _macroReferencePositions;
            std::unordered_map< uIntType, floatVector > _microReferencePositions;

            //Private functions
            errorOut processDomainMassData( const unsigned int &microIncrement, const std::string &domainName,
                                            domainFloatMap &domainMass, domainFloatVectorMap &domainCenterOfMass,
                                            domainFloatVectorMap &domainMomentsOfInertia,
                                            std::unordered_map< uIntType, floatVector > &domainXiVectors );

            //Compute initial values
            errorOut setReferenceStateFromIncrement( const unsigned int &microIncrement, const unsigned int &macroIncrement );

            //Compute the increment's values
            errorOut computeIncrementCentersOfMass( const unsigned int microIncrement, const unsigned int macroIncrement,
                                                    domainFloatMap &freeDomainMass, domainFloatMap &ghostDomainMass,
                                                    domainFloatVectorMap &freeDomainCM, domainFloatVectorMap &ghostDomainCM );

            //Compute the shape functions at the centers of mass
            errorOut buildMacroDomainElement( const unsigned int cellID,
                                              const std::unordered_map< uIntType, floatVector > &nodeLocations,
                                              const std::unordered_map< uIntType, uIntVector > &connectivity,
                                              std::unique_ptr< elib::Element > &element );            

            errorOut buildMacroDomainElement( const unsigned int cellID,
                                              const std::unordered_map< uIntType, floatVector > &nodeReferenceLocations,
                                              const std::unordered_map< uIntType, floatVector > &nodeDisplacements,
                                              const std::unordered_map< uIntType, uIntVector > &connectivity,
                                              std::unique_ptr< elib::Element > &element );            

            errorOut computeShapeFunctionsAtPoint( const unsigned int cellID,
                                                   const std::unordered_map< uIntType, floatVector > &nodeLocations,
                                                   const std::unordered_map< uIntType, uIntVector >  &connectivity,
                                                   const floatVector &point, floatVector &shapeFunction );

            errorOut computeShapeFunctionsAtPoint( const unsigned int cellID,
                                                   const std::unordered_map< uIntType, floatVector > &nodeReferenceLocations,
                                                   const std::unordered_map< uIntType, floatVector > &nodeDisplacements,
                                                   const std::unordered_map< uIntType, uIntVector >  &connectivity,
                                                   const floatVector &point, floatVector &shapeFunction );

            errorOut computeShapeFunctionsAtPoints( const unsigned int cellID,
                                                    const std::unordered_map< uIntType, floatVector > &nodeLocations,
                                                    const std::unordered_map< uIntType, uIntVector > &connectivity,
                                                    const std::unordered_map< uIntType, floatVector > &points,
                                                    std::unordered_map< uIntType, floatVector > &shapeFunctions );

            errorOut computeShapeFunctionsAtPoints( const unsigned int cellID,
                                                    const std::unordered_map< uIntType, floatVector > &nodeReferenceLocations,
                                                    const std::unordered_map< uIntType, floatVector > &nodeDisplacements,
                                                    const std::unordered_map< uIntType, uIntVector > &connectivity,
                                                    const std::unordered_map< uIntType, floatVector > &points,
                                                    std::unordered_map< uIntType, floatVector > &shapeFunctions );

            errorOut computeShapeFunctionGradientsAtPoints( const unsigned int cellID,
                                                            const std::unordered_map< uIntType, floatVector > &nodeReferenceLocations,
                                                            const std::unordered_map< uIntType, floatVector > &nodeDisplacements,
                                                            const std::unordered_map< uIntType, uIntVector > &connectivity,
                                                            const std::unordered_map< uIntType, floatVector > points,
                                                            std::unordered_map< uIntType, floatVector > &shapeFunctionGradients );

            errorOut computeShapeFunctionsAtReferenceCentersOfMass( );

            errorOut computeDomainShapeFunctionInformation( const unsigned int &cellID,
                                                            const std::string &domainName,
                                                            const unsigned int &microIncrement,
                                                            const floatVector &domainCenterOfMass,
                                                            floatVector &domainCenterOfMassShapeFunctionValues,
                                                            std::unordered_map< uIntType, floatVector > &domainMicroPositionShapeFunctionValues );

            errorOut addDomainContributionToInterpolationMatrix( const uIntVector  &domainNodes,
                                                                 const uIntVector  &macroNodes,
                                                                 const std::unordered_map< uIntType, floatVector > &domainReferenceXis,
                                                                 const floatVector &domainCenterOfMassShapeFunctionValues );

            errorOut processDomainReference( const unsigned int &microIncrement,
                                             const std::string &domainName,
                                             const unsigned int cellID, const uIntVector &macroNodes,
                                             domainFloatMap       &referenceMicroDomainMass,
                                             domainFloatVectorMap &referenceMicroDomainCentersOfMass,
                                             domainFloatVectorMap &referenceMicroDomainMomentsOfInertia,
                                             std::unordered_map< uIntType, floatVector > &domainReferenceXiVectors,
                                             floatVector &domainCenterOfMassShapeFunctionValues,
                                             std::unordered_map< uIntType, floatVector > &domainMicroPositionShapeFunctionValues );

            errorOut addDomainToDirectProjectionReferenceValues( const uIntVector &domainNodes,
                                                                 const uIntVector &macroNodes,
                                                                 const std::unordered_map< uIntType, floatVector > &domainReferenceXiVectors,
                                                                 const std::unordered_map< uIntType, floatVector > &domainMicroPositionShapeFunctionValues
                                                               );

            errorOut formTheProjectors( const unsigned int &microIncrement, const unsigned int &macroIncrement );

            errorOut formL2Projectors( );

            errorOut formDirectProjectionProjectors( const unsigned int &microIncrement, const unsigned int &macroIncrement );

            errorOut addDomainContributionToDirectFreeMicroToGhostMacroProjector( const unsigned int &cellIndex,
                                                                                  const unsigned int &cellID,
                                                                                  const unsigned int &microIncrement,
                                                                                  const std::string &domainName,
                                                                                  const uIntVector &macroNodes );
            errorOut projectDegreesOfFreedom( const bool useUpdatedFreeDOF = false );

            errorOut homogenizeMicroScale( const unsigned int &microIncrement );

            errorOut reconstructDomain( const unsigned int &microIncrement, const std::string &microDomainName,
                                        uIntVector &microDomainNodeIds, floatVector &microNodePositions,
                                        std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume );

            errorOut computeDomainVolumeAverages( const uIntType &macroCellID, const std::string &microDomainName,
                                                  const uIntVector &microDomainNodeIDs,
                                                  std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume, 
                                                  const floatVector *microDomainCenterOfMass = NULL );

            errorOut computeDomainSurfaceAverages( const uIntType &macroCellID, const std::string &microDomainName,
                                                   const uIntVector &microDomainNodeIDs,
                                                   const uIntType &microDomainSurfaceDecompositionCount,
                                                   std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume );

            errorOut computeHomogenizedStresses( const uIntType &macroCellName );

            errorOut assembleHomogenizedMatricesAndVectors( );

            errorOut assembleHomogenizedMassMatrix( );

            errorOut assembleHomogenizedInternalForceVector( );

            errorOut assembleHomogenizedExternalForceVector( );

            errorOut assembleFreeMicromorphicMassMatrix( );

            errorOut assembleCouplingMassAndDampingMatrices( );

            errorOut assembleCouplingForceVector( );

            errorOut solveFreeDisplacement( const bool updateGhostDOF );

            errorOut constructKineticEnergyPartitioningCoefficient( const uIntType &macroCellID,
                                                                    const std::unique_ptr< elib::Element > &element,
                                                                    floatVector &res );

            errorOut constructPotentialEnergyPartitioningCoefficient( std::unordered_map< uIntType, floatType > &qes );

            errorOut outputReferenceInformation( );

            errorOut outputHomogenizedResponse( const uIntType collectionNumber = 0 );

            errorOut extractProjectionMatricesFromFile( );

            errorOut writeReferenceMeshDataToFile( const uIntType collectionNumber = 0 );

            errorOut writeUpdatedDOFToFile( const uIntType collectionNumber = 0 );
                                                                   
            //The interpolation matrix
            SparseMatrix _N;

            //Construct the projection matrices for the L2 projection
            Eigen::MatrixXd _L2_BQhatQ;
            Eigen::MatrixXd _L2_BQhatD;
            Eigen::MatrixXd _L2_BDhatQ;
            Eigen::MatrixXd _L2_BDhatD;

            //Construct the projection matrices for the direct projection
            SparseMatrix _DP_BQhatQ;
            SparseMatrix _DP_BQhatD;
            SparseMatrix _DP_BDhatQ;
            SparseMatrix _DP_BDhatD;

            //The homogenized values
            cellDomainFloatMap homogenizedVolumes;
            cellDomainFloatMap homogenizedSurfaceAreas;
            cellDomainFloatMap homogenizedDensities;
            cellDomainFloatVectorMap homogenizedMicroInertias;
            cellDomainFloatVectorMap homogenizedCentersOfMass;
            cellDomainFloatVectorMap homogenizedBodyForces;
            cellDomainFloatVectorMap homogenizedBodyForceCouples;
            cellDomainFloatVectorMap homogenizedAccelerations;
            cellDomainFloatVectorMap homogenizedMicroSpinInertias;
            cellDomainFloatVectorMap homogenizedSymmetricMicroStresses;
            cellDomainFloatVectorMap homogenizedSurfaceRegionAreas;
//            std::unordered_map< uIntType, floatVector > homogenizedSurfaceRegionDensities;
            cellDomainFloatVectorMap homogenizedSurfaceRegionCentersOfMass;
            cellDomainFloatVectorMap homogenizedSurfaceRegionTractions;
            cellDomainFloatVectorMap homogenizedSurfaceRegionCouples;

            //The values at the macro cell quadrature points
            std::unordered_map< uIntType, floatVector > quadraturePointDensities;
            std::unordered_map< uIntType, floatVector > quadraturePointBodyForce;
            std::unordered_map< uIntType, floatVector > quadraturePointAccelerations;
            std::unordered_map< uIntType, floatVector > quadraturePointMicroInertias;
            std::unordered_map< uIntType, floatVector > quadraturePointBodyCouples;
            std::unordered_map< uIntType, floatVector > quadraturePointMicroSpinInertias;
            std::unordered_map< uIntType, floatVector > quadraturePointSymmetricMicroStress;
            std::unordered_map< uIntType, floatVector > quadraturePointCauchyStress;
            std::unordered_map< uIntType, floatVector > quadraturePointHigherOrderStress;

            //The external forces at the nodes
            std::unordered_map< uIntType, floatVector > externalForcesAtNodes;
            std::unordered_map< uIntType, floatVector > externalCouplesAtNodes;

            SparseMatrix homogenizedMassMatrix;
            Eigen::MatrixXd homogenizedFINT;
            Eigen::MatrixXd homogenizedFEXT;

            //The values of the macro-domains
            SparseMatrix freeMicromorphicMassMatrix;
            std::unordered_map< uIntType, floatType > macroKineticPartitioningCoefficient;

            Eigen::MatrixXd _L2_MASS;
            Eigen::MatrixXd _L2_DAMPING;

            SparseMatrix _DP_MASS;
            SparseMatrix _DP_DAMPING;

            Eigen::MatrixXd _FORCE;

            floatVector _updatedFreeMicroDispDOFValues;
            floatVector _updatedFreeMacroDispDOFValues;
            bool _freeDOFValuesUpdated;
            uIntType _currentReferenceOutputIncrement = 0;

    };

    errorOut MADOutlierDetection( const floatVector &x, uIntVector &outliers, const floatType threshold = 10,
                                  const floatType eps = 1e-9 );

    errorOut formMicromorphicElementMassMatrix( const std::unique_ptr< elib::Element > &element,
                                                const floatVector &degreeOfFreedomValues,
                                                const floatVector &momentOfInertia,
                                                const floatVector &density,
                                                const DOFMap *nodeIDToIndex,
                                                std::vector< DOFProjection::T > &coefficients );

    errorOut formMicromorphicElementInternalForceVector( const std::unique_ptr< elib::Element > &element,
                                                         const floatVector &degreeOfFreedomValues,
                                                         const floatVector &cauchyStress,
                                                         const floatVector &symmetricMicroStress,
                                                         const floatVector &higherOrderStress,
                                                         const DOFMap *nodeIDToIndex,
                                                         Eigen::MatrixXd &internalForceVector );

    errorOut computeMicromorphicElementRequiredValues( const std::unique_ptr< elib::Element > &element,
                                                       const elib::quadrature_rule::iterator &qpt,
                                                       const uIntType dim,
                                                       const floatMatrix &reshapedDOFValues,
                                                       const bool useReference,
                                                       floatVector &shapeFunctions,
                                                       floatMatrix &gradShapeFunctions,
                                                       floatVector &deformationGradient,
                                                       floatType &J, floatType &Jxw,
                                                       floatVector &uQpt, floatVector &XiQpt );

    errorOut writeSparseMatrixToXDMF( const SparseMatrix &A, const std::string matrixName,
                                      const std::string &filename, shared_ptr< XdmfDomain > &domain,
                                      shared_ptr< XdmfUnstructuredGrid > &grid );

    errorOut readSparseMatrixFromXDMF( const shared_ptr< XdmfUnstructuredGrid > &grid, const std::string &matrixName, SparseMatrix &A );

    errorOut runOverlapCoupling( const std::string &filename,
                                 DOFMap &microGlobalLocalNodeMap, floatVector &updatedMicroDisplacementDOF,
                                 DOFMap &macroGlobalLocalNodeMap, floatVector &updatedMacroDisplacementDOF
                               );

}

#endif
