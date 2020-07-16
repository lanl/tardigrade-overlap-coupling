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
    typedef inputFileProcessor::stringVector stringVector;
    typedef inputFileProcessor::DOFMap DOFMap;
    typedef DOFProjection::SparseMatrix SparseMatrix;

    class overlapCoupling{
        /*!
         * The implementation of the overlap coupling
         */

        public:

            //Constructors
            overlapCoupling( );

            overlapCoupling( const std::string &configurationFileName );

            errorOut setConfigurationFilename( const std::string &configurationFilename );

            errorOut getConstructorError( );

            errorOut initializeCoupling( );

            errorOut processIncrement( const unsigned int &microIncrement,
                                       const unsigned int &macroIncrement );

            //Access functions
            const floatVector* getReferenceFreeMicroDomainMasses( );
            const floatVector* getReferenceGhostMicroDomainMasses( );
            const floatVector* getReferenceFreeMicroDomainCentersOfMass( );
            const floatVector* getReferenceGhostMicroDomainCentersOfMass( );

            const floatVector* getFreeMicroDomainMasses( );
            const floatVector* getGhostMicroDomainMasses( );
            const floatVector* getFreeMicroDomainCentersOfMass( );
            const floatVector* getGhostMicroDomainCentersOfMass( );

            const floatVector* getReferenceFreeMicroDomainCenterOfMassShapeFunctions( );
            const floatVector* getReferenceGhostMicroDomainCenterOfMassShapeFunctions( );

            const floatVector* getProjectedGhostMacroDisplacement( );
            const floatVector* getProjectedGhostMicroDisplacement( );

        private:

            //Private attributes
            const unsigned int _dim = 3; //The dimension of the problem is hard-coded to be 3D
            errorOut _error;
            inputFileProcessor::inputFileProcessor _inputProcessor;

            //Domain mass properties
            floatVector _referenceFreeMicroDomainMasses;
            floatVector _referenceGhostMicroDomainMasses;
            floatVector _referenceFreeMicroDomainCentersOfMass;
            floatVector _referenceGhostMicroDomainCentersOfMass;

            floatVector _freeMicroDomainMasses;
            floatVector _ghostMicroDomainMasses;
            floatVector _freeMicroDomainCentersOfMass;
            floatVector _ghostMicroDomainCentersOfMass;

            floatVector _referenceFreeMicroDomainCenterOfMassShapeFunctions;
            floatVector _referenceGhostMicroDomainCenterOfMassShapeFunctions;

            floatVector _projected_ghost_macro_displacement;
            floatVector _projected_ghost_micro_displacement;

            floatVector _macroNodeProjectedMass;
            floatVector _macroNodeProjectedMassMomentOfInertia;
            floatVector _macroNodeMassRelativePositionConstant;

            floatVector _macroReferencePositions;
            floatVector _microReferencePositions;

            //Private functions
            errorOut processDomainMassData( const unsigned int &microIncrement, const std::string &domainName,
                                            floatType &domainMass, floatVector &domainCenterOfMass,
                                            floatVector &domainXiVectors );

            //Compute initial values
            errorOut setReferenceStateFromIncrement( const unsigned int &microIncrement, const unsigned int &macroIncrement );

            //Compute the increment's values
            errorOut computeIncrementCentersOfMass( const unsigned int microIncrement, const unsigned int macroIncrement,
                                                    floatVector &freeDomainMass, floatVector &ghostDomainMass,
                                                    floatVector &freeDomainCM, floatVector &ghostDomainCM );

            //Compute the shape functions at the centers of mass
            errorOut computeShapeFunctionsAtPoints( const unsigned int cellID, const floatVector &nodeLocations,
                                                    const uIntVector &connectivity, const uIntVector &connectivityCellIndices,
                                                    const floatVector &points, floatVector &shapeFunctions );

            errorOut computeShapeFunctionsAtPoints( const unsigned int cellID,
                                                    const floatVector &nodeReferenceLocations, const floatVector &nodeDisplacements,
                                                    const uIntVector &connectivity, const uIntVector &connectivityCellIndices,
                                                    const floatVector &points, floatVector &shapeFunctions );

            errorOut computeShapeFunctionsAtReferenceCentersOfMass( );

            errorOut computeDomainShapeFunctionInformation( const unsigned int &cellID,
                                                            const std::string &domainName,
                                                            const unsigned int &microIncrement,
                                                            const floatVector &domainCenterOfMass,
                                                            floatVector &domainCenterOfMassShapeFunctionValues,
                                                            floatVector &domainMicroPositionShapeFunctionValues );

            errorOut addDomainContributionToInterpolationMatrix( const uIntVector  &domainNodes,
                                                                 const uIntVector  &macroNodes,
                                                                 const floatVector &domainReferenceXis,
                                                                 const floatVector &domainCenterOfMassShapeFunctionValues );

            errorOut processDomainReference( const unsigned int &microIncrement,
                                             const unsigned int &domainIndex, const std::string &domainName,
                                             const unsigned int cellID, const uIntVector &macroNodes,
                                             floatType   &referenceMicroDomainMass,
                                             floatVector &referenceMicroDomainCentersOfMass,
                                             floatVector &domainReferenceXiVectors,
                                             floatVector &domainCenterOfMassShapeFunctionValues,
                                             floatVector &domainMicroPositionShapeFunctionValues );

            errorOut addDomainToDirectProjectionReferenceValues( const uIntVector &domainNodes,
                                                                 const uIntVector &macroNodes,
                                                                 const floatVector &domainReferenceXiVectors,
                                                                 const floatVector &domainMicroPositionShapeFunctionValues
                                                               );

            errorOut formTheProjectors( const unsigned int &microIncrement, const unsigned int &macroIncrement );

            errorOut formL2Projectors( );

            errorOut formDirectProjectionProjectors( const unsigned int &microIncrement, const unsigned int &macroIncrement );

            errorOut addDomainContributionToDirectFreeMicroToGhostMacroProjector( const unsigned int &cellIndex,
                                                                                  const unsigned int &cellID,
                                                                                  const unsigned int &microIncrement,
                                                                                  const std::string &domainName,
                                                                                  const uIntVector &macroNodes );
            errorOut projectDegreesOfFreedom( );

            errorOut homogenizeMicroScale( const unsigned int &microIncrement );

            errorOut reconstructDomain( const unsigned int &microIncrement, const std::string &microDomainName,
                                        uIntVector &microDomainNodeIds, floatVector &microNodePositions,
                                        std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume );

            errorOut computeDomainVolumeAverages( const uIntType &macroCellName, const std::string &microDomainName,
                                                  const uIntVector &microDomainNodeIDs,
                                                  std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume );

            errorOut computeDomainSurfaceAverages( const uIntType &macroCellName, const std::string &microDomainName,
                                                   std::shared_ptr< volumeReconstruction::volumeReconstructionBase > &reconstructedVolume );

            errorOut computeHomogenizedStresses( const uIntType &macroCellName );
                                                                   
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
            std::unordered_map< uIntType, floatVector > homogenizedVolumes;
            std::unordered_map< uIntType, floatVector > homogenizedDensities;

    };

}

#endif
