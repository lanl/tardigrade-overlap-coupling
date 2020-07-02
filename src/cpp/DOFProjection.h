/*!============================================================================
|                                DOFProjection.h                              |
|=============================================================================|
| A collection of tools which can be used for degree of freedom projection in |
| a micromorphic context. The techniques are based on those of Wagner and Liu |
| [2003], Kadowaki and Liu [2004] and Klein and Zimmerman [2006] modified by  |
| Regueiro [2012] and Miller [2020].                                          |
=============================================================================*/

#ifndef DOFPROJECTION_H
#define DOFPROJECTION_H

#define USE_EIGEN
#include<vector_tools.h>
#include<error_tools.h>
#include<Eigen/Sparse>
#include<unordered_map>

namespace DOFProjection{

    /*=======================================================
    |                       Typedefs                        |
    =========================================================
    | Typedefs for use in the degree of freedom projections |
    =======================================================*/

    typedef errorTools::Node errorNode; //!Redefinition for the error node
    typedef errorNode* errorOut; //!Redefinition for a pointer to the error node
    typedef double floatType; //!Define the float values type.
    typedef std::vector< floatType > floatVector; //! Define a vector of floats
    typedef std::vector< std::vector< floatType > > floatMatrix; //!Define a matrix of floats
    typedef std::vector< unsigned int > uIntVector; //!Define a vector of unsigned ints

    //Eigen Typedefs
    typedef Eigen::SparseMatrix< floatType > SparseMatrix;
    typedef Eigen::Triplet< floatType > T;

    /*==================================================================
    |                       Projection Functions                       |
    ====================================================================
    | Functions which project the values from the macro to micro-scale |
    ==================================================================*/

    errorOut addMacroDomainDisplacementToMicro( const unsigned int dim,
                                                const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                const floatVector &domainReferenceXis,
                                                const floatVector &domainMacroInterpolationFunctionValues,
                                                const unsigned int &nMacroDOF, const floatVector &macroDOFVector,
                                                const floatVector &microWeights, floatVector &microDisplacements,
                                                const std::unordered_map< unsigned int, unsigned int >* microNodeToLocalIndex = NULL );

    errorOut addMacroDomainDisplacementToMicro( const unsigned int dim, const uIntVector &domainMicroNodeIndices,
                                                const floatVector &u, const floatVector &phi,
                                                const floatVector &domainReferenceXis,
                                                const floatVector &microWeights, floatVector &microDisplacements,
                                                const std::unordered_map< unsigned int, unsigned int >* microNodeToLocalIndex = NULL );

    errorOut addDomainMicroContributionToMacroMass( const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                    const floatVector &microMasses, const floatVector &domainMicroShapeFunctions,
                                                    const floatVector &microWeights,
                                                    floatVector &projectedMicroMasses,
                                                    const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL );

    errorOut addDomainMicroContributionToMacroMass( const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                    const floatVector &microVolumes, const floatVector &microDensities,
                                                    const floatVector &domainMicroShapeFunctions,
                                                    const floatVector &microWeights,
                                                    floatVector &projectedMicroMasses,
                                                    const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL );

    errorOut addDomainMicroContributionToMacroMicroMassMomentOfInertia( const unsigned int &dim,
                                                                        const uIntVector &domainMicroNodeIndices,
                                                                        const uIntVector &domainMacroNodeIndices,
                                                                        const floatVector &domainReferenceXis,
                                                                        const floatVector &microMasses,
                                                                        const floatVector &domainMicroShapeFunctions,
                                                                        const floatVector &microWeights,
                                                                        floatVector &projectedMassMicroMomentOfInertia,
                                                                        const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL
                                                                      );

    errorOut addDomainMassConstant( const unsigned int &dim,
                                    const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                    const floatVector &domainReferenceXis, const floatVector &microMasses,
                                    const floatVector &domainMicroShapeFunctions, const floatVector &microWeights,
                                    floatVector &projectedMassConstant,
                                    const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL
                                  );

    errorOut addDomainMassDisplacement( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                        const floatVector &microMasses, const floatVector &domainMicroShapeFunctions,
                                        const floatVector &microWeights, const floatVector &microDisplacements,
                                        floatVector &projectedMassDisplacement,
                                        const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL
                                      );

    errorOut addDomainMassMicroDisplacementPosition( const unsigned int &dim,
                                                     const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                     const floatVector &domainReferenceXis, const floatVector &microMasses,
                                                     const floatVector &domainMicroShapeFunctions, const floatVector &microWeights,
                                                     const floatVector &microDisplacements,
                                                     floatVector &projectedMassDisplacementPosition,
                                                     const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL
                                                   );

    errorOut addDomainMicroToMacroProjectionTerms( const unsigned int &dim,
                                                   const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                   const floatVector &domainReferenceXis, const floatVector &microMasses,
                                                   const floatVector &domainMicroShapeFunctions, const floatVector &microWeights,
                                                   const floatVector &microDisplacements,
                                                   floatVector &projectedMassMicroMomentOfInertia,
                                                   floatVector &projectedMassConstant,
                                                   floatVector &projectedMassDisplacement,
                                                   floatVector &projectedMassDisplacementPosition,
                                                   const bool computeMassMomentOfInertia = true,
                                                   const bool computeMassConstant = true,
                                                   const bool computeMassMicroDisplacement = true,
                                                   const bool computeMassDisplacementPosition = true,
                                                   const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL
                                                 );

    errorOut addDomainMicroContributionToMacroMicroMassMomentOfInertia( const unsigned int &dim,
                                                                        const uIntVector &domainMicroNodeIndices,
                                                                        const uIntVector &domainMacroNodeIndices,
                                                                        const floatVector &domainReferenceXis,
                                                                        const floatVector &microVolumes,
                                                                        const floatVector &microDensities,
                                                                        const floatVector &domainMicroShapeFunctions,
                                                                        const floatVector &microWeights,
                                                                        floatVector &projectedMassMicroMomentOfInertia,
                                                                        const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL
                                                                      );

    errorOut addDomainMassConstant( const unsigned int &dim,
                                    const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                    const floatVector &domainReferenceXis, const floatVector &microVolumes,
                                    const floatVector &microDensities, const floatVector &domainMicroShapeFunctions,
                                    const floatVector &microWeights, floatVector &projectedMassConstant,
                                    const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL
                                  );

    errorOut addDomainMassDisplacement( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                        const floatVector &microVolumes, const floatVector &microDensities, 
                                        const floatVector &domainMicroShapeFunctions, const floatVector &microWeights,
                                        const floatVector &microDisplacements, floatVector &projectedMassDisplacement,
                                        const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL
                                      );

    errorOut addDomainMassMicroDisplacementPosition( const unsigned int &dim,
                                                     const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                     const floatVector &domainReferenceXis, const floatVector &microVolumes,
                                                     const floatVector &microDensities, const floatVector &domainMicroShapeFunctions,
                                                     const floatVector &microWeights, const floatVector &microDisplacements,
                                                     floatVector &projectedMassDisplacementPosition,
                                                     const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL );

    errorOut addDomainMicroToMacroProjectionTerms( const unsigned int &dim,
                                                   const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                   const floatVector &domainReferenceXis, const floatVector &microVolumes,
                                                   const floatVector &microDensities, const floatVector &domainMicroShapeFunctions,
                                                   const floatVector &microWeights,
                                                   const floatVector &microDisplacements,
                                                   floatVector &projectedMassMicroMomentOfInertia,
                                                   floatVector &projectedMassConstant,
                                                   floatVector &projectedMassDisplacement,
                                                   floatVector &projectedMassDisplacementPosition,
                                                   const bool computeMassMomentOfInertia = true,
                                                   const bool computeMassConstant = true,
                                                   const bool computeMassMicroDisplacement = true,
                                                   const bool computeMassDisplacementPosition = true,
                                                   const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL
                                                 );

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microMasses,
                                        const floatVector &microPositions, const floatVector &microWeights,
                                        floatType &domainMass, floatVector &domainCM );

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microVolumes,
                                        const floatVector &microDensities, const floatVector &microPositions,
                                        const floatVector &microWeights, floatType &domainMass, floatVector &domainCM );

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microVolumes,
                                        const floatVector &microDensities, const floatVector &microReferencePositions,
                                        const floatVector &microDisplacements, const floatVector &microWeights,
                                        floatType &domainMass, floatVector &domainCM );

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microMasses,
                                        const floatVector &microPositions, const floatVector &microWeights,
                                        floatVector &domainCM );

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microVolumes,
                                        const floatVector &microDensities, const floatVector &microPositions,
                                        const floatVector &microWeights, floatVector &domainCM );

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microVolumes,
                                        const floatVector &microDensities, const floatVector &microReferencePositions,
                                        const floatVector &microDisplacements, const floatVector &microWeights,
                                        floatVector &domainCM );

    errorOut computeDomainXis( const unsigned int &dim,
                               const uIntVector &domainMicroNodeIndices, const floatVector &microPositions,
                               const floatVector &domainCM, floatVector &domainXis );

    errorOut computeDomainXis( const unsigned int &dim,
                               const uIntVector &domainMicroNodeIndices, const floatVector &microReferencePositions,
                               const floatVector &microDisplacements, const floatVector &domainCM, floatVector &domainXis );

    /*===================================================
    |                Projection Matrices                |
    =====================================================
    | Functions which construct the projection matrices |
    ===================================================*/

    errorOut formMacroDomainToMicroInterpolationMatrix( const unsigned int &dim,
                                                        const unsigned int &nMicroNodes, const unsigned int &nMacroNodes,
                                                        const uIntVector &domainMicroNodeIndices,
                                                        const uIntVector &domainMacroNodeIndices,
                                                        const floatVector &domainReferenceXis,
                                                        const floatVector &domainMacroInterpolationFunctionValues,
                                                        const floatVector &microWeights, SparseMatrix &domainN,
                                                        const std::unordered_map< unsigned int, unsigned int >* microNodeToLocalIndex = NULL,
                                                        const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL );

    errorOut formMicroDomainToMacroProjectionMatrix( const unsigned int &dim,
                                                     const unsigned int nMicroNodes,
                                                     const unsigned int nMacroNodes,
                                                     const uIntVector  &domainMicroNodeIndices,
                                                     const uIntVector  &domainMacroNodeIndices,
                                                     const floatVector &microVolumes,
                                                     const floatVector &microDensities,
                                                     const floatVector &microWeights,
                                                     const floatVector &domainReferenceXiVectors,
                                                     const floatVector &domainInterpolationFunctionValues,
                                                     const floatVector &domainMacroNodeProjectedMass,
                                                     const floatVector &domainMacroNodeProjectedMassMomentOfInertia,
                                                     const floatVector &domainMacroNodeMassRelativePositionConstant,
                                                     SparseMatrix &projector,
                                                     const std::unordered_map< unsigned int, unsigned int >* microNodeToLocalIndex = NULL,
                                                     const std::unordered_map< unsigned int, unsigned int >* macroNodeToLocalIndex = NULL );
}


#endif
