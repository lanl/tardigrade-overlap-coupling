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
                                                const floatVector &domainMicroWeights, floatVector &microDisplacements );

    errorOut addMacroDomainDisplacementToMicro( const unsigned int dim, const uIntVector &domainMicroNodeIndices,
                                                const floatVector &u, const floatVector &phi,
                                                const floatVector &domainReferenceXis,
                                                const floatVector &domainMicroWeights, floatVector &microDisplacements );

    errorOut addDomainMicroContributionToMacroMass( const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                    const floatVector &microMasses, const floatVector &domainMicroShapeFunctions,
                                                    const floatVector &domainMicroWeights,
                                                    floatVector &projectedMicroMasses );

    errorOut addDomainMicroContributionToMacroMicroMassMomentOfInertia( const unsigned int &dim,
                                                                        const uIntVector &domainMicroNodeIndices,
                                                                        const uIntVector &domainMacroNodeIndices,
                                                                        const floatVector &domainReferenceXis,
                                                                        const floatVector &microMasses,
                                                                        const floatVector &domainMicroShapeFunctions,
                                                                        const floatVector &domainMicroWeights,
                                                                        floatVector &projectedMassMicroMomentOfInertia );

    errorOut addDomainMassConstant( const unsigned int &dim,
                                    const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                    const floatVector &domainReferenceXis, const floatVector &microMasses,
                                    const floatVector &domainMicroShapeFunctions, const floatVector &domainMicroWeights,
                                    floatVector &projectedMassConstant );

    errorOut addDomainMassDisplacement( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                        const floatVector &microMasses, const floatVector &domainMicroShapeFunctions,
                                        const floatVector &domainMicroWeights, const floatVector &microDisplacements,
                                        floatVector &projectedMassDisplacement );

    errorOut addDomainMassMicroDisplacementPosition( const unsigned int &dim,
                                                     const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                     const floatVector &domainReferenceXis, const floatVector &microMasses,
                                                     const floatVector &domainMicroShapeFunctions, const floatVector &domainMicroWeights,
                                                     const floatVector &microDisplacements,
                                                     floatVector &projectedMassDisplacementPosition );

    errorOut addDomainMicroToMacroProjectionTerms( const unsigned int &dim,
                                                   const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                   const floatVector &domainReferenceXis, const floatVector &microMasses,
                                                   const floatVector &domainMicroShapeFunctions, const floatVector &domainMicroWeights,
                                                   const floatVector &microDisplacements,
                                                   floatVector &projectedMassMicroMomentOfInertia,
                                                   floatVector &projectedMassConstant,
                                                   floatVector &projectedMassDisplacement,
                                                   floatVector &projectedMassDisplacementPosition,
                                                   const bool computeMassMomentOfInertia = true,
                                                   const bool computeMassConstant = true,
                                                   const bool computeMassMicroDisplacement = true,
                                                   const bool computeMassDisplacementPosition = true );

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microMasses,
                                        const floatVector &microPositions, const floatVector &domainMicroWeights,
                                        floatType &domainMass, floatVector &domainCM );

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microVolumes,
                                        const floatVector &microDensities, const floatVector &microPositions,
                                        const floatVector &domainMicroWeights, floatType &domainMass, floatVector &domainCM );

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microMasses,
                                        const floatVector &microPositions, const floatVector &domainMicroWeights,
                                        floatVector &domainCM );

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microVolumes,
                                        const floatVector &microDensities, const floatVector &microPositions,
                                        const floatVector &domainMicroWeights, floatVector &domainCM );

    errorOut computeDomainXis( const unsigned int &dim,
                               const uIntVector &domainMicroNodeIndices, const floatVector &microPositions,
                               const floatVector &domainCM, floatVector &domainXis );

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
                                                        const floatVector &domainMicroWeights, SparseMatrix &domainN );
}

#endif
