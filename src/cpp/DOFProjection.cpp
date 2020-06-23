/*!============================================================================
|                               DOFProjection.cpp                             |
|=============================================================================|
| A collection of tools which can be used for degree of freedom projection in |
| a micromorphic context. The techniques are based on those of Wagner and Liu |
| [2003], Kadowaki and Liu [2004] and Klein and Zimmerman [2006] modified by  |
| Regueiro [2012] and Miller [2020].                                          |
=============================================================================*/

#include<DOFProjection.h>

namespace DOFProjection{

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
                                                const floatVector &domainMicroWeights, floatVector &microDisplacements ){

        /*!
         * Add the contribution of a macro domain's deformation to the micro-scale
         *
         * :param const unsigned int dim: The dimension of the problem ( only 3 is tested )
         * :param const uIntVector &domainMicroDOFIndices: The indices of the micro-scale nodes present in the domain
         * :param const uIntVector &domainMacroDOFIndices: The indices of the macro-scale nodes present in the domain
         *     these are ( in a macro-scale FEA context ) the nodes of the micromorphic finite element.
         * :param const floatVector &domainReferenceXis: The Xi vectors of the micro-scale nodes which go from the local
         *     center of mass to the position of the micro-scale node in the reference configuration.
         * :param const floatVector &domainMacroInterpolationFunctionValues: The values of the interpolation functions
         *     at the local center of mass.
         * :param const unsigned int &nMacroDOF: The number of degrees of freedom associated with each macro-scale node.
         *     ( only 12 is tested )
         * :param const floatVector &macroDOFVector: The degree of freedom vector for the macro-scale. The ordering for each node
         *     is assumed to be [ u1, u2, u3, phi11, phi12, phi13, phi21, phi22, phi23, phi31, phi32, phi33, ... ]
         * :param const floatVector &domainMicroWeights: The weight associated with each micro-scale node for this domain. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 projection is being used.
         * :param floatVector &microDisplacements: The displacements of the micro-scale nodes as determined by the macro-scale
         *     projection.
         */

        if ( domainMacroNodeIndices.size() != domainMacroInterpolationFunctionValues.size() ){
            return new errorNode( "projectMacroDomainDisplacementToMicro",
                                  "The macro-scale node indices and the macro-scale interpolation function values must be the same length" );
        }

        //The interpolated degree of freedom vector
        floatVector interpolatedMacroDOF( nMacroDOF, 0 );

        //Compute the interpolated value of u and phi
        for ( unsigned int i = 0; i < domainMacroNodeIndices.size(); i++ ){
            interpolatedMacroDOF += domainMacroInterpolationFunctionValues[ i ]
                                  * floatVector( macroDOFVector.begin() + nMacroDOF * domainMacroNodeIndices[ i ],
                                                 macroDOFVector.begin() + nMacroDOF * domainMacroNodeIndices[ i ] + nMacroDOF );
        }

        //The interpolated macro-displacement
        floatVector u( interpolatedMacroDOF.begin(), interpolatedMacroDOF.begin() + dim );

        //The interpolated micro-displacement ( phi )
        floatVector phi( interpolatedMacroDOF.begin() + dim, interpolatedMacroDOF.begin() + dim + dim * dim );

        //Compute the micro displacements
        errorOut error = addMacroDomainDisplacementToMicro( dim, domainMicroNodeIndices, u, phi, domainReferenceXis,
                                                            domainMicroWeights, microDisplacements );

        if ( error ){
            errorOut result = new errorNode( "projectMacroDomainDisplacementToMicro",
                                             "Error in projection of the macro-displacements to the micro-scale" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut addMacroDomainDisplacementToMicro( const unsigned int dim, const uIntVector &domainMicroNodeIndices,
                                                const floatVector &u, const floatVector &phi,
                                                const floatVector &domainReferenceXis,
                                                const floatVector &domainMicroWeights, floatVector &microDisplacements ){
        /*!
         * Add the contribution of a macro domain's deformation to the micro-scale
         *
         * :param const unsigned int dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The global micro-node indices in the given domain
         * :param const floatVector &u: The macro-displacement at the local center of mass
         * :param const floatVector &phi: The micro-displacement ( in the micro-morphic sense ) at the local center of mass
         * :param const floatVector &domainReferenceXis: The vectors from the local center of mass to the micro-nodes in the domain.
         * :param const floatVector &domainMicroWeights: The weight associated with each micro-scale node for this domain. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &microDisplacements: The displacements of the micro-scale nodes as determined by the macro-scale
         *     projection.
         */

        if ( domainMicroWeights.size() != domainMicroNodeIndices.size() ){
            return new errorNode( "projectMacroDomainDisplacementToMicro",
                                  "The number of micro domain weights is not equal to the number of micro nodes in the domain" );
        }

        if ( domainReferenceXis.size() != dim * domainMicroNodeIndices.size() ){
            return new errorNode( "projectMacroDomainDisplacementToMicro",
                                  "The number of Xi vectors is not equal to the number of micro nodes in the domain" );
        }

        //The current micro node number
        unsigned int m;

        //The current reference micro-position vector
        floatVector Xi;

        //The current value of the micro-degrees of freedom
        floatVector q;

        for ( unsigned int i = 0; i < domainMicroNodeIndices.size(); i++ ){

            //Get the index of the micro-scale node
            m = domainMicroNodeIndices[ i ];

            if ( dim * ( m + 1 ) > microDisplacements.size() ){
                return new errorNode( "projectMacroDomainDisplacementToMicro",
                                      "The micro-displacements vector is too small for the micro-nodes" ); 
            }

            //Get the value of the micro-position
            Xi = floatVector( domainReferenceXis.begin() + dim * i, domainReferenceXis.begin() + dim * i + dim );

            //Compute the value of the micro-node's displacement
            q = u + vectorTools::matrixMultiply( phi, Xi, dim, dim, dim, 1 );

            //Add the contribution to the displacement
            microDisplacements[ dim * m + 0 ] += domainMicroWeights[ i ] * q[ 0 ]; 
            microDisplacements[ dim * m + 1 ] += domainMicroWeights[ i ] * q[ 1 ]; 
            microDisplacements[ dim * m + 2 ] += domainMicroWeights[ i ] * q[ 2 ]; 

        }

        return NULL;
    }

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
                                                        const floatVector &domainMicroWeights, SparseMatrix &domainN ){
        /*!
         * Construct the interpolation matrix for a macro domain overlapping with a 
         * micro-domain
         *
         * Note: It is assumed, that both the macro and micro-scale's have the same dimension. The number of spatial 
         *       variables at the micro-scale is therefore three and the number of macro-scale spatial parameters is
         *       dim + dim * dim.
         *
         * :param const unsigned int &dim: The dimension of the problem.
         * :param const unsigned int &nMicroNodes: The number of micro-nodes in the N matrix
         * :param const unsigned int &nMacroNodes: The number of macro-nodes in the N matrix
         * :param const uIntVector &domainMicroNodeIndices: The global micro-node indices in the given domain
         * :param const uIntVector &domainMacroDOFIndices: The indices of the macro-scale nodes present in the domain
         *     these are ( in a macro-scale FEA context ) the nodes of the micromorphic finite element.
         * :param const floatVector &domainReferenceXis: The vectors from the local center of mass to the micro-nodes in the domain.
         * :param const floatVector &domainMacroInterpolationFunctionValues: The values of the interpolation functions
         *     at the local center of mass.
         * :param const floatVector &domainMicroWeights: The weight associated with each micro-scale node for this domain. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param SparseMatrix &domainN: The macro domain's interpolation function.
         */

        //Error handling
        if ( dim != 3 ){
            return new errorNode( "formMacroDomainToMicroInterpolationMatrix",
                                  "Only 3D domains are currently supported" );
        }

        if ( dim * domainMicroNodeIndices.size() != domainReferenceXis.size() ){
            return new errorNode( "formMacroDomainToMicroInterpolationMatrix",
                                  "The number of micro node indices is not equal to the number of Xi vectors" );
        }

        if ( domainMicroNodeIndices.size() != domainMicroWeights.size() ){
            return new errorNode( "formMacroDomainToMicroInterpolationMatrix",
                                  "The number of micro node indices is not equal to the number of weights" );
        }

        if ( domainMacroNodeIndices.size() != domainMacroInterpolationFunctionValues.size() ){
            return new errorNode( "formMacroDomainToMicroInterpolationMatrix",
                                  "The number of macro indices is not equal to the number of macro interpolation function values" );
        }

        //Set the number of spatial degrees of freedom for the two scales
        const unsigned int nMicroDOF = dim;
        const unsigned int nMacroDOF = dim + dim * dim;

        //Set up the vector of terms for the sparse matrix
        std::vector< T > coefficients;
        coefficients.reserve( nMicroDOF * domainMicroNodeIndices.size() * nMacroDOF * domainMacroNodeIndices.size() );

        //Initialize the row and column indices for the sparse matrix
        unsigned int row0, col0;

        //Initialize the Xi vector
        floatVector Xi;

        //Initialize the value of the weight and shapefunction value
        floatType w, sf;

        for ( unsigned int i = 0; i < domainMicroNodeIndices.size(); i++ ){

            //Set the row index
            row0 = nMicroDOF * domainMicroNodeIndices[ i ];

            //Set the micro-position vector
            Xi = floatVector( domainReferenceXis.begin() + nMicroDOF * i,
                              domainReferenceXis.begin() + nMicroDOF * i + nMicroDOF );

            //Set the value of the weight
            w = domainMicroWeights[ i ];

            for ( unsigned int j = 0; j < domainMacroNodeIndices.size(); j++ ){

                //Set the column index
                col0 = nMacroDOF * domainMacroNodeIndices[ j ];

                //Set the shape function value
                sf = domainMacroInterpolationFunctionValues[ j ];

                //Set the coefficients of the shape function matrix
                coefficients.push_back( T( row0 + 0, col0 +  0, w * sf ) );
                coefficients.push_back( T( row0 + 1, col0 +  1, w * sf ) );
                coefficients.push_back( T( row0 + 2, col0 +  2, w * sf ) );
                coefficients.push_back( T( row0 + 0, col0 +  3, w * sf * Xi[ 0 ] ) );
                coefficients.push_back( T( row0 + 0, col0 +  4, w * sf * Xi[ 1 ] ) );
                coefficients.push_back( T( row0 + 0, col0 +  5, w * sf * Xi[ 2 ] ) );
                coefficients.push_back( T( row0 + 1, col0 +  6, w * sf * Xi[ 0 ] ) );
                coefficients.push_back( T( row0 + 1, col0 +  7, w * sf * Xi[ 1 ] ) );
                coefficients.push_back( T( row0 + 1, col0 +  8, w * sf * Xi[ 2 ] ) );
                coefficients.push_back( T( row0 + 2, col0 +  9, w * sf * Xi[ 0 ] ) );
                coefficients.push_back( T( row0 + 2, col0 + 10, w * sf * Xi[ 1 ] ) );
                coefficients.push_back( T( row0 + 2, col0 + 11, w * sf * Xi[ 2 ] ) );

            }

        }

        //Assemble the sparse matrix
        domainN = SparseMatrix( nMicroDOF * nMicroNodes, nMacroDOF * nMacroNodes );
        domainN.setFromTriplets( coefficients.begin(), coefficients.end() );

        return NULL;
    }

    errorOut addDomainMicroContributionToMacroMass( const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                    const floatVector &microMasses, const floatVector &domainMicroShapeFunctions,
                                                    const floatVector &domainMicroWeights,
                                                    floatVector &projectedMicroMasses ){
        /*!
         * Add the contribution of the micro-nodes' mass to the macro nodes.
         *
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes in the domain
         * :param const floatVector &microMasses: The masses of the micro-nodes
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &domainMicroWeights: The weight associated with each micro-scale node for this domain. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &projectedMicroMasses: The projected micro-masses on all of the macro-scale nodes.
         */

        //Error handling
        if ( domainMacroNodeIndices.size() * domainMicroNodeIndices.size() != domainMicroShapeFunctions.size() ){
            return new errorNode( "addDomainMicroContributionToMacroMass",
                                  "The size of the domain node indices vectors are not consistent with the number of shape functions" );
        }

        if ( domainMicroWeights.size() != domainMicroNodeIndices.size() ){
            return new errorNode( "addDomainMicroContributionToMacroMass", "The size of the domain's micro weights vector is not equal to the number of domain micro node indices" );
        }

        for ( unsigned int i = 0; i < domainMacroNodeIndices.size(); i++ ){
            if ( domainMacroNodeIndices[ i ] >= projectedMicroMasses.size() ){
                return new errorNode( "addDomainMicroContributionToMacroMass",
                                      "The size of the projected micro mass vector is smaller than a macro node requires" );
            }
        }

        //Set the number of macro-nodes in this domain
        unsigned int nMacroNodes = domainMacroNodeIndices.size();

        //Initialize the micro node mass
        floatType mass;

        //Initialize the micro node weights
        floatType weight;

        //Iterate over the micro-node indices
        for ( unsigned int i = 0; i < domainMicroNodeIndices.size(); i++ ){

            //Error handling
            if ( domainMicroNodeIndices[ i ] >= microMasses.size() ){

                return new errorNode( "addDomainMicroContributionToMacroMass",
                                      "The micro-node index is too large for the provided mass and shape function vectors" );

            }

            //Get the micro-node's mass
            mass = microMasses[ domainMicroNodeIndices[ i ] ];

            //Get the micro-node's weight
            weight = domainMicroWeights[ i ];

            //Iterate over the macro-node indices
            for ( unsigned int j = 0; j < domainMacroNodeIndices.size(); j++ ){

                //Add the micro-node contribution to the macro node
                projectedMicroMasses[ domainMacroNodeIndices[ j ] ] += mass * domainMicroShapeFunctions[ i * nMacroNodes + j ] * weight;

            }

        }

        return NULL;
    }
 
    errorOut addDomainMicroContributionToMacroMicroMassMomentOfInertia( const unsigned int &dim,
                                                                        const uIntVector &domainMicroNodeIndices,
                                                                        const uIntVector &domainMacroNodeIndices,
                                                                        const floatVector &domainReferenceXis,
                                                                        const floatVector &microMasses,
                                                                        const floatVector &domainMicroShapeFunctions,
                                                                        const floatVector &domainMicroWeights,
                                                                        floatVector &projectedMassMicroMomentOfInertia ){
        /*!
         * Add the contribution of the micro-nodes in the domain to the macro moment of inertia.
         *
         * :param const unsigned int &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &domainReferenceXis: The micro-position vectors in the domain
         * :param const floatVector &microMasses: The masses of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &domainMicroWeights: The weight associated with each micro-scale node for this domain. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &projectedMassMicroMomentOfInertia: The moments of inertia at the macro-nodes of the domain
         *     as projected from the micro-nodes weighted by the mass.
         */

        floatVector projectedMassConstant;
        floatVector projectedMassDisplacement;
        floatVector projectedMassDisplacementPosition;
        floatVector microDisplacements;
        return addDomainMicroToMacroProjectionTerms( dim, domainMicroNodeIndices, domainMacroNodeIndices,
                                                     domainReferenceXis, microMasses, domainMicroShapeFunctions,
                                                     domainMicroWeights, microDisplacements,
                                                     projectedMassMicroMomentOfInertia,
                                                     projectedMassConstant,
                                                     projectedMassDisplacement,
                                                     projectedMassDisplacementPosition,
                                                     true, false, false, false );

    }

    errorOut addDomainMassConstant( const unsigned int &dim,
                                    const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                    const floatVector &domainReferenceXis, const floatVector &microMasses,
                                    const floatVector &domainMicroShapeFunctions, const floatVector &domainMicroWeights,
                                    floatVector &projectedMassConstant ){
        /*!
         * Add the contributions of the domain to the mass constant
         *
         * :param const unsigned int &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &domainReferenceXis: The micro-position vectors in the domain
         * :param const floatVector &microMasses: The masses of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &domainMicroWeights: The weight associated with each micro-scale node for this domain. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &projectedMassConstant: The projected mass constant at the macro-scale node.
         */

        floatVector projectedMassMicroMomentOfInertia;
        floatVector projectedMassDisplacement;
        floatVector projectedMassDisplacementPosition;
        floatVector microDisplacements;
        return addDomainMicroToMacroProjectionTerms( dim, domainMicroNodeIndices, domainMacroNodeIndices,
                                                     domainReferenceXis, microMasses, domainMicroShapeFunctions,
                                                     domainMicroWeights, microDisplacements,
                                                     projectedMassMicroMomentOfInertia,
                                                     projectedMassConstant,
                                                     projectedMassDisplacement,
                                                     projectedMassDisplacementPosition,
                                                     false, true, false, false );

    }

    errorOut addDomainMassDisplacement( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                        const floatVector &microMasses, const floatVector &domainMicroShapeFunctions,
                                        const floatVector &domainMicroWeights, const floatVector &microDisplacements,
                                        floatVector &projectedMassDisplacement ){
        /*!
         * Add the contributions of the domain to the mass constant
         *
         * :param const unsigned int &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &microMasses: The masses of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &domainMicroWeights: The weight associated with each micro-scale node for this domain. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &microDisplacements: The displacements of the micro-degrees of freedom.
         * :param floatVector &projectedMassDisplacement: The projected mass weighted displacement at the macro-scale node.
         */

        floatVector projectedMassMicroMomentOfInertia;
        floatVector projectedMassConstant;
        floatVector projectedMassDisplacementPosition;
        floatVector domainReferenceXis;
        return addDomainMicroToMacroProjectionTerms( dim, domainMicroNodeIndices, domainMacroNodeIndices,
                                                     domainReferenceXis, microMasses, domainMicroShapeFunctions,
                                                     domainMicroWeights, microDisplacements,
                                                     projectedMassMicroMomentOfInertia,
                                                     projectedMassConstant,
                                                     projectedMassDisplacement,
                                                     projectedMassDisplacementPosition,
                                                     false, false, true, false );

    }

    errorOut addDomainMassMicroDisplacementPosition( const unsigned int &dim,
                                                     const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                     const floatVector &domainReferenceXis, const floatVector &microMasses,
                                                     const floatVector &domainMicroShapeFunctions, const floatVector &domainMicroWeights,
                                                     const floatVector &microDisplacements,
                                                     floatVector &projectedMassDisplacementPosition ){
        /*!
         * Add the contributions of the domain to the mass constant
         *
         * :param const unsigned int &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &domainReferenceXis: The micro-position vectors in the domain
         * :param const floatVector &microMasses: The masses of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &domainMicroWeights: The weight associated with each micro-scale node for this domain. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &microDisplacements: The displacements of the micro-degrees of freedom.
         * :param floatVector &projectedMassDisplacementPosition: The projected mass weighted dyadic product of displacement
         *     and the micro-position at the macro-scale node.
         */

        floatVector projectedMassMicroMomentOfInertia;
        floatVector projectedMassConstant;
        floatVector projectedMassDisplacement;
        return addDomainMicroToMacroProjectionTerms( dim, domainMicroNodeIndices, domainMacroNodeIndices,
                                                     domainReferenceXis, microMasses, domainMicroShapeFunctions,
                                                     domainMicroWeights, microDisplacements,
                                                     projectedMassMicroMomentOfInertia,
                                                     projectedMassConstant,
                                                     projectedMassDisplacement,
                                                     projectedMassDisplacementPosition,
                                                     false, false, false, true );

    }


    errorOut addDomainMicroToMacroProjectionTerms( const unsigned int &dim,
                                                   const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                   const floatVector &domainReferenceXis, const floatVector &microMasses,
                                                   const floatVector &domainMicroShapeFunctions, const floatVector &domainMicroWeights,
                                                   const floatVector &microDisplacements,
                                                   floatVector &projectedMassMicroMomentOfInertia,
                                                   floatVector &projectedMassConstant,
                                                   floatVector &projectedMassDisplacement,
                                                   floatVector &projectedMassDisplacementPosition,
                                                   const bool computeMassMomentOfInertia,
                                                   const bool computeMassConstant,
                                                   const bool computeMassMicroDisplacement,
                                                   const bool computeMassDisplacementPosition ){
        /*!
         * Solve for the terms required to project from the micro-scale to the macro-scale.
         *
         * :param const unsigned int &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &domainReferenceXis: The micro-position vectors in the domain
         * :param const floatVector &microMasses: The masses of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &domainMicroWeights: The weight associated with each micro-scale node for this domain. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &microDisplacements: The displacements of the micro-degrees of freedom.
         * :param floatVector &projectedMassMicroMomentOfInertia: The moments of inertia at the macro-nodes of the domain
         *     as projected from the micro-nodes weighted by the mass.
         * :param floatVector &projectedMassDisplacementPosition: The mass weighted dyadic product of the micro displacement
         *     and the micro position projected to the macro nodes.
         * :param const bool computeMassMomentOfInertia: Boolean for whether the contribution to the mass-weighted moment 
         *     of inertia should be computed.
         * :param const bool computeMassConstant: Boolean for whether the contribution to the mass constant.
         * :param const bool computeMassMicroDisplacement: Boolean for whether the contribution of the mass-weighted
         *     displacement should be calculated.
         */

        //Error handling
        if ( ( computeMassMomentOfInertia ) || ( computeMassConstant ) || ( computeMassDisplacementPosition ) ){
            if ( dim * domainMicroNodeIndices.size() != domainReferenceXis.size() ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The number of micro node indices and the micro position vectors do not have consistent sizes" );
            }
        }

        if ( domainMicroNodeIndices.size() != domainMicroWeights.size() ){
            return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                  "The number of micro node indices and the micro weights are not the same" );
        }

        if ( domainMicroNodeIndices.size() * domainMacroNodeIndices.size() != domainMicroShapeFunctions.size() ){
            return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                  "The number of micro and micro node indices are not consistent with the number of shape functions" );
        }

        for ( unsigned int i = 0; i < domainMacroNodeIndices.size(); i++ ){
            if ( ( projectedMassMicroMomentOfInertia.size() < dim * dim * ( domainMacroNodeIndices[ i ] + 1 ) ) &&
                 ( computeMassMomentOfInertia ) ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The size of the projected micro moment of inertia weighted by the mass is smaller than required for the provided nodes" );
            }

            if ( ( projectedMassConstant.size() < dim * ( domainMacroNodeIndices[ i ] + 1 ) ) &&
                 ( computeMassConstant ) ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The size of the projected mass constant is smaller than required for the provided nodes" );
            }

            if ( ( projectedMassDisplacement.size() < dim * ( domainMacroNodeIndices[ i ] + 1 ) ) &&
                 ( computeMassMicroDisplacement ) ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The size of the projected mass-weighted micro displacement is smaller than required for the provided nodes" );
            }

            if ( ( projectedMassDisplacementPosition.size() < dim * dim * ( domainMacroNodeIndices[ i ] + 1 ) ) &&
                 ( computeMassDisplacementPosition ) ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The size of the projected mass-weighted dyadic product of the micro displacement and the micro position is smaller than required for the provided nodes" );
            }

        }

        if ( ( computeMassMicroDisplacement ) || ( computeMassDisplacementPosition ) ){
            for ( unsigned int i = 0; i < domainMicroNodeIndices.size(); i++ ){
                if ( microDisplacements.size() < dim * ( domainMicroNodeIndices[ i ] + 1 ) ){
                    return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                          "The size of the micro degree of freedom vector is too small for the provided nodes" );
                }
            }
        }

        //Initialize the micro-mass, weight, and shape function value
        floatType mass, weight, sf;

        //Initialize the micro-displacement vector
        floatVector q( dim );

        //Initialize the micro-position vector
        floatVector Xi( dim );

        //Initialize the dyadic product of the micro-position vector
        floatVector XiXi( dim * dim );

        //Loop through the micro nodes
        for ( unsigned int i = 0; i < domainMicroNodeIndices.size(); i++ ){

            //Extract the nodal mass
            mass = microMasses[ domainMicroNodeIndices[ i ] ];

            //Extract the nodal weight
            weight = domainMicroWeights[ i ];

            //Extract the current micro-position vectors
            if ( ( computeMassConstant ) || ( computeMassDisplacementPosition ) ){

                Xi = floatVector( domainReferenceXis.begin() + dim * i, domainReferenceXis.begin() + dim * ( i + 1 ) );

            }

            //Compute the dyadic product of Xi with itself
            if ( computeMassMomentOfInertia ){

                for ( unsigned int j = 0; j < dim; j++ ){
    
                    for ( unsigned int k = 0; k < dim; k++ ){
    
                        XiXi[ dim * j + k ] = domainReferenceXis[ dim * i + j ] * domainReferenceXis[ dim * i + k ];
    
                    }
    
                }
            }

            if ( ( computeMassMicroDisplacement ) || ( computeMassDisplacementPosition ) ){

                q = floatVector( microDisplacements.begin() + dim * domainMicroNodeIndices[ i ],
                                 microDisplacements.begin() + dim * domainMicroNodeIndices[ i ] + dim );

            }

            //Loop through the macro nodes
            for ( unsigned int j = 0; j < domainMacroNodeIndices.size(); j++ ){

                //Set the shape-function value
                sf = domainMicroShapeFunctions[ domainMacroNodeIndices.size() * i + j ];
               
                if ( computeMassMomentOfInertia ){ 

                    for ( unsigned int k = 0; k < dim * dim; k++ ){
    
                        //Add the contribution to the micro-moment of inertia
                        projectedMassMicroMomentOfInertia[ dim * dim * domainMacroNodeIndices[ j ] + k ]
                            += weight * mass * sf * XiXi[ k ];
    
                    }

                }

                if ( computeMassConstant ){ 

                    for ( unsigned int k = 0; k < dim; k++ ){
    
                        //Add the contribution to the micro-moment of inertia
                        projectedMassConstant[ dim * domainMacroNodeIndices[ j ] + k ] += weight * mass * sf * Xi[ k ];
    
                    }

                }

                if ( computeMassMicroDisplacement ) {

                    for ( unsigned int k = 0; k < dim; k++ ){

                        //Add the contribution to the mass weighted micro-displacement
                        projectedMassDisplacement[ dim * domainMacroNodeIndices[ j ] + k ] += weight * mass * sf * q[ k ];

                    }

                }

                if ( computeMassDisplacementPosition ){

                    for ( unsigned int k = 0; k < dim; k++ ){

                        for ( unsigned int l = 0; l < dim; l++ ){
                            projectedMassDisplacementPosition[ dim * dim * domainMacroNodeIndices[ j ] + dim * k + l ]
                                += weight * mass * sf * q[ k ] * Xi[ l ];
                        }

                    }

                }

            }

        }

        return NULL;

    }

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microMasses,
                                        const floatVector &microPositions, const floatVector &domainMicroWeights,
                                        floatVector &domainCM ){
        /*!
         * Compute the center of mass of a micro domain from the masses of the micro-nodes contained
         * within the domain.
         *
         * :param const unsigned int &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const floatVector &microMasses: The masses of the micro nodes.
         * :param const floatVector &microPositions: The positions of the micro-nodes.
         * :param const floatVector &domainMicroWeights: The weight associated with each micro-scale node for this domain. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &domainCG: The center of mass of the domain
         */

        floatType domainMass;
        return computeDomainCenterOfMass( dim, domainMicroNodeIndices, microMasses, microPositions, domainMicroWeights,
                                          domainMass, domainCM );
    }

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microVolumes,
                                        const floatVector &microDensities, const floatVector &microPositions,
                                        const floatVector &domainMicroWeights, floatVector &domainCM ){
        /*!
         * Compute the center of mass of a micro domain from the masses of the micro-nodes contained
         * within the domain.
         *
         * :param const unsigned int &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const floatVector &microVolumes: The volumes of the micro nodes.
         * :param const floatVector &microDensities: The densities of the micro nodes.
         * :param const floatVector &microPositions: The positions of the micro-nodes.
         * :param const floatVector &domainMicroWeights: The weight associated with each micro-scale node for this domain. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &domainCG: The center of mass of the domain
         */

        floatType domainMass;
        return computeDomainCenterOfMass( dim, domainMicroNodeIndices, microVolumes, microDensities, microPositions,
                                          domainMicroWeights, domainMass, domainCM );
    }

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microMasses,
                                        const floatVector &microPositions, const floatVector &domainMicroWeights,
                                        floatType &domainMass, floatVector &domainCM ){
        /*!
         * Compute the center of mass of a micro domain from the masses of the micro-nodes contained
         * within the domain.
         *
         * :param const unsigned int &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const floatVector &microMasses: The masses of the micro nodes.
         * :param const floatVector &microPositions: The positions of the micro-nodes.
         * :param const floatVector &domainMicroWeights: The weight associated with each micro-scale node for this domain. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &domainMass: The mass of the domain
         * :param floatVector &domainCM: The center of mass of the domain
         */

        if ( domainMicroNodeIndices.size() != domainMicroWeights.size() ){
            return new errorNode( "computeDomainCenterOfMass",
                                  "The number of node indices and the weights must have the same length" );
        }

        for ( unsigned int i = 0; i < domainMicroNodeIndices.size(); i++ ){
            if ( microPositions.size() < dim * domainMicroNodeIndices[ i ] + dim ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-positions vector is not consistent with the micro indices" );
            }

            if ( microMasses.size() < domainMicroNodeIndices[ i ] ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-masses vector is not consistent with the micro indices" );
            }
        }

        //Initialize the domain mass
        domainMass = 0;

        //Initialize the center of mass vector
        domainCM = floatVector( dim, 0 );

        for ( unsigned int i = 0; i < domainMicroNodeIndices.size(); i++ ){


            //Add to the domain's mass
            domainMass += microMasses[ domainMicroNodeIndices[ i ] ] * domainMicroWeights[ i ];

            //Add to the domain's mass weighted position
            domainCM += microMasses[ domainMicroNodeIndices[ i ] ] * domainMicroWeights[ i ]
                      * floatVector( microPositions.begin() + dim * domainMicroNodeIndices[ i ],
                                     microPositions.begin() + dim * ( domainMicroNodeIndices[ i ] + 1 ) );

        }

        //Normalize the center of mass by the domain's mass
        domainCM /= domainMass;

        return NULL;

    }

    errorOut computeDomainCenterOfMass( const unsigned int &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microVolumes,
                                        const floatVector &microDensities, const floatVector &microPositions,
                                        const floatVector &domainMicroWeights, floatType &domainMass, floatVector &domainCM ){
        /*!
         * Compute the center of mass of a micro domain from the masses of the micro-nodes contained
         * within the domain.
         *
         * :param const unsigned int &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const floatVector &microVolumes: The volumes of the micro nodes.
         * :param const floatVector &microDensities: The densities of the micro nodes.
         * :param const floatVector &microPositions: The positions of the micro-nodes.
         * :param const floatVector &domainMicroWeights: The weight associated with each micro-scale node for this domain. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &domainMass: The mass of the domain
         * :param floatVector &domainCM: The center of mass of the domain
         */

        if ( domainMicroNodeIndices.size() != domainMicroWeights.size() ){
            return new errorNode( "computeDomainCenterOfMass",
                                  "The number of node indices and the weights must have the same length" );
        }

        for ( unsigned int i = 0; i < domainMicroNodeIndices.size(); i++ ){
            if ( microPositions.size() < dim * domainMicroNodeIndices[ i ] + dim ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-positions vector is not consistent with the micro indices" );
            }

            if ( microVolumes.size() < domainMicroNodeIndices[ i ] ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-volumes vector is not consistent with the micro indices" );
            }

            if ( microDensities.size() < domainMicroNodeIndices[ i ] ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-densities vector is not consistent with the micro indices" );
            }
        }

        //Initialize the domain mass
        domainMass = 0;

        //Initialize the center of mass vector
        domainCM = floatVector( dim, 0 );

        for ( unsigned int i = 0; i < domainMicroNodeIndices.size(); i++ ){


            //Add to the domain's mass
            domainMass += microVolumes[ domainMicroNodeIndices[ i ] ] * microDensities[ domainMicroNodeIndices[ i ] ]
                        * domainMicroWeights[ i ];

            //Add to the domain's mass weighted position
            domainCM += microVolumes[ domainMicroNodeIndices[ i ] ] * microDensities[ domainMicroNodeIndices[ i ] ]
                      * domainMicroWeights[ i ]
                      * floatVector( microPositions.begin() + dim * domainMicroNodeIndices[ i ],
                                     microPositions.begin() + dim * ( domainMicroNodeIndices[ i ] + 1 ) );

        }

        //Normalize the center of mass by the domain's mass
        domainCM /= domainMass;

        return NULL;

    }

    errorOut computeDomainXis( const unsigned int &dim,
                               const uIntVector &domainMicroNodeIndices, const floatVector &microPositions,
                               const floatVector &domainCM, floatVector &domainXis ){
        /*
         * Compute the relative position vector between the center of mass of a micro domain and the 
         * micro position.
         *
         * :param const unsigned int &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const floatVector &microPositions: The positions of the micro-nodes.
         * :param floatVector &domainCM: The center of mass of the domain
         * :param floatVector &domainXis: The relative positions of the micro nodes.
         */

        //Error Handling
        if ( domainCM.size() != dim ){
            return new errorNode( "computeDomainXis",
                                  "The center of mass is not consistent with the dimension" );
        }

        for ( unsigned int i = 0; i < domainMicroNodeIndices.size(); i++ ){
            if ( microPositions.size() < dim * domainMicroNodeIndices[ i ] + dim ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-positions vector is not consistent with the micro indices" );
            }
        }

        //Resize the Xi vector
        domainXis.resize( dim * domainMicroNodeIndices.size() );

        for ( unsigned int i = 0; i < domainMicroNodeIndices.size(); i++ ){

            for ( unsigned int j = 0; j < dim; j++ ){
           
                //Compute the relative position vector 
                domainXis[ dim * i + j ] = microPositions[ dim * domainMicroNodeIndices[ i ] + j ] - domainCM[ j ];

            }

        }

        return NULL;

    }

}
