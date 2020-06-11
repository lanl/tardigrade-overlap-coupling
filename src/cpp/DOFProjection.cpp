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
                                                const floatVector &referenceXis,
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
         * :param const floatVector &referenceXis: The Xi vectors of the micro-scale nodes which go from the local
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
        errorOut error = addMacroDomainDisplacementToMicro( dim, domainMicroNodeIndices, u, phi, referenceXis,
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
                                                const floatVector &referenceXis,
                                                const floatVector &domainMicroWeights, floatVector &microDisplacements ){
        /*!
         * Add the contribution of a macro domain's deformation to the micro-scale
         *
         * :param const unsigned int dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The global micro-node indices in the given domain
         * :param const floatVector &u: The macro-displacement at the local center of mass
         * :param const floatVector &phi: The micro-displacement ( in the micro-morphic sense ) at the local center of mass
         * :param const floatVector &referenceXis: The vectors from the local center of mass to the micro-nodes in the domain.
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

        if ( referenceXis.size() != dim * domainMicroNodeIndices.size() ){
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
            Xi = floatVector( referenceXis.begin() + dim * i, referenceXis.begin() + dim * i + dim );

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
                                                        const floatVector &referenceXis,
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
         * :param const floatVector &referenceXis: The vectors from the local center of mass to the micro-nodes in the domain.
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

        if ( dim * domainMicroNodeIndices.size() != referenceXis.size() ){
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
            Xi = floatVector( referenceXis.begin() + nMicroDOF * i,
                              referenceXis.begin() + nMicroDOF * i + nMicroDOF );

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
         * :param floatVector &projectedMicroMasses: The projected micro-masses on all of the macro-scale nodes.
         */

        //Error handling
        if ( domainMacroNodeIndices.size() * domainMicroNodeIndices.size() != domainMicroShapeFunctions.size() ){
            return new errorNode( "addDomainMicroContributionToMacroMass",
                                  "The size of the domain node indices vectors are not consistent with the number of shape functions" );
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

        //Iterate over the micro-node indices
        for ( unsigned int i = 0; i < domainMicroNodeIndices.size(); i++ ){

            //Error handling
            if ( domainMicroNodeIndices[ i ] >= microMasses.size() ){
                return new errorNode( "addDomainMicroContributionToMacroMass",
                                      "The micro-node index is too large for the provided mass and shape function vectors" );
            }

            //Get the index of the micro node
            mass = microMasses[ domainMicroNodeIndices[ i ] ];

            //Iterate over the macro-node indices
            for ( unsigned int j = 0; j < domainMacroNodeIndices.size(); j++ ){

                //Add the micro-node contribution to the macro node
                projectedMicroMasses[ domainMacroNodeIndices[ j ] ] += mass * domainMicroShapeFunctions[ i * nMacroNodes ];

            }

        }

        return NULL;
    }
}
