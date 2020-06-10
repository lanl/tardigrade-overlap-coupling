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

    errorOut addMacroDomainDisplacementToMicro( const unsigned int dim,
                                                const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                const floatVector &referenceXis,
                                                const floatVector &domainMacroInterpolationFunctionValues,
                                                const unsigned int &nMacroDOF, const floatVector &MacroDOFVector,
                                                const floatVector &domainMicroWeights, floatVector &microDisplacements ){

        /*!
         * Add the contribution of a macro domain's deformation to the micro-scale
         *
         * :param const unsigned int dim: The dimension of the problem ( only 3 is tested )
         * :param const uIntVector &domainMicroDOFIndices: The indices of the micro-scale nodes present in the domain
         * :param const uIntVector &domainMacroDOFIndices: THe indices of the macro-scale nodes present in the domain
         *     these are ( in a macro-scale FEA context ) the nodes of the micromorphic finite element.
         * :param const floatVector &referenceXis: The Xi vectors of the micro-scale nodes which go from the local
         *     center of mass to the position of the micro-scale node in the reference configuration.
         * :param const floatVector &domainMacroInterpolationFunctionValues: The values of the interpolation functions
         *     at the local center of mass.
         * :param const unsigned int &nMacroDOF: The number of degrees of freedom associated with each macro-scale node.
         *     ( only 12 is tested )
         * :param const floatVector &MacroDOFVector: The degree of freedom vector for the macro-scale. The ordering for each node
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
                                  * floatVector( MacroDOFVector.begin() + nMacroDOF * domainMacroNodeIndices[ i ],
                                                 MacroDOFVector.begin() + nMacroDOF * domainMacroNodeIndices[ i ] + nMacroDOF );
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
         *       if the minimum L2 projection is being used.
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
}
