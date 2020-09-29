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

    errorOut addMacroDomainDisplacementToMicro( const uIntType dim,
                                                const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                const floatVector &domainReferenceXis,
                                                const floatVector &domainMacroInterpolationFunctionValues,
                                                const uIntType &nMacroDOF, const floatVector &macroDOFVector,
                                                const floatVector &microWeights, floatVector &microDisplacements,
                                                const std::unordered_map< uIntType, uIntType > *microNodeToLocalIndex ){

        /*!
         * Add the contribution of a macro domain's deformation to the micro-scale
         *
         * :param const uIntType dim: The dimension of the problem ( only 3 is tested )
         * :param const uIntVector &domainMicroDOFIndices: The indices of the micro-scale nodes present in the domain
         * :param const uIntVector &domainMacroDOFIndices: The indices of the macro-scale nodes present in the domain
         *     these are ( in a macro-scale FEA context ) the nodes of the micromorphic finite element.
         * :param const floatVector &domainReferenceXis: The Xi vectors of the micro-scale nodes which go from the local
         *     center of mass to the position of the micro-scale node in the reference configuration.
         * :param const floatVector &domainMacroInterpolationFunctionValues: The values of the interpolation functions
         *     at the local center of mass.
         * :param const uIntType &nMacroDOF: The number of degrees of freedom associated with each macro-scale node.
         *     ( only 12 is tested )
         * :param const floatVector &macroDOFVector: The degree of freedom vector for the macro-scale. The ordering for each node
         *     is assumed to be [ u1, u2, u3, phi11, phi12, phi13, phi21, phi22, phi23, phi31, phi32, phi33, ... ]
         * :param const floatVector &microWeights: The weight associated with each micro-scale node. 
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 projection is being used.
         * :param floatVector &microDisplacements: The displacements of the micro-scale nodes as determined by the macro-scale
         *     projection.
         * :param std::unordered_map< uIntType, uIntType > *microNodeToLocalIndex: A map from the micro node index 
         *     to the indices to be used in the output vector. The micro nodes which are either free or ghost may not be all
         *     of the micro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the macro index values
         *     will be used.
         *
         *     Note: If a node is not located in microNodeToLocalIndex ( and microNodeToLocalIndex is not NULL ) this node
         *           will be skipped.
         */

        if ( domainMacroNodeIndices.size() != domainMacroInterpolationFunctionValues.size() ){
            return new errorNode( "addMacroDomainDisplacementToMicro",
                                  "The macro-scale node indices and the macro-scale interpolation function values must be the same length" );
        }

        //The interpolated degree of freedom vector
        floatVector interpolatedMacroDOF( nMacroDOF, 0 );

        //Compute the interpolated value of u and phi
        for ( uIntType i = 0; i < domainMacroNodeIndices.size(); i++ ){
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
                                                            microWeights, microDisplacements, microNodeToLocalIndex );

        if ( error ){
            errorOut result = new errorNode( "addMacroDomainDisplacementToMicro",
                                             "Error in projection of the macro-displacements to the micro-scale" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut addMacroDomainDisplacementToMicro( const uIntType dim, const uIntVector &domainMicroNodeIndices,
                                                const floatVector &u, const floatVector &phi,
                                                const floatVector &domainReferenceXis,
                                                const floatVector &microWeights, floatVector &microDisplacements,
                                                const std::unordered_map< uIntType, uIntType > *microNodeToLocalIndex ){
        /*!
         * Add the contribution of a macro domain's deformation to the micro-scale
         *
         * :param const uIntType dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The global micro-node indices in the given domain
         * :param const floatVector &u: The macro-displacement at the local center of mass
         * :param const floatVector &phi: The micro-displacement ( in the micro-morphic sense ) at the local center of mass
         * :param const floatVector &domainReferenceXis: The vectors from the local center of mass to the micro-nodes in the domain.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &microDisplacements: The displacements of the micro-scale nodes as determined by the macro-scale
         *     projection.
         * :param std::unordered_map< uIntType, uIntType > *microNodeToLocalIndex: A map from the micro node index 
         *     to the indices to be used in the output vector. The micro nodes which are either free or ghost may not be all
         *     of the micro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the macro index values
         *     will be used.
         *
         *     Note: If a node is not located in microNodeToLocalIndex ( and microNodeToLocalIndex is not NULL ) this node
         *           will be skipped.
         */

        if ( !microNodeToLocalIndex ){
            if ( dim * microWeights.size() != microDisplacements.size() ){
                return new errorNode( "addMacroDomainDisplacementToMicro",
                                      "The number of micro domain weights is not consistent with the number of micro displacements" );
            }
        }

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){

            if ( microWeights.size() <= domainMicroNodeIndices[ i ] ){
                return new errorNode( "addMacroDomainDisplacementToMicro",
                                      "A micro-node index is larger than the micro-weights vector" );
            }

        }

        if ( domainReferenceXis.size() != dim * domainMicroNodeIndices.size() ){
            return new errorNode( "addMacroDomainDisplacementToMicro",
                                  "The number of Xi vectors is not equal to the number of micro nodes in the domain" );
        }

        //The current micro node number
        uIntType m;

        //The current output index
        uIntType o;

        //The current reference micro-position vector
        floatVector Xi;

        //The current value of the micro-degrees of freedom
        floatVector q;

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){

            //Get the index of the micro-scale node and the output index
            if ( microNodeToLocalIndex ){

                auto indx = microNodeToLocalIndex->find( domainMicroNodeIndices[ i ] );

                if ( indx == microNodeToLocalIndex->end( ) ){
                    continue;
//                    return new errorNode( "addMacroDomainDisplacementToMicro",
//                                          "Micro node index " + std::to_string( domainMicroNodeIndices[ i ] )
//                                          + " is not found in microNodeToLocalIndex" );
                }

                if ( indx->second >= microDisplacements.size( ) ){
                    return new errorNode( "addMacroDomainDisplacementToMicro",
                                          "Micro node index " + std::to_string( domainMicroNodeIndices[ i ] )
                                          + " is greater than the size of the micro displacements" );
                }

                o = indx->second;
                m = domainMicroNodeIndices[ i ];

            }
            else{
                o = domainMicroNodeIndices[ i ];
                m = domainMicroNodeIndices[ i ];
            }

            if ( dim * ( o + 1 ) > microDisplacements.size() ){
                return new errorNode( "addMacroDomainDisplacementToMicro",
                                      "The micro-displacements vector is too small for the micro-nodes" ); 
            }

            //Get the value of the micro-position
            Xi = floatVector( domainReferenceXis.begin() + dim * i, domainReferenceXis.begin() + dim * i + dim );

            //Compute the value of the micro-node's displacement
            q = u + vectorTools::matrixMultiply( phi, Xi, dim, dim, dim, 1 );

            //Add the contribution to the displacement
            microDisplacements[ dim * o + 0 ] += microWeights[ m ] * q[ 0 ]; 
            microDisplacements[ dim * o + 1 ] += microWeights[ m ] * q[ 1 ]; 
            microDisplacements[ dim * o + 2 ] += microWeights[ m ] * q[ 2 ]; 

        }

        return NULL;
    }

    /*===================================================
    |                Projection Matrices                |
    =====================================================
    | Functions which construct the projection matrices |
    ===================================================*/

    errorOut formMacroDomainToMicroInterpolationMatrix( const uIntType &dim,
                                                        const uIntType &nMicroNodes, const uIntType &nMacroNodes,
                                                        const uIntVector &domainMicroNodeIndices,
                                                        const uIntVector &domainMacroNodeIndices,
                                                        const floatVector &domainReferenceXis,
                                                        const floatVector &domainMacroInterpolationFunctionValues,
                                                        const floatVector &microWeights, SparseMatrix &domainN,
                                                        const std::unordered_map< uIntType, uIntType >* microNodeToLocalIndex,
                                                        const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex ){
        /*!
         * Construct the interpolation matrix for a macro domain overlapping with a 
         * micro-domain
         *
         * Note: It is assumed, that both the macro and micro-scale's have the same dimension. The number of spatial 
         *       variables at the micro-scale is therefore three and the number of macro-scale spatial parameters is
         *       dim + dim * dim.
         *
         * :param const uIntType &dim: The dimension of the problem.
         * :param const uIntType &nMicroNodes: The number of micro-nodes in the N matrix
         * :param const uIntType &nMacroNodes: The number of macro-nodes in the N matrix
         * :param const uIntVector &domainMicroNodeIndices: The global micro-node indices in the given domain
         * :param const uIntVector &domainMacroDOFIndices: The indices of the macro-scale nodes present in the domain
         *     these are ( in a macro-scale FEA context ) the nodes of the micromorphic finite element.
         * :param const floatVector &domainReferenceXis: The vectors from the local center of mass to the micro-nodes in the domain.
         * :param const floatVector &domainMacroInterpolationFunctionValues: The values of the interpolation functions
         *     at the local center of mass.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param SparseMatrix &domainN: The macro domain's interpolation function.
         * :param std::unordered_map< uIntType, uIntType > *microNodeToLocalIndex: A map from the micro node index 
         *     to the indices to be used in the output vector. The micro nodes which are either free or ghost may not be all
         *     of the micro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         *
         *     Note: If a node is not located in microNodeToLocalIndex ( and microNodeToLocalIndex is not NULL ) this node
         *           will be skipped.
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
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

        if ( !microNodeToLocalIndex ){
            if ( nMicroNodes != microWeights.size() ){
                return new errorNode( "formMacroDomainToMicroInterpolationMatrix",
                                      "The number of micro nodes is not equal to the number of weights" );
            }
        }

        if ( domainMacroNodeIndices.size() != domainMacroInterpolationFunctionValues.size() ){
            return new errorNode( "formMacroDomainToMicroInterpolationMatrix",
                                  "The number of macro indices is not equal to the number of macro interpolation function values" );
        }

        //Set the number of spatial degrees of freedom for the two scales
        const uIntType nMicroDOF = dim;
        const uIntType nMacroDOF = dim + dim * dim;

        //Set up the vector of terms for the sparse matrix
        std::vector< T > coefficients;
        coefficients.reserve( nMicroDOF * domainMicroNodeIndices.size() * nMacroDOF * domainMacroNodeIndices.size() );

        //Initialize the row and column indices for the sparse matrix
        uIntType row0, col0;

        //Initialize the Xi vector
        floatVector Xi;

        //Initialize the value of the weight and shapefunction value
        floatType w, sf;

        //Initialize the global micro node index
        uIntType m;

        //Initialize the global macro node index
        uIntType n;

        //Initialize the output index for the micro-nodes
        uIntType o;

        //Initialize the output index for the macro-nodes
        uIntType p;

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){

            //Set the global micro node index
            m = domainMicroNodeIndices[ i ];

            if ( m >= microWeights.size( ) ){
                return new errorNode( "formMacroDomainToMicroInterpolationMatrix",
                                      "The number of micro-weights is smaller than required for micro-node " + std::to_string( m ) );
            }

            //Set the row index
            if ( microNodeToLocalIndex ){

                auto indx = microNodeToLocalIndex->find( m );

                if ( indx == microNodeToLocalIndex->end( ) ){

                    continue;

                }

                o = indx->second;
            }
            else{
                o = domainMicroNodeIndices[ i ];
            }

            row0 = nMicroDOF * o;

            //Set the micro-position vector
            Xi = floatVector( domainReferenceXis.begin() + nMicroDOF * i,
                              domainReferenceXis.begin() + nMicroDOF * i + nMicroDOF );

            //Set the value of the weight
            w = microWeights[ m ];

            for ( uIntType j = 0; j < domainMacroNodeIndices.size(); j++ ){

                //Set the global macro node index
                n = domainMacroNodeIndices[ j ];

                //Set the column index
                if ( macroNodeToLocalIndex ){
    
                    auto indx = macroNodeToLocalIndex->find( n );
    
                    if ( indx == macroNodeToLocalIndex->end( ) ){
    
                        return new errorNode( "formMacroDomaintoMicroInterpolationMatrix",
                                              "The macro node " + std::to_string( n ) + " is not found in the macro node to local index map" );
    
                    }
    
                    p = indx->second;
    
                }
                else{
                    p = domainMacroNodeIndices[ j ];
                }

                col0 = nMacroDOF * p;

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

    errorOut formMacroDomainToMicroInterpolationMatrix( const uIntType &dim,
                                                        const uIntType &nMicroNodes, const uIntType &nMacroNodes,
                                                        const uIntVector &domainMicroNodeIndices,
                                                        const uIntVector &domainMacroNodeIndices,
                                                        const std::unordered_map< uIntType, floatVector > &domainReferenceXis,
                                                        const floatVector &domainMacroInterpolationFunctionValues,
                                                        const std::unordered_map< uIntType, floatType > &microWeights,
                                                        SparseMatrix &domainN,
                                                        const std::unordered_map< uIntType, uIntType >* microNodeToLocalIndex,
                                                        const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex ){
        /*!
         * Construct the interpolation matrix for a macro domain overlapping with a 
         * micro-domain
         *
         * Note: It is assumed, that both the macro and micro-scale's have the same dimension. The number of spatial 
         *       variables at the micro-scale is therefore three and the number of macro-scale spatial parameters is
         *       dim + dim * dim.
         *
         * :param const uIntType &dim: The dimension of the problem.
         * :param const uIntType &nMicroNodes: The number of micro-nodes in the N matrix
         * :param const uIntType &nMacroNodes: The number of macro-nodes in the N matrix
         * :param const uIntVector &domainMicroNodeIndices: The global micro-node indices in the given domain
         * :param const uIntVector &domainMacroDOFIndices: The indices of the macro-scale nodes present in the domain
         *     these are ( in a macro-scale FEA context ) the nodes of the micromorphic finite element.
         * :param const std::unordered_map< uIntType, floatType > &domainReferenceXis: The vectors from the local center of
         *     mass to the micro-nodes in the domain.
         * :param const floatVector &domainMacroInterpolationFunctionValues: The values of the interpolation functions
         *     at the local center of mass.
         * :param const std::unordered_map< uIntType, floatType >: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param SparseMatrix &domainN: The macro domain's interpolation function.
         * :param std::unordered_map< uIntType, uIntType > *microNodeToLocalIndex: A map from the micro node index 
         *     to the indices to be used in the output vector. The micro nodes which are either free or ghost may not be all
         *     of the micro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         *
         *     Note: If a node is not located in microNodeToLocalIndex ( and microNodeToLocalIndex is not NULL ) this node
         *           will be skipped.
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         */

        //Error handling
        if ( dim != 3 ){
            return new errorNode( "formMacroDomainToMicroInterpolationMatrix",
                                  "Only 3D domains are currently supported" );
        }

        if ( !microNodeToLocalIndex ){
            if ( nMicroNodes != microWeights.size() ){
                return new errorNode( "formMacroDomainToMicroInterpolationMatrix",
                                      "The number of micro nodes is not equal to the number of weights" );
            }
        }

        if ( domainMacroNodeIndices.size() != domainMacroInterpolationFunctionValues.size() ){
            return new errorNode( "formMacroDomainToMicroInterpolationMatrix",
                                  "The number of macro indices is not equal to the number of macro interpolation function values" );
        }

        //Set the number of spatial degrees of freedom for the two scales
        const uIntType nMicroDOF = dim;
        const uIntType nMacroDOF = dim + dim * dim;

        //Set up the vector of terms for the sparse matrix
        std::vector< T > coefficients;
        coefficients.reserve( nMicroDOF * domainMicroNodeIndices.size() * nMacroDOF * domainMacroNodeIndices.size() );

        //Initialize the row and column indices for the sparse matrix
        uIntType row0, col0;

        //Initialize the Xi vector
        floatVector Xi;

        //Initialize the value of the weight and shapefunction value
        floatType w, sf;

        //Initialize the global micro node index
        uIntType m;

        //Initialize the global macro node index
        uIntType n;

        //Initialize the output index for the micro-nodes
        uIntType o;

        //Initialize the output index for the macro-nodes
        uIntType p;

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){

            //Set the global micro node index
            m = domainMicroNodeIndices[ i ];

            auto microWeight = microWeights.find( m );

            if ( microWeight == microWeights.end( ) ){
                return new errorNode( "formMacroDomainToMicroInterpolationMatrix",
                                      "The micro node " + std::to_string( m ) + " was not found in the micro weight map" );
            }

            auto referenceXi = domainReferenceXis.find( m );

            if ( referenceXi == domainReferenceXis.end( ) ){
                return new errorNode( "formMacroDomainToMicroInterpolationMatrix",
                                      "The micro node " + std::to_string( m ) + " was not found in the reference Xi vector map" );
            }

            //Set the row index
            if ( microNodeToLocalIndex ){

                auto indx = microNodeToLocalIndex->find( m );

                if ( indx == microNodeToLocalIndex->end( ) ){

                    continue;

                }

                o = indx->second;
            }
            else{
                o = domainMicroNodeIndices[ i ];
            }

            row0 = nMicroDOF * o;

            //Set the micro-position vector
            Xi = referenceXi->second;

            //Set the value of the weight
            w = microWeight->second;

            for ( uIntType j = 0; j < domainMacroNodeIndices.size(); j++ ){

                //Set the global macro node index
                n = domainMacroNodeIndices[ j ];

                //Set the column index
                if ( macroNodeToLocalIndex ){
    
                    auto indx = macroNodeToLocalIndex->find( n );
    
                    if ( indx == macroNodeToLocalIndex->end( ) ){
    
                        return new errorNode( "formMacroDomaintoMicroInterpolationMatrix",
                                              "The macro node " + std::to_string( n ) + " is not found in the macro node to local index map" );
    
                    }
    
                    p = indx->second;
    
                }
                else{
                    p = domainMacroNodeIndices[ j ];
                }

                col0 = nMacroDOF * p;

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
                                                    const floatVector &microWeights,
                                                    floatVector &projectedMicroMasses,
                                                    const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex ){
        /*!
         * Add the contribution of the micro-nodes' mass to the macro nodes.
         *
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes in the domain
         * :param const floatVector &microMasses: The masses of the micro-nodes
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &projectedMicroMasses: The projected micro-masses on all of the macro-scale nodes.
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         */

        //Error handling
        if ( domainMacroNodeIndices.size() * domainMicroNodeIndices.size() != domainMicroShapeFunctions.size() ){
            return new errorNode( "addDomainMicroContributionToMacroMass",
                                  "The size of the domain node indices vectors are not consistent with the number of shape functions" );
        }

        if ( microWeights.size() != microMasses.size() ){
            return new errorNode( "addDomainMicroContributionToMacroMass", "The size of the domain's micro weights vector is not equal to the number of micro masses" );
        }

        for ( uIntType i = 0; i < domainMacroNodeIndices.size(); i++ ){
            if ( !macroNodeToLocalIndex ){
                if ( domainMacroNodeIndices[ i ] >= projectedMicroMasses.size() ){
                    return new errorNode( "addDomainMicroContributionToMacroMass",
                                          "The size of the projected micro mass vector is smaller than a macro node requires" );
                }
            }
            else{
                if ( macroNodeToLocalIndex->find( domainMacroNodeIndices[ i ] ) == macroNodeToLocalIndex->end( ) ){
                    return new errorNode( "addDomainMicroContributionToMacroMass",
                                          std::to_string( domainMacroNodeIndices[ i ] ) + " is not found in the gobal to local macro id map" );
                }
            }
        }

        //Set the number of macro-nodes in this domain
        uIntType nMacroNodes = domainMacroNodeIndices.size();

        //Initialize the micro node mass
        floatType mass;

        //Initialize the micro node weights
        floatType weight;

        //Initialize the global micro id
        uIntType m;

        //Initialize the global macro id
        uIntType n;

        //Initialize the local macro id
        uIntType p;

        //Iterate over the micro-node indices
        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){

            //Set the global micro id
            m = domainMicroNodeIndices[ i ];

            //Error handling
            if ( m >= microWeights.size() ){

                return new errorNode( "addDomainMicroContributionToMacroMass",
                                      "The size of the micro weights vector is smaller than a micro node requires" );

            }
            if ( m >= microMasses.size() ){

                return new errorNode( "addDomainMicroContributionToMacroMass",
                                      "The micro-node index is too large for the provided mass and shape function vectors" );

            }

            //Get the micro-node's mass
            mass = microMasses[ m ];

            //Get the micro-node's weight
            weight = microWeights[ m ];

            //Iterate over the macro-node indices
            for ( uIntType j = 0; j < domainMacroNodeIndices.size(); j++ ){

                //Set the global macro id
                n = domainMacroNodeIndices[ j ];

                //Find the local macro node index
                if ( macroNodeToLocalIndex ){

                    auto indx = macroNodeToLocalIndex->find( n );
                    p = indx->second;

                }
                else{

                    p = n;

                }

                //Add the micro-node contribution to the macro node
                projectedMicroMasses[ p ] += mass * domainMicroShapeFunctions[ i * nMacroNodes + j ] * weight;

            }

        }

        return NULL;
    }

    errorOut addDomainMicroContributionToMacroMass( const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                    const floatVector &microVolumes, const floatVector &microDensities,
                                                    const floatVector &domainMicroShapeFunctions,
                                                    const floatVector &microWeights,
                                                    floatVector &projectedMicroMasses,
                                                    const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex ){
        /*!
         * Add the contribution of the micro-nodes' mass to the macro nodes.
         *
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes in the domain
         * :param const floatVector &microVolumes: The volumes of the micro-nodes
         * :param const floatVector &microDensities: The densities of the micro-nodes
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &projectedMicroMasses: The projected micro-masses on all of the macro-scale nodes.
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         */

        //Error handling
        if ( domainMacroNodeIndices.size() * domainMicroNodeIndices.size() != domainMicroShapeFunctions.size() ){
            return new errorNode( "addDomainMicroContributionToMacroMass",
                                  "The size of the domain node indices vectors are not consistent with the number of shape functions" );
        }

        if ( microWeights.size() != microVolumes.size() ){
            return new errorNode( "addDomainMicroContributionToMacroMass", "The size of the domain's micro weights vector is not equal to the number of micro volumes" );
        }

        if ( microWeights.size() != microDensities.size() ){
            return new errorNode( "addDomainMicroContributionToMacroMass", "The size of the domain's micro weights vector is not equal to the number of micro densities" );
        }

        for ( uIntType i = 0; i < domainMacroNodeIndices.size(); i++ ){
            if ( !macroNodeToLocalIndex ){
                if ( domainMacroNodeIndices[ i ] >= projectedMicroMasses.size() ){
                    return new errorNode( "addDomainMicroContributionToMacroMass",
                                          "The size of the projected micro mass vector is smaller than a macro node requires" );
                }
            }
            else{
                if ( macroNodeToLocalIndex->find( domainMacroNodeIndices[ i ] ) == macroNodeToLocalIndex->end( ) ){
                    return new errorNode( "addDomainMicroContributionToMacroMass",
                                          std::to_string( domainMacroNodeIndices[ i ] ) + " is not found in the gobal to local macro id map" );
                }
            }
        }

        //Set the number of macro-nodes in this domain
        uIntType nMacroNodes = domainMacroNodeIndices.size();

        //Initialize the micro node mass
        floatType mass;

        //Initialize the micro node weights
        floatType weight;

        //Initialize the global micro id
        uIntType m;

        //Initialize the global macro id
        uIntType n;

        //Initialize the local macro id
        uIntType p;

        //Iterate over the micro-node indices
        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){

            //Set the global micro id
            m = domainMicroNodeIndices[ i ];

            //Error handling
            if ( m >= microWeights.size() ){

                    return new errorNode( "addDomainMicroContributionToMacroMass",
                                      "The size of the micro weights vector is smaller than a micro node requires" );

            }

            if ( m >= microVolumes.size() ){

                return new errorNode( "addDomainMicroContributionToMacroMass",
                                      "The micro-node index is too large for the provided volume vector" );

            }

            if ( m >= microDensities.size() ){

                return new errorNode( "addDomainMicroContributionToMacroMass",
                                      "The micro-node index is too large for the provided density vector" );

            }

            //Get the micro-node's mass
            mass = microVolumes[ m ] * microDensities[ m ];

            //Get the micro-node's weight
            weight = microWeights[ m ];

            //Iterate over the macro-node indices
            for ( uIntType j = 0; j < domainMacroNodeIndices.size(); j++ ){

                //Set the global macro id
                n = domainMacroNodeIndices[ j ];

                //Find the local macro node index
                if ( macroNodeToLocalIndex ){

                    auto indx = macroNodeToLocalIndex->find( n );
                    p = indx->second;

                }
                else{

                    p = n;

                }

                //Add the micro-node contribution to the macro node
                projectedMicroMasses[ p ] += mass * domainMicroShapeFunctions[ i * nMacroNodes + j ] * weight;

            }

        }

        return NULL;
    }
 
    errorOut addDomainMicroContributionToMacroMicroMassMomentOfInertia( const uIntType &dim,
                                                                        const uIntVector &domainMicroNodeIndices,
                                                                        const uIntVector &domainMacroNodeIndices,
                                                                        const floatVector &domainReferenceXis,
                                                                        const floatVector &microMasses,
                                                                        const floatVector &domainMicroShapeFunctions,
                                                                        const floatVector &microWeights,
                                                                        floatVector &projectedMassMicroMomentOfInertia,
                                                                        const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex
                                                                      ){
        /*!
         * Add the contribution of the micro-nodes in the domain to the macro moment of inertia.
         *
         * :param const uIntType &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &domainReferenceXis: The micro-position vectors in the domain
         * :param const floatVector &microMasses: The masses of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &projectedMassMicroMomentOfInertia: The moments of inertia at the macro-nodes of the domain
         *     as projected from the micro-nodes weighted by the mass.
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         */

        floatVector projectedMassConstant;
        floatVector projectedMassDisplacement;
        floatVector projectedMassDisplacementPosition;
        floatVector microDisplacements;
        return addDomainMicroToMacroProjectionTerms( dim, domainMicroNodeIndices, domainMacroNodeIndices,
                                                     domainReferenceXis, microMasses, domainMicroShapeFunctions,
                                                     microWeights, microDisplacements,
                                                     projectedMassMicroMomentOfInertia,
                                                     projectedMassConstant,
                                                     projectedMassDisplacement,
                                                     projectedMassDisplacementPosition,
                                                     true, false, false, false,
                                                     macroNodeToLocalIndex
                                                    );

    }

    errorOut addDomainMassConstant( const uIntType &dim,
                                    const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                    const floatVector &domainReferenceXis, const floatVector &microMasses,
                                    const floatVector &domainMicroShapeFunctions, const floatVector &microWeights,
                                    floatVector &projectedMassConstant,
                                    const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex
                                  ){
        /*!
         * Add the contributions of the domain to the mass constant
         *
         * :param const uIntType &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &domainReferenceXis: The micro-position vectors in the domain
         * :param const floatVector &microMasses: The masses of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &projectedMassConstant: The projected mass constant at the macro-scale node.
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         */

        floatVector projectedMassMicroMomentOfInertia;
        floatVector projectedMassDisplacement;
        floatVector projectedMassDisplacementPosition;
        floatVector microDisplacements;
        return addDomainMicroToMacroProjectionTerms( dim, domainMicroNodeIndices, domainMacroNodeIndices,
                                                     domainReferenceXis, microMasses, domainMicroShapeFunctions,
                                                     microWeights, microDisplacements,
                                                     projectedMassMicroMomentOfInertia,
                                                     projectedMassConstant,
                                                     projectedMassDisplacement,
                                                     projectedMassDisplacementPosition,
                                                     false, true, false, false,
                                                     macroNodeToLocalIndex
                                                   );

    }

    errorOut addDomainMassDisplacement( const uIntType &dim,
                                        const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                        const floatVector &microMasses, const floatVector &domainMicroShapeFunctions,
                                        const floatVector &microWeights, const floatVector &microDisplacements,
                                        floatVector &projectedMassDisplacement,
                                        const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex
                                      ){
        /*!
         * Add the contributions of the domain to the mass constant
         *
         * :param const uIntType &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &microMasses: The masses of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &microDisplacements: The displacements of the micro-degrees of freedom.
         * :param floatVector &projectedMassDisplacement: The projected mass weighted displacement at the macro-scale node.
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         */

        floatVector projectedMassMicroMomentOfInertia;
        floatVector projectedMassConstant;
        floatVector projectedMassDisplacementPosition;
        floatVector domainReferenceXis;
        return addDomainMicroToMacroProjectionTerms( dim, domainMicroNodeIndices, domainMacroNodeIndices,
                                                     domainReferenceXis, microMasses, domainMicroShapeFunctions,
                                                     microWeights, microDisplacements,
                                                     projectedMassMicroMomentOfInertia,
                                                     projectedMassConstant,
                                                     projectedMassDisplacement,
                                                     projectedMassDisplacementPosition,
                                                     false, false, true, false,
                                                     macroNodeToLocalIndex
                                                   );

    }

    errorOut addDomainMassMicroDisplacementPosition( const uIntType &dim,
                                                     const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                     const floatVector &domainReferenceXis, const floatVector &microMasses,
                                                     const floatVector &domainMicroShapeFunctions, const floatVector &microWeights,
                                                     const floatVector &microDisplacements,
                                                     floatVector &projectedMassDisplacementPosition,
                                                     const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex
                                                   ){
        /*!
         * Add the contributions of the domain to the mass constant
         *
         * :param const uIntType &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &domainReferenceXis: The micro-position vectors in the domain
         * :param const floatVector &microMasses: The masses of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &microDisplacements: The displacements of the micro-degrees of freedom.
         * :param floatVector &projectedMassDisplacementPosition: The projected mass weighted dyadic product of displacement
         *     and the micro-position at the macro-scale node.
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         */

        floatVector projectedMassMicroMomentOfInertia;
        floatVector projectedMassConstant;
        floatVector projectedMassDisplacement;
        return addDomainMicroToMacroProjectionTerms( dim, domainMicroNodeIndices, domainMacroNodeIndices,
                                                     domainReferenceXis, microMasses, domainMicroShapeFunctions,
                                                     microWeights, microDisplacements,
                                                     projectedMassMicroMomentOfInertia,
                                                     projectedMassConstant,
                                                     projectedMassDisplacement,
                                                     projectedMassDisplacementPosition,
                                                     false, false, false, true,
                                                     macroNodeToLocalIndex
                                                   );

    }


    errorOut addDomainMicroToMacroProjectionTerms( const uIntType &dim,
                                                   const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                   const floatVector &domainReferenceXis, const floatVector &microMasses,
                                                   const floatVector &domainMicroShapeFunctions, const floatVector &microWeights,
                                                   const floatVector &microDisplacements,
                                                   floatVector &projectedMassMicroMomentOfInertia,
                                                   floatVector &projectedMassConstant,
                                                   floatVector &projectedMassDisplacement,
                                                   floatVector &projectedMassDisplacementPosition,
                                                   const bool computeMassMomentOfInertia,
                                                   const bool computeMassConstant,
                                                   const bool computeMassMicroDisplacement,
                                                   const bool computeMassDisplacementPosition,
                                                   const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex
                                                 ){
        /*!
         * Solve for the terms required to project from the micro-scale to the macro-scale.
         *
         * :param const uIntType &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &domainReferenceXis: The micro-position vectors in the domain
         * :param const floatVector &microMasses: The masses of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
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
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         */

        //Error handling
        if ( ( computeMassMomentOfInertia ) || ( computeMassConstant ) || ( computeMassDisplacementPosition ) ){
            if ( dim * domainMicroNodeIndices.size() != domainReferenceXis.size() ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The number of micro node indices and the micro position vectors do not have consistent sizes" );
            }
        }

        if ( microWeights.size( ) != microMasses.size( ) ){
            return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                  "The micro weight and micro mass vectors are not consistent in size" );
        }

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){
            if ( domainMicroNodeIndices[ i ] >= microWeights.size() ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The number of micro node weights is smaller than the micro indices requires" );
            }
        }

        if ( domainMicroNodeIndices.size() * domainMacroNodeIndices.size() != domainMicroShapeFunctions.size() ){
            return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                  "The number of micro and micro node indices are not consistent with the number of shape functions" );
        }

        for ( uIntType i = 0; i < domainMacroNodeIndices.size(); i++ ){
            uIntType n = domainMacroNodeIndices[ i ];
            uIntType p = n;

            if ( macroNodeToLocalIndex ){

                auto indx = macroNodeToLocalIndex->find( n );

                if ( indx == macroNodeToLocalIndex->end( ) ){
                    return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                          "Macro node " + std::to_string( n ) + " was not found in macroNodeToLocalIndex" );
                }
                else{
                    p = indx->second;
                }

            }

            if ( ( projectedMassMicroMomentOfInertia.size() < dim * dim * ( p + 1 ) ) &&
                 ( computeMassMomentOfInertia ) ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The size of the projected micro moment of inertia weighted by the mass is smaller than required for the provided nodes" );
            }

            if ( ( projectedMassConstant.size() < dim * ( p + 1 ) ) &&
                 ( computeMassConstant ) ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The size of the projected mass constant is smaller than required for the provided nodes" );
            }

            if ( ( projectedMassDisplacement.size() < dim * ( p + 1 ) ) &&
                 ( computeMassMicroDisplacement ) ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The size of the projected mass-weighted micro displacement is smaller than required for the provided nodes" );
            }

            if ( ( projectedMassDisplacementPosition.size() < dim * dim * ( p + 1 ) ) &&
                 ( computeMassDisplacementPosition ) ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The size of the projected mass-weighted dyadic product of the micro displacement and the micro position is smaller than required for the provided nodes" );
            }

        }

        if ( ( computeMassMicroDisplacement ) || ( computeMassDisplacementPosition ) ){
            for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){
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

        //Initialize the micro-node global id
        uIntType m;

        //Initialize the macro-node global id
        uIntType n;

        //Initialize the macro-node local id
        uIntType p;

        //Loop through the micro nodes
        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){

            //Set the micro-node global id
            m = domainMicroNodeIndices[ i ];

            if ( m >= microMasses.size( ) ){

                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The micro node index " + std::to_string( m ) + " is too large for the micro-mass vector" );

            }

            //Extract the nodal mass
            mass = microMasses[ m ];

            //Extract the nodal weight
            weight = microWeights[ m ];

            //Extract the current micro-position vectors
            if ( ( computeMassConstant ) || ( computeMassDisplacementPosition ) ){

                Xi = floatVector( domainReferenceXis.begin() + dim * i, domainReferenceXis.begin() + dim * ( i + 1 ) );

            }

            //Compute the dyadic product of Xi with itself
            if ( computeMassMomentOfInertia ){

                for ( uIntType j = 0; j < dim; j++ ){
    
                    for ( uIntType k = 0; k < dim; k++ ){
    
                        XiXi[ dim * j + k ] = domainReferenceXis[ dim * i + j ] * domainReferenceXis[ dim * i + k ];
    
                    }
    
                }
            }

            if ( ( computeMassMicroDisplacement ) || ( computeMassDisplacementPosition ) ){

                q = floatVector( microDisplacements.begin() + dim * m,
                                 microDisplacements.begin() + dim * m + dim );

            }

            //Loop through the macro nodes
            for ( uIntType j = 0; j < domainMacroNodeIndices.size(); j++ ){

                //Set the macro-node global id
                n = domainMacroNodeIndices[ j ];
                p = n;

                if ( macroNodeToLocalIndex ){
                    auto indx = macroNodeToLocalIndex->find( n );
                    p = indx->second;
                }


                //Set the shape-function value
                sf = domainMicroShapeFunctions[ domainMacroNodeIndices.size() * i + j ];
               
                if ( computeMassMomentOfInertia ){ 

                    for ( uIntType k = 0; k < dim * dim; k++ ){
    
                        //Add the contribution to the micro-moment of inertia
                        projectedMassMicroMomentOfInertia[ dim * dim * p + k ]
                            += weight * mass * sf * XiXi[ k ];
    
                    }

                }

                if ( computeMassConstant ){ 

                    for ( uIntType k = 0; k < dim; k++ ){
    
                        //Add the contribution to the micro-moment of inertia
                        projectedMassConstant[ dim * p + k ] += weight * mass * sf * Xi[ k ];
    
                    }

                }

                if ( computeMassMicroDisplacement ) {

                    for ( uIntType k = 0; k < dim; k++ ){

                        //Add the contribution to the mass weighted micro-displacement
                        projectedMassDisplacement[ dim * p + k ] += weight * mass * sf * q[ k ];

                    }

                }

                if ( computeMassDisplacementPosition ){

                    for ( uIntType k = 0; k < dim; k++ ){

                        for ( uIntType l = 0; l < dim; l++ ){
                            projectedMassDisplacementPosition[ dim * dim * p + dim * k + l ]
                                += weight * mass * sf * q[ k ] * Xi[ l ];
                        }

                    }

                }

            }

        }

        return NULL;

    }

    errorOut addDomainMicroContributionToMacroMicroMassMomentOfInertia( const uIntType &dim,
                                                                        const uIntVector &domainMicroNodeIndices,
                                                                        const uIntVector &domainMacroNodeIndices,
                                                                        const floatVector &domainReferenceXis,
                                                                        const floatVector &microVolumes,
                                                                        const floatVector &microDensities,
                                                                        const floatVector &domainMicroShapeFunctions,
                                                                        const floatVector &microWeights,
                                                                        floatVector &projectedMassMicroMomentOfInertia,
                                                                        const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex
                                                                      ){
        /*!
         * Add the contribution of the micro-nodes in the domain to the macro moment of inertia.
         *
         * :param const uIntType &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &domainReferenceXis: The micro-position vectors in the domain
         * :param const floatVector &microVolumes: The volumes of the micro-nodes.
         * :param const floatVector &microDensities: The densities of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &projectedMassMicroMomentOfInertia: The moments of inertia at the macro-nodes of the domain
         *     as projected from the micro-nodes weighted by the mass.
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         */

        floatVector projectedMassConstant;
        floatVector projectedMassDisplacement;
        floatVector projectedMassDisplacementPosition;
        floatVector microDisplacements;
        return addDomainMicroToMacroProjectionTerms( dim, domainMicroNodeIndices, domainMacroNodeIndices,
                                                     domainReferenceXis, microVolumes, microDensities,
                                                     domainMicroShapeFunctions,
                                                     microWeights, microDisplacements,
                                                     projectedMassMicroMomentOfInertia,
                                                     projectedMassConstant,
                                                     projectedMassDisplacement,
                                                     projectedMassDisplacementPosition,
                                                     true, false, false, false,
                                                     macroNodeToLocalIndex );

    }

    errorOut addDomainMassConstant( const uIntType &dim,
                                    const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                    const floatVector &domainReferenceXis, const floatVector &microVolumes,
                                    const floatVector &microDensities, const floatVector &domainMicroShapeFunctions,
                                    const floatVector &microWeights,
                                    floatVector &projectedMassConstant,
                                    const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex ){
        /*!
         * Add the contributions of the domain to the mass constant
         *
         * :param const uIntType &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &domainReferenceXis: The micro-position vectors in the domain
         * :param const floatVector &microVolumes: The volumes of the micro-nodes.
         * :param const floatVector &microDensities: The densities of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &projectedMassConstant: The projected mass constant at the macro-scale node.
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         */

        floatVector projectedMassMicroMomentOfInertia;
        floatVector projectedMassDisplacement;
        floatVector projectedMassDisplacementPosition;
        floatVector microDisplacements;
        return addDomainMicroToMacroProjectionTerms( dim, domainMicroNodeIndices, domainMacroNodeIndices,
                                                     domainReferenceXis, microVolumes, microDensities,
                                                     domainMicroShapeFunctions,
                                                     microWeights, microDisplacements,
                                                     projectedMassMicroMomentOfInertia,
                                                     projectedMassConstant,
                                                     projectedMassDisplacement,
                                                     projectedMassDisplacementPosition,
                                                     false, true, false, false,
                                                     macroNodeToLocalIndex );

    }

    errorOut addDomainMassDisplacement( const uIntType &dim,
                                        const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                        const floatVector &microVolumes, const floatVector &microDensities,
                                        const floatVector &domainMicroShapeFunctions,
                                        const floatVector &microWeights, const floatVector &microDisplacements,
                                        floatVector &projectedMassDisplacement,
                                        const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex ){
        /*!
         * Add the contributions of the domain to the mass constant
         *
         * :param const uIntType &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &microVolumes: The volumes of the micro-nodes.
         * :param const floatVector &microDensities: The densities of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &microDisplacements: The displacements of the micro-degrees of freedom.
         * :param floatVector &projectedMassDisplacement: The projected mass weighted displacement at the macro-scale node.
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         */

        floatVector projectedMassMicroMomentOfInertia;
        floatVector projectedMassConstant;
        floatVector projectedMassDisplacementPosition;
        floatVector domainReferenceXis;
        return addDomainMicroToMacroProjectionTerms( dim, domainMicroNodeIndices, domainMacroNodeIndices,
                                                     domainReferenceXis, microVolumes, microDensities,
                                                     domainMicroShapeFunctions,
                                                     microWeights, microDisplacements,
                                                     projectedMassMicroMomentOfInertia,
                                                     projectedMassConstant,
                                                     projectedMassDisplacement,
                                                     projectedMassDisplacementPosition,
                                                     false, false, true, false,
                                                     macroNodeToLocalIndex );

    }

    errorOut addDomainMassMicroDisplacementPosition( const uIntType &dim,
                                                     const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                     const floatVector &domainReferenceXis, const floatVector &microVolumes,
                                                     const floatVector &microDensities,
                                                     const floatVector &domainMicroShapeFunctions, const floatVector &microWeights,
                                                     const floatVector &microDisplacements,
                                                     floatVector &projectedMassDisplacementPosition,
                                                     const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex ){
        /*!
         * Add the contributions of the domain to the mass constant
         *
         * :param const uIntType &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &domainReferenceXis: The micro-position vectors in the domain
         * :param const floatVector &microVolumes: The volumes of the micro-nodes.
         * :param const floatVector &microDensities: The densities of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &microDisplacements: The displacements of the micro-degrees of freedom.
         * :param floatVector &projectedMassDisplacementPosition: The projected mass weighted dyadic product of displacement
         *     and the micro-position at the macro-scale node.
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         */

        floatVector projectedMassMicroMomentOfInertia;
        floatVector projectedMassConstant;
        floatVector projectedMassDisplacement;
        return addDomainMicroToMacroProjectionTerms( dim, domainMicroNodeIndices, domainMacroNodeIndices,
                                                     domainReferenceXis, microVolumes, microDensities, domainMicroShapeFunctions,
                                                     microWeights, microDisplacements,
                                                     projectedMassMicroMomentOfInertia,
                                                     projectedMassConstant,
                                                     projectedMassDisplacement,
                                                     projectedMassDisplacementPosition,
                                                     false, false, false, true,
                                                     macroNodeToLocalIndex );

    }

    errorOut addDomainMicroToMacroProjectionTerms( const uIntType &dim,
                                                   const uIntVector &domainMicroNodeIndices, const uIntVector &domainMacroNodeIndices,
                                                   const floatVector &domainReferenceXis, const floatVector &microVolumes,
                                                   const floatVector &microDensities, const floatVector &domainMicroShapeFunctions,
                                                   const floatVector &microWeights,
                                                   const floatVector &microDisplacements,
                                                   floatVector &projectedMassMicroMomentOfInertia,
                                                   floatVector &projectedMassConstant,
                                                   floatVector &projectedMassDisplacement,
                                                   floatVector &projectedMassDisplacementPosition,
                                                   const bool computeMassMomentOfInertia,
                                                   const bool computeMassConstant,
                                                   const bool computeMassMicroDisplacement,
                                                   const bool computeMassDisplacementPosition,
                                                   const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex
                                                 ){
        /*!
         * Solve for the terms required to project from the micro-scale to the macro-scale.
         *
         * :param const uIntType &dim: The dimension of the problem.
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const uIntVector &domainMacroNodeIndices: The indices of the macro-nodes associated with the domain
         * :param const floatVector &domainReferenceXis: The micro-position vectors in the domain
         * :param const floatVector &microVolumes: The volumes of the micro-nodes.
         * :param const floatVector &microDensities: The densities of the micro-nodes.
         * :param const floatVector &domainMicroShapeFunctions: The shape functions of the macro interpolation functions
         *     at the micro nodes. Organized as [ N_11, N_12, N_13, ... N_21, N_22, ... ] where the first index is the 
         *     micro-node number and the second is the macro-node number.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
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
         * :param std::unordered_map< uIntType, uIntType > *macroNodeToLocalIndex: A map from the macro node index 
         *     to the indices to be used in the output vector. The macro nodes which are either free or ghost may not be all
         *     of the macro-scale nodes so the projection matrices would include large zero regions and be ordered in a less
         *     than optimal way. We can use this to define the mapping better. This defaults to NULL so the global ID index values
         *     will be used.
         */

        //Error handling
        if ( ( computeMassMomentOfInertia ) || ( computeMassConstant ) || ( computeMassDisplacementPosition ) ){
            if ( dim * domainMicroNodeIndices.size() != domainReferenceXis.size() ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The number of micro node indices and the micro position vectors do not have consistent sizes" );
            }
        }

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){
            if ( domainMicroNodeIndices[ i ] >= microWeights.size() ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The number of micro node weights is smaller than the micro indices requires" );
            }
        }

        if ( microWeights.size( ) != microVolumes.size( ) ){
            return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                  "The micro weight and micro volume vectors are not consistent in size" );
        }

        if ( microVolumes.size( ) != microDensities.size( ) ){
            return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                  "The micro density and micro volume vectors are not consistent in size" );
        }

        if ( domainMicroNodeIndices.size() * domainMacroNodeIndices.size() != domainMicroShapeFunctions.size() ){
            return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                  "The number of micro and micro node indices are not consistent with the number of shape functions" );
        }

        for ( uIntType i = 0; i < domainMacroNodeIndices.size(); i++ ){
            uIntType n = domainMacroNodeIndices[ i ];
            uIntType p = n;

            if ( macroNodeToLocalIndex ){

                auto indx = macroNodeToLocalIndex->find( n );

                if ( indx == macroNodeToLocalIndex->end( ) ){
                    return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                          "Macro node " + std::to_string( n ) + " was not found in macroNodeToLocalIndex" );
                }
                else{
                    p = indx->second;
                }

            }

            if ( ( projectedMassMicroMomentOfInertia.size() < dim * dim * ( p + 1 ) ) &&
                 ( computeMassMomentOfInertia ) ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The size of the projected micro moment of inertia weighted by the mass is smaller than required for the provided nodes" );
            }

            if ( ( projectedMassConstant.size() < dim * ( p + 1 ) ) &&
                 ( computeMassConstant ) ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The size of the projected mass constant is smaller than required for the provided nodes" );
            }

            if ( ( projectedMassDisplacement.size() < dim * ( p + 1 ) ) &&
                 ( computeMassMicroDisplacement ) ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The size of the projected mass-weighted micro displacement is smaller than required for the provided nodes" );
            }

            if ( ( projectedMassDisplacementPosition.size() < dim * dim * ( p + 1 ) ) &&
                 ( computeMassDisplacementPosition ) ){
                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The size of the projected mass-weighted dyadic product of the micro displacement and the micro position is smaller than required for the provided nodes" );
            }

        }

        if ( ( computeMassMicroDisplacement ) || ( computeMassDisplacementPosition ) ){
            for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){
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

        //Initialize the micro-node global id
        uIntType m;

        //Initialize the macro-node global id
        uIntType n;

        //Initialize the macro-node local id
        uIntType p;

        //Loop through the micro nodes
        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){

            //Set the micro-node global id
            m = domainMicroNodeIndices[ i ];

            if ( m >= microVolumes.size( ) ){

                return new errorNode( "addDomainMicroToMacroProjectionTerms",
                                      "The micro node index " + std::to_string( m )
                                      + " is too large for the micro-density and volume vectors" );

            }

            //Extract the nodal mass
            mass = microVolumes[ m ] * microDensities[ m ];

            //Extract the nodal weight
            weight = microWeights[ m ];

            //Extract the current micro-position vectors
            if ( ( computeMassConstant ) || ( computeMassDisplacementPosition ) ){

                Xi = floatVector( domainReferenceXis.begin() + dim * i, domainReferenceXis.begin() + dim * ( i + 1 ) );

            }

            //Compute the dyadic product of Xi with itself
            if ( computeMassMomentOfInertia ){

                for ( uIntType j = 0; j < dim; j++ ){
    
                    for ( uIntType k = 0; k < dim; k++ ){
    
                        XiXi[ dim * j + k ] = domainReferenceXis[ dim * i + j ] * domainReferenceXis[ dim * i + k ];
    
                    }
    
                }
            }

            if ( ( computeMassMicroDisplacement ) || ( computeMassDisplacementPosition ) ){

                q = floatVector( microDisplacements.begin() + dim * m,
                                 microDisplacements.begin() + dim * m + dim );

            }

            //Loop through the macro nodes
            for ( uIntType j = 0; j < domainMacroNodeIndices.size(); j++ ){

                //Set the macro-node global id
                n = domainMacroNodeIndices[ j ];
                p = n;

                if ( macroNodeToLocalIndex ){
                    auto indx = macroNodeToLocalIndex->find( n );
                    p = indx->second;
                }

                //Set the shape-function value
                sf = domainMicroShapeFunctions[ domainMacroNodeIndices.size() * i + j ];
               
                if ( computeMassMomentOfInertia ){ 

                    for ( uIntType k = 0; k < dim * dim; k++ ){
    
                        //Add the contribution to the micro-moment of inertia
                        projectedMassMicroMomentOfInertia[ dim * dim * p + k ]
                            += weight * mass * sf * XiXi[ k ];
    
                    }

                }

                if ( computeMassConstant ){ 

                    for ( uIntType k = 0; k < dim; k++ ){
    
                        //Add the contribution to the micro-moment of inertia
                        projectedMassConstant[ dim * p + k ] += weight * mass * sf * Xi[ k ];
    
                    }

                }

                if ( computeMassMicroDisplacement ) {

                    for ( uIntType k = 0; k < dim; k++ ){

                        //Add the contribution to the mass weighted micro-displacement
                        projectedMassDisplacement[ dim * p + k ] += weight * mass * sf * q[ k ];

                    }

                }

                if ( computeMassDisplacementPosition ){

                    for ( uIntType k = 0; k < dim; k++ ){

                        for ( uIntType l = 0; l < dim; l++ ){
                            projectedMassDisplacementPosition[ dim * dim * p + dim * k + l ]
                                += weight * mass * sf * q[ k ] * Xi[ l ];
                        }

                    }

                }

            }

        }

        return NULL;

    }

    errorOut computeDomainCenterOfMass( const uIntType &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microMasses,
                                        const floatVector &microPositions, const floatVector &microWeights,
                                        floatVector &domainCM ){
        /*!
         * Compute the center of mass of a micro domain from the masses of the micro-nodes contained
         * within the domain.
         *
         * :param const uIntType &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const floatVector &microMasses: The masses of the micro nodes.
         * :param const floatVector &microPositions: The positions of the micro-nodes.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &domainCG: The center of mass of the domain
         */

        floatType domainMass;
        return computeDomainCenterOfMass( dim, domainMicroNodeIndices, microMasses, microPositions, microWeights,
                                          domainMass, domainCM );
    }

    errorOut computeDomainCenterOfMass( const uIntType &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microVolumes,
                                        const floatVector &microDensities, const floatVector &microPositions,
                                        const floatVector &microWeights, floatVector &domainCM ){
        /*!
         * Compute the center of mass of a micro domain from the masses of the micro-nodes contained
         * within the domain.
         *
         * :param const uIntType &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const floatVector &microVolumes: The volumes of the micro nodes.
         * :param const floatVector &microDensities: The densities of the micro nodes.
         * :param const floatVector &microPositions: The positions of the micro-nodes.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &domainCG: The center of mass of the domain
         */

        floatType domainMass;
        return computeDomainCenterOfMass( dim, domainMicroNodeIndices, microVolumes, microDensities, microPositions,
                                          microWeights, domainMass, domainCM );
    }

    errorOut computeDomainCenterOfMass( const uIntType &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microVolumes,
                                        const floatVector &microDensities, const floatVector &microReferencePositions,
                                        const floatVector &microDisplacements, const floatVector &microWeights,
                                        floatVector &domainCM ){
        /*!
         * Compute the center of mass of a micro domain from the masses of the micro-nodes contained
         * within the domain.
         *
         * :param const uIntType &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const floatVector &microVolumes: The volumes of the micro nodes.
         * :param const floatVector &microDensities: The densities of the micro nodes.
         * :param const floatVector &microReferencePositions: The reference positions of the micro-nodes.
         * :param const floatVector &microDisplacements: The displacements of the micro-nodes.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &domainCG: The center of mass of the domain
         */

        floatType domainMass;
        return computeDomainCenterOfMass( dim, domainMicroNodeIndices, microVolumes, microDensities, microReferencePositions,
                                          microDisplacements, microWeights, domainMass, domainCM );
    }

    errorOut computeDomainCenterOfMass( const uIntType &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microMasses,
                                        const floatVector &microPositions, const floatVector &microWeights,
                                        floatType &domainMass, floatVector &domainCM ){
        /*!
         * Compute the center of mass of a micro domain from the masses of the micro-nodes contained
         * within the domain.
         *
         * :param const uIntType &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const floatVector &microMasses: The masses of the micro nodes.
         * :param const floatVector &microPositions: The positions of the micro-nodes.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &domainMass: The mass of the domain
         * :param floatVector &domainCM: The center of mass of the domain
         */

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){
            if ( microPositions.size() <= dim * domainMicroNodeIndices[ i ] + dim ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-positions vector is not consistent with the micro indices" );
            }

            if ( microMasses.size() <= domainMicroNodeIndices[ i ] ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-masses vector is not consistent with the micro indices" );
            }

            if ( microWeights.size( ) <= domainMicroNodeIndices[ i ] ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "the size of the micro-weights vector is not consistent with the micro indices" );
            }
        }

        //Initialize the domain mass
        domainMass = 0;

        //Initialize the center of mass vector
        domainCM = floatVector( dim, 0 );

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){


            //Add to the domain's mass
            domainMass += microMasses[ domainMicroNodeIndices[ i ] ] * microWeights[ domainMicroNodeIndices[ i ] ];

            //Add to the domain's mass weighted position
            domainCM += microMasses[ domainMicroNodeIndices[ i ] ] * microWeights[ domainMicroNodeIndices[ i ] ]
                      * floatVector( microPositions.begin() + dim * domainMicroNodeIndices[ i ],
                                     microPositions.begin() + dim * ( domainMicroNodeIndices[ i ] + 1 ) );

        }

        //Normalize the center of mass by the domain's mass
        domainCM /= domainMass;

        return NULL;

    }

    errorOut computeDomainCenterOfMass( const uIntType &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microVolumes,
                                        const floatVector &microDensities, const floatVector &microPositions,
                                        const floatVector &microWeights, floatType &domainMass, floatVector &domainCM ){
        /*!
         * Compute the center of mass of a micro domain from the masses of the micro-nodes contained
         * within the domain.
         *
         * :param const uIntType &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const floatVector &microVolumes: The volumes of the micro nodes.
         * :param const floatVector &microDensities: The densities of the micro nodes.
         * :param const floatVector &microPositions: The positions of the micro-nodes.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &domainMass: The mass of the domain
         * :param floatVector &domainCM: The center of mass of the domain
         */

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){
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

            if ( microWeights.size( ) <= domainMicroNodeIndices[ i ] ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "the size of the micro-weights vector is not consistent with the micro indices" );
            }
        }

        //Initialize the domain mass
        domainMass = 0;

        //Initialize the center of mass vector
        domainCM = floatVector( dim, 0 );

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){


            //Add to the domain's mass
            domainMass += microVolumes[ domainMicroNodeIndices[ i ] ] * microDensities[ domainMicroNodeIndices[ i ] ]
                        * microWeights[ domainMicroNodeIndices[ i ] ];

            //Add to the domain's mass weighted position
            domainCM += microVolumes[ domainMicroNodeIndices[ i ] ] * microDensities[ domainMicroNodeIndices[ i ] ]
                      * microWeights[ domainMicroNodeIndices[ i ] ]
                      * floatVector( microPositions.begin() + dim * domainMicroNodeIndices[ i ],
                                     microPositions.begin() + dim * ( domainMicroNodeIndices[ i ] + 1 ) );

        }

        //Normalize the center of mass by the domain's mass
        domainCM /= domainMass;

        return NULL;

    }

    errorOut computeDomainCenterOfMass( const uIntType &dim,
                                        const uIntVector &domainMicroNodeIndices, const floatVector &microVolumes,
                                        const floatVector &microDensities, const floatVector &microReferencePositions,
                                        const floatVector &microDisplacements, const floatVector &microWeights,
                                        floatType &domainMass, floatVector &domainCM ){
        /*!
         * Compute the center of mass of a micro domain from the masses of the micro-nodes contained
         * within the domain.
         *
         * :param const uIntType &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const floatVector &microVolumes: The volumes of the micro nodes.
         * :param const floatVector &microDensities: The densities of the micro nodes.
         * :param const floatVector &microReferencePositions: The reference positions of the micro-nodes.
         * :param const floatVector &microDisplacements: The displacements of the micro-nodes relative to their
         *     reference positions.
         * :param const floatVector &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &domainMass: The mass of the domain
         * :param floatVector &domainCM: The center of mass of the domain
         */

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){
            if ( microReferencePositions.size() < dim * domainMicroNodeIndices[ i ] + dim ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-reference positions vector is not consistent with the micro indices" );
            }

            if ( microDisplacements.size() < dim * domainMicroNodeIndices[ i ] + dim ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-displacements vector is not consistent with the micro indices" );
            }

            if ( microVolumes.size() < domainMicroNodeIndices[ i ] ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-volumes vector is not consistent with the micro indices" );
            }

            if ( microDensities.size() < domainMicroNodeIndices[ i ] ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-densities vector is not consistent with the micro indices" );
            }

            if ( microWeights.size( ) <= domainMicroNodeIndices[ i ] ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "the size of the micro-weights vector is not consistent with the micro indices" );
            }
        }

        //Initialize the domain mass
        domainMass = 0;

        //Initialize the center of mass vector
        domainCM = floatVector( dim, 0 );

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){


            //Add to the domain's mass
            domainMass += microVolumes[ domainMicroNodeIndices[ i ] ] * microDensities[ domainMicroNodeIndices[ i ] ]
                        * microWeights[ domainMicroNodeIndices[ i ] ];

            //Add to the domain's mass weighted position
            domainCM += microVolumes[ domainMicroNodeIndices[ i ] ] * microDensities[ domainMicroNodeIndices[ i ] ]
                      * microWeights[ domainMicroNodeIndices[ i ] ]
                      * ( floatVector( microReferencePositions.begin() + dim * domainMicroNodeIndices[ i ],
                                       microReferencePositions.begin() + dim * ( domainMicroNodeIndices[ i ] + 1 ) )
                        + floatVector( microDisplacements.begin() + dim * domainMicroNodeIndices[ i ],
                                       microDisplacements.begin() + dim * ( domainMicroNodeIndices[ i ] + 1 ) )
                        );


        }

        //Normalize the center of mass by the domain's mass
        domainCM /= domainMass;

        return NULL;

    }

    errorOut computeDomainCenterOfMass( const uIntType &dim,
                                        const uIntVector &domainMicroNodeIndices,
                                        const std::unordered_map< uIntType, floatType > &microVolumes,
                                        const std::unordered_map< uIntType, floatType > &microDensities,
                                        const std::unordered_map< uIntType, floatVector > &microReferencePositions,
                                        const std::unordered_map< uIntType, floatVector > &microDisplacements,
                                        const std::unordered_map< uIntType, floatType > &microWeights,
                                        floatType &domainMass, floatVector &domainCM ){
        /*!
         * Compute the center of mass of a micro domain from the masses of the micro-nodes contained
         * within the domain.
         *
         * :param const uIntType &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const std::unordered_map< uIntType, floatType > &microVolumes: The volumes of the micro nodes.
         * :param const std::unordered_map< uIntType, floatType > &microDensities: The densities of the micro nodes.
         * :param const std::unordered_map< uIntType, floatType > &microReferencePositions: The reference positions of the micro-nodes.
         * :param const std::unordered_map< uIntType, floatType > &microDisplacements: The displacements of the micro-nodes
         *     relative to their reference positions.
         * :param const std::unordered_map< uIntType, floatType > &microWeights: The weight associated with each micro-scale node.
         *     This is important for two cases:
         *     - Nodes which are shared between macro-scale domains. ( we don't want to double count )
         *     - Weighting the influence of nodes if nodes which have no mass are being used. This may be important
         *       if the minimum L2 norm projection is being used.
         * :param floatVector &domainMass: The mass of the domain
         * :param floatVector &domainCM: The center of mass of the domain
         */

        //Initialize the domain mass
        domainMass = 0;

        //Initialize the center of mass vector
        domainCM = floatVector( dim, 0 );

        for ( auto index = domainMicroNodeIndices.begin( ); index != domainMicroNodeIndices.end( ); index++ ){

            auto microVolume = microVolumes.find( *index );

            if ( microVolume == microVolumes.end( ) ){

                return new errorNode( "computeDomainCenterOfMass",
                                      "The micro index " + std::to_string( *index ) + " was not found in the micro volume map" );

            }

            auto microDensity = microDensities.find( *index );

            if ( microDensity == microDensities.end( ) ){

                return new errorNode( "computeDomainCenterOfMass",
                                      "The micro index " + std::to_string( *index ) + " was not found in the micro density map" );

            }

            auto microWeight = microWeights.find( *index );

            if ( microWeight == microWeights.end( ) ){

                return new errorNode( "computeDomainCenterOfMass",
                                      "The micro index " + std::to_string( *index ) + " was not found in the micro weight map" );

            }

            auto microReferencePosition = microReferencePositions.find( *index );

            if ( microReferencePosition == microReferencePositions.end( ) ){

                return new errorNode( "computeDomainCenterOfMass",
                                      "The micro index " + std::to_string( *index ) + " was not found in the micro reference position map" );

            }

            auto microDisplacement = microDisplacements.find( *index );

            if ( microDisplacement == microDisplacements.end( ) ){

                return new errorNode( "computeDomainCenterOfMass",
                                      "The micro index " + std::to_string( *index ) + " was not found in the micro displacement map" );

            }

            //Add to the domain's mass
            domainMass += microVolume->second * microDensity->second * microWeight->second;

            //Add to the domain's mass weighted position
            domainCM += microVolume->second * microDensity->second * microWeight->second
                      * ( microReferencePosition->second + microDisplacement->second );

        }

        //Normalize the center of mass by the domain's mass
        domainCM /= domainMass;

        return NULL;

    }

    errorOut computeDomainXis( const uIntType &dim,
                               const uIntVector &domainMicroNodeIndices, const floatVector &microPositions,
                               const floatVector &domainCM, floatVector &domainXis ){
        /*
         * Compute the relative position vector between the center of mass of a micro domain and the 
         * micro position.
         *
         * :param const uIntType &dim: The dimension of the problem
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

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){
            if ( microPositions.size() < dim * domainMicroNodeIndices[ i ] + dim ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-positions vector is not consistent with the micro indices" );
            }
        }

        //Resize the Xi vector
        domainXis.resize( dim * domainMicroNodeIndices.size() );

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){

            for ( uIntType j = 0; j < dim; j++ ){
           
                //Compute the relative position vector 
                domainXis[ dim * i + j ] = microPositions[ dim * domainMicroNodeIndices[ i ] + j ] - domainCM[ j ];

            }

        }

        return NULL;

    }

    errorOut computeDomainXis( const uIntType &dim,
                               const uIntVector &domainMicroNodeIndices, const floatVector &microReferencePositions,
                               const floatVector &microDisplacements, const floatVector &domainCM, floatVector &domainXis ){
        /*
         * Compute the relative position vector between the center of mass of a micro domain and the 
         * micro position.
         *
         * :param const uIntType &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const floatVector &microReferencePositions: The reference positions of the micro-nodes.
         * :param const floatVector &microDisplacements: The displacements of the micro-nodes.
         * :param floatVector &domainCM: The center of mass of the domain
         * :param floatVector &domainXis: The relative positions of the micro nodes.
         */

        //Error Handling
        if ( domainCM.size() != dim ){
            return new errorNode( "computeDomainXis",
                                  "The center of mass is not consistent with the dimension" );
        }

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){
            if ( microReferencePositions.size() < dim * domainMicroNodeIndices[ i ] + dim ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-reference positions vector is not consistent with the micro indices" );
            }
            if ( microDisplacements.size() < dim * domainMicroNodeIndices[ i ] + dim ){
                return new errorNode( "computeDomainCenterOfMass",
                                      "The size of the micro-displacements vector is not consistent with the micro indices" );
            }
        }

        //Resize the Xi vector
        domainXis.resize( dim * domainMicroNodeIndices.size() );

        for ( uIntType i = 0; i < domainMicroNodeIndices.size(); i++ ){

            for ( uIntType j = 0; j < dim; j++ ){

                //Compute the relative position vector 
                domainXis[ dim * i + j ] = (  microReferencePositions[ dim * domainMicroNodeIndices[ i ] + j ]
                                            + microDisplacements[ dim * domainMicroNodeIndices[ i ] + j ] ) - domainCM[ j ];

            }

        }

        return NULL;

    }

    errorOut computeDomainXis( const uIntType &dim,
                               const uIntVector &domainMicroNodeIndices,
                               const std::unordered_map< uIntType, floatVector > &microReferencePositions,
                               const std::unordered_map< uIntType, floatVector > &microDisplacements,
                               const floatVector &domainCM,
                               std::unordered_map< uIntType, floatVector > &domainXis ){
        /*
         * Compute the relative position vector between the center of mass of a micro domain and the 
         * micro position.
         *
         * :param const uIntType &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The indices of the micro-nodes in the domain.
         * :param const std::unordered_map< uIntType, floatVector > &microReferencePositions: The reference positions of the micro-nodes.
         * :param const std::unordered_map< uIntType, floatVector > &microDisplacements: The displacements of the micro-nodes.
         * :param floatVector &domainCM: The center of mass of the domain
         * :param std::unordered_map< uIntType, floatVector > &domainXis: The relative positions of the micro nodes.
         */

        //Error Handling
        if ( domainCM.size() != dim ){
            return new errorNode( "computeDomainXis",
                                  "The center of mass is not consistent with the dimension" );
        }

        //Resize the Xi map
        domainXis.reserve( domainMicroNodeIndices.size() );

        for ( auto index = domainMicroNodeIndices.begin( ); index != domainMicroNodeIndices.end( ); index++ ){

            auto microReferencePosition = microReferencePositions.find( *index );

            if ( microReferencePosition == microReferencePositions.end( ) ){

                return new errorNode( "computeDomainXis", "Micro node " + std::to_string( *index ) + " was not found in the micro reference positions map" );

            }

            auto microDisplacement = microDisplacements.find( *index );

            if ( microDisplacement == microDisplacements.end( ) ){

                return new errorNode( "computeDomainXis", "Micro node " + std::to_string( *index ) + " was not found in the micro displacements map" );

            }

            //Compute the relative position vector 
            domainXis.emplace( *index, ( microReferencePosition->second + microDisplacement->second ) - domainCM );

        }

        return NULL;

    }


    errorOut formMicroDomainToMacroProjectionMatrix( const uIntType &dim,
                                                     const uIntType nMicroNodes,
                                                     const uIntType nMacroNodes,
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
                                                     const std::unordered_map< uIntType, uIntType >* microNodeToLocalIndex,
                                                     const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex ){

        /*!
         * Form the micro to macro projection matrix due to the current domain
         *
         * :param const uIntType &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The global micro node IDs in the domain
         * :param const uIntVector &domainMacroNodeIndices: The global macro node IDs in the domain
         * :param const floatVector &microVolumes: The volumes of the micro-nodes ( global vector )
         * :param const floatVector &microDensities: The densities of the micro-nodes ( global vector )
         * :param const floatVector &microWeights: The weighting values of the micro-nodes ( global vector )
         * :param const floatVector &domainReferenceXiVectors: The relative position vectors of the nodes in the domain
         * :param const floatVector &domainInterpolationFunctionValues: The values of the macro interpolation functions
         *     at each of the micro nodes organized via [ N_11, N_12, ... ] where the first index is the micro node number
         *     and the second is the macro node number in the domain.
         * :param const floatVector &domainMacroNodeProjectedMass: The masses of the macro nodes derived from the
         *     projection of the micro masses in the domain.
         * :param const floatVector &domainMacroNodeProjectedMassMomentOfInertia: The mass weighted moment of inertia of the macro nodes
         *     derived from the projection of the micro masses and relative positions. [ I1_11, I1_12, I1_13, I1_21, ... ] in the domain
         * :param const floatVEctor &domainMacroNodeMassRelativePositionConstant: The mass-weighted relative position vector constant
         *     at each of the macro nodes derived from the projection of the micro masses and relative positions in the domain.
         *     [ C1_1, C1_2, C1_3, C2_1, ... ]
         * :param SparseMatrix &projector: The contribution of the current domain to the projector
         * :param const std::unordered_map< uIntType, uIntType >* microNodeToLocalIndex: The map from the global micro node IDs
         *     to the local micro node ids
         * :param const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex: The mapp from the global macro node IDs
         *     to the local macro node ids
         */

        //Error handling
        if ( dim != 3 ){
            return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                  "Only 3D domains are currently supported" );
        }

        if ( dim * domainMicroNodeIndices.size() != domainReferenceXiVectors.size() ){
            return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                  "The number of micro node indices is not equal to the number of Xi vectors" );
        }

        if ( microWeights.size( ) != microDensities.size( ) ){
            return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                  "The micro weight and micro density vectors are of inconsistent sizes" );
        }

        if ( microWeights.size( ) != microVolumes.size( ) ){
            return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                  "The micro weight and micro volume vectors are of inconsistent sizes" );
        }

        if ( !microNodeToLocalIndex ){
            if ( nMicroNodes != microWeights.size() ){
                return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                      "The number of micro nodes is not equal to the number of weights" );
            }
        }

        if ( dim * dim * domainMacroNodeProjectedMass.size( ) != domainMacroNodeProjectedMassMomentOfInertia.size( ) ){
            return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                  "The macro node projected mass and macro node projected mass moment of inertia vectors are not of consistent sizes" );
        }

        if ( dim * domainMacroNodeProjectedMass.size( ) != domainMacroNodeMassRelativePositionConstant.size( ) ){
            return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                  "The macro node projected mass and macro node mass weighted relative position constant vectors are not of consistent sizes" );
        }

        //Set the number of spatial degrees of freedom for the two scales
        const uIntType nMicroDOF = dim;
        const uIntType nMacroDOF = dim + dim * dim;

        //Set up the vector of terms for the sparse matrix
        std::vector< T > coefficients;
        coefficients.reserve( nMicroDOF * domainMicroNodeIndices.size() * nMacroDOF * domainMacroNodeIndices.size() );

        //Initialize the row and column indices for the sparse matrix
        uIntType row0, col0;

        //Initialize the macro node mass
        floatType macroNodeMass;

        //Initialize the macro inverse mass moment of inertia
        floatVector inverseMacroMassMomentOfInertia;

        //Initialize the macro mass weighted relative position constant
        floatVector C;

        //Initialize the micro mass
        floatType microMass;

        //Initialize the Xi vector
        floatVector Xi;

        //Initialize the value of the weight and shapefunction value
        floatType w, sf;

        //Initialize terms to simplify the expressions
        floatType weightedMassTerm;
        floatVector positionTerm;

        //Initialize the global micro node index
        uIntType m;

        //Initialize the global macro node index
        uIntType n;

        //Initialize the output index for the micro-nodes
        uIntType o;

        //Initialize the output index for the macro-nodes
        uIntType p;

        //Loop over the macro nodes
        for ( uIntType i = 0; i < domainMacroNodeIndices.size( ); i++ ){

            //Set the global macro node index
            n = domainMacroNodeIndices[ i ];

            //Set the column index
            if ( macroNodeToLocalIndex ){

                auto indx = macroNodeToLocalIndex->find( n );

                if ( indx == macroNodeToLocalIndex->end( ) ){

                    return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                          "The macro node " + std::to_string( n ) + " is not found in the macro node to local index map" );

                }

                p = indx->second;

            }
            else{
                p = domainMacroNodeIndices[ i ];
            }

            row0 = nMacroDOF * p;

            if ( i > domainMacroNodeProjectedMass.size( ) ){

                return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                      "The macro node " + std::to_string( n ) + " is too large for the macro node projected mass vector" );

            }

            //Set the macro node mass
            macroNodeMass = domainMacroNodeProjectedMass[ i ];

            //Set the inverse macro node mass moment of inertia
            inverseMacroMassMomentOfInertia
                = vectorTools::inverse( floatVector( domainMacroNodeProjectedMassMomentOfInertia.begin( ) + dim * dim * i,
                                                     domainMacroNodeProjectedMassMomentOfInertia.begin( ) + dim * dim * ( i + 1 ) ),
                                        dim, dim );

            //Set the mass weighted relative position constant
            C = floatVector ( domainMacroNodeMassRelativePositionConstant.begin( ) + dim * i,
                              domainMacroNodeMassRelativePositionConstant.begin( ) + dim * ( i + 1 ) );

            //Loop over the micro nodes

            for ( uIntType j = 0; j < domainMicroNodeIndices.size( ); j++ ){

                //Set the global micro node index
                m = domainMicroNodeIndices[ j ];

                if ( m >= microWeights.size( ) ){
                    return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                          "The number of micro-weights is smaller than required for micro-node " + std::to_string( m ) );
                }
    
                //Set the row index
                if ( microNodeToLocalIndex ){
    
                    auto indx = microNodeToLocalIndex->find( m );
    
                    if ( indx == microNodeToLocalIndex->end( ) ){
    
                        continue;
    
                    }
   
                    o = indx->second;
                }
                else{
                    o = domainMicroNodeIndices[ j ];
                }
    
                col0 = nMicroDOF * o;

                if ( m >= microWeights.size( ) ){
                    return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                          "The number of micro-weights is smaller than required for micro-node " + std::to_string( m ) );
                }

                if ( domainMacroNodeIndices.size( ) * j + i >= domainInterpolationFunctionValues.size( ) ){
                    return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                          "The number of micro shape functions in the domain is msaller than required for micro-node "
                                          + std::to_string( m ) + " and macro node " + std::to_string( n ) );
                }

                //Set the micro mass
                microMass = microDensities[ m ] * microVolumes[ m ];

                //Set the micro weight
                w = microWeights[ m ];

                //Set the shape function
                sf = domainInterpolationFunctionValues[ domainMacroNodeIndices.size( ) * j + i ];

                //Extract the micro-position
                Xi = floatVector( domainReferenceXiVectors.begin( ) + dim * j,
                                  domainReferenceXiVectors.begin( ) + dim * ( j + 1 ) );

                //Compute the weighted mass term
                weightedMassTerm = microMass * w * sf;

                //Compute the position term
                positionTerm = weightedMassTerm
                             * vectorTools::matrixMultiply( Xi - C / macroNodeMass, inverseMacroMassMomentOfInertia, 1, 3, 3, 3 );

                //Add the matrix coefficients
                
                //Macro displacements
                coefficients.push_back( T( row0 + 0, col0 + 0, ( weightedMassTerm / macroNodeMass ) ) );
                coefficients.push_back( T( row0 + 1, col0 + 1, ( weightedMassTerm / macroNodeMass ) ) );
                coefficients.push_back( T( row0 + 2, col0 + 2, ( weightedMassTerm / macroNodeMass ) ) );

                //Micro displacement ( phi )
                coefficients.push_back( T( row0 +  3, col0 + 0, positionTerm[ 0 ] ) ); 
                coefficients.push_back( T( row0 +  4, col0 + 0, positionTerm[ 1 ] ) ); 
                coefficients.push_back( T( row0 +  5, col0 + 0, positionTerm[ 2 ] ) ); 
                coefficients.push_back( T( row0 +  6, col0 + 1, positionTerm[ 0 ] ) ); 
                coefficients.push_back( T( row0 +  7, col0 + 1, positionTerm[ 1 ] ) ); 
                coefficients.push_back( T( row0 +  8, col0 + 1, positionTerm[ 2 ] ) ); 
                coefficients.push_back( T( row0 +  9, col0 + 2, positionTerm[ 0 ] ) ); 
                coefficients.push_back( T( row0 + 10, col0 + 2, positionTerm[ 1 ] ) ); 
                coefficients.push_back( T( row0 + 11, col0 + 2, positionTerm[ 2 ] ) );

            }

        }

        //Assemble the sparse matrix
        projector = SparseMatrix( nMacroDOF * nMacroNodes, nMicroDOF * nMicroNodes );
        projector.setFromTriplets( coefficients.begin(), coefficients.end() );

        return NULL;
    }

    errorOut formMicroDomainToMacroProjectionMatrix( const uIntType &dim,
                                                     const uIntType nMicroNodes,
                                                     const uIntType nMacroNodes,
                                                     const uIntVector  &domainMicroNodeIndices,
                                                     const uIntVector  &domainMacroNodeIndices,
                                                     const std::unordered_map< uIntType, floatType > &microVolumes,
                                                     const std::unordered_map< uIntType, floatType > &microDensities,
                                                     const std::unordered_map< uIntType, floatType > &microWeights,
                                                     const std::unordered_map< uIntType, floatVector > &domainReferenceXiVectors,
                                                     const std::unordered_map< uIntType, floatVector > &domainInterpolationFunctionValues,
                                                     const std::unordered_map< uIntType, floatType > &domainMacroNodeProjectedMass,
                                                     const std::unordered_map< uIntType, floatVector > &domainMacroNodeProjectedMassMomentOfInertia,
                                                     const std::unordered_map< uIntType, floatVector > &domainMacroNodeMassRelativePositionConstant,
                                                     SparseMatrix &projector,
                                                     const std::unordered_map< uIntType, uIntType >* microNodeToLocalIndex,
                                                     const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex ){

        /*!
         * Form the micro to macro projection matrix due to the current domain
         *
         * :param const uIntType &dim: The dimension of the problem
         * :param const uIntVector &domainMicroNodeIndices: The global micro node IDs in the domain
         * :param const uIntVector &domainMacroNodeIndices: The global macro node IDs in the domain
         * :param const std::unordered_map< uIntType, floatType > &microVolumes: The volumes of the micro-nodes
         * :param const std::unordered_map< uIntType, floatType > &microDensities: The densities of the micro-nodes
         * :param const std::unordered_map< uIntType, floatType > &microWeights: The weighting values of the micro-nodes
         * :param const std::unordered_map< uIntType, floatType > &domainReferenceXiVectors: The relative position 
         *     vectors of the nodes in the domain
         * :param const std::unordered_map< uIntType, floatVector > &domainInterpolationFunctionValues: The values of the 
         *     macro interpolation functions at each of the micro nodes organized via [ N_11, N_12, ... ] where the first 
         *     index is the micro node number and the second is the macro node number in the domain.
         * :param const std::unordered_map< uIntType, floatVector > &domainMacroNodeProjectedMass: The masses of the macro
         *     nodes derived from the projection of the micro masses in the domain.
         * :param const std::unordered_map< uIntType, floatVector > &domainMacroNodeProjectedMassMomentOfInertia: The mass
         *     weighted moment of inertia of the macro nodes derived from the projection of the micro masses and relative
         *     positions. [ I1_11, I1_12, I1_13, I1_21, ... ] in the domain
         * :param const std::unordered_map< uIntType, floatVector > &domainMacroNodeMassRelativePositionConstant: The 
         *     mass-weighted relative position vector constant at each of the macro nodes derived from the projection of 
         *     the micro masses and relative positions in the domain. [ C1_1, C1_2, C1_3, C2_1, ... ]
         * :param SparseMatrix &projector: The contribution of the current domain to the projector
         * :param const std::unordered_map< uIntType, uIntType >* microNodeToLocalIndex: The map from the global micro node IDs
         *     to the local micro node ids
         * :param const std::unordered_map< uIntType, uIntType >* macroNodeToLocalIndex: The mapp from the global macro node IDs
         *     to the local macro node ids
         */

        //Error handling
        if ( dim != 3 ){
            return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                  "Only 3D domains are currently supported" );
        }

        if ( !microNodeToLocalIndex ){
            if ( nMicroNodes != microWeights.size() ){
                return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                      "The number of micro nodes is not equal to the number of weights" );
            }
        }

        //Set the number of spatial degrees of freedom for the two scales
        const uIntType nMicroDOF = dim;
        const uIntType nMacroDOF = dim + dim * dim;

        //Set up the vector of terms for the sparse matrix
        std::vector< T > coefficients;
        coefficients.reserve( nMicroDOF * domainMicroNodeIndices.size() * nMacroDOF * domainMacroNodeIndices.size() );

        //Initialize the row and column indices for the sparse matrix
        uIntType row0, col0;

        //Initialize the macro inverse mass moment of inertia
        floatVector inverseMacroMassMomentOfInertia;

        //Initialize the micro mass
        floatType microMass;

        //Initialize the Xi vector
        floatVector Xi;

        //Initialize the value of the weight and shapefunction value
        floatType w, sf;

        //Initialize terms to simplify the expressions
        floatType weightedMassTerm;
        floatVector positionTerm;

        //Initialize the global micro node index
        uIntType m;

        //Initialize the global macro node index
        uIntType n;

        //Initialize the output index for the micro-nodes
        uIntType o;

        //Initialize the output index for the macro-nodes
        uIntType p;

        //Loop over the macro nodes
        for ( uIntType i = 0; i < domainMacroNodeIndices.size( ); i++ ){

            //Set the global macro node index
            n = domainMacroNodeIndices[ i ];

            //Set the column index
            if ( macroNodeToLocalIndex ){

                auto indx = macroNodeToLocalIndex->find( n );

                if ( indx == macroNodeToLocalIndex->end( ) ){

                    return new errorNode( "formMicroDomaintoMacroProjectionMatrix",
                                          "The macro node " + std::to_string( n ) + " is not found in the macro node to local index map" );

                }

                p = indx->second;

            }
            else{
                p = domainMacroNodeIndices[ i ];
            }

            row0 = nMacroDOF * p;

            if ( i > domainMacroNodeProjectedMass.size( ) ){

                return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                      "The macro node " + std::to_string( n ) + " is too large for the macro node projected mass vector" );

            }

            //Set the macro node mass
            auto macroNodeMass = domainMacroNodeProjectedMass.find( n );

            if ( macroNodeMass == domainMacroNodeProjectedMass.end( ) ){

                return new errorNode( "formMacroDomainToMacroProjectionMatrix",
                                      "The macro node " + std::to_string( n ) + " is not found in the macro node mass map" );

            }

            //Get the macro node mass moment of inertia
            auto macroNodeMassMomentOfInertia = domainMacroNodeProjectedMassMomentOfInertia.find( n );

            if ( macroNodeMassMomentOfInertia == domainMacroNodeProjectedMassMomentOfInertia.end( ) ){

                return new errorNode( "formMacroDomainToMacroProjectionMatrix",
                                      "The macro node " + std::to_string( n ) + " is not found in the macro node mass moment of inertia map" );

            }


            //Set the inverse macro node mass moment of inertia
            inverseMacroMassMomentOfInertia = vectorTools::inverse( macroNodeMassMomentOfInertia->second, dim, dim );

            //Set the mass weighted relative position constant
            auto C = domainMacroNodeMassRelativePositionConstant.find( n );

            if ( C == domainMacroNodeMassRelativePositionConstant.end( ) ){

                return new errorNode( "formMacroDomainToMacroProjectionMatrix",
                                      "The macro node " + std::to_string( n ) + " is not found in the macro node relative position constant" );

            }

            //Loop over the micro nodes

            for ( uIntType j = 0; j < domainMicroNodeIndices.size( ); j++ ){

                //Set the global micro node index
                m = domainMicroNodeIndices[ j ];

                auto microWeight = microWeights.find( m );
                if ( microWeight == microWeights.end( ) ){
                    return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                          "The micro node " + std::to_string( m ) + " was not found in the micro weight map" );
                }
    
                auto microVolume = microVolumes.find( m );
                if ( microVolume == microVolumes.end( ) ){
                    return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                          "The micro node " + std::to_string( m ) + " was not found in the micro volume map" );
                }
    
                auto microDensity = microDensities.find( m );
                if ( microDensity == microDensities.end( ) ){
                    return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                          "The micro node " + std::to_string( m ) + " was not found in the micro density map" );
                }

                auto domainReferenceXi = domainReferenceXiVectors.find( m );
                if ( domainReferenceXi == domainReferenceXiVectors.end( ) ){
                    return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                          "The micro node " + std::to_string( m ) + " was not found in the reference Xi vector map" );
                }

                auto shapefunctions = domainInterpolationFunctionValues.find( m );
                if ( shapefunctions == domainInterpolationFunctionValues.end( ) ){
                    return new errorNode( "formMicroDomainToMacroProjectionMatrix",
                                          "The micro node " + std::to_string( m ) + " was not found in the interpolation function map" );
                }
    
                //Set the row index
                if ( microNodeToLocalIndex ){
    
                    auto indx = microNodeToLocalIndex->find( m );
    
                    if ( indx == microNodeToLocalIndex->end( ) ){
    
                        continue;
    
                    }
   
                    o = indx->second;
                }
                else{
                    o = domainMicroNodeIndices[ j ];
                }
    
                col0 = nMicroDOF * o;

//                if ( domainMacroNodeIndices.size( ) * j + i >= domainInterpolationFunctionValues.size( ) ){
//                    return new errorNode( "formMicroDomainToMacroProjectionMatrix",
//                                          "The number of micro shape functions in the domain is smaller than required for micro-node "
//                                          + std::to_string( m ) + " and macro node " + std::to_string( n ) );
//                }

                //Set the micro mass
                microMass = microDensity->second * microVolume->second;

                //Set the micro weight
                w = microWeight->second;

                //Set the shape function
                sf = shapefunctions->second[ i ];

                //Extract the micro-position
                Xi = domainReferenceXi->second;

                //Compute the weighted mass term
                weightedMassTerm = microMass * w * sf;

                //Compute the position term
                positionTerm = weightedMassTerm
                             * vectorTools::matrixMultiply( Xi - C->second / macroNodeMass->second, inverseMacroMassMomentOfInertia, 1, 3, 3, 3 );

                //Add the matrix coefficients
                
                //Macro displacements
                coefficients.push_back( T( row0 + 0, col0 + 0, ( weightedMassTerm / macroNodeMass->second ) ) );
                coefficients.push_back( T( row0 + 1, col0 + 1, ( weightedMassTerm / macroNodeMass->second ) ) );
                coefficients.push_back( T( row0 + 2, col0 + 2, ( weightedMassTerm / macroNodeMass->second ) ) );

                //Micro displacement ( phi )
                coefficients.push_back( T( row0 +  3, col0 + 0, positionTerm[ 0 ] ) ); 
                coefficients.push_back( T( row0 +  4, col0 + 0, positionTerm[ 1 ] ) ); 
                coefficients.push_back( T( row0 +  5, col0 + 0, positionTerm[ 2 ] ) ); 
                coefficients.push_back( T( row0 +  6, col0 + 1, positionTerm[ 0 ] ) ); 
                coefficients.push_back( T( row0 +  7, col0 + 1, positionTerm[ 1 ] ) ); 
                coefficients.push_back( T( row0 +  8, col0 + 1, positionTerm[ 2 ] ) ); 
                coefficients.push_back( T( row0 +  9, col0 + 2, positionTerm[ 0 ] ) ); 
                coefficients.push_back( T( row0 + 10, col0 + 2, positionTerm[ 1 ] ) ); 
                coefficients.push_back( T( row0 + 11, col0 + 2, positionTerm[ 2 ] ) );

            }

        }

        //Assemble the sparse matrix
        projector = SparseMatrix( nMacroDOF * nMacroNodes, nMicroDOF * nMicroNodes );
        projector.setFromTriplets( coefficients.begin(), coefficients.end() );

        return NULL;
    }
}
