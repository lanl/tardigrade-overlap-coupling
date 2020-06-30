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

namespace overlapCoupling{

    /*=========================================================================
    |                                Typedefs                                 |
    =========================================================================*/

    typedef inputFileProcessor::errorNode errorNode;
    typedef inputFileProcessor::errorOut errorOut;
    typedef inputFileProcessor::floatType floatType;
    typedef inputFileProcessor::floatVector floatVector;
    typedef inputFileProcessor::floatMatrix floatMatrix;
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

            const errorOut getConstructorError( );

            errorOut initializeCoupling( );

            errorOut processIncrement( const unsigned int &increment );

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

            //Private functions
            errorOut processDomainMassData( const unsigned int &increment, const std::string &domainName,
                                            floatType &domainMass, floatVector &domainCenterOfMass,
                                            floatVector &domainXiVectors );

            //Compute initial values
            errorOut setReferenceStateFromIncrement( const unsigned int &increment );

            //Compute the increment's values
            errorOut computeIncrementCentersOfMass( const unsigned int increment,
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
                                                            const unsigned int &increment,
                                                            const floatVector &domainCenterOfMass,
                                                            floatVector &domainCenterOfMassShapeFunctionValues,
                                                            floatVector &domainMicroPositionShapeFunctionValues );

            errorOut addDomainContributionToInterpolationMatrix( const uIntVector  &domainNodes,
                                                                 const uIntVector  &macroNodes,
                                                                 const floatVector &domainReferenceXis,
                                                                 const floatVector &domainCenterOfMassShapeFunctionValues );
                                                                   
            //Construct the interpolation matrix
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

    };

}

#endif
