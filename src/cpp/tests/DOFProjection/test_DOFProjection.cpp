//!The test file for DOFProjection.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>

#include<DOFProjection.h>

typedef DOFProjection::errorNode errorNode; //!Redefinition for the error node
typedef DOFProjection::errorOut errorOut; //!Redefinition for a pointer to the error node
typedef DOFProjection::floatType floatType; //!Define the float values type.
typedef DOFProjection::floatVector floatVector; //! Define a vector of floats
typedef DOFProjection::floatMatrix floatMatrix; //!Define a matrix of floats
typedef DOFProjection::uIntVector uIntVector; //!Define a vector of unsigned ints


int test_addMacroDomainDisplacementToMicro( std::ofstream &results ){
    /*!
     * Test the projection of the macro-domain's displacement to the 
     * micro-scale when u and phi at a given point are known.
     */

    const unsigned int dim = 3;

    const unsigned int nMicroNodes = 100;

    const uIntVector domainMicroNodeIndices = { 53, 28, 63, 97, 93, 90,  8,  5,  0, 62 };

    const floatVector u = { 0.4802733 , 0.63413557, 0.47580155 };

    const floatVector phi = { 0.24395441, 0.46860497, 0.43078742,
                              0.61868352, 0.46794329, 0.66017423,
                              0.58630018, 0.55379286, 0.50449636 };

    const floatVector referenceXis = { -0.02920635,  0.39712726, -0.83686303,  0.73820473, -0.13378864,
                                       -0.01133987, -0.00851906, -0.25855584,  0.84425732,  0.68255644,
                                        0.31105184, -0.0746299 ,  0.13002262, -0.15216899, -0.42357609,
                                       -0.47203856,  0.38147646,  0.66567306,  0.18478316, -0.06695484,
                                        0.34731997, -0.21704129,  0.40420874,  0.93139529, -0.81898025,
                                       -0.88374973,  0.17603484,  0.50234751,  0.02263478,  0.84100238 };

    const floatVector domainMicroWeights = { 0.3039641 , 0.49300273, 0.97936034, 0.32350827, 0.18956717,
                                             0.30522911, 0.34411193, 0.67953029, 0.053815  , 0.80660376 };

    const floatVector answer = { -0.00311139, -0.00914237, -0.02179357,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.69174329,  0.88602979,  0.70826331,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.22146963,
                                  0.32567389,  0.24854602,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.2922434 ,
                                  0.50323692,  0.40860017,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.09080442,
                                  0.07581585,  0.07793936,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.78702201,  1.21856017,  0.97368818,  0.7058543 ,
                                  1.0432464 ,  0.73799228,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.25353635,  0.29303907,  0.22774176,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.04894899,
                                  0.06895278,  0.04816355,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                  0.        ,  0.24599464,  0.37291053,  0.32693496,  0.        ,
                                  0.        ,  0.        ,  0.        ,  0.        ,  0.        };

    floatVector result( dim * nMicroNodes );

    errorOut error = DOFProjection::addMacroDomainDisplacementToMicro( dim, domainMicroNodeIndices, u, phi, referenceXis,
                                                                       domainMicroWeights, result );

    if ( error ){
        error->print();
        results << "test_addMacroDomainDisplacementToMicro & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_addMacroDomainDisplacementToMicro (test 1) & False\n";
        return 1;
    }

    results << "test_addMicroDomainDisplacementToMicro & True\n";
    return 0;

}

int main(){
    /*!
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or False 
    if the test passes or fails respectively.
    */

    //Open the results file
    std::ofstream results;
    results.open("results.tex");

    test_addMacroDomainDisplacementToMicro( results );

    //Close the results file
    results.close();
}
