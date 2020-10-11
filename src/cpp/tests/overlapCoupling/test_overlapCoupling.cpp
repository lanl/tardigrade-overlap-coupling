//!The test file for dataFileInterface.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>
#include <boost/algorithm/string.hpp>
#define USE_EIGEN
#include<vector_tools.h>

#include<overlapCoupling.h>
#include<generateXDMFData.h>

typedef overlapCoupling::errorNode errorNode; //!Redefinition for the error node
typedef overlapCoupling::errorOut errorOut; //!Redefinition for a pointer to the error node
typedef overlapCoupling::uIntType uIntType; //!Define a type for an unsigned int
typedef overlapCoupling::uIntVector uIntVector; //!Define a type for a vector of unsigned ints
typedef overlapCoupling::floatType floatType; //!Define the float values type.
typedef overlapCoupling::floatVector floatVector; //! Define a vector of floats
typedef overlapCoupling::floatMatrix floatMatrix; //!Define a matrix of floats
typedef overlapCoupling::uIntVector uIntVector; //!Define a vector of unsigned ints
typedef overlapCoupling::SparseMatrix SparseMatrix; //!Define a sparse matrix
typedef DOFProjection::T T; //!Define the triplet

typedef overlapCoupling::domainFloatMap domainFloatMap; //!Define a map from a micro domain to a float property
typedef overlapCoupling::domainFloatVectorMap domainFloatVectorMap; //!Define a map from a micro domain to a float vector property
typedef overlapCoupling::cellDomainUIntVectorMap cellDomainUIntVectorMap; //!Define a map from a macro cell to a micro domain intvector property
typedef overlapCoupling::cellDomainFloatMap cellDomainFloatMap; //!Define a map from a macro cell to a micro domain float property
typedef overlapCoupling::cellDomainFloatVectorMap cellDomainFloatVectorMap; //!Define a map from a macro cell to a micro domain float vector property

errorOut _createXDMFDatafiles( ){

    fileGenerator::fileGenerator fg( "macroscale.yaml" );
    if ( fg.build( ) ){

        fg.getError( )->print( );
        errorOut result = new errorNode( "_createXDMFDatafiles", "Error in creation of the macroscale datafile" );
        return result;

    }

    fg = fileGenerator::fileGenerator( "microscale.yaml" );
    if ( fg.build( ) ){

        fg.getError( )->print( );
        errorOut result = new errorNode( "_createXDMFDatafiles", "Error in creation of the microscale datafile" );
        return result;

    }

    return NULL;

}

errorOut readMatrixFromFile( const std::string filename, floatVector &data, Eigen::MatrixXd &matrix ){
    /*!
     * Read a matrix of doubles from a file
     */

    std::ifstream file( filename );
    std::string line;
    std::vector< std::string > splitLine;

    unsigned int rows = 0;
    unsigned int cols = 0;

    if ( file.is_open( ) ){

        while ( std::getline( file, line ) ){


            boost::algorithm::split( splitLine, line, boost::is_any_of( "," ) );
    
            if ( ( splitLine.size( ) != rows ) && ( rows != 0 ) ){
    
                return new errorNode( "readMatrixFromFile", "The matrix is not a consistent matrix" );
    
            }
            else{

                rows = splitLine.size( );

            }
    
            for ( auto v = splitLine.begin( ); v != splitLine.end( ); v++ ){
    
                data.push_back( std::stod( *v ) );
    
            }
    
            cols++;
    
        }

    }
    else{

        return new errorNode( "readMatrixFromFile", "Can't open " + filename );

    }

    if ( cols == 0 ){

        return new errorNode( "readMatrixFromFile", "There are no columns in the matrix" );

    }

    if ( rows == 0 ){

        return new errorNode( "readMatrixFromFile", "there are no rows in the matrix" );

    }

    matrix = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>( data.data(), rows, cols );

    return NULL;

}

template< class T1, class T2 >
int _compare_domainMaps( std::ofstream &results,
                         const std::unordered_map< T1, T2 > &answer,
                         const std::unordered_map< T1, T2 > &result,
                         const std::string testName,
                         uIntType &testNum, const floatType tolr = 1e-6,
                         const floatType tola = 1e-6 ){
    /*!
     * Compare domain maps to eachother
     *
     * :param std::ofstream &results: The output file
     * :param const std::unordered_map< T1, T2 > &answer: The answer map
     * :param const std::unordered_map< T1, T2 > &result: The result map
     * :param const std::string testName: The name of the test
     * :param uIntType &testNum: The number of the test
     * :param floatType &tolr: The relative tolerance
     * :param flaotType &tola: The absolute tolerance
     */

    for ( auto a_domain = answer.begin( ); a_domain != answer.end( ); a_domain++ ){

        auto r_domain = result.find( a_domain->first );

        if ( r_domain == result.end( ) ){

            std::string outstr = "test_";
            outstr += testName;
            outstr += " (test ";
            outstr += std::to_string( testNum + 1 );
            outstr += ") & False\n";
            results << outstr;
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( r_domain->second, a_domain->second, tolr, tola ) ){

            std::cout << r_domain->first << "\n";
            std::cout << a_domain->first << "\n";
            std::string outstr = "test_";
            outstr += testName;
            outstr += " (test ";
            outstr += std::to_string( testNum + 2 );
            outstr += ") & False\n";
            results << outstr;
            return 1;

        }

    }

    testNum += 2;

    return 0;

}

template< class T >
int _compare_cellDomainMaps( std::ofstream &results,
                             const std::unordered_map< uIntType, std::unordered_map< std::string, T > > &answer,
                             const std::unordered_map< uIntType, std::unordered_map< std::string, T > > &result,
                             const std::string testName,
                             uIntType &testNum,
                             const floatType tolr = 1e-6, const floatType tola = 1e-6 ){
    /*!
     * Compare cell domain maps to eachother
     *
     * :param std::ofstream &results: The output file
     * :param const std::unordered_map< uIntType, std::unordered_map< std::string, T > > &answer: The answer map
     * :param const std::unordered_map< uIntType, std::unordered_map< std::string, T > > &result: The result map
     * :param const std::string testName: The name of the test
     * :param uIntType &testNum: The number of the test
     * :param floatType &tolr: The relative tolerance
     * :param floatType &tola: The absolute tolerance
     */

    uIntType tmp;
    for ( auto a_cell = answer.begin( ); a_cell != answer.end( ); a_cell++ ){

        auto r_cell = result.find( a_cell->first );

        if ( r_cell == result.end( ) ){

            std::string outstr = "test_";
            outstr += testName;
            outstr += " (test ";
            outstr += std::to_string( testNum + 1 );
            outstr += ") & False\n";
            results << outstr;
            return 1;

        }

        tmp = testNum + 1;
        if ( _compare_domainMaps( results, a_cell->second, r_cell->second, testName, tmp, tolr, tola ) ){

            return 1;

        }

    }

    testNum = tmp;

    return 0;

}

template< class T >
int _compare_cellDomainPointMaps( std::ofstream &results,
                                  const std::unordered_map< uIntType, std::unordered_map< std::string, std::unordered_map< uIntType, T > > > &answer,
                                  const std::unordered_map< uIntType, std::unordered_map< std::string, std::unordered_map< uIntType, T > > > &result,
                                  const std::string testName,
                                  uIntType &testNum ){
    /*!
     * Compare cell domain maps of points to eachother
     *
     * :param std::ofstream &results: The output file
     * :param const std::unordered_map< uIntType, std::unordered_map< std::string, std::unordered_map< uIntType, T > > > &answer: The answer map
     * :param const std::unordered_map< uIntType, std::unordered_map< std::string, std::unordered_map< uIntType, T > > > &result: The result map
     * :param const std::string testName: The name of the test
     * :param uIntType &testNum: The number of the test
     */

    for ( auto a_cell = answer.begin( ); a_cell != answer.end( ); a_cell++ ){

        auto r_cell = result.find( a_cell->first );

        if ( r_cell == result.end( ) ){

            std::string outstr = "test_";
            outstr += testName;
            outstr += " (test ";
            outstr += std::to_string( testNum + 1 );
            outstr += ") & False\n";
            results << outstr;
            return 1;

        }

        for ( auto a_domain = a_cell->second.begin( ); a_domain != a_cell->second.end( ); a_domain++ ){

            auto r_domain = r_cell->second.find( a_domain->first );

            if ( r_domain == r_cell->second.end( ) ){

                std::string outstr = "test_";
                outstr += testName;
                outstr += " (test ";
                outstr += std::to_string( testNum + 2 );
                outstr += ") & False\n";
                results << outstr;

            }

            for ( auto a_point = a_domain->second.begin( ); a_point != a_domain->second.end( ); a_point++ ){

                auto r_point = r_domain->second.find( a_point->first );

                if ( r_point == r_domain->second.end( ) ){

                    std::string outstr = "test_";
                    outstr += testName;
                    outstr += " (test ";
                    outstr += std::to_string( testNum + 3 );
                    outstr += ") & False\n";
                    results << outstr;

                }

                if ( !vectorTools::fuzzyEquals( r_point->second, a_point->second ) ){

                    std::string outstr = "test_";
                    outstr += testName;
                    outstr += " (test ";
                    outstr += std::to_string( testNum + 4 );
                    outstr += ") & False\n";
                    results << outstr;
                    return 1;

                }

            }

        }

    }

    testNum += 4;

    return 0;

}


int test_overlapCoupling_constructor( std::ostream &results ){
    /*!
     * Test the constructor to make sure that the code generates properly
     *
     * :param std::ofstream &results: The output file
     */

    std::string filename = "testConfig_averaged_l2_projection.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_constructor & False\n";
        return 1;
    }

    results << "test_overlapCoupling_constructor & True\n";
    return 0;
}

int test_overlapCoupling_initializeCoupling_l2_projection( std::ofstream &results ){
    /*!
     * Test the initialization of the coupling for the l2_projection
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string filename = "testConfig_l2_projection.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_initializeCoupling_l2_projection & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_initializeCoupling_l2_projection & False\n";
        return 1;
    }

    if ( !std::ifstream( "reference_information.xdmf" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_l2_projection (test 1) & False\n";
        return 1;

    }

    if ( !std::ifstream( "reference_information.h5" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_l2_projection (test 2) & False\n";
        return 1;

    }

    if ( !std::ifstream( "homogenized_response.xdmf" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_l2_projection (test 3) & False\n";
        return 1;

    }

    if ( !std::ifstream( "homogenized_response.h5" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_l2_projection (test 4) & False\n";
        return 1;

    }

    std::string testName = "overlapCoupling_initializeCoupling";
    uIntType testNum = 4;

    cellDomainFloatMap domainMassAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1" , 0.25 },
                    { "ghost_nodeset_volume_2" , 0.25 },
                    { "ghost_nodeset_volume_3" , 0.25 },
                    { "ghost_nodeset_volume_4" , 0.25 },
                    { "ghost_nodeset_volume_5" , 0.25 },
                    { "ghost_nodeset_volume_6" , 0.25 },
                    { "ghost_nodeset_volume_7" , 0.25 },
                    { "ghost_nodeset_volume_8" , 0.25 }
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1" , 0.25 },
                    { "free_nodeset_volume_2" , 0.25 },
                    { "free_nodeset_volume_3" , 0.25 },
                    { "free_nodeset_volume_4" , 0.25 },
                    { "free_nodeset_volume_5" , 0.375 },
                    { "free_nodeset_volume_6" , 0.375 },
                    { "free_nodeset_volume_7" , 0.375 },
                    { "free_nodeset_volume_8" , 0.375 }
                }
            }
        };

    if ( _compare_cellDomainMaps( results, domainMassAnswer, oc._test_domainMass,
                                  testName, testNum ) ){
        return 1;
    }

    cellDomainFloatVectorMap domainCOMAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { 0.250000, 0.250000, 0.250000 } },
                    { "ghost_nodeset_volume_2", { 0.750000, 0.250000, 0.250000 } },
                    { "ghost_nodeset_volume_3", { 0.750000, 0.750000, 0.250000 } },
                    { "ghost_nodeset_volume_4", { 0.250000, 0.750000, 0.250000 } },
                    { "ghost_nodeset_volume_5", { 0.250000, 0.250000, 0.750000 } },
                    { "ghost_nodeset_volume_6", { 0.750000, 0.250000, 0.750000 } },
                    { "ghost_nodeset_volume_7", { 0.750000, 0.750000, 0.750000 } },
                    { "ghost_nodeset_volume_8", { 0.250000, 0.750000, 0.750000 } }
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { 0.250000, 0.250000, 1.250000 } },
                    { "free_nodeset_volume_2", { 0.750000, 0.250000, 1.250000 } },
                    { "free_nodeset_volume_3", { 0.750000, 0.750000, 1.250000 } },
                    { "free_nodeset_volume_4", { 0.250000, 0.750000, 1.250000 } },
                    { "free_nodeset_volume_5", { 0.250000, 0.250000, 1.833333 } },
                    { "free_nodeset_volume_6", { 0.750000, 0.250000, 1.833333 } },
                    { "free_nodeset_volume_7", { 0.750000, 0.750000, 1.833333 } },
                    { "free_nodeset_volume_8", { 0.250000, 0.750000, 1.833333 } }
                }
            }
        };

    if ( _compare_cellDomainMaps( results, domainCOMAnswer, oc._test_domainCOM,
                                  testName, testNum ) ){
        return 1;
    }

    cellDomainFloatVectorMap freeDomainMomentsOfInertiaAnswer
        =
        {
            { 2,
                {
                    { "free_nodeset_volume_1", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "free_nodeset_volume_2", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "free_nodeset_volume_3", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "free_nodeset_volume_4", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "free_nodeset_volume_5", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.055556 } },
                    { "free_nodeset_volume_6", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.055556 } },
                    { "free_nodeset_volume_7", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.055556 } },
                    { "free_nodeset_volume_8", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.055556 } }
                }
            }
        };

    if ( _compare_cellDomainMaps( results, freeDomainMomentsOfInertiaAnswer, *oc.getReferenceFreeMicroDomainMomentsOfInertia( ),
                                  testName, testNum ) ){
        return 1;
    }

    cellDomainFloatVectorMap ghostDomainMomentsOfInertiaAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_2", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_3", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_4", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_5", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_6", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_7", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_8", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } }
                }
            },
        };

    if ( _compare_cellDomainMaps( results, ghostDomainMomentsOfInertiaAnswer, *oc.getReferenceGhostMicroDomainMomentsOfInertia( ),
                                  testName, testNum ) ){
        return 1;
    }

    std::unordered_map< uIntType, std::unordered_map< std::string, std::unordered_map< uIntType, floatVector > > > domainXiAnswer
        =
        {
            {
                1,
                {
                    {
                        "ghost_nodeset_volume_1",
                        {
                            { 24, { -0.250000, -0.250000, -0.250000,  } },
                            { 39, { -0.250000, -0.250000, +0.250000,  } },
                            { 40, { +0.250000, -0.250000, -0.250000,  } },
                            { 57, { +0.250000, -0.250000, +0.250000,  } },
                            { 44, { -0.250000, +0.250000, -0.250000,  } },
                            { 58, { -0.250000, +0.250000, +0.250000,  } },
                            { 29, { +0.250000, +0.250000, -0.250000,  } },
                            { 59, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_2",
                        {
                            { 40, { -0.250000, -0.250000, -0.250000,  } },
                            { 57, { -0.250000, -0.250000, +0.250000,  } },
                            { 11, { +0.250000, -0.250000, -0.250000,  } },
                            {  0, { +0.250000, -0.250000, +0.250000,  } },
                            { 29, { -0.250000, +0.250000, -0.250000,  } },
                            { 59, { -0.250000, +0.250000, +0.250000,  } },
                            { 20, { +0.250000, +0.250000, -0.250000,  } },
                            { 60, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_3",
                        {
                            { 29, { -0.250000, -0.250000, -0.250000,  } },
                            { 59, { -0.250000, -0.250000, +0.250000,  } },
                            { 20, { +0.250000, -0.250000, -0.250000,  } },
                            { 60, { +0.250000, -0.250000, +0.250000,  } },
                            { 47, { -0.250000, +0.250000, -0.250000,  } },
                            { 49, { -0.250000, +0.250000, +0.250000,  } },
                            { 17, { +0.250000, +0.250000, -0.250000,  } },
                            { 38, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_4",
                        {
                            { 44, { -0.250000, -0.250000, -0.250000,  } },
                            { 58, { -0.250000, -0.250000, +0.250000,  } },
                            { 29, { +0.250000, -0.250000, -0.250000,  } },
                            { 59, { +0.250000, -0.250000, +0.250000,  } },
                            { 14, { -0.250000, +0.250000, -0.250000,  } },
                            { 55, { -0.250000, +0.250000, +0.250000,  } },
                            { 47, { +0.250000, +0.250000, -0.250000,  } },
                            { 49, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_5",
                        {
                            { 39, { -0.250000, -0.250000, -0.250000,  } },
                            { 15, { -0.250000, -0.250000, +0.250000,  } },
                            { 57, { +0.250000, -0.250000, -0.250000,  } },
                            { 13, { +0.250000, -0.250000, +0.250000,  } },
                            { 58, { -0.250000, +0.250000, -0.250000,  } },
                            { 53, { -0.250000, +0.250000, +0.250000,  } },
                            { 59, { +0.250000, +0.250000, -0.250000,  } },
                            { 37, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_6",
                        {
                            { 57, { -0.250000, -0.250000, -0.250000,  } },
                            { 13, { -0.250000, -0.250000, +0.250000,  } },
                            {  0, { +0.250000, -0.250000, -0.250000,  } },
                            {  5, { +0.250000, -0.250000, +0.250000,  } },
                            { 59, { -0.250000, +0.250000, -0.250000,  } },
                            { 37, { -0.250000, +0.250000, +0.250000,  } },
                            { 60, { +0.250000, +0.250000, -0.250000,  } },
                            {  3, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_7",
                        {
                            { 59, { -0.250000, -0.250000, -0.250000,  } },
                            { 37, { -0.250000, -0.250000, +0.250000,  } },
                            { 60, { +0.250000, -0.250000, -0.250000,  } },
                            {  3, { +0.250000, -0.250000, +0.250000,  } },
                            { 49, { -0.250000, +0.250000, -0.250000,  } },
                            { 32, { -0.250000, +0.250000, +0.250000,  } },
                            { 38, { +0.250000, +0.250000, -0.250000,  } },
                            { 34, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_8",
                        {
                            { 58, { -0.250000, -0.250000, -0.250000,  } },
                            { 53, { -0.250000, -0.250000, +0.250000,  } },
                            { 59, { +0.250000, -0.250000, -0.250000,  } },
                            { 37, { +0.250000, -0.250000, +0.250000,  } },
                            { 55, { -0.250000, +0.250000, -0.250000,  } },
                            { 25, { -0.250000, +0.250000, +0.250000,  } },
                            { 49, { +0.250000, +0.250000, -0.250000,  } },
                            { 32, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                }
            },
            {
                2,
                {
                    {
                        "free_nodeset_volume_1",
                        {
                            { 15, { -0.250000, -0.250000, -0.250000,  } },
                            { 31, { -0.250000, -0.250000, +0.250000,  } },
                            { 13, { +0.250000, -0.250000, -0.250000,  } },
                            { 26, { +0.250000, -0.250000, +0.250000,  } },
                            { 53, { -0.250000, +0.250000, -0.250000,  } },
                            { 21, { -0.250000, +0.250000, +0.250000,  } },
                            { 37, { +0.250000, +0.250000, -0.250000,  } },
                            { 48, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_2",
                        {
                            { 13, { -0.250000, -0.250000, -0.250000,  } },
                            { 26, { -0.250000, -0.250000, +0.250000,  } },
                            {  5, { +0.250000, -0.250000, -0.250000,  } },
                            { 10, { +0.250000, -0.250000, +0.250000,  } },
                            { 37, { -0.250000, +0.250000, -0.250000,  } },
                            { 48, { -0.250000, +0.250000, +0.250000,  } },
                            {  3, { +0.250000, +0.250000, -0.250000,  } },
                            {  4, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_3",
                        {
                            { 37, { -0.250000, -0.250000, -0.250000,  } },
                            { 48, { -0.250000, -0.250000, +0.250000,  } },
                            {  3, { +0.250000, -0.250000, -0.250000,  } },
                            {  4, { +0.250000, -0.250000, +0.250000,  } },
                            { 32, { -0.250000, +0.250000, -0.250000,  } },
                            { 33, { -0.250000, +0.250000, +0.250000,  } },
                            { 34, { +0.250000, +0.250000, -0.250000,  } },
                            { 28, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_4",
                        {
                            { 53, { -0.250000, -0.250000, -0.250000,  } },
                            { 21, { -0.250000, -0.250000, +0.250000,  } },
                            { 37, { +0.250000, -0.250000, -0.250000,  } },
                            { 48, { +0.250000, -0.250000, +0.250000,  } },
                            { 25, { -0.250000, +0.250000, -0.250000,  } },
                            { 50, { -0.250000, +0.250000, +0.250000,  } },
                            { 32, { +0.250000, +0.250000, -0.250000,  } },
                            { 33, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_5",
                        {
                            { 31, { -0.250000, -0.250000, -0.333333,  } },
                            { 43, { -0.250000, -0.250000, +0.166667,  } },
                            { 26, { +0.250000, -0.250000, -0.333333,  } },
                            { 27, { +0.250000, -0.250000, +0.166667,  } },
                            { 21, { -0.250000, +0.250000, -0.333333,  } },
                            {  1, { -0.250000, +0.250000, +0.166667,  } },
                            { 48, { +0.250000, +0.250000, -0.333333,  } },
                            {  7, { +0.250000, +0.250000, +0.166667,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_6",
                        {
                            { 26, { -0.250000, -0.250000, -0.333333,  } },
                            { 27, { -0.250000, -0.250000, +0.166667,  } },
                            { 10, { +0.250000, -0.250000, -0.333333,  } },
                            { 30, { +0.250000, -0.250000, +0.166667,  } },
                            { 48, { -0.250000, +0.250000, -0.333333,  } },
                            {  7, { -0.250000, +0.250000, +0.166667,  } },
                            {  4, { +0.250000, +0.250000, -0.333333,  } },
                            { 16, { +0.250000, +0.250000, +0.166667,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_7",
                        {
                            { 48, { -0.250000, -0.250000, -0.333333,  } },
                            {  7, { -0.250000, -0.250000, +0.166667,  } },
                            {  4, { +0.250000, -0.250000, -0.333333,  } },
                            { 16, { +0.250000, -0.250000, +0.166667,  } },
                            { 33, { -0.250000, +0.250000, -0.333333,  } },
                            { 22, { -0.250000, +0.250000, +0.166667,  } },
                            { 28, { +0.250000, +0.250000, -0.333333,  } },
                            {  2, { +0.250000, +0.250000, +0.166667,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_8",
                        {
                            { 21, { -0.250000, -0.250000, -0.333333,  } },
                            {  1, { -0.250000, -0.250000, +0.166667,  } },
                            { 48, { +0.250000, -0.250000, -0.333333,  } },
                            {  7, { +0.250000, -0.250000, +0.166667,  } },
                            { 50, { -0.250000, +0.250000, -0.333333,  } },
                            { 46, { -0.250000, +0.250000, +0.166667,  } },
                            { 33, { +0.250000, +0.250000, -0.333333,  } },
                            { 22, { +0.250000, +0.250000, +0.166667,  } },
                        }
                    },
                }
            }
        };

    if ( _compare_cellDomainPointMaps( results, domainXiAnswer, oc._test_domainXi,
                                       testName, testNum ) ){
        return 1;
    }

//    std::cout << "shape functions\n";
    cellDomainFloatVectorMap domainCOMSFAnswer
        =
        {
            {
                1,
                {
                    {
                        "ghost_nodeset_volume_1",
                        { +0.421875, +0.140625, +0.046875, +0.140625, +0.140625, +0.046875, +0.015625, +0.046875 }
                    },
                    {
                        "ghost_nodeset_volume_2",
                        { +0.140625, +0.421875, +0.140625, +0.046875, +0.046875, +0.140625, +0.046875, +0.015625 }
                    },
                    {
                        "ghost_nodeset_volume_3",
                        { +0.046875, +0.140625, +0.421875, +0.140625, +0.015625, +0.046875, +0.140625, +0.046875 }
                    },
                    {
                        "ghost_nodeset_volume_4",
                        { +0.140625, +0.046875, +0.140625, +0.421875, +0.046875, +0.015625, +0.046875, +0.140625 }
                    },
                    {
                        "ghost_nodeset_volume_5",
                        { +0.140625, +0.046875, +0.015625, +0.046875, +0.421875, +0.140625, +0.046875, +0.140625 }
                    },
                    {
                        "ghost_nodeset_volume_6",
                        { +0.046875, +0.140625, +0.046875, +0.015625, +0.140625, +0.421875, +0.140625, +0.046875 }
                    },
                    {
                        "ghost_nodeset_volume_7",
                        { +0.015625, +0.046875, +0.140625, +0.046875, +0.046875, +0.140625, +0.421875, +0.140625 }
                    },
                    {
                        "ghost_nodeset_volume_8",
                        { +0.046875, +0.015625, +0.046875, +0.140625, +0.140625, +0.046875, +0.140625, +0.421875 }
                    },
                }
            },
            {
                2,
                {
                    {
                        "free_nodeset_volume_1",
                        { +0.421875, +0.140625, +0.046875, +0.140625, +0.140625, +0.046875, +0.015625, +0.046875 }
                    },
                    {
                        "free_nodeset_volume_2",
                        { +0.140625, +0.421875, +0.140625, +0.046875, +0.046875, +0.140625, +0.046875, +0.015625 }
                    },
                    {
                        "free_nodeset_volume_3",
                        { +0.046875, +0.140625, +0.421875, +0.140625, +0.015625, +0.046875, +0.140625, +0.046875 }
                    },
                    {
                        "free_nodeset_volume_4",
                        { +0.140625, +0.046875, +0.140625, +0.421875, +0.046875, +0.015625, +0.046875, +0.140625 }
                    },
                    {
                        "free_nodeset_volume_5",
                        { +0.093750, +0.031250, +0.010417, +0.031250, +0.468750, +0.156250, +0.052083, +0.156250 }
                    },
                    {
                        "free_nodeset_volume_6",
                        { +0.031250, +0.093750, +0.031250, +0.010417, +0.156250, +0.468750, +0.156250, +0.052083 }
                    },
                    {
                        "free_nodeset_volume_7",
                        { +0.010417, +0.031250, +0.093750, +0.031250, +0.052083, +0.156250, +0.468750, +0.156250 }
                    },
                    {
                        "free_nodeset_volume_8",
                        { +0.031250, +0.010417, +0.031250, +0.093750, +0.156250, +0.052083, +0.156250, +0.468750 }
                    },
                }
            }
        };

    if ( _compare_cellDomainMaps( results, domainCOMSFAnswer, *oc.getReferenceCellDomainCenterOfMassShapeFunctions( ),
                                  testName, testNum ) ){
        return 1;
    }

    std::string xdmf_filename = "reference_information.xdmf";
    shared_ptr< XdmfReader > reader = XdmfReader::New( );
    shared_ptr< XdmfDomain > _readDomain = shared_dynamic_cast< XdmfDomain >( reader->read( xdmf_filename ) );
    shared_ptr< XdmfUnstructuredGrid > _readGrid = _readDomain->getUnstructuredGrid( 0 );

    SparseMatrix N;
    overlapCoupling::readSparseMatrixFromXDMF( _readGrid, "N", N );

    Eigen::MatrixXd A( N.rows( ), 1 );
    A << -0.416617, -0.311112, +0.013488, -0.337511, -0.627147, +0.058891,
         -0.307855, -0.308153, -0.002976, -0.315971, -0.555133, -0.095850,
         -0.207090, -0.086065, +0.005617, -0.256976, -0.393331, +0.003183,
         -0.220261, +0.018260, -0.024482, -0.247338, -0.246172, -0.070223,
         -0.136143, -0.383353, -0.008253, -0.127478, -0.638260, -0.265770,
         -0.217825, -0.039818, +0.083267, -0.096157, -0.273991, -0.095080,
         +0.017088, +0.100889, +0.122720, -0.111532, -0.125117, +0.091137,
         -0.062550, +0.096134, +0.286642, +0.065641, -0.082587, +0.172500,
         +0.064989, -0.141004, +0.223305, -0.172570, -0.362464, +0.122049,
         -0.442175, -0.696125, +0.145688, -0.397721, -0.669501, -0.096621,
         -0.415235, -0.551976, -0.020117, -0.289330, -0.434870, -0.099156,
         -0.135692, -0.839601, -0.370336, +0.044901, -0.488363, -0.194583,
         -0.186298, -0.299832, +0.008876, +0.273937, -0.235734, +0.065720,
         -0.447454, -0.508404, -0.049338, -0.700900, -0.081796, +0.229061,
         -0.597343, -0.063885, +0.150099, -0.382640, -0.187373, +0.006893,
         -0.351521, -0.173152, +0.065005, -0.538410, -0.236101, +0.172070,
         -0.408010, -0.053649, +0.099167, -0.389010, -0.237217, -0.062455,
         -0.305336, -0.042800, -0.015233, -0.181560, -0.265425, -0.181198,
         -0.159504, -0.290039, +0.007613, -0.405291, -0.308462, -0.148630,
         -0.307276, -0.139741, +0.025647, -0.319713, -0.217274, +0.032810,
         -0.141414, +0.016770, +0.075093, -0.410465, -0.234103, -0.012248,
         -0.215067, -0.014769, +0.152827, -0.443142, -0.368229, +0.340488,
         -0.223185, -0.159651, +0.280218;

    Eigen::MatrixXd macroD( N.cols( ), 1 );

    macroD << -0.942534, +0.179256, +0.819716, +0.453604, +0.857718, +0.104167,
              -0.531297, -0.616251, +0.726625, +0.713301, -0.561171, -0.036437,
              +0.226544, -0.764067, -0.567154, -0.083834, -0.760801, -0.184202,
              +0.099935, -0.981089, -0.640083, +0.471241, +0.284384, +0.911188,
              -0.612098, -0.151590, -0.359352, -0.498748, +0.681872, +0.931696,
              -0.130505, +0.258422, +0.598219, +0.449634, +0.437597, +0.189190,
              -0.725022, -0.415684, +0.225260, +0.777793, -0.316170, -0.697904,
              +0.760474, -0.172924, +0.469180, -0.923765, +0.554894, -0.436341,
              -0.584481, -0.417923, -0.484523, -0.042049, +0.580823, -0.183014,
              -0.286460, +0.753883, -0.669810, +0.192213, -0.784086, -0.479125,
              -0.102530, -0.289361, +0.034742, +0.471416, -0.674051, +0.672879,
              -0.177298, +0.925295, -0.369792, -0.364725, -0.197006, -0.405645,
              -0.253373, +0.669836, +0.545734, -0.563213, +0.781067, -0.720527,
              -0.803555, -0.153175, +0.275870, +0.938778, +0.031605, +0.964556,
              +0.484850, +0.129173, +0.201998, -0.189893, +0.740700, -0.353216,
              +0.770499, -0.982987, -0.968853, -0.971307, +0.447054, -0.359427,
              -0.567646, -0.661847, +0.304800, -0.719283, -0.056555, -0.766686,
              +0.119887, +0.525746, +0.649757, -0.457472, -0.811261, +0.059171,
              +0.029514, -0.785691, -0.929496, +0.266773, +0.672461, +0.022964,
              -0.971627, +0.648887, -0.556750, -0.568914, +0.036906, +0.370040,
              +0.488436, -0.216337, +0.139964, +0.762068, +0.872199, +0.852070,
              +0.211857, +0.395989, +0.230612, +0.163435, +0.272140, +0.406392,
              -0.962130, -0.700213, +0.079188, -0.223784, -0.630731, -0.245312,
              +0.465122, +0.385622, -0.372280, +0.479940, -0.277007, -0.881563;

    Eigen::MatrixXd R = N * macroD;

    if ( ( A - R ).norm( ) > ( 1e-6 * A.norm( ) ) + 1e-6 ){

        results << "test_overlapCoupling_initializeCoupling_l2_projection (" + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    //Check the center of mass interpolation matrix and the center of mass projector
    SparseMatrix centerOfMassInterpolator;
    error = overlapCoupling::readSparseMatrixFromXDMF( _readGrid, "centerOfMassInterpolator", centerOfMassInterpolator );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_initializeCoupling_l2_projection & False\n";
        return 1;
    }

    Eigen::MatrixXd DX( 12, 1 );
    DX << -1.00911786,  1.51428288,  1.75159184, -0.77596151, -0.13860077,
          -1.30538174, -1.11042458, -0.86808735,  0.47158175, -1.21084958,
           1.4369616 , -0.41944997;

    Eigen::MatrixXd PA( 16, 1 );
    PA << -0.40381902, -0.85441036, -0.64758047, -0.59150002, -0.05132026,
          -0.46167635,  0.26993602, -0.09548531, -0.37851177,  0.45119406,
           0.53527639, -0.35298579, -0.4960972 , -0.53141861, -0.51544242,
          -0.65371376;
    
    R = centerOfMassInterpolator * DX;

    if ( ( R - PA ).norm( ) > 1e-6 * ( PA.norm( ) + 1 ) ){

        results << "test_overlapCoupling_initializeCoupling_l2_projection (" + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;


    Eigen::MatrixXd centerOfMassProjector;
    error = overlapCoupling::readDenseMatrixFromXDMF( _readGrid, "centerOfMassProjector", centerOfMassProjector );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_initializeCoupling_l2_projection & False\n";
        return 1;
    }

    R = centerOfMassProjector * PA;

    if ( ( R - DX ).norm( ) > 1e-6 * ( DX.norm( ) + 1 ) ){

        results << "test_overlapCoupling_initializeCoupling_l2_projection (" + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    //Check the projection matrices
    Eigen::MatrixXd BDhatQ;
    overlapCoupling::readDenseMatrixFromXDMF( _readGrid, "BDhatQ", BDhatQ );

    Eigen::MatrixXd Q( 81, 1 );
    Q << 1.6082461 ,  0.23123014,  0.62393882,  1.32988565, -1.20046325,
        -1.49098297,  1.08575643, -0.27084579, -0.45887108,  1.13457348,
         1.14212648, -1.34876558, -1.53954667, -0.6699138 ,  0.4062938 ,
         1.51120934,  0.45950889, -0.3039844 ,  1.8313851 ,  1.41633787,
         1.0965811 ,  1.50251364, -1.68657903, -1.87216511,  0.82496983,
         0.21188063,  1.42106996,  1.81642989, -0.1000955 ,  0.19266961,
         0.93810141,  0.15452743,  0.98045664,  0.3140218 , -1.29539698,
         1.0298772 ,  1.79294532,  1.51096488,  1.42206134, -0.7942898 ,
        -1.56131436,  1.62426425,  1.67991981, -0.33085656, -1.8824174 ,
        -1.98883142, -1.86904329, -1.5365518 ,  1.39131847, -0.47607648,
         0.00974553, -0.15420091, -0.6692329 , -0.29326975, -1.78084752,
         1.97746862, -0.418282  , -1.04194253,  0.15101235,  1.55810889,
         0.29150197, -0.99929398, -0.4581576 ,  1.09085781, -0.59822029,
        -0.22436283, -0.34358714,  0.15518958,  1.67276323, -0.94694814,
         1.11237832,  0.39840522, -1.04803035,  0.15294796, -0.5688733 ,
        -0.3469194 ,  0.02140078, -1.85645887, -0.78465718,  1.49107402,
         1.9616645;

    Eigen::MatrixXd D( 48, 1 );
    D << -0.24194266,  1.25961845, -0.87935036, -1.71921134,  1.70558356,
          0.75569485, -1.69431444,  0.7158976 ,  0.8212172 , -1.45008094,
          1.56941873,  1.78945147, -1.65800529,  0.34847407, -0.42676962,
         -0.19490982, -0.01828974,  1.7880325 ,  0.32964821, -1.07369484,
          0.46494527, -1.86369121, -1.56866323,  0.00889209,  0.16946288,
         -1.94731671, -1.81322178,  1.28646336,  0.85564197,  0.28811254,
         -0.46973343,  0.14448512, -1.03384903,  0.15534826, -0.77913744,
          1.22798127,  0.06452942,  0.09612534,  1.43803989, -0.57649306,
         -1.68445039, -0.46275924,  1.60444853,  1.23426519, -1.0681013 ,
          0.60927561, -0.21281336, -1.07731193;

    Eigen::MatrixXd DhatAnswer1( 96, 1 );
    DhatAnswer1 <<   0.05449762,   0.37789505,  -0.56974931,  -0.41063288,
                   -21.40331577,   3.37113424,  -2.11587749,   6.80072956,
                    -2.86345496, -22.90626947,   0.77140762,  -2.83040308,
                    -0.09033652,  -0.28276276,   3.13348275,   0.33821581,
                    12.32822362,   6.53850242,   1.78554859,  -2.10831916,
                     3.58290492,  24.75788551,   7.0184034 ,  -5.20611124,
                     6.08524375,   2.98060686,  -1.57083513,  -5.37640008,
                    -9.24043348,  -9.96016475, -13.94944368,   3.74000397,
                    -8.85502491, -28.6549696 ,  -9.37056234,   5.37402281,
                    -0.50808992,   0.65670541,   3.15621434,   8.67306691,
                    21.122022  ,   1.76581633,  15.1113944 ,  -6.66132438,
                     1.141037  ,  26.29144487,   1.0915742 ,  -2.41849736,
                     0.47916177,  -0.36041937,  -4.64044921,   2.36811439,
                    22.30927678, -11.52665762,   1.89119019, -15.37170053,
                     9.76587169,   9.66743994,  -8.01008613,   8.45287362,
                     2.54158223,   1.30638226,   3.12256768,  -1.33690416,
                    -9.43707769, -21.85436559,  -1.05778937,   6.30992196,
                   -11.98007612,  -5.78593149,   2.0443998 ,  18.6658619 ,
                    -5.98023142,  -4.05929903,  -3.3194877 ,  -0.31886817,
                     5.17617087,  32.11206995,  10.0659197 ,  -8.99276261,
                    29.14955761,   8.00567045,  -5.26542749, -19.44939161,
                     1.5113475 ,  -0.50557262,   2.1030087 ,  -3.42692129,
                   -21.79318391,  -5.60993344, -11.84278291,  15.29912391,
                    -4.16057157, -10.71691865,  11.38181508,   8.71119504;

    Eigen::MatrixXd DhatResult = BDhatQ * Q;

    if ( ( DhatAnswer1 - DhatResult ).norm( ) > 1e-6 * ( DhatAnswer1.norm( ) + 1 ) ){

        results << "test_overlapCoupling_initializeCoupling_l2_projection (test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    Eigen::MatrixXd BDhatD;
    overlapCoupling::readDenseMatrixFromXDMF( _readGrid, "BDhatD", BDhatD );

    Eigen::MatrixXd DhatAnswer2( 96, 1 );
    DhatAnswer2 <<  -3.32011991e-02, -2.11061647e-01, -1.17399639e-02,  1.94653084e-01,
                    -2.41843825e-01,  3.75291840e-02,  1.94530314e-01,  9.40304243e-02,
                     1.40679257e-01,  4.50060810e-03, -1.75849596e-01,  1.75234795e-03,
                     2.83948570e-01, -4.56688945e-02,  1.31477626e-01, -3.60781989e-02,
                    -2.20871708e-02, -1.89508586e-01, -1.11833938e-01,  1.86761177e-01,
                     4.22868547e-02,  6.71081869e-02,  9.98195528e-02, -6.93222860e-02,
                    -6.89716326e-02,  3.32587595e-01,  2.91737120e-01, -1.64667501e-01,
                    -1.54372931e-01,  6.81047326e-02,  2.38097094e-01,  2.36706741e-03,
                    -2.22646593e-01,  1.51187099e-01, -1.96898056e-02, -1.94319067e-01,
                     2.57124178e-02,  8.13097088e-03, -1.75590644e-01,  1.17325476e-01,
                     2.71300634e-01, -2.53644156e-02, -7.58687823e-02,  1.55658846e-02,
                    -4.63795374e-04,  8.24767835e-02,  9.39242564e-02,  1.28992144e-01,
                     1.48786626e-02,  9.67932754e-02,  4.84984114e-03, -7.93761706e-02,
                     8.98382034e-02, -1.22352792e-01, -7.64555286e-02, -3.99504712e-02,
                    -4.38378182e-01, -4.38274389e-03,  6.94649168e-02, -3.53031624e-03,
                    -1.28944873e-01,  2.15844277e-02, -5.98084940e-02,  7.46440297e-03,
                     1.53788719e-02,  5.91006022e-01,  3.88511047e-02, -7.80383390e-02,
                    -1.35658625e-01, -2.79464237e-02, -4.07822162e-02,  2.09860973e-01,
                     3.21476740e-02, -1.51767721e-01, -1.32977635e-01,  6.41124026e-02,
                     6.51674014e-02, -2.19689372e-01, -9.60395007e-02, -8.63773552e-03,
                     6.94361272e-01, -6.60845391e-02,  4.19764574e-03,  6.05726668e-01,
                    -1.20339547e-02, -3.94504797e-03,  7.96393755e-02, -4.20215882e-02,
                    -1.03294512e-01,  8.17987982e-02,  2.21281639e-02, -1.04186905e-02,
                    -3.25827746e-04, -4.02239661e-02, -3.20701496e-02, -4.06341701e-01;

    DhatResult = BDhatD * D;

    if ( ( DhatAnswer2 - DhatResult ).norm( ) > 1e-6 * ( DhatAnswer2.norm( ) + 1 ) ){

        results << "test_overlapCoupling_initializeCoupling_l2_projection (test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    Eigen::MatrixXd BQhatQ;
    overlapCoupling::readDenseMatrixFromXDMF( _readGrid, "BQhatQ", BQhatQ );

    Eigen::MatrixXd QhatAnswer1( 54, 1 );
    QhatAnswer1 <<  0.22416624, -0.02292138,  0.59147205,  0.62525827, -0.12574234,
                    1.01507914,  0.16866512,  0.14663159, -0.00666511,  0.505569  ,
                    0.24985976, -0.17859599, -0.51723019,  0.18822046,  0.4243965 ,
                   -0.94079483,  0.29488332,  0.7270837 ,  0.25971782,  0.59860762,
                    0.52114014,  0.54623701,  1.08792558,  0.96288982, -0.07125665,
                    0.04498831,  0.77091726,  0.01703845,  0.08306936,  1.37916797,
                    0.6210215 ,  0.22276315,  0.39031317,  1.20198021,  0.30850459,
                    0.74355473,  1.00566679,  0.81265456,  0.54679542,  1.89669757,
                    1.45013323,  1.10007569,  0.96672725,  0.46433354, -0.75352109,
                    1.69377718,  0.66153093, -1.41851883,  0.39773832, -0.14029142,
                    0.01625116,  0.80588194, -0.36379848, -0.0430513 ;

    Eigen::MatrixXd QhatResult = BQhatQ * Q;

    if ( ( QhatAnswer1 - QhatResult ).norm( ) > 1e-6 * ( QhatAnswer1.norm( ) + 1 ) ){

        results << "test_overlapCoupling_initializeCoupling_l2_projection (test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    Eigen::MatrixXd BQhatD;
    overlapCoupling::readDenseMatrixFromXDMF( _readGrid, "BQhatD", BQhatD );

    Eigen::MatrixXd QhatAnswer2( 54, 1 );
    QhatAnswer2 << -0.37582848,  0.44064848, -0.37554179, -0.11016934,  0.29553775,
                   -0.08169968, -0.8273574 ,  0.16322941, -0.47229178, -0.33408832,
                    0.12543549, -0.16411263, -0.01849311,  0.22514316,  0.00885351,
                    0.05018321,  0.10252856,  0.07026189, -0.42812998, -0.03899359,
                   -0.3483852 , -0.1551361 , -0.05436877, -0.12348079, -0.992427  ,
                    0.06374769, -0.64447523, -0.38733   ,  0.06135132, -0.29147012,
                   -0.53623405, -0.30760085, -0.78285188, -0.18079705, -0.21392971,
                   -0.3632372 , -0.25288531, -0.0866296 , -0.31284907, -0.10963659,
                   -0.14206748, -0.13550156, -0.1142353 , -0.49321385, -0.89513038,
                    0.0053624 , -0.37854635, -0.41945456, -0.07492904,  0.13307628,
                    0.19040976, -0.03599619, -0.01693391,  0.10136851;

    QhatResult = BQhatD * D;

    if ( ( QhatAnswer2 - QhatResult ).norm( ) > 1e-6 * ( QhatAnswer2.norm( ) + 1 ) ){

        results << "test_overlapCoupling_initializeCoupling_l2_projection (test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_initializeCoupling_l2_projection & True\n";
    return 0;
}

int test_overlapCoupling_initializeCoupling_averaged_l2_projection( std::ofstream &results ){
    /*!
     * Test the initialization of the coupling for the averaged_l2_projection
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string filename = "testConfig_averaged_l2_projection.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection & False\n";
        return 1;
    }

    if ( !std::ifstream( "reference_information.xdmf" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection (test 1) & False\n";
        return 1;

    }

    if ( !std::ifstream( "reference_information.h5" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection (test 2) & False\n";
        return 1;

    }

    if ( !std::ifstream( "homogenized_response.xdmf" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection (test 3) & False\n";
        return 1;

    }

    if ( !std::ifstream( "homogenized_response.h5" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection (test 4) & False\n";
        return 1;

    }

    std::string testName = "overlapCoupling_initializeCoupling";
    uIntType testNum = 4;

    cellDomainFloatMap domainMassAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1" , 0.25 },
                    { "ghost_nodeset_volume_2" , 0.25 },
                    { "ghost_nodeset_volume_3" , 0.25 },
                    { "ghost_nodeset_volume_4" , 0.25 },
                    { "ghost_nodeset_volume_5" , 0.25 },
                    { "ghost_nodeset_volume_6" , 0.25 },
                    { "ghost_nodeset_volume_7" , 0.25 },
                    { "ghost_nodeset_volume_8" , 0.25 }
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1" , 0.25 },
                    { "free_nodeset_volume_2" , 0.25 },
                    { "free_nodeset_volume_3" , 0.25 },
                    { "free_nodeset_volume_4" , 0.25 },
                    { "free_nodeset_volume_5" , 0.375 },
                    { "free_nodeset_volume_6" , 0.375 },
                    { "free_nodeset_volume_7" , 0.375 },
                    { "free_nodeset_volume_8" , 0.375 }
                }
            }
        };

    if ( _compare_cellDomainMaps( results, domainMassAnswer, oc._test_domainMass,
                                  testName, testNum ) ){
        return 1;
    }

    cellDomainFloatVectorMap domainCOMAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { 0.250000, 0.250000, 0.250000 } },
                    { "ghost_nodeset_volume_2", { 0.750000, 0.250000, 0.250000 } },
                    { "ghost_nodeset_volume_3", { 0.750000, 0.750000, 0.250000 } },
                    { "ghost_nodeset_volume_4", { 0.250000, 0.750000, 0.250000 } },
                    { "ghost_nodeset_volume_5", { 0.250000, 0.250000, 0.750000 } },
                    { "ghost_nodeset_volume_6", { 0.750000, 0.250000, 0.750000 } },
                    { "ghost_nodeset_volume_7", { 0.750000, 0.750000, 0.750000 } },
                    { "ghost_nodeset_volume_8", { 0.250000, 0.750000, 0.750000 } }
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { 0.250000, 0.250000, 1.250000 } },
                    { "free_nodeset_volume_2", { 0.750000, 0.250000, 1.250000 } },
                    { "free_nodeset_volume_3", { 0.750000, 0.750000, 1.250000 } },
                    { "free_nodeset_volume_4", { 0.250000, 0.750000, 1.250000 } },
                    { "free_nodeset_volume_5", { 0.250000, 0.250000, 1.833333 } },
                    { "free_nodeset_volume_6", { 0.750000, 0.250000, 1.833333 } },
                    { "free_nodeset_volume_7", { 0.750000, 0.750000, 1.833333 } },
                    { "free_nodeset_volume_8", { 0.250000, 0.750000, 1.833333 } }
                }
            }
        };

    if ( _compare_cellDomainMaps( results, domainCOMAnswer, oc._test_domainCOM,
                                  testName, testNum ) ){
        return 1;
    }

    cellDomainFloatVectorMap freeDomainMomentsOfInertiaAnswer
        =
        {
            { 2,
                {
                    { "free_nodeset_volume_1", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "free_nodeset_volume_2", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "free_nodeset_volume_3", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "free_nodeset_volume_4", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "free_nodeset_volume_5", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.055556 } },
                    { "free_nodeset_volume_6", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.055556 } },
                    { "free_nodeset_volume_7", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.055556 } },
                    { "free_nodeset_volume_8", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.055556 } }
                }
            }
        };

    if ( _compare_cellDomainMaps( results, freeDomainMomentsOfInertiaAnswer, *oc.getReferenceFreeMicroDomainMomentsOfInertia( ),
                                  testName, testNum ) ){
        return 1;
    }

    cellDomainFloatVectorMap ghostDomainMomentsOfInertiaAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_2", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_3", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_4", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_5", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_6", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_7", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "ghost_nodeset_volume_8", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } }
                }
            },
        };

    if ( _compare_cellDomainMaps( results, ghostDomainMomentsOfInertiaAnswer, *oc.getReferenceGhostMicroDomainMomentsOfInertia( ),
                                  testName, testNum ) ){
        return 1;
    }

    std::unordered_map< uIntType, std::unordered_map< std::string, std::unordered_map< uIntType, floatVector > > > domainXiAnswer
        =
        {
            {
                1,
                {
                    {
                        "ghost_nodeset_volume_1",
                        {
                            { 24, { -0.250000, -0.250000, -0.250000,  } },
                            { 39, { -0.250000, -0.250000, +0.250000,  } },
                            { 40, { +0.250000, -0.250000, -0.250000,  } },
                            { 57, { +0.250000, -0.250000, +0.250000,  } },
                            { 44, { -0.250000, +0.250000, -0.250000,  } },
                            { 58, { -0.250000, +0.250000, +0.250000,  } },
                            { 29, { +0.250000, +0.250000, -0.250000,  } },
                            { 59, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_2",
                        {
                            { 40, { -0.250000, -0.250000, -0.250000,  } },
                            { 57, { -0.250000, -0.250000, +0.250000,  } },
                            { 11, { +0.250000, -0.250000, -0.250000,  } },
                            {  0, { +0.250000, -0.250000, +0.250000,  } },
                            { 29, { -0.250000, +0.250000, -0.250000,  } },
                            { 59, { -0.250000, +0.250000, +0.250000,  } },
                            { 20, { +0.250000, +0.250000, -0.250000,  } },
                            { 60, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_3",
                        {
                            { 29, { -0.250000, -0.250000, -0.250000,  } },
                            { 59, { -0.250000, -0.250000, +0.250000,  } },
                            { 20, { +0.250000, -0.250000, -0.250000,  } },
                            { 60, { +0.250000, -0.250000, +0.250000,  } },
                            { 47, { -0.250000, +0.250000, -0.250000,  } },
                            { 49, { -0.250000, +0.250000, +0.250000,  } },
                            { 17, { +0.250000, +0.250000, -0.250000,  } },
                            { 38, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_4",
                        {
                            { 44, { -0.250000, -0.250000, -0.250000,  } },
                            { 58, { -0.250000, -0.250000, +0.250000,  } },
                            { 29, { +0.250000, -0.250000, -0.250000,  } },
                            { 59, { +0.250000, -0.250000, +0.250000,  } },
                            { 14, { -0.250000, +0.250000, -0.250000,  } },
                            { 55, { -0.250000, +0.250000, +0.250000,  } },
                            { 47, { +0.250000, +0.250000, -0.250000,  } },
                            { 49, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_5",
                        {
                            { 39, { -0.250000, -0.250000, -0.250000,  } },
                            { 15, { -0.250000, -0.250000, +0.250000,  } },
                            { 57, { +0.250000, -0.250000, -0.250000,  } },
                            { 13, { +0.250000, -0.250000, +0.250000,  } },
                            { 58, { -0.250000, +0.250000, -0.250000,  } },
                            { 53, { -0.250000, +0.250000, +0.250000,  } },
                            { 59, { +0.250000, +0.250000, -0.250000,  } },
                            { 37, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_6",
                        {
                            { 57, { -0.250000, -0.250000, -0.250000,  } },
                            { 13, { -0.250000, -0.250000, +0.250000,  } },
                            {  0, { +0.250000, -0.250000, -0.250000,  } },
                            {  5, { +0.250000, -0.250000, +0.250000,  } },
                            { 59, { -0.250000, +0.250000, -0.250000,  } },
                            { 37, { -0.250000, +0.250000, +0.250000,  } },
                            { 60, { +0.250000, +0.250000, -0.250000,  } },
                            {  3, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_7",
                        {
                            { 59, { -0.250000, -0.250000, -0.250000,  } },
                            { 37, { -0.250000, -0.250000, +0.250000,  } },
                            { 60, { +0.250000, -0.250000, -0.250000,  } },
                            {  3, { +0.250000, -0.250000, +0.250000,  } },
                            { 49, { -0.250000, +0.250000, -0.250000,  } },
                            { 32, { -0.250000, +0.250000, +0.250000,  } },
                            { 38, { +0.250000, +0.250000, -0.250000,  } },
                            { 34, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "ghost_nodeset_volume_8",
                        {
                            { 58, { -0.250000, -0.250000, -0.250000,  } },
                            { 53, { -0.250000, -0.250000, +0.250000,  } },
                            { 59, { +0.250000, -0.250000, -0.250000,  } },
                            { 37, { +0.250000, -0.250000, +0.250000,  } },
                            { 55, { -0.250000, +0.250000, -0.250000,  } },
                            { 25, { -0.250000, +0.250000, +0.250000,  } },
                            { 49, { +0.250000, +0.250000, -0.250000,  } },
                            { 32, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                }
            },
            {
                2,
                {
                    {
                        "free_nodeset_volume_1",
                        {
                            { 15, { -0.250000, -0.250000, -0.250000,  } },
                            { 31, { -0.250000, -0.250000, +0.250000,  } },
                            { 13, { +0.250000, -0.250000, -0.250000,  } },
                            { 26, { +0.250000, -0.250000, +0.250000,  } },
                            { 53, { -0.250000, +0.250000, -0.250000,  } },
                            { 21, { -0.250000, +0.250000, +0.250000,  } },
                            { 37, { +0.250000, +0.250000, -0.250000,  } },
                            { 48, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_2",
                        {
                            { 13, { -0.250000, -0.250000, -0.250000,  } },
                            { 26, { -0.250000, -0.250000, +0.250000,  } },
                            {  5, { +0.250000, -0.250000, -0.250000,  } },
                            { 10, { +0.250000, -0.250000, +0.250000,  } },
                            { 37, { -0.250000, +0.250000, -0.250000,  } },
                            { 48, { -0.250000, +0.250000, +0.250000,  } },
                            {  3, { +0.250000, +0.250000, -0.250000,  } },
                            {  4, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_3",
                        {
                            { 37, { -0.250000, -0.250000, -0.250000,  } },
                            { 48, { -0.250000, -0.250000, +0.250000,  } },
                            {  3, { +0.250000, -0.250000, -0.250000,  } },
                            {  4, { +0.250000, -0.250000, +0.250000,  } },
                            { 32, { -0.250000, +0.250000, -0.250000,  } },
                            { 33, { -0.250000, +0.250000, +0.250000,  } },
                            { 34, { +0.250000, +0.250000, -0.250000,  } },
                            { 28, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_4",
                        {
                            { 53, { -0.250000, -0.250000, -0.250000,  } },
                            { 21, { -0.250000, -0.250000, +0.250000,  } },
                            { 37, { +0.250000, -0.250000, -0.250000,  } },
                            { 48, { +0.250000, -0.250000, +0.250000,  } },
                            { 25, { -0.250000, +0.250000, -0.250000,  } },
                            { 50, { -0.250000, +0.250000, +0.250000,  } },
                            { 32, { +0.250000, +0.250000, -0.250000,  } },
                            { 33, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_5",
                        {
                            { 31, { -0.250000, -0.250000, -0.333333,  } },
                            { 43, { -0.250000, -0.250000, +0.166667,  } },
                            { 26, { +0.250000, -0.250000, -0.333333,  } },
                            { 27, { +0.250000, -0.250000, +0.166667,  } },
                            { 21, { -0.250000, +0.250000, -0.333333,  } },
                            {  1, { -0.250000, +0.250000, +0.166667,  } },
                            { 48, { +0.250000, +0.250000, -0.333333,  } },
                            {  7, { +0.250000, +0.250000, +0.166667,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_6",
                        {
                            { 26, { -0.250000, -0.250000, -0.333333,  } },
                            { 27, { -0.250000, -0.250000, +0.166667,  } },
                            { 10, { +0.250000, -0.250000, -0.333333,  } },
                            { 30, { +0.250000, -0.250000, +0.166667,  } },
                            { 48, { -0.250000, +0.250000, -0.333333,  } },
                            {  7, { -0.250000, +0.250000, +0.166667,  } },
                            {  4, { +0.250000, +0.250000, -0.333333,  } },
                            { 16, { +0.250000, +0.250000, +0.166667,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_7",
                        {
                            { 48, { -0.250000, -0.250000, -0.333333,  } },
                            {  7, { -0.250000, -0.250000, +0.166667,  } },
                            {  4, { +0.250000, -0.250000, -0.333333,  } },
                            { 16, { +0.250000, -0.250000, +0.166667,  } },
                            { 33, { -0.250000, +0.250000, -0.333333,  } },
                            { 22, { -0.250000, +0.250000, +0.166667,  } },
                            { 28, { +0.250000, +0.250000, -0.333333,  } },
                            {  2, { +0.250000, +0.250000, +0.166667,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_8",
                        {
                            { 21, { -0.250000, -0.250000, -0.333333,  } },
                            {  1, { -0.250000, -0.250000, +0.166667,  } },
                            { 48, { +0.250000, -0.250000, -0.333333,  } },
                            {  7, { +0.250000, -0.250000, +0.166667,  } },
                            { 50, { -0.250000, +0.250000, -0.333333,  } },
                            { 46, { -0.250000, +0.250000, +0.166667,  } },
                            { 33, { +0.250000, +0.250000, -0.333333,  } },
                            { 22, { +0.250000, +0.250000, +0.166667,  } },
                        }
                    },
                }
            }
        };

    if ( _compare_cellDomainPointMaps( results, domainXiAnswer, oc._test_domainXi,
                                       testName, testNum ) ){
        return 1;
    }

//    std::cout << "shape functions\n";
    cellDomainFloatVectorMap domainCOMSFAnswer
        =
        {
            {
                1,
                {
                    {
                        "ghost_nodeset_volume_1",
                        { +0.421875, +0.140625, +0.046875, +0.140625, +0.140625, +0.046875, +0.015625, +0.046875 }
                    },
                    {
                        "ghost_nodeset_volume_2",
                        { +0.140625, +0.421875, +0.140625, +0.046875, +0.046875, +0.140625, +0.046875, +0.015625 }
                    },
                    {
                        "ghost_nodeset_volume_3",
                        { +0.046875, +0.140625, +0.421875, +0.140625, +0.015625, +0.046875, +0.140625, +0.046875 }
                    },
                    {
                        "ghost_nodeset_volume_4",
                        { +0.140625, +0.046875, +0.140625, +0.421875, +0.046875, +0.015625, +0.046875, +0.140625 }
                    },
                    {
                        "ghost_nodeset_volume_5",
                        { +0.140625, +0.046875, +0.015625, +0.046875, +0.421875, +0.140625, +0.046875, +0.140625 }
                    },
                    {
                        "ghost_nodeset_volume_6",
                        { +0.046875, +0.140625, +0.046875, +0.015625, +0.140625, +0.421875, +0.140625, +0.046875 }
                    },
                    {
                        "ghost_nodeset_volume_7",
                        { +0.015625, +0.046875, +0.140625, +0.046875, +0.046875, +0.140625, +0.421875, +0.140625 }
                    },
                    {
                        "ghost_nodeset_volume_8",
                        { +0.046875, +0.015625, +0.046875, +0.140625, +0.140625, +0.046875, +0.140625, +0.421875 }
                    },
                }
            },
            {
                2,
                {
                    {
                        "free_nodeset_volume_1",
                        { +0.421875, +0.140625, +0.046875, +0.140625, +0.140625, +0.046875, +0.015625, +0.046875 }
                    },
                    {
                        "free_nodeset_volume_2",
                        { +0.140625, +0.421875, +0.140625, +0.046875, +0.046875, +0.140625, +0.046875, +0.015625 }
                    },
                    {
                        "free_nodeset_volume_3",
                        { +0.046875, +0.140625, +0.421875, +0.140625, +0.015625, +0.046875, +0.140625, +0.046875 }
                    },
                    {
                        "free_nodeset_volume_4",
                        { +0.140625, +0.046875, +0.140625, +0.421875, +0.046875, +0.015625, +0.046875, +0.140625 }
                    },
                    {
                        "free_nodeset_volume_5",
                        { +0.093750, +0.031250, +0.010417, +0.031250, +0.468750, +0.156250, +0.052083, +0.156250 }
                    },
                    {
                        "free_nodeset_volume_6",
                        { +0.031250, +0.093750, +0.031250, +0.010417, +0.156250, +0.468750, +0.156250, +0.052083 }
                    },
                    {
                        "free_nodeset_volume_7",
                        { +0.010417, +0.031250, +0.093750, +0.031250, +0.052083, +0.156250, +0.468750, +0.156250 }
                    },
                    {
                        "free_nodeset_volume_8",
                        { +0.031250, +0.010417, +0.031250, +0.093750, +0.156250, +0.052083, +0.156250, +0.468750 }
                    },
                }
            }
        };

    if ( _compare_cellDomainMaps( results, domainCOMSFAnswer, *oc.getReferenceCellDomainCenterOfMassShapeFunctions( ),
                                  testName, testNum ) ){
        return 1;
    }

//    std::unordered_map< uIntType, std::unordered_map< std::string, std::unordered_map< uIntType, floatVector > > > domainPointSFAnswer
//        =
//        {
//            {
//                1,
//                {
//                    {
//                        "ghost_nodeset_volume_1",
//                        {
//                            { 24, { +1.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 39, { +0.500000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000,  } },
//                            { 40, { +0.500000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 57, { +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000,  } },
//                            { 44, { +0.500000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 58, { +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000,  } },
//                            { 29, { +0.250000, +0.250000, +0.250000, +0.250000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 59, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                        }
//                    },
//                    {
//                        "ghost_nodeset_volume_2",
//                        {
//                            { 40, { +0.500000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 57, { +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000,  } },
//                            { 11, { +0.000000, +1.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            {  0, { +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000,  } },
//                            { 29, { +0.250000, +0.250000, +0.250000, +0.250000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 59, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            { 20, { +0.000000, +0.500000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 60, { +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000,  } },
//                        }
//                    },
//                    {
//                        "ghost_nodeset_volume_3",
//                        {
//                            { 29, { +0.250000, +0.250000, +0.250000, +0.250000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 59, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            { 20, { +0.000000, +0.500000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 60, { +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000,  } },
//                            { 47, { +0.000000, +0.000000, +0.500000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 49, { +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000,  } },
//                            { 17, { +0.000000, +0.000000, +1.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 38, { +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000,  } },
//                        }
//                    },
//                    {
//                        "ghost_nodeset_volume_4",
//                        {
//                            { 44, { +0.500000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 58, { +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000,  } },
//                            { 29, { +0.250000, +0.250000, +0.250000, +0.250000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 59, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            { 14, { +0.000000, +0.000000, +0.000000, +1.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 55, { +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.500000,  } },
//                            { 47, { +0.000000, +0.000000, +0.500000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 49, { +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000,  } },
//                        }
//                    },
//                    {
//                        "ghost_nodeset_volume_5",
//                        {
//                            { 39, { +0.500000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000,  } },
//                            { 15, { +0.000000, +0.000000, +0.000000, +0.000000, +1.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 57, { +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000,  } },
//                            { 13, { +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.500000, +0.000000, +0.000000,  } },
//                            { 58, { +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000,  } },
//                            { 53, { +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.500000,  } },
//                            { 59, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            { 37, { +0.000000, +0.000000, +0.000000, +0.000000, +0.250000, +0.250000, +0.250000, +0.250000,  } },
//                        }
//                    },
//                    {
//                        "ghost_nodeset_volume_6",
//                        {
//                            { 57, { +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000,  } },
//                            { 13, { +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.500000, +0.000000, +0.000000,  } },
//                            {  0, { +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000,  } },
//                            {  5, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +1.000000, +0.000000, +0.000000,  } },
//                            { 59, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            { 37, { +0.000000, +0.000000, +0.000000, +0.000000, +0.250000, +0.250000, +0.250000, +0.250000,  } },
//                            { 60, { +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000,  } },
//                            {  3, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.500000, +0.000000,  } },
//                        }
//                    },
//                    {
//                        "ghost_nodeset_volume_7",
//                        {
//                            { 59, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            { 37, { +0.000000, +0.000000, +0.000000, +0.000000, +0.250000, +0.250000, +0.250000, +0.250000,  } },
//                            { 60, { +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000,  } },
//                            {  3, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.500000, +0.000000,  } },
//                            { 49, { +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000,  } },
//                            { 32, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.500000,  } },
//                            { 38, { +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000,  } },
//                            { 34, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +1.000000, +0.000000,  } },
//                        }
//                    },
//                    {
//                        "ghost_nodeset_volume_8",
//                        {
//                            { 58, { +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000,  } },
//                            { 53, { +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.500000,  } },
//                            { 59, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            { 37, { +0.000000, +0.000000, +0.000000, +0.000000, +0.250000, +0.250000, +0.250000, +0.250000,  } },
//                            { 55, { +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.500000,  } },
//                            { 25, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +1.000000,  } },
//                            { 49, { +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000,  } },
//                            { 32, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.500000,  } },
//                        }
//                    },
//                }
//            },
//            {
//                2,
//                {
//                    {
//                        "free_nodeset_volume_1",
//                        {
//                            { 15, { +1.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 31, { +0.500000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000,  } },
//                            { 13, { +0.500000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 26, { +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000,  } },
//                            { 53, { +0.500000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 21, { +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000,  } },
//                            { 37, { +0.250000, +0.250000, +0.250000, +0.250000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 48, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                        }
//                    },
//                    {
//                        "free_nodeset_volume_2",
//                        {
//                            { 13, { +0.500000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 26, { +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000,  } },
//                            {  5, { +0.000000, +1.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 10, { +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000,  } },
//                            { 37, { +0.250000, +0.250000, +0.250000, +0.250000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 48, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            {  3, { +0.000000, +0.500000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            {  4, { +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000,  } },
//                        }
//                    },
//                    {
//                        "free_nodeset_volume_3",
//                        {
//                            { 37, { +0.250000, +0.250000, +0.250000, +0.250000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 48, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            {  3, { +0.000000, +0.500000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            {  4, { +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000,  } },
//                            { 32, { +0.000000, +0.000000, +0.500000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 33, { +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000,  } },
//                            { 34, { +0.000000, +0.000000, +1.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 28, { +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000,  } },
//                        }
//                    },
//                    {
//                        "free_nodeset_volume_4",
//                        {
//                            { 53, { +0.500000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 21, { +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000,  } },
//                            { 37, { +0.250000, +0.250000, +0.250000, +0.250000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 48, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            { 25, { +0.000000, +0.000000, +0.000000, +1.000000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 50, { +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.500000,  } },
//                            { 32, { +0.000000, +0.000000, +0.500000, +0.500000, +0.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 33, { +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000,  } },
//                        }
//                    },
//                    {
//                        "free_nodeset_volume_5",
//                        {
//                            { 31, { +0.500000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000,  } },
//                            { 43, { +0.000000, +0.000000, +0.000000, +0.000000, +1.000000, +0.000000, +0.000000, +0.000000,  } },
//                            { 26, { +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000,  } },
//                            { 27, { +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.500000, +0.000000, +0.000000,  } },
//                            { 21, { +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000,  } },
//                            {  1, { +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.500000,  } },
//                            { 48, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            {  7, { +0.000000, +0.000000, +0.000000, +0.000000, +0.250000, +0.250000, +0.250000, +0.250000,  } },
//                        }
//                    },
//                    {
//                        "free_nodeset_volume_6",
//                        {
//                            { 26, { +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000,  } },
//                            { 27, { +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.500000, +0.000000, +0.000000,  } },
//                            { 10, { +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000,  } },
//                            { 30, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +1.000000, +0.000000, +0.000000,  } },
//                            { 48, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            {  7, { +0.000000, +0.000000, +0.000000, +0.000000, +0.250000, +0.250000, +0.250000, +0.250000,  } },
//                            {  4, { +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000,  } },
//                            { 16, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.500000, +0.000000,  } },
//                        }
//                    },
//                    {
//                        "free_nodeset_volume_7",
//                        {
//                            { 48, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            {  7, { +0.000000, +0.000000, +0.000000, +0.000000, +0.250000, +0.250000, +0.250000, +0.250000,  } },
//                            {  4, { +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000,  } },
//                            { 16, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.500000, +0.000000,  } },
//                            { 33, { +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000,  } },
//                            { 22, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.500000,  } },
//                            { 28, { +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000,  } },
//                            {  2, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +1.000000, +0.000000,  } },
//                        }
//                    },
//                    {
//                        "free_nodeset_volume_8",
//                        {
//                            { 21, { +0.250000, +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000,  } },
//                            {  1, { +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.500000,  } },
//                            { 48, { +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000, +0.125000,  } },
//                            {  7, { +0.000000, +0.000000, +0.000000, +0.000000, +0.250000, +0.250000, +0.250000, +0.250000,  } },
//                            { 50, { +0.000000, +0.000000, +0.000000, +0.500000, +0.000000, +0.000000, +0.000000, +0.500000,  } },
//                            { 46, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +1.000000,  } },
//                            { 33, { +0.000000, +0.000000, +0.250000, +0.250000, +0.000000, +0.000000, +0.250000, +0.250000,  } },
//                            { 22, { +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.000000, +0.500000, +0.500000,  } },
//                        }
//                    },
//                }
//            }
//        };
//
//    if ( _compare_cellDomainPointMaps( results, domainPointSFAnswer, oc._test_domainMUP,
//                                       testName, testNum ) ){
//        return 1;
//    }

    std::string xdmf_filename = "reference_information.xdmf";
    shared_ptr< XdmfReader > reader = XdmfReader::New( );
    shared_ptr< XdmfDomain > _readDomain = shared_dynamic_cast< XdmfDomain >( reader->read( xdmf_filename ) );
    shared_ptr< XdmfUnstructuredGrid > _readGrid = _readDomain->getUnstructuredGrid( 0 );

    SparseMatrix N;
    overlapCoupling::readSparseMatrixFromXDMF( _readGrid, "N", N );

    Eigen::MatrixXd A( N.rows( ), 1 );
    A << -0.416617, -0.311112, +0.013488, -0.337511, -0.627147, +0.058891,
         -0.307855, -0.308153, -0.002976, -0.315971, -0.555133, -0.095850,
         -0.207090, -0.086065, +0.005617, -0.256976, -0.393331, +0.003183,
         -0.220261, +0.018260, -0.024482, -0.247338, -0.246172, -0.070223,
         -0.136143, -0.383353, -0.008253, -0.127478, -0.638260, -0.265770,
         -0.217825, -0.039818, +0.083267, -0.096157, -0.273991, -0.095080,
         +0.017088, +0.100889, +0.122720, -0.111532, -0.125117, +0.091137,
         -0.062550, +0.096134, +0.286642, +0.065641, -0.082587, +0.172500,
         +0.064989, -0.141004, +0.223305, -0.172570, -0.362464, +0.122049,
         -0.442175, -0.696125, +0.145688, -0.397721, -0.669501, -0.096621,
         -0.415235, -0.551976, -0.020117, -0.289330, -0.434870, -0.099156,
         -0.135692, -0.839601, -0.370336, +0.044901, -0.488363, -0.194583,
         -0.186298, -0.299832, +0.008876, +0.273937, -0.235734, +0.065720,
         -0.447454, -0.508404, -0.049338, -0.700900, -0.081796, +0.229061,
         -0.597343, -0.063885, +0.150099, -0.382640, -0.187373, +0.006893,
         -0.351521, -0.173152, +0.065005, -0.538410, -0.236101, +0.172070,
         -0.408010, -0.053649, +0.099167, -0.389010, -0.237217, -0.062455,
         -0.305336, -0.042800, -0.015233, -0.181560, -0.265425, -0.181198,
         -0.159504, -0.290039, +0.007613, -0.405291, -0.308462, -0.148630,
         -0.307276, -0.139741, +0.025647, -0.319713, -0.217274, +0.032810,
         -0.141414, +0.016770, +0.075093, -0.410465, -0.234103, -0.012248,
         -0.215067, -0.014769, +0.152827, -0.443142, -0.368229, +0.340488,
         -0.223185, -0.159651, +0.280218;

    Eigen::MatrixXd macroD( N.cols( ), 1 );

    macroD << -0.942534, +0.179256, +0.819716, +0.453604, +0.857718, +0.104167,
              -0.531297, -0.616251, +0.726625, +0.713301, -0.561171, -0.036437,
              +0.226544, -0.764067, -0.567154, -0.083834, -0.760801, -0.184202,
              +0.099935, -0.981089, -0.640083, +0.471241, +0.284384, +0.911188,
              -0.612098, -0.151590, -0.359352, -0.498748, +0.681872, +0.931696,
              -0.130505, +0.258422, +0.598219, +0.449634, +0.437597, +0.189190,
              -0.725022, -0.415684, +0.225260, +0.777793, -0.316170, -0.697904,
              +0.760474, -0.172924, +0.469180, -0.923765, +0.554894, -0.436341,
              -0.584481, -0.417923, -0.484523, -0.042049, +0.580823, -0.183014,
              -0.286460, +0.753883, -0.669810, +0.192213, -0.784086, -0.479125,
              -0.102530, -0.289361, +0.034742, +0.471416, -0.674051, +0.672879,
              -0.177298, +0.925295, -0.369792, -0.364725, -0.197006, -0.405645,
              -0.253373, +0.669836, +0.545734, -0.563213, +0.781067, -0.720527,
              -0.803555, -0.153175, +0.275870, +0.938778, +0.031605, +0.964556,
              +0.484850, +0.129173, +0.201998, -0.189893, +0.740700, -0.353216,
              +0.770499, -0.982987, -0.968853, -0.971307, +0.447054, -0.359427,
              -0.567646, -0.661847, +0.304800, -0.719283, -0.056555, -0.766686,
              +0.119887, +0.525746, +0.649757, -0.457472, -0.811261, +0.059171,
              +0.029514, -0.785691, -0.929496, +0.266773, +0.672461, +0.022964,
              -0.971627, +0.648887, -0.556750, -0.568914, +0.036906, +0.370040,
              +0.488436, -0.216337, +0.139964, +0.762068, +0.872199, +0.852070,
              +0.211857, +0.395989, +0.230612, +0.163435, +0.272140, +0.406392,
              -0.962130, -0.700213, +0.079188, -0.223784, -0.630731, -0.245312,
              +0.465122, +0.385622, -0.372280, +0.479940, -0.277007, -0.881563;

    Eigen::MatrixXd R = N * macroD;

    if ( ( A - R ).norm( ) > ( 1e-6 * A.norm( ) ) + 1e-6 ){

        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection (" + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    //Check the center of mass interpolation matrix and the center of mass projector
    SparseMatrix centerOfMassInterpolator;
    error = overlapCoupling::readSparseMatrixFromXDMF( _readGrid, "centerOfMassInterpolator", centerOfMassInterpolator );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection & False\n";
        return 1;
    }

    Eigen::MatrixXd DX( 12, 1 );
    DX << -1.00911786,  1.51428288,  1.75159184, -0.77596151, -0.13860077,
          -1.30538174, -1.11042458, -0.86808735,  0.47158175, -1.21084958,
           1.4369616 , -0.41944997;

    Eigen::MatrixXd PA( 16, 1 );
    PA << -0.40381902, -0.85441036, -0.64758047, -0.59150002, -0.05132026,
          -0.46167635,  0.26993602, -0.09548531, -0.37851177,  0.45119406,
           0.53527639, -0.35298579, -0.4960972 , -0.53141861, -0.51544242,
          -0.65371376;
    
    R = centerOfMassInterpolator * DX;

    if ( ( R - PA ).norm( ) > 1e-6 * ( PA.norm( ) + 1 ) ){

        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection (" + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;


    Eigen::MatrixXd centerOfMassProjector;
    error = overlapCoupling::readDenseMatrixFromXDMF( _readGrid, "centerOfMassProjector", centerOfMassProjector );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection & False\n";
        return 1;
    }

    R = centerOfMassProjector * PA;

    if ( ( R - DX ).norm( ) > 1e-6 * ( DX.norm( ) + 1 ) ){

        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection (" + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

//    std::cout << "checking the projection constants\n";
//    //Check the direct projection constants
//    std::cout << "projected mass\n";
//    std::unordered_map< uIntType, floatType > macroNodeProjectedMassAnswer
//        =
//        {
//            {  5, 0.250000 },
//            {  9, 0.250000 },
//            {  8, 0.250000 },
//            { 11, 0.250000 },
//            {  3, 0.500000 },
//            {  1, 0.500000 },
//            {  6, 0.500000 },
//            { 15, 0.500000 },
//            { 12, 0.375000 },
//            {  2, 0.375000 },
//            { 13, 0.375000 },
//            { 14, 0.375000 }
//        };
//
//    if ( _compare_domainMaps( results, macroNodeProjectedMassAnswer, *oc.getMacroNodeProjectedMass( ),
//                              testName, testNum ) ){
//        return 1;
//    }
//
//    //Check the projection constant
//    std::cout << "projected constant\n";
//    std::unordered_map< uIntType, floatVector > macroNodeProjectedConstantAnswer
//        =
//        {
//            {  5, { -0.031250, -0.031250, -0.031250 } },
//            {  9, { +0.031250, -0.031250, -0.031250 } },
//            {  8, { +0.031250, +0.031250, -0.031250 } },
//            { 11, { -0.031250, +0.031250, -0.031250 } },
//            {  3, { -0.062500, -0.062500, -0.005208 } },
//            {  1, { +0.062500, -0.062500, -0.005208 } },
//            {  6, { +0.062500, +0.062500, -0.005208 } },
//            { 15, { -0.062500, +0.062500, -0.005208 } },
//            { 12, { -0.046875, -0.046875, +0.036458 } },
//            {  2, { +0.046875, -0.046875, +0.036458 } },
//            { 13, { +0.046875, +0.046875, +0.036458 } },
//            { 14, { -0.046875, +0.046875, +0.036458 } }
//        };
//
//    if ( _compare_domainMaps( results, macroNodeProjectedConstantAnswer, *oc.getMacroNodeMassRelativePositionConstant( ),
//                              testName, testNum ) ){
//        return 1;
//    }
//
//    //Check the mass-weighted moment of inertia
//    std::cout << "projected mass moment of inertia\n";
//    std::unordered_map< uIntType, floatVector > macroNodeProjectedMassMomentOfInertiaAnswer
//        =
//        {
//            {  5, { +0.015625, +0.003906, +0.003906, +0.003906, +0.015625, +0.003906, +0.003906, +0.003906, +0.015625 } },
//            {  9, { +0.015625, -0.003906, -0.003906, -0.003906, +0.015625, +0.003906, -0.003906, +0.003906, +0.015625 } },
//            {  8, { +0.015625, +0.003906, -0.003906, +0.003906, +0.015625, -0.003906, -0.003906, -0.003906, +0.015625 } },
//            { 11, { +0.015625, -0.003906, +0.003906, -0.003906, +0.015625, -0.003906, +0.003906, -0.003906, +0.015625 } },
//            {  3, { +0.031250, +0.007812, +0.000651, +0.007812, +0.031250, +0.000651, +0.000651, +0.000651, +0.034288 } },
//            {  1, { +0.031250, -0.007812, -0.000651, -0.007812, +0.031250, +0.000651, -0.000651, +0.000651, +0.034288 } },
//            {  6, { +0.031250, +0.007812, -0.000651, +0.007812, +0.031250, -0.000651, -0.000651, -0.000651, +0.034288 } },
//            { 15, { +0.031250, -0.007812, +0.000651, -0.007812, +0.031250, -0.000651, +0.000651, -0.000651, +0.034288 } },
//            { 12, { +0.023438, +0.005859, -0.004557, +0.005859, +0.023438, -0.004557, -0.004557, -0.004557, +0.017795 } },
//            {  2, { +0.023438, -0.005859, +0.004557, -0.005859, +0.023438, -0.004557, +0.004557, -0.004557, +0.017795 } },
//            { 13, { +0.023438, +0.005859, +0.004557, +0.005859, +0.023438, +0.004557, +0.004557, +0.004557, +0.017795 } },
//            { 14, { +0.023438, -0.005859, -0.004557, -0.005859, +0.023438, +0.004557, -0.004557, +0.004557, +0.017795 } }
//        };
//
//    if ( _compare_domainMaps( results, macroNodeProjectedMassMomentOfInertiaAnswer, *oc.getMacroNodeProjectedMassMomentOfInertia( ),
//                              testName, testNum ) ){
//        return 1;
//    }

    //Check the projection matrices
    Eigen::MatrixXd BDhatQ;
    overlapCoupling::readDenseMatrixFromXDMF( _readGrid, "BDhatQ", BDhatQ );

    Eigen::MatrixXd Q( 81, 1 );
    Q << 1.6082461 ,  0.23123014,  0.62393882,  1.32988565, -1.20046325,
        -1.49098297,  1.08575643, -0.27084579, -0.45887108,  1.13457348,
         1.14212648, -1.34876558, -1.53954667, -0.6699138 ,  0.4062938 ,
         1.51120934,  0.45950889, -0.3039844 ,  1.8313851 ,  1.41633787,
         1.0965811 ,  1.50251364, -1.68657903, -1.87216511,  0.82496983,
         0.21188063,  1.42106996,  1.81642989, -0.1000955 ,  0.19266961,
         0.93810141,  0.15452743,  0.98045664,  0.3140218 , -1.29539698,
         1.0298772 ,  1.79294532,  1.51096488,  1.42206134, -0.7942898 ,
        -1.56131436,  1.62426425,  1.67991981, -0.33085656, -1.8824174 ,
        -1.98883142, -1.86904329, -1.5365518 ,  1.39131847, -0.47607648,
         0.00974553, -0.15420091, -0.6692329 , -0.29326975, -1.78084752,
         1.97746862, -0.418282  , -1.04194253,  0.15101235,  1.55810889,
         0.29150197, -0.99929398, -0.4581576 ,  1.09085781, -0.59822029,
        -0.22436283, -0.34358714,  0.15518958,  1.67276323, -0.94694814,
         1.11237832,  0.39840522, -1.04803035,  0.15294796, -0.5688733 ,
        -0.3469194 ,  0.02140078, -1.85645887, -0.78465718,  1.49107402,
         1.9616645;

    Eigen::MatrixXd D( 48, 1 );
    D << -0.24194266,  1.25961845, -0.87935036, -1.71921134,  1.70558356,
          0.75569485, -1.69431444,  0.7158976 ,  0.8212172 , -1.45008094,
          1.56941873,  1.78945147, -1.65800529,  0.34847407, -0.42676962,
         -0.19490982, -0.01828974,  1.7880325 ,  0.32964821, -1.07369484,
          0.46494527, -1.86369121, -1.56866323,  0.00889209,  0.16946288,
         -1.94731671, -1.81322178,  1.28646336,  0.85564197,  0.28811254,
         -0.46973343,  0.14448512, -1.03384903,  0.15534826, -0.77913744,
          1.22798127,  0.06452942,  0.09612534,  1.43803989, -0.57649306,
         -1.68445039, -0.46275924,  1.60444853,  1.23426519, -1.0681013 ,
          0.60927561, -0.21281336, -1.07731193;

    Eigen::MatrixXd DhatAnswer1( 96, 1 );
    DhatAnswer1 << 1.00753524,  -0.0959485 ,  -0.66625132,   2.44514371,
                  -6.26144891,   3.6251003 ,   1.76774608,   2.15231818,
                  -0.29617276,  -6.37111855,   1.42668662,  -4.2850262 ,
                   1.44621962,   0.23464651,   0.55786313,  -1.56777495,
                   2.62357923,   3.85130026,  -1.10655606,  -0.5841102 ,
                   0.12926275,   8.81934745,   3.05590507,  -1.35489717,
                   1.20926748,  -0.13574076,   0.16402413,  -3.63973743,
                  -3.09748352,  -2.91303167,  -6.06361495,  -0.15666434,
                  -5.47178654,  -9.22359596,  -3.84358307,   0.51947356,
                   0.53768312,   0.16217523,   0.89769605,   4.98290614,
                   5.32174467,   1.07605055,   5.46778297,  -1.63607072,
                  -1.20007404,   7.75625213,   1.70120811,   0.39469567,
                   0.16140048,   0.10547403,  -1.44956365,   2.275294  ,
                  11.04128873,  -5.85944414,  -1.40348803,  -8.45642463,
                   0.72614375,   2.71306206,  -3.90484329,   4.05830282,
                   0.44395985,   0.19595009,   1.60498907,  -1.71017525,
                  -0.28317691,  -5.91725792,   1.58444014,  -0.30921696,
                   0.04658451,   3.3960434 ,  -2.78475037,   4.94243764,
                  -1.30186081,  -1.01481543,  -1.67203878,  -2.9131604 ,
                  -3.20846441,   1.05931864,   4.7149403 ,  -2.11231406,
                   6.40715802,  -0.03511209,  -3.76930535,  -4.26232784,
                   0.30771234,  -0.43832265,   0.76778491,  -0.3059859 ,
                 -10.748665  ,  -1.36924381,  -5.86792585,   7.36883127,
                   0.8484395 ,  -4.8445353 ,   8.33954043,   0.36725083;

    Eigen::MatrixXd DhatResult = BDhatQ * Q;

    if ( ( DhatAnswer1 - DhatResult ).norm( ) > 1e-6 * ( DhatAnswer1.norm( ) + 1 ) ){

        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection (test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    Eigen::MatrixXd BDhatD;
    overlapCoupling::readDenseMatrixFromXDMF( _readGrid, "BDhatD", BDhatD );

    Eigen::MatrixXd DhatAnswer2( 96, 1 );
    DhatAnswer2 <<  0.01117404, -0.06213769, -0.00323124,  0.14326677, -0.11300732,
                   -0.00151848,  0.10364709,  0.10441377,  0.00844412,  0.04219827,
                   -0.0810997 ,  0.00043911,  0.080921  , -0.00415823,  0.05832703,
                   -0.00377284, -0.08803151, -0.01099665,  0.01231181,  0.14838209,
                    0.00056508,  0.08091828,  0.0346689 , -0.00792627,  0.00166544,
                    0.09655481,  0.08629784, -0.08395274, -0.07047961, -0.00022632,
                    0.13392157,  0.053044  , -0.01312119,  0.13885125,  0.02127272,
                   -0.01172733,  0.00436115,  0.01059559, -0.03227229,  0.07856132,
                    0.09938155, -0.00059265,  0.03799686,  0.04105279, -0.00143987,
                    0.098289  ,  0.02301762,  0.0043856 , -0.00229423,  0.01275796,
                    0.00066343, -0.02941519,  0.02320239, -0.00944099, -0.02128057,
                   -0.02143798,  0.05250038, -0.00866405,  0.01665119,  0.00273009,
                   -0.0166145 ,  0.00085376, -0.01197557,  0.00077463,  0.01807442,
                   -0.06837048, -0.00252783, -0.03046545,  0.00351331, -0.01661395,
                   -0.00711813, -0.04928074, -0.00034194, -0.0198244 , -0.01771846,
                    0.01723697,  0.0144707 , -0.00140714, -0.02749645, -0.01089086,
                   -0.08157954, -0.0285086 , -0.00436766, -0.07291339, -0.00089542,
                   -0.00217546,  0.00662607, -0.01613002, -0.02040478, -0.00368476,
                   -0.00780142, -0.00842886, -0.00895226, -0.02018046, -0.00472592,
                    0.02726698;

    DhatResult = BDhatD * D;

    if ( ( DhatAnswer2 - DhatResult ).norm( ) > 1e-6 * ( DhatAnswer2.norm( ) + 1 ) ){

        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection (test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    Eigen::MatrixXd BQhatQ;
    overlapCoupling::readDenseMatrixFromXDMF( _readGrid, "BQhatQ", BQhatQ );

    Eigen::MatrixXd QhatAnswer1( 54, 1 );
    QhatAnswer1 <<  0.10712051, -0.0912848 ,  0.11184787,  0.38804931, -0.22690467,
                    0.06382711,  0.27234912,  0.1258435 , -0.05194806,  0.70557497,
                    0.19565114, -0.22893988, -0.13646913,  0.01526603,  0.13639149,
                   -0.15542481, -0.04627468,  0.17507661,  0.2273456 ,  0.25597847,
                    0.19916026,  0.54280712,  0.40510115,  0.32448061,  0.11474337,
                    0.0280608 ,  0.23099898,  0.37743189, -0.01161507,  0.37177906,
                    0.196299  ,  0.04651278,  0.20601083,  0.45131641, -0.04387927,
                    0.36204822,  0.36553366,  0.28108397,  0.24323837,  0.74642245,
                    0.40449223,  0.46384067,  0.19066721,  0.02580771, -0.13810142,
                    0.3508261 , -0.15445757, -0.28593084,  0.07351183, -0.04908502,
                    0.0659988 ,  0.20824225, -0.20744848,  0.09645349;

    Eigen::MatrixXd QhatResult = BQhatQ * Q;

    if ( ( QhatAnswer1 - QhatResult ).norm( ) > 1e-6 * ( QhatAnswer1.norm( ) + 1 ) ){

        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection (test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    Eigen::MatrixXd BQhatD;
    overlapCoupling::readDenseMatrixFromXDMF( _readGrid, "BQhatD", BQhatD );

    Eigen::MatrixXd QhatAnswer2( 54, 1 );
    QhatAnswer2 <<  -0.38024752,  0.46469597, -0.37771336, -0.11826487,  0.33930194,
                    -0.08611647, -0.84585307,  0.17462444, -0.48222434, -0.36815128,
                     0.14577774, -0.18211656, -0.01682891,  0.23384865,  0.01568758,
                     0.05403782,  0.11812839,  0.08288495, -0.43687213, -0.04302016,
                    -0.35728098, -0.17112188, -0.0618814 , -0.13942753, -1.0192094 ,
                     0.06916286, -0.66372833, -0.43578057,  0.07161681, -0.32618026,
                    -0.54903347, -0.32445362, -0.80910362, -0.20392506, -0.24474305,
                    -0.41100596, -0.25659012, -0.10282278, -0.32256879, -0.11697756,
                    -0.17092501, -0.15311257, -0.11072741, -0.52792405, -0.92599509,
                     0.01220563, -0.44161753, -0.47551061, -0.0789383 ,  0.12853926,
                     0.20014234, -0.04370486, -0.0252995 ,  0.11881716;

    QhatResult = BQhatD * D;

    if ( ( QhatAnswer2 - QhatResult ).norm( ) > 1e-6 * ( QhatAnswer2.norm( ) + 1 ) ){

        results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection (test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_initializeCoupling_averaged_l2_projection & True\n";
    return 0;
}

int test_overlapCoupling_getReferenceFreeMicroDomainMasses( std::ofstream &results ){
    /*!
     * Test the extraction of the reference free micro-domain mass
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string testName = "overlapCoupling_getReferenceFreeMicroDomainMasses";
    uIntType testNum = 0;

    std::string filename = "testConfig_averaged_l2_projection.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_getReferenceFreeMicroDomainMasses & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_getReferenceFreeMicroDomainMasses & False\n";
        return 1;
    }

    const cellDomainFloatMap referenceFreeMicroDomainMassesAnswer
        =
        { { 1,
              { 
                  { "free_nodeset_volume_1", 0.25 },
                  { "free_nodeset_volume_2", 0.25 },
                  { "free_nodeset_volume_3", 0.25 },
                  { "free_nodeset_volume_4", 0.25 },
                  { "free_nodeset_volume_5", 0.25 },
                  { "free_nodeset_volume_6", 0.25 },
                  { "free_nodeset_volume_7", 0.25 },
                  { "free_nodeset_volume_8", 0.25 },
              }
          }
        };

    const cellDomainFloatMap *referenceFreeMicroDomainMassesResult = oc.getReferenceFreeMicroDomainMasses( );

    if ( _compare_cellDomainMaps( results, referenceFreeMicroDomainMassesAnswer, *referenceFreeMicroDomainMassesResult,
                                  testName, testNum ) ){
        return 1;
    }

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_getReferenceFreeMicroDomainMasses & True\n";
    return 0;
}

int test_overlapCoupling_getReferenceGhostMicroDomainMasses( std::ofstream &results ){
    /*!
     * Test the extraction of the reference ghost micro-domain mass
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string testName = "overlapCoupling_getReferenceGhostMicroDomainMasses";
    uIntType testNum = 0;

    std::string filename = "testConfig_averaged_l2_projection.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_getReferenceGhostMicroDomainMasses & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_getReferenceGhostMicroDomainMasses & False\n";
        return 1;
    }

    const cellDomainFloatMap referenceGhostMicroDomainMassesAnswer
        =
        { { 2,
              { 
                  { "ghost_nodeset_volume_1", 0.25 },
                  { "ghost_nodeset_volume_2", 0.25 },
                  { "ghost_nodeset_volume_3", 0.25 },
                  { "ghost_nodeset_volume_4", 0.25 },
                  { "ghost_nodeset_volume_5", 0.25 },
                  { "ghost_nodeset_volume_6", 0.25 },
                  { "ghost_nodeset_volume_7", 0.25 },
                  { "ghost_nodeset_volume_8", 0.25 },
              }
          }
        };

    const cellDomainFloatMap *referenceGhostMicroDomainMassesResult = oc.getReferenceGhostMicroDomainMasses( );

    if ( _compare_cellDomainMaps( results, referenceGhostMicroDomainMassesAnswer, *referenceGhostMicroDomainMassesResult,
                                  testName, testNum ) ){
        return 1;
    }

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_getReferenceGhostMicroDomainMasses & True\n";
    return 0;
}
int test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass( std::ofstream &results ){
    /*!
     * Test the extraction of the reference free micro-domain centers of mass
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string testName = "overlapCoupling_getReferenceFreeMicroDomainCentersOfMass";
    uIntType testNum = 0;

    std::string filename = "testConfig_averaged_l2_projection.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass & False\n";
        return 1;
    }

    const cellDomainFloatVectorMap referenceFreeMicroDomainCentersOfMassAnswer
        =
        {
            { 1,
                {
                    { "free_nodeset_volume_1", { 0.75, 0.25, 2.75 } },
                    { "free_nodeset_volume_2", { 0.75, 0.25, 2.25 } },
                    { "free_nodeset_volume_3", { 0.25, 0.25, 2.75 } },
                    { "free_nodeset_volume_4", { 0.25, 0.25, 2.25 } },
                    { "free_nodeset_volume_5", { 0.75, 0.75, 2.75 } },
                    { "free_nodeset_volume_6", { 0.75, 0.75, 2.25 } },
                    { "free_nodeset_volume_7", { 0.25, 0.75, 2.75 } },
                    { "free_nodeset_volume_8", { 0.25, 0.75, 2.25 } }
                }
            }
        };

    const cellDomainFloatVectorMap *referenceFreeMicroDomainCentersOfMassResult = oc.getReferenceFreeMicroDomainCentersOfMass( );

    if ( _compare_cellDomainMaps( results, referenceFreeMicroDomainCentersOfMassAnswer, *referenceFreeMicroDomainCentersOfMassResult,
                                  testName, testNum ) ){
        return 1;
    }

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass & True\n";
    return 0;
}

int test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass( std::ofstream &results ){
    /*!
     * Test the extraction of the reference ghost micro-domain centers of mass
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string testName = "overlapCoupling_getReferenceGhostMicroDomainCentersOfMass";
    uIntType testNum = 0;

    std::string filename = "testConfig_averaged_l2_projection.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass & False\n";
        return 1;
    }

    const cellDomainFloatVectorMap referenceFreeMicroDomainCentersOfMassAnswer
        =
        {
            { 1,
                {
                    { "free_nodeset_volume_1", { 0.75, 0.25, 2.75 } },
                    { "free_nodeset_volume_2", { 0.75, 0.25, 2.25 } },
                    { "free_nodeset_volume_3", { 0.25, 0.25, 2.75 } },
                    { "free_nodeset_volume_4", { 0.25, 0.25, 2.25 } },
                    { "free_nodeset_volume_5", { 0.75, 0.75, 2.75 } },
                    { "free_nodeset_volume_6", { 0.75, 0.75, 2.25 } },
                    { "free_nodeset_volume_7", { 0.25, 0.75, 2.75 } },
                    { "free_nodeset_volume_8", { 0.25, 0.75, 2.25 } }
                }
            }
        };

    const cellDomainFloatVectorMap *referenceFreeMicroDomainCentersOfMassResult = oc.getReferenceFreeMicroDomainCentersOfMass( );

    if ( _compare_cellDomainMaps( results, referenceFreeMicroDomainCentersOfMassAnswer, *referenceFreeMicroDomainCentersOfMassResult,
                                  testName, testNum ) ){
        return 1;
    }

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass & True\n";
    return 0;
}

//int test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions( std::ofstream &results ){
//    /*!
//     * Test of the retrieval of the reference free micro domain centers of mass shape functions
//     */
//
//    remove( "reference_information.xdmf" );
//    remove( "reference_information.h5" );
//
//    remove( "homogenized_response.xdmf" );
//    remove( "homogenized_response.h5" );
//
//    std::string filename = "testConfig_averaged_l2_projection.yaml";
//    overlapCoupling::overlapCoupling oc( filename );
//
//    if ( oc.getConstructorError( ) ){
//        oc.getConstructorError( )->print( );
//        results << "test_overlapCoupling_ & False\n";
//        return 1;
//    }
//
//    errorOut error = oc.initializeCoupling( );
//
//    if ( error ){
//        error->print( );
//        results << "test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions & False\n";
//        return 1;
//    }
//
//    const floatVector referenceFreeMicroDomainCenterOfMassShapeFunctionsAnswer =
//        {
//            0.140625, 0.046875, 0.140625, 0.421875, 0.046875, 0.015625, 0.046875, 0.140625,
//            0.046875, 0.015625, 0.046875, 0.140625, 0.140625, 0.046875, 0.140625, 0.421875,
//            0.421875, 0.140625, 0.046875, 0.140625, 0.140625, 0.046875, 0.015625, 0.046875,
//            0.140625, 0.046875, 0.015625, 0.046875, 0.421875, 0.140625, 0.046875, 0.140625,
//            0.046875, 0.140625, 0.421875, 0.140625, 0.015625, 0.046875, 0.140625, 0.046875,
//            0.015625, 0.046875, 0.140625, 0.046875, 0.046875, 0.140625, 0.421875, 0.140625,
//            0.140625, 0.421875, 0.140625, 0.046875, 0.046875, 0.140625, 0.046875, 0.015625,
//            0.046875, 0.140625, 0.046875, 0.015625, 0.140625, 0.421875, 0.140625, 0.046875
//        };
//
//    const floatVector *referenceFreeMicroDomainCenterOfMassShapeFunctionsResult =
//        oc.getReferenceFreeMicroDomainCenterOfMassShapeFunctions( );
//
//    if ( !vectorTools::fuzzyEquals( referenceFreeMicroDomainCenterOfMassShapeFunctionsAnswer,
//                                   *referenceFreeMicroDomainCenterOfMassShapeFunctionsResult ) ){
//        vectorTools::print( *referenceFreeMicroDomainCenterOfMassShapeFunctionsResult );
//        results << "test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions (test 1) & False\n";
//        return 1;
//    } 
//
//    remove( "reference_information.xdmf" );
//    remove( "reference_information.h5" );
//
//    remove( "homogenized_response.xdmf" );
//    remove( "homogenized_response.h5" );
//
//    results << "test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions & True\n";
//    return 0;
//
//}
//
//int test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions( std::ofstream &results ){
//    /*!
//     * Test of the retrieval of the reference ghost micro domain centers of mass shape functions
//     */
//
//    remove( "reference_information.xdmf" );
//    remove( "reference_information.h5" );
//
//    remove( "homogenized_response.xdmf" );
//    remove( "homogenized_response.h5" );
//
//    std::string filename = "../testFiles/testConfig_averaged_l2_projection.yaml";
//    overlapCoupling::overlapCoupling oc( filename );
//
//    if ( oc.getConstructorError( ) ){
//        oc.getConstructorError( )->print( );
//        results << "test_overlapCoupling_ & False\n";
//        return 1;
//    }
//
//    errorOut error = oc.initializeCoupling( );
//
//    if ( error ){
//        error->print( );
//        results << "test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions & False\n";
//        return 1;
//    }
//
//    const floatVector referenceGhostMicroDomainCenterOfMassShapeFunctionsAnswer =
//        {
//            0.046875, 0.015625, 0.046875, 0.140625, 0.140625, 0.046875, 0.140625, 0.421875,
//            0.015625, 0.046875, 0.140625, 0.046875, 0.046875, 0.140625, 0.421875, 0.140625,
//            0.140625, 0.046875, 0.015625, 0.046875, 0.421875, 0.140625, 0.046875, 0.140625,
//            0.046875, 0.140625, 0.046875, 0.015625, 0.140625, 0.421875, 0.140625, 0.046875,
//            0.140625, 0.046875, 0.140625, 0.421875, 0.046875, 0.015625, 0.046875, 0.140625,
//            0.046875, 0.140625, 0.421875, 0.140625, 0.015625, 0.046875, 0.140625, 0.046875,
//            0.421875, 0.140625, 0.046875, 0.140625, 0.140625, 0.046875, 0.015625, 0.046875,
//            0.140625, 0.421875, 0.140625, 0.046875, 0.046875, 0.140625, 0.046875, 0.015625
//        };
//
//    const floatVector *referenceGhostMicroDomainCenterOfMassShapeFunctionsResult =
//        oc.getReferenceGhostMicroDomainCenterOfMassShapeFunctions( );
//
//    if ( !vectorTools::fuzzyEquals( referenceGhostMicroDomainCenterOfMassShapeFunctionsAnswer,
//                                   *referenceGhostMicroDomainCenterOfMassShapeFunctionsResult ) ){
//        vectorTools::print( *referenceGhostMicroDomainCenterOfMassShapeFunctionsResult );
//        results << "test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions (test 1) & False\n";
//        return 1;
//    } 
//
//    remove( "reference_information.xdmf" );
//    remove( "reference_information.h5" );
//
//    remove( "homogenized_response.xdmf" );
//    remove( "homogenized_response.h5" );
//
//    results << "test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions & True\n";
//    return 0;
//
//}

int test_overlapCoupling_processIncrement( std::ofstream &results ){
    /*!
     * Test the initialization of the coupling
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    remove( "macroscale_dof.xdmf" );
    remove( "macroscale_dof.h5" );

    remove( "microscale_dof.xdmf" );
    remove( "microscale_dof.h5" );

    std::string testName = "overlapCoupling_processIncrement";
    uIntType testNum = 0;

    std::string filename = "testConfig_averaged_l2_projection.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_processIncrement & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_processIncrement & False\n";
        return 1;
    }

    error = oc.processIncrement( 1, 1 );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_processIncrement & False\n";
        return 1;
    }

    //Check the mass properties
    std::unordered_map< std::string, floatType > freeDomainMassAnswer
        =
        {
            { "free_nodeset_volume_1", 0.250000 },
            { "free_nodeset_volume_2", 0.250000 },
            { "free_nodeset_volume_3", 0.250000 },
            { "free_nodeset_volume_4", 0.250000 },
            { "free_nodeset_volume_5", 0.375000 },
            { "free_nodeset_volume_6", 0.375000 },
            { "free_nodeset_volume_7", 0.375000 },
            { "free_nodeset_volume_8", 0.375000 },
        };

    std::unordered_map< std::string, floatType > ghostDomainMassAnswer
        =
        {
            { "ghost_nodeset_volume_1", 0.250000 },
            { "ghost_nodeset_volume_2", 0.250000 },
            { "ghost_nodeset_volume_3", 0.250000 },
            { "ghost_nodeset_volume_4", 0.250000 },
            { "ghost_nodeset_volume_5", 0.250000 },
            { "ghost_nodeset_volume_6", 0.250000 },
            { "ghost_nodeset_volume_7", 0.250000 },
            { "ghost_nodeset_volume_8", 0.250000 },
        };

    const std::unordered_map< std::string, floatType > *freeDomainMassResult = oc.getFreeMicroDomainMasses( );

    if( _compare_domainMaps( results, freeDomainMassAnswer, *freeDomainMassResult, testName, testNum ) ){

        return 1;

    }

    const std::unordered_map< std::string, floatType > *ghostDomainMassResult = oc.getGhostMicroDomainMasses( );

    if( _compare_domainMaps( results, ghostDomainMassAnswer, *ghostDomainMassResult, testName, testNum ) ){

        return 1;

    }

    std::unordered_map< std::string, floatVector > freeDomainCenterOfMassAnswer
        =
        {
            { "free_nodeset_volume_1", { 0.250000, 0.250000, 1.251000 } },
            { "free_nodeset_volume_2", { 0.750000, 0.250000, 1.251000 } },
            { "free_nodeset_volume_3", { 0.750000, 0.750000, 1.251000 } },
            { "free_nodeset_volume_4", { 0.250000, 0.750000, 1.251000 } },
            { "free_nodeset_volume_5", { 0.250000, 0.250000, 1.834333 } },
            { "free_nodeset_volume_6", { 0.750000, 0.250000, 1.834333 } },
            { "free_nodeset_volume_7", { 0.750000, 0.750000, 1.834333 } },
            { "free_nodeset_volume_8", { 0.250000, 0.750000, 1.834333 } },
        };
        
    std::unordered_map< std::string, floatVector > ghostDomainCenterOfMassAnswer
        =
        {
            { "ghost_nodeset_volume_1", { 0.250000, 0.250000, 0.251000 } },
            { "ghost_nodeset_volume_2", { 0.750000, 0.250000, 0.251000 } },
            { "ghost_nodeset_volume_3", { 0.750000, 0.750000, 0.251000 } },
            { "ghost_nodeset_volume_4", { 0.250000, 0.750000, 0.251000 } },
            { "ghost_nodeset_volume_5", { 0.250000, 0.250000, 0.751000 } },
            { "ghost_nodeset_volume_6", { 0.750000, 0.250000, 0.751000 } },
            { "ghost_nodeset_volume_7", { 0.750000, 0.750000, 0.751000 } },
            { "ghost_nodeset_volume_8", { 0.250000, 0.750000, 0.751000 } },
        };
        
    const std::unordered_map< std::string, floatVector > *freeDomainCenterOfMassResult = oc.getFreeMicroDomainCentersOfMass( );

    if( _compare_domainMaps( results, freeDomainCenterOfMassAnswer, *freeDomainCenterOfMassResult, testName, testNum ) ){

        return 1;

    }

    const std::unordered_map< std::string, floatVector > *ghostDomainCenterOfMassResult = oc.getGhostMicroDomainCentersOfMass( );

    if( _compare_domainMaps( results, ghostDomainCenterOfMassAnswer, *ghostDomainCenterOfMassResult, testName, testNum ) ){

        return 1;

    }

    //Test the initial projected displacements
    floatVector DhatAnswer =
        {
            0.00000000e+00,  0.00000000e+00,  9.73677106e-04,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00, -9.83716477e-20,  3.38388162e-19,  1.40010799e-03,
            0.00000000e+00,  0.00000000e+00,  9.73677106e-04,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00,  3.68675543e-19,  1.09620245e-19,  1.40010799e-03,
            0.00000000e+00,  0.00000000e+00,  9.73677106e-04,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00, -8.68532583e-20,  1.23294745e-19,  1.40010799e-03,
            0.00000000e+00,  0.00000000e+00,  9.73677106e-04,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00,  5.49603859e-20,  2.56747739e-19,  1.40010799e-03,
            0.00000000e+00,  0.00000000e+00,  1.03199244e-03,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00, -6.29037035e-21,  2.11103544e-19, -6.86285097e-04,
            0.00000000e+00,  0.00000000e+00,  1.03199244e-03,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00, -8.34160044e-19, -4.64924223e-19, -6.86285097e-04,
            0.00000000e+00,  0.00000000e+00,  1.03199244e-03,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00, -1.11438260e-19,  3.29118343e-19, -6.86285097e-04,
            0.00000000e+00,  0.00000000e+00,  1.03199244e-03,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00, -2.21234058e-19, -2.19907477e-19, -6.86285097e-04
        };

    if ( !vectorTools::fuzzyEquals( DhatAnswer, oc._test_initial_projected_ghost_macro_displacement ) ){

        results << testName + "(test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum++;

    floatVector QhatAnswer =
        {
            0.        ,  0.        , -0.00059409,  0.        ,  0.        ,
           -0.00010067,  0.        ,  0.        , -0.00059409,  0.        ,
            0.        , -0.00010067,  0.        ,  0.        , -0.00059409,
            0.        ,  0.        , -0.00010067,  0.        ,  0.        ,
           -0.00059409,  0.        ,  0.        , -0.00010067,  0.        ,
            0.        , -0.00059409,  0.        ,  0.        , -0.00010067,
            0.        ,  0.        , -0.00059409,  0.        ,  0.        ,
           -0.00010067,  0.        ,  0.        , -0.00059409,  0.        ,
            0.        , -0.00010067,  0.        ,  0.        , -0.00059409,
            0.        ,  0.        , -0.00010067,  0.        ,  0.        ,
           -0.00059409,  0.        ,  0.        , -0.00010067
        };

    if ( !vectorTools::fuzzyEquals( DhatAnswer, oc._test_initial_projected_ghost_macro_displacement ) ){

        results << testName + "(test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum++;

    const cellDomainFloatMap homogenizedVolumesAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", 0.125000 },
                    { "ghost_nodeset_volume_2", 0.125000 },
                    { "ghost_nodeset_volume_3", 0.125000 },
                    { "ghost_nodeset_volume_4", 0.125000 },
                    { "ghost_nodeset_volume_5", 0.125000 },
                    { "ghost_nodeset_volume_6", 0.125000 },
                    { "ghost_nodeset_volume_7", 0.125000 },
                    { "ghost_nodeset_volume_8", 0.125000 },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", 0.125000 },
                    { "free_nodeset_volume_2", 0.125000 },
                    { "free_nodeset_volume_3", 0.125000 },
                    { "free_nodeset_volume_4", 0.125000 },
                    { "free_nodeset_volume_5", 0.125000 },
                    { "free_nodeset_volume_6", 0.125000 },
                    { "free_nodeset_volume_7", 0.125000 },
                    { "free_nodeset_volume_8", 0.125000 },
                }
            }
        };

    //Note higher tolerance because it's an approximate volume reconstruction value
    if( _compare_cellDomainMaps( results, homogenizedVolumesAnswer, *oc.getHomogenizedVolumes( ), testName, testNum, 1e-6, 1e-3 ) ){

        return 1;

    }

    const cellDomainFloatMap homogenizedDensitiesAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", 2. },
                    { "ghost_nodeset_volume_2", 2. },
                    { "ghost_nodeset_volume_3", 2. },
                    { "ghost_nodeset_volume_4", 2. },
                    { "ghost_nodeset_volume_5", 2. },
                    { "ghost_nodeset_volume_6", 2. },
                    { "ghost_nodeset_volume_7", 2. },
                    { "ghost_nodeset_volume_8", 2. },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", 2. },
                    { "free_nodeset_volume_2", 2. },
                    { "free_nodeset_volume_3", 2. },
                    { "free_nodeset_volume_4", 2. },
                    { "free_nodeset_volume_5", 2. },
                    { "free_nodeset_volume_6", 2. },
                    { "free_nodeset_volume_7", 2. },
                    { "free_nodeset_volume_8", 2. },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedDensitiesAnswer, *oc.getHomogenizedDensities( ), testName, testNum ) ){

        return 1;

    }

    const cellDomainFloatVectorMap homogenizedSymmetricMicroStressesAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "ghost_nodeset_volume_2", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "ghost_nodeset_volume_3", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "ghost_nodeset_volume_4", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "ghost_nodeset_volume_5", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "ghost_nodeset_volume_6", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "ghost_nodeset_volume_7", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "ghost_nodeset_volume_8", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "free_nodeset_volume_2", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "free_nodeset_volume_3", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "free_nodeset_volume_4", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "free_nodeset_volume_5", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "free_nodeset_volume_6", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "free_nodeset_volume_7", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                    { "free_nodeset_volume_8", { 1.000000, 2.000000, 3.000000, 4.000000, 5.000000, 6.000000, 7.000000, 8.000000, 9.000000 } },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedSymmetricMicroStressesAnswer,
                                 *oc.getHomogenizedSymmetricMicroStresses( ), testName, testNum ) ){

        return 1;

    }

    const cellDomainFloatVectorMap homogenizedCentersOfMassAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { 0.250000, 0.250000, 0.251000 } },
                    { "ghost_nodeset_volume_2", { 0.750000, 0.250000, 0.251000 } },
                    { "ghost_nodeset_volume_3", { 0.750000, 0.750000, 0.251000 } },
                    { "ghost_nodeset_volume_4", { 0.250000, 0.750000, 0.251000 } },
                    { "ghost_nodeset_volume_5", { 0.250000, 0.250000, 0.751000 } },
                    { "ghost_nodeset_volume_6", { 0.750000, 0.250000, 0.751000 } },
                    { "ghost_nodeset_volume_7", { 0.750000, 0.750000, 0.751000 } },
                    { "ghost_nodeset_volume_8", { 0.250000, 0.750000, 0.751000 } },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { 0.250000, 0.250000, 1.251000 } },
                    { "free_nodeset_volume_2", { 0.750000, 0.250000, 1.251000 } },
                    { "free_nodeset_volume_3", { 0.750000, 0.750000, 1.251000 } },
                    { "free_nodeset_volume_4", { 0.250000, 0.750000, 1.251000 } },
                    { "free_nodeset_volume_5", { 0.250000, 0.250000, 1.751000 } },
                    { "free_nodeset_volume_6", { 0.750000, 0.250000, 1.751000 } },
                    { "free_nodeset_volume_7", { 0.750000, 0.750000, 1.751000 } },
                    { "free_nodeset_volume_8", { 0.250000, 0.750000, 1.751000 } },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedCentersOfMassAnswer,
                                 *oc.getHomogenizedCentersOfMass( ), testName, testNum ) ){

        return 1;

    }

    const cellDomainFloatVectorMap homogenizedBodyForcesAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { -1.000000, -2.000000, -3.000000 } },
                    { "ghost_nodeset_volume_2", { -1.000000, -2.000000, -3.000000 } },
                    { "ghost_nodeset_volume_3", { -1.000000, -2.000000, -3.000000 } },
                    { "ghost_nodeset_volume_4", { -1.000000, -2.000000, -3.000000 } },
                    { "ghost_nodeset_volume_5", { -1.000000, -2.000000, -3.000000 } },
                    { "ghost_nodeset_volume_6", { -1.000000, -2.000000, -3.000000 } },
                    { "ghost_nodeset_volume_7", { -1.000000, -2.000000, -3.000000 } },
                    { "ghost_nodeset_volume_8", { -1.000000, -2.000000, -3.000000 } },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { -1.000000, -2.000000, -3.000000 } },
                    { "free_nodeset_volume_2", { -1.000000, -2.000000, -3.000000 } },
                    { "free_nodeset_volume_3", { -1.000000, -2.000000, -3.000000 } },
                    { "free_nodeset_volume_4", { -1.000000, -2.000000, -3.000000 } },
                    { "free_nodeset_volume_5", { -1.000000, -2.000000, -3.000000 } },
                    { "free_nodeset_volume_6", { -1.000000, -2.000000, -3.000000 } },
                    { "free_nodeset_volume_7", { -1.000000, -2.000000, -3.000000 } },
                    { "free_nodeset_volume_8", { -1.000000, -2.000000, -3.000000 } },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedBodyForcesAnswer,
                                 *oc.getHomogenizedBodyForces( ), testName, testNum ) ){

        return 1;

    }

    const cellDomainFloatVectorMap homogenizedAccelerationsAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { 0.000000, 0.000000, 0.003000 } },
                    { "ghost_nodeset_volume_2", { 0.000000, 0.000000, 0.003000 } },
                    { "ghost_nodeset_volume_3", { 0.000000, 0.000000, 0.003000 } },
                    { "ghost_nodeset_volume_4", { 0.000000, 0.000000, 0.003000 } },
                    { "ghost_nodeset_volume_5", { 0.000000, 0.000000, 0.003000 } },
                    { "ghost_nodeset_volume_6", { 0.000000, 0.000000, 0.003000 } },
                    { "ghost_nodeset_volume_7", { 0.000000, 0.000000, 0.003000 } },
                    { "ghost_nodeset_volume_8", { 0.000000, 0.000000, 0.003000 } }
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { 0.000000, 0.000000, 0.003000 } },
                    { "free_nodeset_volume_2", { 0.000000, 0.000000, 0.003000 } },
                    { "free_nodeset_volume_3", { 0.000000, 0.000000, 0.003000 } },
                    { "free_nodeset_volume_4", { 0.000000, 0.000000, 0.003000 } },
                    { "free_nodeset_volume_5", { 0.000000, 0.000000, 0.003000 } },
                    { "free_nodeset_volume_6", { 0.000000, 0.000000, 0.003000 } },
                    { "free_nodeset_volume_7", { 0.000000, 0.000000, 0.003000 } },
                    { "free_nodeset_volume_8", { 0.000000, 0.000000, 0.003000 } }
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedAccelerationsAnswer,
                                 *oc.getHomogenizedAccelerations( ), testName, testNum ) ){

        return 1;

    }

    //NOTE: The different free nodesets have come directly from the code output. I believe they are correct but they
    //      were not computed outside of the code.
    const cellDomainFloatVectorMap homogenizedMicroInertiasAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "ghost_nodeset_volume_2", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "ghost_nodeset_volume_3", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "ghost_nodeset_volume_4", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "ghost_nodeset_volume_5", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "ghost_nodeset_volume_6", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "ghost_nodeset_volume_7", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "ghost_nodeset_volume_8", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "free_nodeset_volume_2", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "free_nodeset_volume_3", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "free_nodeset_volume_4", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "free_nodeset_volume_5", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "free_nodeset_volume_6", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "free_nodeset_volume_7", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                    { "free_nodeset_volume_8", { 0.062500, 0.000000, 0.000000, 0.000000, 0.062500, 0.000000, 0.000000, 0.000000, 0.062500 } },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedMicroInertiasAnswer,
                                 *oc.getHomogenizedMicroInertias( ), testName, testNum ) ){

        return 1;

    }

    const cellDomainFloatVectorMap homogenizedBodyForceCouplesAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_2", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_3", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_4", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_5", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_6", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_7", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_8", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_2", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_3", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_4", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_5", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_6", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_7", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_8", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedBodyForceCouplesAnswer,
                                 *oc.getHomogenizedBodyForceCouples( ), testName, testNum ) ){

        return 1;

    }

    const cellDomainFloatVectorMap homogenizedMicroSpinInertiasAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_2", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_3", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_4", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_5", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_6", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_7", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "ghost_nodeset_volume_8", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_2", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_3", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_4", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_5", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_6", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_7", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                    { "free_nodeset_volume_8", { 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedMicroSpinInertiasAnswer,
                                 *oc.getHomogenizedMicroSpinInertias( ), testName, testNum ) ){

        return 1;

    }

    const cellDomainFloatMap homogenizedSurfaceAreasAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", 1.5 },
                    { "ghost_nodeset_volume_2", 1.5 },
                    { "ghost_nodeset_volume_3", 1.5 },
                    { "ghost_nodeset_volume_4", 1.5 },
                    { "ghost_nodeset_volume_5", 1.5 },
                    { "ghost_nodeset_volume_6", 1.5 },
                    { "ghost_nodeset_volume_7", 1.5 },
                    { "ghost_nodeset_volume_8", 1.5 },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", 1.5 },
                    { "free_nodeset_volume_2", 1.5 },
                    { "free_nodeset_volume_3", 1.5 },
                    { "free_nodeset_volume_4", 1.5 },
                    { "free_nodeset_volume_5", 1.5 },
                    { "free_nodeset_volume_6", 1.5 },
                    { "free_nodeset_volume_7", 1.5 },
                    { "free_nodeset_volume_8", 1.5 },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedSurfaceAreasAnswer,
                                 *oc.getHomogenizedSurfaceAreas( ), testName, testNum, 1e-6, 1e-2 ) ){

        return 1;

    }

    const cellDomainUIntVectorMap cellDomainMacroSurfacesAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { 0, 2, 4 } },
                    { "ghost_nodeset_volume_2", { 1, 2, 4 } },
                    { "ghost_nodeset_volume_3", { 1, 3, 4 } },
                    { "ghost_nodeset_volume_4", { 0, 3, 4 } },
                    { "ghost_nodeset_volume_5", { 0, 2, 5 } },
                    { "ghost_nodeset_volume_6", { 1, 2, 5 } },
                    { "ghost_nodeset_volume_7", { 1, 3, 5 } },
                    { "ghost_nodeset_volume_8", { 0, 3, 5 } },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { 0, 2, 4 } },
                    { "free_nodeset_volume_2", { 1, 2, 4 } },
                    { "free_nodeset_volume_3", { 1, 3, 4 } },
                    { "free_nodeset_volume_4", { 0, 3, 4 } },
                    { "free_nodeset_volume_5", { 0, 2, 5 } },
                    { "free_nodeset_volume_6", { 1, 2, 5 } },
                    { "free_nodeset_volume_7", { 1, 3, 5 } },
                    { "free_nodeset_volume_8", { 0, 3, 5 } },
                }
            }
        };

    if( _compare_cellDomainMaps( results, cellDomainMacroSurfacesAnswer,
                                 *oc.getCellDomainMacroSurfaces( ), testName, testNum ) ){

        return 1;

    }


    const cellDomainFloatVectorMap homogenizedSurfaceRegionAreasAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { 0.25, 0.00, 0.25, 0.00, 0.25, 0.00 } },
                    { "ghost_nodeset_volume_2", { 0.00, 0.25, 0.25, 0.00, 0.25, 0.00 } },
                    { "ghost_nodeset_volume_3", { 0.00, 0.25, 0.00, 0.25, 0.25, 0.00 } },
                    { "ghost_nodeset_volume_4", { 0.25, 0.00, 0.00, 0.25, 0.25, 0.00 } },
                    { "ghost_nodeset_volume_5", { 0.25, 0.00, 0.25, 0.00, 0.00, 0.25 } },
                    { "ghost_nodeset_volume_6", { 0.00, 0.25, 0.25, 0.00, 0.00, 0.25 } },
                    { "ghost_nodeset_volume_7", { 0.00, 0.25, 0.00, 0.25, 0.00, 0.25 } },
                    { "ghost_nodeset_volume_8", { 0.25, 0.00, 0.00, 0.25, 0.00, 0.25 } },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { 0.25, 0.00, 0.25, 0.00, 0.25, 0.00 } },
                    { "free_nodeset_volume_2", { 0.00, 0.25, 0.25, 0.00, 0.25, 0.00 } },
                    { "free_nodeset_volume_3", { 0.00, 0.25, 0.00, 0.25, 0.25, 0.00 } },
                    { "free_nodeset_volume_4", { 0.25, 0.00, 0.00, 0.25, 0.25, 0.00 } },
                    { "free_nodeset_volume_5", { 0.25, 0.00, 0.25, 0.00, 0.00, 0.25 } },
                    { "free_nodeset_volume_6", { 0.00, 0.25, 0.25, 0.00, 0.00, 0.25 } },
                    { "free_nodeset_volume_7", { 0.00, 0.25, 0.00, 0.25, 0.00, 0.25 } },
                    { "free_nodeset_volume_8", { 0.25, 0.00, 0.00, 0.25, 0.00, 0.25 } },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedSurfaceRegionAreasAnswer,
                                 *oc.getHomogenizedSurfaceRegionAreas( ), testName, testNum, 1e-6, 1e-2 ) ){

        return 1;

    }

    const cellDomainFloatVectorMap homogenizedSurfaceRegionCentersOfMassAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { 0.00, 0.25, 0.251,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 0.00, 0.251,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 0.25, 0.001,
                                                  0.00, 0.00, 0.000 } },
                    { "ghost_nodeset_volume_2", { 0.00, 0.00, 0.000,
                                                  1.00, 0.25, 0.251,
                                                  0.75, 0.00, 0.251,
                                                  0.00, 0.00, 0.000,
                                                  0.75, 0.25, 0.001,
                                                  0.00, 0.00, 0.000 } },
                    { "ghost_nodeset_volume_3", { 0.00, 0.00, 0.000,
                                                  1.00, 0.75, 0.251,
                                                  0.00, 0.00, 0.000,
                                                  0.75, 1.00, 0.251,
                                                  0.75, 0.75, 0.001,
                                                  0.00, 0.00, 0.000 } },
                    { "ghost_nodeset_volume_4", { 0.00, 0.75, 0.251,
                                                  0.00, 0.00, 0.000,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 1.00, 0.251,
                                                  0.25, 0.75, 0.001,
                                                  0.00, 0.00, 0.000 } },
                    { "ghost_nodeset_volume_5", { 0.00, 0.25, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 0.00, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 0.25, 1.001 } },
                    { "ghost_nodeset_volume_6", { 0.00, 0.00, 0.000,
                                                  1.00, 0.25, 0.751,
                                                  0.75, 0.00, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.00, 0.00, 0.000,
                                                  0.75, 0.25, 1.001 } },
                    { "ghost_nodeset_volume_7", { 0.00, 0.00, 0.000,
                                                  1.00, 0.75, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.75, 1.00, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.75, 0.75, 1.001 } },
                    { "ghost_nodeset_volume_8", { 0.00, 0.75, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 1.00, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 0.75, 1.001 } },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { 0.00, 0.25, 1.251,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 0.00, 1.251,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 0.25, 1.001,
                                                 0.00, 0.00, 0.000 } },
                    { "free_nodeset_volume_2", { 0.00, 0.00, 0.000,
                                                 1.00, 0.25, 1.251,
                                                 0.75, 0.00, 1.251,
                                                 0.00, 0.00, 0.000,
                                                 0.75, 0.25, 1.001,
                                                 0.00, 0.00, 0.000 } },
                    { "free_nodeset_volume_3", { 0.00, 0.00, 0.000,
                                                 1.00, 0.75, 1.251,
                                                 0.00, 0.00, 0.000,
                                                 0.75, 1.00, 1.251,
                                                 0.75, 0.75, 1.001,
                                                 0.00, 0.00, 0.000 } },
                    { "free_nodeset_volume_4", { 0.00, 0.75, 1.251,
                                                 0.00, 0.00, 0.000,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 1.00, 1.251,
                                                 0.25, 0.75, 1.001,
                                                 0.00, 0.00, 0.000 } },
                    { "free_nodeset_volume_5", { 0.00, 0.25, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 0.00, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 0.25, 2.001 } },
                    { "free_nodeset_volume_6", { 0.00, 0.00, 0.000,
                                                 1.00, 0.25, 1.751,
                                                 0.75, 0.00, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.00, 0.00, 0.000,
                                                 0.75, 0.25, 2.001 } },
                    { "free_nodeset_volume_7", { 0.00, 0.00, 0.000,
                                                 1.00, 0.75, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.75, 1.00, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.75, 0.75, 2.001 } },
                    { "free_nodeset_volume_8", { 0.00, 0.75, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 1.00, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 0.75, 2.001 } },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedSurfaceRegionCentersOfMassAnswer,
                                 *oc.getHomogenizedSurfaceRegionCentersOfMass( ), testName, testNum, 1e-6, 1e-6 ) ){

        return 1;

    }

    const cellDomainFloatVectorMap homogenizedSurfaceRegionProjectedLocalCentersOfMassAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { -1.00, -0.50, -0.496,
                                                   0.00,  0.00,  0.000,
                                                  -0.50, -1.00, -0.496,
                                                   0.00,  0.00,  0.000,
                                                  -0.50, -0.50, -1.000,
                                                   0.00,  0.00,  0.000 } },
                    { "ghost_nodeset_volume_2", {  0.00,  0.00,  0.000,
                                                   1.00, -0.50, -0.496,
                                                   0.50, -1.00, -0.496,
                                                   0.00,  0.00,  0.000,
                                                   0.50, -0.50, -1.000,
                                                   0.00,  0.00,  0.000 } },
                    { "ghost_nodeset_volume_3", {  0.00,  0.00,  0.000,
                                                   1.00,  0.50, -0.496,
                                                   0.00,  0.00,  0.000,
                                                   0.50,  1.00, -0.496,
                                                   0.50,  0.50, -1.000,
                                                   0.00,  0.00,  0.000 } },
                    { "ghost_nodeset_volume_4", { -1.00,  0.50, -0.496,
                                                   0.00,  0.00,  0.000,
                                                   0.00,  0.00,  0.000,
                                                  -0.50,  1.00, -0.496,
                                                  -0.50,  0.50, -1.000,
                                                   0.00,  0.00,  0.000 } },
                    { "ghost_nodeset_volume_5", { -1.00, -0.50,  0.504,
                                                   0.00,  0.00,  0.000,
                                                  -0.50, -1.00,  0.504,
                                                   0.00,  0.00,  0.000,
                                                   0.00,  0.00,  0.000,
                                                  -0.50, -0.50,  1.000 } },
                    { "ghost_nodeset_volume_6", {  0.00,  0.00,  0.000,
                                                   1.00, -0.50,  0.504,
                                                   0.50, -1.00,  0.504,
                                                   0.00,  0.00,  0.000,
                                                   0.00,  0.00,  0.000,
                                                   0.50, -0.50,  1.000 } },
                    { "ghost_nodeset_volume_7", {  0.00,  0.00,  0.000,
                                                   1.00,  0.50,  0.504,
                                                   0.00,  0.00,  0.000,
                                                   0.50,  1.00,  0.504,
                                                   0.00,  0.00,  0.000,
                                                   0.50,  0.50,  1.000 } },
                    { "ghost_nodeset_volume_8", { -1.00,  0.50,  0.504,
                                                   0.00,  0.00,  0.000,
                                                   0.00,  0.00,  0.000,
                                                  -0.50,  1.00,  0.504,
                                                   0.00,  0.00,  0.000,
                                                  -0.50,  0.50,  1.000 } },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { -1.00, -0.50, -0.496,
                                                  0.00,  0.00,  0.000,
                                                 -0.50, -1.00, -0.496,
                                                  0.00,  0.00,  0.000,
                                                 -0.50, -0.50, -1.000,
                                                  0.00,  0.00,  0.000 } },
                    { "free_nodeset_volume_2", {  0.00,  0.00,  0.000,
                                                  1.00, -0.50, -0.496,
                                                  0.50, -1.00, -0.496,
                                                  0.00,  0.00,  0.000,
                                                  0.50, -0.50, -1.000,
                                                  0.00,  0.00,  0.000 } },
                    { "free_nodeset_volume_3", {  0.00,  0.00,  0.000,
                                                  1.00,  0.50, -0.496,
                                                  0.00,  0.00,  0.000,
                                                  0.50,  1.00, -0.496,
                                                  0.50,  0.50, -1.000,
                                                  0.00,  0.00,  0.000 } },
                    { "free_nodeset_volume_4", { -1.00,  0.50, -0.496,
                                                  0.00,  0.00,  0.000,
                                                  0.00,  0.00,  0.000,
                                                 -0.50,  1.00, -0.496,
                                                 -0.50,  0.50, -1.000,
                                                  0.00,  0.00,  0.000 } },
                    { "free_nodeset_volume_5", { -1.00, -0.50,  0.504,
                                                  0.00,  0.00,  0.000,
                                                 -0.50, -1.00,  0.504,
                                                  0.00,  0.00,  0.000,
                                                  0.00,  0.00,  0.000,
                                                 -0.50, -0.50,  1.000 } },
                    { "free_nodeset_volume_6", {  0.00,  0.00,  0.000,
                                                  1.00, -0.50,  0.504,
                                                  0.50, -1.00,  0.504,
                                                  0.00,  0.00,  0.000,
                                                  0.00,  0.00,  0.000,
                                                  0.50, -0.50,  1.000 } },
                    { "free_nodeset_volume_7", {  0.00,  0.00,  0.000,
                                                  1.00,  0.50,  0.504,
                                                  0.00,  0.00,  0.000,
                                                  0.50,  1.00,  0.504,
                                                  0.00,  0.00,  0.000,
                                                  0.50,  0.50,  1.000 } },
                    { "free_nodeset_volume_8", { -1.00,  0.50,  0.504,
                                                  0.00,  0.00,  0.000,
                                                  0.00,  0.00,  0.000,
                                                 -0.50,  1.00,  0.504,
                                                  0.00,  0.00,  0.000,
                                                 -0.50,  0.50,  1.000 } },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedSurfaceRegionProjectedLocalCentersOfMassAnswer,
                                 *oc.getHomogenizedSurfaceRegionProjectedLocalCentersOfMass( ), testName, testNum, 1e-6, 1e-3 ) ){

        return 1;

    }

    const cellDomainFloatVectorMap homogenizedSurfaceRegionProjectedCentersOfMassAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { 0.00, 0.25, 0.251,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 0.00, 0.251,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 0.25,-0.001,
                                                  0.00, 0.00, 0.000 } },
                    { "ghost_nodeset_volume_2", { 0.00, 0.00, 0.000,
                                                  1.00, 0.25, 0.251,
                                                  0.75, 0.00, 0.251,
                                                  0.00, 0.00, 0.000,
                                                  0.75, 0.25,-0.001,
                                                  0.00, 0.00, 0.000 } },
                    { "ghost_nodeset_volume_3", { 0.00, 0.00, 0.000,
                                                  1.00, 0.75, 0.251,
                                                  0.00, 0.00, 0.000,
                                                  0.75, 1.00, 0.251,
                                                  0.75, 0.75,-0.001,
                                                  0.00, 0.00, 0.000 } },
                    { "ghost_nodeset_volume_4", { 0.00, 0.75, 0.251,
                                                  0.00, 0.00, 0.000,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 1.00, 0.251,
                                                  0.25, 0.75,-0.001,
                                                  0.00, 0.00, 0.000 } },
                    { "ghost_nodeset_volume_5", { 0.00, 0.25, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 0.00, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 0.25, 0.999 } },
                    { "ghost_nodeset_volume_6", { 0.00, 0.00, 0.000,
                                                  1.00, 0.25, 0.751,
                                                  0.75, 0.00, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.00, 0.00, 0.000,
                                                  0.75, 0.25, 0.999 } },
                    { "ghost_nodeset_volume_7", { 0.00, 0.00, 0.000,
                                                  1.00, 0.75, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.75, 1.00, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.75, 0.75, 0.999 } },
                    { "ghost_nodeset_volume_8", { 0.00, 0.75, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 1.00, 0.751,
                                                  0.00, 0.00, 0.000,
                                                  0.25, 0.75, 0.999 } },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { 0.00, 0.25, 1.251,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 0.00, 1.251,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 0.25, 0.999,
                                                 0.00, 0.00, 0.000 } },
                    { "free_nodeset_volume_2", { 0.00, 0.00, 0.000,
                                                 1.00, 0.25, 1.251,
                                                 0.75, 0.00, 1.251,
                                                 0.00, 0.00, 0.000,
                                                 0.75, 0.25, 0.999,
                                                 0.00, 0.00, 0.000 } },
                    { "free_nodeset_volume_3", { 0.00, 0.00, 0.000,
                                                 1.00, 0.75, 1.251,
                                                 0.00, 0.00, 0.000,
                                                 0.75, 1.00, 1.251,
                                                 0.75, 0.75, 0.999,
                                                 0.00, 0.00, 0.000 } },
                    { "free_nodeset_volume_4", { 0.00, 0.75, 1.251,
                                                 0.00, 0.00, 0.000,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 1.00, 1.251,
                                                 0.25, 0.75, 0.999,
                                                 0.00, 0.00, 0.000 } },
                    { "free_nodeset_volume_5", { 0.00, 0.25, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 0.00, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 0.25, 1.999 } },
                    { "free_nodeset_volume_6", { 0.00, 0.00, 0.000,
                                                 1.00, 0.25, 1.751,
                                                 0.75, 0.00, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.00, 0.00, 0.000,
                                                 0.75, 0.25, 1.999 } },
                    { "free_nodeset_volume_7", { 0.00, 0.00, 0.000,
                                                 1.00, 0.75, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.75, 1.00, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.75, 0.75, 1.999 } },
                    { "free_nodeset_volume_8", { 0.00, 0.75, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 1.00, 1.751,
                                                 0.00, 0.00, 0.000,
                                                 0.25, 0.75, 1.999 } },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedSurfaceRegionProjectedCentersOfMassAnswer,
                                 *oc.getHomogenizedSurfaceRegionProjectedCentersOfMass( ), testName, testNum, 1e-6, 1e-6 ) ){

        return 1;

    }

    const cellDomainFloatVectorMap homogenizedSurfaceRegionTractionsAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", { -1, -2, -3,
                                                   0,  0,  0,
                                                  -4, -5, -6,
                                                   0,  0,  0,
                                                  -7, -8, -9,
                                                   0,  0,  0 } },
                    { "ghost_nodeset_volume_2", {  0,  0,  0,
                                                   1,  2,  3,
                                                  -4, -5, -6,
                                                   0,  0,  0,
                                                  -7, -8, -9,
                                                   0,  0,  0 } },
                    { "ghost_nodeset_volume_3", {  0,  0,  0,
                                                   1,  2,  3,
                                                   0,  0,  0,
                                                   4,  5,  6,
                                                  -7, -8, -9,
                                                   0,  0,  0 } },
                    { "ghost_nodeset_volume_4", { -1, -2, -3,
                                                   0,  0,  0,
                                                   0,  0,  0,
                                                   4,  5,  6,
                                                  -7, -8, -9,
                                                   0,  0,  0 } },
                    { "ghost_nodeset_volume_5", { -1, -2, -3,
                                                   0,  0,  0,
                                                  -4, -5, -6,
                                                   0,  0,  0,
                                                   0,  0,  0,
                                                   7,  8,  9 } },
                    { "ghost_nodeset_volume_6", {  0,  0,  0,
                                                   1,  2,  3,
                                                  -4, -5, -6,
                                                   0,  0,  0,
                                                   0,  0,  0,
                                                   7,  8,  9 } },
                    { "ghost_nodeset_volume_7", {  0,  0,  0,
                                                   1,  2,  3,
                                                   0,  0,  0,
                                                   4,  5,  6,
                                                   0,  0,  0,
                                                   7,  8,  9 } },
                    { "ghost_nodeset_volume_8", { -1, -2, -3,
                                                   0,  0,  0,
                                                   0,  0,  0,
                                                   4,  5,  6,
                                                   0,  0,  0,
                                                   7,  8,  9 } },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", { -1, -2, -3,
                                                  0,  0,  0,
                                                 -4, -5, -6,
                                                  0,  0,  0,
                                                 -7, -8, -9,
                                                  0,  0,  0 } },
                    { "free_nodeset_volume_2", {  0,  0,  0,
                                                  1,  2,  3,
                                                 -4, -5, -6,
                                                  0,  0,  0,
                                                 -7, -8, -9,
                                                  0,  0,  0 } },
                    { "free_nodeset_volume_3", {  0,  0,  0,
                                                  1,  2,  3,
                                                  0,  0,  0,
                                                  4,  5,  6,
                                                 -7, -8, -9,
                                                  0,  0,  0 } },
                    { "free_nodeset_volume_4", { -1, -2, -3,
                                                  0,  0,  0,
                                                  0,  0,  0,
                                                  4,  5,  6,
                                                 -7, -8, -9,
                                                  0,  0,  0 } },
                    { "free_nodeset_volume_5", { -1, -2, -3,
                                                  0,  0,  0,
                                                 -4, -5, -6,
                                                  0,  0,  0,
                                                  0,  0,  0,
                                                  7,  8,  9 } },
                    { "free_nodeset_volume_6", {  0,  0,  0,
                                                  1,  2,  3,
                                                 -4, -5, -6,
                                                  0,  0,  0,
                                                  0,  0,  0,
                                                  7,  8,  9 } },
                    { "free_nodeset_volume_7", {  0,  0,  0,
                                                  1,  2,  3,
                                                  0,  0,  0,
                                                  4,  5,  6,
                                                  0,  0,  0,
                                                  7,  8,  9 } },
                    { "free_nodeset_volume_8", { -1, -2, -3,
                                                  0,  0,  0,
                                                  0,  0,  0,
                                                  4,  5,  6,
                                                  0,  0,  0,
                                                  7,  8,  9 } },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedSurfaceRegionTractionsAnswer,
                                 *oc.getHomogenizedSurfaceRegionTractions( ), testName, testNum, 1e-6, 1e-6 ) ){

        return 1;

    }

    const cellDomainFloatVectorMap homogenizedSurfaceRegionCouplesAnswer
        =
        {
            { 1,
                {
                    { "ghost_nodeset_volume_1", floatVector( 6 * 9, 0 ) },
                    { "ghost_nodeset_volume_2", floatVector( 6 * 9, 0 ) },
                    { "ghost_nodeset_volume_3", floatVector( 6 * 9, 0 ) },
                    { "ghost_nodeset_volume_4", floatVector( 6 * 9, 0 ) },
                    { "ghost_nodeset_volume_5", floatVector( 6 * 9, 0 ) },
                    { "ghost_nodeset_volume_6", floatVector( 6 * 9, 0 ) },
                    { "ghost_nodeset_volume_7", floatVector( 6 * 9, 0 ) },
                    { "ghost_nodeset_volume_8", floatVector( 6 * 9, 0 ) },
                }
            },
            { 2,
                {
                    { "free_nodeset_volume_1", floatVector( 6 * 9, 0 ) },
                    { "free_nodeset_volume_2", floatVector( 6 * 9, 0 ) },
                    { "free_nodeset_volume_3", floatVector( 6 * 9, 0 ) },
                    { "free_nodeset_volume_4", floatVector( 6 * 9, 0 ) },
                    { "free_nodeset_volume_5", floatVector( 6 * 9, 0 ) },
                    { "free_nodeset_volume_6", floatVector( 6 * 9, 0 ) },
                    { "free_nodeset_volume_7", floatVector( 6 * 9, 0 ) },
                    { "free_nodeset_volume_8", floatVector( 6 * 9, 0 ) },
                }
            }
        };

    if( _compare_cellDomainMaps( results, homogenizedSurfaceRegionCouplesAnswer,
                                 *oc.getHomogenizedSurfaceRegionCouples( ), testName, testNum, 1e-6, 1e-2 ) ){

        return 1;

    }

    floatVector elementNodalVolumesAnswer( 8, .125 );

    for ( auto c = oc._test_elementNodalVolumes.begin( ); c != oc._test_elementNodalVolumes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( elementNodalVolumesAnswer, c->second, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    floatVector volumeAtNodesAnswer( 8, .125 );

    for ( auto c = oc._test_volumeAtNodes.begin( ); c != oc._test_volumeAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( volumeAtNodesAnswer, c->second, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    floatVector densityAtNodesAnswer( 8, 2 );

    for ( auto c = oc._test_densityAtNodes.begin( ); c != oc._test_densityAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( densityAtNodesAnswer, c->second, 1e-6, 1e-2 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    floatVector tmp = { -1, -2, -3 };
    floatMatrix bodyForceAtNodesAnswer( 8, tmp );

    for ( auto c = oc._test_bodyForceAtNodes.begin( ); c != oc._test_bodyForceAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( bodyForceAtNodesAnswer, c->second, 1e-6, 1e-6 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    tmp = { 0., 0., 3e-3 };
    floatMatrix accelerationAtNodesAnswer( 8, tmp );

    for ( auto c = oc._test_accelerationAtNodes.begin( ); c != oc._test_accelerationAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( accelerationAtNodesAnswer, c->second, 1e-6, 1e-6 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    tmp = { 0.0625, 0, 0, 0, 0.0625, 0, 0, 0, 0.0625 };

    floatMatrix microInertiaAtNodesAnswer( 8, tmp );

    for ( auto c = oc._test_microInertiaAtNodes.begin( ); c != oc._test_microInertiaAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( microInertiaAtNodesAnswer, c->second, 1e-6, 1e-6 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    tmp = floatVector( 9, 0 );
    floatMatrix bodyCoupleAtNodesAnswer( 8, tmp );

    for ( auto c = oc._test_bodyCoupleAtNodes.begin( ); c != oc._test_bodyCoupleAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( bodyCoupleAtNodesAnswer, c->second, 1e-6, 1e-6 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    tmp = floatVector( 9, 0 );
    floatMatrix microSpinInertiaAtNodesAnswer( 8, tmp );

    for ( auto c = oc._test_microSpinInertiaAtNodes.begin( ); c != oc._test_microSpinInertiaAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( microSpinInertiaAtNodesAnswer, c->second, 1e-6, 1e-6 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    tmp = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    floatMatrix symmetricMicroStressAtNodesAnswer( 8, tmp );

    for ( auto c = oc._test_symmetricMicroStressAtNodes.begin( ); c != oc._test_symmetricMicroStressAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( symmetricMicroStressAtNodesAnswer, c->second, 1e-6, 1e-6 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    //Check the force vectors
    floatVector forceVectorsAnswer
        = { -3.25   , -4.25   , -5.25075, -2.75   , -3.25   , -3.75075,
            -0.75   , -0.75   , -0.75075, -1.25   , -1.75   , -2.25075,
             0.25   , -0.25   , -0.75075,  0.75   ,  0.75   ,  0.74925,
             2.75   ,  3.25   ,  3.74925,  2.25   ,  2.25   ,  2.24925 };

    floatVector extForceVectorsAnswer
        = { -3.25   , -4.25   , -5.25, -2.75   , -3.25   , -3.75,
            -0.75   , -0.75   , -0.75, -1.25   , -1.75   , -2.25,
             0.25   , -0.25   , -0.75,  0.75   ,  0.75   ,  0.75,
             2.75   ,  3.25   ,  3.74,  2.25   ,  2.25   ,  2.25 };

    floatVector coupleVectorsAnswer
        = { -0.125, -0.25 , -0.375, -0.5  , -0.625, -0.75 , -0.875, -1.   ,
            -1.125, -0.125, -0.25 , -0.375, -0.5  , -0.625, -0.75 , -0.875,
            -1.   , -1.125, -0.125, -0.25 , -0.375, -0.5  , -0.625, -0.75 ,
            -0.875, -1.   , -1.125, -0.125, -0.25 , -0.375, -0.5  , -0.625,
            -0.75 , -0.875, -1.   , -1.125, -0.125, -0.25 , -0.375, -0.5  ,
            -0.625, -0.75 , -0.875, -1.   , -1.125, -0.125, -0.25 , -0.375,
            -0.5  , -0.625, -0.75 , -0.875, -1.   , -1.125, -0.125, -0.25 ,
            -0.375, -0.5  , -0.625, -0.75 , -0.875, -1.   , -1.125, -0.125,
            -0.25 , -0.375, -0.5  , -0.625, -0.75 , -0.875, -1.   , -1.125 };

    for ( auto c = oc.getExternalForcesAtNodes( )->begin( ); c != oc.getExternalForcesAtNodes( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( extForceVectorsAnswer, c->second, 1e-6, 1e-1 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    for ( auto c = oc.getExternalCouplesAtNodes( )->begin( ); c != oc.getExternalCouplesAtNodes( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( floatVector( coupleVectorsAnswer.size( ), 0 ), c->second, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    for ( auto c = oc._test_cellLinearMomentumRHS.begin( ); c != oc._test_cellLinearMomentumRHS.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( forceVectorsAnswer, c->second, 1e-6, 1e-1 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    for ( auto c = oc._test_cellFirstMomentRHS.begin( ); c != oc._test_cellFirstMomentRHS.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( coupleVectorsAnswer, c->second, 1e-6, 1e-1 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    Eigen::MatrixXd LHSX( 288, 1 );
    LHSX <<  1.6082461 ,  0.23123014,  0.62393882,  1.32988565, -1.20046325,
            -1.49098297,  1.08575643, -0.27084579, -0.45887108,  1.13457348,
             1.14212648, -1.34876558, -1.53954667, -0.6699138 ,  0.4062938 ,
             1.51120934,  0.45950889, -0.3039844 ,  1.8313851 ,  1.41633787,
             1.0965811 ,  1.50251364, -1.68657903, -1.87216511,  0.82496983,
             0.21188063,  1.42106996,  1.81642989, -0.1000955 ,  0.19266961,
             0.93810141,  0.15452743,  0.98045664,  0.3140218 , -1.29539698,
             1.0298772 ,  1.79294532,  1.51096488,  1.42206134, -0.7942898 ,
            -1.56131436,  1.62426425,  1.67991981, -0.33085656, -1.8824174 ,
            -1.98883142, -1.86904329, -1.5365518 ,  1.39131847, -0.47607648,
             0.00974553, -0.15420091, -0.6692329 , -0.29326975, -1.78084752,
             1.97746862, -0.418282  , -1.04194253,  0.15101235,  1.55810889,
             0.29150197, -0.99929398, -0.4581576 ,  1.09085781, -0.59822029,
            -0.22436283, -0.34358714,  0.15518958,  1.67276323, -0.94694814,
             1.11237832,  0.39840522, -1.04803035,  0.15294796, -0.5688733 ,
            -0.3469194 ,  0.02140078, -1.85645887, -0.78465718,  1.49107402,
             1.9616645 , -0.24194266,  1.25961845, -0.87935036, -1.71921134,
             1.70558356,  0.75569485, -1.69431444,  0.7158976 ,  0.8212172 ,
            -1.45008094,  1.56941873,  1.78945147, -1.65800529,  0.34847407,
            -0.42676962, -0.19490982, -0.01828974,  1.7880325 ,  0.32964821,
            -1.07369484,  0.46494527, -1.86369121, -1.56866323,  0.00889209,
             0.16946288, -1.94731671, -1.81322178,  1.28646336,  0.85564197,
             0.28811254, -0.46973343,  0.14448512, -1.03384903,  0.15534826,
            -0.77913744,  1.22798127,  0.06452942,  0.09612534,  1.43803989,
            -0.57649306, -1.68445039, -0.46275924,  1.60444853,  1.23426519,
            -1.0681013 ,  0.60927561, -0.21281336, -1.07731193, -0.55479226,
            -0.6091404 , -0.23743334,  1.59429283, -0.82822957,  0.32881152,
             0.76887587, -0.80735223, -0.81406656,  1.37861004, -1.44708557,
            -1.52856592,  1.32006201,  0.69897149, -0.12453674, -1.4602061 ,
            -1.15926572, -1.82017397,  0.5993131 ,  1.68027963,  1.12213658,
            -1.48578834, -1.59165138, -1.99922335,  0.5415541 ,  1.43641856,
             1.26490651,  1.78234528, -1.94824744, -0.747312  ,  1.29833448,
            -1.20067926, -1.68068102,  1.36716021,  1.80866173, -0.0364494 ,
             1.06392003, -1.46910731, -0.17134657, -0.02810908,  1.24846583,
             0.59492076,  0.78519705, -1.95111884, -1.14141891, -1.62851376,
             0.83826821,  1.68252774,  1.0556338 ,  0.58519686, -0.02516275,
            -0.72941457, -1.32498254,  0.46014727,  0.48855993,  1.70462867,
            -0.68164314,  1.97719623,  0.1438778 ,  0.76477815, -0.3680267 ,
             1.2548148 , -0.03005103,  0.60236049, -1.97593119,  1.86036645,
             1.82280531,  1.33938005, -0.50673755, -1.78690982,  1.35691525,
            -1.43122857, -1.67233715, -0.52498148, -0.12109349, -0.98761515,
             0.68838949, -0.85784641,  0.41203733, -0.25579901, -0.76895987,
            -0.9219585 , -0.61476178, -1.71490687, -0.43732827,  0.29157371,
            -0.06878555,  0.66687259,  1.01339296,  0.25225556,  1.75900635,
            -0.48708225,  1.23372708, -1.13049026, -1.95865261, -0.33059471,
             0.20845165, -0.20076705,  0.96219912, -0.8753787 ,  0.65915043,
            -1.89939549,  1.02776798,  1.17042635,  0.65582458, -1.25221368,
             1.96969231,  0.82878036,  0.39866316,  0.91485086, -0.31437875,
            -0.27374567,  0.08063585,  1.01416116, -1.09062704, -0.24350887,
            -1.76250544, -0.91778241, -0.94918421, -0.49795888, -1.58898471,
            -0.92796655,  0.04034382, -0.56300233,  1.2696228 , -1.30373122,
            -1.0500919 , -1.60771873, -1.70433096,  1.07673513,  1.94832553,
             0.71579226, -0.47986902,  1.02920333, -0.67345196, -1.88666695,
            -1.29978151,  0.30628494, -0.58037746, -0.12743149,  0.96912632,
             1.26811874, -1.15144818, -1.98883533,  0.32852892,  1.83145568,
            -1.89243164,  1.80968967,  1.75802041, -0.16549997, -1.21975459,
             1.47819122,  0.19478798,  1.51252287,  1.23971307, -0.74991309,
            -0.13283132, -0.13583932, -0.24156527;

    Eigen::MatrixXd LHSA( 96, 1 );
    LHSA << -0.67148377,  0.08644475,  0.1876113 ,  0.05316974,  0.40034055,
             0.10525142,  0.2532617 ,  0.01081782, -0.08153879, -0.21072644,
            -0.27487082, -0.14290544,  0.22260118,  0.0548315 , -0.39934768,
             0.21540622,  0.1000854 , -0.25457557,  0.16625719, -0.04739766,
             0.20677403, -0.02848581, -0.33025154,  0.37873074, -0.01076609,
            -0.46845323, -0.02160498,  0.41184923,  0.28058437,  0.32835265,
            -0.01610606, -0.20314495, -0.13646886, -0.23956204, -0.16970513,
            -0.16122469, -0.13167802,  0.14602455,  0.15379931, -0.10467893,
            -0.26175818, -0.04269552,  0.25807113,  0.12635382,  0.03454933,
             0.05360817, -0.26872938, -0.20561066,  0.17163892, -0.20123063,
             0.06966427,  0.06336308,  0.00530582, -0.07900626,  0.04271629,
             0.10468538,  0.15505774, -0.33192724,  0.24736403,  0.13116234,
            -0.42203307,  0.21572573,  0.1599863 , -0.3424212 ,  0.14623527,
             0.16944365, -0.38643524, -0.23316837,  0.10107278, -0.11397124,
            -0.14235407,  0.22499911, -0.36225596,  0.17564057, -0.38667163,
             0.48963247,  0.32025427,  0.10998442, -0.19407242,  0.34018225,
            -0.40723539, -0.21322009, -0.14876766, -0.09710544, -0.02182683,
            -0.20298848,  0.03559792, -0.02912419, -0.08736181, -0.32624219,
             0.07755545,  0.20602909,  0.10546668,  0.22379181,  0.17361177,
            -0.19989885;

    std::cout << "testing LHS matrix\n";
    for ( auto LHS = oc._test_stressProjectionLHS.begin( ); LHS != oc._test_stressProjectionLHS.end( ); LHS++ ){

        if ( ( LHS->second * LHSX - LHSA ).norm( ) > 1e-6 * ( LHSA.norm( ) + 1 ) ){

            std::cout << LHS->second * LHSX << "\n";

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    remove( "macroscale_dof.xdmf" );
    remove( "macroscale_dof.h5" );

    remove( "microscale_dof.xdmf" );
    remove( "microscale_dof.h5" );

    results << "test_overlapCoupling_processIncrement & True\n";
    return 0;
}

int test_overlapCoupling_processLastIncrements( std::ofstream &results ){
    /*!
     * Test processing the last increments
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    remove( "macroscale_dof.xdmf" );
    remove( "macroscale_dof.h5" );

    remove( "microscale_dof.xdmf" );
    remove( "microscale_dof.h5" );

    std::string filename = "testConfig_averaged_l2_projection.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_processLastIncrements & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_processLastIncrements & False\n";
        return 1;
    }

    error = oc.processLastIncrements( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_processLastIncrements & False\n";
        return 1;
    }

    const floatVector *projectedGhostMacroDisplacement = oc.getProjectedGhostMacroDisplacement( );
    const floatVector *projectedGhostMicroDisplacement = oc.getProjectedGhostMicroDisplacement( );

    if ( projectedGhostMacroDisplacement->size( ) == 0 ){
        results << "test_overlapCoupling_processLastIncrements (test 1) & False\n";
        return 1;
    }

    if ( projectedGhostMicroDisplacement->size( ) == 0 ){
        results << "test_overlapCoupling_processLastIncrements (test 2) & False\n";
        return 1;
    }

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    remove( "macroscale_dof.xdmf" );
    remove( "macroscale_dof.h5" );

    remove( "microscale_dof.xdmf" );
    remove( "microscale_dof.h5" );

    results << "test_overlapCoupling_processLastIncrements & True\n";
    return 0;
}

int test_MADOutlierDetection( std::ofstream &results ){
    /*!
     * Test the computation of outliers using the maximum absolution deviation 
     * detection metric.
     * 
     * :param std::ofstream &results: The output file.
     */

    floatVector x = { 0.70154526, 0.00265005, 0.29766985, 0.0570927 , 0.12136678 };
    uIntVector outliers;

    errorOut error = overlapCoupling::MADOutlierDetection( x, outliers, 5 );

    if ( error ){

        error->print( );
        results << "test_MADOutlierDetection & False\n";
        return 1;

    }

    if ( outliers.size() > 0 ){
        results << "test_MADOutlierDetection (test 1) & False\n";
        return 1;
    }

    overlapCoupling::MADOutlierDetection( x, outliers, 4 );

    if ( !vectorTools::fuzzyEquals( outliers, { 0 } ) ){

        results << "test_MADOutlierDetection (test 1) & False\n";
        return 1;

    }

    results << "test_MADOutlierDetection & True\n";
    return 0;
}

int test_formMicromorphicElementMassMatrix( std::ofstream &results ){
    /*!
     * Test the formation of the mass matrix for a single micromorphic element.
     *
     * :param std::ofstream &results: The output file
     */

    std::unique_ptr< elib::Element > element;

    floatMatrix reference_nodes =
        {
            { 0, 0, 0 },
            { 1, 0, 0 },
            { 1, 1, 0 },
            { 0, 1, 0 },
            { 0, 0, 1 },
            { 1, 0, 1 },
            { 1, 1, 1 },
            { 0, 1, 1 }
        };

    floatVector degreeOfFreedomValues =
        {
            -0.04066559,  0.0390943 , -0.00232655, -0.03264781,  0.02111065,
             0.00079745, -0.00409562,  0.03073762,  0.01799572, -0.02048283,
             0.00914281, -0.04062647,  0.02985775, -0.02334215,  0.01175122,
            -0.02032996,  0.0489607 , -0.01085405, -0.00867754,  0.00475848,
             0.04859963,  0.00442464,  0.0211827 ,  0.04006121, -0.03019419,
            -0.04428046,  0.04353155, -0.03999713,  0.03122015,  0.01574864,
            -0.0025961 , -0.0232323 , -0.03535742, -0.00621033, -0.00485358,
             0.00917265, -0.03230599,  0.0228427 , -0.00198795,  0.04252169,
             0.04223892,  0.01351459,  0.01295129, -0.04424972,  0.0322117 ,
            -0.03558341, -0.01547058,  0.01382653, -0.0379305 ,  0.04802716,
            -0.03400389,  0.0348995 , -0.02256884, -0.04231619,  0.00654967,
             0.01205778, -0.01045064, -0.01303744, -0.00900963, -0.01305943,
            -0.04141735, -0.02224897,  0.03300644, -0.03487759, -0.04996436,
            -0.00305152, -0.00471929, -0.04222132, -0.03447746, -0.03298123,
            -0.02164259, -0.03886441,  0.02980088,  0.00453406,  0.03759968,
             0.00967625,  0.04231608, -0.001698  ,  0.01340712, -0.03988135,
            -0.04814852,  0.0196056 , -0.01305974, -0.0333702 , -0.03428086,
            -0.04671768,  0.02380142, -0.02290295,  0.02739323,  0.00481225,
             0.0424606 , -0.00402194, -0.00948801,  0.0395755 ,  0.02496865,
            -0.04143367
        };

    floatMatrix nodeDisplacement;
    nodeDisplacement.reserve( 8 );
    for ( unsigned int i = 0; i < 8; i++ ){

        nodeDisplacement.push_back( floatVector( degreeOfFreedomValues.begin( ) + 12 * i,
                                                 degreeOfFreedomValues.begin( ) + 12 * i + 3 ) );

    }

    auto qrule = elib::default_qrules.find( "Hex8" );
    element = elib::build_element_from_string( "Hex8", { 10, 7, 3, 9, 1, 8, 13, 4 }, reference_nodes, qrule->second );
    element->update_node_positions( nodeDisplacement );

    overlapCoupling::DOFMap nodeIDToIndex =
        {
            { 10, 0 },
            {  7, 2 },
            {  3, 4 },
            {  9, 3 },
            {  1, 6 },
            {  8, 1 },
            { 13, 7 },
            {  4, 5 }
        };

    floatVector density =
        { 693.53490713, 1765.4802207 ,   91.36052518,  666.64526727,
          51.16415254,  398.63874113,  702.24020488, 1190.92397094 };

    floatVector momentOfInertia =
        {
            0.44595488,  0.13676299, -0.2525482 ,  0.13676299, -0.19746052,
           -0.30931581, -0.2525482 , -0.30931581,  0.17223508, -0.31771288,
            0.3408996 ,  0.11782056,  0.3408996 ,  0.38466266,  0.1586499 ,
            0.11782056,  0.1586499 ,  0.19104204,  0.11479176,  0.17524618,
            0.04975924,  0.17524618, -0.17921051,  0.15571073,  0.04975924,
            0.15571073, -0.0283223 , -0.28271601,  0.18506962,  0.26725315,
            0.18506962,  0.17534291, -0.12381001,  0.26725315, -0.12381001,
            0.33227658,  0.32864468, -0.16755297, -0.25824399, -0.16755297,
            0.10339622, -0.40619467, -0.25824399, -0.40619467, -0.19037855,
            0.27080415,  0.35140531, -0.00926281,  0.35140531, -0.1955969 ,
           -0.0311474 , -0.00926281, -0.0311474 , -0.31288414,  0.41621005,
            0.10722768,  0.16218443,  0.10722768, -0.00765149,  0.06275192,
            0.16218443,  0.06275192, -0.31463761,  0.49814788, -0.14514796,
            0.24525217, -0.14514796, -0.46891292,  0.10017765,  0.24525217,
            0.10017765, -0.21952903
        };

    floatVector answerData;
    Eigen::MatrixXd answer;

    errorOut error = readMatrixFromFile( "mass_matrix_answer.csv", answerData, answer );

    if ( error ){

        error->print( );
        results << "test_formMicromorphicElementMassMatrix & False\n";
        return 1;

    }

    std::vector< DOFProjection::T > coefficients;

    error = overlapCoupling::formMicromorphicElementMassMatrix( element, degreeOfFreedomValues,
                                                                momentOfInertia, density, &nodeIDToIndex, coefficients );
    DOFProjection::SparseMatrix result( 8 * 12, 8 * 12 );
    result.setFromTriplets( coefficients.begin( ), coefficients.end( ) );

    if ( error ){

        error->print( );
        results << "test_formMicromorphicElementMassMatrix & False\n";
        return 1;

    }

    if ( !answer.isApprox( result.toDense( ), 1e-5 ) ){

        results << "test_formMicromorphicElementMassMatrix (test 1) & False\n";
        return 1;

    }

    results << "test_formMicromorphicElementMassMatrix & True\n";
    return 0;
}

int test_computeMicromorphicElementRequiredValues( std::ofstream &results ){
    /*!
     * Test the computation of the default required values from the micromorphic element
     *
     * :param std::ofstream &results: The output file
     */

    std::unique_ptr< elib::Element > element;

    floatMatrix reference_nodes =
        {
            { 0, 0, 0 },
            { 1, 0, 0 },
            { 1, 1, 0 },
            { 0, 1, 0 },
            { 0, 0, 1 },
            { 1, 0, 1 },
            { 1, 1, 1 },
            { 0, 1, 1 }
        };

    floatVector degreeOfFreedomValues =
        {
            -0.1       , -0.1       , -0.1       , -0.03264781,  0.02111065,
             0.00079745, -0.00409562,  0.03073762,  0.01799572, -0.02048283,
             0.00914281, -0.04062647,  0.1       , -0.1       , -0.1       ,
            -0.02032996,  0.0489607 , -0.01085405, -0.00867754,  0.00475848,
             0.04859963,  0.00442464,  0.0211827 ,  0.04006121,  0.1       ,
             0.1       , -0.1       , -0.03999713,  0.03122015,  0.01574864,
            -0.0025961 , -0.0232323 , -0.03535742, -0.00621033, -0.00485358,
             0.00917265, -0.1       ,  0.1       , -0.1       ,  0.04252169,
             0.04223892,  0.01351459,  0.01295129, -0.04424972,  0.0322117 ,
            -0.03558341, -0.01547058,  0.01382653, -0.1       , -0.1       ,
             0.1       ,  0.0348995 , -0.02256884, -0.04231619,  0.00654967,
             0.01205778, -0.01045064, -0.01303744, -0.00900963, -0.01305943,
             0.1       , -0.1       ,  0.1       , -0.03487759, -0.04996436,
            -0.00305152, -0.00471929, -0.04222132, -0.03447746, -0.03298123,
            -0.02164259, -0.03886441,  0.1       ,  0.1       ,  0.1       ,
             0.00967625,  0.04231608, -0.001698  ,  0.01340712, -0.03988135,
            -0.04814852,  0.0196056 , -0.01305974, -0.0333702 , -0.1       ,
             0.1       ,  0.1       , -0.02290295,  0.02739323,  0.00481225,
             0.0424606 , -0.00402194, -0.00948801,  0.0395755 ,  0.02496865,
            -0.04143367
        };

    floatMatrix nodeDisplacement;
    nodeDisplacement.reserve( 8 );
    for ( unsigned int i = 0; i < 8; i++ ){

        nodeDisplacement.push_back( floatVector( degreeOfFreedomValues.begin( ) + 12 * i,
                                                 degreeOfFreedomValues.begin( ) + 12 * i + 3 ) );

    }

    auto qrule = elib::default_qrules.find( "Hex8" );
    element = elib::build_element_from_string( "Hex8", { 10, 7, 3, 9, 1, 8, 13, 4 }, reference_nodes, qrule->second );
    element->update_node_positions( nodeDisplacement );

    floatMatrix reshapedDOFValues = vectorTools::inflate( degreeOfFreedomValues, 8, 12 );
    floatVector uQptResult, XiQptResult, shapeFunctionsResult, deformationGradientResult;
    floatMatrix gradShapeFunctionsResult;
    floatType JResult, JxwResult;

    floatType JAnswer = 1.728;
    floatType JxwReferenceAnswer = 0.125;
    floatType JxwCurrentAnswer = 0.216;
    floatType tmp = 0.0577350269;

    floatMatrix shapeFunctionAnswer
        {
            { 0.490563,   0.131446,   0.0352208,  0.131446,   0.131446,   0.0352208,  0.00943739, 0.0352208 }, 
            { 0.131446,   0.490563,   0.131446,   0.0352208,  0.0352208,  0.131446,   0.0352208,  0.00943739 },
            { 0.0352208,  0.131446,   0.490563,   0.131446,   0.00943739, 0.0352208,  0.131446,   0.0352208 },
            { 0.131446,   0.0352208,  0.131446,   0.490563,   0.0352208,  0.00943739, 0.0352208,  0.131446 },
            { 0.131446,   0.0352208,  0.00943739, 0.0352208,  0.490563,   0.131446,   0.0352208,  0.131446 },
            { 0.0352208,  0.131446,   0.0352208,  0.00943739, 0.131446,   0.490563,   0.131446,   0.0352208 },
            { 0.00943739, 0.0352208,  0.131446,   0.0352208,  0.0352208,  0.131446,   0.490563,   0.131446 },
            { 0.0352208,  0.00943739, 0.0352208,  0.131446,   0.131446,   0.0352208,  0.131446,   0.490563 },
        };

    floatMatrix gradShapeFunctionsReferenceAnswer
        {
            { -0.622008, -0.622008, -0.622008, 0.622008, -0.166667, -0.166667, 0.166667, 0.166667, -0.0446582, -0.166667, 0.622008, -0.166667, -0.166667, -0.166667, 0.622008, 0.166667, -0.0446582, 0.166667, 0.0446582, 0.0446582, 0.0446582, -0.0446582, 0.166667, 0.166667 },
            { -0.622008, -0.166667, -0.166667, 0.622008, -0.622008, -0.622008, 0.166667, 0.622008, -0.166667, -0.166667, 0.166667, -0.0446582, -0.166667, -0.0446582, 0.166667, 0.166667, -0.166667, 0.622008, 0.0446582, 0.166667, 0.166667, -0.0446582, 0.0446582, 0.0446582 },
            { -0.166667, -0.166667, -0.0446582, 0.166667, -0.622008, -0.166667, 0.622008, 0.622008, -0.622008, -0.622008, 0.166667, -0.166667, -0.0446582, -0.0446582, 0.0446582, 0.0446582, -0.166667, 0.166667, 0.166667, 0.166667, 0.622008, -0.166667, 0.0446582, 0.166667 },
            { -0.166667, -0.622008, -0.166667, 0.166667, -0.166667, -0.0446582, 0.622008, 0.166667, -0.166667, -0.622008, 0.622008, -0.622008, -0.0446582, -0.166667, 0.166667, 0.0446582, -0.0446582, 0.0446582, 0.166667, 0.0446582, 0.166667, -0.166667, 0.166667, 0.622008 },
            { -0.166667, -0.166667, -0.622008, 0.166667, -0.0446582, -0.166667, 0.0446582, 0.0446582, -0.0446582, -0.0446582, 0.166667, -0.166667, -0.622008, -0.622008, 0.622008, 0.622008, -0.166667, 0.166667, 0.166667, 0.166667, 0.0446582, -0.166667, 0.622008, 0.166667 },
            { -0.166667, -0.0446582, -0.166667, 0.166667, -0.166667, -0.622008, 0.0446582, 0.166667, -0.166667, -0.0446582, 0.0446582, -0.0446582, -0.622008, -0.166667, 0.166667, 0.622008, -0.622008, 0.622008, 0.166667, 0.622008, 0.166667, -0.166667, 0.166667, 0.0446582 },
            { -0.0446582, -0.0446582, -0.0446582, 0.0446582, -0.166667, -0.166667, 0.166667, 0.166667, -0.622008, -0.166667, 0.0446582, -0.166667, -0.166667, -0.166667, 0.0446582, 0.166667, -0.622008, 0.166667, 0.622008, 0.622008, 0.622008, -0.622008, 0.166667, 0.166667 },
            { -0.0446582, -0.166667, -0.166667, 0.0446582, -0.0446582, -0.0446582, 0.166667, 0.0446582, -0.166667, -0.166667, 0.166667, -0.622008, -0.166667, -0.622008, 0.166667, 0.166667, -0.166667, 0.0446582, 0.622008, 0.166667, 0.166667, -0.622008, 0.622008, 0.622008 }
        };

    floatMatrix gradShapeFunctionsCurrentAnswer
        {
            { -0.51834, -0.51834, -0.51834, 0.51834, -0.138889, -0.138889, 0.138889, 0.138889, -0.0372152, -0.138889, 0.51834, -0.138889, -0.138889, -0.138889, 0.51834, 0.138889, -0.0372152, 0.138889, 0.0372152, 0.0372152, 0.0372152, -0.0372152, 0.138889, 0.138889 },
            { -0.51834, -0.138889, -0.138889, 0.51834, -0.51834, -0.51834, 0.138889, 0.51834, -0.138889, -0.138889, 0.138889, -0.0372152, -0.138889, -0.0372152, 0.138889, 0.138889, -0.138889, 0.51834, 0.0372152, 0.138889, 0.138889, -0.0372152, 0.0372152, 0.0372152 },
            { -0.138889, -0.138889, -0.0372152, 0.138889, -0.51834, -0.138889, 0.51834, 0.51834, -0.51834, -0.51834, 0.138889, -0.138889, -0.0372152, -0.0372152, 0.0372152, 0.0372152, -0.138889, 0.138889, 0.138889, 0.138889, 0.51834, -0.138889, 0.0372152, 0.138889 },
            { -0.138889, -0.51834, -0.138889, 0.138889, -0.138889, -0.0372152, 0.51834, 0.138889, -0.138889, -0.51834, 0.51834, -0.51834, -0.0372152, -0.138889, 0.138889, 0.0372152, -0.0372152, 0.0372152, 0.138889, 0.0372152, 0.138889, -0.138889, 0.138889, 0.51834 },
            { -0.138889, -0.138889, -0.51834, 0.138889, -0.0372152, -0.138889, 0.0372152, 0.0372152, -0.0372152, -0.0372152, 0.138889, -0.138889, -0.51834, -0.51834, 0.51834, 0.51834, -0.138889, 0.138889, 0.138889, 0.138889, 0.0372152, -0.138889, 0.51834, 0.138889 },
            { -0.138889, -0.0372152, -0.138889, 0.138889, -0.138889, -0.51834, 0.0372152, 0.138889, -0.138889, -0.0372152, 0.0372152, -0.0372152, -0.51834, -0.138889, 0.138889, 0.51834, -0.51834, 0.51834, 0.138889, 0.51834, 0.138889, -0.138889, 0.138889, 0.0372152 },
            { -0.0372152, -0.0372152, -0.0372152, 0.0372152, -0.138889, -0.138889, 0.138889, 0.138889, -0.51834, -0.138889, 0.0372152, -0.138889, -0.138889, -0.138889, 0.0372152, 0.138889, -0.51834, 0.138889, 0.51834, 0.51834, 0.51834, -0.51834, 0.138889, 0.138889 },
            { -0.0372152, -0.138889, -0.138889, 0.0372152, -0.0372152, -0.0372152, 0.138889, 0.0372152, -0.138889, -0.138889, 0.138889, -0.51834, -0.138889, -0.51834, 0.138889, 0.138889, -0.138889, 0.0372152, 0.51834, 0.138889, 0.138889, -0.51834, 0.51834, 0.51834 },
        };

    floatVector deformationGradientAnswer = { 1.2, 0.0, 0.0,
                                              0.0, 1.2, 0.0,
                                              0.0, 0.0, 1.2 };

    floatMatrix uQptAnswer =
        {
            { -tmp, -tmp, -tmp },
            {  tmp, -tmp, -tmp },
            {  tmp,  tmp, -tmp },
            { -tmp,  tmp, -tmp },
            { -tmp, -tmp,  tmp },
            {  tmp, -tmp,  tmp },
            {  tmp,  tmp,  tmp },
            { -tmp,  tmp,  tmp },
        };

    floatMatrix XiQptAnswer =
        {
            {  9.88136118e-01,  2.00813257e-02, -4.22070732e-03,  7.77908101e-04,  1.00864936e+00,  1.48284518e-02, -1.56589841e-02,  3.87462985e-03,  9.82616852e-01 },
            {  9.78744967e-01,  2.67710907e-02, -4.57959127e-03,  -4.19704309e-03,  9.95194617e-01,  1.60081991e-02, -6.32184786e-03,  7.02386931e-03,  1.00887017e+00 },
            { 9.81712207e-01,  3.26011268e-02,  7.54296900e-03,   2.29737042e-03,  9.77737555e-01, -1.40649506e-02, -5.17738385e-03, -2.99266479e-03,  1.00281440e+00 },
            {  1.00882498e+00,  3.31486739e-02,  7.47592280e-03,  1.14078262e-02,  9.77539770e-01,  1.15950117e-02, -1.56866541e-02, -3.97893908e-03,  9.96610936e-01 },
            {  1.00597891e+00, -6.26624679e-03, -2.02399680e-02,  8.23385344e-03,  1.00086213e+00, -7.72359125e-03, -8.68680450e-03, -3.08540289e-03,  9.78507688e-01 },
            {  9.83113378e-01, -1.22727139e-02, -7.82937249e-03,  5.49537395e-04,  9.75961052e-01, -1.88694411e-02, -1.46165030e-02, -9.84911837e-03,  9.77660786e-01 },
            {  9.93597014e-01,  2.45120434e-02,  7.93844761e-05,  1.15392913e-02,  9.70627254e-01, -3.13985669e-02,  7.91831580e-03, -6.63726565e-03,  9.75335260e-01 },
            {  9.96234380e-01,  2.21313473e-02, -1.27548343e-03,  2.46714378e-02,  9.87375422e-01, -9.49015591e-03,  1.35403120e-02,  6.90293147e-03,  9.73290037e-01 }
        };

    errorOut error = NULL;
    for ( auto qpt = element->qrule.begin( ); qpt != element->qrule.end( ); qpt++ ){

        error = overlapCoupling::computeMicromorphicElementRequiredValues( element, qpt, 3, reshapedDOFValues, true,
                                                                           shapeFunctionsResult, gradShapeFunctionsResult,
                                                                           deformationGradientResult,
                                                                           JResult, JxwResult, uQptResult, XiQptResult );

        if ( error ){
        
            error->print( );
            results << "test_computeMicromorphicElementRequiredValues & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( shapeFunctionAnswer[ qpt - element->qrule.begin( ) ], shapeFunctionsResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 1) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( gradShapeFunctionsReferenceAnswer[ qpt - element->qrule.begin( ) ],
                                        vectorTools::appendVectors( gradShapeFunctionsResult ) ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 2) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( deformationGradientAnswer, deformationGradientResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 3) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( JResult, JAnswer ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 4) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( JxwResult, JxwReferenceAnswer ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 5) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( uQptAnswer[ qpt - element->qrule.begin( ) ], uQptResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 6) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( XiQptAnswer[ qpt - element->qrule.begin( ) ], XiQptResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 7) & False\n";
            return 1;

        }

        error = overlapCoupling::computeMicromorphicElementRequiredValues( element, qpt, 3, reshapedDOFValues, false,
                                                                           shapeFunctionsResult, gradShapeFunctionsResult,
                                                                           deformationGradientResult,
                                                                           JResult, JxwResult, uQptResult, XiQptResult );

        if ( error ){
        
            error->print( );
            results << "test_computeMicromorphicElementRequiredValues & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( shapeFunctionAnswer[ qpt - element->qrule.begin( ) ], shapeFunctionsResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 8) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( gradShapeFunctionsCurrentAnswer[ qpt - element->qrule.begin( ) ],
                                        vectorTools::appendVectors( gradShapeFunctionsResult ) ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 9) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( deformationGradientAnswer, deformationGradientResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 10) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( JResult, JAnswer ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 11) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( JxwResult, JxwCurrentAnswer ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 12) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( uQptAnswer[ qpt - element->qrule.begin( ) ], uQptResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 13) & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( XiQptAnswer[ qpt - element->qrule.begin( ) ], XiQptResult ) ){

            results << "test_computeMicromorphicElementRequiredValues (test 14) & False\n";
            return 1;

        }

    }

    results << "test_computeMicromorphicElementRequiredValues & True\n";
    return 0;

}

int test_computeMicromorphicElementInternalForceVector( std::ofstream &results ){
    /*!
     * Test the computation of the micromorphic internal force vector
     *
     * :param std::ofstream &results: The output file
     */

    std::unique_ptr< elib::Element > element;

    floatMatrix reference_nodes =
        {
            { 0, 0, 0 },
            { 1, 0, 0 },
            { 1, 1, 0 },
            { 0, 1, 0 },
            { 0, 0, 1 },
            { 1, 0, 1 },
            { 1, 1, 1 },
            { 0, 1, 1 }
        };

    floatVector degreeOfFreedomValues =
        {
            -0.04066559,  0.0390943 , -0.00232655, -0.03264781,  0.02111065,
             0.00079745, -0.00409562,  0.03073762,  0.01799572, -0.02048283,
             0.00914281, -0.04062647,  0.02985775, -0.02334215,  0.01175122,
            -0.02032996,  0.0489607 , -0.01085405, -0.00867754,  0.00475848,
             0.04859963,  0.00442464,  0.0211827 ,  0.04006121, -0.03019419,
            -0.04428046,  0.04353155, -0.03999713,  0.03122015,  0.01574864,
            -0.0025961 , -0.0232323 , -0.03535742, -0.00621033, -0.00485358,
             0.00917265, -0.03230599,  0.0228427 , -0.00198795,  0.04252169,
             0.04223892,  0.01351459,  0.01295129, -0.04424972,  0.0322117 ,
            -0.03558341, -0.01547058,  0.01382653, -0.0379305 ,  0.04802716,
            -0.03400389,  0.0348995 , -0.02256884, -0.04231619,  0.00654967,
             0.01205778, -0.01045064, -0.01303744, -0.00900963, -0.01305943,
            -0.04141735, -0.02224897,  0.03300644, -0.03487759, -0.04996436,
            -0.00305152, -0.00471929, -0.04222132, -0.03447746, -0.03298123,
            -0.02164259, -0.03886441,  0.02980088,  0.00453406,  0.03759968,
             0.00967625,  0.04231608, -0.001698  ,  0.01340712, -0.03988135,
            -0.04814852,  0.0196056 , -0.01305974, -0.0333702 , -0.03428086,
            -0.04671768,  0.02380142, -0.02290295,  0.02739323,  0.00481225,
             0.0424606 , -0.00402194, -0.00948801,  0.0395755 ,  0.02496865,
            -0.04143367
        };

    floatMatrix nodeDisplacement;
    nodeDisplacement.reserve( 8 );
    for ( unsigned int i = 0; i < 8; i++ ){

        nodeDisplacement.push_back( floatVector( degreeOfFreedomValues.begin( ) + 12 * i,
                                                 degreeOfFreedomValues.begin( ) + 12 * i + 3 ) );

    }

    auto qrule = elib::default_qrules.find( "Hex8" );
    element = elib::build_element_from_string( "Hex8", { 10, 7, 3, 9, 1, 8, 13, 4 }, reference_nodes, qrule->second );
    element->update_node_positions( nodeDisplacement );

    overlapCoupling::DOFMap nodeIDToIndex =
        {
            { 10, 0 },
            {  7, 2 },
            {  3, 4 },
            {  9, 3 },
            {  1, 6 },
            {  8, 1 },
            { 13, 7 },
            {  4, 5 }
        };

    floatVector cauchyStress =
        {
            -0.45969764,  1.14822033, -1.36295921, -0.5000321 , -1.42325377,
            -0.90204189, -0.44663969,  0.28638596,  1.33786465,  1.16783921,
            -0.98509841, -0.34726126, -0.59238156, -1.46301493, -0.53019081,
            -0.39287991,  0.51670525, -0.95313864,  0.81999322,  0.5800708 ,
             1.2254044 ,  1.15398799, -0.20563341,  0.12685255,  1.15753279,
             0.57312611,  0.34437528,  1.01836067, -0.59267923,  0.03311643,
            -0.53763152, -0.0644638 ,  0.89123469,  0.99872816, -0.08496691,
            -0.84814803,  0.39031588,  0.44741445,  0.72010183,  0.52602873,
            -0.33027138,  1.15610447, -0.41258865,  0.99682975,  0.98593404,
             0.40067375, -1.3895129 , -1.40599158,  0.31018865, -1.30229246,
            -0.15995105, -1.13487559, -0.57113566,  0.81241244,  1.15210152,
             1.13622788,  0.95633037, -0.58679069,  0.90543097, -1.19180477,
            -1.09231537, -0.93865242,  1.24863015, -0.2600532 ,  1.11112689,
             0.90341929, -0.02295448,  0.42265386, -0.1380203 , -0.04614231,
            -0.94391283,  1.49444364
        };

    floatVector symmetricMicroStress =
        {
            0.10117987,  1.37582714, -0.97206764, -1.40673876,  0.47691406,
            0.09568587,  0.12415184, -0.6585871 , -0.63620366, -1.00348131,
           -0.05910748, -1.22504474,  0.95867523,  0.02521843,  0.62481998,
           -0.33498708,  0.72102278, -0.39559034, -1.37320129, -1.32955232,
            0.5171092 ,  1.13703248, -1.35457667,  0.53214959, -0.02469652,
            0.70267216,  0.03722338, -1.31772695, -0.21335467,  0.81220832,
           -0.98121639, -0.887004  , -0.46146348,  1.48095631,  0.90664032,
            0.9653973 , -0.68733076, -0.46874014,  0.59222427,  0.26895499,
           -1.4239325 , -0.74967125, -0.64881934, -1.10759892, -1.33468558,
           -0.34723361, -0.74556723, -0.01161972, -0.93394231,  0.40266023,
            0.39595606,  0.10547423,  1.43015732,  1.07281307, -0.01329458,
            1.27192279,  0.03659976,  1.46190328,  0.55927909, -0.78215141,
           -0.51928948, -1.18045737,  0.03809132, -0.60385602, -0.00554656,
            0.19209442,  1.02908222,  0.82778923, -1.38683747, -0.25118759,
            0.28963075, -0.85042722
        };

    floatVector higherOrderStress =
        {
            0.39229651, -1.27188227, -0.75057189, -1.36700197, -1.15241392,
           -0.83755545, -0.21498426,  1.02131234, -1.24288282, -0.76007716,
           -0.92450045, -1.04121465,  1.14264551, -0.94492367,  1.02005295,
            0.13905738, -1.08039572, -1.27358044,  0.7775386 , -1.49857357,
            0.68073573, -0.9120657 , -1.196205  , -0.27346835, -0.86683224,
            0.11928698, -1.40885654, -0.25191655, -1.21035354, -0.13516862,
            0.42551714, -1.367081  ,  1.41603026, -0.29139481,  0.66615639,
            1.45493764, -1.00795599, -1.08847277, -0.40932263,  1.37012762,
            1.3090307 , -1.0598115 ,  1.42470902, -1.10917598, -0.15906364,
            1.17699552, -0.74454968,  0.37418966,  0.5419793 ,  0.52773656,
            0.53722857,  0.60228796,  0.15701482, -1.40756007, -0.84260176,
            0.46713455, -0.9109566 , -0.4587495 , -0.27443065,  1.0929497 ,
           -0.56260743, -0.37428778,  1.1547629 , -1.26992331,  0.02180244,
           -0.99848932, -1.10047995, -0.73077205, -1.48752174,  0.36737644,
           -0.58238242,  1.01795009, -0.66254986, -1.10319654,  0.66144406,
            0.0965511 , -0.49245866, -0.35159224, -0.35986417,  0.52077983,
            0.37951793,  1.02467304,  1.12864357, -0.07594849,  1.31868883,
            0.58521494, -0.03745289, -0.31158896, -0.17807641, -0.91852362,
            0.03084061,  0.92307018,  0.0202181 ,  1.29077688,  1.30206557,
            0.0921027 ,  1.27553382,  0.34262177, -1.48643965,  0.61201724,
            1.49347647,  0.64445638,  0.28889739, -1.43068598,  0.03570288,
           -0.08223876, -0.44111068, -1.09584718, -0.19206467,  0.46059502,
           -0.40809829,  1.08533245,  0.51658791,  0.7905599 , -0.96952095,
            1.23554638,  0.87045523, -0.84229178, -1.44913607, -0.18983496,
            0.10093825, -1.1185102 ,  0.44234559,  0.60141309, -0.89067448,
            0.31832896, -0.55737341,  1.12297842, -1.23493082, -1.41979439,
            0.18887879, -0.70628017,  0.14308575, -0.39197145, -0.5690881 ,
            0.57602306,  0.99915749,  0.86749036, -0.43421454,  1.43900802,
            0.29587861,  1.48878773, -0.08768882, -0.14582274,  1.1333723 ,
            0.24309208,  0.23585109,  0.81010256, -0.27329013, -0.23483857,
           -0.25537111,  1.31371579, -0.36190462, -0.67439408,  1.15628843,
           -1.09453436, -0.16136072, -0.76631781, -1.3900636 ,  0.51565904,
            0.51245367,  1.41764668, -0.58441173,  1.34112171,  1.31309759,
           -0.39597866, -0.14284022,  1.41928088, -0.10742838,  0.9797636 ,
            1.0535908 , -0.63374945, -1.13702384,  0.82982881, -0.50849777,
           -0.00453001,  1.19397083,  0.44670022,  1.28185961, -0.45946417,
            1.09819494,  0.75844966,  0.30450907, -1.46525729, -1.07253823,
            1.31673755,  0.56305741,  0.7642357 , -0.42590002,  0.30657098,
           -1.2917611 ,  0.94766539,  0.89394059, -0.30803095, -0.16673604,
           -1.22619946, -0.46748525, -1.229485  , -0.73181422,  1.213918  ,
           -0.17310261,  0.28313866, -1.43347574,  0.99141871, -0.79494399,
           -0.0173043 ,  0.80792229, -0.79151497,  1.02627251, -1.09927188,
            0.0545878 ,  0.99842225,  0.3788543 , -0.01641077, -0.60538728,
            0.91039346
        };


    floatVector answer =
        {
            0.02739017,  0.08038737, -0.08497879,  0.00601456,  0.08946546,
            0.27811672,  0.01057751,  0.38096787, -0.13814203, -0.04606262,
            0.00897415,  0.18528464, -0.03866123,  0.07169987,  0.13703867,
           -0.00299104,  0.18306985,  0.04642505, -0.03120374,  0.18365814,
            0.08300798,  0.0682409 ,  0.07193369, -0.0009522 ,  0.2231008 ,
            0.16225866, -0.13121907, -0.03476567,  0.01903948,  0.12197495,
           -0.03377065,  0.0310203 ,  0.33407419, -0.33467015,  0.26452521,
            0.06854623, -0.27872297, -0.1310955 , -0.22403472, -0.34850733,
           -0.29972842, -0.022039  ,  0.18845184,  0.06034523, -0.01170348,
            0.2754791 , -0.12785453, -0.138628  ,  0.22756336, -0.09145956,
           -0.11057031, -0.47555049, -0.00796385, -0.26483904,  0.15577536,
            0.10029588, -0.04724648, -0.0130081 , -0.15922499,  0.03647805,
           -0.0681886 , -0.25736625,  0.03305347, -0.23813408,  0.16353921,
           -0.21829181, -0.1012466 , -0.2838822 ,  0.09664438,  0.08055681,
           -0.25891319, -0.2240477 , -0.19875004,  0.1424041 ,  0.19372511,
            0.05390138, -0.03151748, -0.05120333, -0.36223144,  0.00697449,
           -0.431687  , -0.0595205 , -0.2103422 , -0.43627243,  0.10626851,
            0.0231713 ,  0.18698564, -0.15977423,  0.12031186,  0.28120184,
           -0.01767001, -0.20123309,  0.29285515, -0.07107898,  0.31030247,
           -0.16779977
        };

    Eigen::MatrixXd internalForceVectorAnswer = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned>( answer.data(), 12 * 8, 1 );
    Eigen::MatrixXd internalForceVectorResult = Eigen::MatrixXd::Zero( 12 * 8, 1 );

    errorOut error = overlapCoupling::formMicromorphicElementInternalForceVector( element, degreeOfFreedomValues,
                                                                                  cauchyStress, symmetricMicroStress, higherOrderStress,
                                                                                  &nodeIDToIndex, internalForceVectorResult );

    if ( error ){

        error->print( );
        results << "test_computeMicromorphicElementInternalForceVector & False\n";
        return 1;

    }

    if ( !internalForceVectorAnswer.isApprox( internalForceVectorResult, 1e-5 ) ){

        results << "test_computeMicromorphicElementInternalForceVector (test 1) & False\n";
        return 1;

    }

    results << "test_computeMicromorphicElementInternalForceVector & True\n";
    return 0;
}

int test_readWriteSparseMatrixToXDMF( std::ofstream &results ){
    /*!
     * Test reading and writing a sparse matrix to XDMF file
     *
     * :param std::ofstream &results: The output file
     */

    shared_ptr< XdmfDomain > domain = XdmfDomain::New( );
    shared_ptr< XdmfUnstructuredGrid > grid = XdmfUnstructuredGrid::New( );

    std::string filename = "test_output_file";
    std::string h5_filename = filename + ".h5";
    std::string xdmf_filename = filename + ".xdmf";

    remove( h5_filename.c_str( ) );
    remove( xdmf_filename.c_str( ) );

    SparseMatrix A1( 3, 4 );
    std::vector< T > triplets;
    triplets.reserve( 7 );
    triplets.push_back( T( 0, 0, 1.0 ) );
    triplets.push_back( T( 0, 3, 1.5 ) );
    triplets.push_back( T( 2, 1, 7.0 ) );
    triplets.push_back( T( 1, 2, 5.0 ) );
    triplets.push_back( T( 0, 1, 2.0 ) );
    triplets.push_back( T( 1, 0, 1.6 ) );
    triplets.push_back( T( 0, 2, 3.0 ) );

    A1.setFromTriplets( triplets.begin( ), triplets.end( ) );

    domain->insert( grid );

    std::string matrixName = "A_MATRIX";
    errorOut error = overlapCoupling::writeSparseMatrixToXDMF( A1, matrixName, filename, domain, grid );

    if ( error ){

        error->print( );
        results << "test_readWriteSparseMatrixToXDMF & False\n";
        return 1;

    }

    shared_ptr< XdmfReader > reader = XdmfReader::New( );
    shared_ptr< XdmfDomain > _readDomain = shared_dynamic_cast< XdmfDomain >( reader->read( "test_output_file.xdmf" ) );
    shared_ptr< XdmfUnstructuredGrid > _readGrid = _readDomain->getUnstructuredGrid( 0 );

    SparseMatrix A1_result;
    error = overlapCoupling::readSparseMatrixFromXDMF( _readGrid, matrixName, A1_result );
    
    if ( error ){

        error->print( );
        results << "test_readWriteSparseMatrixToXDMF & False\n";
        return 1;

    }

    if ( !A1.isApprox( A1_result ) ){

        results << "test_readWriteSparseMatrixToXDMF (test 1) & False\n";
        return 1;

    }

    remove( h5_filename.c_str( ) );
    remove( xdmf_filename.c_str( ) );

    results << "test_readWriteSparseMatrixToXDMF & True\n";
    return 0;

}

int test_readWriteDenseMatrixToXDMF( std::ofstream &results ){
    /*!
     * Test reading and writing a dense matrix to XDMF file
     *
     * :param std::ofstream &results: The output file
     */

    shared_ptr< XdmfDomain > domain = XdmfDomain::New( );
    shared_ptr< XdmfUnstructuredGrid > grid = XdmfUnstructuredGrid::New( );

    std::string filename = "test_output_file";
    std::string h5_filename = filename + ".h5";
    std::string xdmf_filename = filename + ".xdmf";

    remove( h5_filename.c_str( ) );
    remove( xdmf_filename.c_str( ) );

    Eigen::MatrixXd A( 3, 4 );
    A << 1,  2,  3,  4,
         5,  6,  7,  8,
         9, 10, 11, 12;

    domain->insert( grid );

    std::string matrixName = "A_MATRIX";
    errorOut error = overlapCoupling::writeDenseMatrixToXDMF( A, matrixName, filename, domain, grid );

    if ( error ){

        error->print( );
        results << "test_readWriteDenseMatrixToXDMF & False\n";
        return 1;

    }

    shared_ptr< XdmfReader > reader = XdmfReader::New( );
    shared_ptr< XdmfDomain > _readDomain = shared_dynamic_cast< XdmfDomain >( reader->read( "test_output_file.xdmf" ) );
    shared_ptr< XdmfUnstructuredGrid > _readGrid = _readDomain->getUnstructuredGrid( 0 );

    Eigen::MatrixXd result;
    error = overlapCoupling::readDenseMatrixFromXDMF( _readGrid, matrixName, result );

    if ( error ){

        error->print( );
        results << "test_readWriteDenseMatrixToXDMF & False\n";
        return 1;

    }

    if ( !A.isApprox( result ) ){

        results << "test_readWriteDenseMatrixToXDMF (test 1) & False\n";
        return 1;

    }

    remove( h5_filename.c_str( ) );
    remove( xdmf_filename.c_str( ) );

    domain = XdmfDomain::New( );
    grid = XdmfUnstructuredGrid::New( );
    A = Eigen::MatrixXd( 100, 200 );

    uIntType index = 0;
    for ( uIntType i = 0; i < 100; i++ ){
        for ( uIntType j = 0; j < 200; j++, index++ ){
            A( i, j ) = index;
        }
    }

    domain->insert( grid );

    matrixName = "A_MATRIX";
    error = overlapCoupling::writeDenseMatrixToXDMF( A, matrixName, filename, domain, grid );

    if ( error ){

        error->print( );
        results << "test_readWriteDenseMatrixToXDMF & False\n";
        return 1;

    }

    reader = XdmfReader::New( );
    _readDomain = shared_dynamic_cast< XdmfDomain >( reader->read( "test_output_file.xdmf" ) );
    _readGrid = _readDomain->getUnstructuredGrid( 0 );

    error = overlapCoupling::readDenseMatrixFromXDMF( _readGrid, matrixName, result );

    if ( error ){

        error->print( );
        results << "test_readWriteDenseMatrixToXDMF & False\n";
        return 1;

    }

    if ( !A.isApprox( result ) ){

        results << "test_readWriteDenseMatrixToXDMF (test 2) & False\n";
        return 1;

    }
    
    remove( h5_filename.c_str( ) );
    remove( xdmf_filename.c_str( ) );

    results << "test_readWriteDenseMatrixToXDMF & True\n";
    return 0;

}


int main(){
    /*!
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or False 
    if the test passes or fails respectively.
    */

    //Run the setup
    std::unique_ptr< errorNode > error;
    error.reset( _createXDMFDatafiles( ) );
    if ( error ){

        error->print( );
        return 1;

    }

    //Open the results file
    std::ofstream results;
    results.open("results.tex");

    test_overlapCoupling_constructor( results );
    test_overlapCoupling_initializeCoupling_l2_projection( results );
    test_overlapCoupling_initializeCoupling_averaged_l2_projection( results );
    test_overlapCoupling_processIncrement( results );
//    test_overlapCoupling_processLastIncrements( results );
//    test_overlapCoupling_getReferenceFreeMicroDomainMasses( results );
//    test_overlapCoupling_getReferenceGhostMicroDomainMasses( results );
//    test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass( results );
//    test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass( results );
////    test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions( results );
////    test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions( results );
//    test_MADOutlierDetection( results );
//    test_formMicromorphicElementMassMatrix( results );
//    test_computeMicromorphicElementRequiredValues( results );
//    test_computeMicromorphicElementInternalForceVector( results );
    test_readWriteSparseMatrixToXDMF( results );
    test_readWriteDenseMatrixToXDMF( results );

    //Close the results file
    results.close();
}
