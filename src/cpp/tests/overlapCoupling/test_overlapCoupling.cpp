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

    fg = fileGenerator::fileGenerator( "macroscale_Arlequin.yaml" );
    if ( fg.build( ) ){

        fg.getError( )->print( );
        errorOut result = new errorNode( "_createXDMFDatafiles", "Error in creation of the macroscale Arlequin datafile" );
        return result;

    }

    fg = fileGenerator::fileGenerator( "macroscale_realistic.yaml" );
    if ( fg.build( ) ){

        fg.getError( )->print( );
        errorOut result = new errorNode( "_createXDMFDatafiles", "Error in creation of the macroscale realistic datafile" );
        return result;

    }

    fg = fileGenerator::fileGenerator( "microscale_realistic.yaml" );
    if ( fg.build( ) ){

        fg.getError( )->print( );
        errorOut result = new errorNode( "_createXDMFDatafiles", "Error in creation of the microscale realistic datafile" );
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

    std::string testName = "overlapCoupling_initializeCoupling_l2_projection";
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

    std::string testName = "overlapCoupling_initializeCoupling_averaged_l2_projection";
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
    SparseMatrix BDhatQ;
    overlapCoupling::readSparseMatrixFromXDMF( _readGrid, "BDhatQ", BDhatQ );

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

    SparseMatrix BDhatD;
    overlapCoupling::readSparseMatrixFromXDMF( _readGrid, "BDhatD", BDhatD );

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

    SparseMatrix BQhatQ;
    overlapCoupling::readSparseMatrixFromXDMF( _readGrid, "BQhatQ", BQhatQ );

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

    SparseMatrix BQhatD;
    overlapCoupling::readSparseMatrixFromXDMF( _readGrid, "BQhatD", BQhatD );

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

int test_overlapCoupling_initializeCoupling_Arlequin( std::ofstream &results ){
    /*!
     * Test the initialization of the coupling for the Arlequin method
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string filename = "testConfig_Arlequin.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_initializeCoupling_Arlequin & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_initializeCoupling_Arlequin & False\n";
        return 1;
    }

    if ( !std::ifstream( "reference_information.xdmf" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_Arlequin (test 1) & False\n";
        return 1;

    }

    if ( !std::ifstream( "reference_information.h5" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_Arlequin (test 2) & False\n";
        return 1;

    }

    if ( !std::ifstream( "homogenized_response.xdmf" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_Arlequin (test 3) & False\n";
        return 1;

    }

    if ( !std::ifstream( "homogenized_response.h5" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_Arlequin (test 4) & False\n";
        return 1;

    }

    std::string testName = "overlapCoupling_initializeCoupling_Arlequin";
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

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_overlapCoupling_initializeCoupling_Arlequin & True\n";
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

    remove( "macroscale_out.xdmf" );
    remove( "macroscale_out.h5" );

    remove( "microscale_out.xdmf" );
    remove( "microscale_out.h5" );

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
        = { -0.125, -0.5  , -0.875, -0.25 , -0.625, -1.   , -0.375, -0.75 ,
            -1.125, -0.125, -0.5  , -0.875, -0.25 , -0.625, -1.   , -0.375,
            -0.75 , -1.125, -0.125, -0.5  , -0.875, -0.25 , -0.625, -1.   ,
            -0.375, -0.75 , -1.125, -0.125, -0.5  , -0.875, -0.25 , -0.625,
            -1.   , -0.375, -0.75 , -1.125, -0.125, -0.5  , -0.875, -0.25 ,
            -0.625, -1.   , -0.375, -0.75 , -1.125, -0.125, -0.5  , -0.875,
            -0.25 , -0.625, -1.   , -0.375, -0.75 , -1.125, -0.125, -0.5  ,
            -0.875, -0.25 , -0.625, -1.   , -0.375, -0.75 , -1.125, -0.125,
            -0.5  , -0.875, -0.25 , -0.625, -1.   , -0.375, -0.75 , -1.125 };

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

    for ( auto LHS = oc._test_stressProjectionLHS.begin( ); LHS != oc._test_stressProjectionLHS.end( ); LHS++ ){

        if ( ( LHS->second * LHSX - LHSA ).norm( ) > 1e-6 * ( LHSA.norm( ) + 1 ) ){

            std::cout << LHS->second * LHSX << "\n";

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    floatVector cauchyQptAnswer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    floatVector cauchyAnswer = vectorTools::appendVectors( floatMatrix( 8, cauchyQptAnswer ) );
    for ( auto c = oc.getQuadraturePointCauchyStress( )->begin( ); c != oc.getQuadraturePointCauchyStress( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, cauchyAnswer, 1e-6, 1e-1 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    floatVector higherOrderAnswer( 27 * 8, 0 );
    for ( auto c = oc.getQuadraturePointHigherOrderStress( )->begin( ); c != oc.getQuadraturePointHigherOrderStress( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, higherOrderAnswer, 1e-6, 1e-1 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    floatVector densityAnswer( 8, 2 );
    for ( auto c = oc.getQuadraturePointDensities( )->begin( ); c != oc.getQuadraturePointDensities( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, densityAnswer, 1e-6, 1e-2 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    floatVector bodyForceQptAnswer = { -1, -2, -3 };
    floatVector bodyForceAnswer = vectorTools::appendVectors( floatMatrix( 8, bodyForceQptAnswer ) );
    for ( auto c = oc.getQuadraturePointBodyForce( )->begin( ); c != oc.getQuadraturePointBodyForce( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, bodyForceAnswer, 1e-6, 1e-2 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    floatVector accelerationQptAnswer = { 0., 0., 3e-3 };
    floatVector accelerationAnswer = vectorTools::appendVectors( floatMatrix( 8, accelerationQptAnswer ) );
    for ( auto c = oc.getQuadraturePointAccelerations( )->begin( ); c != oc.getQuadraturePointAccelerations( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, accelerationAnswer, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    floatVector microInertiaQptAnswer = { 0.0625, 0., 0., 0., 0.0625, 0., 0., 0., 0.0625 };
    floatVector microInertiasAnswer = vectorTools::appendVectors( floatMatrix( 8, microInertiaQptAnswer ) );
    for ( auto c = oc.getQuadraturePointMicroInertias( )->begin( ); c != oc.getQuadraturePointMicroInertias( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, microInertiasAnswer, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    floatVector bodyCouplesAnswer( 9 * 8, 0 );
    for ( auto c = oc.getQuadraturePointBodyCouples( )->begin( ); c != oc.getQuadraturePointBodyCouples( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, bodyCouplesAnswer, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    floatVector microSpinInertiasAnswer( 9 * 8, 0 );
    for ( auto c = oc.getQuadraturePointMicroSpinInertias( )->begin( ); c != oc.getQuadraturePointMicroSpinInertias( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, microSpinInertiasAnswer, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    floatVector symmetricMicroStressQptAnswer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    floatVector symmetricMicroStressAnswer = vectorTools::appendVectors( floatMatrix( 8, symmetricMicroStressQptAnswer ) );
    for ( auto c = oc.getQuadraturePointSymmetricMicroStress( )->begin( ); c != oc.getQuadraturePointSymmetricMicroStress( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, symmetricMicroStressAnswer, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }

    Eigen::MatrixXd homogenizedFEXTAnswer( 144, 1 );
    homogenizedFEXTAnswer << -3.25, -4.25, -5.25,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,
                              0.  ,  0.  ,  0.  , -2.75, -3.25, -3.75,  0.  ,  0.  ,  0.  ,
                              0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  , -0.75, -0.75, -0.75,
                              0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,
                             -1.25, -1.75, -2.25,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,
                              0.  ,  0.  ,  0.  , -3.  , -4.5 , -6.  ,  0.  ,  0.  ,  0.  ,
                              0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  , -2.  , -2.5 , -3.  ,
                              0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,
                              2.  ,  2.5 ,  3.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,
                              0.  ,  0.  ,  0.  ,  1.  ,  0.5 ,  0.  ,  0.  ,  0.  ,  0.  ,
                              0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.25, -0.25, -0.75,
                              0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,
                              0.75,  0.75,  0.75,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,
                              0.  ,  0.  ,  0.  ,  2.75,  3.25,  3.75,  0.  ,  0.  ,  0.  ,
                              0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  2.25,  2.25,  2.25,
                              0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ;

    if ( ( homogenizedFEXTAnswer - *oc.getHomogenizedFEXT( ) ).norm( ) > 1e-2 * ( homogenizedFEXTAnswer.norm( ) + 1 ) ){

        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;

    }
    testNum++;

    Eigen::MatrixXd homogenizedFINTAnswer( 144, 1 );
    homogenizedFINTAnswer << -3.00, -3.75, -4.50,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                             -2.50, -2.75, -3.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                             -0.50, -0.25,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                             -1.00, -1.25, -1.50,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                             -2.50, -3.50, -4.50,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                             -1.50, -1.50, -1.50,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              2.50,  3.50,  4.50,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              1.50,  1.50,  1.50,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              0.50,  0.25,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              1.00,  1.25,  1.50,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              3.00,  3.75,  4.50,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              2.50,  2.75,  3.00,  0.00,
                              0.00,  0.00,  0.00,  0.00,
                              0.00,  0.00,  0.00,  0.00;

    if ( ( homogenizedFINTAnswer - *oc.getHomogenizedFINT( ) ).norm( ) > 1e-2 * ( homogenizedFINTAnswer.norm( ) + 1 ) ){

        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;

    }
    testNum++;

    Eigen::MatrixXd MX( 144, 1 );
    MX <<  1.03221787, -0.92246333,  1.27200681, -0.48761604,  0.37037086,
           1.02800141, -0.06719396,  0.25023267, -0.50356931,  1.98702766,
          -1.33517121, -0.98618331,  1.67520885, -0.32002969,  0.31204602,
          -0.62283931,  1.88485606,  0.14203533,  1.8679691 ,  1.95733768,
           0.48543641,  1.22018616, -0.31776795,  0.21502873,  0.64847517,
          -0.44455271, -1.59968356, -1.38163596, -0.46210322,  0.28330518,
          -1.53446563,  0.52655343, -0.40943491, -1.28177123,  0.31722293,
          -1.20203397,  1.8838875 , -0.41063676,  0.98500843, -1.91833504,
           1.65884954,  0.65857992,  0.75360733, -1.04602099, -1.92897066,
          -0.23344673,  1.47255343,  1.69531131, -1.27639949, -1.16911678,
           1.59589692,  0.30165247,  0.56896295,  1.63962054, -0.89920825,
          -1.63118444,  0.02242142,  0.72002951, -1.01593182, -1.68371068,
          -1.54483822,  1.87738139,  1.61092459, -1.03754536,  1.10254212,
          -1.24979886,  0.42068153,  1.80115139,  0.58122646,  0.69030709,
          -0.30644823,  1.61575065,  1.76422882,  0.49346777,  1.9656111 ,
           1.09841017, -1.43482989, -1.24436724,  0.40463851, -1.94135932,
          -0.93089867, -0.83737637,  0.19009932, -0.73989132, -0.41730156,
           0.56866435, -0.24432187,  0.74592874, -0.29527219,  0.79572626,
          -1.16371926,  1.48413216,  0.16221102,  1.33208692, -0.57881395,
           0.75027824,  0.37030067, -0.27326494, -1.1946677 , -0.54679603,
           1.99430325, -1.57736132, -0.04530493, -1.92386778, -1.63076862,
           0.54026597,  1.14049803, -1.01269005, -0.12859646,  0.95372803,
           0.04832678, -0.76899178, -0.80743886, -0.33441147, -0.89729875,
          -1.57117314,  0.08745574, -1.80049292,  0.92430717,  0.85418227,
           0.91247804,  0.46793921,  1.71737995, -1.50876716,  0.05699081,
          -1.66527938,  0.93377481, -1.53093827,  0.81617988, -1.67771678,
           0.27021483,  1.96160743, -1.88528753,  1.40697224,  0.29062036,
           1.80679358, -0.47801091,  1.32833655,  0.21266734, -0.58380344,
           0.14261189,  0.01538911, -0.66847627,  0.69894409;

    Eigen::MatrixXd MA( 144, 1 );
    MA <<  1.53013240e-01, -9.70591273e-02,  2.15252865e-01, -8.74254100e-03,
           1.08039031e-02,  9.49019236e-03,  1.27404678e-03,  2.78050058e-03,
          -5.77308992e-03,  1.35226657e-02, -6.40736937e-03, -3.12201996e-03,
           1.69177983e-01, -1.96544754e-02,  1.52574728e-01, -9.78039351e-03,
           1.18129517e-02,  2.48015410e-03,  5.54062296e-03,  1.05425081e-02,
          -1.71082442e-03,  9.24435367e-03, -4.12332277e-03, -7.38391094e-04,
           2.16156099e-01, -2.43256333e-02,  6.59843492e-02, -1.04635169e-02,
           4.43456111e-03,  1.89751832e-03, -1.55702366e-03,  3.20144958e-03,
          -7.10130397e-03, -5.31160114e-04,  1.42391437e-03, -2.23292618e-03,
           2.12091704e-01, -6.10414527e-02,  1.38426904e-01, -1.11823334e-02,
           8.60144588e-03,  7.82497060e-03, -1.07938453e-03, -4.36201176e-04,
          -1.08216397e-02,  5.31121836e-03,  1.62021374e-03,  2.89892269e-03,
          -1.74084656e-01,  1.23787924e-02,  4.44709934e-01, -1.82155625e-03,
           1.35372354e-02,  1.02045982e-02, -9.19553019e-03, -1.02826186e-02,
          -4.61509962e-03,  1.73050212e-02, -1.15431895e-02, -6.52759960e-03,
          -9.60369460e-02,  2.60129749e-01,  5.22192232e-01, -9.22743394e-03,
           1.09619679e-02, -9.43119594e-03,  8.92387796e-04,  4.24418192e-04,
          -1.55860336e-04,  6.96227600e-03, -5.67669229e-03,  7.90324480e-03,
           1.68163899e-01,  2.30974708e-01,  4.86910323e-01,  4.89276667e-04,
          -5.18664979e-03, -1.12686056e-02, -9.88594225e-04, -1.12747180e-02,
          -6.82377576e-03, -3.03595997e-03, -1.81014178e-03,  4.81567162e-03,
          -1.37554394e-02,  1.24789998e-01,  3.41477515e-01,  5.38052796e-03,
           1.43447515e-03,  8.55733369e-03, -1.01180268e-02, -5.27708505e-03,
          -7.16161152e-03,  1.20479851e-02, -7.22517313e-03,  4.19344403e-03,
          -9.75354004e-02,  8.24227048e-02,  5.84778798e-02, -8.78997982e-04,
           7.74431155e-03, -4.37950188e-03, -3.42111426e-03, -1.67640365e-02,
          -5.69899113e-03, -4.98792314e-05,  2.91874144e-03, -4.09997849e-04,
          -4.89743303e-02,  1.66197757e-01,  1.51676819e-01, -6.57680975e-03,
           1.83611457e-03, -9.49198148e-03, -2.12344241e-03, -1.50535683e-02,
          -9.27848881e-04, -8.71779329e-03,  4.77088219e-03,  8.32933609e-03,
           1.70474627e-02,  1.69783618e-01,  2.30528143e-01, -2.83591561e-03,
          -2.46127477e-03, -9.69169617e-03,  2.24209101e-03, -1.59381960e-02,
           1.14228823e-03, -1.04550379e-02,  1.99072884e-03,  1.15559992e-02,
          -1.15247816e-01,  1.55024739e-01,  1.13602949e-01,  5.46319865e-03,
           5.53200681e-04, -1.67292005e-04, -1.02008927e-03, -1.21755818e-02,
          -1.46379182e-03, -1.29844626e-03, -1.23247302e-03,  6.28775835e-03;

    const SparseMatrix *MASS = oc.getHomogenizedMassMatrix( );

    if ( ( ( *MASS ) * MX - MA ).norm( ) > 1e-3 * ( MA.norm( ) + 1 ) ){

        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;

    }
    testNum++;

    MA = Eigen::MatrixXd::Zero( 144, 1 );

    MASS = oc.getFreeMicromorphicMassMatrix( );

    if ( ( ( *MASS ) * MX - MA ).norm( ) > 1e-3 * ( MA.norm( ) + 1 ) ){

        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;

    }
    testNum++;

    //Test of the micro mass matrices
    Eigen::MatrixXd MQX;
    MQX = Eigen::MatrixXd::Ones( 81, 1 );
    Eigen::MatrixXd MQA( 81, 1 );
    MQA << 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.0625 ,
           0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.0625 ,
           0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.125  , 0.125  , 0.125  ,
           0.125  , 0.125  , 0.125  , 0.03125, 0.03125, 0.03125, 0.03125,
           0.03125, 0.03125, 0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.0625 ,
           0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.0625 ,
           0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125,
           0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125, 0.03125,
           0.03125, 0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.0625 ,
           0.125  , 0.125  , 0.125  , 0.03125, 0.03125, 0.03125, 0.0625 ,
           0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.0625 , 0.03125, 0.03125,
           0.03125, 0.03125, 0.03125, 0.03125;

    if ( ( oc._test_MQ * MQX - MQA ).norm( ) > 1e-6 * ( MQA.norm( ) + 1 ) ){

        std::cout << oc._test_MQ << "\n";
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;

    }
    testNum++;

    Eigen::MatrixXd MQhatX;
    MQhatX = Eigen::MatrixXd::Ones( 54, 1 );
    Eigen::MatrixXd MQhatA( 54, 1 );
    MQhatA << 0.015625, 0.015625, 0.015625, 0.03125 , 0.03125 , 0.03125 ,
              0.03125 , 0.03125 , 0.03125 , 0.0625  , 0.0625  , 0.0625  ,
              0.03125 , 0.03125 , 0.03125 , 0.0625  , 0.0625  , 0.0625  ,
              0.0625  , 0.0625  , 0.0625  , 0.125   , 0.125   , 0.125   ,
              0.015625, 0.015625, 0.015625, 0.03125 , 0.03125 , 0.03125 ,
              0.03125 , 0.03125 , 0.03125 , 0.0625  , 0.0625  , 0.0625  ,
              0.03125 , 0.03125 , 0.03125 , 0.0625  , 0.0625  , 0.0625  ,
              0.015625, 0.015625, 0.015625, 0.03125 , 0.03125 , 0.03125 ,
              0.015625, 0.015625, 0.015625, 0.03125 , 0.03125 , 0.03125;

    if ( ( oc._test_MQhat * MQhatX - MQhatA ).norm( ) > 1e-6 * ( MQhatA.norm( ) + 1 ) ){

        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;

    }
    testNum++;

    Eigen::MatrixXd blockQ( 81, 1 );
    blockQ <<  1.0011287 , -0.48217873,  1.0658671 , -0.59590148,  1.86397241,
               0.83663432,  1.35532451, -1.12671449,  0.29344055,  0.95589525,
              -0.23928356,  1.42599166,  0.1912066 ,  1.06965081,  1.01712103,
              -1.18132365, -1.21585256,  1.20078413,  0.33294429,  0.37493735,
               0.25446677, -1.2453524 , -0.99254126,  0.07586945,  0.14820318,
               0.14313054,  0.6756101 ,  1.52370397, -0.02636185,  1.45696374,
               1.36324949, -0.33198489,  1.51940143, -0.76151492, -0.40468061,
              -0.18816881, -0.71433009,  0.60194972, -1.65484134, -0.06390784,
               0.38962463,  1.7471683 ,  0.01489609, -0.24756449,  0.85431249,
               0.60728149, -1.08881333,  1.20554868,  1.05935139,  0.45570988,
               0.55916691, -0.37152396, -1.81866729,  1.62818807, -0.83506505,
              -0.41756133,  0.68115436, -0.69577408,  0.00356817,  0.10222119,
               1.48197214,  0.82315746,  1.52177125,  1.07312683, -1.81523853,
               1.64862571,  1.0919167 , -1.5641087 ,  1.77508949, -0.6794552 ,
               1.49221457, -1.56626164, -0.14480897,  1.92560457,  0.21194309,
               1.57820827,  0.16747939, -0.5245916 ,  1.52191223,  1.11487778,
               0.0923271 ;

    Eigen::MatrixXd blockD( 48, 1 );
    blockD << 1.12887639, -1.87755231,  1.65442528,  0.04342622,  0.60115804,
             -0.43565936,  1.08036051, -1.33575861, -1.44532869,  1.73659983,
              1.1254876 ,  0.45722009,  0.11743958, -1.98988391,  1.79062452,
             -1.0165452 , -1.62741096,  0.82120009, -0.45170424,  1.4438855 ,
              1.96121981,  0.7152962 , -1.57750846, -1.48106602, -0.79377388,
              0.63141761, -1.48792601,  1.38464132, -0.71136603, -0.23120376,
             -1.85495407, -0.68366174, -0.88327987,  1.95952626,  0.58011715,
             -1.50747095,  0.94703545,  1.57112632,  0.45692654, -1.50705853,
             -0.34519697,  1.93386415,  0.62898342,  0.10554378,  1.56442787,
              0.185254  ,  0.15645488, -1.3053099;

    Eigen::MatrixXd MQQA( 81,1 );
    MQQA <<  0.05971404, -0.00634264,  0.07567346, -0.02312268,  0.07946925,
             0.05725601,  0.17421939, -0.12251207,  0.05275722,  0.10101681,
            -0.05784871,  0.13633191, -0.01855223,  0.09252345,  0.14506602,
            -0.09723474, -0.0908336 ,  0.16430976,  0.0058038 ,  0.0446031 ,
             0.02926295, -0.208282  , -0.19572726,  0.07179328,  0.04300342,
             0.00201703,  0.06173752,  0.07983716, -0.00548067,  0.0837787 ,
             0.12478498, -0.04690799,  0.17904863, -0.04181574, -0.0482353 ,
             0.02175389, -0.08616034,  0.06672298, -0.15381323, -0.04211104,
             0.03117406,  0.12550784,  0.01070461, -0.02753107,  0.06434679,
             0.02785674, -0.04336496,  0.06063639,  0.04895516,  0.00655007,
             0.04899355, -0.01006352, -0.06625401,  0.08254097, -0.04220474,
             0.00683894,  0.03022279, -0.04268879, -0.03011562,  0.04085348,
             0.134798  ,  0.03287278,  0.14974673,  0.14406718, -0.29309928,
             0.26318819,  0.04883612, -0.05268398,  0.07437171, -0.05680173,
             0.10015744, -0.1452058 ,  0.01769978,  0.16191392,  0.04986132,
             0.07071237,  0.02726372, -0.03155772,  0.0730217 ,  0.06246673,
             0.01150372;

    Eigen::MatrixXd MQDA( 81, 1 );
    MQDA << 4.25412368e-03, -4.83706217e-03,  4.65503691e-03, -6.24424127e-04,
            6.91429725e-04, -9.97329797e-04,  4.95674094e-03, -9.85705697e-03,
            1.15710905e-02, -1.34479347e-03,  1.69563166e-03, -1.79963383e-03,
            6.60369722e-03, -3.70328248e-03,  6.11858194e-03, -1.29577009e-03,
            9.02223570e-04, -1.68768615e-03,  4.49150134e-03, -6.45745867e-03,
            1.30548871e-02, -2.79573235e-03,  3.11459917e-03, -3.40666392e-03,
            1.21354535e-03, -6.86671168e-03,  7.01574954e-03, -5.31369041e-04,
            4.00832953e-04, -7.53312353e-04, -3.45792566e-04, -6.62695535e-03,
            7.66242056e-03, -9.28281130e-04,  9.58984158e-04, -1.48397680e-03,
           -4.94468731e-04,  2.73782707e-03,  2.63000058e-03, -1.46039863e-03,
            1.20479098e-03, -1.23607110e-03, -1.39850209e-03, -4.27338525e-04,
            1.51432490e-03, -3.56001707e-04,  3.63604163e-04, -4.68562575e-04,
            2.15950859e-03,  1.13910331e-03,  1.74209513e-03, -7.21716094e-04,
            1.91164357e-04, -5.81499228e-04, -6.84854880e-04,  8.62630704e-04,
           -5.65334209e-04, -6.90947683e-04,  1.41227146e-03, -1.64549941e-03,
           -1.19867050e-03,  4.57421389e-04, -5.17248011e-04, -2.59175561e-04,
           -6.96433922e-06, -8.97018306e-04, -1.48703333e-04,  1.19718488e-03,
           -1.08793337e-03,  3.11091313e-04,  9.13367827e-04, -6.38088181e-04,
            4.42170477e-04, -1.18380857e-03,  3.40715374e-04,  3.84123841e-04,
           -1.11352132e-05,  1.10195132e-04, -4.27746456e-04, -4.42463925e-04,
           -2.00294777e-05;

    Eigen::MatrixXd MDQA( 48, 1 );
    MDQA <<  1.32144420e-02,  4.96725596e-03,  2.44421827e-03,  1.85640648e-05,
            -2.03399171e-03, -1.23911385e-03, -4.87928752e-04,  1.42563921e-03,
             1.15926905e-05, -8.41448936e-04, -1.02217996e-03, -1.86229352e-03,
             1.39302514e-02,  2.48464733e-03,  4.38463273e-03,  3.39340638e-04,
            -1.66853177e-03, -1.69299470e-03, -7.53375562e-04,  1.02011287e-03,
             4.88358839e-04,  1.81165617e-03, -6.64278407e-04, -1.86779630e-03,
             8.75502505e-03,  4.91124795e-03,  7.67494716e-04,  5.65920359e-04,
            -9.19081417e-04, -6.24080196e-04, -1.18815244e-03,  1.93187435e-04,
             9.52826063e-05,  2.24408941e-03, -1.14429060e-03, -1.23232108e-03,
             7.81347114e-03,  8.34081386e-03, -1.96538764e-03, -9.51434045e-05,
            -6.66493724e-04, -2.77222083e-04, -5.26630514e-04,  2.61139735e-04,
            -3.98920678e-04, -8.77648238e-04, -1.18262299e-03, -1.18055443e-03;

    Eigen::MatrixXd MDDA( 48, 1 );
    MDDA <<  8.22523748e-02, -1.06405556e-01,  1.36382744e-01, -5.20077281e-03,
            -3.20320083e-03,  5.56583707e-04,  1.22964613e-03,  3.87050016e-03,
             3.98982808e-03,  7.17758133e-03, -2.38869519e-03, -9.26531978e-03,
             3.59178836e-02, -1.22479549e-01,  1.25187239e-01, -4.13342988e-03,
            -6.42792210e-03,  1.74559291e-03, -3.82161886e-03,  6.78972927e-03,
             7.24077649e-03,  7.20461385e-03, -6.32056251e-03, -1.23396557e-02,
             8.22906150e-03, -2.02034579e-02,  3.59646849e-02, -2.32339095e-03,
            -5.81384762e-03,  2.85308645e-03, -5.82937262e-03,  4.61477744e-03,
             3.61520050e-03,  7.06669287e-03, -3.62735911e-03, -1.20293241e-02,
             6.25874855e-02,  3.78425533e-03,  6.98456179e-02, -6.16914072e-03,
            -3.89857351e-03,  3.65918980e-03, -1.08849087e-04,  4.20479456e-03,
             3.94116334e-03,  5.07410030e-03, -1.89383480e-03, -1.03132406e-02;

    if ( ( oc._test_dense_MQQ * blockQ - MQQA ).norm( ) > 1e-4 * ( MQQA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( oc._test_dense_MQD * blockD - MQDA ).norm( ) > 1e-4 * ( MQDA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( oc._test_dense_MDQ * blockQ - MDQA ).norm( ) > 1e-4 * ( MDQA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( oc._test_dense_MDD * blockD - MDDA ).norm( ) > 2e-4 * ( MDDA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    Eigen::MatrixXd CQQA( 81,1 );
    CQQA <<  0.00817607,  0.00026187,  0.01118441, -0.0029048 ,  0.01009808,
             0.00867662,  0.02503721, -0.01703899,  0.00853178,  0.01395977,
            -0.00991062,  0.01848636, -0.00484049,  0.01124219,  0.02175531,
            -0.01191313, -0.01070831,  0.02508456, -0.00298069,  0.00364059,
             0.00313561, -0.02575939, -0.02690848,  0.01385988,  0.00740482,
             0.00014861,  0.00959539,  0.01106456, -0.00096424,  0.0120508 ,
             0.01582144, -0.0071143 ,  0.02503935, -0.00366011, -0.00706069,
             0.00515751, -0.01254446,  0.00919911, -0.01951676, -0.0078921 ,
             0.00372743,  0.01460382,  0.00191605, -0.00455106,  0.00976626,
             0.00363611, -0.00518773,  0.00822522,  0.00634638, -0.00026499,
             0.00782357, -0.00090503, -0.00755122,  0.01132992, -0.0057281 ,
             0.00266321,  0.00402134, -0.00397984, -0.00613005,  0.00752077,
             0.01765461,  0.00152601,  0.02060092,  0.01531913, -0.03582454,
             0.03187299,  0.00646829, -0.00568165,  0.00943298, -0.00703331,
             0.01067067, -0.01902579,  0.00438981,  0.02040663,  0.00844939,
             0.00923888,  0.00489003, -0.0045931 ,  0.0098756 ,  0.00902198,
             0.00205837;

    Eigen::MatrixXd CQDA( 81, 1 );
    CQDA <<  7.08435052e-04, -8.28217337e-04,  6.91061938e-04,  8.84113776e-05,
            -1.13071924e-04,  2.05012990e-05,  8.73844896e-04, -1.54683275e-03,
             1.74025113e-03,  6.77726809e-05, -1.29313339e-04,  1.29454831e-04,
             1.23361666e-03, -5.62100292e-04,  8.19698219e-04,  1.74883709e-04,
            -1.15117596e-05, -5.73697310e-05,  8.78176724e-04, -8.04815747e-04,
             1.60569448e-03, -4.72642683e-05,  2.21518912e-04, -1.77468385e-04,
             2.38427320e-04, -1.06090086e-03,  1.05153572e-03,  9.56055494e-06,
            -1.57807034e-04,  1.09924032e-04, -2.26881561e-05, -9.72262416e-04,
             9.22780541e-04, -8.45253782e-05, -6.87021928e-05, -6.35261303e-05,
            -1.17430104e-06,  6.17353765e-04,  8.13631327e-05, -1.17314225e-04,
             2.99272862e-04, -2.17621089e-04, -2.19869257e-04, -5.80074161e-05,
             5.30721557e-05, -7.70269352e-05,  2.84536745e-05, -9.82483087e-05,
             4.78429263e-04,  2.88099663e-04,  1.62728726e-04,  6.71360572e-05,
             1.10651942e-04, -6.37707573e-05, -9.81595666e-05,  1.17051162e-04,
            -8.02186930e-05, -1.11324363e-04,  1.99268594e-04, -2.20398299e-04,
            -1.75873394e-04,  6.40286869e-05, -7.58495621e-05, -8.46155980e-05,
             3.54407103e-05, -1.33201470e-04, -2.82644221e-05,  1.53000242e-04,
            -1.40664856e-04,  2.24464965e-05,  1.22278456e-04, -8.56381698e-05,
             2.78474029e-05, -1.38048189e-04,  4.25457654e-05,  4.21814196e-05,
            -3.96202904e-07,  1.74257592e-05, -6.80456903e-05, -5.75683631e-05,
            -2.68100544e-06;

    Eigen::MatrixXd CDQA( 48, 1 );
    CDQA <<  1.68643108e-03,  3.67726167e-04,  1.21931810e-03,  4.55228169e-05,
            -3.79650513e-04, -3.26646668e-05, -1.22082884e-04,  1.87921502e-04,
            -3.10905881e-05, -1.83562202e-04, -1.59292727e-04,  5.74947149e-05,
             1.98848996e-03, -7.66647939e-05,  1.40699441e-03,  1.05506627e-04,
            -3.31294035e-04, -2.04332640e-05, -1.00112596e-04,  1.37823563e-04,
            -3.24464978e-05,  2.77400361e-04, -1.20067259e-04,  5.53531546e-05,
             9.92563336e-04,  2.84409532e-04,  7.99249106e-04,  1.49730304e-04,
            -1.66669280e-04, -3.31428120e-05, -1.85516179e-04,  4.27136002e-05,
            -4.21505552e-05,  3.51167769e-04, -1.83805395e-04,  5.73928002e-05,
             7.14964222e-04,  8.21268716e-04,  5.13488747e-04, -1.09307470e-05,
            -1.06082914e-04, -4.43179315e-05, -8.29134130e-05,  3.88497717e-05,
            -4.30952353e-05, -2.08287589e-04, -1.93621947e-04,  5.94514343e-05;

    Eigen::MatrixXd CDDA( 48, 1 );
    CDDA <<  1.35554199e-02, -1.77165212e-02,  2.24811070e-02, -7.13310416e-04,
            -4.47176794e-04,  2.38525760e-04,  2.85607075e-04,  2.02807421e-04,
             4.38552065e-04,  1.33700625e-03, -9.81203748e-05, -1.24770181e-03,
             5.46775199e-03, -2.03940081e-02,  1.98101172e-02, -5.63914978e-04,
            -1.04718954e-03,  3.85271442e-04, -5.35345793e-04,  7.71628745e-04,
             1.01828223e-03,  1.32348839e-03, -7.77638768e-04, -1.79189088e-03,
             8.62935199e-04, -2.17211855e-03,  3.75629490e-03, -1.93782350e-04,
            -9.34690595e-04,  5.16307249e-04, -9.23302071e-04,  3.92608258e-04,
             4.87375207e-04,  1.36526785e-03, -3.16601297e-04, -1.85375758e-03,
             1.04228243e-02,  2.05783831e-03,  1.06161247e-02, -8.48826106e-04,
            -5.72986996e-04,  7.65721748e-04,  1.94506180e-05,  2.80450124e-04,
             6.02876147e-04,  1.02106916e-03, -3.71641385e-05, -1.54751459e-03;

    if ( ( oc._test_dense_CQQ * blockQ - CQQA ).norm( ) > 1e-4 * ( CQQA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( oc._test_dense_CQD * blockD - CQDA ).norm( ) > 1e-4 * ( CQDA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( oc._test_dense_CDQ * blockQ - CDQA ).norm( ) > 1e-4 * ( CDQA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( oc._test_dense_CDD * blockD - CDDA ).norm( ) > 1e-4 * ( CDDA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    Eigen::MatrixXd QD( blockQ.size( ) + blockD.size( ), 1 );
    QD << blockQ, blockD;

    Eigen::MatrixXd L2MASSA( blockQ.size( ) + blockD.size( ), 1 );
    L2MASSA << MQQA + MQDA, MDQA + MDDA;

    Eigen::MatrixXd L2DAMPINGA( blockQ.size( ) + blockD.size( ), 1 );
    L2DAMPINGA << CQQA + CQDA, CDQA + CDDA;

    if ( ( ( *oc.getMass( ) ) * QD - L2MASSA ).norm( )  > 2e-4 * ( L2MASSA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( ( *oc.getDamping( ) ) * QD - L2DAMPINGA ).norm( )  > 1e-4 * ( L2DAMPINGA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    Eigen::MatrixXd FintQA( 81, 1 );
    FintQA << 3.75780853e-01,  1.06556668e+00, -1.16829679e-01,  6.22032962e-01,
              1.23021166e+00, -5.99072345e-01,  8.46752728e-01,  5.85982914e-01,
              7.15721545e-02,  1.09300484e+00,  7.50627895e-01, -4.10670511e-01,
              2.73121662e-01,  1.48027906e+00, -3.78334926e-01,  5.19373771e-01,
              1.64492404e+00, -8.60577592e-01,  7.44093537e-01,  1.00069530e+00,
             -1.89933092e-01,  9.90345645e-01,  1.16534028e+00, -6.72175758e-01,
              1.31772460e+00,  1.06399153e-01,  2.59973988e-01,  1.56397671e+00,
              2.71044134e-01, -2.22268678e-01,  1.21506541e+00,  5.21111539e-01,
             -1.53125819e-03,  1.46131752e+00,  6.85756520e-01, -4.83773924e-01,
              6.41434346e-01,  1.41540769e+00, -4.51438338e-01,  8.87686454e-01,
              1.58005267e+00, -9.33681004e-01,  1.11240622e+00,  9.35823925e-01,
             -2.63036505e-01,  1.35865833e+00,  1.10046891e+00, -7.45279171e-01,
              1.70462471e-01,  1.89499145e+00, -6.39840172e-01,  4.16714580e-01,
              2.05963643e+00, -1.12208284e+00,  8.68285070e-01,  1.39485664e+00,
             -1.08131501e+00,  1.33925694e+00,  9.15272876e-01, -8.92913177e-01,
              7.65625879e-01,  1.80956902e+00, -1.34282026e+00,  1.23659775e+00,
              1.32998526e+00, -1.15441842e+00,  1.81022882e+00,  4.35689115e-01,
             -7.04511344e-01,  1.70756963e+00,  8.50401502e-01, -9.66016590e-01,
              1.13393856e+00,  1.74469765e+00, -1.41592367e+00,  1.60491044e+00,
              1.26511389e+00, -1.22752184e+00,  6.62966688e-01,  2.22428141e+00,
             -1.60432550e+00;

    Eigen::MatrixXd FintQhatA( 54, 1 );
    FintQhatA << -0.11672336,  0.73627671,  0.84765565,  0.12952874,  0.90092169,
                  0.36541299,  0.35424851,  0.25669295,  1.03605749,  0.60050062,
                  0.42133793,  0.55381482, -0.21938255,  1.1509891 ,  0.58615041,
                  0.02686955,  1.31563408,  0.10390774,  0.25158932,  0.67140534,
                  0.77455224,  0.49784143,  0.83605032,  0.29230957,  0.82522038,
                 -0.22289081,  1.22445932,  1.07147249, -0.05824583,  0.74221665,
                  0.72256119,  0.19182158,  0.96295407,  0.9688133 ,  0.35646656,
                  0.48071141,  0.14893013,  1.08611772,  0.51304699,  0.39518224,
                  1.25076271,  0.03080433,  0.619902  ,  0.60653396,  0.70144883,
                  0.86615411,  0.77117894,  0.21920616, -0.32204175,  1.56570149,
                  0.32464516, -0.07578964,  1.73034647, -0.15759751;

    Eigen::MatrixXd FextQA( 81,1 );
    FextQA << 0.65061584, -1.77812354, -0.40151144,  1.05334699, -2.09057326,
             -0.25620483,  0.26440628, -2.03750988, -0.29634106,  0.66713743,
             -2.3499596 , -0.15103445,  1.11996119, -1.50248801,  0.0487764 ,
              1.52269233, -1.81493773,  0.19408302,  0.73375163, -1.76187435,
              0.15394679,  1.13648278, -2.07432407,  0.2992534 , -0.12180327,
             -2.29689622, -0.19117068,  0.28092787, -2.60934594, -0.04586407,
              0.34754207, -2.02126069,  0.25911717,  0.75027322, -2.33371041,
              0.40442378,  1.20309698, -1.48623882,  0.60423463,  1.60582812,
             -1.79868854,  0.74954125,  0.81688742, -1.74562517,  0.70940502,
              1.21961857, -2.05807488,  0.85471163,  1.58930653, -1.22685248,
              0.49906425,  1.99203768, -1.5393022 ,  0.64437087,  1.45607813,
             -2.40302297, -0.11089821,  1.06986858, -2.66240931, -0.00572783,
              1.92542348, -2.12738744,  0.33938963,  1.53921392, -2.38677378,
              0.44456002,  0.68365902, -2.92179565,  0.09944255,  1.15300437,
             -2.64616012,  0.5497304 ,  2.00855927, -2.11113825,  0.89484786,
              1.62234971, -2.37052459,  1.00001825,  2.39476883, -1.85175191,
              0.78967748;

    Eigen::MatrixXd FextQhatA( 54, 1 );
    FextQhatA << -0.15484645, -1.15322411, -0.69212468,  0.24788469, -1.46567383,
                 -0.54681806, -0.54105601, -1.41261045, -0.58695429, -0.13832486,
                 -1.72506017, -0.44164768,  0.31449889, -0.87758858, -0.24183683,
                  0.71723004, -1.1900383 , -0.09653021, -0.07171066, -1.13697492,
                 -0.13666645,  0.33102048, -1.44942464,  0.00864017, -0.92726557,
                 -1.6719968 , -0.48178391, -0.52453442, -1.98444651, -0.3364773 ,
                 -0.45792022, -1.39636127, -0.03149606, -0.05518907, -1.70881098,
                  0.11381055,  0.39763468, -0.8613394 ,  0.3136214 ,  0.80036583,
                 -1.17378911,  0.45892802,  0.01142513, -1.12072574,  0.41879178,
                  0.41415627, -1.43317545,  0.5640984 ,  0.78384424, -0.60195305,
                  0.20845102,  1.18657539, -0.91440277,  0.35375764;

    if ( ( oc._test_FintQ - FintQA ).norm( ) > 1e-6 * ( FintQA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( oc._test_FintQhat - FintQhatA ).norm( ) > 1e-6 * ( FintQhatA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( oc._test_FextQ - FextQA ).norm( ) > 1e-6 * ( FextQA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( oc._test_FextQhat - FextQhatA ).norm( ) > 1e-6 * ( FextQhatA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    Eigen::MatrixXd FintDA( 48, 1 );
    FintDA << -1.50000000e+00, -1.87500000e+00, -2.25000000e+00,  0.00000000e+00,
               0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
               0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
              -1.25000000e+00, -1.37500000e+00, -1.50000000e+00,  0.00000000e+00,
               0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
               0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
              -2.50000000e-01, -1.25000000e-01,  2.56739074e-16,  0.00000000e+00,
               0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
               0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
              -5.00000000e-01, -6.25000000e-01, -7.50000000e-01,  0.00000000e+00,
               0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
               0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00;

    Eigen::MatrixXd FintDhatA( 96, 1 );
    FintDhatA << -1.25000000e+00, -1.75000000e+00, -2.25000000e+00,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                 -7.50000000e-01, -7.50000000e-01, -7.50000000e-01,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  1.25000000e+00,  1.75000000e+00,  2.25000000e+00,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  7.50000000e-01,  7.50000000e-01,  7.50000000e-01,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  2.50000000e-01,  1.25000000e-01, -2.22044605e-16,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  5.00000000e-01,  6.25000000e-01,  7.50000000e-01,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  1.50000000e+00,  1.87500000e+00,  2.25000000e+00,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  1.25000000e+00,  1.37500000e+00,  1.50000000e+00,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00;

    Eigen::MatrixXd FextDA( 48, 1 );
    FextDA << -1.625, -2.125, -2.625,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
               0.   ,  0.   ,  0.   ,  0.   , -1.375, -1.625, -1.875,  0.   ,
               0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
              -0.375, -0.375, -0.375,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
               0.   ,  0.   ,  0.   ,  0.   , -0.625, -0.875, -1.125,  0.   ,
               0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ;

    Eigen::MatrixXd FextDhatA( 96,1 );
    FextDhatA << -1.5  , -2.25 , -3.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                  0.   ,  0.   ,  0.   ,  0.   , -1.   , -1.25 , -1.5  ,  0.   ,
                  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                  1.   ,  1.25 ,  1.5  ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                  0.   ,  0.   ,  0.   ,  0.   ,  0.5  ,  0.25 ,  0.   ,  0.   ,
                  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                  0.125, -0.125, -0.375,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                  0.   ,  0.   ,  0.   ,  0.   ,  0.375,  0.375,  0.375,  0.   ,
                  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                  1.375,  1.625,  1.875,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
                  0.   ,  0.   ,  0.   ,  0.   ,  1.125,  1.125,  1.125,  0.   ,
                  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ;

    if ( ( oc._test_FintD - FintDA ).norm( ) > 1e-2 * ( FintDA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( oc._test_FintDhat - FintDhatA ).norm( ) > 1e-2 * ( FintDhatA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( oc._test_FextD - FextDA ).norm( ) > 1e-2 * ( FextDA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( oc._test_FextDhat - FextDhatA ).norm( ) > 1e-2 * ( FextDhatA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    Eigen::MatrixXd FORCEA( 129, 1 );
    FORCEA <<  7.41572446e-01, -3.37882955e+00, -6.64396274e-02,  6.75585483e-01,
              -3.15872885e+00,  7.06987166e-01, -3.96516507e-01, -2.97209062e+00,
              -1.06275977e-01, -1.43272296e-01, -2.63537046e+00,  7.57288269e-01,
               1.53352065e+00, -4.07299625e+00,  7.49033914e-01,  1.05028350e+00,
              -3.40772068e+00,  1.25376301e+00,  2.41003195e-02, -3.54414684e+00,
               2.93374201e-01, -1.07155659e-01, -2.87156906e+00,  7.65913394e-01,
              -1.96731252e+00, -3.05311249e+00, -8.05579873e-01, -1.34683094e+00,
              -2.92313814e+00,  1.45398784e-01, -1.69987143e+00, -3.98269697e+00,
              -5.17145693e-01, -1.08579328e+00, -3.42687090e+00,  3.15924535e-01,
               2.30165721e-01, -5.08360261e+00,  3.38164198e-01,  1.07762519e-01,
              -4.19922112e+00,  8.12399276e-01, -9.25144345e-01, -4.28607508e+00,
               3.23580687e-01, -5.84448865e-01, -3.85991038e+00,  9.65458262e-01,
               1.78374062e+00, -4.61179214e+00,  1.06272093e+00,  1.43796755e+00,
              -4.09550109e+00,  1.52704664e+00,  5.61451546e-01, -3.61323670e+00,
               1.10742672e+00, -3.01573985e-01, -3.53987249e+00,  9.05627735e-01,
               7.99527252e-01, -4.04978573e+00,  1.43548718e+00, -4.09266167e-01,
              -4.59617147e+00,  6.65651287e-01, -1.08136102e+00, -3.33136346e+00,
               7.67655549e-01, -8.68931626e-01, -3.90146570e+00,  9.12970047e-01,
               2.32169612e-01, -4.41137894e+00,  1.44282950e+00, -2.74914751e-01,
              -3.89830852e+00,  1.70769408e+00,  1.36789782e+00, -4.18018175e+00,
               2.04746526e+00, -6.97556127e-01, -2.21689262e+00, -1.94240251e+00,
              -1.66270848e-01,  6.59146956e-02,  1.82958084e-02,  2.77961583e-01,
              -8.40178319e-02,  2.13630595e-01, -8.22957140e-03,  1.73876884e-01,
               3.26770915e-01, -1.53165889e+00, -1.36753800e+00, -2.27606679e+00,
              -2.50780532e-01,  1.01700599e-01,  1.18845411e-01,  1.46715729e-01,
              -1.20458093e-01,  1.11242408e-01, -1.58602567e-01,  1.88192240e-01,
               3.66993550e-01, -1.34984694e+00, -2.10200382e+00, -1.81293457e+00,
              -2.42980168e-01, -1.07946245e-02,  9.69283041e-02,  1.15204596e-01,
              -2.46774819e-01,  1.99780955e-01, -1.38732587e-01,  4.33738684e-02,
               3.11163786e-01, -5.15744178e-01, -2.95135844e+00, -1.47927030e+00,
              -1.74071212e-01,  2.49912790e-02, -3.62129833e-03,  3.09472716e-01,
              -2.83215080e-01,  3.02169142e-01, -2.80995510e-02,  5.76892239e-02,
               2.70941152e-01;

    if ( ( *oc.getFORCE( ) - FORCEA ).norm( ) > 1e-2 * ( FORCEA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }

    Eigen::MatrixXd DOFTA( 129, 1 );
    DOFTA = Eigen::MatrixXd::Zero( 129, 1 );

    Eigen::MatrixXd DotDOFTA( 129, 1 );
    DotDOFTA <<  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
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
                 0.        ,  0.        ,  0.        , -0.        ,  0.33024635,
                -0.55968633, -0.72735076,  1.24046985,  1.13823231,  0.27563121,
                 1.15197975,  0.13506319, -0.23919711,  0.        ,  0.        ,
                -0.        , -1.38293804,  0.08902328, -0.66596965, -0.33235995,
                 0.6922183 ,  1.04286003,  0.26561463, -1.83426801,  0.38520627,
                 0.        ,  0.        , -0.        , -3.15697149,  1.9007297 ,
                -2.5389396 , -1.22303943,  1.0373448 ,  2.06632253,  0.63118904,
                -2.03612508,  0.88366176,  0.        ,  0.        , -0.        ,
                -1.44378711,  1.25202009, -2.6003207 ,  0.34979037,  1.48335881,
                 1.29909371,  1.51755416, -0.06679388,  0.25925838;

    Eigen::MatrixXd DotDotDOFTA( 129, 1 );
    DotDotDOFTA <<  0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
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
                    0.        ,  0.        ,  0.        , -0.        , -0.89718641,
                    0.42644159, -1.67788156, -1.19866342,  0.03566113,  1.00466326,
                    0.47680901, -0.31530177,  1.38300667,  0.        ,  0.        ,
                   -0.        , -0.1124791 ,  2.26606881, -2.59540061,  0.77169844,
                    1.52953487,  0.74842145,  2.34750761, -2.18910047,  3.20318666,
                    0.        ,  0.        , -0.        ,  0.61956589,  0.48823537,
                   -0.63236534,  2.25807607,  0.46795257,  0.13635096,  3.35049415,
                   -0.20548571,  3.85886988,  0.        ,  0.        , -0.        ,
                   -0.16514143, -1.35139186,  0.28515371,  0.28771421, -1.02592116,
                    0.39259277,  1.47979555,  1.668313  ,  2.03868988;

    if ( !DOFTA.isApprox( oc._test_DOF_t ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( DotDOFTA - oc._test_DotDOF_t ).norm( ) > 1e-6 * ( DotDOFTA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( DotDotDOFTA - oc._test_DotDotDOF_t ).norm( ) > 1e-6 * ( DotDotDOFTA.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    Eigen::MatrixXd DotDotDOFTP1A( 129, 1 );
    DotDotDOFTP1A << 1.03430057e+01, -4.82813189e+01, -1.11898072e+00,  9.03736566e+00,
                    -4.63731522e+01,  1.21740548e+01, -1.97876122e+00, -1.99671032e+01,
                     6.41970944e-02, -2.09745568e-01, -1.55646078e+01,  6.50581602e+00,
                     1.00855313e+01, -2.59295640e+01,  5.06581401e+00,  5.04582239e+00,
                    -2.02594855e+01,  8.41424597e+00,  5.13431613e-01, -1.94754691e+01,
                     2.49816817e+00, -5.14771466e-01, -1.48675425e+01,  4.44886800e+00,
                    -2.71979247e+01, -4.50305995e+01, -1.07226573e+01, -1.74243964e+01,
                    -4.29323702e+01,  4.78676493e+00, -1.02163985e+01, -2.69745851e+01,
                    -2.91654779e+00, -5.90033161e+00, -2.07779072e+01,  2.52357188e+00,
                     1.84732426e+00, -3.29371480e+01,  2.08488097e+00, -6.45468790e-01,
                    -2.54729113e+01,  4.43176896e+00, -1.17897823e+01, -6.06724600e+01,
                     5.66129401e+00, -7.49074520e+00, -5.53330367e+01,  1.56069825e+01,
                     2.57511845e+01, -6.39231446e+01,  1.52645086e+01,  1.89711177e+01,
                    -5.87737224e+01,  2.29929915e+01,  1.09351248e+01, -7.15733152e+01,
                     2.11969429e+01, -3.08908106e+00, -3.23611178e+01,  7.29757610e+00,
                     7.28987487e+00, -3.73298598e+01,  1.23093909e+01, -1.66556850e+00,
                    -2.72176845e+01,  4.57909428e+00, -2.13494012e+01, -6.57375962e+01,
                     1.45478537e+01, -8.51578856e+00, -3.58016517e+01,  7.29005404e+00,
                     1.86250767e+00, -4.07705119e+01,  1.23016510e+01, -5.44518885e+00,
                    -7.68907928e+01,  3.32090276e+01,  2.68393530e+01, -8.27264966e+01,
                     3.98579143e+01, -7.80470522e-01, -1.28365059e+01, -1.20856277e+01,
                    -2.02038719e+00, -5.54418324e+00, -6.69569577e+00,  2.43883196e+01,
                     2.77756586e+01,  1.01376942e+01,  2.55757186e+01,  2.37956704e+01,
                     3.31379321e+01, -1.52331473e+01,  3.99479178e+00, -1.80066878e+01,
                    -3.62368979e+01,  3.88325488e+01,  1.55737751e+01,  2.03688451e+01,
                    -1.92952995e+01, -2.20356705e+01, -4.18797125e+01,  4.27848551e+01,
                     4.25235872e+01, -1.22756695e+01, -9.48124613e+00, -9.44039734e+00,
                    -2.71537739e+01, -3.15581059e+01,  1.17918589e+01, -2.03818665e+01,
                    -1.71368574e+01,  8.64168121e-01, -1.78139513e+01, -1.90200165e+01,
                     2.75584151e+01,  2.17597215e+00, -2.63135333e+01, -3.50618702e+00,
                    -1.05887469e+01,  1.33862379e+01, -1.04764274e+01,  6.51953027e+01,
                    -6.41010276e+01,  3.30386652e+01,  1.15959259e+00, -1.05257549e+00,
                     1.81577115e+01;

    Eigen::MatrixXd DOFTP1A( 129, 1 );
    DOFTP1A <<  2.58575143e+00, -1.20703297e+01, -2.79245181e-01,  2.25934142e+00,
               -1.15932881e+01,  3.04401370e+00, -4.94690306e-01, -4.99177579e+00,
                1.65492736e-02, -5.24363920e-02, -3.89115194e+00,  1.62695400e+00,
                2.52138282e+00, -6.48239100e+00,  1.26695350e+00,  1.26145560e+00,
               -5.06487138e+00,  2.10406149e+00,  1.28357903e-01, -4.86886727e+00,
                6.25042041e-01, -1.28692867e-01, -3.71688562e+00,  1.11271700e+00,
               -6.79948119e+00, -1.12576499e+01, -2.68016433e+00, -4.35609910e+00,
               -1.07330925e+01,  1.19719123e+00, -2.55409962e+00, -6.74364627e+00,
               -7.28636947e-01, -1.47508290e+00, -5.19447680e+00,  6.31392970e-01,
                4.61831065e-01, -8.23428699e+00,  5.21720244e-01, -1.61367198e-01,
               -6.36822783e+00,  1.10844224e+00, -2.94744558e+00, -1.51681150e+01,
                1.41582350e+00, -1.87268630e+00, -1.38332592e+01,  3.90224563e+00,
                6.43779614e+00, -1.59807862e+01,  3.81662716e+00,  4.74277942e+00,
               -1.46934306e+01,  5.74874788e+00,  2.73378120e+00, -1.78933288e+01,
                5.29973573e+00, -7.72270266e-01, -8.09027946e+00,  1.82489403e+00,
                1.82246872e+00, -9.33246495e+00,  3.07784773e+00, -4.16392124e-01,
               -6.80442114e+00,  1.14527357e+00, -5.33735029e+00, -1.64343991e+01,
                3.63746343e+00, -2.12894714e+00, -8.95041292e+00,  1.82301351e+00,
                4.65626917e-01, -1.01926280e+01,  3.07591274e+00, -1.36129721e+00,
               -1.92226982e+01,  8.30275689e+00,  6.70983826e+00, -2.06816242e+01,
                9.96497856e+00, -1.95117630e-01, -3.20912648e+00, -3.02140693e+00,
               -3.99147053e-01, -1.83912174e+00, -2.82074509e+00,  7.03788389e+00,
                8.09106225e+00,  3.06122057e+00,  7.66511164e+00,  6.00515535e+00,
                8.39103759e+00, -3.80828682e+00,  9.98697945e-01, -4.50167194e+00,
               -1.04702823e+01,  1.03636777e+01,  2.57862396e+00,  4.95277593e+00,
               -3.74922286e+00, -4.27895224e+00, -9.61743660e+00,  8.31467065e+00,
                1.18168997e+01, -3.06891737e+00, -2.37031153e+00, -2.36009933e+00,
               -9.79052348e+00, -5.86673795e+00,  2.50933789e-01, -5.75398704e+00,
               -3.12988140e+00,  2.31645230e+00, -2.98467524e+00, -6.84250062e+00,
                8.73798302e+00,  5.43993037e-01, -6.57838333e+00, -8.76546754e-01,
               -4.13225918e+00,  4.26073160e+00, -5.14813913e+00,  1.67205446e+01,
               -1.47983784e+01,  9.65690819e+00,  2.17740120e+00,  8.71405012e-02,
                5.30835874e+00;

    if ( ( DotDotDOFTP1A - oc._test_DotDotDOF_tp1 ).norm( ) > 1e-2 * ( DotDotDOFTP1A.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( ( DOFTP1A - oc._test_DOF_tp1 ).norm( ) > 1e-2 * ( DOFTP1A.norm( ) + 1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    //Check the projected values of the ghost displacement
    const floatVector projectedGhostMicroDisplacementAnswer
        =
        {
              0.57409872,  -5.99059648,  -5.16487288,   0.67580229,
             -5.99320405,  -1.79300314,  -1.44547133,  -1.95639062,
             -4.44631554,  -1.21323214,  -2.51187232,  -1.22145146,
              1.2378831 ,  -5.55784005,  -3.03467742,   1.11118926,
             -4.62340866,  -0.19554806,  -0.31464079,  -1.64478451,
             -2.6333718 ,  -0.37313869,  -1.24955668,  -0.05957304,
             -5.48313625,  -0.97685356,  -6.75852263,  -4.89647693,
             -3.34426222,  -3.14832106,  -3.37070592,  -1.79922411,
             -4.31787507,  -2.89837418,  -2.94074004,  -1.28963773,
             -0.17848315,  -5.19257368,  -2.78177396,  -0.25607076,
             -4.92226921,  -0.18598817,  -3.74635013,  -5.28059651,
             -4.30853555,  -3.21445839,  -6.54075802,  -1.08172494,
              2.40039622, -10.18487233,  -2.39586573,   2.41471218,
             -9.12012943,   0.47609048
        };

    const floatVector projectedGhostMacroDisplacementAnswer
        =
        {
              2.03991416,  -4.44139433,   0.65063145,  -0.25566626,
             -5.15643919,   2.345488  ,  13.00243902,  12.88600197,
             -5.54187157,   4.03360134,   0.69049268,   2.98534445,
             -3.37023733,  -4.6460331 ,  -1.19770455, -10.56463672,
              9.14816409,  -3.45182511, -13.41171656,   9.75022929,
             -5.21748719,  -7.73027334,   1.91047105,   0.59754516,
             -1.72825907,  -7.49368327,   0.238362  ,  -0.70707224,
             -5.86420758,  -2.52872347, -15.22437972, -15.44552964,
             -9.26663334,   4.36752067,   0.96166206,   1.19374423,
              3.67073789,  -7.30268861,   2.04677339, -10.09092168,
              8.41808665,   3.27015416,  14.84239041, -18.60859054,
             -9.58911703,  -7.98434343,   2.1017912 ,   3.58635051,
              1.67330269,  -9.41617329,   2.73781623,  -2.11523401,
             -8.30742041,  -1.52278702,  22.97603963,  23.15610625,
             -7.63282171,  -5.36148545,  -9.21353538,   0.52935468,
             -4.00189592,  -8.49770197,   1.22715036,  -9.23516322,
              8.6092407 ,  -0.6200974 , -21.13909699,  18.65505205,
             -6.15000311,   2.34015371,  -7.51825517,   2.32144923,
             -1.7859679 , -11.08466576,   3.71156221,   1.96299978,
             -4.17738467,   1.3768405 , -24.03039625, -23.82897963,
             -3.34513811,  16.3265205 ,  12.48707889,   5.48973865,
              3.89153568, -12.0003216 ,   5.2302379 , -13.31800695,
             12.74388638,   0.48353066,  25.86170794, -28.32440287,
             -4.8164715 , -19.36387188,  14.19837873,   3.73194476
        };
    
    const floatVector *projectedGhostMicroDisplacementResult = oc.getProjectedGhostMicroDisplacement( );
    const floatVector *projectedGhostMacroDisplacementResult = oc.getProjectedGhostMacroDisplacement( );

    if ( !vectorTools::fuzzyEquals( *projectedGhostMicroDisplacementResult, projectedGhostMicroDisplacementAnswer, 1e-6, 1e-2 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( !vectorTools::fuzzyEquals( *projectedGhostMacroDisplacementResult, projectedGhostMacroDisplacementAnswer, 1e-6, 1e-1 ) ){
        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;
    }
    testNum++;

    //check the outputted DOF
    shared_ptr< XdmfReader > reader = XdmfReader::New( );
    shared_ptr< XdmfDomain > _readDomain = shared_dynamic_cast< XdmfDomain >( reader->read( "macroscale_dof.xdmf" ) );
    shared_ptr< XdmfUnstructuredGrid > _readGrid = _readDomain->getGridCollection( 0 )->getUnstructuredGrid( 0 );
    shared_ptr< XdmfAttribute > nodeIDs = _readGrid->getAttribute( "node_ids" );
    shared_ptr< XdmfAttribute > updatedDOF = _readGrid->getAttribute( "updated_DOF" );

    nodeIDs->read( );
    updatedDOF->read( );

    uIntVector macroNodeIdsAnswer
       = { 5,  9,  8, 11,  3,  1,  6, 15, 12,  2, 13, 14 };

    if ( nodeIDs->getSize( ) != macroNodeIdsAnswer.size( ) ){

        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;

    }

    for ( uIntType i = 0; i < macroNodeIdsAnswer.size( ); i++ ){

        if ( !vectorTools::fuzzyEquals( macroNodeIdsAnswer[ i ], nodeIDs->getValue< uIntType >( i ) ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    if ( updatedDOF->getSize( ) != oc.getUpdatedMacroDisplacementDOF( ).size( ) ){

        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;

    }
    testNum++;

    for ( uIntType i = 0; i < oc.getUpdatedMacroDisplacementDOF( ).size( ); i++ ){

        if ( !vectorTools::fuzzyEquals( oc.getUpdatedMacroDisplacementDOF( )[ i ], updatedDOF->getValue< floatType >( i ) ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;
    
    reader = XdmfReader::New( );
    _readDomain = shared_dynamic_cast< XdmfDomain >( reader->read( "microscale_dof.xdmf" ) );
    _readGrid = _readDomain->getGridCollection( 0 )->getUnstructuredGrid( 0 );
    nodeIDs = _readGrid->getAttribute( "node_ids" );
    updatedDOF = _readGrid->getAttribute( "updated_DOF" );

    nodeIDs->read( );
    updatedDOF->read( );

    uIntVector microNodeIdsAnswer
       = { 15, 31, 13, 26, 53, 21, 37, 48,  5, 10,  3,  4, 32, 33, 34, 28, 25,
           50, 43, 27,  1,  7, 30, 16, 22,  2, 46, 24, 39, 40, 57, 44, 58, 29,
           59, 11,  0, 20, 60, 47, 49, 17, 38, 14, 55 };

    if ( nodeIDs->getSize( ) != microNodeIdsAnswer.size( ) ){

        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;

    }

    for ( uIntType i = 0; i < microNodeIdsAnswer.size( ); i++ ){

        if ( !vectorTools::fuzzyEquals( microNodeIdsAnswer[ i ], nodeIDs->getValue< uIntType >( i ) ) ){

            results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
            return 1;

        }

    }
    testNum++;

    if ( updatedDOF->getSize( ) != oc.getUpdatedMicroDisplacementDOF( ).size( ) ){

        results << "test_overlapCoupling_processIncrement (test " + std::to_string( testNum ) + ") & False\n";
        return 1;

    }
    testNum++;

    for ( uIntType i = 0; i < oc.getUpdatedMicroDisplacementDOF( ).size( ); i++ ){

        if ( !vectorTools::fuzzyEquals( oc.getUpdatedMicroDisplacementDOF( )[ i ], updatedDOF->getValue< floatType >( i ) ) ){

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

int test_overlapCoupling_processIncrement_Arlequin( std::ofstream &results ){
    /*!
     * Test processing an increment using the Arlequin method
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    remove( "macroscale_out.xdmf" );
    remove( "macroscale_out.h5" );

    remove( "microscale_out.xdmf" );
    remove( "microscale_out.h5" );

    std::string testName = "overlapCoupling_processIncrement_Arlequin";
    uIntType testNum = 0;

    std::string filename = "testConfig_Arlequin.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_processIncrement_Arlequin & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_processIncrement_Arlequin & False\n";
        return 1;
    }

    error = oc.processIncrement( 1, 1 );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_processIncrement_Arlequin & False\n";
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

    const std::unordered_map< uIntType, floatType > arlequinMicroWeightingFactorsAnswer
        =
        {
            { 15, 0.50 },
            { 31, 0.75 },
            { 13, 0.50 },
            { 26, 0.75 },
            { 53, 0.50 },
            { 21, 0.75 },
            { 37, 0.50 },
            { 48, 0.75 },
            {  5, 0.50 },
            { 10, 0.75 },
            {  3, 0.50 },
            {  4, 0.75 },
            { 32, 0.50 },
            { 33, 0.75 },
            { 34, 0.50 },
            { 28, 0.75 },
            { 25, 0.50 },
            { 50, 0.75 },
            { 43, 1.00 },
            { 27, 1.00 },
            {  1, 1.00 },
            {  7, 1.00 },
            { 30, 1.00 },
            { 16, 1.00 },
            { 22, 1.00 },
            {  2, 1.00 },
            { 46, 1.00 },
            { 24, 0.00 },
            { 39, 0.25 },
            { 40, 0.00 },
            { 57, 0.25 },
            { 44, 0.00 },
            { 58, 0.25 },
            { 29, 0.00 },
            { 59, 0.25 },
            { 11, 0.00 },
            {  0, 0.25 },
            { 20, 0.00 },
            { 60, 0.25 },
            { 47, 0.00 },
            { 49, 0.25 },
            { 17, 0.00 },
            { 38, 0.25 },
            { 14, 0.00 },
            { 55, 0.25 }
        };

    const std::unordered_map< uIntType, floatType > *arlequinMicroWeightingFactorsResult = oc.getArlequinMicroWeightingFactors( );

    for( auto a = arlequinMicroWeightingFactorsAnswer.begin( ); a != arlequinMicroWeightingFactorsAnswer.end( ); a++ ){

        auto r = arlequinMicroWeightingFactorsResult->find( a->first );

        if ( r == arlequinMicroWeightingFactorsResult->end( ) ){

            results << testName + "(test " + std::to_string( testNum + 1 ) + ") & False\n";
            return 1;

        }

        if ( !vectorTools::fuzzyEquals( r->second, a->second, 1e-2, 1e-2 ) ){

            results << testName + "(test " + std::to_string( testNum + 2 ) + ") & False\n";
            return 1;

        }

    }
    testNum += 2;

    //Test the micromorphic mass matrix
    Eigen::MatrixXd MDX( 144, 1 );
    MDX <<  1.03221787, -0.92246333,  1.27200681, -0.48761604,  0.37037086,
            1.02800141, -0.06719396,  0.25023267, -0.50356931,  1.98702766,
           -1.33517121, -0.98618331,  1.67520885, -0.32002969,  0.31204602,
           -0.62283931,  1.88485606,  0.14203533,  1.8679691 ,  1.95733768,
            0.48543641,  1.22018616, -0.31776795,  0.21502873,  0.64847517,
           -0.44455271, -1.59968356, -1.38163596, -0.46210322,  0.28330518,
           -1.53446563,  0.52655343, -0.40943491, -1.28177123,  0.31722293,
           -1.20203397,  1.8838875 , -0.41063676,  0.98500843, -1.91833504,
            1.65884954,  0.65857992,  0.75360733, -1.04602099, -1.92897066,
           -0.23344673,  1.47255343,  1.69531131, -1.27639949, -1.16911678,
            1.59589692,  0.30165247,  0.56896295,  1.63962054, -0.89920825,
           -1.63118444,  0.02242142,  0.72002951, -1.01593182, -1.68371068,
           -1.54483822,  1.87738139,  1.61092459, -1.03754536,  1.10254212,
           -1.24979886,  0.42068153,  1.80115139,  0.58122646,  0.69030709,
           -0.30644823,  1.61575065,  1.76422882,  0.49346777,  1.9656111 ,
            1.09841017, -1.43482989, -1.24436724,  0.40463851, -1.94135932,
           -0.93089867, -0.83737637,  0.19009932, -0.73989132, -0.41730156,
            0.56866435, -0.24432187,  0.74592874, -0.29527219,  0.79572626,
           -1.16371926,  1.48413216,  0.16221102,  1.33208692, -0.57881395,
            0.75027824,  0.37030067, -0.27326494, -1.1946677 , -0.54679603,
            1.99430325, -1.57736132, -0.04530493, -1.92386778, -1.63076862,
            0.54026597,  1.14049803, -1.01269005, -0.12859646,  0.95372803,
            0.04832678, -0.76899178, -0.80743886, -0.33441147, -0.89729875,
           -1.57117314,  0.08745574, -1.80049292,  0.92430717,  0.85418227,
            0.91247804,  0.46793921,  1.71737995, -1.50876716,  0.05699081,
           -1.66527938,  0.93377481, -1.53093827,  0.81617988, -1.67771678,
            0.27021483,  1.96160743, -1.88528753,  1.40697224,  0.29062036,
            1.80679358, -0.47801091,  1.32833655,  0.21266734, -0.58380344,
            0.14261189,  0.01538911, -0.66847627,  0.69894409;

    Eigen::MatrixXd MDA( 144, 1 );
    MDA << 1.08489179e-02, -1.13861182e-02,  3.98692508e-02,  1.06097249e-03,
           1.88010929e-03,  2.21167184e-03, -2.86397165e-04, -4.64097754e-04,
          -5.67150313e-04, -3.29316751e-04, -5.64138010e-04, -5.10352321e-04,
           1.47813607e-02,  5.29863271e-03,  3.45748689e-02,  3.58082359e-04,
           7.05686896e-04,  7.77609290e-04,  4.85896197e-04,  9.98473185e-04,
           1.30778511e-03,  2.22201075e-05,  3.19662655e-05,  1.52122868e-04,
           2.93257542e-02,  3.52143598e-03,  2.23798973e-02, -9.93283414e-05,
          -1.87789594e-04, -3.13534150e-04, -4.82010813e-04, -7.95566194e-04,
          -9.72435761e-04, -1.15549063e-04, -1.87197163e-04, -2.20190601e-04,
           2.40849889e-02, -4.44084777e-03,  2.63668358e-02,  7.00832876e-04,
           1.21276710e-03,  1.38365198e-03, -7.59445653e-04, -1.29438086e-03,
          -1.60287425e-03,  3.07143521e-04,  5.58621363e-04,  7.71196193e-04,
          -1.18298271e-01,  2.85797547e-02,  2.02797288e-01,  3.67363688e-03,
           6.80774499e-03,  8.38041496e-03, -4.35798160e-03, -8.34823423e-03,
          -1.07580395e-02, -1.76534630e-03, -3.06261660e-03, -3.15450490e-03,
          -7.52329741e-02,  1.53249940e-01,  2.60984102e-01, -2.29899017e-03,
          -3.60452083e-03, -4.59441176e-03, -1.02599409e-03, -2.05514754e-03,
          -2.61051612e-03,  2.08499778e-03,  3.35674414e-03,  4.09184613e-03,
           5.92431469e-02,  1.39702483e-01,  2.63982000e-01, -4.31218420e-03,
          -7.49155922e-03, -8.98555220e-03, -3.71534663e-03, -7.04933329e-03,
          -8.71795523e-03,  1.27348919e-03,  1.88649770e-03,  1.97252326e-03,
          -4.77133248e-02,  8.93492564e-02,  1.67641969e-01,  2.24610053e-03,
           3.96713993e-03,  5.20900316e-03, -3.27966713e-03, -6.17318798e-03,
          -8.02407771e-03,  7.71352997e-04,  1.13958931e-03,  1.73966242e-03,
          -7.68778268e-02,  7.12389826e-02,  3.83230834e-02, -2.35151526e-05,
           6.16467448e-04,  8.65037513e-04, -7.48452232e-03, -1.41644118e-02,
          -1.77166101e-02,  7.92633619e-04,  1.61093729e-03,  1.98144227e-03,
          -3.63230377e-02,  1.37528206e-01,  1.17214862e-01, -4.31931644e-03,
          -7.35031740e-03, -9.32231083e-03, -4.95122959e-03, -9.76493068e-03,
          -1.23175796e-02,  3.56386417e-03,  5.98995118e-03,  6.57098690e-03,
           1.27596196e-02,  1.41852943e-01,  1.87456864e-01, -4.97049416e-03,
          -8.63435902e-03, -1.06883844e-02, -3.51502106e-03, -7.21540259e-03,
          -8.90435329e-03,  3.97121612e-03,  6.32700458e-03,  6.76549938e-03,
          -9.81698600e-02,  1.32295264e-01,  9.03577259e-02,  7.41369217e-04,
           1.53759949e-03,  2.29039869e-03, -4.09621493e-03, -7.99445372e-03,
          -9.98475098e-03,  2.15206617e-03,  3.42194005e-03,  3.93265861e-03;

    if( ( ( *oc.getMD( ) ) * MDX - MDA ).norm() > 1e-6*( MDA.norm() + 1 ) ){
        std::cerr << "  Failure in test of the micromorphic mass matrix\n";
        results << testName + "(test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;
    }
    testNum++;

    floatVector FQAnswer = {  2.74285316e-01, -2.83800284e+00, -2.84112402e-01,  2.14794384e-01,
                             -1.65375089e+00,  1.70748023e-01, -5.81181752e-01, -2.61824581e+00,
                             -3.67177391e-01, -2.12081969e-01, -1.54409257e+00,  1.29298760e-01,
                              8.45145845e-01, -2.97680154e+00,  4.26257107e-01,  4.99652644e-01,
                             -1.72301116e+00,  5.25220984e-01, -1.03212230e-02, -2.75704452e+00,
                              3.43192118e-01,  7.27762913e-02, -1.61335285e+00,  4.83771721e-01,
                             -1.43664882e+00, -2.39848879e+00, -4.50242380e-01, -6.38958321e-01,
                             -1.43443426e+00,  8.78494969e-02, -8.65788291e-01, -2.53728749e+00,
                              2.60127129e-01, -3.54100061e-01, -1.50369453e+00,  4.42322458e-01,
                              5.60539306e-01, -2.89584322e+00,  1.05356163e+00,  3.57634551e-01,
                             -1.68261312e+00,  8.38244682e-01, -2.94927762e-01, -2.67608619e+00,
                              9.70496637e-01, -6.92418012e-02, -1.57295481e+00,  7.96795419e-01,
                              1.41600637e+00, -3.11560024e+00,  1.13662661e+00,  7.84510904e-01,
                             -1.79227144e+00,  8.79693945e-01,  1.17558612e-03, -7.59575922e-03,
                              1.94083360e-03, -5.38776738e-04, -7.15536438e-03,  1.77437069e-03,
                              2.31959520e-03, -7.87391293e-03,  3.36441978e-03,  6.05232338e-04,
                             -7.43351809e-03,  3.19795688e-03, -2.25313960e-03, -6.71496954e-03,
                              1.60790779e-03, -1.10913052e-03, -6.99312325e-03,  3.03149398e-03,
                              1.74924141e-03, -7.71167181e-03,  4.62154307e-03,  3.48785530e-05,
                             -7.27127696e-03,  4.45508017e-03,  3.46360427e-03, -8.15206665e-03,
                              4.78800597e-03, -7.61699343e-02, -3.77522265e+00, -3.07648110e+00,
                              1.77297210e-01, -3.54516009e+00, -1.36652211e+00, -1.78881843e+00,
                             -3.33526821e+00, -3.24277754e+00, -1.10676057e+00, -3.21530436e+00,
                             -1.49120282e+00,  1.06669513e+00, -4.05309821e+00, -1.65431849e+00,
                              1.03416001e+00, -3.75349722e+00, -3.00256052e-01, -6.45953365e-01,
                             -3.61314377e+00, -1.82061493e+00, -2.49897775e-01, -3.42364149e+00,
                             -4.24936767e-01, -3.50146693e+00, -2.89531376e+00, -3.40907398e+00,
                             -2.39081836e+00, -2.88544862e+00, -1.61588354e+00, -2.35860186e+00,
                             -3.17318932e+00, -1.98691137e+00, -1.53395556e+00, -3.09378575e+00,
                             -5.49617482e-01,  4.96911702e-01, -3.89101932e+00, -3.98452331e-01,
                              6.06965023e-01, -3.63197862e+00,  6.41329288e-01, -1.21573680e+00,
                             -3.45106488e+00, -5.64748772e-01, -6.77092760e-01, -3.30212288e+00,
                              5.16648573e-01,  2.20956020e+00, -4.33097377e+00, -2.32155890e-01,
                              1.89102281e+00, -3.96183436e+00,  7.66010003e-01 };

    if ( !vectorTools::fuzzyEquals( *oc.getFQ( ), FQAnswer ) ){
        std::cerr << "  Failure in test of the micro force vector\n";
        results << testName + "(test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;
    }
    testNum++;

    floatVector FDAnswer = {  3.23858499e-03,  2.09589199e-03, -3.42551959e-03, -3.40171968e-03,
                              3.12532604e-04, -2.87943629e-03,  2.39812580e-03, -2.91374371e-03,
                             -4.16822570e-04,  5.35491008e-05, -3.96354509e-03, -2.76773985e-03,
                              4.56237215e-03,  1.64920222e-03, -2.36794154e-03, -2.07443744e-03,
                             -3.16741627e-03, -4.24670327e-03,  2.36664420e-03, -6.30484437e-03,
                              1.83546267e-03, -1.60978102e-03, -4.99350093e-03, -4.14636879e-03,
                              8.70950927e-03,  2.61210642e-03, -2.34107920e-03, -2.03278932e-04,
                             -6.02159803e-03, -5.88398646e-04,  3.68715026e-03, -1.84694936e-03,
                             -2.51209428e-03, -4.53739709e-03, -7.76773475e-03, -6.08135464e-03,
                              7.38572212e-03,  3.05879619e-03, -3.39865725e-03, -1.53056117e-03,
                             -2.54164916e-03,  7.78868338e-04,  3.71863185e-03,  1.54415130e-03,
                             -4.76437952e-03, -2.87406697e-03, -6.73777891e-03, -4.70272571e-03,
                              1.17803670e-01,  1.03645637e+00, -3.75514827e+00, -2.61559903e+00,
                              1.09281113e+00, -1.65637048e+00,  1.49324606e+00, -1.19123623e+00,
                              1.89878202e+00,  1.56943401e+00, -3.33115710e+00, -2.03239898e+00,
                              7.79697248e-01,  8.13111484e-01, -3.22635925e+00, -1.94335947e+00,
                             -6.24491430e-01, -2.31939443e+00,  1.50411732e+00, -2.84236792e+00,
                              3.07419831e+00,  7.21786772e-01, -3.89924128e+00, -2.79304431e+00,
                              2.85326581e+00,  1.29456359e+00, -3.21292808e+00, -1.03859819e+00,
                             -2.10093603e+00, -5.58563342e-01,  2.18687582e+00, -5.78294023e-01,
                              9.40840646e-01, -6.74153864e-01, -5.16200793e+00, -3.60197241e+00,
                              2.19137223e+00,  1.51790847e+00, -3.74171710e+00, -1.71083775e+00,
                             -3.83633468e-01,  1.04460614e-01,  2.17600456e+00,  1.07283767e+00,
                             -2.34575649e-01,  1.73493375e-01, -4.59392375e+00, -2.84132708e+00,
                             -2.76460294e+00,  2.04788356e+00, -1.15834784e+01, -6.55090245e+00,
                              4.91390770e+00, -2.64290383e+00,  2.78468843e+00, -3.28819690e+00,
                              6.18444966e+00,  5.92737766e+00, -9.85263495e+00, -6.00962178e+00,
                             -1.44213958e+00,  1.60164048e+00, -1.05269579e+01, -5.20764660e+00,
                              1.48305678e+00, -3.96733517e+00,  2.80678435e+00, -6.58653187e+00,
                              8.53362602e+00,  4.23355318e+00, -1.09884158e+01, -7.53039669e+00,
                              2.70085041e+00,  2.56358177e+00, -1.05001225e+01, -3.40036800e+00,
                             -1.46757525e+00, -4.50157755e-01,  4.17125309e+00, -2.06241705e+00,
                              4.27174721e+00,  1.44542050e+00, -1.35096706e+01, -9.14439978e+00,
                              1.37838704e+00,  3.00982485e+00, -1.15566429e+01, -4.74362385e+00,
                              1.96327566e+00,  8.74273580e-01,  4.14915716e+00,  1.23591792e+00,
                              1.92257085e+00,  3.13924498e+00, -1.23738898e+01, -7.62362488e+00 };

    if ( !vectorTools::fuzzyEquals( *oc.getFD( ), FDAnswer ) ){
        std::cerr << "  Failure in test of the macro force vector\n";
        results << testName + "(test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;
    }
    testNum++;

    Eigen::VectorXd DotDotD_tp1Answer( 144 );
    DotDotD_tp1Answer <<  6.47725051e+01,  3.96167211e+01, -1.52219309e+02,  3.49828544e+03,
                         -9.56889125e+03,  6.41508339e+03,  7.38300065e+03, -1.69679779e+04,
                          1.04894435e+04, -4.75042811e+03,  1.17360164e+04, -7.15883929e+03,
                         -3.04104225e+02, -7.02518555e+02,  3.26087481e+02,  4.15753775e+03,
                          1.69849567e+04, -1.60211921e+04,  9.00104701e+03,  3.06945396e+04,
                         -2.93274706e+04, -6.57258930e+03, -2.06237217e+04,  1.99529539e+04,
                          2.99731585e+02,  6.51154751e+02, -4.38410746e+02, -1.77195145e+04,
                          1.50592160e+04, -1.68208571e+03, -3.32152402e+04,  2.72394312e+04,
                         -2.80367616e+03,  2.28495408e+04, -1.86452065e+04,  2.19386989e+03,
                         -3.38042111e+02, -5.19155629e+02,  4.08688594e+02, -1.64730383e+04,
                         -7.24392811e+03,  1.11730902e+04, -2.98321664e+04, -1.42437346e+04,
                          2.14193781e+04,  2.05987523e+04,  9.47333329e+03, -1.45552959e+04,
                          3.60212753e+02,  6.40479061e+02, -4.31737390e+02,  1.74608371e+04,
                         -9.38640882e+03, -9.83367542e+02,  3.14789317e+04, -1.79860458e+04,
                         -4.95435652e+02, -2.12955837e+04,  1.19920434e+04,  5.44589333e+02,
                         -2.04374618e+02, -2.01586974e+02,  2.15687473e+02,  1.74180563e+04,
                          4.66347486e+03, -1.05841953e+04,  3.09981330e+04,  8.51417718e+03,
                         -1.85873381e+04, -2.08479961e+04, -5.56113771e+03,  1.23791532e+04,
                          8.10313459e+01,  7.66301810e+01,  5.59501232e+00, -5.51917388e+03,
                          4.66394940e+03, -1.19827304e+03, -1.00904431e+04,  9.70849377e+03,
                         -3.45809125e+03,  6.11826938e+03, -6.57375168e+03,  2.71746601e+03,
                         -1.45193217e+02, -3.30912252e+02,  1.39986437e+02, -6.79163366e+03,
                         -1.08093764e+04,  1.08980493e+04, -1.30379307e+04, -2.07051379e+04,
                          2.10166556e+04,  8.25392643e+03,  1.39711701e+04, -1.39221648e+04,
                         -3.94689916e+02, -7.13387345e+02,  4.30127612e+02, -2.11086519e+04,
                          1.74522179e+04, -5.00813769e+03, -3.58654387e+04,  2.80652006e+04,
                         -6.27723312e+03,  2.47842579e+04, -1.93914760e+04,  4.55197000e+03,
                          5.01028980e+02,  2.09403793e+02, -8.95200046e+02, -2.28065596e+04,
                         -1.28517095e+04,  1.82010614e+04, -3.94992681e+04, -2.56453215e+04,
                          3.55082492e+04,  2.69973816e+04,  1.38378617e+04, -2.21241427e+04,
                         -7.51933587e+01, -3.74725013e+02,  2.86743108e+02,  1.49759639e+04,
                         -1.19782365e+04,  1.98123118e+03,  2.91732809e+04, -2.09496449e+04,
                          2.00217151e+03, -1.69856368e+04,  1.11616136e+04, -1.36687532e+03,
                         -1.40281659e+01,  7.61912447e+02,  8.36658984e+01,  1.36710990e+04,
                          1.45369219e+04, -1.47592681e+04,  2.78478839e+04,  2.62533716e+04,
                         -2.86241850e+04, -1.59574336e+04, -1.77479864e+04,  1.78663100e+04;

    if ( ( DotDotD_tp1Answer - oc._test_DotDotD_tp1 ).norm( ) > 1e-6 * ( DotDotD_tp1Answer.norm( ) + 1 ) ){
        std::cerr << "  Failure in test of the macro DOF acceleration vector\n";
        results << testName + "(test " + std::to_string( testNum + 1 ) + ") & False\n";
    }
    testNum++;

    Eigen::VectorXd DotDotQ_tp1Answer( 135 );
    DotDotQ_tp1Answer <<  -5.42317374e+02, -8.87978563e+02,  6.02246559e+02,  8.05196770e+00,
                           2.28246530e+01, -4.07172779e+01,  1.36472442e+02,  3.25534813e+02,
                          -2.02051872e+02, -4.90427924e+02, -7.45249477e+02,  5.85895747e+02,
                          -3.36031350e+02, -6.60480560e+02,  4.41804920e+02,  6.30752312e+02,
                           8.17339522e+02, -6.70130094e+02,  1.95903117e+01, -2.21992399e+01,
                           1.60806144e+01,  1.49025841e+01, -1.13517759e+02,  4.44358064e+01,
                           8.43488439e+02,  1.54519860e+03, -1.02522128e+03, -1.08591925e+03,
                          -1.79547983e+03,  1.37933838e+03,  4.73660711e+02,  7.82606111e+02,
                          -5.43107209e+02, -5.58626056e+02, -1.06043467e+03,  7.42099698e+02,
                          -2.54948926e+02, -4.79160907e+02,  3.04540740e+02,  3.01854915e+02,
                           6.26197279e+02, -5.11283320e+02,  9.29086054e+00,  9.67907405e+01,
                          -1.31516928e+02, -8.06697301e+01,  1.15002832e+02, -1.41340196e+02,
                          -3.50517419e+02, -7.28214902e+02,  4.92542350e+02,  8.66033172e+02,
                           1.38742868e+03, -1.08186007e+03,  7.04828408e+01,  4.19681210e+02,
                          -2.85847815e+02,  1.36245235e+02,  8.24758933e+02, -2.55924928e+02,
                          -1.78243882e+01, -4.90809831e+02,  4.39507675e+01,  7.25718102e+01,
                           4.52319628e+01, -4.84770022e+01,  1.50721584e+01,  7.45679705e+02,
                           5.75682850e+01,  1.46127366e+02,  3.91121159e+02, -5.30623538e+01,
                          -2.32059527e+02, -4.65580640e+02,  7.18216422e+01,  2.48383906e+02,
                           6.84974234e+02, -5.06992581e+02, -5.59248775e+02, -1.51228376e+03,
                           5.42750375e+02, -3.32392474e+01, -1.60212358e+02, -7.95895051e+00,
                           7.88241643e+01,  5.62859010e+01, -1.03400148e+02,  3.24467886e+01,
                           7.01443867e+01, -1.22332823e+02, -1.44926046e+02, -2.91118809e+02,
                           1.62370790e+02, -1.75714124e+02, -3.93604597e+02,  2.04671394e+02,
                           2.43661349e+02,  4.00889431e+02, -2.97660713e+02, -3.99560174e+00,
                          -3.82273454e+00, -2.30121401e+01,  1.00084661e+01,  1.45722730e+00,
                          -1.63929252e+00,  1.48771439e+02,  3.91894133e+02, -4.00250065e+02,
                          -3.11587068e+02, -5.46598859e+02,  3.02454620e+02,  1.49581147e+02,
                           2.96491597e+02, -2.63682726e+02, -1.85007975e+02, -3.55201375e+02,
                           2.15248302e+02, -4.60965414e+01, -1.93877965e+02,  8.01653671e+01,
                           1.02355818e+02,  1.53135458e+02, -1.22708780e+02, -3.22781506e+01,
                          -1.43509153e+02,  2.88016631e+01, -2.09367969e+02, -3.75512420e+02,
                           2.10028132e+02, -1.46837095e+02, -5.14574842e+02,  2.69098337e+02,
                           4.34263298e+02,  6.75405879e+02, -4.88026236e+02;

    if ( ( DotDotQ_tp1Answer - oc._test_DotDotQ_tp1 ).norm( ) > 1e-6 * ( DotDotQ_tp1Answer.norm( ) + 1 ) ){
        std::cerr << "  Failure in test of the micro DOF acceleration vector\n";
        results << testName + "(test " + std::to_string( testNum + 1 ) + ") & False\n";
    }
    testNum++;

    const floatVector projectedGhostMacroDisplacementAnswer
        = {  9.00531882e+01,  1.60119765e+02, -1.07934347e+02,  4.36450048e+03,
            -2.34697440e+03, -2.46512192e+02,  7.87237621e+03, -4.49425561e+03,
            -1.22815410e+02, -5.32349948e+03,  2.99904211e+03,  1.35445863e+02,
            -5.10936546e+01, -5.03967436e+01,  5.39218681e+01,  4.35228829e+03,
             1.16660514e+03, -2.64688712e+03,  7.75109630e+03,  2.13072760e+03,
            -4.64508785e+03, -5.21202127e+03, -1.39169095e+03,  3.09516629e+03,
             2.02578365e+01,  1.91575453e+01,  1.39875308e+00, -1.38361028e+03,
             1.16809102e+03, -3.01788778e+02, -2.52156683e+03,  2.42938648e+03,
            -8.61905697e+02,  1.53016142e+03, -1.64455039e+03,  6.80406856e+02,
            -3.62983043e+01, -8.27280631e+01,  3.49966093e+01, -1.70020822e+03,
            -2.70134905e+03,  2.72245980e+03, -3.25735848e+03, -5.17394888e+03,
             5.25607784e+03,  2.06449437e+03,  3.49411783e+03, -3.48058029e+03,
            -9.86724790e+01, -1.78346836e+02,  1.07531903e+02, -5.27868650e+03,
             4.36276316e+03, -1.25222821e+03, -8.96201392e+03,  7.01966470e+03,
            -1.56774807e+03,  6.19558618e+03, -4.84586272e+03,  1.13648300e+03,
             1.25257245e+02,  5.23509481e+01, -2.23800012e+02, -5.70468042e+03,
            -3.21211007e+03,  4.54990356e+03, -9.87155150e+03, -6.40803835e+03,
             8.87932568e+03,  6.74844842e+03,  3.45903392e+03, -5.53146572e+03,
            -1.87983397e+01, -9.36812533e+01,  7.16857770e+01,  3.73935943e+03,
            -2.99237455e+03,  4.93563794e+02,  7.29606665e+03, -5.23403948e+03,
             5.03676700e+02, -4.24668986e+03,  2.79026594e+03, -3.41486500e+02,
            -3.50704148e+00,  1.90478112e+02,  2.09164746e+01,  3.41466020e+03,
             3.63530641e+03, -3.69139303e+03,  6.96579765e+03,  6.56678719e+03,
            -7.15361559e+03, -3.98922037e+03, -4.43469627e+03,  4.46573039e+03 };

    const floatVector projectedGhostMicroDisplacementAnswer
        = {  -8.30981185,  -40.0530895 ,   -1.98973763,   19.70604107,
             14.07147525,  -25.85003704,    8.11169714,   17.53609668,
            -30.58320567,  -36.23151157,  -72.7797022 ,   40.59269753,
            -43.92853106,  -98.40114933,   51.16784844,   60.91533714,
            100.22235775,  -74.41517823,   -0.99890044,   -0.95568364,
             -5.75303502,    2.50211652,    0.36430683,   -0.40982313,
             37.19285967,   97.97353321, -100.06251628,  -77.89676692,
           -136.64971468,   75.61365503,   37.39528666,   74.12289921,
            -65.92068155,  -46.25199379,  -88.8003438 ,   53.81207544,
            -11.52413535,  -48.46949126,   20.04134176,   25.58895446,
             38.28386461,  -30.677195  ,   -8.06953765,  -35.87728815,
              7.20041577,  -52.34199235,  -93.87810502,   52.50703309,
            -36.70927366, -128.64371042,   67.27458419,  108.56582459,
            168.85146987, -122.00655909 };

    const floatVector updatedFreeMacroDisplacementDOFAnswer
        = {    16.19312628,     9.90418028,   -38.05482716,   874.67731082,
            -2392.67588963,  1602.62402521,  1846.69096728, -4240.84733416,
             2622.88766205, -1186.33584674,  2934.06034855, -1789.60326676,
              -76.02605634,  -175.62963868,    81.52187017,  1037.97337869,
             4246.8947068 , -4006.61283849,  2250.12231719,  7674.70950782,
            -7330.63769453, -1642.29483288, -5158.31196797,  4989.42447017,
               74.93289631,   162.7886878 ,  -109.60268638, -4432.88069876,
             3766.82679763,  -423.21845776, -8304.46856992,  6811.01212449,
             -698.81862998,  5713.85401173, -4663.38911845,   550.31585138,
              -84.51052768,  -129.7889072 ,   102.17214839, -4119.74465185,
            -1810.0678566 ,  2790.74351923, -7457.61989196, -3559.70677784,
             5356.24176504,  5151.57558245,  2368.68360682, -3638.05504579 };

    const floatVector updatedFreeMicroDisplacementDOFAnswer
        = { -135.57934359, -221.99464073,  150.56163986,    2.01299192,
               5.70616325,  -10.17931948,   34.1181106 ,   81.3837033 ,
             -50.51296811, -122.60698097, -186.31236929,  146.47393687,
             -84.00783759, -165.12014012,  110.45122995,  157.68807789,
             204.33488041, -167.5325235 ,    4.89757791,   -5.54980998,
               4.02015361,    3.72564603,  -28.37943972,   11.10895159,
             210.87210984,  386.29964943, -256.30532076, -271.47981315,
            -448.86995839,  344.83459508,  118.4151777 ,  195.65152778,
            -135.77680227, -139.65651395, -265.1086687 ,  185.5249244 ,
             -63.7372315 , -119.79022676,   76.13518493,   75.46372876,
             156.54931982, -127.82082997,    2.32271513,   24.19768513,
             -32.87923212,  -20.16743253,   28.75070811,  -35.33504911,
             -87.62935463, -182.05372553,  123.1355876 ,  216.50829307,
             346.85717037, -270.46501725,   17.62071021,  104.92030246,
             -71.46195368,   34.06130872,  206.18973321,  -63.98123197,
              -4.45609704, -122.70245777,   10.98769187,   18.14295254,
              11.30799069,  -12.11925054,    3.76803961,  186.41992623,
              14.39207125,   36.53184147,   97.78028966,  -13.26558845,
             -58.01488172, -116.39516007,   17.95541056,   62.09597653,
             171.24355844, -126.74814524, -139.81219366, -378.07093983,
             135.68759386 };

    const floatVector* projectedGhostMacroDisplacementResult = oc.getProjectedGhostMacroDisplacement( );
    const floatVector* projectedGhostMicroDisplacementResult = oc.getProjectedGhostMicroDisplacement( );

    const floatVector* updatedFreeMacroDisplacementDOFResult = oc.getUpdatedFreeMacroDispDOFValues( );
    const floatVector* updatedFreeMicroDisplacementDOFResult = oc.getUpdatedFreeMicroDispDOFValues( );

    if ( !vectorTools::fuzzyEquals( *updatedFreeMacroDisplacementDOFResult, updatedFreeMacroDisplacementDOFAnswer ) ){
        std::cerr << "  Failure in test of the updated macro DOF vector\n";
        results << testName + "(test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( !vectorTools::fuzzyEquals( *updatedFreeMicroDisplacementDOFResult, updatedFreeMicroDisplacementDOFAnswer ) ){
        std::cerr << "  Failure in test of the updated micro DOF vector\n";
        results << testName + "(test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( !vectorTools::fuzzyEquals( *projectedGhostMacroDisplacementResult, projectedGhostMacroDisplacementAnswer ) ){
        std::cerr << "  Failure in test of the updated macro DOF vector\n";
        results << testName + "(test " + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;
    }
    testNum++;

    if ( !vectorTools::fuzzyEquals( *projectedGhostMicroDisplacementResult, projectedGhostMicroDisplacementAnswer ) ){
        std::cerr << "  Failure in test of the updated micro DOF vector\n";
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

        std::cerr << "  Failure in test of the homogenized volumes\n";
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

        std::cerr << "  Failure in test of the homogenized densities\n";
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

        std::cerr << "  Failure in test of the homogenized symmetric micro stress\n";
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

        std::cerr << "  Failure in the homogenized centers of mass\n";
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

        std::cerr << "  Failure in the homogenized body forces\n";
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

        std::cerr << "  Failure in the homogenized accelerations\n";
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

        std::cerr << "  Failure in the homogenized micro inertias\n";
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

        std::cerr << "  Failure in the homogenized body force couples\n";
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

        std::cerr << "  Failure in the homogenized micro spin inertias\n";
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

        std::cerr << "  Failure in the homogenized surface areas\n";
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

        std::cerr << "  Failure in the cell domain macro surfaces\n";
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

        std::cerr << "  Failure in the homogenized surface region areas\n";
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

        std::cerr << "  Failure in the homogenized surface region centers of mass\n";
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

        std::cerr << "  Failure in the homogenized surface region projected local centers of mass\n";
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

        std::cerr << "  Failure in homogenized surface region projected centers of mass\n";
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

        std::cerr << "  Failure in homogenization of the surface region tractions\n";
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

        std::cerr << "  Failure in the homogenization of the surface region couples\n";
        return 1;

    }

    floatVector elementNodalVolumesAnswer( 8, .125 );

    for ( auto c = oc._test_elementNodalVolumes.begin( ); c != oc._test_elementNodalVolumes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( elementNodalVolumesAnswer, c->second, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in element nodal volumes\n";
            return 1;

        }

    }
    testNum++;

    floatVector volumeAtNodesAnswer( 8, .125 );

    for ( auto c = oc._test_volumeAtNodes.begin( ); c != oc._test_volumeAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( volumeAtNodesAnswer, c->second, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in volume at nodes\n";
            return 1;

        }

    }
    testNum++;

    floatVector densityAtNodesAnswer( 8, 2 );

    for ( auto c = oc._test_densityAtNodes.begin( ); c != oc._test_densityAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( densityAtNodesAnswer, c->second, 1e-6, 1e-2 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in density at nodes\n";
            return 1;

        }

    }
    testNum++;

    floatVector tmp = { -1, -2, -3 };
    floatMatrix bodyForceAtNodesAnswer( 8, tmp );

    for ( auto c = oc._test_bodyForceAtNodes.begin( ); c != oc._test_bodyForceAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( bodyForceAtNodesAnswer, c->second, 1e-6, 1e-6 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in body force at nodes\n";
            return 1;

        }

    }
    testNum++;

    tmp = { 0., 0., 3e-3 };
    floatMatrix accelerationAtNodesAnswer( 8, tmp );

    for ( auto c = oc._test_accelerationAtNodes.begin( ); c != oc._test_accelerationAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( accelerationAtNodesAnswer, c->second, 1e-6, 1e-6 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in acceleration at nodes\n";
            return 1;

        }

    }
    testNum++;

    tmp = { 0.0625, 0, 0, 0, 0.0625, 0, 0, 0, 0.0625 };

    floatMatrix microInertiaAtNodesAnswer( 8, tmp );

    for ( auto c = oc._test_microInertiaAtNodes.begin( ); c != oc._test_microInertiaAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( microInertiaAtNodesAnswer, c->second, 1e-6, 1e-6 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in micro inertas at nodes\n";
            return 1;

        }

    }
    testNum++;

    tmp = floatVector( 9, 0 );
    floatMatrix bodyCoupleAtNodesAnswer( 8, tmp );

    for ( auto c = oc._test_bodyCoupleAtNodes.begin( ); c != oc._test_bodyCoupleAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( bodyCoupleAtNodesAnswer, c->second, 1e-6, 1e-6 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in body couple at nodes\n";
            return 1;

        }

    }
    testNum++;

    tmp = floatVector( 9, 0 );
    floatMatrix microSpinInertiaAtNodesAnswer( 8, tmp );

    for ( auto c = oc._test_microSpinInertiaAtNodes.begin( ); c != oc._test_microSpinInertiaAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( microSpinInertiaAtNodesAnswer, c->second, 1e-6, 1e-6 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in micro spin inerta at nodes\n";
            return 1;

        }

    }
    testNum++;

    tmp = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    floatMatrix symmetricMicroStressAtNodesAnswer( 8, tmp );

    for ( auto c = oc._test_symmetricMicroStressAtNodes.begin( ); c != oc._test_symmetricMicroStressAtNodes.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( symmetricMicroStressAtNodesAnswer, c->second, 1e-6, 1e-6 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in symmetric micro stress at nodes\n";
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
        = { -0.125, -0.5  , -0.875, -0.25 , -0.625, -1.   , -0.375, -0.75 ,
            -1.125, -0.125, -0.5  , -0.875, -0.25 , -0.625, -1.   , -0.375,
            -0.75 , -1.125, -0.125, -0.5  , -0.875, -0.25 , -0.625, -1.   ,
            -0.375, -0.75 , -1.125, -0.125, -0.5  , -0.875, -0.25 , -0.625,
            -1.   , -0.375, -0.75 , -1.125, -0.125, -0.5  , -0.875, -0.25 ,
            -0.625, -1.   , -0.375, -0.75 , -1.125, -0.125, -0.5  , -0.875,
            -0.25 , -0.625, -1.   , -0.375, -0.75 , -1.125, -0.125, -0.5  ,
            -0.875, -0.25 , -0.625, -1.   , -0.375, -0.75 , -1.125, -0.125,
            -0.5  , -0.875, -0.25 , -0.625, -1.   , -0.375, -0.75 , -1.125 };

    for ( auto c = oc.getExternalForcesAtNodes( )->begin( ); c != oc.getExternalForcesAtNodes( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( extForceVectorsAnswer, c->second, 1e-6, 1e-1 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in external forces at nodes answer\n";
            return 1;

        }

    }
    testNum++;

    for ( auto c = oc.getExternalCouplesAtNodes( )->begin( ); c != oc.getExternalCouplesAtNodes( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( floatVector( coupleVectorsAnswer.size( ), 0 ), c->second, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in external couples at nodes answer\n";
            return 1;

        }

    }
    testNum++;

    for ( auto c = oc._test_cellLinearMomentumRHS.begin( ); c != oc._test_cellLinearMomentumRHS.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( forceVectorsAnswer, c->second, 1e-6, 1e-1 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in linear momentum RHS\n";
            return 1;

        }

    }
    testNum++;

    for ( auto c = oc._test_cellFirstMomentRHS.begin( ); c != oc._test_cellFirstMomentRHS.end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( coupleVectorsAnswer, c->second, 1e-6, 1e-1 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in first moment RHS\n";
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

    for ( auto LHS = oc._test_stressProjectionLHS.begin( ); LHS != oc._test_stressProjectionLHS.end( ); LHS++ ){

        if ( ( LHS->second * LHSX - LHSA ).norm( ) > 1e-6 * ( LHSA.norm( ) + 1 ) ){

            std::cout << LHS->second * LHSX << "\n";

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in stress projection LHS\n";
            return 1;

        }

    }
    testNum++;

    floatVector cauchyQptAnswer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    floatVector cauchyAnswer = vectorTools::appendVectors( floatMatrix( 8, cauchyQptAnswer ) );
    for ( auto c = oc.getQuadraturePointCauchyStress( )->begin( ); c != oc.getQuadraturePointCauchyStress( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, cauchyAnswer, 1e-6, 1e-1 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in cauchy qpt\n";
            return 1;

        }

    }
    testNum++;

    floatVector higherOrderAnswer( 27 * 8, 0 );
    for ( auto c = oc.getQuadraturePointHigherOrderStress( )->begin( ); c != oc.getQuadraturePointHigherOrderStress( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, higherOrderAnswer, 1e-6, 1e-1 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in higher order qpt\n";
            return 1;

        }

    }
    testNum++;

    floatVector densityAnswer( 8, 2 );
    for ( auto c = oc.getQuadraturePointDensities( )->begin( ); c != oc.getQuadraturePointDensities( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, densityAnswer, 1e-6, 1e-2 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in density qpt\n";
            return 1;

        }

    }
    testNum++;

    floatVector bodyForceQptAnswer = { -1, -2, -3 };
    floatVector bodyForceAnswer = vectorTools::appendVectors( floatMatrix( 8, bodyForceQptAnswer ) );
    for ( auto c = oc.getQuadraturePointBodyForce( )->begin( ); c != oc.getQuadraturePointBodyForce( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, bodyForceAnswer, 1e-6, 1e-2 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in body force qt\n";
            return 1;

        }

    }
    testNum++;

    floatVector accelerationQptAnswer = { 0., 0., 3e-3 };
    floatVector accelerationAnswer = vectorTools::appendVectors( floatMatrix( 8, accelerationQptAnswer ) );
    for ( auto c = oc.getQuadraturePointAccelerations( )->begin( ); c != oc.getQuadraturePointAccelerations( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, accelerationAnswer, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in acceleration qpt\n";
            return 1;

        }

    }
    testNum++;

    floatVector microInertiaQptAnswer = { 0.0625, 0., 0., 0., 0.0625, 0., 0., 0., 0.0625 };
    floatVector microInertiasAnswer = vectorTools::appendVectors( floatMatrix( 8, microInertiaQptAnswer ) );
    for ( auto c = oc.getQuadraturePointMicroInertias( )->begin( ); c != oc.getQuadraturePointMicroInertias( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, microInertiasAnswer, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in micro inertia qpt\n";
            return 1;

        }

    }
    testNum++;

    floatVector bodyCouplesAnswer( 9 * 8, 0 );
    for ( auto c = oc.getQuadraturePointBodyCouples( )->begin( ); c != oc.getQuadraturePointBodyCouples( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, bodyCouplesAnswer, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in body couples qpt\n";
            return 1;

        }

    }
    testNum++;

    floatVector microSpinInertiasAnswer( 9 * 8, 0 );
    for ( auto c = oc.getQuadraturePointMicroSpinInertias( )->begin( ); c != oc.getQuadraturePointMicroSpinInertias( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, microSpinInertiasAnswer, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in micro spin inertias qpt\n";
            return 1;

        }

    }
    testNum++;

    floatVector symmetricMicroStressQptAnswer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    floatVector symmetricMicroStressAnswer = vectorTools::appendVectors( floatMatrix( 8, symmetricMicroStressQptAnswer ) );
    for ( auto c = oc.getQuadraturePointSymmetricMicroStress( )->begin( ); c != oc.getQuadraturePointSymmetricMicroStress( )->end( ); c++ ){

        if ( !vectorTools::fuzzyEquals( c->second, symmetricMicroStressAnswer, 1e-6, 1e-3 ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in symmetric micro stress qpt\n";
            return 1;

        }

    }

    //check the outputted DOF
    shared_ptr< XdmfReader > reader = XdmfReader::New( );
    shared_ptr< XdmfDomain > _readDomain = shared_dynamic_cast< XdmfDomain >( reader->read( "macroscale_dof.xdmf" ) );
    shared_ptr< XdmfUnstructuredGrid > _readGrid = _readDomain->getGridCollection( 0 )->getUnstructuredGrid( 0 );
    shared_ptr< XdmfAttribute > nodeIDs = _readGrid->getAttribute( "node_ids" );
    shared_ptr< XdmfAttribute > updatedDOF = _readGrid->getAttribute( "updated_DOF" );

    nodeIDs->read( );
    updatedDOF->read( );

    uIntVector macroNodeIdsAnswer
       = { 5,  9,  8, 11,  3,  1,  6, 15, 12,  2, 13, 14 };

    if ( nodeIDs->getSize( ) != macroNodeIdsAnswer.size( ) ){

        results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
        std::cerr << "  Failure in macro node ids ( length )\n";
        return 1;

    }

    for ( uIntType i = 0; i < macroNodeIdsAnswer.size( ); i++ ){

        if ( !vectorTools::fuzzyEquals( macroNodeIdsAnswer[ i ], nodeIDs->getValue< uIntType >( i ) ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in macro node ids ( values )\n";
            return 1;

        }

    }
    testNum++;

    if ( updatedDOF->getSize( ) != oc.getUpdatedMacroDisplacementDOF( ).size( ) ){

        results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
        std::cerr << "  Failure in updated macro displacement DOF ( length )\n";
        return 1;

    }
    testNum++;

    for ( uIntType i = 0; i < oc.getUpdatedMacroDisplacementDOF( ).size( ); i++ ){

        if ( !vectorTools::fuzzyEquals( oc.getUpdatedMacroDisplacementDOF( )[ i ], updatedDOF->getValue< floatType >( i ) ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in updated macro displacement DOF ( value )\n";
            return 1;

        }

    }
    testNum++;
    
    reader = XdmfReader::New( );
    _readDomain = shared_dynamic_cast< XdmfDomain >( reader->read( "microscale_dof.xdmf" ) );
    _readGrid = _readDomain->getGridCollection( 0 )->getUnstructuredGrid( 0 );
    nodeIDs = _readGrid->getAttribute( "node_ids" );
    updatedDOF = _readGrid->getAttribute( "updated_DOF" );

    nodeIDs->read( );
    updatedDOF->read( );

    uIntVector microNodeIdsAnswer
       = { 15, 31, 13, 26, 53, 21, 37, 48,  5, 10,  3,  4, 32, 33, 34, 28, 25,
           50, 43, 27,  1,  7, 30, 16, 22,  2, 46, 24, 39, 40, 57, 44, 58, 29,
           59, 11,  0, 20, 60, 47, 49, 17, 38, 14, 55 };

    if ( nodeIDs->getSize( ) != microNodeIdsAnswer.size( ) ){

        results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
        std::cerr << "  Failure in updated micro displacement ids ( length )\n";
        return 1;

    }

    for ( uIntType i = 0; i < microNodeIdsAnswer.size( ); i++ ){

        if ( !vectorTools::fuzzyEquals( microNodeIdsAnswer[ i ], nodeIDs->getValue< uIntType >( i ) ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in updated micro displacement ids ( value )\n";
            return 1;

        }

    }
    testNum++;

    if ( updatedDOF->getSize( ) != oc.getUpdatedMicroDisplacementDOF( ).size( ) ){

        results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
        std::cerr << "  Failure in updated micro displacement DOF ( length )\n";
        return 1;

    }
    testNum++;

    for ( uIntType i = 0; i < oc.getUpdatedMicroDisplacementDOF( ).size( ); i++ ){

        if ( !vectorTools::fuzzyEquals( oc.getUpdatedMicroDisplacementDOF( )[ i ], updatedDOF->getValue< floatType >( i ) ) ){

            results << "test_overlapCoupling_processIncrement_Arlequin (test " + std::to_string( testNum ) + ") & False\n";
            std::cerr << "  Failure in updated micro displacement DOF ( value )\n";
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

    results << "test_overlapCoupling_processIncrement_Arlequin & True\n";
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

int test_overlapCoupling_runOverlapCoupling_Arlequin_realistic( std::ofstream &results ){
    /*!
     * Test the function that runs the overlap coupling for a more realistic macro and
     * micro-scale using the Arlequin method
     *
     * :param std::ofstream &results: The output file
     */

    remove( "macroscale_dof.xdmf" );
    remove( "macroscale_dof.h5" );

    remove( "microscale_dof.xdmf" );
    remove( "microscale_dof.h5" );

    std::string filename = "testConfig_Arlequin_realistic.yaml";

    std::unordered_map< uIntType, uIntType > microGlobalLocalNodeMap;
    std::unordered_map< uIntType, uIntType > macroGlobalLocalNodeMap;

    floatVector updatedMicroDisplacementDOFResult;
    floatVector updatedMacroDisplacementDOFResult;

    floatVector microAugmentedLagrangianForceResult;
    floatVector macroAugmentedLagrangianForceResult;

    errorOut error = overlapCoupling::runOverlapCoupling( filename,
                                                          microGlobalLocalNodeMap, updatedMicroDisplacementDOFResult,
                                                          microAugmentedLagrangianForceResult,
                                                          macroGlobalLocalNodeMap, updatedMacroDisplacementDOFResult,
                                                          macroAugmentedLagrangianForceResult
                                                        );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_runOverlapCoupling_Arlequin_realistic & False\n";
        return 1;
    }

    std::cout << "updated micro displacement DOF Result:\n";
    vectorTools::print(updatedMicroDisplacementDOFResult);
    std::cout << "\n";
    std::cout << "updated macro displacement DOF Result:\n";
    vectorTools::print(updatedMacroDisplacementDOFResult);

    remove( "macroscale_dof.xdmf" );
    remove( "macroscale_dof.h5" );

    remove( "microscale_dof.xdmf" );
    remove( "microscale_dof.h5" );

    results << "test_overlapCoupling_runOverlapCoupling_Arlequin_realistic & True\n";
    return 0;
}

int test_overlapCoupling_initializeCoupling_Arlequin_reconstructed( std::ofstream &results ){
    /*!
     * Test the initialization of the coupling for the Arlequin method when the reconstructed
     * mass centers are used.
     *
     * :param std::ofstream &results: The output file.
     */

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    std::string filename = "testConfig_Arlequin_reconstructed.yaml";
    overlapCoupling::overlapCoupling oc( filename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_initializeCoupling_Arlequin_reconstructed & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_initializeCoupling_Arlequin_reconstructed & False\n";
        return 1;
    }

    if ( !std::ifstream( "reference_information.xdmf" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_Arlequin_reconstructed (test 1) & False\n";
        return 1;

    }

    if ( !std::ifstream( "reference_information.h5" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_Arlequin_reconstructed (test 2) & False\n";
        return 1;

    }

    if ( !std::ifstream( "homogenized_response.xdmf" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_Arlequin_reconstructed (test 3) & False\n";
        return 1;

    }

    if ( !std::ifstream( "homogenized_response.h5" ).good( ) ){

        results << "test_overlapCoupling_initializeCoupling_Arlequin_reconstructed (test 4) & False\n";
        return 1;

    }

    std::string testName = "overlapCoupling_initializeCoupling_Arlequin_reconstructed";
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
                    { "free_nodeset_volume_5" , 0.25 },
                    { "free_nodeset_volume_6" , 0.25 },
                    { "free_nodeset_volume_7" , 0.25 },
                    { "free_nodeset_volume_8" , 0.25 }
                }
            }
        };

    if ( _compare_cellDomainMaps( results, domainMassAnswer, oc._test_domainMass,
                                  testName, testNum, 1e-3, 1e-3 ) ){
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
                    { "free_nodeset_volume_5", { 0.250000, 0.250000, 1.750000 } },
                    { "free_nodeset_volume_6", { 0.750000, 0.250000, 1.750000 } },
                    { "free_nodeset_volume_7", { 0.750000, 0.750000, 1.750000 } },
                    { "free_nodeset_volume_8", { 0.250000, 0.750000, 1.750000 } }
                }
            }
        };

    if ( _compare_cellDomainMaps( results, domainCOMAnswer, oc._test_domainCOM,
                                  testName, testNum, 1e-3, 1e-3 ) ){
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
                    { "free_nodeset_volume_5", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "free_nodeset_volume_6", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "free_nodeset_volume_7", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } },
                    { "free_nodeset_volume_8", { +0.062500, +0.000000, +0.000000, +0.000000, +0.062500, +0.000000, +0.000000, +0.000000, +0.062500 } }
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
                            { 31, { -0.250000, -0.250000, -0.250000,  } },
                            { 43, { -0.250000, -0.250000, +0.250000,  } },
                            { 26, { +0.250000, -0.250000, -0.250000,  } },
                            { 27, { +0.250000, -0.250000, +0.250000,  } },
                            { 21, { -0.250000, +0.250000, -0.250000,  } },
                            {  1, { -0.250000, +0.250000, +0.250000,  } },
                            { 48, { +0.250000, +0.250000, -0.250000,  } },
                            {  7, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_6",
                        {
                            { 26, { -0.250000, -0.250000, -0.250000,  } },
                            { 27, { -0.250000, -0.250000, +0.250000,  } },
                            { 10, { +0.250000, -0.250000, -0.250000,  } },
                            { 30, { +0.250000, -0.250000, +0.250000,  } },
                            { 48, { -0.250000, +0.250000, -0.250000,  } },
                            {  7, { -0.250000, +0.250000, +0.250000,  } },
                            {  4, { +0.250000, +0.250000, -0.250000,  } },
                            { 16, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_7",
                        {
                            { 48, { -0.250000, -0.250000, -0.250000,  } },
                            {  7, { -0.250000, -0.250000, +0.250000,  } },
                            {  4, { +0.250000, -0.250000, -0.250000,  } },
                            { 16, { +0.250000, -0.250000, +0.250000,  } },
                            { 33, { -0.250000, +0.250000, -0.250000,  } },
                            { 22, { -0.250000, +0.250000, +0.250000,  } },
                            { 28, { +0.250000, +0.250000, -0.250000,  } },
                            {  2, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                    {
                        "free_nodeset_volume_8",
                        {
                            { 21, { -0.250000, -0.250000, -0.250000,  } },
                            {  1, { -0.250000, -0.250000, +0.250000,  } },
                            { 48, { +0.250000, -0.250000, -0.250000,  } },
                            {  7, { +0.250000, -0.250000, +0.250000,  } },
                            { 50, { -0.250000, +0.250000, -0.250000,  } },
                            { 46, { -0.250000, +0.250000, +0.250000,  } },
                            { 33, { +0.250000, +0.250000, -0.250000,  } },
                            { 22, { +0.250000, +0.250000, +0.250000,  } },
                        }
                    },
                }
            }
        };

    if ( _compare_cellDomainPointMaps( results, domainXiAnswer, oc._test_domainXi,
                                       testName, testNum ) ){
        return 1;
    }

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
                        { +0.140625, +0.046875, +0.015625, +0.046875, +0.421875, +0.140625, +0.046875, +0.140625 }
                    },
                    {
                        "free_nodeset_volume_6",
                        { +0.046875, +0.140625, +0.046875, +0.015625, +0.140625, +0.421875, +0.140625, +0.046875 }
                    },
                    {
                        "free_nodeset_volume_7",
                        { +0.015625, +0.046875, +0.140625, +0.046875, +0.046875, +0.140625, +0.421875, +0.140625 }
                    },
                    {
                        "free_nodeset_volume_8",
                        { +0.046875, +0.015625, +0.046875, +0.140625, +0.140625, +0.046875, +0.140625, +0.421875 }
                    }
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
    A << 3.01086178e-01,  3.73083127e-01, -2.38128356e-02,  2.59667767e-01,
         3.81825184e-01,  2.47669483e-03,  4.57671591e-01,  4.62375983e-01,
         2.14785351e-01,  4.38818433e-01,  4.17093025e-01,  2.23207977e-01,
         5.51383155e-01,  5.19387478e-01,  1.77949499e-01,  5.03685954e-01,
         5.01863971e-01,  2.48959916e-01,  6.53689216e-01,  6.08389210e-01,
         3.25284090e-01,  6.18922811e-01,  5.46460887e-01,  3.60351963e-01,
         6.24685518e-01,  5.76649979e-01,  4.44867399e-01,  6.35748959e-01,
         4.90362784e-01,  4.21448837e-01,  7.99115468e-01,  6.75194985e-01,
         4.66079526e-01,  7.87324774e-01,  5.85351468e-01,  4.57872552e-01,
         8.05310891e-01,  8.62517011e-01,  4.62486565e-01,  8.02914650e-01,
         7.74076152e-01,  5.59111811e-01,  9.37733511e-01,  8.68813552e-01,
         5.50358901e-01,  9.68673837e-01,  7.58454610e-01,  5.99295180e-01,
         7.48700137e-01,  7.86847416e-01,  3.70052056e-01,  7.25705814e-01,
         7.40283104e-01,  5.13675948e-01,  3.92083649e-01,  4.48868899e-01,
         1.97336971e-01,  6.10254590e-01,  4.77522623e-01,  3.99710841e-01,
         5.98624443e-01,  5.62847803e-01,  4.67857453e-01,  7.62302834e-01,
         6.06590588e-01,  5.40978522e-01,  8.52098910e-01,  5.40964575e-01,
         5.80758732e-01,  9.80779591e-01,  6.49959451e-01,  6.04506237e-01,
         9.56672848e-01,  8.11585854e-01,  7.63645571e-01,  1.19021276e+00,
         8.12162668e-01,  7.51572977e-01,  8.09056288e-01,  7.75472967e-01,
         7.77857436e-01,  2.30540336e-01,  6.64703123e-04, -6.27165920e-02,
         3.15448579e-01,  2.43688429e-01,  5.54817991e-02,  3.78963633e-01,
         1.28344015e-01,  6.84807385e-02,  4.64507330e-01,  3.64486259e-01,
         2.25097533e-01,  3.99475103e-01,  2.66251345e-01,  1.64820730e-01,
         5.15440047e-01,  4.70753468e-01,  2.69561731e-01,  5.21078101e-01,
         3.48998935e-01,  2.56955259e-01,  6.34620679e-01,  5.58838000e-01,
         3.75343353e-01,  5.42313346e-01,  1.93404759e-01,  2.69184455e-01,
         6.25514633e-01,  4.58348143e-01,  4.32777754e-01,  6.47009543e-01,
         3.64429355e-01,  4.09158872e-01,  7.73319801e-01,  5.96587366e-01,
         5.13695143e-01,  6.46419743e-01,  5.70863574e-01,  4.36057418e-01,
         7.69225089e-01,  8.01625423e-01,  5.26690549e-01,  7.44822756e-01,
         5.40980183e-01,  5.30784472e-01,  8.92474568e-01,  7.81152975e-01,
         6.05166075e-01,  5.41747205e-01,  5.28731195e-01,  3.91962147e-01,
         6.73064042e-01,  7.48363484e-01,  4.75290871e-01;

    Eigen::MatrixXd macroD( N.cols( ), 1 );

    macroD << 0.69646919, 0.28613933, 0.22685145, 0.55131477, 0.71946897,
              0.42310646, 0.9807642 , 0.68482974, 0.4809319 , 0.39211752,
              0.34317802, 0.72904971, 0.43857224, 0.0596779 , 0.39804426,
              0.73799541, 0.18249173, 0.17545176, 0.53155137, 0.53182759,
              0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338,
              0.32295891, 0.36178866, 0.22826323, 0.29371405, 0.63097612,
              0.09210494, 0.43370117, 0.43086276, 0.4936851 , 0.42583029,
              0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668,
              0.62395295, 0.1156184 , 0.31728548, 0.41482621, 0.86630916,
              0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
              0.12062867, 0.8263408 , 0.60306013, 0.54506801, 0.34276383,
              0.30412079, 0.41702221, 0.68130077, 0.87545684, 0.51042234,
              0.66931378, 0.58593655, 0.6249035 , 0.67468905, 0.84234244,
              0.08319499, 0.76368284, 0.24366637, 0.19422296, 0.57245696,
              0.09571252, 0.88532683, 0.62724897, 0.72341636, 0.01612921,
              0.59443188, 0.55678519, 0.15895964, 0.15307052, 0.69552953,
              0.31876643, 0.6919703 , 0.55438325, 0.38895057, 0.92513249,
              0.84167   , 0.35739757, 0.04359146, 0.30476807, 0.39818568,
              0.70495883, 0.99535848, 0.35591487, 0.76254781, 0.59317692,
              0.6917018 , 0.15112745, 0.39887629, 0.2408559 , 0.34345601,
              0.51312815, 0.66662455, 0.10590849, 0.13089495, 0.32198061,
              0.66156434, 0.84650623, 0.55325734, 0.85445249, 0.38483781,
              0.3167879 , 0.35426468, 0.17108183, 0.82911263, 0.33867085,
              0.55237008, 0.57855147, 0.52153306, 0.00268806, 0.98834542,
              0.90534158, 0.20763586, 0.29248941, 0.52001015, 0.90191137,
              0.98363088, 0.25754206, 0.56435904, 0.80696868, 0.39437005,
              0.73107304, 0.16106901, 0.60069857, 0.86586446, 0.98352161,
              0.07936579, 0.42834727, 0.20454286, 0.45063649, 0.54776357,
              0.09332671, 0.29686078, 0.92758424, 0.56900373;

    Eigen::MatrixXd R = N * macroD;

    if ( ( A - R ).norm( ) > ( 1e-6 * A.norm( ) ) + 1e-6 ){

        results << "test_" + testName +" (" + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    //Check the center of mass interpolation matrix and the center of mass projector
    SparseMatrix centerOfMassInterpolator;
    error = overlapCoupling::readSparseMatrixFromXDMF( _readGrid, "centerOfMassInterpolator", centerOfMassInterpolator );

    if ( error ){
        error->print( );
        results << "test_" + testName + " & False\n";
        return 1;
    }

    Eigen::MatrixXd DX( 12, 1 );
    DX << -1.00911786,  1.51428288,  1.75159184, -0.77596151, -0.13860077,
          -1.30538174, -1.11042458, -0.86808735,  0.47158175, -1.21084958,
           1.4369616 , -0.41944997;

    Eigen::MatrixXd PA( 16, 1 );
    PA << -0.40381903, -0.85441036, -0.64758047, -0.59150002, -0.10167723,
          -0.51778121,  0.13886224, -0.16634456, -0.37851177,  0.45119406,
           0.53527639, -0.35298579, -0.49609721, -0.53141861, -0.51544242,
          -0.65371376;
    
    R = centerOfMassInterpolator * DX;

    if ( ( R - PA ).norm( ) > 1e-6 * ( PA.norm( ) + 1 ) ){

        results << "test_" + testName + " (" + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    Eigen::MatrixXd centerOfMassProjector;
    error = overlapCoupling::readDenseMatrixFromXDMF( _readGrid, "centerOfMassProjector", centerOfMassProjector );

    if ( error ){
        error->print( );
        results << "test_" + testName + " & False\n";
        return 1;
    }

    R = centerOfMassProjector * PA;

    if ( ( R - DX ).norm( ) > 1e-6 * ( DX.norm( ) + 1 ) ){

        results << "test_" + testName + " (" + std::to_string( testNum + 1 ) + ") & False\n";
        return 1;

    }
    testNum += 1;

    remove( "reference_information.xdmf" );
    remove( "reference_information.h5" );

    remove( "homogenized_response.xdmf" );
    remove( "homogenized_response.h5" );

    results << "test_" + testName + " & True\n";
    return 0;
}

int temp_processBigFile( ){


    std::string filename = "testConfig.yaml";
    overlapCoupling::DOFMap d1, d2;
    floatVector f1, f2, f3, f4;
    errorOut error = overlapCoupling::runOverlapCoupling( filename, d1, f1, f2, d2, f3, f4 );
    if ( error ){
        error->print( );
    }

}

int test_overlapCoupling_filterMode( std::ofstream &results ){
    /*!
     * Test the overlap coupling used in filter mode
     * 
     * :param std::ofstream &results: The output file
     */

    std::string inputFilename = "filterConfig.yaml";

    overlapCoupling::overlapCoupling oc( inputFilename );

    if ( oc.getConstructorError( ) ){
        oc.getConstructorError( )->print( );
        results << "test_overlapCoupling_filterMode & False\n";
        return 1;
    }

    errorOut error = oc.initializeCoupling( );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_filterMode & False\n";
        return 1;
    }

    error = oc.processIncrement( 10, 0 );

    if ( error ){
        error->print( );
        results << "test_overlapCoupling_filterMode & False\n";
        return 1;
    }

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

//    test_overlapCoupling_constructor( results );
//    test_overlapCoupling_initializeCoupling_l2_projection( results );
//    test_overlapCoupling_initializeCoupling_averaged_l2_projection( results );
//    test_overlapCoupling_initializeCoupling_Arlequin( results );
//    test_overlapCoupling_initializeCoupling_Arlequin_reconstructed( results );
//    test_overlapCoupling_processIncrement( results );
//    test_overlapCoupling_processIncrement_Arlequin( results );
//    test_overlapCoupling_processLastIncrements( results );
////    test_overlapCoupling_getReferenceFreeMicroDomainMasses( results );
////    test_overlapCoupling_getReferenceGhostMicroDomainMasses( results );
////    test_overlapCoupling_getReferenceFreeMicroDomainCentersOfMass( results );
////    test_overlapCoupling_getReferenceGhostMicroDomainCentersOfMass( results );
//////    test_overlapCoupling_getReferenceFreeMicroDomainCenterOfMassShapeFunctions( results );
//////    test_overlapCoupling_getReferenceGhostMicroDomainCenterOfMassShapeFunctions( results );
//    test_overlapCoupling_runOverlapCoupling_Arlequin_realistic( results );
    test_overlapCoupling_filterMode( results );
//    test_MADOutlierDetection( results );
//    test_formMicromorphicElementMassMatrix( results );
//    test_computeMicromorphicElementRequiredValues( results );
//    test_computeMicromorphicElementInternalForceVector( results );
//    test_readWriteSparseMatrixToXDMF( results );
//    test_readWriteDenseMatrixToXDMF( results );
//
//    temp_processBigFile( );

    //Close the results file
    results.close();
}
