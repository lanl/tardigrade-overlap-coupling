/*
A test file for the assembly library
*/

#include<iostream>
#include<fstream>
#include<assembly.h>

#define BOOST_TEST_MODULE test_assembly
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE( testRead_connectivity_data ){
    /*
     * Test reading in connectivity data from an example file.
    */

    std::string input_filename = "assembly_connectivity.txt";
    assembly::node_map nodes;
    assembly::element_map elements;
    assembly::qrule_map qrules;

    assembly::read_connectivity_data(input_filename, nodes, elements, qrules);

    BOOST_CHECK( nodes.size() == 8);

    BOOST_CHECK( elements.size() == 1 );

    BOOST_CHECK( qrules.size() == 1 );

    for (auto it = qrules.begin(); it != qrules.end(); it++){
        BOOST_CHECK( std::strcmp( it->first.c_str( ), "Hex8" ) == 0 );
    }

}
