/*
A test file for the assembly library
*/

#include<iostream>
#include<fstream>
#include<assembly.h>

int test_read_connectivity_data(std::ofstream &results){
    /*
     * Test reading in connectivity data from an example file.
    */

    std::string input_filename = "connectivity.txt";
    assembly::node_map nodes;
    assembly::element_map elements;
    assembly::qrule_map qrules;

    assembly::read_connectivity_data(input_filename, nodes, elements, qrules);

    std::cout << "Nodes\n";
    assembly::print_node_map(nodes);
    std::cout << "Elements\n";
    assembly::print_element_map(elements);
    std::cout << "Qrules\n";
    assembly::print_qrule_map(qrules);

    if (nodes.size() != 8){
        results << "test_read_connectivity_data (test 1) & False\n";
    }

    if (elements.size() != 1){
        results << "test_read_connectivity_data (test 2) & False\n";
    }

    results << "test_read_connectivity_data & True\n";
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

    //Test reading connectivity data from a text file
    test_read_connectivity_data(results);

    //Close the results file
    results.close();

    return 0;
}
