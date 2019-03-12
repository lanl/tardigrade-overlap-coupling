//!The test file for overlap_coupling.cpp

#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>

#include "overlap_coupling.h"

template<typename T>
void print_vector(std::vector< T > vector){
    /*!
    Print the vector to the terminal
    */
    for (unsigned int i=0; i<vector.size(); i++){
        std::cout << vector[i] << " ";
    }
    std::cout << "\n";
    return;
}

template<typename T>
void print_matrix(std::vector< std::vector< T > > matrix){
    /*!
    Print the matrix to the terminal.
    */
    for (unsigned int i=0; i<matrix.size(); i++){
        print_vector(matrix[i]);
    }
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

    overlap::ParsedData data = overlap::read_data_from_file("overlap.txt");
    print_matrix(data.global_nodes);
    print_matrix(data.local_nodes);
    print_vector(data.node_numbers);
    print_vector(data.volumes);
    print_vector(data.densities);
    print_matrix(data.coordinates);

    //Close the results file
    results.close();
}
