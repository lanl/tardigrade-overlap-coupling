/*!
===============================================================================
|                                 assembly.h                                  |
===============================================================================
| A collection of utilities to assist in assembling a finite element problem. |
===============================================================================
*/

#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include<string>
#include<cstring>
#include<vector>
#include<map>
#include<iostream>
#include<fstream>
#include<algorithm>

namespace assembly{

        typedef std::map< unsigned int, std::vector< double > > node_map;
	typedef std::map< std::string, std::map< unsigned int, std::vector< unsigned int > > > element_map;
	typedef std::map< std::string, std::vector< std::pair< std::vector< double >, double > > > qrule_map;

	int read_connectivity_data(const std::string &input_filename, node_map &nodes, element_map &connectivity, qrule_map &qrules);
	int read_past_header(std::ifstream &file);
	int split_string(const std::string &line, const std::string &delimiter, std::vector< std::string > &parsed_line);
        int parsed_line_to_data(const std::vector< std::string > &parsed_line, node_map &nodes, element_map &elements, qrule_map &qrules);
        int remove_blank_spaces(std::string &str);

	int print_node_map(const node_map &nodes);
	int print_element_map(const element_map &elements);
	int print_qrule_map(const qrule_map &qrules);
}

#endif
