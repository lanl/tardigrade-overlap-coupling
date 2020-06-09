/*!
===============================================================================
|                               assembly.cpp                                  |
===============================================================================
| A collection of utilities to assist in assembling a finite element problem. |
===============================================================================
*/

#include "assembly.h"

namespace assembly{

        int split_string(const std::string &line, const std::string &delimiter, std::vector< std::string > &parsed_line){
            /*!
            Split the provided string at the given delimiter
            */

            size_t last, next;
            last = next = 0;
            parsed_line.resize(0);

            while ((next = line.find(delimiter, last)) != std::string::npos){
                parsed_line.push_back(line.substr(last, next - last));
                last = next + 1;
            }
	    if(last<line.size()){
                parsed_line.push_back(line.substr(last));
	    }
            return 0;
        }

        int read_past_header(std::ifstream &file){
            /*!
            Read past the header information

            :param std::ifstream file: The input file
            */

            std::string line;

            while (std::getline(file, line)){
                if (std::strcmp(line.c_str(), "BEGIN DATA") == 0){
                    return 0;
                }
            }
            return 1;
        }

        int parsed_line_to_data(const std::vector< std::string > &parsed_line, node_map &nodes, element_map &elements, qrule_map &qrules){
            /*!
             * Convert a parsed line into data
             *
             * :param const std::vector< std::string > &parsed_line: A parsed line from the input file
	     * :param node_map &nodes: The coordinates of the nodes which make up the elements.
	     * :param element_map &elements: The connectivity of each of the nodes to form elements.
	     * :param qrule_map &qrules: The quadrature rule for the elements.
            */

            //Add a node to the nodes map
            if (std::strcmp(parsed_line[0].c_str(), "N")==0){
                    //The line is a node definition
		    unsigned int node_id = std::stoul(parsed_line[1]);
		    std::vector< double > line_data(parsed_line.size()-2);
		    for (unsigned int i=2; i<parsed_line.size(); i++){
                            line_data[i-2] = std::stod(parsed_line[i]);
		    }

		    //Check if the node was already defined
		    node_map::iterator it = nodes.find(node_id);
		    if (it != nodes.end()){
                            std::cout << "Error: node " << node_id << " is defined twice.\n";
		            return 1;
		    }
		    else{
			    nodes.emplace(node_id, line_data);
		    }

		    return 0;
            }
            else if (std::strcmp(parsed_line[0].c_str(), "E")==0){
		    //The line is an element definition

		    unsigned int elid = std::stoul(parsed_line[2]);
                    std::vector< unsigned int > line_data(parsed_line.size()-3);
		    for (unsigned int i=3; i<parsed_line.size(); i++){
                        line_data[i-3] = std::stoul(parsed_line[i]);
		    }

                    std::map<unsigned int, std::vector< unsigned int > > emap;
		    emap.emplace(elid, line_data);

		    element_map::iterator it = elements.find(parsed_line[1]);

		    if (it == elements.end()){
			    std::cout << "Adding element type: " << parsed_line[1] << "\n";
			    
			    elements.emplace(parsed_line[1], emap);
		    }
		    else{
			    std::map< unsigned int, std::vector< unsigned int > >::iterator itid = it->second.find(elid);

			    if (itid != it->second.end()){
                                std::cout << "Error: element " << elid << " already defined.\n";
				return 1;
			    }
			    it->second.emplace( elid, line_data);
		    }

		    return 0;

	    }
	    else if (std::strcmp(parsed_line[0].c_str(), "Q")==0){
                    //The line is an element type quadrature rule

                    std::vector< std::pair< std::vector< double >, double > > el_qrule;
		    unsigned int dim = std::stoul(parsed_line[2]);
		    std::vector< double > pt(dim, 0);
		    double w;
		    for (unsigned int i=3; i<parsed_line.size(); i+=(dim+1)){
                        for (unsigned int j=0; j<dim; j++){
                            pt[j] = std::stod(parsed_line[i+j]);
			}
			w = std::stod(parsed_line[i+dim]);
			el_qrule.push_back(std::pair< std::vector< double > , double >(pt, w));
		    }

		    qrule_map::iterator it = qrules.find(parsed_line[0]);

		    if (it != qrules.end()){
                            std::cout << "Error: quadrature rule for " << parsed_line[0] << " already defined.\n";
			    return 1;
		    }
		    else{
                        qrules.emplace(parsed_line[1], el_qrule);
		    }
	    }

            return 1;
        }

	int read_connectivity_data(const std::string &input_filename, node_map &nodes, element_map &elements, qrule_map &qrules){
            /*!
	     * Read in the connectivity data from a text file. This will be the location of the nodes and the connectivity of the 
	     * nodes to form elements.
	     *
	     * :param const std::string &input_filename: The name of the input file
	     * :param node_map &nodes: The coordinates of the nodes which make up the elements.
	     * :param element_map &elements: The connectivity of each of the nodes to form elements.
	     * :param qrule_map &qrules: The quadrature rules associated with the elements.
	    */

	//Open the file
        std::ifstream file;
	file.open(input_filename);

        if (!file.is_open()){

		std::cout << "Error: file failed to open.\n";
		return 1;
	}

	//Skip the header
	int rph_retvalue = read_past_header(file);
	if (rph_retvalue > 0){
                std::cout << "Error: problem skipping header\n";
		return 1;
	}

        //Read in the lines
	std::string line;
	std::vector< std::string > split_line;
	while(std::getline(file, line)){
                remove_blank_spaces(line);
		if (line.size() > 0){

			int ss_retvalue = split_string(line, ",", split_line);
			if (ss_retvalue>0){
                            return 1;
			}
			int pld_retvalue = parsed_line_to_data(split_line, nodes, elements, qrules);
			if (pld_retvalue>0){
                            return 1;
                        }
		}
            
	}

        return 0;
}
    int print_node_map(const node_map &nodes){
        /*!
	 * Print the node map to the terminal
	 *
	 * :param const node_map &nodes: The node_map to be printed
	*/

	for (auto it=nodes.begin(); it!=nodes.end(); it++){
            std::cout << it->first << ": ";
	    for (unsigned int i=0; i<it->second.size(); i++){
                std::cout << it->second[i] << " ";
	    }
	    std::cout << "\n";
	}
	return 0;
    }

    int print_element_map(const element_map &elements){
        /*!
	 * Print the element map to the terminal
	 *
	 * :param const element_map &elements: The element_map to be printed
	*/

	for (auto it=elements.begin(); it!=elements.end(); it++){
            std::cout << it->first << "\n";
	    for (auto elit=it->second.begin(); elit!=it->second.end(); elit++){
		    std::cout << "   " << elit->first << ": ";
		    for (unsigned int i=0; i<elit->second.size(); i++){
                        std::cout << elit->second[i] << " ";
		    }
		    std::cout << "\n";

	    }
	}
	return 0;
    }

    int print_qrule_map(const qrule_map &qrules){
        /*!
	 * Print the quadrature rules to the terminal
	 *
	 * :param const qrule_map &qrules: The qrule_map to the printed
	*/

	for (auto it=qrules.begin(); it!=qrules.end(); it++){
            std::cout << it->first << "\n";
	    for (unsigned int i=0; i<it->second.size(); i++){
                std::cout << "   ";
                for (unsigned int j=0; j<it->second[i].first.size(); j++){
                    std::cout << it->second[i].first[j] << " ";
		}
		std::cout << "(" << it->second[i].second << ")\n";
	    }
	}
	return 0;
    }

    int remove_blank_spaces(std::string &str){
        /*!
	 * Remove blank spaces from a string
	 *
	 * :param std::string &str: The string to modify
	*/

	str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
	return 0;
    }
}
