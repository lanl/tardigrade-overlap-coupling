/*!
===============================================================================
|                                 filter.cpp                                  |
===============================================================================
| An implementation of the micromorphic filter using the overlap coupling     |
| library.                                                                    |
===============================================================================
*/

#include "filter.h"

namespace filter{

    int open_format1_file(const std::string &fn, std::ifstream &file){
        /*!
        Open a file in format 1
    
        :param std::string fn: The filename
        */
    
        file = std::ifstream(fn);
        if(file.is_open()){
            return 0;
        }
        else{
            std::cout << "Error: cannot open file: " << fn << "\n";
            return 1;
        }
    }
    
    int open_input_file(const std::string &fn, const int format, std::ifstream &file){
        /*!
        Open an input file of the provided format.
    
        :param const std::string &fn: The input filename
        :param const int format: The file format
        :param std::ifstream &file: The opened file.
        */
    
        if (format == 1){
            return open_format1_file(fn, file);
        }
    
        return 1;
    }
    
    int read_past_header(std::ifstream &file, input_format &mp_format, input_format &dof_format, const int format){
        /*!
         * Read past the file header for the given format
         * 
         * :param const std::ifstream &file: The input file
         * :param input_format &mp_format: The format of a material point line
         * :param input_format &dof_format: The format of a dof point line
         * :param int format: The format of the file
         */
    
        if (format==1){
            return read_past_header_format1(file, mp_format, dof_format);
        }
    
        return 1;
    }

// FOLLOWING LINES FROM https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring

    // trim from start (in place)
    static inline void ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
            return !std::isspace(ch);
        }));
    }

    // trim from end (in place)
    static inline void rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
            return !std::isspace(ch);
        }).base(), s.end());
    }

    // trim from both ends (in place)
    static inline void trim(std::string &s) {
        ltrim(s);
        rtrim(s);
    }   
// END LINES FROM https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring    

    int set_format(std::string &line, input_format &format){
        /*!
         * Set the format from the input line
         * 
         * :param std::string &line: A line from the input file
         * :param input_format &format: The format to be specified
         */
         
         std::vector< std::string > sline;
         split_string(line, ",", sline);
         std::vector< unsigned int > vals(2, 0);

         for (unsigned int i=1; i<sline.size(); i+=3){
             vals[0] = std::stol(sline[i+1]);
             vals[1] = std::stol(sline[i+2]);
             format.emplace(sline[i], vals);
         }
         return 0;
    }

    int read_past_header_format1(std::ifstream &file, input_format &mp_format, input_format &dof_format){
        /*!
         * Read past the header information for format 1
         * 
         * :param std::ifstream file: The input file
         * :param input_format &mp_format: The format of a material point line
         * :param input_format &dof_format: The format of a dof point line
         */
    
        std::string line;
   
        std::cout << "File header:\n\n";
 
        while (std::getline(file, line)){
            if (line.find("*MPFORMAT") != std::string::npos){
                set_format(line, mp_format);
            }
            else if (line.find("*DOFFORMAT") != std::string::npos){
                set_format(line, dof_format);
            }

            else if (std::strcmp(line.c_str(), "BEGIN DATA") == 0){
                std::cout << "\n";
                return 0;
            }
            else{
                std::cout << line << "\n";
            }
        }
        return 1;
    }
    
    int read_timestep(std::ifstream &file, const int format, std::ofstream &output_file, elib::vecOfvec &data){
        /*!
        Read in the information for a single timestep in the provided format.
    
        :param std::ifstream &file: The input file to read
        :param const int format: The input file format
        :param std::ofstream &output_file: The output file
        :param elib::vecOfvec &data: The parsed data. Note: All data is converted to double
        */
    
        double time;
        int result;
        if (format == 1){
            result = read_timestep_format1(file, time, data);
        }
        else{
            result = 1;
        }

        output_file << "*TIMESTEP, " << time << "\n";
        return result;
    }
    
    int find_current_time_format1(std::ifstream &file, elib::vecOfvec &data, double &time){
        /*!
        Read in the current time from the input file.
    
        :param std::ifstream &file: The input file to read
        :param elib::vecOfvec &data: The parsed data. Note: All data is converted to double.
        */
    
        std::string line;
        std::size_t found;
        std::string time_indicator = "t = ";
    
        while (std::getline(file, line)){
            //Search for a line that begins with 't = '
            found = line.find(time_indicator);
            if (found!=std::string::npos){
                time = std::stod(line.substr(found+4, line.size()-(found+4)));
                std::cout << "Retrieving data from timestep: " << time << "\n";
                return 0;
            }
        }
        return 1;
    }
    
    int split_string(const std::string &line, const std::string &delimiter, std::vector< std::string > &parsed_line){
        /*!
        Split the provided string at the given delimiter
        */
    
        size_t last, next;
        last = next = 0;
        parsed_line.resize(0);
        std::string substr;
    
        while ((next = line.find(delimiter, last)) != std::string::npos){
            substr = line.substr(last, next - last);
            trim(substr);
            parsed_line.push_back(substr);//line.substr(last, next - last));
            last = next + 1;
        }
        parsed_line.push_back(line.substr(last));
        return 0;
    }
    
    int parsed_line_to_data(const std::vector< std::string > &parsed_line, elib::vec &line_data){
        /*!
        Convert a parsed line into data
        
        :param const std::vector< std::string > &parsed_line: A parsed line from the input file
        :param elib::vec &line_data: The double representation of the data line
        */
    
        //Convert the string to data
        line_data.resize(parsed_line.size());
        if (std::strcmp(parsed_line[0].c_str(), "MP")==0){
            //Assign a material point a value of 1
            line_data[0] = 1;
        }
        else if (std::strcmp(parsed_line[0].c_str(), "DOFP")==0){
            //Assign a data point a value of 2
            line_data[0] = 2;
        }
        for (unsigned int i=1; i<parsed_line.size(); i++){
            line_data[i] = std::stod(parsed_line[i]);
        }
        return 0;
    }
    
    int read_timestep_data_format1(std::ifstream &file, elib::vecOfvec &data){
        /*!
        Read in the data for the current timestep using format 1
    
        :param std::ifstream &file: The input file to read
        :param elib::vecOfvec &data: The parsed data. Note: All data is converted to double.
        */
    
        std::string line;
        std::string time_indicator = "t = ";
        std::string delimiter = ",";
    
        std::streampos oldpos = file.tellg();
        std::vector< std::string > parsed_line;
        size_t found;
    
        elib::vec line_data;
    
    
        while (std::getline(file, line)){
            if (line.size() > 0){
    
                //Check if a new timestep has been found
                found = line.find(time_indicator);
                if (found!=std::string::npos){
                    file.seekg(oldpos);
                    return 0;
                }
    
                //Split the line at the commas
                split_string(line, delimiter, parsed_line);
    
                //Convert the parsed line to a vector
                parsed_line_to_data(parsed_line, line_data);
                data.push_back(line_data);
            }
        }
    
        return 0;    
    }
    
    int read_timestep_format1(std::ifstream &file, double &time, elib::vecOfvec &data){
        /*!
        Read in the information for a single timestep.
    
        :param std::ifstream &file: The input file to read
        :param elib::vecOfvec &data: The parsed data. Note: All data is converted to double.
        */
    
        //Get the current time
        int get_time_result = find_current_time_format1(file, data, time);
        if (get_time_result>0){
            std::cout << "Error reading in new timestep\n";
            return 1;
        }
    
        //Read in the data
        return read_timestep_data_format1(file, data);
    
    }

    int build_filters(const assembly::node_map &nodes, const assembly::element_map elements,
                      const assembly::qrule_map &qrules, const unsigned int num_macro_dof, filter_map &filters){
        /*!
        * Build the filters for processing data.
        *
        * :param const assembly::node_map &nodes: The nodes of the micromorphic filter.
        * :param const assembly::element_map &elements: The elements of the micromorphic filter.
        * :param const assembly::qrule_map &qrules: The quadrature rules for the elements.
        * :param const unsigned int num_macro_dof: The number of macro-scale degrees of freedom.
        * :param filter_map &filters: Return map of the element ids to their corresponding sub-filters.
        */

        //Clear the filters
        filters.clear();

        //Iterate over the element types.
        for (auto _element_types = elements.begin(); _element_types!=elements.end(); _element_types++){

            //Find the quadrature rule
            auto qrule = qrules.find(_element_types->first);
            if (qrule == qrules.end()){
                std::cout << "Error: quadrature rule for " << _element_types->first << " not found.\n";
                return 1;
            }

            //Loop over the elements. This could be parallelized but it probably won't make a difference.
            for (auto _element = _element_types->second.begin(); _element!=_element_types->second.end(); _element++){

                auto elid = filters.find(_element->first);
                if (elid != filters.end()){
                    std::cout << "Error: filters must be made of elements with unique ids.\n"
                              << "       " << _element->first << " is already used.\n";
                    return 1;
                }

                //Create a vector of the filter's nodes
                elib::vecOfvec element_nodes(_element->second.size());
                for (unsigned int i=0; i<element_nodes.size(); i++){
                    auto point = nodes.find(_element->second[i]);
                    if (point == nodes.end()){
                        std::cout << "Error: node " << _element->second[i] << " not found.\n";
                        return 1;
                    }
                    element_nodes[i] = point->second;
                }

                //Initialize the filter
                filters.emplace(_element->first, overlap::MicromorphicFilter(_element->first, _element_types->first,
                                                                             _element->second, element_nodes, qrule->second,
                                                                             num_macro_dof));

            }
        }
        return 0;
    }

    int populate_filters(const elib::vecOfvec &data, const input_format &mp_format, const input_format &dof_format, 
                        const assembly::node_map &nodes,
                        const assembly::element_map &elements, const assembly::qrule_map &qrules,
                        const bool update_shapefunction, const bool shared_dof_material,
                        const unsigned int num_macro_dof,
                        std::map<unsigned int, unsigned int > &micro_node_to_row,
                        std::map< unsigned int, unsigned int > &micro_node_elcount, filter_map &filters){
        /*!
        * Populate the micromorphic sub-filters using the provided data.
        *
        * :param const elib::vecOfvec &data: The DNS data at the current timestep.
        * :param const input_format &mp_format: The format of a material point dataline
        * :param const input_format &dof_format: The format of a degree of freedom dataline
        * :param const assembly::node_map &nodes: The nodes of the micromorphic filter.
        * :param const assembly::element_map &elements: The elements of the micromorphic filter.
        * :param const assembly::qrule_map &qrules: The quadrature rules for the elements.
        * :param const bool update_shapefunction: Flag indicating of the shapefunction matrix is to be updated.
        * :param const bool shared_dof_material: Flag indicating if the degrees of freedom and material quantities are known at 
        *                                       the same points.
        * :param const unsigned int num_macro_dof: The number of degrees of freedom to macro-scale has.
        * :param std::map< unsigned int, unsigned int > &micro_node_to_row: The mapping from the micro node number to the 
        *                                                                   row of the shape function matrix.
        * :param std::map< unsigned int, unsigned int > &micro_node_elcount: The number of filters a node is contained within.
        * :param filter_map &filters: Return map of the element ids to their corresponding sub-filters.
        */

        //Construct the filters if they haven't been formed yet
        if (filters.size() == 0){
            std::cout << " Filter list unpopulated. Initial construction of filters occuring\n";
            int bf_result = build_filters(nodes, elements, qrules, num_macro_dof, filters);
            if (bf_result>0){
                return 1;
            }
        }
        else{
            for (auto filter=filters.begin(); filter!=filters.end(); filter++){
                filter->second.clear_microscale();
            }
        }

        //Erase micro_node_to_row
        if (update_shapefunction){
            micro_node_to_row.clear();
	    micro_node_elcount.clear();
        }

        //This probably could/should be parallelized
        unsigned int index = 0;
        std::vector< double > position;
        unsigned int nodeid;

        for (auto datapoint=data.begin(); datapoint!=data.end(); datapoint++){
            unsigned int containing_filters = 0;

            //Determine whether the point is a material point or dof point
            int pointtype = (int)((*datapoint)[0]+0.5);

	    //Iterate over the macro-scale filters
            for (auto filter = filters.begin(); filter!=filters.end(); filter++){
                bool iscontained = false;
                
                if (pointtype==1){
                    get_position(*datapoint, mp_format, nodeid, position);
                    iscontained = filter->second.add_micro_material_point(nodeid, position);

		}
		else if (pointtype==2){
                    get_position(*datapoint, dof_format, nodeid, position);
                    iscontained = filter->second.add_micro_dof_point(nodeid, position);

		}

                if (iscontained){
                    containing_filters += 1;
                }
            }

            //Giving the current point a mapping to the shape-function matrix if required.
            if ((update_shapefunction) && (shared_dof_material) && (containing_filters>0) && (pointtype==1)){
                //Update in the case of a material point and shared dof-material points (i.e. MPM)
                micro_node_to_row.emplace((*datapoint)[1], index);

    		//Adjust the weighting on the point if it is contained in elements.
                //Note that this should only happen if a point is on the boundary between 
                //multiple filter domains.

                if (containing_filters > 1){
                   micro_node_elcount.emplace((unsigned int)((*datapoint)[1]+0.5), containing_filters);
                }

            }
            else if ((update_shapefunction) && (!shared_dof_material) && (containing_filters>0) && (pointtype==2)){
                //Update in the case of a dof point and non-shared dof-material points (i.e. FEA)
                micro_node_to_row.emplace((*datapoint)[1], index);

		//Adjust the weighting on the point if it is contained in elements.
                //Note that this should only happen if a point is on the boundary between 
                //multiple filter domains.

                if (containing_filters > 1){
                    micro_node_elcount.emplace((unsigned int)((*datapoint)[1]+0.5), containing_filters);
                }

            }
            //Increment index if it was contained
            if (containing_filters>0){
                index++;
            }
        }

        return 0;
    }

    int construct_micro_displacement_vector_from_positions(const elib::vecOfvec &data, const input_format &mp_format, 
                                                           const input_format &dof_format, const uint_to_vec &reference_coordinates,
                                                           const bool shared_dof_material,
                                                           const unsigned int num_micro_dof, const uint_map &micro_node_to_row,
                                                           elib::vec &micro_displacement_vector){

        /*!
         * Construct the micro displacement vector by computing the difference between the current position and some 
         * reference position.
         * 
         * :param const elib::vecOfvec &data: A collection of micro-scale data
         * :param const input_format &mp_format: The format of a material point dataline
         * :param const input_format &dof_format: The format of a degree of freedom dataline
         * :param const uint_to_vec &reference_coordinates: The reference coordinates for the data
         * :param const bool shared_dof_material: Whether the degrees of freedom are located at the material points or not. 
         * :param const unsigned int num_micro_dof: The number of micro-scale degrees of freedom per node
         * :param const uint_map &micro_node_to_row: The map from the micro node to the row of the shape-function matrix
         * :param elib::vec &micro_displacement_vector: The returned micro-scale displacement vector
         */

        micro_displacement_vector.resize(micro_node_to_row.size()*num_micro_dof);
//        std::cout << "reference coordinates:\n"; print(reference_coordinates);

        std::vector< double > pi;
        bool dof_point;
        int nodetype;
        unsigned int nodeid;
        int gp_result;

        for (auto dataline=data.begin(); dataline!=data.end(); dataline++){
            //Compute the difference between the current coordinates and the 
            //reference
            dof_point = false;
            nodetype = (int)((*dataline)[0]+0.5);
            gp_result = 0;

            //Extract the point information if required
            if ((nodetype==1) && (shared_dof_material)){
                dof_point = true;
                gp_result = get_position(*dataline, mp_format, nodeid, pi);
            }
            else if ((nodetype==2) && (!shared_dof_material)){
                dof_point = true;
                gp_result = get_position(*dataline, dof_format, nodeid, pi);
            }

            if (gp_result > 0){
                return gp_result;
            }

            //Check if the current node is located in the referencecoordinates map (it better be!)
            if (dof_point){
                auto ref = reference_coordinates.find(nodeid);
                if (ref == reference_coordinates.end()){
                    std::cerr << "Error: node " << nodeid << " not found in reference coordinates.\n";
                    std::cerr << "       it is currently required that nodes cannot be deleted.\n";
                    return 1;
                }
                //Assign the value to the micro-displacement (dof) vector
                auto it = micro_node_to_row.find(nodeid);
                if (it != micro_node_to_row.end()){
                    for (unsigned int i=0; i<pi.size(); i++){
                        micro_displacement_vector[num_micro_dof*it->second+i] = pi[i] - ref->second[i];
                    }
                }
            }
        }
        return 0;   
        
    }

    int assign_dof_information_to_filters(const assembly::node_map &nodes, const assembly::element_map &elements,
                                          const uint_map &macro_node_to_col, const unsigned int num_macro_dof,
                                          const std::vector< double > &macro_displacement, filter_map &filters){
        /*!
         * Assign the computed macro-scale degree of freedom values to the nodes
         * 
         */

        //Assign the dof information to the filter and update their current nodal positions.
        for (auto filter=filters.begin(); filter!=filters.end(); filter++){

            //Get the elements of the given type
            auto elementtype=elements.find(filter->second.element_type());
            if (elementtype == elements.end()){
                std::cout << "Error: cant find filter element type in elements\n";
                return 1;
            }

            //Get the element nodal definition
            auto element = elementtype->second.find(filter->second.id());
            if (element == elementtype->second.end()){
                std::cout << "Error: filter not found in element list\n";
                return 1;
            }
                
            //Iterate through the element's nodes
            unsigned int index=0;
//            std::cout << "macro_node_to_col:\n";
//            print(macro_node_to_col);
//            std::cout << "macro_displacement:\n";
//            for (unsigned int i=0; i<8; i++){
//                for (unsigned int j=0; j<num_macro_dof; j++){
//                    std::cout << macro_displacement[num_macro_dof*i + j] << ", ";
//                }
//                std::cout << "\n";
//            }
//            std::cout << "element nodes:\n";
            for (auto n=element->second.begin(); n!=element->second.end(); n++){

                auto dofptr = macro_node_to_col.find(*n);
                if (dofptr == macro_node_to_col.end()){
                    std::cerr << "Error: filter node " << *n << "not found in macro_node_to_col map\n";
                    return 1;
                }

                filter->second.update_dof_values(index, 
                    std::vector< double >(macro_displacement.begin() + num_macro_dof*dofptr->second,
                                          macro_displacement.begin() + num_macro_dof*dofptr->second + num_macro_dof));

                //We assume that the first dof values are the macro displacements
                unsigned int uenp_result = filter->second.update_element_node_position(index);
//                unsigned int uenp_result = filter->second.update_element_node_positions(index,
//                    std::vector< double > (&macro_displacement[num_macro_dof*dofptr->second],
//                                           &macro_displacement[num_macro_dof*dofptr->second+filter->second.dim()])
//                );
                if (uenp_result > 0){
                    return 1;
                }

                index++;
            }
            //filter->second.print();
        }
        return 0;
    }

    int assemble_micro_density(const elib::vecOfvec &data, const input_format mp_format, std::map< unsigned int, double > &micro_density){
        /*!
         * Assemble the micro-density map for use in homogenizing the density
         */
        
        auto idit = mp_format.find("ID");
        auto densityit = mp_format.find("DENSITY");

        if ((idit == mp_format.end()) || (densityit == mp_format.end())){
            std::cout << "MPFORMAT:\n"; print(mp_format);
            std::cout << "Error: density or id not defined in *MPFORMAT\n";
            return 1;
        }

        for (auto dataline = data.begin(); dataline!=data.end(); dataline++){

            if ((unsigned int)((*dataline)[0]+0.5) == 1){
                micro_density.emplace((*dataline)[ idit->second[0] ], (*dataline)[ densityit->second[0] ]);
            }
        }
        return 0;
    }

    int get_position(const elib::vec &dataline, const input_format &format, unsigned int &node_id, elib::vec &position){
        /*!
         * Get the position from the provided dataline given the format
         * 
         * :param const elib::vec &dataline: The line of data as read from the input file
         * :param const input_format &format: The formatting of the line
         * :param unsigned int &node_id: The id associated with the position.
         * :param elib::vec &position: The returned position vector
         */
        
        auto idit = format.find("ID");
        auto pit = format.find("POSITION");

        if ((idit == format.end()) || (pit == format.end())){
            std::cerr << "Error: ID or POSITION not found in format\n";
            return 1;
        }

        node_id = dataline[idit->second[0]];
        position.resize(pit->second[1]);
        for (unsigned int i=0; i<pit->second[1]; i++){
            position[i] = dataline[pit->second[0]+i];
        }
        return 0;
    }

    int populate_reference_coordinates(const elib::vecOfvec &data, const bool shared_dof_material, const input_format &mp_format,
                                       const input_format &dof_format, uint_to_vec &reference_coordinates){
        /*!
         * Populate the reference coordinates map
         * 
         * :param const elib::vecOfvec &data: The dns data to be parsed.
         * :param const input_format &mp_format: The format of material points.
         * :param const input_format &dof_format: The format of degree of freedom points
         * :param uint_to_vec &reference_coordinates: The reference coordinates map
         */
 
        //Populate the reference coordinate map
        reference_coordinates.clear();

        std::vector< double > pi;
        bool dof_point;
        int nodetype;
        unsigned int nodeid;
        int gp_result;

        for (auto dataline=data.begin(); dataline!=data.end(); dataline++){
//            std::vector< double > pi(num_micro_dof);
            dof_point = false;
            gp_result = 0;
            nodetype = (int)((*dataline)[0]+0.5);

            //Extract the point information if required
            if ((nodetype==1) && (shared_dof_material)){
                dof_point = true;
                gp_result = get_position(*dataline, mp_format, nodeid, pi);
            }
            else if ((nodetype==2) && (!shared_dof_material)){
                dof_point = true;
                gp_result = get_position(*dataline, dof_format, nodeid, pi);
            }

            if (gp_result > 0){
                return gp_result;
            }

            //Store the reference coordinates of the point if it is a DOF point
            if (dof_point){
                reference_coordinates.emplace(nodeid, pi);
            }
        }
        return 0;
    }

    int process_timestep_totalLagrangian(const elib::vecOfvec &data, const input_format &mp_format, const input_format &dof_format, 
                                         const assembly::node_map &nodes,
                                         const assembly::element_map &elements, const assembly::qrule_map &qrules,
                                         const bool shared_dof_material,
                                         std::map< unsigned int, unsigned int > &macro_node_to_col,
                                         std::map< unsigned int, unsigned int > &micro_node_to_row,
                                         std::map< unsigned int, unsigned int > &micro_node_elcount,
					 std::map< unsigned int, elib::vec > &reference_coordinates,
                                         overlap::SpMat &shapefunctions, overlap::QRsolver &dof_solver, filter_map &filters,
                                         std::ofstream &output_file,
                                         const unsigned int num_macro_dof, const unsigned int num_micro_dof){
        /*!
        * Process the current timestep using a Total-Lagrangian approach.
        * If filters is empty, it is assumed that this is the first increment and they will be 
        * populated along with the shape-function matrix.
        *
        * Note that, even in a total-Lagrangian framework, we still must re-populate the material 
        * point integrators every timestep but we can leave the shape-function matrix in tact. This 
        * is because we do not, generally, know how the micro-volumes volumetrically deform but we 
        * can choose to use the points originally contained within the filter as the points that 
        * determine the filter's motion. This is what is meant by Total-Lagrangian.
        * 
        * :param const elib::vecOfvec &data: The DNS data at the current timestep.
        * :param const input_format &mp_format: The format of a material point dataline
        * :param const input_format &dof_format: The format of a degree of freedom dataline
        * :param const assembly::node_map &nodes: The nodes of the micromorphic filter.
        * :param const assembly::element_map &elements: The elements of the micromorphic filter.
        * :param const assembly::qrule_map &qrules: The quadrature rules for the elements.
        * :param const bool shared_dof_material: Whether the dof and material information is co-located.
        * :param std::map< unsigned int, unsigned int > &macro_node_to_col: The map from macro nodes to the corresponding 
	*     column of the matrix (normalized by the num_macro_dof)
        * :param std::map< unsigned int, unsigned int > &micro_node_to_row: The map from micro nodes to the corresponding
	*     row of the matrix (normalized by the num_micro_dof)
        * :param std::map< unsigned int, unsigned int > &micro_node_elcount: The number of macro elements (filters) 
	*     containing the given micro node. Only occurs for micro nodes on the boundaries between two filters and 
	*     is only populated for nodes greater than 1.
        * :param  std::map< unsigned int, elib::vec > &reference_coordinates: The reference coordinates used to compute 
	*     the change in position of the dof points.
        * :param overlap::SpMat &shapefunctions: The shapefunction matrix.
        * :param overlap::QRsolver &dof_solver: The degree of freedom solver object.
        * :param filter_map &filters: Return map of the element ids to their corresponding sub-filters.
        * :param std::ofstream &output_file: The output file to write the filter results to
        * :param const unsigned int num_macro_dof: The number of macro-scale degrees of freedom (default 12)
        * :param const unsigned int num_micro_dof: The number of micro-scale degrees of freedom (default 3)
        */

        //Check if the filters have been populated yet. If not, we will re-compute the shape-function matrix.
        bool populated_filters = false;
//        std::vector< unsigned int > macro_node_ids(nodes.size());
//        std::vector< double > u_450(3);
        if (filters.size() == 0){
            populated_filters = true;

            macro_node_to_col.clear();

            unsigned int index = 0;
            for (auto it=nodes.begin(); it!=nodes.end(); it++){
                macro_node_to_col.emplace(it->first, index);
//                macro_node_ids[index] = it->first;
                index++;
            }
        }
	else{
            //Compute the macro displacement
            std::cout << " Computing the macro-displacement\n";
	    std::vector< double > macro_displacement(shapefunctions.cols());
	    std::vector< double > micro_displacement(shapefunctions.rows());

            int cmdvfp_result = construct_micro_displacement_vector_from_positions(data, mp_format, dof_format, reference_coordinates,
                                                                                   shared_dof_material, num_micro_dof,
                                                                                   micro_node_to_row,
                                                                                   micro_displacement);

            if (cmdvfp_result > 0){
                return 1;
            }

	    //Solve for the macro dof
	    Eigen::Map< overlap::EigVec > b(micro_displacement.data(), micro_displacement.size(), 1);
	    Eigen::Map< overlap::EigVec > x(macro_displacement.data(), macro_displacement.size(), 1);
	    x = dof_solver.solve(b);

            std::cout << " Assigning macro-dof values to the filter\n";
            assign_dof_information_to_filters(nodes, elements, macro_node_to_col, num_macro_dof,
                                              macro_displacement, filters);

	}

        //Populate the filters
        int pf_result = populate_filters(data, mp_format, dof_format, nodes, elements, qrules,
                                         populated_filters, shared_dof_material, num_macro_dof,
                                         micro_node_to_row, micro_node_elcount, filters);

        if (pf_result > 0){
            std::cout << "Error in population of filters.\n";
            return 1;
        }

        //Perform the voronoi cell decomposition (maybe should be parallelized)
        std::cout << "Constructing Integrators\n";
        for (auto filter = filters.begin(); filter!=filters.end(); filter++){
            filter->second.construct_integrators();
        }
        
        //Compute the mass and volume properties of the filters
        std::cout << "Computing Mass Properties\n";
        std::map< unsigned int, double > micro_density;
        assemble_micro_density(data, mp_format, micro_density);

        for (auto filter = filters.begin(); filter!=filters.end(); filter++){
            std::cout << "filter: " << filter->first << "\n";
            filter->second.compute_mass_properties(micro_density);
        }

        //Construct the shapefunction matrix if required
        if (populated_filters){
            std::cout << "Computing the Shape-Function matrix\n";
            std::cout << " Formulating the shape-function matrix terms\n";
            std::vector< overlap::T > tripletList;
            for (auto filter = filters.begin(); filter!=filters.end(); filter++){

                //Compute the contribution of the current filter to the shape-function matrix
                const std::vector< unsigned int > *macro_node_ids = filter->second.get_element_global_node_ids();
//                std::cout << "macro_node_ids:\n";
//                for (unsigned int i=0; i<macro_node_ids->size(); i++){
//                    std::cout << (*macro_node_ids)[i] << "\n";
//                }

                filter->second.add_shapefunction_matrix_contribution(macro_node_to_col, micro_node_to_row, *macro_node_ids,
                                                                     micro_node_elcount, num_macro_dof, num_micro_dof, 
                                                                     data.size(), tripletList);
            }

            //Construct the shapefunction matrix
            std::cout << " Constructing the shape-function matrix from " << tripletList.size() << " terms.\n";
            std::cout << "  rows, cols: " << num_micro_dof*micro_node_to_row.size() << ", " << num_macro_dof*macro_node_to_col.size() << "\n";
            shapefunctions = overlap::SpMat(num_micro_dof*micro_node_to_row.size(), num_macro_dof*macro_node_to_col.size());
            shapefunctions.setFromTriplets(tripletList.begin(), tripletList.end());

            //Construct the DOF solver
            std::cout << " Performing QR decomposition\n";
            dof_solver.compute(shapefunctions);
            
            if (dof_solver.info()!=Eigen::Success){
                std::cout << "Error: Failure in QR decomposition.\n";
                return 1;
            }

	    //Populate the reference coordinate map
            populate_reference_coordinates(data, shared_dof_material, mp_format, dof_format, reference_coordinates);
        }

        for (auto filter = filters.begin(); filter!=filters.end(); filter++){
            filter->second.write_to_file(output_file);
        }

        return 0;
    }

    
    int process_timestep(const elib::vecOfvec &data, const input_format &mp_format, 
                         const input_format &dof_format, const assembly::node_map &nodes,
                         const assembly::element_map &elements, const assembly::qrule_map &qrules,
                         const unsigned int mode, const bool shared_dof_material,
                         std::map< unsigned int, unsigned int > &macro_node_to_col,
                         std::map< unsigned int, unsigned int > &micro_node_to_row,
                         std::map< unsigned int, unsigned int > &micro_node_elcount,
			 std::map< unsigned int, elib::vec > &reference_coordinates,
                         overlap::SpMat &shapefunctions, overlap::QRsolver &dof_solver, filter_map &filters,
                         std::ofstream &output_file){
        /*!
        * Process the current timestep and compute the macro-scale stress and deformation quantities
        * 
        * :param const elib::vecOfvec &data: The DNS data at the current timestep.
        * :param const input_format &mp_format: The format of a material point dataline
        * :param const input_format &dof_format: The format of a degree of freedom dataline
        * :param const assembly::node_map &nodes: The nodes of the micromorphic filter.
        * :param const assembly::element_map &elements: The elements of the micromorphic filter.
        * :param const assembly::qrule_map &qrules: The quadrature rules for the elements.
        * :param const unsigned int mode: The mode that the filter is operating in.
        *     Options:
        *         0 Total Lagrangian. If filters is empty, it will be assumed that 
        *         the current timestep is the reference configuration.
        * :param bool shared_dof_material: Flag which indicates if the dof and material information are 
        *                                  co-located.
        * :param std::map< unsigned int, unsigned int > &macro_node_to_col: The map from macro nodes to the corresponding 
	*     column of the matrix (normalized by the num_macro_dof)
        * :param std::map< unsigned int, unsigned int > &micro_node_to_row: The map from micro nodes to the corresponding
	*     row of the matrix (normalized by the num_micro_dof)
        * :param std::map< unsigned int, unsigned int > &micro_node_elcount: The number of macro elements (filters) 
	*     containing the given micro node. Only occurs for micro nodes on the boundaries between two filters and 
	*     is only populated for nodes greater than 1.
        * :param  std::map< unsigned int, elib::vec > &reference_coordinates: The reference coordinates used to compute 
	*     the change in position of the dof points.
        * :param overlap::SpMat &shapefunctions: The shapefunction matrix.
        * :param overlap::QRsolver &dof_solver: The solver for the degrees of freedom.
        * :param filter_map &filters: Return map of the element ids to their corresponding sub-filters.
        * :param std::ofstream &output_file: The file to write the output data to
        */
    
        if (mode == 0){
            return process_timestep_totalLagrangian(data, mp_format, dof_format, nodes, elements, qrules,
                                                    shared_dof_material, macro_node_to_col, micro_node_to_row, 
						    micro_node_elcount, reference_coordinates, 
						    shapefunctions, dof_solver, filters, output_file);
        }
        return 1;
    }

    int print(const std::map< unsigned int, unsigned int > &map){
        /*!
         * print the map to the terminal
         */

         for (auto it=map.begin(); it!=map.end(); it++){
             std::cerr << it->first << ": " << it->second << "\n";
         }
         return 0;
    }

    int print(const uint_to_vec &map){
        /*!
         * print the map to the terminal
         */

         for (auto it=map.begin(); it!=map.end(); it++){
             std::cerr << it->first << ": "; elib::print(it->second);
         }
         return 0;
    }

    int print(const input_format &format){
        /*!
         * print the format to the terminal
         */
        
        for (auto it=format.begin(); it!=format.end(); it++){
            std::cerr << it->first << ": " << it->second[0] << ", " << it->second[1] << "\n";
        }
        return 0;
    }

}

int main(int argc, char **argv){
    /*!
    * Main
    * 
    * Read in an input file and write out the filtered values.
    * 
    * format should be:
    * ./filter input_filename filter_filename output_filename
    */

    std::string input_fn;
    std::string filter_fn;
    std::string output_fn;

    std::ifstream input_file;
    std::ofstream output_file;

    elib::vecOfvec data;

    std::cout << "\n\n";
    std::cout << "###########################\n";
    std::cout << "### MICROMORPHIC FILTER ###\n";
    std::cout << "###########################\n";
    std::cout << "\n";
    std::cout << "  author: Nathan Miller\n";
    std::cout << "  email: nathanm@lanl.gov\n\n";

    int format=1; //Hard-coded to format 1
    int mode=0; //Hard-coded to total Lagrangian
    bool shared_dof_material = 1; //Hard-coded to co-located information

    if (argc != 4){
        std::cout << "argc: " << argc << "\n";
        std::cout << "Error: require three filenames. A dns data filename to read, a filter definition, and a filename to write.\n";
    }
    else{
        input_fn  = argv[1];
	filter_fn = argv[2];
        output_fn = argv[3];

        std::cout << "Opening filter definition file.\n";

        //Open the filter definition file
	assembly::node_map nodes;
	assembly::element_map elements;
	assembly::qrule_map qrules;
        overlap::SpMat shapefunctions;
        overlap::QRsolver dof_solver;
        
        filter::input_format mp_format;
        filter::input_format dof_format;
        filter::filter_map filters;
	filter::uint_map macro_node_to_col;
	filter::uint_map micro_node_to_row;
	filter::uint_map micro_node_elcount;
	std::map< unsigned int, elib::vec > reference_coordinates;

        int connresult = assembly::read_connectivity_data(filter_fn, nodes, elements, qrules);
	if (connresult > 1){
            std::cout << "Error in constructing filter\n";
	    return 1;
	}

        output_file.open(output_fn);

        output_file << "*INPUT_FILE, " << input_fn << "\n";
        output_file << "*FILTER_CONFIGURATION_FILE, " << filter_fn << "\n";

        //Open the dns file
        std::cout << "Opening micro-scale data file.\n";
        int openresult = filter::open_input_file(input_fn, format, input_file);
        if (openresult>0){
            std::cout << "Error in opening the DNS input file\n";
            return 1;
        }
        
        //Read past the header
        int headerresult = filter::read_past_header(input_file, mp_format, dof_format, format);
        if (headerresult>0){
            std::cout << "Error in skipping the header.\n Check the input file format.\n";
            return 1;
        }

        //Read the timeseries
        while (!input_file.eof()){
            data.resize(0);
            int filterresult = filter::read_timestep(input_file, format, output_file, data);
            if (filterresult>0){
                std::cout << "Error reading timestep.\n";
                return 1;
            }

            std::cout << "Initializing filters\n";
            int pt_result = filter::process_timestep(data, mp_format, dof_format, nodes, elements, qrules, 
                                                     mode, shared_dof_material,
                                                     macro_node_to_col, micro_node_to_row, micro_node_elcount, 
						     reference_coordinates,
                                                     shapefunctions, dof_solver, filters, output_file);
            if (pt_result > 0){
                std::cout << "Error in processing timestep\n";
                return 1;
            }
            else{
                std::cout << "Timestep processing successful\n\n";
            }
        }
        output_file.close();
    }
    std::cout << "Processing of input file " << input_fn << " completed.\n";
    std::cout << "Output written to " << output_fn << "\n";
    return 0;
}
