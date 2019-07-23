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
    
    int read_past_header(std::ifstream &file, const int format){
        /*!
        Read past the file header for the given format
    
        :param const std::ifstream &file: The input file
        :param int format: The format of the file
        */
    
        if (format==1){
            return read_past_header_format1(file);
        }
    
        return 1;
    }
    
    int read_past_header_format1(std::ifstream &file){
        /*!
        Read past the header information for format 1
    
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
    
    int read_timestep(std::ifstream &file, const int format, elib::vecOfvec &data){
        /*!
        Read in the information for a single timestep in the provided format.
    
        :param std::ifstream &file: The input file to read
        :param const int format: The input file format
        :param elib::vecOfvec &data: The parsed data. Note: All data is converted to double
        */
    
        if (format == 1){
            return read_timestep_format1(file, data);
        }
    
        return 1;
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
    
        while ((next = line.find(delimiter, last)) != std::string::npos){
            parsed_line.push_back(line.substr(last, next - last));
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
    
    int read_timestep_format1(std::ifstream &file, elib::vecOfvec &data){
        /*!
        Read in the information for a single timestep.
    
        :param std::ifstream &file: The input file to read
        :param elib::vecOfvec &data: The parsed data. Note: All data is converted to double.
        */
    
        double time;
    
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
                      const assembly::qrule_map &qrules, filter_map &filters){
        /*!
        * Build the filters for processing data.
        *
        * :param const assembly::node_map &nodes: The nodes of the micromorphic filter.
        * :param const assembly::element_map &elements: The elements of the micromorphic filter.
        * :param const assembly::qrule_map &qrules: The quadrature rules for the elements.
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


//                //Construct the element
//                elib::Element* _current_element = elib::build_element_from_string(_element_types->first,
//                                                                                                  element_nodes,
//                                                                                                  qrule->second);

                //Initialize the filter
                filters.emplace(_element->first, overlap::MicromorphicFilter(_element_types->first, element_nodes, qrule->second));

            }
        }
        return 0;
    }

    int populate_filters(const elib::vecOfvec &data, const assembly::node_map &nodes,
                        const assembly::element_map &elements, const assembly::qrule_map &qrules,
                        const bool update_shapefunction, const bool shared_dof_material,
                        std::map<unsigned int, unsigned int > &micro_node_to_row,
                        std::map< unsigned int, unsigned int > &micro_node_elcount, filter_map &filters){
        /*!
        * Populate the micromorphic sub-filters using the provided data.
        *
        * :param const elib::vecOfvec &data: The DNS data at the current timestep.
        * :param const assembly::node_map &nodes: The nodes of the micromorphic filter.
        * :param const assembly::element_map &elements: The elements of the micromorphic filter.
        * :param const assembly::qrule_map &qrules: The quadrature rules for the elements.
        * :param const bool update_shapefunction: Flag indicating of the shapefunction matrix is to be updated.
        * :param const bool shared_dof_material: Flag indicating if the degrees of freedom and material quantities are known at 
        *                                       the same points.
        * :param std::map< unsigned int, unsigned int > &micro_node_to_row: The mapping from the micro node number to the 
        *                                                                   row of the shape function matrix.
        * :param std::map< unsigned int, unsigned int > &micro_node_elcount: The number of filters a node is contained within.
        * :param filter_map &filters: Return map of the element ids to their corresponding sub-filters.
        */

        //Construct the filters if they haven't been formed yet
        int bf_result = build_filters(nodes, elements, qrules, filters);
        if (bf_result>0){
            return 1;
        }

        //Erase micro_node_to_row
        if (update_shapefunction){
            micro_node_to_row.clear();
        }

        //This probably could/should be parallelized
        unsigned int index = 0;
        for (auto datapoint=data.begin(); datapoint!=data.end(); datapoint++){
            unsigned int containing_filters = 0;

            //Determine whether the point is a material point or dof point
            int pointtype = (int)((*datapoint)[0]+0.5);

	    //Iterate over the macro-scale filters
            for (auto filter = filters.begin(); filter!=filters.end(); filter++){
                bool iscontained = false;
                if (pointtype==1){
		    iscontained = filter->second.add_micro_material_point((*datapoint)[1],
				                            std::vector< double >((*datapoint).begin()+2, (*datapoint).begin()+5));
		}
		else if (pointtype==2){
                    iscontained = filter->second.add_micro_dof_point((*datapoint)[1],
                                                            std::vector< double >((*datapoint).begin()+2, (*datapoint).begin()+5));
		}

                if (iscontained){
                    containing_filters += 1;
                }
            }

            //Adjust the weighting on the point if it is contained in elements.
            //Note that this should only happen if a point is on the boundary between 
            //multiple filter domains.

            if (containing_filters > 1){
                micro_node_elcount.emplace((unsigned int)((*datapoint)[0]+0.5), containing_filters);
            }

            //Giving the current point a mapping to the shape-function matrix if required.
            if ((update_shapefunction) && (shared_dof_material) && (containing_filters>0) && (pointtype==1)){
                //Update in the case of a material point and shared dof-material points (i.e. MPM)
                micro_node_to_row.emplace((*datapoint)[1], index);
            }
            else if ((update_shapefunction) && (!shared_dof_material) && (containing_filters>0) && (pointtype==2)){
                //Update in the case of a dof point and non-shared dof-material points (i.e. FEA)
                micro_node_to_row.emplace((*datapoint)[1], index);
            }
            //Increment index if it was contained
            if (containing_filters>0){
                index++;
            }
        }

        return 0;
    }

    int process_timestep_totalLagrangian(const elib::vecOfvec &data, const assembly::node_map &nodes,
                                         const assembly::element_map &elements, const assembly::qrule_map &qrules,
                                         const bool shared_dof_material,
                                         overlap::SpMat &shapefunctions, overlap::QRsolver &dof_solver, filter_map &filters,
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
        * :param const assembly::node_map &nodes: The nodes of the micromorphic filter.
        * :param const assembly::element_map &elements: The elements of the micromorphic filter.
        * :param const assembly::qrule_map &qrules: The quadrature rules for the elements.
        * :param const bool shared_dof_material: Whether the dof and material information is co-located.
        * :param overlap::SpMat &shapefunctions: The shapefunction matrix.
        * :param overlap::QRsolver &dof_solver: The degree of freedom solver object.
        * :param filter_map &filters: Return map of the element ids to their corresponding sub-filters.
        * :param const unsigned int num_macro_dof: The number of macro-scale degrees of freedom (default 12)
        * :param const unsigned int num_micro_dof: The number of micro-scale degrees of freedom (default 3)
        */

        //Check if the filters have been populated yet. If not, we will re-compute the shape-function matrix.
        bool populated_filters = false;
        std::vector< unsigned int > macro_node_ids(nodes.size());
        std::map< unsigned int, unsigned int > macro_node_to_col;
        std::map< unsigned int, unsigned int > micro_node_to_row;
        if (filters.size() != elements.size()){
            populated_filters = true;
            unsigned int index = 0;
            for (auto it=nodes.begin(); it!=nodes.end(); it++){
                macro_node_to_col.emplace(it->first, index);
                macro_node_ids[index] = it->first;
                index++;
            }
        }

        //Populate the filters
        std::map< unsigned int, unsigned int > micro_node_elcount;
        int pf_result = populate_filters(data, nodes, elements, qrules,
                                         populated_filters, shared_dof_material,
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
        std::map< unsigned int, elib::vec> micro_position;
        for (auto dataline = data.begin(); dataline!=data.end(); dataline++){
            micro_density.emplace((*dataline)[1], (*dataline)[5]);
        }
        for (auto filter = filters.begin(); filter!=filters.end(); filter++){
            filter->second.compute_mass_properties(micro_density);
        }
        
        //Construct the shapefunction matrix if required
        if (populated_filters){
            std::cout << "Computing the Shape-Function matrix\n";
            std::cout << " Formulating the shape-function matrix terms\n";
            std::vector< overlap::T > tripletList;
            for (auto filter = filters.begin(); filter!=filters.end(); filter++){

                //Compute the contribution of the current filter to the shape-function matrix
                filter->second.add_shapefunction_matrix_contribution(macro_node_to_col, micro_node_to_row, macro_node_ids,
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
            
        }

        return 0;
    }

    
    int process_timestep(const elib::vecOfvec &data, const assembly::node_map &nodes,
                         const assembly::element_map &elements, const assembly::qrule_map &qrules,
                         const unsigned int mode, const bool shared_dof_material,
                         overlap::SpMat &shapefunctions, overlap::QRsolver &dof_solver, filter_map &filters){
        /*!
        * Process the current timestep and compute the macro-scale stress and deformation quantities
        * 
        * :param const elib::vecOfvec &data: The DNS data at the current timestep.
        * :param const assembly::node_map &nodes: The nodes of the micromorphic filter.
        * :param const assembly::element_map &elements: The elements of the micromorphic filter.
        * :param const assembly::qrule_map &qrules: The quadrature rules for the elements.
        * :param const unsigned int mode: The mode that the filter is operating in.
        *     Options:
        *         0 Total Lagrangian. If filters is empty, it will be assumed that 
        *         the current timestep is the reference configuration.
        * :param bool shared_dof_material: Flag which indicates if the dof and material information are 
        *                                  co-located.
        * :param overlap::SpMat &shapefunctions: The shapefunction matrix.
        * :param overlap::QRsolver &dof_solver: The solver for the degrees of freedom.
        * :param filter_map &filters: Return map of the element ids to their corresponding sub-filters.
        */
    
        if (mode == 0){
            return process_timestep_totalLagrangian(data, nodes, elements, qrules,
                                                    shared_dof_material, shapefunctions, dof_solver, filters);
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

    elib::vecOfvec data;

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

        //Open the filter definition file
	assembly::node_map nodes;
	assembly::element_map elements;
	assembly::qrule_map qrules;
        overlap::SpMat shapefunctions;
        overlap::QRsolver dof_solver;
        filter::filter_map filters;

        int connresult = assembly::read_connectivity_data(filter_fn, nodes, elements, qrules);
	if (connresult > 1){
            std::cout << "Error in constructing filter\n";
	    return 1;
	}

        //Open the dns file
        int openresult = filter::open_input_file(input_fn, format, input_file);
        if (openresult>0){
            std::cout << "Error in opening the DNS input file\n";
            return 1;
        }
        
        //Read past the header
        int headerresult = filter::read_past_header(input_file, format);
        if (headerresult>0){
            std::cout << "Error in skipping the header.\n Check the input file format.\n";
            return 1;
        }

        //Read the timeseries
        while (!input_file.eof()){
            data.resize(0);
            int filterresult = filter::read_timestep(input_file, format, data);
            if (filterresult>0){
                std::cout << "Error reading timestep.\n";
                return 1;
            }

            std::cout << "Initializing filters\n";
            int pt_result = filter::process_timestep(data, nodes, elements, qrules, mode, shared_dof_material,
                                                     shapefunctions, dof_solver, filters);
            if (pt_result > 0){
                std::cout << "Error in processing timestep\n";
                return 1;
            }
        }
    }
    return 0;
}
