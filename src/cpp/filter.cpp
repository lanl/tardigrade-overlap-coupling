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

}

int main(int argc, char **argv){
    /*!
    Main

    Read in an input file and write out the filtered values.

    format should be:
    ./filter input_filename output_filename
    */

    std::string input_fn;
    std::string filter_fn;
    std::string output_fn;

    std::ifstream input_file;

    elib::vecOfvec data;

    int format=1; //Hard-coded to format 1

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
        int connresult = assembly::read_connectivity_data(filter_fn, nodes, elements, qrules);
	if (connresult > 1){
            std::cout << "Error in constructing filter\n";
	    return 1;
	}
	assembly::print_qrule_map(qrules);

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
            }
            std::cout << "data.size(): " << data.size() << "\n";
        }
    }
    return 0;
}
