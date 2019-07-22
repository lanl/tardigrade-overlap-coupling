/*!
===============================================================================
|                                   filter.h                                  |
===============================================================================
| An implementation of the micromorphic filter using the overlap coupling     |
| library.                                                                    |
===============================================================================
*/

#ifndef MICROMORPHIC_FILTER_H
#define MICROMORPHIC_FILTER_H
#include<iostream>
#include<string>
#include "occonfiguration.h"
#define OVERLAP_LIBCOMPILE
#include "overlap_coupling.h"
#include "assembly.h"

namespace filter{

    //!Typedefs
//    typedef std::map< unsigned int, overlap::MicromorphicFilter > filter_map;
    typedef std::map< unsigned int, overlap::MicromorphicFilter > filter_map;

    //!File IO
    int open_input_file(const std::string &fn, const int format, std::ifstream &file);
    int read_past_header(std::ifstream &file, const int format);
    int read_timestep(std::ifstream &file, const int format, elib::vecOfvec &data);
    int split_string(const std::string &line, const std::string &delimiter, std::vector< std::string > &parsed_line);
    int parsed_line_to_data(const std::vector< std::string > &parsed_line, elib::vec &line_data);

    //!File IO for a text based file conforming to standard 1
    int open_format1_file(const std::string &fn, std::ifstream &file);
    int read_past_header_format1(std::ifstream &file);
    int read_timestep_format1(std::ifstream &file, elib::vecOfvec &data);
    int find_current_time_format1(std::ifstream &file, elib::vecOfvec &data, double &time);

    //!Processing commands
    int build_filters(const assembly::node_map &nodes, const assembly::element_map elements, 
                      const assembly::qrule_map &qrules, filter_map &filters);
    int populate_filters(const elib::vecOfvec &data, const assembly::node_map &nodes,
                         const assembly::element_map &elements, const assembly::qrule_map &qrules,
                         std::map< unsigned int, double > &weights, filter_map &filters);
    int process_timestep_totalLagrangian(const elib::vecOfvec &data, const assembly::node_map &nodes,
                                         const assembly::element_map &elements, const assembly::qrule_map &qrules,
                                         overlap::SpMat &shapefunctions, filter_map &filters);
    int process_timestep(const elib::vecOfvec &data, const assembly::node_map &nodes,
                         const assembly::element_map &elements, const assembly::qrule_map &qrules,
                         const unsigned int mode, overlap::SpMat &shapefunctions, filter_map &filters);
}

int main(int argc, char **argv);
#endif
