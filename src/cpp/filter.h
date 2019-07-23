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
                         const bool update_shapefunction, const bool shared_dof_material,
                         std::map< unsigned int, unsigned int> &micro_node_to_row,
                         std::map< unsigned int, unsigned int > &micro_node_elcount, filter_map &filters); 
    int process_timestep_totalLagrangian(const elib::vecOfvec &data, const assembly::node_map &nodes,
                                         const assembly::element_map &elements, const assembly::qrule_map &qrules,
                                         const bool shared_dof_material,
                                         overlap::SpMat &shapefunctions, overlap::QRsolver &dof_solver, filter_map &filters,
                                         const unsigned int num_macro_dof=12, const unsigned int num_micro_dof=3);
    int process_timestep(const elib::vecOfvec &data, const assembly::node_map &nodes,
                         const assembly::element_map &elements, const assembly::qrule_map &qrules,
                         const unsigned int mode, const bool shared_dof_material, 
                         overlap::SpMat &shapefunctions, overlap::QRsolver &dof_solver, filter_map &filters);

    int print(const std::map< unsigned int, unsigned int > &map);
}

int main(int argc, char **argv);
#endif
