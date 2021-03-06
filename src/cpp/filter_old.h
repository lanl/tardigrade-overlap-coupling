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
    typedef std::map< unsigned int, overlap::MicromorphicFilter > filter_map;
    typedef std::map< unsigned int, unsigned int > uint_map;
    typedef std::map< unsigned int, elib::vec > uint_to_vec;
    typedef std::map< std::string, std::vector< unsigned int > > input_format;

    //!File IO
    int open_input_file(const std::string &fn, const int format, std::ifstream &file);
    int read_past_header(std::ifstream &file, input_format &mp_format, input_format &dof_format, const int format);
    int set_format(const std::string &line, input_format &format);
    int read_timestep(std::ifstream &file, const int format, std::ofstream &output_file, elib::vecOfvec &data);
    int split_string(const std::string &line, const std::string &delimiter, std::vector< std::string > &parsed_line);
    int parsed_line_to_data(const std::vector< std::string > &parsed_line, elib::vec &line_data);

    //!File IO for a text based file conforming to standard 1
    int open_format1_file(const std::string &fn, std::ifstream &file);
    int read_past_header_format1(std::ifstream &file, input_format &mp_format, input_format &dof_format);
    int read_timestep_format1(std::ifstream &file, double &time, elib::vecOfvec &data);
    int find_current_time_format1(std::ifstream &file, elib::vecOfvec &data, double &time);

    //!Data collection commands
    int assemble_micro_density(const elib::vecOfvec &data, const input_format &mp_format, std::map< unsigned int, double > &micro_density);
    int assemble_micro_stress(const elib::vecOfvec &data, const input_format &mp_format, std::map< unsigned int, double > &micro_stress);

    int get_position(const elib::vec &dataline, const input_format &format, unsigned int &node_id, elib::vec &position);

    int populate_reference_coordinates(const elib::vecOfvec &data, const bool shared_dof_material, const input_format &mp_format,
                                       const input_format &dof_format, uint_to_vec &reference_coordinates);

    //!Processing commands
    int build_filters(const assembly::node_map &nodes, const assembly::element_map elements, 
                      const assembly::qrule_map &qrules, const unsigned int num_macro_dof, filter_map &filters,
                      const bool shared_dof_material);
    int populate_filters(const elib::vecOfvec &data, const input_format &mp_format, const input_format &dof_format, 
                         const assembly::node_map &nodes,
                         const assembly::element_map &elements, const assembly::qrule_map &qrules,
                         const bool update_shapefunction, const bool shared_dof_material,
                         const unsigned int num_macro_dof,
                         uint_map &micro_node_to_row,
                         uint_map &micro_node_elcount, filter_map &filters);
    int construct_micro_displacement_vector_from_positions(const elib::vecOfvec &data, const input_format &mp_format,
                                                           const input_format &dof_format, const uint_to_vec &reference_coordinates,
                                                           const bool shared_dof_material,
                                                           const unsigned int num_micro_dof, const uint_map &micro_node_to_row,
                                                           elib::vec &micro_displacement_vector);
    int assign_dof_information_to_filters(const assembly::node_map &nodes, const assembly::element_map &elements,
                                          const uint_map &macro_node_to_col, const unsigned int num_macro_dof,
                                          const std::vector< double > &macro_displacement, filter_map &filters); 
    int process_timestep_totalLagrangian(const elib::vecOfvec &data, const input_format &mp_format, const input_format &dof_format,
                                         const assembly::node_map &nodes,
                                         const assembly::element_map &elements, const assembly::qrule_map &qrules,
                                         const bool shared_dof_material,
                                         uint_map &macro_node_to_col,
                                         uint_map &micro_node_to_row,
                                         uint_map &micro_node_elcount,
                                         uint_to_vec &reference_coordinates,
                                         overlap::SpMat &shapefunctions, overlap::QRsolver &dof_solver, filter_map &filters,
                                         std::ofstream &output_file,
                                         const unsigned int num_macro_dof=12, const unsigned int num_micro_dof=3);
    int process_timestep(const elib::vecOfvec &data, const input_format &mp_format, const input_format &dof_format, 
                         const assembly::node_map &nodes,
                         const assembly::element_map &elements, const assembly::qrule_map &qrules,
                         const unsigned int mode, const bool shared_dof_material,
                         uint_map &macro_node_to_col,
                         uint_map &micro_node_to_row,
                         uint_map &micro_node_elcount,
                         uint_to_vec &reference_coordinates,
                         overlap::SpMat &shapefunctions, overlap::QRsolver &dof_solver, filter_map &filters,
                         std::ofstream &output_file);

    int print(const uint_map &map);
    int print(const uint_to_vec &map);
    int print(const input_format &format);
}

int main(int argc, char **argv);
#endif
