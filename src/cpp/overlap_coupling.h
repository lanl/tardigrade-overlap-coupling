/*!
===============================================================================
|                              overlap_coupling.h                             |
===============================================================================
| The header file for the overlap coupling library. This library provides the |
| classes, functions, and methods required to compute the required terms for  |
| the micro/meso-scale to macro-scale coupling following the micromorphic     |
| continuum mechanics framework.                                              |
===============================================================================
*/

#ifndef OVERLAP_COUPLING_H
#define OVERLAP_COUPLING_H

#include<stdio.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<math.h>
#include<assert.h>
#include<string.h>

#ifdef __GNUC__
//Avoid warnings from Voro++ root code
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include "voro++.hh"

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __GNUC__
//Avoid warnings from Eigen/Sparse root code
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#endif

#ifdef OVERLAP_LIBCOMPILE
#include<Eigen/Dense>
#include<Eigen/Sparse>
#endif

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif


//#include<Eigen/Dense>

#if CONVEXLIB == QUICKHULL
    #include "3d-quickhull/quickhull.h"
    typedef qh_vertex_t vertex_t;
    typedef qh_mesh_t mesh_t;
#elif CONVEXLIB == CONVHULL_3D
    #include "conhull_3d/convhull_3d.h"
    typedef ch_vertex vertex_t;
    typedef std::pair< std::vector< int >, std::vector< vertex_t > > mesh_t;
#elif CONVEXLIB == AKUUKKA
    #include "quickhull/QuickHull.hpp"
    typedef double FloatType;
    typedef quickhull::Vector3<FloatType> vertex_t;
    typedef quickhull::ConvexHull<FloatType> mesh_t;
#else
    #error CONVEXLIB must be defined. If defined, check that the value is supported.
#endif

#include "element.h"

namespace overlap{

    //!===
    //! | Classes
    //!===

    //Forward declarations
    class OverlapCoupling;
    class ParsedData;
    class MicroPoint;

    //Typedefs
    typedef std::vector< std::vector< double > > vecOfvec;
    typedef std::map< std::vector< double >, std::vector< double > > planeMap;
    typedef std::map< unsigned int, MicroPoint > integrateMap;
    typedef std::vector< std::map< unsigned int, double > > scalar_surface_map;
    typedef std::vector< std::map< unsigned int, std::vector< double > > > vector_surface_map;

    //Ignore Eigen definitions when not being compiled into the library
    #ifdef OVERLAP_LIBCOMPILE
        typedef Eigen::SparseMatrix< FloatType > SpMat;
        typedef Eigen::Matrix< FloatType, Eigen::Dynamic, 1 > EigVec;
        typedef Eigen::SparseVector< FloatType > SpEigVec;
        typedef Eigen::Triplet< FloatType > T;
        typedef Eigen::SparseQR< SpMat, Eigen::COLAMDOrdering<int> > QRsolver;
        typedef Eigen::Ref< const SpMat > SpRef;
    #endif

    class OverlapCoupling{
        /*!
        The primary class for performing the overlap coupling operation. This 
        object will return the integration weights for volume integrals, area  
        weighted normal vectors for surface integrals, and also constructs 
        the shape function matrix which can be used to construct the projectors.
        */

        public:
            OverlapCoupling();
            OverlapCoupling(const vecOfvec &local_coordinates, const vecOfvec &gauss_points);
            void initialize(const vecOfvec &local_coordinates, const vecOfvec &gauss_points);

            //! > Interface to 3D-quickhull
            vertex_t map_vector_to_quickhull(const std::vector< double > &vector) const;
            std::vector< double > map_quickhull_to_vector(const vertex_t &vertex) const;
            void map_vectors_to_quickhull(const std::map< unsigned int, std::vector< FloatType > > &vectors, std::vector< vertex_t > &vertices) const;
            void map_quickhull_to_vectors(const std::vector< vertex_t > &vertices, vecOfvec &vectors) const;
#if CONVEXLIB != AKUUKKA
            void extract_mesh_info(const mesh_t &mesh, vecOfvec &normals, vecOfvec &points) const;
#else
            void extract_mesh_info(mesh_t &mesh, vecOfvec &normals, vecOfvec &points) const;//, vecOfvec &facepoints) const;
#endif
            void compute_node_bounds(const std::map< unsigned int, std::vector< FloatType > > &coordinates, planeMap &planes, 
                std::vector< double > &xbnds, std::vector< double > &ybnds, std::vector< double > &zbnds,
                const double tolr=1e-6, const double tola=1e-6);

            void compute_dns_bounds(const std::map< unsigned int, std::vector< FloatType > > &dns_coordinates, bool use_dns_bounds);

            void compute_weights(const std::map< unsigned int, std::vector< FloatType > > &positions,
                                 std::vector< integrateMap > &points, bool use_dns_bounds=true);

            //! > Interface to defined quantities
            const planeMap* get_element_planes() const;
            const vecOfvec* get_element_bounds() const;
            
            const planeMap* get_dns_planes() const;
            const vecOfvec* get_dns_bounds() const;
            const std::vector< MicroPoint >* get_gauss_domains() const;

            void print_element() const;

        protected:
            vecOfvec local_coordinates;
            vecOfvec gauss_points;

            planeMap element_planes;
            planeMap dns_planes;
            vecOfvec element_bounds;
            vecOfvec dns_bounds;

            std::vector< MicroPoint > gauss_domains;
            std::map< unsigned int, FloatType> boundary_node_volumes;
//            voro::container element_container;
//            voro::container dns_container;

            void compute_element_bounds();
            planeMap compute_unique_planes(const vecOfvec &normals, const vecOfvec &points, const double tolr=1e-6, const double tola=1e-6) const;
            void construct_gauss_domains();


    };

    class ParsedData{
        /*!
        Class which stores data read from an input file. Used for testing.
        */

        public:
            vecOfvec global_nodes;
            vecOfvec local_nodes;
            vecOfvec local_gpts;
            std::vector< unsigned int > node_numbers;
            std::vector< double > volumes;
            std::vector< double > densities;
            vecOfvec coordinates;

            //! > Constructors
            ParsedData(){}

            ParsedData(vecOfvec _global_nodes, vecOfvec _local_nodes, vecOfvec _local_gpts,
                       std::vector< unsigned int > _node_numbers, std::vector< double > _volumes,
                       std::vector< double > _densities, vecOfvec _coordinates){
                global_nodes = _global_nodes;
                local_nodes = _local_nodes;
                local_gpts = _local_gpts;
                node_numbers = _node_numbers;
                volumes = _volumes;
                densities = _densities;
                coordinates = _coordinates;
            }
    };

   class MicroPoint{
        /*!
        Class which stores micro-point information.
        */

        public:
            //! > Constructors
            MicroPoint(){}
            MicroPoint( double ivolume, std::vector< double > icoordinates,
                        std::vector< double > pcoordinates,
                        std::vector< int > iplanes, std::vector< double > iareas,
                        vecOfvec inormals, vecOfvec iface_centroids){
                          /*!
                          Constructor for MicroPoint;

                          :param double ivolume: The volume of the voronoi cell
                          :param std::vector< double > icoordinates: The coordinates of the centroid of the voronoi cell
                          :param std::vector< double > pcoordinates: The coordinates of the particle which defined the cell
                          :param std::vector< int > iplanes: The exterior planes cutting the cell
                          :param std::vector< double > iareas: Areas of the surfaces corresponding to iplanes
                          :param vecOfvec inormals: The normals corresponding to iplanes
                          :param vecOfvec iface_centroids: The centroids of the faces corresponding to iplanes
                          */

                          volume = ivolume;

                          coordinates.reserve(icoordinates.size());
                          for (unsigned int i=0; i<icoordinates.size(); i++){coordinates.push_back(icoordinates[i]);}

                          particle_coordinates.reserve(pcoordinates.size());
                          for (unsigned int i=0; i<pcoordinates.size(); i++){particle_coordinates.push_back(pcoordinates[i]);}

                          planes.reserve(iplanes.size());
                          for (unsigned int i=0; i<iplanes.size(); i++){planes.push_back(iplanes[i]);}

//                          areas.reserve(iareas.size());
//                          for (unsigned int i=0; i<iareas.size(); i++){areas.push_back(iareas[i]);}
//
//                          normals.resize(inormals.size());
//                          for (unsigned int i=0; i<inormals.size(); i++){
//                              normals[i].reserve(inormals[i].size());
//                              for (unsigned int j=0; j<inormals[i].size(); j++){
//                                  normals[i].push_back(inormals[i][j]);
//                              }
//                          }
                          das.resize(inormals.size());
                          for (unsigned int i=0; i<inormals.size(); i++){
                              das[i].reserve(inormals[i].size());
                              for (unsigned int j=0; j<inormals[i].size(); j++){
                                  das[i].push_back(iareas[i]*inormals[i][j]);
                              }
                          }

                          face_centroids.resize(iface_centroids.size());
                          for (unsigned int i=0; i<iface_centroids.size(); i++){
                              face_centroids[i].reserve(iface_centroids[i].size());
                              for (unsigned int j=0; j<iface_centroids[i].size(); j++){
                                  face_centroids[i].push_back(iface_centroids[i][j]);
                              }
                          }

                      }
//            MicroPoint(const MicroPoint &p){
//                /*!
//                Copy constructor
//                :param MicroPoint p: The point to copy from
//                */
//
//                volume = p.volume;
//                coordinates = p.coordinates;
//
//                planes.reserve(p.planes.size());
//                for (unsigned int i=0; i<p.planes.size(); i++){
//                    planes.push_back( p.planes[i]);
//                }
//
//                areas.reserve(p.areas.size());
//                for (unsigned int i=0; i<p.areas.size(); i++){
//                    areas.push_back(p.areas[i]);
//                }
//
//                normals.reserve(p.normals.size());
//                for (unsigned int i=0; i<p.normals.size(); i++){
//                    normals.push_back(p.normals[i]);
//                }
//
//                face_centroids.reserve(p.face_centroids.size());
//                for (unsigned int i=0; i<face_centroids.size(); i++){
//                    face_centroids.push_back(p.face_centroids[i]);
//                }
//            }

            //! > Methods
            void print() const;
            std::vector< double > normal(unsigned int i) const;
            double area(unsigned int i) const;

            //! > Attributes
            double volume;
            std::vector< double > coordinates;
            std::vector< double > particle_coordinates;
            std::vector< int > planes;
//            std::vector< double > areas;
//            vecOfvec normals;
            vecOfvec das;
            vecOfvec face_centroids;
            FloatType weight = 1.;
    };

    class Projector{
        /*!
        A class which projects degrees of freedom and forces from the macro to micro scales and 
        vice versa. This involves constructing the shape function matrix and forming the linear 
        solvers which perform the projectors.

        The current form makes the assumption that no free DNS (micro-scale) nodes are located 
        in the domain of influence of micromorphic (macro-scale) nodes. This allows for a 
        significant reduction in the expense of the projection.
        */
        public:
            Projector();
            Projector(unsigned int _num_macro_dof, unsigned int _num_micro_dof,
                      unsigned int _num_macro_ghost, unsigned int _num_macro_free,
                      unsigned int _num_micro_ghost, unsigned int _num_micro_free);
            Projector(const Projector &p);

            void add_shapefunction_terms(const std::map< unsigned int, unsigned int >* macro_node_to_col_map,
                                         const std::map< unsigned int, unsigned int >* micro_node_to_row_map,
                                         const std::vector< unsigned int > &macro_node_ids,
                                         const std::vector< FloatType > &cg,
                                         const vecOfvec &psis,
                                         const integrateMap &dns_weights,
                                         const std::map< unsigned int, unsigned int>* micro_node_elcount,
                                         bool share_ghost_free_boundary_nodes,
                                         bool macro_elem_is_ghost,
                                         unsigned int num_micro_free);//,
//                                         unsigned int n_macro_dof,
//                                         unsigned int n_micro_dof);

            void initialize(unsigned int _num_macro_dof, unsigned int _num_micro_dof,
                            unsigned int _num_macro_ghost, unsigned int _num_macro_free,
                            unsigned int _num_micro_ghost, unsigned int _num_micro_free);
            void form_shapefunction_matrix(unsigned int nrows, unsigned int ncols);
            int form_BDhQsolver();
            int form_NQDh_PR_transpose_solver();

            void solve_BDhQ(const std::vector< FloatType > &Qvec, std::vector< FloatType > &Dhvec) const;
            void solve_BDhQtranspose(const std::vector< FloatType > &bvec, std::vector< FloatType > &xvec) const;
            void solve_BQhDtranspose(const std::vector< FloatType > &bvec, std::vector< FloatType > &xvec) const;
            void solve_BQhQtranspose(const std::vector< FloatType > &bvec, std::vector< FloatType > &xvec) const;
            void project_dof(const std::vector< double > &D, const std::vector<double > &Q,
                             std::vector< double > &Dh, std::vector< double > &Qh) const;
            int _run_tests(bool solve_for_projectors);

        protected:
            unsigned int num_macro_dof, num_macro_ghost, num_macro_free;
            unsigned int num_micro_dof, num_micro_ghost, num_micro_free;

            //!Variables and functions not to be exposed to non-library code
//            #ifdef OVERLAP_LIBCOMPILE
                //!Variables and attributes
//            Eigen::SparseMatrix< FloatType > shapefunction; //!The shape-function matrix. Formed as a sparse matrix.
            #ifdef OVERLAP_LIBCOMPILE
                SpMat shapefunction; //!The shape-function matrix. Formed as a sparse matrix.
                QRsolver BDhQsolver; //!The solver for the BDhQ projection operation
                QRsolver NQDh_PR_transpose_solver; //!The solver for PR.T x = b of the QR decomposition of NQDh
                std::vector< T > tripletList; //!A list of triplets which will be used to construct the shapefunction matrix
            #endif
    };

    class MicromorphicFilter{
        /*!
        A class which implements the micromorphic filter using other 
        overlap-coupling classes rather than a stand-alone library. The filter 
        tries to implement the approach as generally as possible though it 
        does currently assume the finite-element framework for the filter.
        */

        public:
            MicromorphicFilter(){};
//            MicromorphicFilter(elib::Element* element, bool _shared_dof_material = true);
            MicromorphicFilter(const unsigned int id, const std::string &element_type, 
                               const std::vector< unsigned int > &global_node_ids, const elib::vecOfvec &nodes,
                               const elib::quadrature_rule &qrule, const unsigned int num_macro_dof, bool _shared_dof_material = true, 
                               bool _use_dns_bounds = false);
            
            //Point loading and integration domain construction
            bool add_micro_dof_point(const unsigned int &id, const elib::vec &coordinates, FloatType tol=1e-3);
            bool add_micro_material_point(const unsigned int &id, const elib::vec &coordinates, FloatType tol=1e-1);

            int construct_integrators(bool update_shapefunction_matrix = false);

            //Compute mass properties
            int compute_mass_properties(const std::map< unsigned int, double > &micro_density);
            int get_cg_phis(elib::vecOfvec &cg_phis);

            //Construct the contributions to the shape-function matrix of the filter
            int add_shapefunction_matrix_contribution(const std::map< unsigned int, unsigned int > &macro_node_to_col,
                                                      const std::map< unsigned int, unsigned int > &micro_node_to_row,
                                                      const std::vector< unsigned int > &macro_node_ids,
                                                      const std::map< unsigned int, unsigned int > &micro_node_elcount,
                                                      const unsigned int num_macro_dof, const unsigned int num_micro_dof,
                                                      const unsigned int num_micro_free,
                                                      std::vector< T > &tripletList);
            const unsigned int id();
            const unsigned int dim();
            const std::vector< unsigned int >* get_element_global_node_ids();
            unsigned int get_dns_point_gauss_domain(const unsigned int dns_id);
            const std::vector< FloatType >* get_center_of_mass(const unsigned int &gp_id);

            //Element query/setting tools
            const std::string element_type();
            int update_element_node_position(const unsigned int n);
            int update_element_node_position(const unsigned int n, const elib::vec &displacement);
            int update_element_node_positions(const elib::vecOfvec &displacements);
            int update_dof_values(const unsigned int n, const std::vector< FloatType > &new_dof_values);
            int update_dof_values(const vecOfvec &new_dof_values);

            //Display tools
	    int print(const bool show_microscale_info = false);
            int print_mass_properties();

            //File IO tools
            int write_to_file(std::ofstream &file);

            //Re-initialize
            int clear_microscale();

        protected:
            unsigned int filter_id;
	    std::unique_ptr<elib::Element> element;
            unsigned int filter_dim;
            bool shared_dof_material;
            bool use_dns_bounds;
            OverlapCoupling material_overlap;
            OverlapCoupling dof_overlap;

            //ID number vectors
//            std::vector< unsigned int > dof_id_numbers;
//            std::vector< unsigned int > material_id_numbers;

            //Local coordinates
//            elib::vecOfvec micro_dof_local_coordinates;
//            elib::vecOfvec micro_material_local_coordinates;
            std::map< unsigned int, std::vector< FloatType > > micro_dof_local_coordinates;
            std::map< unsigned int, std::vector< FloatType > > micro_material_local_coordinates;

            //Integrators
            std::vector< integrateMap > dof_weights;
            std::vector< integrateMap > material_weights;
            int construct_dof_point_integrator();
            int construct_material_point_integrator();

            //Mass/geometric properties
            elib::vec volume;
            scalar_surface_map surface_area;
            elib::vec density;
            elib::vecOfvec local_center_of_mass;
            elib::vecOfvec center_of_mass;

            int compute_volume();
            int compute_surface_area();
            int compute_density(const std::map< unsigned int, double > &micro_density);
            int compute_centers_of_mass(const std::map< unsigned int, double > &micro_density);

            //Degree of freedom properties
            vecOfvec dof_values;

            //Stress properties
            int compute_symmetric_microstress(const std::map< unsigned int, std::vector< double > > &micro_cauchy);
            
            vecOfvec symmetric_microstress;

    };
    
    //!===
    //! | Functions
    //!===

    ParsedData read_data_from_file(std::string filename);

    std::vector< std::string > split(std::string s, char delimiter);

    template<typename Out>
    void split(const std::string &s, char delimiter, Out result);

    double dot(const std::vector< double > &a, const std::vector< double > &b);

    std::vector< double > subtract(const std::vector< double > &a, const std::vector< double > &b);

    bool point_on_surface(const std::vector< double > &p, const std::vector< double > &n, const std::vector< double > &a);

    std::vector< double > cross(const std::vector< double > &a, const std::vector< double > &b);

    bool fuzzy_equals(const double a, const double b, const double tolr=1e-6, const double tola=1e-6);

    bool fuzzy_equals(const std::vector< double > &a, const std::vector< double > &b, const double tolr=1e-6, const double tola=1e-6);

    bool compare_vector_directions(const std::vector< double > &v1, const std::vector< double > &v2, const double tolr=1e-6, const double tola=1e-6);

    int normal_from_vertices(const vertex_t &p1, const vertex_t &p2, const vertex_t &p3, std::vector< double > &normal, double tolr=1e-6, double tola=1e-6);

    void print_vertex(const vertex_t &vertex);
    void print_vector(const std::vector< FloatType > &vector);
    void print_matrix(const std::vector< std::vector< FloatType > > &matrix);
    void print_planeMap(const planeMap &planes);
    void print_coordinateMap(const std::map< unsigned int, std::vector< FloatType > > &coordinates);    

    void add_planes_to_container(std::vector< voro::wall_plane > &planes, voro::container &container);

    voro::container* construct_container(const std::map< unsigned int, std::vector< FloatType > > &point_coords,
                                         const vecOfvec &bounds, std::vector< voro::wall_plane> &planes, double expand=1);

    void evaluate_container_information(const std::map< unsigned int, std::vector< FloatType > > &positions,
                                        const std::map< int, std::pair< std::vector< FloatType >, std::vector< FloatType > > > &bounding_faces,
                                        voro::container *container, integrateMap &points, std::map< unsigned int,
                                        FloatType > &boundary_node_volumes);

    void find_face_centroid(const std::vector< int > &face_vertices, const std::vector< double > &vertices, const int &index, std::vector< double > &centroid);

    void map_planes_to_voro(const planeMap &planes, std::vector< voro::wall_plane > &vplanes, int j=0);
    void map_domain_to_voro(const MicroPoint &domains, std::vector< voro::wall_plane > &vplanes);
    void apply_nansons_relation(const std::vector< double > &N, const double &JdA, const vecOfvec &Finv, std::vector< double > &nda);
    void perform_volume_integration( const std::map< unsigned int, double > &values, const std::vector< integrateMap > &weights, std::vector< double > &result);
    void perform_volume_integration( const std::map< unsigned int, std::vector< double > > &values, const std::vector< integrateMap > &weights, std::vector< std::vector< double > > &result);
    void perform_position_weighted_volume_integration( const std::map< unsigned int, double > &values, const std::vector< integrateMap > &weights, std::vector< std::vector< double > > &result);
    void compute_surface_area(const std::vector< integrateMap > &weights, scalar_surface_map &surface_areas);
    void perform_surface_integration( const std::map< unsigned int, double > &values, const std::vector< integrateMap > &weights, std::vector< std::map< unsigned int, double > > &result);
    void perform_surface_integration( const std::map< unsigned int, std::vector< double > > &values, const std::vector< integrateMap > &weights, std::vector< std::map< unsigned int, std::vector< double > > > &result);
    #ifdef OVERLAP_LIBCOMPILE
        void construct_triplet_list(const std::map< unsigned int, unsigned int >* macro_node_to_row_map,
                                    const std::map< unsigned int, unsigned int >* dns_node_to_col_map,
                                    const std::vector< unsigned int > &macro_node_ids,
                                    const std::vector< FloatType > &cg, const vecOfvec &psis,
                                    const integrateMap &dns_weights,
                                    const std::map< unsigned int, unsigned int >* micro_node_elcount,
                                    bool share_ghost_free_boundary_nodes,
                                    bool macro_elem_is_ghost,
                                    unsigned int num_micro_free,
                                    std::vector< T > &tripletList,
                                    unsigned int n_macro_dof=12, unsigned int n_micro_dof=3);
        void solve_for_projector(const SpMat &A, const SpMat &B, SpMat &X);
        SpMat form_sparsematrix(const std::vector< T > &tripletList, unsigned int nrows, unsigned int ncols, const bool &ignore_dup);
        SpMat extract_block(const SpMat &A, unsigned int start_row, unsigned int start_col, unsigned int nrows, unsigned int ncols);
    #endif
//    QRsolver form_solver(SpMat &A);
}

#endif
