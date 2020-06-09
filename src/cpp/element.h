/*!
===============================================================================
|                                  element.h                                  |
===============================================================================
| A collection of finite elements which can be used in various projects.      |
===============================================================================
*/

#ifndef ELEMENT_H
#define ELEMENT_H
#include<vector>
#include<string>
#include<math.h>
#include<memory>
#include<iostream>
#include<Eigen/Dense>

namespace elib{

    typedef std::vector< unsigned int > uivec;
    typedef std::vector< double > vec;
    typedef std::vector< uivec > vecOfuivec;
    typedef std::vector< vec > vecOfvec;
    typedef std::vector<std::pair< vec, double> > quadrature_rule;

    class Element{
        /*!
        The base finite element class
        */

        public:
            std::string name; //!The name of the element
            std::vector< unsigned int > global_node_ids; //!The global id numbers of the nodes
            vecOfvec nodes; //!The global coordinates of the nodes
            vecOfvec reference_nodes; //!The global reference coordinates of the nodes
            quadrature_rule qrule; //!The quadrature rule of the element
            vecOfvec local_node_coordinates; //!The local coordinates of the nodes
            vecOfvec bounding_box; //!The bounding box of the element

            Element(){}
            Element(const std::vector< unsigned int > &global_node_ids, const vecOfvec &nodes, const quadrature_rule &qrule);
	    virtual ~Element() = default;

            void interpolate(const vec &nodal_values, const vec &local_coordinates,
                             double &value);

            void interpolate(const vecOfvec &nodal_values, const vec &local_coordinates,
                             vec &value);

            void get_local_gradient(const vec &nodal_values, const vec &local_coordinates,
			            vec &value);

            void get_local_gradient(const vecOfvec &nodal_values, const vec &local_coordinates,
                                    vecOfvec &value);

            void get_global_gradient(const vec  &nodal_values, const vec &local_coordinates, const vecOfvec &coords,
                                     vec &value);

            void get_global_gradient(const vecOfvec  &nodal_values, const vec &local_coordinates, const vecOfvec &coords,
                                     vecOfvec &value);

            void get_global_gradient(const vec  &nodal_values, const vec &local_coordinates,
                                     vec &value);

            void get_global_gradient(const vecOfvec  &nodal_values, const vec &local_coordinates,
                                     vecOfvec &value);

            void get_global_shapefunction_gradients(const vec &local_coordinates, vecOfvec &dNdx, bool use_reference = false);

            void get_jacobian(const vec &local_coordinates, const vecOfvec &reference_coordinates,
                              vecOfvec &jacobian);

            void estimate_local_coordinates(const vec &global_coordinates, vec &local_coordinates, double tolr=1e-9, double tola=1e-9);

            int compute_local_coordinates(const vec &global_coordinates, vec &local_coordinates,
                                          double tolr=1e-9, double tola=1e-9, unsigned int maxiter=20, unsigned int maxls=5);

            virtual void get_shape_functions(const vec &local_coordinates, vec &result){} //Must be over-ridden
            virtual void get_local_grad_shape_functions(const vec &local_coordinates, vecOfvec &result){} //Must be over-ridden
            virtual bool local_point_inside(const vec &local_coordinates, const double tol=1e-9){return false;} //Must be over-ridden

            bool bounding_box_contains_point(const vec &x);

            bool contains_point(const vec &x);

            int update_node_position(const unsigned int n, const vec &displacement, const bool bounding_box_update=true);

            int update_node_positions(const vecOfvec &displacements);

            int update_bounding_box();

            const std::vector< unsigned int > *get_global_node_ids();
    };

    class Hex8 : public Element{
        /*!
        An 8 noded hex element.
        */

        public:

            Hex8(const std::vector< unsigned int > &global_node_ids, 
                 const vecOfvec &nodes, const quadrature_rule &qrule) : Element(global_node_ids, nodes, qrule){
                name = "Hex8";
                local_node_coordinates = {{-1, -1, -1},
                                          { 1, -1, -1},
                                          { 1,  1, -1},
                                          {-1,  1, -1},
                                          {-1, -1,  1},
                                          { 1, -1,  1},
                                          { 1,  1,  1},
                                          {-1,  1,  1}};
            }

            void get_shape_functions(const vec &local_coordinates, vec &result);
            void get_local_grad_shape_functions(const vec &local_coordinates, vecOfvec &result);
            bool local_point_inside(const vec &local_coordinates, const double tol=1e-8);

    };

    //Functions

    void invert(const vecOfvec &A, vecOfvec &Ainv);
    void solve(const vecOfvec &A, const vec &b, vec &x, int mode=1);
    void print(const vec &a);
    void print(const uivec &a);
    void print(const vecOfvec &A);
    void print(const vecOfuivec &A);
    void print(const quadrature_rule &qrule);
    void print(const Element &element);

    std::unique_ptr<Element> build_element_from_string(const std::string &elname, const std::vector< unsigned int > &global_node_is, 
                                                       const vecOfvec &nodes, const quadrature_rule &qrule);
    void determinant_3x3(const vecOfvec &A, double &d);
}
#endif
