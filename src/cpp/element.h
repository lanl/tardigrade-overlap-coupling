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

    typedef std::vector< double > vec;
    typedef std::vector< vec > vecOfvec;
    typedef std::vector<std::pair< vec, double> > quadrature_rule;

    class Element{
        /*!
        The base finite element class
        */

        public:
            std::string name; //!The name of the element
            vecOfvec nodes; //!The global coordinates of the nodes
            quadrature_rule qrule; //!The quadrature rule of the element
            vecOfvec local_node_coordinates; //!The local coordinates of the nodes
            vecOfvec bounding_box; //!The bounding box of the element

            Element(){}
            Element(const vecOfvec &nodes, const quadrature_rule &qrule);

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

            void get_global_shapefunction_gradient(const vec &local_coordinates, const int& node,
                                                   vec &value);

            void get_jacobian(const vec &local_coordinates, const vecOfvec &reference_coordinates,
                              vecOfvec &jacobian);

            void compute_local_coordinates(const vec &global_coordinates, vec &local_coordinates,
                                           double tolr=1e-9, double tola=1e-9, unsigned int maxiter=20, unsigned int maxls=5);

            virtual void get_shape_functions(const vec &local_coordinates, vec &result){} //Must be over-ridden
            virtual void get_local_grad_shape_functions(const vec &local_coordiantes, vecOfvec &result){} //Must be over-ridden
            virtual bool local_point_inside(const vec &local_coordinates){return false;} //Must be over-ridden

            bool bounding_box_contains_point(const vec &x);

            bool contains_point(const vec &x);
    };

    class Hex8 : public Element{
        /*!
        An 8 noded hex element.
        */

        public:

            Hex8(const vecOfvec &nodes, const quadrature_rule &qrule) : Element(nodes, qrule){
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
            bool local_point_inside(const vec &local_coordinates);
    };

    //Functions

    void invert(const vecOfvec &A, vecOfvec &Ainv);
    void solve(const vecOfvec &A, const vec &b, vec &x);
    void print(const vec &a);
    void print(const vecOfvec &A);
    void print(const quadrature_rule &qrule);
    void print(const Element &element);

    std::unique_ptr<Element> build_element_from_string(const std::string &elname, const vecOfvec &nodes, const quadrature_rule &qrule);
}
#endif
