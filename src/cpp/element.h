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
#include<map>
#include<Eigen/Dense>
#include<error_tools.h>
#include<vector_tools.h>

namespace elib{

    typedef errorTools::Node errorNode;
    typedef errorNode* errorOut;

    typedef std::vector< unsigned int > uivec;
    typedef std::vector< double > vec;
    typedef std::vector< uivec > vecOfuivec;
    typedef std::vector< vec > vecOfvec;
    typedef std::vector<std::pair< vec, double> > quadrature_rule;

    //Map of currently implemented elements to the number of faces and the number of nodes on each face
    const std::map< std::string, std::pair< unsigned int, std::vector< unsigned int > > > elementRegistry =
        {
            { "Hex8", { 6, { 4, 4, 4, 4, 4, 4 } } },
        };

    const std::map< std::string, unsigned int > elementNameToXDMFType =
        {
            { "Hex8", 9 },
        };

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

            errorOut interpolate(const vec &nodal_values, const vec &local_coordinates,
                                 double &value);

            errorOut interpolate(const vecOfvec &nodal_values, const vec &local_coordinates,
                                 vec &value);

            errorOut get_local_gradient(const vec &nodal_values, const vec &local_coordinates,
	    	                        vec &value);

            errorOut get_local_gradient(const vecOfvec &nodal_values, const vec &local_coordinates,
                                        vecOfvec &value);

            errorOut get_global_gradient(const vec  &nodal_values, const vec &local_coordinates, const vecOfvec &coords,
                                         vec &value);

            errorOut get_global_gradient(const vecOfvec  &nodal_values, const vec &local_coordinates, const vecOfvec &coords,
                                         vecOfvec &value);

            errorOut get_global_gradient(const vec  &nodal_values, const vec &local_coordinates,
                                         vec &value);

            errorOut get_global_gradient(const vecOfvec  &nodal_values, const vec &local_coordinates,
                                         vecOfvec &value);

            errorOut get_global_shapefunction_gradients(const vec &local_coordinates, vecOfvec &dNdx, bool use_reference = false);

            errorOut get_jacobian(const vec &local_coordinates, const vecOfvec &reference_coordinates,
                                  vecOfvec &jacobian);

            errorOut estimate_local_coordinates(const vec &global_coordinates, vec &local_coordinates,
                                                double tolr=1e-9, double tola=1e-9);

            errorOut compute_local_coordinates(const vec &global_coordinates, vec &local_coordinates,
                                              double tolr=1e-9, double tola=1e-9, unsigned int maxiter=20, unsigned int maxls=5);

            virtual errorOut get_shape_functions(const vec &local_coordinates, vec &result){
                return new errorNode( "get_shape_functions", "Not implemented" );
            } //Must be over-ridden
            virtual errorOut get_local_grad_shape_functions(const vec &local_coordinates, vecOfvec &result){
                return new errorNode( "get_local_grad_shape_functions", "Not implemented" );
            } //Must be over-ridden
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

            errorOut get_shape_functions(const vec &local_coordinates, vec &result);
            errorOut get_local_grad_shape_functions(const vec &local_coordinates, vecOfvec &result);
            bool local_point_inside(const vec &local_coordinates, const double tol=1e-8);

    };

    //Functions

    errorOut invert(const vecOfvec &A, vecOfvec &Ainv);
    errorOut solve(const vecOfvec &A, const vec &b, vec &x, int mode=1);
    void print(const vec &a);
    void print(const uivec &a);
    void print(const vecOfvec &A);
    void print(const vecOfuivec &A);
    void print(const quadrature_rule &qrule);
    void print(const Element &element);

    std::unique_ptr<Element> build_element_from_string(const std::string &elname, const std::vector< unsigned int > &global_node_is, 
                                                       const vecOfvec &nodes, const quadrature_rule &qrule);
    void determinant_3x3(const vecOfvec &A, double &d);

    errorOut getPolyhedralCellEquivalentElementType( const unsigned int &index0, const uivec &connectivity,
                                                     unsigned int &XDMFCellType, std::string &elementName,
                                                     unsigned int &deltaIndex );

    errorOut getPolyhedralCellEquivalentElementType( const unsigned int &index0, const uivec &connectivity,
                                                     unsigned int &XDMFCellType, std::string &elementName,
                                                     unsigned int &deltaIndex,
                                                     unsigned int &nFaces, uivec &nNodesOnFace, uivec &nodeIndexArrays );
}
#endif
