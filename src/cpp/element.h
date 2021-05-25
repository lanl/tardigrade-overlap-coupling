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

    typedef unsigned int uitype;
    typedef std::vector< uitype > uivec;
    typedef std::vector< double > vec;
    typedef std::vector< uivec > vecOfuivec;
    typedef std::vector< vec > vecOfvec;
    typedef std::vector<std::pair< vec, double> > quadrature_rule;

    //Map of currently implemented elements to the number of faces and the number of nodes on each face
    const std::map< std::string, std::pair< uitype, std::vector< uitype > > > elementRegistry =
        {
            { "Hex8", { 6, { 4, 4, 4, 4, 4, 4 } } },
            { "Quad4", {4, { 2, 2, 2, 2} } }
        };

    const std::map< std::string, uitype > elementNameToXDMFType =
        {
            { "Hex8", 9 },
            { "Quad4", 5 },
        };

    const std::map< uitype, std::string > XDMFTypeToElementName =
        {
            { 9, "Hex8" },
            { 5, "Quad4" },
        };

    const std::map< uitype, uitype > XDMFTypeToNodeCount =
        {
            {  1,  1 }, //Polyvertex
            {  2,  0 }, //Polyline ( special case )
            {  3,  0 }, //Polygon ( special case )
            {  4,  3 }, //Triangle
            {  5,  4 }, //Quadrilateral
            {  6,  4 }, //Tetrahedron
            {  7,  5 }, //Pyramid
            {  8,  6 }, //Wedge
            {  9,  8 }, //Hexahedron
            { 16,  0 }, //Polyhedron ( special case )
        };

    class Element{
        /*!
        The base finite element class
        */

        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            std::string name; //!The name of the element
            std::vector< uitype > global_node_ids; //!The global id numbers of the nodes
            vecOfvec nodes; //!The global coordinates of the nodes
            vecOfvec reference_nodes; //!The global reference coordinates of the nodes
            quadrature_rule qrule; //!The quadrature rule of the element
            vecOfvec local_node_coordinates; //!The local coordinates of the nodes
            vecOfvec bounding_box; //!The bounding box of the element
            vecOfuivec local_surface_node_ids; //!The local ids of the local surface nodes
            vecOfvec local_surface_points; //!The local coordinates of points on the element surface
            vecOfvec local_surface_normals; //!The normal vectors of points on the element surface
            std::vector< quadrature_rule > surface_quadrature_rules; //!The quadrature rules for the surfaces
            uivec surface_fixed_dimension; //!The dimension fixed on each surface
            std::vector< std::string > surface_element_names; //!The names of the elements that make up the surface

            Element(){}
            Element(const std::vector< uitype > &global_node_ids, const vecOfvec &nodes, const quadrature_rule &qrule);
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
                                               double tolr=1e-9, double tola=1e-9, uitype maxiter=20, uitype maxls=5);

            virtual errorOut get_shape_functions(const vec &local_coordinates, vec &result) = 0;
//            {
//                unusedArgs( local_coordinates, result );
//                return new errorNode( "get_shape_functions", "Not implemented" );
//            } //Must be over-ridden
            virtual errorOut get_local_grad_shape_functions(const vec &local_coordinates, vecOfvec &result) = 0;
//            {
//                unusedArgs( local_coordinates, result );
//                return new errorNode( "get_local_grad_shape_functions", "Not implemented" );
//            } //Must be over-ridden
            virtual bool local_point_inside(const vec &local_coordinates, const double tol=1e-9) = 0;
//            {
//                unusedArgs( local_coordinates, tol );
//                return false;
//            } //Must be over-ridden

            bool bounding_box_contains_point(const vec &x);

            bool contains_point( const vec &x, const double tol = 1e-8 );

            int update_node_position(const uitype n, const vec &displacement, const bool bounding_box_update=true);

            int update_node_positions(const vecOfvec &displacements);

            int update_bounding_box();

            const std::vector< uitype > *get_global_node_ids();

            virtual bool point_on_surface( const vec &x, std::vector< uitype > &surf, const double tol );

            virtual bool local_point_on_surface( const vec &ix, std::vector< uitype > &surf, const double tol );

            virtual errorOut transform_local_vector( const vec &xi, const vec &local_vector, vec &global_vector,
                                                     const bool &useCurrent = true );
    };

    const double sqrt3 = std::sqrt( 3. );
    const quadrature_rule Hex8_default_qrule = { { { -1 / sqrt3, -1 / sqrt3, -1 / sqrt3 }, 1. }, 
                                                 { {  1 / sqrt3, -1 / sqrt3, -1 / sqrt3 }, 1. },
                                                 { {  1 / sqrt3,  1 / sqrt3, -1 / sqrt3 }, 1. },
                                                 { { -1 / sqrt3,  1 / sqrt3, -1 / sqrt3 }, 1. },
                                                 { { -1 / sqrt3, -1 / sqrt3,  1 / sqrt3 }, 1. },
                                                 { {  1 / sqrt3, -1 / sqrt3,  1 / sqrt3 }, 1. },
                                                 { {  1 / sqrt3,  1 / sqrt3,  1 / sqrt3 }, 1. },
                                                 { { -1 / sqrt3,  1 / sqrt3,  1 / sqrt3 }, 1. }  };

    const quadrature_rule Quad4_default_qrule = { { { -1 / sqrt3, -1 / sqrt3 }, 1.},
                                                  { {  1 / sqrt3, -1 / sqrt3 }, 1.},
                                                  { {  1 / sqrt3,  1 / sqrt3 }, 1.},
                                                  { { -1 / sqrt3,  1 / sqrt3 }, 1.} };

    const quadrature_rule Bar2_default_qrule = { { { -1 / sqrt3 }, 1. },
                                                 { {  1 / sqrt3 }, 1. } };

    class Hex8 : public Element{
        /*!
        An 8 noded hex element.
        */

        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            Hex8(const std::vector< uitype > &global_node_ids, 
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
                local_surface_points = { { -1,  0,  0 },
                                         {  1,  0,  0 },
                                         {  0, -1,  0 },
                                         {  0,  1,  0 },
                                         {  0,  0, -1 },
                                         {  0,  0,  1 } };
                local_surface_normals = { { -1,  0,  0 },
                                          {  1,  0,  0 },
                                          {  0, -1,  0 },
                                          {  0,  1,  0 },
                                          {  0,  0, -1 },
                                          {  0,  0,  1 } };

                local_surface_node_ids = { { 0, 4, 7, 3 },
                                           { 1, 2, 6, 5 },
                                           { 0, 1, 5, 4 },
                                           { 2, 3, 7, 6 },
                                           { 3, 2, 1, 0 },
                                           { 4, 5, 6, 7 } };

                surface_quadrature_rules = { { { { -1, -1 / sqrt3, -1 / sqrt3 }, 1.},
                                               { { -1,  1 / sqrt3, -1 / sqrt3 }, 1.},
                                               { { -1,  1 / sqrt3,  1 / sqrt3 }, 1.},
                                               { { -1, -1 / sqrt3,  1 / sqrt3 }, 1.} },
                                             { { {  1, -1 / sqrt3, -1 / sqrt3 }, 1.},
                                               { {  1,  1 / sqrt3, -1 / sqrt3 }, 1.},
                                               { {  1,  1 / sqrt3,  1 / sqrt3 }, 1.},
                                               { {  1, -1 / sqrt3,  1 / sqrt3 }, 1.} },
                                             { { { -1 / sqrt3, -1, -1 / sqrt3 }, 1.},
                                               { {  1 / sqrt3, -1, -1 / sqrt3 }, 1.},
                                               { {  1 / sqrt3, -1,  1 / sqrt3 }, 1.},
                                               { { -1 / sqrt3, -1,  1 / sqrt3 }, 1.} },
                                             { { { -1 / sqrt3,  1, -1 / sqrt3 }, 1.},
                                               { {  1 / sqrt3,  1, -1 / sqrt3 }, 1.},
                                               { {  1 / sqrt3,  1,  1 / sqrt3 }, 1.},
                                               { { -1 / sqrt3,  1,  1 / sqrt3 }, 1.} },
                                             { { { -1 / sqrt3, -1 / sqrt3, -1 }, 1.},
                                               { {  1 / sqrt3, -1 / sqrt3, -1 }, 1.},
                                               { {  1 / sqrt3,  1 / sqrt3, -1 }, 1.},
                                               { { -1 / sqrt3,  1 / sqrt3, -1 }, 1.} },
                                             { { { -1 / sqrt3, -1 / sqrt3,  1 }, 1.},
                                               { {  1 / sqrt3, -1 / sqrt3,  1 }, 1.},
                                               { {  1 / sqrt3,  1 / sqrt3,  1 }, 1.},
                                               { { -1 / sqrt3,  1 / sqrt3,  1 }, 1.} } };

            }

            errorOut get_shape_functions(const vec &local_coordinates, vec &result);
            errorOut get_local_grad_shape_functions(const vec &local_coordinates, vecOfvec &result);
            bool local_point_inside(const vec &local_coordinates, const double tol=1e-8);

    };

    class Quad4 : public Element{
        /*!
         * A 4 noded quad element
         */

        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            Quad4(const std::vector< uitype > &global_node_ids,
                  const vecOfvec &nodes, const quadrature_rule &qrule) : Element(global_node_ids, nodes, qrule){
                name = "Quad4";
                local_node_coordinates = { { -1, -1 },
                                           {  1, -1 },
                                           {  1,  1 },
                                           { -1,  1 } };

                local_surface_points = { { -1,  0 },
                                         {  1,  0 },
                                         {  0, -1 },
                                         {  0,  1 } };

                local_surface_normals = { { -1,  0 },
                                          {  1,  0 },
                                          {  0, -1 },
                                          {  0,  1 } };

                local_surface_node_ids = { { 4, 0 },
                                           { 1, 2 },
                                           { 0, 1 },
                                           { 2, 3 } };

                surface_quadrature_rules = { { { { -1, -1 / sqrt3 }, 1 },
                                               { { -1,  1 / sqrt3 }, 1 } },
                                             { { {  1, -1 / sqrt3 }, 1 },
                                               { {  1,  1 / sqrt3 }, 1 } },
                                             { { { -1 / sqrt3, -1 }, 1 },
                                               { {  1 / sqrt3, -1 }, 1 } },
                                             { { { -1 / sqrt3,  1 }, 1 },
                                               { {  1 / sqrt3,  1 }, 1 } } };
 
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

    std::unique_ptr<Element> build_element_from_string(const std::string &elname, const std::vector< uitype > &global_node_is, 
                                                       const vecOfvec &nodes, const quadrature_rule &qrule);
    void determinant_3x3(const vecOfvec &A, double &d);

    errorOut getPolyhedralCellEquivalentElementType( const uitype &index0, const uivec &connectivity,
                                                     uitype &XDMFCellType, std::string &elementName,
                                                     uitype &deltaIndex );

    errorOut getPolyhedralCellEquivalentElementType( const uitype &index0, const uivec &connectivity,
                                                     uitype &XDMFCellType, std::string &elementName,
                                                     uitype &deltaIndex,
                                                     uitype &nFaces, uivec &nNodesOnFace, uivec &nodeIndexArrays );

    const std::map< std::string, quadrature_rule > default_qrules =
        {
            { "Hex8", Hex8_default_qrule },
        };
}
#endif
