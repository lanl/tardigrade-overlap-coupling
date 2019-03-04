/*!
===============================================================================
|                            overlap_coupling.cpp                             |
===============================================================================
| Header file for the overlap coupling classes and functions. These will      |
| compute the required weights and other values for the multi-scale overlap   |
| coupling. The current strategy is to only explicitly support a linear hex   |
| element. This is not considered to be a major restriction as the            |
| micromorphic continuum is relatively costly. Furthermore, hex elements are  |
| generally preferred over tetrahedral elements for mechanics applications.   |
===============================================================================
*/

#include "overlap_coupling.h"

namespace overlap{

    //!===
    //! | Function definitions
    //!===

//   template< typename T1, typename T2 >
    bool fuzzy_compare(const double &a, const double &b, double tolr, double tola){
        /*!
        Compare two values using a tolerance.
        */

        double tol = fmin(tolr*fabs(a) + tola, tolr*fabs(b) + tola);
        return tol > fabs(a - b);
    }


    //!===
    //! | Class definitions
    //!===

    //! > BaseElement

    BaseElement::BaseElement(const std::vector< std::vector< double > > &_global_nodes){
        /*!
        Generic constructor for an isoparametric element.
        */

        //Add the global nodes
        global_coordinates.resize(_global_nodes.size());
        for (unsigned int i=0; i<_global_nodes.size(); i++){
            global_coordinates[i] = Vector(_global_nodes[i]);
        }
    }

    BaseElement::BaseElement(const std::vector< Vector > & _global_nodes){
        /*!
        Generic constructor for an isoparametric element.
        */

        global_coordinates.resize(_global_nodes.size());
        for (unsigned int i=0; i<_global_nodes.size(); i++){
            global_coordinates[i] = _global_nodes[i];
        }
    }

    void BaseElement::interpolate(const std::vector< Vector > &nodal_values, const Vector &Position, Vector &result){
        /*!
        Interpolate a nodally valued function at the provided position.
        */

        if (nodal_values.size() != local_coordinates.size()){
            std::cout << "Error: nodal_values must have the same number of values as there are local coordinates.\n";
            assert(1==0);
        }

        result = shape_function(0, Position)*nodal_values[0];

        for (unsigned int i=1; i<nodal_values.size(); i++){
            result += shape_function(i, Position)*nodal_values[i];
        }
        return;
    }

    void BaseElement::local_gradient(const std::vector< Vector > &nodal_values, const Vector &Position, std::vector< Vector > &result){
        /*!
        Compute the gradient of a nodally valued function w.r.t. the local coordinates at the provided position.

        The first index is the component of v and the second is the local coordinate.
        */

        if (nodal_values.size() != local_coordinates.size()){
            std::cout << "Error: nodal_values must have the same number of values as there are local coordinates.\n";
            assert(1==0);
        }

        result = nodal_values[0].dyadic_product(grad_shape_function(0, Position));
        std::vector< Vector > temp;

        for (unsigned int n=1; n<nodal_values.size(); n++){
            temp = nodal_values[n].dyadic_product(grad_shape_function(n, Position));
            for (unsigned int i=0; i<temp.size(); i++){
                result[i] += temp[i];
            }
        }
    }

    void BaseElement::print() const{
        /*!
        Print output related to the element.
        */

        std::cout << "Global coordinates:\n";
        for (unsigned int i=0; i<global_coordinates.size(); i++){
            std::cout << "node " << i << ": ", global_coordinates[i].print();
        }

        std::cout << "\nlocal coordinates:\n";
        for (unsigned int i=0; i<local_coordinates.size(); i++){
            std::cout << "node " << i << ": ", local_coordinates[i].print();
        }

        std::cout << "\nquadrature points:\n";
        for (unsigned int i=0; i<gauss_points.size(); i++){
            std::cout << "node " << i << ": ", gauss_points[i].print();
        }
    }

    Vector BaseElement::get_local_coordinates(const int &n) const{
        /*!
        Get the local coordinates of node n
        */
        return this->local_coordinates[n];
    }

    //! > Hex8

    void Hex8::initialize(){
        /*!
        Initialize the Hex8 Element
        */

        //! Populate the local coordinates
        local_coordinates.resize(0);
        //! Node 1
        std::vector< double > vec(3, -1);
        local_coordinates.push_back(Vector(vec));

        //! Node 2
        vec[0] = 1;
        local_coordinates.push_back(Vector(vec));

        //! Node 3
        vec[1] = 1;
        local_coordinates.push_back(Vector(vec));

        //! Node 4
        vec[0] = -1;
        local_coordinates.push_back(Vector(vec));

        //! Node 5
        vec[1] = -1;
        vec[2] = 1;
        local_coordinates.push_back(Vector(vec));

        //! Node 6
        vec[0] = 1;
        local_coordinates.push_back(Vector(vec));

        //! Node 7
        vec[1] = 1;
        local_coordinates.push_back(Vector(vec));

        //! Node 8
        vec[0] = -1;
        local_coordinates.push_back(Vector(vec));
        
        //! Populate the gauss point locations.
        double factor = 1./sqrt(3);

        for (unsigned int i=0; i<local_coordinates.size(); i++){
            gauss_points.push_back(local_coordinates[i]*factor);
        }

        //! Populate the gauss points weights
        gauss_weights = std::vector<double>(local_coordinates.size(), 1.);
    }

    double Hex8::shape_function(const int &node, const Vector &Position){
        /*!
        Compute the shape function for the given node at the provided local position.
        */

        return (1 + (local_coordinates[node]*Position)).product()/8;
    }

    Vector Hex8::grad_shape_function(const int &node, const Vector &Position){
        /*!
        Compute the gradient of the shape function for a given node at the provided local position.
        */

        //!Compute the individual components of the shape function product
        Vector tmp = 1 + local_coordinates[node]*Position ;

        //!Assemble the gradient terms
        std::vector< double > vec(3, 0);
        vec[0] = local_coordinates[node](0)*tmp(1)*tmp(2);
        vec[1] = tmp(0)*local_coordinates[node](1)*tmp(2);
        vec[2] = tmp(0)*tmp(1)*local_coordinates[node](2);

        Vector result = Vector(vec)/8;
        return result; 

    }
}
