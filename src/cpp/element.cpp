/*!
===============================================================================
|                                 element.cpp                                 |
===============================================================================
| A collection of finite elements which can be used in various projects.      |
===============================================================================
*/

#include "element.h"
#include<iostream>

namespace elib{
    
    Element::Element(vecOfvec _nodes, quadrature_rule _qrule){
        /*!
        The constructor for the element

        :param const vecOfvec &nodes: The global coordinates of the nodes
        :param const quadrature_rule &qrule: The quadrature rule of the element
        */

        nodes = _nodes;
        qrule = _qrule;
        
        bounding_box.resize(2);
        bounding_box[0] = nodes[0];
        bounding_box[1] = nodes[0];

        for (unsigned int n=0; n<nodes.size(); n++){
            for (unsigned int i=0; i<nodes[n].size(); i++){
                bounding_box[0][i] = std::min(bounding_box[0][i], nodes[n][i]);
                bounding_box[1][i] = std::max(bounding_box[1][i], nodes[n][i]);
            }
        }
    }

    void Element::interpolate(const vec &nodal_values, const vec &local_coordinates,
                              double &value){
        /*!
        Interpolate the vector of values to the provided local coordinates.

        :param const vec &nodal_values: The values to be interpolated (nnodes, value at node n)
        :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
        :param double &value: The interpolated value
        */

        vec shape_functions;
        get_shape_functions(local_coordinates, shape_functions);

        value = 0;
        for (unsigned int n=0; n<nodes.size(); n++){
            value += shape_functions[n]*nodal_values[n];
        }
        return;
    }

    void Element::interpolate(const vecOfvec &nodal_values, const vec &local_coordinates,
                              vec &value){
        /*!
        Interpolate the vector of values to the provided local coordinates.

        :param const vecOfvec &nodal_values: The values to be interpolated (nnodes, value at node n)
        :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
        :param vec &value: The interpolated value
        */
        
        vec shape_functions;
        get_shape_functions(local_coordinates, shape_functions);

        value.resize(nodal_values[0].size());
        for (unsigned int i=0; i<value.size(); i++){
            value[i] = 0;
        }

        for (unsigned int n=0; n<nodes.size(); n++){
            for (unsigned int v=0; v<nodal_values[n].size(); v++){
                value[v] += shape_functions[n]*nodal_values[n][v];
            }
        }
        return;
    }

    void Element::get_local_gradient(const vec &nodal_values, const vec &local_coordinates,
		                     vec &value){
        /*!
        Compute the gradient of the vector of values at the provided local coordinates.

        :param const vec &nodal_values: The values to be interpolated (nnodes,)
        :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
        :param vecOfvec &value: The global gradient of the nodal_values
        */

	vecOfvec local_grad_shape_functions;
	get_local_grad_shape_functions(local_coordinates, local_grad_shape_functions);

	value = vec(local_grad_shape_functions[0].size(), 0);

	for (unsigned int n=0; n<nodes.size(); n++){
            for (unsigned int i=0; i<local_grad_shape_functions[n].size(); i++){
                value[i] += nodal_values[n]*local_grad_shape_functions[n][i];
            }
	}
	return;
    }

    void Element::get_local_gradient(const vecOfvec  &nodal_values, const vec &local_coordinates,
                                     vecOfvec &value){
        /*!
        Compute the gradient of the vector of values at the provided local coordinates.

        :param const vecOfvec &nodal_values: The values to be interpolated (nnodes, value at node n)
        :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
        :param vecOfvec &value: The global gradient of the nodal_values
        */

        vecOfvec local_grad_shape_functions;
        get_local_grad_shape_functions(local_coordinates, local_grad_shape_functions);

        value.resize(nodal_values[0].size());
        for (unsigned int v=0; v<value.size(); v++){
            value[v].resize(local_node_coordinates[0].size());
            for (unsigned int w=0; w<value[v].size(); w++){
                value[v][w] = 0.;
            }
        }

        for (unsigned int n=0; n<nodes.size(); n++){
            for (unsigned int i=0; i<nodal_values[n].size(); i++){
                for (unsigned int j=0; j<local_grad_shape_functions[n].size(); j++){
                    value[i][j] += nodal_values[n][i]*local_grad_shape_functions[n][j];
                }
            }
        }
        return;
    }

    void Element::get_global_gradient(const vec  &nodal_values, const vec &local_coordinates,
                                      const vecOfvec &coords, vec &value){
        /*!
        Compute the gradient of the values with respect to the provided coordinates.

        :param const vec &nodal_values: The values to be interpolated (nnodes, value at node n)
        :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
        :param const vecOfvec &coords: The coordinates of the nodes (nnodes, n dim)
        :param vecOfvec &value: The global gradient of the nodal_values
        */

        vec local_gradient;
        vecOfvec dxdxi, dxidx;
        get_local_gradient(nodal_values, local_coordinates, local_gradient);
        get_local_gradient(coords, local_coordinates, dxdxi);
        
        invert(dxdxi, dxidx);

        value.resize(dxidx[0].size());
        for (unsigned int i=0; i<local_gradient.size(); i++){
            value[i] = 0;
            for (unsigned int j=0; j<dxidx.size(); j++){
                    value[i] += local_gradient[j]*dxidx[j][i];
            }
        }

        return;
    }

    void Element::get_global_gradient(const vecOfvec  &nodal_values, const vec &local_coordinates,
                                      const vecOfvec &coords, vecOfvec &value){
        /*!
        Compute the gradient of the values with respect to the provided coordinates.

        :param const vecOfvec &nodal_values: The values to be interpolated (nnodes, value at node n)
        :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
        :param const vecOfvec &coords: The coordinates of the nodes (nnodes, n dim)
        :param vecOfvec &value: The global gradient of the nodal_values
        */

        vecOfvec local_gradient, dxdxi, dxidx;
        get_local_gradient(nodal_values, local_coordinates, local_gradient);
        get_local_gradient(coords, local_coordinates, dxdxi);
        invert(dxdxi, dxidx);

        value.resize(nodal_values[0].size());
        for (unsigned int v=0; v<value.size(); v++){
            value[v].resize(coords[0].size());
            for (unsigned int w=0; w<value[v].size(); w++){
                value[v][w] = 0.;
            }
        }

        for (unsigned int i=0; i<local_gradient.size(); i++){
            for (unsigned int k=0; k<dxidx.size(); k++){
                for (unsigned int j=0; j<dxidx[k].size(); j++){
                    value[i][j] += local_gradient[i][k]*dxidx[k][j];
                }
            }
        }

        return;
    }

    void Element::get_global_gradient(const vec &nodal_values, const vec &local_coordinates,
                                      vec &value){
        /*!
        Compute the gradient of the values with respect to the nodal coordinates.

        :param const vecOfvec &nodal_values: The values to be interpolated (nnodes, value at node n)
        :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
        :param vecOfvec &value: The global gradient of the nodal_values
        */

        get_global_gradient(nodal_values, local_coordinates, nodes, value);
        return;
    }

    void Element::get_global_gradient(const vecOfvec  &nodal_values, const vec &local_coordinates,
                                      vecOfvec &value){
        /*!
        Compute the gradient of the values with respect to the nodal coordinates.

        :param const vecOfvec &nodal_values: The values to be interpolated (nnodes, value at node n)
        :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
        :param vecOfvec &value: The global gradient of the nodal_values
        */

        get_global_gradient(nodal_values, local_coordinates, nodes, value);
        return;
    }

    void Element::get_jacobian(const vec &local_coordinates, const vecOfvec &reference_coordinates,
                               vecOfvec &jacobian){
        /*!
        Compute the jacobian matrix (dxdX) of the element at the given local coordinates.
        :param vec local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
        :param vecOfvec reference_coordinates: The reference coordinates of the element (nnodes, ndim)
        */

        get_global_gradient(nodes, local_coordinates, reference_coordinates, jacobian);
    }

    void Element::compute_local_coordinates(const vec &global_coordinates, vec &local_coordinates,
                                            double tolr, double tola, unsigned int maxiter, unsigned int maxls){
        /*!
        Compute the local coordinates given the global coordinates. Does this via Newton iteration.

        :param const vec &global_coordinates: The global coordinates to compute the local coordinates at (ndim, )
        :param const vec &local_coordinates: The local coordinates (n local dim, )
        :param double tolr: The relative tolerance
        :param double tola: The absolute tolerance
        :param int maxiter: The maximum number of iterations
        */

        // Set the initial iterate
        vec xi(local_node_coordinates[0].size(), 0);

        // Compute the initial result
        vec x(global_coordinates.size(), 0);
        interpolate(nodes, xi, x);

        // Set the initial residual vector
        vec R(global_coordinates.size(), 0);
        double R0, Rp, Rnorm;
        R0 = 0;
        for (unsigned int i=0; i<x.size(); i++){
            R[i] = global_coordinates[i] - x[i];
            R0 += pow(R[i], 2);
        }
        R0 = Rp = Rnorm = sqrt(R0);

        //Set the tolerance
        double tol = tolr*R0 + tola;

        //Begin the Newton iteration
        unsigned int niter = 0;
        unsigned int nls = 0;
        double lambda = 1.;
        vec dxi;
        vecOfvec J;
        while ((niter < maxiter) && (Rnorm > tol)){
            get_local_gradient(nodes, xi, J);
            solve(J, R, dxi);
            for (unsigned int i=0; i<xi.size(); i++){
                xi[i] += dxi[i];
            }
            interpolate(nodes, xi, x);
            Rnorm = 0;
            for (unsigned int i=0; i<x.size(); i++){
               R[i] = global_coordinates[i] - x[i]; 
               Rnorm += pow(R[i], 2);
            }
            Rnorm = sqrt(Rnorm);

            //Line search
            nls = 0;
            while (Rnorm >= Rp){
                std::cout << "in line search\n";

                lambda *= 0.5;

                for (unsigned int i=0; i<x.size(); i++){
                    xi[i] -= dxi[i];
                    dxi[i] *= lambda;
                    xi[i] += dxi[i];
                }

                interpolate(nodes, xi, x);
                Rnorm = 0;
                for (unsigned int i=0; i<x.size(); i++){
                   R[i] = global_coordinates[i] - x[i];
                   Rnorm += pow(R[i], 2);
                }
                Rnorm = sqrt(Rnorm);

                nls += 1;
                if (nls > maxls){
                    break;
                }
            }
            Rp = Rnorm;

            niter += 1;
        }
        if (Rnorm > tol){
            std::cout << "Error in Newton-Raphson solve\n";
            assert(1==-1);
        }
        else{
            local_coordinates = xi;
        }
    }

    bool Element::bounding_box_contains_point(const vec &x){
        /*!
        Determines if a point is contained within the element's bounding box

        :param const vec &x: The point in global coordinates.
        */

        for (unsigned int i=0; i<bounding_box[0].size(); i++){
            if ((bounding_box[0][i]>x[i]) || bounding_box[1][i]<x[i]){
                return false;
            } 
        }
        return true;
    }

    bool Element::contains_point(const vec &x){
        /*!
        Determines if a point is contained within the element.

        Note: This requires a Newton-Raphson solve. If you need the local coordinates you 
              might want to break this function up into its components. A rough check is 
              available in bounding_box_contains_point.

        :param const vec &x: The point in global coordinates.
        */

        vec xi;
        compute_local_coordinates(x, xi);
        return local_point_inside(xi);
    }

    void Hex8::get_shape_functions(const vec &local_coordinates, vec &result){
        /*!
        Compute the shape functions for a Hex8 element.

        :param const vec &local_coordinates: The local coordinates (n local dim, )
        :param vec &result: The vector of shape function values
        */

        result.resize(local_node_coordinates.size());

        for (unsigned int n=0; n<local_node_coordinates.size(); n++){
            result[n] = 0.125*(1 + local_node_coordinates[n][0]*local_coordinates[0])*
                              (1 + local_node_coordinates[n][1]*local_coordinates[1])*
                              (1 + local_node_coordinates[n][2]*local_coordinates[2]);
        }
    }

    void Hex8::get_local_grad_shape_functions(const vec &local_coordinates, vecOfvec &result){
        /*!
        Compute the local gradients of the shape functions for a Hex8 element.

        :param const vec &local_coordinates: The local coordinates (n local dim, )
        :param vecOfvec &result: The gradients of the shape functions w.r.t. the local coordinates (n nodes, n local dim)
        */

        result.resize(local_node_coordinates.size());
        for (unsigned int n=0; n<local_node_coordinates.size(); n++){
            result[n] = {0.125*local_node_coordinates[n][0]*(1 + local_node_coordinates[n][1]*local_coordinates[1])*(1 + local_node_coordinates[n][2]*local_coordinates[2]),
                         0.125*(1 + local_node_coordinates[n][0]*local_coordinates[0])*local_node_coordinates[n][1]*(1 + local_node_coordinates[n][2]*local_coordinates[2]),
                         0.125*(1 + local_node_coordinates[n][0]*local_coordinates[0])*(1 + local_node_coordinates[n][1]*local_coordinates[1])*local_node_coordinates[n][2]};
        }
    }

    bool Hex8::local_point_inside(const vec &local_coordinates){
        /*!

        Determine whether local coordinates are inside of the element or not

        :param const vec &local_coordinates: The local coordinates (n local dim, )
        */

        for (unsigned int i=0; i<local_coordinates.size(); i++){
            if (abs(local_coordinates[i])>1){
                return false;
            }
        }

        return true;
    }

    //Functions

    void invert(const vecOfvec &A, vecOfvec &Ainv){
        /*!
        Invert the matrix A. Should only be used in very special circumstances.

        :param const vecOfvec &A: The matrix to be inverted
        :param vecOvec &Ainv: The inverse
        */

        //Transfer the values to an Eigen matrix.
        Eigen::MatrixXd _A(A.size(), A[0].size());
        for (unsigned int i=0; i<A.size(); i++){
            for (unsigned int j=0; j<A.size(); j++){
                _A(i, j) = A[i][j];
            }
        }

        //Compute the inverse
        Eigen::MatrixXd _Ainv = _A.inverse();

        //Transfer the values
        Ainv.resize(A.size());
        for (unsigned int i=0; i<_Ainv.rows(); i++){
            Ainv[i].resize(_Ainv.cols());
            for (unsigned int j=0; j<_Ainv.cols(); j++){
                Ainv[i][j] = _Ainv(i, j);
            }
        }

    }

    void solve(const vecOfvec &A, const vec &b, vec &x){
        /*!
        Solve an equation of the form Ax = b

        :param const vecOfvec &A: The A matrix
        :param const vec &b: The b matrix
        :param vec &x: The solution vector
        */

        //Transfer the values of A to an Eigen matrix.
        Eigen::MatrixXd _A(A.size(), A[0].size());
        for (unsigned int i=0; i<A.size(); i++){
            for (unsigned int j=0; j<A.size(); j++){
                _A(i, j) = A[i][j];
            }
        }

        //Create the map of A and b so we can use Eigen.
        Eigen::Map<const Eigen::VectorXd> _b(b.data(), b.size(), 1);

        //Resize x
        x.resize(A[0].size());

        //Create the map for x
        Eigen::Map< Eigen::VectorXd> _x(x.data(), x.size(), 1);

        //Solve the matrix equation
        _x = _A.partialPivLu().solve(_b);
    }

    void print(const vec &a){
        /*!
        Print the vector to the terminal
        */

        for (unsigned int i=0; i<a.size(); i++){
            std::cout << a[i] << " ";
        }
        std::cout << "\n";
    }

    void print(const vecOfvec &A){
        /*!
        Print the matrix to the terminal
        */

        for (unsigned int i=0; i<A.size(); i++){
            print(A[i]);
        }
    }

}

