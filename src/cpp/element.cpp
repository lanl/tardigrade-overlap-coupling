/*!
===============================================================================
|                                 element.cpp                                 |
===============================================================================
| A collection of finite elements which can be used in various projects.      |
===============================================================================
*/

#include "element.h"

namespace elib{
    
    Element::Element(const std::vector< unsigned int > &_global_node_ids, const vecOfvec &_nodes, const quadrature_rule &_qrule){
        /*!
         * The constructor for the element
         * 
         * :param std::vector< unsigned int > global_node_ids: The global id numbers of the element's nodes
         * :param const vecOfvec &nodes: The global coordinates of the nodes
         * :param const quadrature_rule &qrule: The quadrature rule of the element
         */

        global_node_ids = _global_node_ids;
        nodes = _nodes;
        reference_nodes = _nodes;
        qrule = _qrule;
        
        bounding_box.resize(2);
        update_bounding_box();
//        bounding_box[0] = nodes[0];
//        bounding_box[1] = nodes[0];
//
//        for (unsigned int n=0; n<nodes.size(); n++){
//            for (unsigned int i=0; i<nodes[n].size(); i++){
//                bounding_box[0][i] = std::min(bounding_box[0][i], nodes[n][i]);
//                bounding_box[1][i] = std::max(bounding_box[1][i], nodes[n][i]);
//            }
//        }
    }

    errorOut Element::interpolate(const vec &nodal_values, const vec &local_coordinates,
                                  double &value){
        /*!
         * Interpolate the vector of values to the provided local coordinates.
         * 
         * :param const vec &nodal_values: The values to be interpolated (nnodes, value at node n)
         * :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
         * :param double &value: The interpolated value
         */

        vec shape_functions;
        errorOut error = get_shape_functions(local_coordinates, shape_functions);
        if ( error ){
            errorOut result = new errorNode( "interpolate", "Error in get_shape_functions" );
            result->addNext( error );
            return result;
        }

        value = 0;
        for (unsigned int n=0; n<nodes.size(); n++){
            value += shape_functions[n]*nodal_values[n];
        }
        return NULL;
    }

    errorOut Element::interpolate(const vecOfvec &nodal_values, const vec &local_coordinates,
                              vec &value){
        /*!
         * Interpolate the vector of values to the provided local coordinates.
         * 
         * :param const vecOfvec &nodal_values: The values to be interpolated (nnodes, value at node n)
         * :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
         * :param vec &value: The interpolated value
         */
         
        vec shape_functions;
        errorOut error = get_shape_functions(local_coordinates, shape_functions);
        if ( error ){
            errorOut result = new errorNode( "interpolate", "Error in get_shape_functions" );
            result->addNext( error );
            return result;
        }

        value.resize(nodal_values[0].size());
        for (unsigned int i=0; i<value.size(); i++){
            value[i] = 0;
        }

        for (unsigned int n=0; n<nodes.size(); n++){
            for (unsigned int v=0; v<nodal_values[n].size(); v++){
                value[v] += shape_functions[n]*nodal_values[n][v];
            }
        }
        return NULL;
    }

    errorOut Element::get_local_gradient(const vec &nodal_values, const vec &local_coordinates,
	                                 vec &value){
        /*!
         * Compute the gradient of the vector of values at the provided local coordinates.
         * 
         * :param const vec &nodal_values: The values to be interpolated (nnodes,)
         * :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
         * :param vecOfvec &value: The global gradient of the nodal_values
         */

	vecOfvec local_grad_shape_functions;
	errorOut error = get_local_grad_shape_functions(local_coordinates, local_grad_shape_functions);
        if ( error ){
            errorOut result = new errorNode( "get_local_gradient", "Error in get_local_grad_shape_functions" );
            result->addNext( error );
            return result;
        }

	value = vec(local_grad_shape_functions[0].size(), 0);

	for (unsigned int n=0; n<nodes.size(); n++){
            for (unsigned int i=0; i<local_grad_shape_functions[n].size(); i++){
                value[i] += nodal_values[n]*local_grad_shape_functions[n][i];
            }
	}
	return NULL;
    }

    errorOut Element::get_local_gradient(const vecOfvec  &nodal_values, const vec &local_coordinates,
                                         vecOfvec &value){
        /*!
         * Compute the gradient of the vector of values at the provided local coordinates.
         * 
         * :param const vecOfvec &nodal_values: The values to be interpolated (nnodes, value at node n)
         * :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
         * :param vecOfvec &value: The global gradient of the nodal_values
         */

        vecOfvec local_grad_shape_functions;
        errorOut error = get_local_grad_shape_functions(local_coordinates, local_grad_shape_functions);
        if ( error ){
            errorOut result = new errorNode( "get_local_gradient", "Error in get_local_grad_shape_functions" );
            result->addNext( error );
            return result;
        }

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
        return NULL;
    }

    errorOut Element::get_global_gradient(const vec  &nodal_values, const vec &local_coordinates,
                                          const vecOfvec &coords, vec &value){
        /*!
         * Compute the gradient of the values with respect to the provided coordinates.
         * 
         * :param const vec &nodal_values: The values to be interpolated (nnodes, value at node n)
         * :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
         * :param const vecOfvec &coords: The coordinates of the nodes (nnodes, n dim)
         * :param vecOfvec &value: The global gradient of the nodal_values
         */

        vec local_gradient;
        vecOfvec dxdxi, dxidx;
        errorOut error = get_local_gradient(nodal_values, local_coordinates, local_gradient);
        if ( error ){
            errorOut result = new errorNode( "get_global_gradient", "Error in getting the local gradient of the reference coordinates" );
            result->addNext( error );
            return result;
        }

        error = get_local_gradient(coords, local_coordinates, dxdxi);
        if ( error ){
            errorOut result = new errorNode( "get_global_gradient", "Error in getting the local gradient of the current coordinates" );
            result->addNext( error );
            return result;
        }
        
        error = invert(dxdxi, dxidx);
        if ( error ){
            errorOut result = new errorNode( "get_global_gradient",
                                             "Error in computing the inverse of the current local gradient of the current coordinates" );
            result->addNext( error );
            return result;
        }

        value.resize(dxidx[0].size());
        for (unsigned int i=0; i<local_gradient.size(); i++){
            value[i] = 0;
            for (unsigned int j=0; j<dxidx.size(); j++){
                    value[i] += local_gradient[j]*dxidx[j][i];
            }
        }

        return NULL;
    }

    errorOut Element::get_global_gradient(const vecOfvec  &nodal_values, const vec &local_coordinates,
                                          const vecOfvec &coords, vecOfvec &value){
        /*!
         * Compute the gradient of the values with respect to the provided coordinates.
         * 
         * :param const vecOfvec &nodal_values: The values to be interpolated (nnodes, value at node n)
         * :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
         * :param const vecOfvec &coords: The coordinates of the nodes (nnodes, n dim)
         * :param vecOfvec &value: The global gradient of the nodal_values
         */

        vecOfvec local_gradient, dxdxi, dxidx;
        errorOut error = get_local_gradient(nodal_values, local_coordinates, local_gradient);
        if ( error ){
            errorOut result = new errorNode( "get_global_gradient", "Error in getting the local gradient of the reference coordinates" );
            result->addNext( error );
            return result;
        }

        error = get_local_gradient(coords, local_coordinates, dxdxi);
        if ( error ){
            errorOut result = new errorNode( "get_global_gradient", "Error in getting the local gradient of the current coordinates" );
            result->addNext( error );
            return result;
        }

        error = invert(dxdxi, dxidx);
        if ( error ){
            errorOut result = new errorNode( "get_global_gradient",
                                             "Error in computing the inverse of the current local gradient of the current coordinates" );
            result->addNext( error );
            return result;
        }

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

        return NULL;
    }

    errorOut Element::get_global_shapefunction_gradients(const vec &local_coordinates, vecOfvec &dNdx, bool use_reference){
        /*!
         * Compute the gradient of the shape functions w.r.t. the global coordinates.
         * 
         * :param const vec &local_coordinates: The local coordinates at which to 
         *     compute the shape function gradients.
         * :param vecOfvec &dNdx: The gradients of the shape functions in the global
         *     reference frame.
         * :param bool &use_reference: Boolean indicating if the reference coordinates should be used
         */

        vecOfvec dNdxi, dxdxi, dxidx;
        errorOut error = get_local_grad_shape_functions(local_coordinates, dNdxi);
        if ( error ){
            errorOut result = new errorNode( "get_global_gradient",
                                             "Error in computing the inverse of the current local gradient of the shape functions" );
            result->addNext( error );
            return result;
        }

        if (use_reference){
            error = get_local_gradient(reference_nodes, local_coordinates, dxdxi);
            if ( error ){
                errorOut result = new errorNode( "get_global_gradient",
                                                 "Error in computing the inverse of the current local gradient of the shape functions w.r.t. the reference configuration" );
                result->addNext( error );
                return result;
            }
        }

        else{
            error = get_local_gradient(nodes, local_coordinates, dxdxi);
            if ( error ){
                errorOut result = new errorNode( "get_global_gradient",
                                                 "Error in computing the inverse of the current local gradient of the shape functions w.r.t. the current configuration" );
                result->addNext( error );
                return result;
            }
        }

        error = invert(dxdxi, dxidx);
        if ( error ){
            errorOut result = new errorNode( "get_global_gradient",
                                             "Error in computing the inverse of the current local gradient of the current coordinates" );
            result->addNext( error );
            return result;
        }

        dNdx.resize(nodes.size());
        for (unsigned int n=0; n<nodes.size(); n++){
            dNdx[n] = std::vector< double >(nodes[n].size(), 0);
            for (unsigned int i=0; i<nodes[n].size(); i++){
                for (unsigned int j=0; j<dxidx.size(); j++){
                    dNdx[n][i] += dNdxi[n][j]*dxidx[j][i];
                }
            }
        }
        return NULL; 
    }

    errorOut Element::get_global_gradient(const vec &nodal_values, const vec &local_coordinates,
                                          vec &value){
        /*!
         * Compute the gradient of the values with respect to the nodal coordinates.
         * 
         * :param const vecOfvec &nodal_values: The values to be interpolated (nnodes, value at node n)
         * :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
         * :param vecOfvec &value: The global gradient of the nodal_values
         */

        return get_global_gradient(nodal_values, local_coordinates, nodes, value);
    }

    errorOut Element::get_global_gradient(const vecOfvec  &nodal_values, const vec &local_coordinates,
                                          vecOfvec &value){
        /*!
         * Compute the gradient of the values with respect to the nodal coordinates.
         * 
         * :param const vecOfvec &nodal_values: The values to be interpolated (nnodes, value at node n)
         * :param const vec &local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
         * :param vecOfvec &value: The global gradient of the nodal_values
         */

        return get_global_gradient(nodal_values, local_coordinates, nodes, value);
    }

    errorOut Element::get_jacobian(const vec &local_coordinates, const vecOfvec &reference_coordinates,
                                   vecOfvec &jacobian){
        /*!
         * Compute the jacobian matrix (dxdX) of the element at the given local coordinates.
         * :param vec local_coordinates: The local coordinates at which to interpolate the values. (nnodes, n local dim)
         * :param vecOfvec reference_coordinates: The reference coordinates of the element (nnodes, ndim)
         */

        return get_global_gradient(nodes, local_coordinates, reference_coordinates, jacobian);
    }

    errorOut Element::estimate_local_coordinates(const vec &global_coordinates, vec &local_coordinates, double tolr, double tola){
        /*!
         * Estimate the local coordinates of a globally defined node
         * 
         * :param const vec &global_coordinates: The global coordinates of the node in question
         * :param vec &local_coordinates: The estimate of the local coordinates
         * :param double tolr: The relative tolerance
         * :param double tola: The absolute tolerance
         */

        //Initialize the distances vector        
        vec distance(nodes.size(), 0);

        //Initialize the sum distances value
        double sum_distance = 0;
        double sum_inv_distance = 0;

        //Compute the distances from the point to the nodes
        unsigned int index=0;
        for (auto node = nodes.begin(); node != nodes.end(); node++){

            //Make sure the current node and global point have the same global dimension
            if ((*node).size() != global_coordinates.size()){
                return new errorNode( "estimate_local_coordinates", "Error: point and node have different global dimensions" );
            }

            //Compute the distance
            for (unsigned int i=0; i<(*node).size(); i++){
                distance[index] += pow((*node)[i] - global_coordinates[i], 2);
            }
            distance[index] = std::sqrt(distance[index]);

            //Update the sum of the distances
            sum_distance += distance[index];

            index++;
            
        }

        //Make sure none of the distances are too small
        index = 0;

        //Set the tolerance
        double tol = tolr*sum_distance + tola;
        for (auto d=distance.begin(); d!=distance.end(); d++){
            if ((*d) < tol){
                local_coordinates = local_node_coordinates[index];
                return NULL;
            }
            sum_inv_distance += 1/(*d);
            index++;
        }

        //Use the distances and the sum of the distances to estimate the local coordinates
        index = 0;
        local_coordinates = vec(local_node_coordinates[0].size(), 0);
        for (auto node = local_node_coordinates.begin(); node != local_node_coordinates.end(); node++){
            if (local_coordinates.size() != (*node).size()){
                return new errorNode( "estimate_local_coordinates", "Error: local node coordinates have different local dimensions" );
            }

            //Compute the weighted distance
            for (unsigned int i=0; i<(*node).size(); i++){
                local_coordinates[i] += (*node)[i]*(1/distance[index])/sum_inv_distance;
            }

            //Increment the index
            index++;
        }
        return NULL;
    }

    errorOut Element::compute_local_coordinates(const vec &global_coordinates, vec &local_coordinates,
                                                double tolr, double tola, unsigned int maxiter, unsigned int maxls){
        /*!
         * Compute the local coordinates given the global coordinates. Does this via Newton iteration.
         * 
         * :param const vec &global_coordinates: The global coordinates to compute the local coordinates at (ndim, )
         * :param const vec &local_coordinates: The local coordinates (n local dim, )
         * :param double tolr: The relative tolerance
         * :param double tola: The absolute tolerance
         * :param int maxiter: The maximum number of iterations
         */

        // Set the initial iterate
        vec xi(local_node_coordinates[0].size(), 0);
        errorOut error = estimate_local_coordinates(global_coordinates, xi);
        if ( error ){
            errorOut result = new errorNode( "compute_local_coordinates", "Error in estimation of the local coordinates" );
            result->addNext( error );
            return result;
        }

        // Compute the initial result
        vec x(global_coordinates.size(), 0);

        error = interpolate(nodes, xi, x);
        if ( error ){
            errorOut result = new errorNode( "compute_local_coordinates", "Error in interpolation of the nodes in initialization" );
            result->addNext( error );
            return result;
        }

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
            error = get_local_gradient(nodes, xi, J);
            if ( error ){
                errorOut result = new errorNode( "compute_local_coordinates", "Error in computation of the local gradient in non-linear solve" );
                result->addNext( error );
                return result;
            }

/*            //Estimate the gradient using finite differences
            vecOfvec Jest(R.size(), vec(R.size(), 0));
            double eps = 1e-6;
            vec delta(xi.size(), 0);
            for (unsigned int i=0; i<xi.size(); i++){
                delta = xi;
                delta[i] = xi[i]*(1 + eps);
                interpolate(nodes, delta, x);
                for (unsigned int j=0; j<x.size(); j++){
                    Jest[j][i] = ((global_coordinates[j] - x[j]) - R[j])/(xi[i]*eps);
                }
                
            }

            std::cout << "################################\n";
            std::cout << "Jest:\n"; print(Jest);
            std::cout << "J:\n"; print(J);
*/
            error = solve(J, R, dxi, 2);
            if ( error ){
                errorOut result = new errorNode( "compute_local_coordinates", "Error in non-linear solve" );
                result->addNext( error );
                return result;
            }

            for (unsigned int i=0; i<xi.size(); i++){
                xi[i] += dxi[i];
            }
            error = interpolate(nodes, xi, x);
            if ( error ){
                errorOut result = new errorNode( "compute_local_coordinates", "Error in interpolation in non-linear solve" );
                result->addNext( error );
                return result;
            }

            Rnorm = 0;
            for (unsigned int i=0; i<x.size(); i++){
               R[i] = global_coordinates[i] - x[i]; 
               Rnorm += pow(R[i], 2);
            }
            Rnorm = sqrt(Rnorm);

            //Line search
            nls = 0;
            lambda = 1.;
            while (Rnorm >= Rp){
//                std::cout << "global coordinates: "; print(global_coordinates);
//                std::cout << "nodes:\n"; print(nodes);
//                std::cout << "x: "; print(x);
//                std::cout << "Rp: " << Rp << "\n";
//                std::cout << "Rnorm: " << Rnorm << "\n";
//                std::cout << "lambda: " << lambda << "\n";
//                std::cout << "xi: "; print(xi);
//                std::cout << "dxi: "; print(dxi);

//                std::cout << "in line search iteration " << nls << "\n";
//                std::cout << " global coordinates: "; print(global_coordinates);
//                std::cout << " nodes:\n"; print(nodes);
//                std::cout << " x: "; print(x);
//                std::cout << " Rp: " << Rp << "\n";
//                std::cout << " Rnorm: " << Rnorm << "\n";
//                std::cout << " lambda: " << lambda << "\n";
//                std::cout << " xi: "; print(xi);
//                std::cout << " dxi: "; print(dxi);

                lambda *= 0.5;

                for (unsigned int i=0; i<x.size(); i++){
                    xi[i] -= dxi[i];
                    dxi[i] *= lambda;
                    xi[i] += dxi[i];
                }

                error = interpolate(nodes, xi, x);
                if ( error ){
                    errorOut result = new errorNode( "compute_local_coordinates", "Error in interpolation in line search" );
                    result->addNext( error );
                    return result;
                }

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

//                std::cout << "Rnorm new: " << Rnorm << "\n";
//                std::cout << "xi: "; print(xi);
//                std::cout << "dxi: "; print(dxi);
//                std::cout << "x: "; print(x);
//                assert(1==0);
            }

            if (nls > maxls){
                return new errorNode( "compute_local_coordinates", "Failure in line search" );
            }

            Rp = Rnorm;

            niter += 1;
        }
        if (Rnorm > tol){
            return new errorNode( "compute_local_coordinates", "Newton-raphson did not converge" );
        }
        else{
            local_coordinates = xi;
        }
        return NULL;
    }

    bool Element::bounding_box_contains_point(const vec &x){
        /*!
         * Determines if a point is contained within the element's bounding box
         * 
         * :param const vec &x: The point in global coordinates.
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
         * Determines if a point is contained within the element.
         * 
         * Note: This requires a Newton-Raphson solve. If you need the local coordinates you 
         *       might want to break this function up into its components. A rough check is 
         *       available in bounding_box_contains_point.
         *
         * :param const vec &x: The point in global coordinates.
         */

        vec xi;
        errorOut error = compute_local_coordinates(x, xi);

        //We assume that if the local coordinates cannot be computed the point must be outside of the element
        if ( error ){
            return false;
        }
        return local_point_inside(xi);
    }

    int Element::update_node_position(const unsigned int n, const vec &displacement, const bool bounding_box_update){
        /*!
         * Update the nodal position of node n in the element
         * 
         * :param const unsigned int n: The local node number
         * :param const elib::vec &displacement: The displacement of the node from the reference state
         * :param const bool bounding_box_update: Boolean to indicate if the bounding box update should be calculated.
         */
        
        if (reference_nodes[n].size() != displacement.size()){
            std::cerr << "Error: local node " << n << " has a dimension of " << reference_nodes[n].size() << ".\n";
            std::cerr << "       the nodal displacement has a dimension of " << displacement.size() << ".\n";
            return 1;
        }
        for (unsigned int i=0; i<reference_nodes[n].size(); i++){
            nodes[n][i] = reference_nodes[n][i] + displacement[i];
        }

        if (bounding_box_update){
            update_bounding_box();
        }

        return 0;
    }

    int Element::update_node_positions(const vecOfvec &displacements){
        /*!
         * Update the nodal positions of the elements to reflect a movement of the filter.
         * 
         * :param elib::vecOfvec &displacements: The displacements at the nodes from the reference coordinates
         */
        
        if (nodes.size() != displacements.size()){
            std::cerr << "Error: " << displacements.size() << " nodal displacements provided to an element which has " << reference_nodes.size() << "nodes.\n";
            return 1;
        }
        for (unsigned int n=0; n<nodes.size(); n++){
            int uenp_result = update_node_position(n, displacements[n], false);
            if (uenp_result > 0){
                return uenp_result;
            }
        }

        update_bounding_box();
        return 0;
    }

    int Element::update_bounding_box(){
        /*!
         * Update the element's bounding box.
         */

        
        bounding_box[0] = nodes[0];
        bounding_box[1] = nodes[0];

        for (unsigned int n=1; n<nodes.size(); n++){
            for (unsigned int i=0; i<nodes[n].size(); i++){
                bounding_box[0][i] = std::min(bounding_box[0][i], nodes[n][i]);
                bounding_box[1][i] = std::max(bounding_box[1][i], nodes[n][i]);
            }
        }
        return 0;
    }

    const std::vector< unsigned int >* Element::get_global_node_ids(){
        /*!
         * Return a constant pointer to the global node ids
         */
        return &global_node_ids;
    }

    errorOut Hex8::get_shape_functions(const vec &local_coordinates, vec &result){
        /*!
         * Compute the shape functions for a Hex8 element.
         * 
         * :param const vec &local_coordinates: The local coordinates (n local dim, )
         * :param vec &result: The vector of shape function values
         */

        result.resize(local_node_coordinates.size());

        for (unsigned int n=0; n<local_node_coordinates.size(); n++){
            result[n] = 0.125*(1 + local_node_coordinates[n][0]*local_coordinates[0])*
                              (1 + local_node_coordinates[n][1]*local_coordinates[1])*
                              (1 + local_node_coordinates[n][2]*local_coordinates[2]);
        }
        return NULL;
    }

    errorOut Hex8::get_local_grad_shape_functions(const vec &local_coordinates, vecOfvec &result){
        /*!
         * Compute the local gradients of the shape functions for a Hex8 element.
         *
         * :param const vec &local_coordinates: The local coordinates (n local dim, )
         * :param vecOfvec &result: The gradients of the shape functions w.r.t. the local coordinates (n nodes, n local dim)
         */

        result.resize(local_node_coordinates.size());
        for (unsigned int n=0; n<local_node_coordinates.size(); n++){
            result[n] = {0.125*local_node_coordinates[n][0]*(1 + local_node_coordinates[n][1]*local_coordinates[1])*(1 + local_node_coordinates[n][2]*local_coordinates[2]),
                         0.125*(1 + local_node_coordinates[n][0]*local_coordinates[0])*local_node_coordinates[n][1]*(1 + local_node_coordinates[n][2]*local_coordinates[2]),
                         0.125*(1 + local_node_coordinates[n][0]*local_coordinates[0])*(1 + local_node_coordinates[n][1]*local_coordinates[1])*local_node_coordinates[n][2]};
        }

        return NULL;
    }

    bool Hex8::local_point_inside(const vec &local_coordinates, const double tol){
        /*!
         * Determine whether local coordinates are inside of the element or not
         * 
         * :param const vec &local_coordinates: The local coordinates (n local dim, )
         * :param double tol: The tolerance
         */

        for (unsigned int i=0; i<local_coordinates.size(); i++){
            if ((abs(local_coordinates[i]) - 1)>tol){
                return false;
            }
        }

        return true;
    }

    //Functions

    errorOut invert(const vecOfvec &A, vecOfvec &Ainv){
        /*!
         * Invert the matrix A. Should only be used in very special circumstances.
         * 
         * :param const vecOfvec &A: The matrix to be inverted
         * :param vecOvec &Ainv: The inverse
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

        return NULL;

    }

    errorOut solve(const vecOfvec &A, const vec &b, vec &x, int mode){
        /*!
         * Solve an equation of the form Ax = b
         *
         * :param const vecOfvec &A: The A matrix
         * :param const vec &b: The b matrix
         * :param vec &x: The solution vector
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
        if (mode==1){
            _x = _A.partialPivLu().solve(_b);
        }
        else if (mode==2){
            _x = _A.fullPivLu().solve(_b);
        }
        else if (mode==3){
            _x = _A.colPivHouseholderQr().solve(_b);
        }
        else{
            _x = _A.fullPivHouseholderQr().solve(_b);
        }
        return NULL;
    }

    void print(const uivec &a){
        /*!
        * Print the vector to the terminal
        */

        for (unsigned int i=0; i<a.size(); i++){
            std::cout << a[i] << " ";
        }
        std::cout << "\n";
    }

    void print(const vec &a){
        /*!
        * Print the vector to the terminal
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

    void print(const vecOfuivec &A){
        /*!
        * Print the matrix to the terminal
        */

        for (unsigned int i=0; i<A.size(); i++){
            print(A[i]);
        }
    }


    void print(const quadrature_rule &qrule){
        /*!
        * Print the quadrature rule to the terminal.
        */

        for (unsigned int i=0; i<qrule.size(); i++){
            for (unsigned int j=0; j<qrule[i].first.size(); j++){
                std::cout << qrule[i].first[j] << " ";
            }
            std::cout << "(" << qrule[i].second << ")\n";
        }

    }

    void print(const Element &element){
        /*!
        * Print the element information to the terminal.
        */

        std::cout << "Element of type: " << element.name << "\n";
        std::cout << "\nglobal nodes:\n";
        print(element.nodes);
        std::cout << "\nglobal reference nodes:\n";
        print(element.reference_nodes);
        std::cout << "\nlocal nodes:\n";
        print(element.local_node_coordinates);
        std::cout << "\nquadrature rule:\n";
        print(element.qrule);
	std::cout << "\nbounding box:\n";
	print(element.bounding_box);
    }

    std::unique_ptr<Element> build_element_from_string(const std::string &eltype, const std::vector< unsigned int > &global_node_ids, 
                                                       const vecOfvec &nodes, const quadrature_rule &qrule){
	    /*
	     * Build an element from the element name, the nodes, and the quadrature rule
	     *
	     * TODO: Make the list of elements automatically populate.
	     *
	     * :param std::string &eltype: The name of the element
             * :param std::vector< unsigned int > &global_node_ids: The id numbers of the global nodes
	     * :param const vecOfvec &nodes: The element's nodes
	     * :param const quadrature_rule &qrule: The quadrature rule of the element
             */

	    if (std::strcmp(eltype.c_str(), "Hex8")==0){
                return std::unique_ptr<Element>(new Hex8(global_node_ids, nodes, qrule));
	    }
	    return NULL;
    }

    void determinant_3x3(const vecOfvec &A, double &d){
        /*
         * Compute the determinant of matrix A. Assumes A is a 3x3 matrix.
         * 
         * :param const vecOfvec &A: The matrix 
         * :param double d: The determinant of A
         */

         d  = 0;
         d += A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]);
         d -= A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]);
         d += A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
    }

}

