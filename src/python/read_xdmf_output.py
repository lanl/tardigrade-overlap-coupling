import os

import numpy as np
import matplotlib.pyplot as plt
import h5py
import xml.parsers.expat
import xml.etree

def get_attribute_data(xml_node, path='.'):
    
    if ("Attribute" != xml_node.tag):
        raise ValueError("XML node is not of type 'Attribute'")
        
    value = None
        
    for child in xml_node:
                            
        if "DataItem" == child.tag:
                                
            value = read_data(child, path=path)
                
            break
            
    return value

def get_set(xml_node, path='.'):
    
    if ("Set" != xml_node.tag):
        raise ValueError("XML node is not of type 'Set'")
        
    for child in xml_node:
        
        if "DataItem" == child.tag:
            
            value = read_data(child, path=path)
            
            break     
        
    return value

def get_geometry(xml_node, path='.'):
    
    if ("Geometry" != xml_node.tag):
        raise ValueError( "XML node is not of type 'Geometry'")
        
    for child in xml_node:
        
        if "DataItem" == child.tag:
            
            value = read_data(child, path=path)
            
            break
            
    return value

def get_topology(xml_node, path='.'):
    
    if ("Topology" != xml_node.tag):
        raise ValueError( "XML node is not of type 'Topology'")
        
    for child in xml_node:
        
        if "DataItem" == child.tag:
            
            value = read_data(child, path=path)
            
            break
            
    return value

def read_data(xml_node, path='.'):
    
    if ("DataItem" != xml_node.tag):
        raise ValueError("XML node is not of type 'DataItem'")
        
    shape = [int(v) for v in xml_node.attrib["Dimensions"].split()]
    
    if xml_node.attrib["DataType"] == "Float":
        dtype = float
    else:
        dtype = int
        
    if xml_node.attrib["Format"] == "XML":
        value = np.hstack([np.array(line.split()).astype(dtype) for line in xml_node.text.split("\n")]).reshape(shape)
        
    else:
        
        filename, keyname = xml_node.text.split(':')
        
        with h5py.File(os.path.join(path, filename), 'r') as h5file:
            
            value = np.array(h5file[keyname]).astype(dtype).reshape(shape)
        
    return value

def parse_xdmf_output(filename):
    """
    Parse an XDMF filter output file into a python dictionary object

    :param str filename: The filename to be processed
    """

    tree = xml.etree.ElementTree.parse(filename)
    root = tree.getroot()
    
    data = {}
    geometry = {}
    topology = {}
    
    for domain in root:
    
        for collection in domain:
            if collection.tag != "Grid":
                continue
    
            for grid in collection:
                if grid.tag != "Grid":
                    continue
    
                for child in grid:
    
                    if ("Attribute" == child.tag):
    
                        name = child.attrib["Name"]
    
                        value = get_attribute_data(child)
    
                        if name in data.keys():
                            data[name].append(value)
                        else:
                            data.update({name:[value]})
    
                    if ("Geometry" == child.tag):
    
                        name = "geometry"
    
                        value = get_geometry(child)
    
                        if name in geometry.keys():
                            geometry[name].append(value)
                        else:
                            geometry.update({name:[value]})
    
                    if ("Topology" == child.tag):
    
                        name = "topology"
    
                        value = get_topology(child)
    
                        if name in topology.keys():
                            topology[name].append(value)
                        else:
                            topology.update({name:[value]})
    
    return dict([(name, np.vstack(data[name])) for name in data.keys()]), geometry, topology

def construct_degrees_of_freedom(data, dim = 3):
    """
    Construct the degrees of freedom from the data

    :param dict data: The output dictionary produced by parse_xdmf_output
    :param int dim: The spatial dimension of the filter
    """

    ninc = data['density_0'].shape[0]

    nel = data['density_0'].shape[1]

    nqp = len([key for key in data.keys() if 'density' in key])

    position = {}

    grad_u = {}

    phi = {}

    grad_phi = {}

    for qp in range(nqp):

        # Collect the position of the outputs
        values = np.zeros((ninc, dim, nel))

        for i in range(dim):

            root_string = f'position_{i+1}_{qp}'

            values[:,i,:] = data[root_string]

        position.update({qp:np.copy(values)})

        # Collect the gradient of the macro deformation
        values = np.zeros((ninc, nel, dim, dim))

        for i in range(dim):
            for j in range(dim):
                root_string = f'DOF_{i+1},{j+1}_{qp}'

                values[:,:,i,j] = data[root_string]

        grad_u.update({qp:np.copy(values)})

        # Collect the micro deformation
        values = np.zeros((ninc, nel, dim, dim))

        for i in range(dim):
            for j in range(dim):
                I = dim * i + j

                root_string = f'DOF_{I+4}_{qp}'

                values[:,:,i,j] = data[root_string]

        phi.update({qp:np.copy(values)})

        # Collect the gradient of the micro deformation
        values = np.zeros((ninc, nel, dim, dim, dim))

        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    I = dim * i + j

                    root_string = f'DOF_{I+4},{k+1}_{qp}'

                    values[:,:,i,j,k] = data[root_string]

        grad_phi.update({qp:np.copy(values)})

    return position, grad_u, phi, grad_phi

def collect_stresses(data, dim=3):
    """
    Collect all of the stresses

    :param dict data: The output dictionary produced by parse_xdmf_output
    :param int dim: The spatial dimension of the filter
    """

    ninc = data['density_0'].shape[0]

    nel = data['density_0'].shape[1]

    nqp = len([key for key in data.keys() if 'density' in key])

    cauchy_stress = {}

    symmetric_micro_stress = {}

    higher_order_stress = {}

    for qp in range(nqp):

        # Collect the Cauchy stress
        values = np.zeros((ninc, nel, dim, dim))

        for i in range(dim):
            for j in range(dim):
                root_string = f'cauchy_stress_{i+1}{j+1}_{qp}'

                values[:,:,i,j] = data[root_string]

        cauchy_stress.update({qp:np.copy(values)})

        # Collect the symmetric micro stress
        values = np.zeros((ninc, nel, dim, dim))

        for i in range(dim):
            for j in range(dim):
                root_string = f'symmetric_micro_stress_{i+1}{j+1}_{qp}'

                values[:,:,i,j] = data[root_string]

        symmetric_micro_stress.update({qp:np.copy(values)})

        values = np.zeros((ninc, nel, dim, dim, dim))

        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    root_string = f'higher_order_stress_{i+1}{j+1}{k+1}_{qp}'

                    values[:,:,i,j,k] = data[root_string]

        higher_order_stress.update({qp:np.copy(values)})

    return cauchy_stress, symmetric_micro_stress, higher_order_stress

def get_reference_configuration_stresses(data, dim=3):
    """
    Collect all of the stresses in the reference configuration

    :param dict data: The output dictionary produced by parse_xdmf_output
    :param int dim: The spatial dimension of the filter
    """

    position, grad_u, phi, grad_phi = construct_degrees_of_freedom(data, dim=dim)

    cauchy_stress, symmetric_micro_stress, higher_order_stress = collect_stresses(data, dim=dim)

    PK2 = {}

    SIGMA = {}

    M = {}

    ninc = data['density_0'].shape[0]

    nel = data['density_0'].shape[1]

    nqp = len([key for key in data.keys() if 'density' in key])

    for qp in range(nqp):

        PK2.update({qp:np.zeros((ninc, nel, dim, dim))})

        SIGMA.update({qp:np.zeros((ninc, nel, dim, dim))})

        M.update({qp:np.zeros((ninc, nel, dim, dim, dim))})

        for inc in range(ninc):

            for el in range(nel):

                # Construct the deformation gradient

                F = grad_u[qp][inc,el,:,:] + np.eye(3)

                Finv = np.linalg.inv(F)

                J = np.linalg.det(F)

                # Construct the micro-displacement

                chi = phi[qp][inc,el,:,:] + np.eye(3)

                chiinv = np.linalg.inv(chi)

                PK2[qp][inc,el,:,:] = J * np.einsum("Kk,kl,Ll->KL", Finv, cauchy_stress[qp][inc,el,:,:], Finv)

                SIGMA[qp][inc,el,:,:] = J * np.einsum("Kk,kl,Ll->KL", Finv, symmetric_micro_stress[qp][inc,el,:,:], Finv)

                M[qp][inc,el,:,:,:] = J * np.einsum("Kk,Ll,Mm,klm->KLM", Finv, Finv, chiinv, higher_order_stress[qp][inc,el,:,:,:])

    return PK2, SIGMA, M

def compute_2order_J2yield(F, A):
    """
    Compute the second order J2 flow yield condition

    :param np.ndarray F: The deformation gradient
    :param np.ndarray A: The second order tensor stress metric
    """

    C = np.dot(F.T, F)

    p = np.einsum('AB,AB', C, A) / 3

    devA = A - p * np.linalg.inv(C)

    return np.sqrt(np.einsum('AB,AB', devA, devA))

def compute_3order_J2yield(F, A):
    """
    Compute the third order J2 flow yield condition

    :param np.ndarray F: The deformation gradient
    :param np.ndarray A: The third order tensor stress metric
    """

    C = np.dot(F.T, F)

    p = np.einsum('AB,ABK', C, A) / 3

    devA = A - np.einsum('IJ,K', np.linalg.inv(C), p)

    G = np.zeros(3)

    for i in range(3):
        for j in range(3):
            for k in range(3):
                G[k] += devA[i,j,k] * devA[i,j,k]

    return np.sqrt(G)

def compute_strains(data, dim=3):
    """
    Compute the strain metrics

    :param dict data: The output dictionary produced by parse_xdmf_output
    :param int dim: The spatial dimension of the filter
    """

    position, grad_u, phi, grad_phi = construct_degrees_of_freedom(data, dim=dim)

    GreenLagrangeStrain = {}

    MicroGreenLagrangeStrain = {}

    Gamma = {}

    ninc = data['density_0'].shape[0]

    nel = data['density_0'].shape[1]

    nqp = len([key for key in data.keys() if 'density' in key])

    for qp in range(nqp):

        GreenLagrangeStrain.update({qp:np.zeros((ninc,nel,dim,dim))})

        MicroGreenLagrangeStrain.update({qp:np.zeros((ninc,nel,dim,dim))})

        Gamma.update({qp:np.zeros((ninc,nel,dim,dim,dim))})

        for inc in range(ninc):

            for e in range(nel):

                F = grad_u[qp][inc,e,:,:] + np.eye(3)

                chi = phi[qp][inc,e,:,:] + np.eye(3)

                grad_chi = grad_phi[qp][inc,e,:,:,:]

                GreenLagrangeStrain[qp][inc,e,:,:] = 0.5 * (np.dot(F.T, F) - np.eye(3))

                MicroGreenLagrangeStrain[qp][inc,e,:,:] = np.dot(F.T, chi) - np.eye(3)

                Gamma[qp][inc,e,:,:,:] = np.einsum('iI,iJK->IJK', F, grad_chi)

    return GreenLagrangeStrain, MicroGreenLagrangeStrain, Gamma
