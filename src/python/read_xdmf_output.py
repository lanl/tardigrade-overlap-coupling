import os

import numpy as np
import matplotlib.pyplot as plt
import h5py
import xml.parsers.expat
import xml.etree

import micromorphic


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

                    if ("Time" == child.tag):
 
                        name = "time"

                        value = float(child.attrib['Value'])

                        if name in data.keys():
                            data[name].append(value)
                        else:
                            data.update({name:[value]})
    
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

def evaluate_model(data, parameters, model_name, parameters_to_fparams, nsdvs, maxinc=None, dim=3):
    """
    Evaluate the model given the parameters
    
    :param dict data: The data dictionary
    :param np.ndarray parameters: The array of parameters
    :param str model_name: The name of the model
    :param func parameters_to_fparams: A function that converts the parameters vector to the fparams vector required
        for the function
    :param int nsdvs: The number of solution dependant state variables
    :param int maxinc: The maximum increment to evaluate
    :param int dim: The spatial dimension of the problem
    """
    
    ninc = data['density_0'].shape[0]

    nel = data['density_0'].shape[1]

    nqp = len([key for key in data.keys() if 'density' in key])

    position, grad_u, phi, grad_phi = construct_degrees_of_freedom(data, dim = dim)

    PK2_sim   = dict([(qp,np.zeros((ninc,nel,dim * dim))) for qp in range(nqp)])
    SIGMA_sim = dict([(qp,np.zeros((ninc,nel,dim * dim))) for qp in range(nqp)])
    M_sim     = dict([(qp,np.zeros((ninc,nel,dim * dim * dim))) for qp in range(nqp)])
    SDVS_sim  = dict([(qp,np.zeros((ninc,nel,nsdvs))) for qp in range(nqp)])

    keys = ['errorCode', 'PK2', 'SIGMA', 'M', 'SDVS',\
            'DPK2Dgrad_u', 'DPK2Dphi', 'DPK2Dgrad_phi',\
            'DSIGMADgrad_u', 'DSIGMADphi', 'DSIGMADgrad_phi',\
            'DMDgrad_u', 'DMDphi', 'DMDgrad_phi',\
            'ADD_TERMS', 'ADD_JACOBIANS', 'output_message']
    
    time = data['time'][1:]
    
    if maxinc is None:
        maxinc = ninc

    tp = 0
        
    for e in range(nel):
    
        for qp in range(nqp):
        
            for i in range(maxinc):

                # Map the parameters vector to the function parameters
                fparams = parameters_to_fparams(parameters)
            
                sp = 0
                ds = 1.
                
                if (i == 0):
                    previous_SDVS_s = np.zeros(nsdvs)
                else:
                    previous_SDVS_s = np.copy(SDVS_sim[qp][i-1,e,:])
            
                while (sp < 1.0):
                
                    s = sp + ds
                
                    time_1     = time[i]
                    grad_u_1   = grad_u[qp][i, e, :, :]
                    phi_1      = phi[qp][i, e, :, :]
                    grad_phi_1 = grad_phi[qp][i, e, :, :, :]
                    
                    if (i == 0):
                        time_0     = 0
                        grad_u_0   = np.zeros((3,3))
                        phi_0      = np.zeros((3,3))
                        grad_phi_0 = np.zeros((3,3,3))
                        
                    else:
                        time_0     = time[i-1]
                        grad_u_0   = grad_u[qp][i-1, e, :, :]
                        phi_0      = phi[qp][i-1, e, :, :]
                        grad_phi_0 = grad_phi[qp][i-1, e, :, :]
                        
                    t                = (time_1 - time_0) * s + time_0
                    current_grad_u   = (grad_u_1 - grad_u_0) * s + grad_u_0
                    current_phi      = (phi_1 - phi_0) * s + phi_0
                    current_grad_phi = (grad_phi_1 - grad_phi_0) * s + grad_phi_0
                    
                    tp                = (time_1 - time_0) * sp + time_0
                    previous_grad_u   = (grad_u_1 - grad_u_0) * sp + grad_u_0
                    previous_phi      = (phi_1 - phi_0) * sp + phi_0
                    previous_grad_phi = (grad_phi_1 - grad_phi_0) * sp + grad_phi_0

                    current_phi = current_phi.flatten()
                    previous_phi = previous_phi.flatten()
                    
                    current_grad_phi = current_grad_phi.reshape((dim * dim, dim))
                    previous_grad_phi = previous_grad_phi.reshape((dim * dim, dim))

                    #TODO: add dof and add grad dof not currently used
                    current_ADD_DOF = np.zeros((1))
                    current_ADD_grad_DOF = np.zeros((1,3))

                    previous_ADD_DOF = np.zeros((1))
                    previous_ADD_grad_DOF = np.zeros((1,3))
                    
                    # Evaluate the model
                    values = micromorphic.evaluate_model(model_name, np.array([t, t - tp]), fparams,\
                                                         current_grad_u, current_phi, current_grad_phi,
                                                         previous_grad_u, previous_phi, previous_grad_phi,
                                                         previous_SDVS_s,\
                                                         current_ADD_DOF, current_ADD_grad_DOF,
                                                         previous_ADD_DOF, previous_ADD_grad_DOF)
                    
                    results = dict(zip(keys, values))
                    
                    if (results['errorCode'] == 1):
                        print("error")
                        ds = 0.5 * ds
                        nsubiter += 1
                        
                        if (nsubiter > maxsubiter):
                            break
                            
                    elif (results['errorCode'] == 2):
                        errormessage = f"evaluate_model return error code {results['errorCode']}\n\n"
                        errormessage += results['output_message'].decode("utf-8")
                        raise IOError(errormessage)
                        
                    else:
                        sp += ds
                        nsubiter = 0
                        
                        if np.isclose(sp, 1):
                            ds = 1
                        else:
                            ds = 1 - sp
                            
                        previous_SDVS_s = np.copy(results['SDVS'])
                        
                if (results['errorCode'] != 0):
                    errormessage = f"evaluate_model returned error code {results['errorCode']}\n\n"
                    errormessage += results['output_message'].decode('utf-8')
                    print(parameters, 'fail')
                    
                    return np.nan
                
                PK2_sim[qp][i,e,:]   = results['PK2']
                SIGMA_sim[qp][i,e,:] = results['SIGMA']
                M_sim[qp][i,e,:]     = results['M']
                SDVS                 = results['SDVS']
                SDVS_sim[qp][i,e,:]  = results['SDVS']
                
    return PK2_sim, SIGMA_sim, M_sim, SDVS_sim
