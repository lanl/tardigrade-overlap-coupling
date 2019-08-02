import numpy as np

def read_values(datafile):
    """
    Read in values from a string until the keyword character is found.

    :param file datafile: The datafile
    """
    
    last_pos = datafile.tell()
    line = datafile.readline().strip()
    while (len(line)==0):
        line = datafile.readline()
            
    sline = [v.strip() for v in line.split(',') if len(v.strip())>0]
    kwdata = np.array([float(s) for s in sline])
        
    last_pos = datafile.tell()
    line = datafile.readline()
    while not ("*" in line):
        if (len(line)==0):
            last_pos = datafile.tell()
            line = datafile.readline().strip()
            continue
        sline = [v.strip() for v in line.split(',') if len(v.strip())>0]
        kwdata = np.vstack([kwdata, [float(s) for s in sline]])
        x = datafile.tell()
        line = datafile.readline().strip()
    datafile.seek(x)
    return kwdata

class MicromorphicFilterData(object):
    """
    A class which stores micromorphic filter data
    """
    def __init__(self):
        self.nodes = [] #The filter's nodes
        self.dof_values = [] #The filter's degree of freedom values
        self.gauss_point_info = GaussPointInformation() #The filter's values at the gauss points
        
    def __repr__(self):
        
        out_str  = "MicromorphicFilterData:\n"
        out_str += " nodes:\n"
        for n in self.nodes:
            out_str += "  "
            for ni in n:
                out_str += " {0:+2.4f}".format(ni)
            out_str += "\n"
            
        out_str += " dof values:\n"
        for d in self.dof_values:
            out_str += "  "
            for di in d:
                out_str += " {0:+2.4f}".format(di)
            out_str += "\n"
            
        out_str += self.gauss_point_info.__repr__(offset=1)
        


        return out_str
        
class GaussPointInformation(object):
    """
    A class which stores gauss point information
    """
    def __init__(self):
        self.volume = [] #The volume of the gauss domain
        self.density = [] #The density at the gauss point
        self.local_mass_center = [] #The local mass center associated with the gauss point
        self.global_mass_center = [] #The global mass center associated with the gauss point
        self.surface_area = [] #The surface area of the faces associated with the gauss point
        
    def __repr__(self, offset=0):
        out_str = " "*offset + "Gauss Point Information\n"
        
        out_str += " "*offset + " Volumes:  "
        
        index = 0
        for v in self.volume:
            if (index > 8):
                out_str += "\n" + " "*offset + "          "
            
            out_str += " {0:2.4f}".format(v)
            index += 1
        out_str += "\n"
        
        out_str += " "*offset + " Densities:"
        index = 0
        for d in self.density:
            if (index > 8):
                out_str += "\n" + " "*offset + "            "
            
            out_str += " {0:2.4f}".format(d)
            index += 1
        out_str += "\n"
        
        out_str += " "*offset + " Local Mass Centers:\n"
        for lmc in self.local_mass_center:
            out_str += " "*offset + "  "
            for lmci in lmc:
                out_str += " {0:+1.4f}".format(lmci)
            out_str += "\n"
            
        out_str += " "*offset + " Global Mass Centers:\n"
        for gmc in self.global_mass_center:
            out_str += " "*offset + "  "
            for gmci in gmc:
                out_str += " {0:+1.4f}".format(gmci)
            out_str += "\n"

        out_str += " "*offset + " Surface Areas:\n";
        for sa in self.surface_area:
            for key in sa.keys():
                out_str += " "*offset + "  "
                out_str += "{0:3d} : ".format(key)
                out_str += "{0:1.4f}\n".format(sa[key])
            out_str += "\n"
        
        return out_str
        
def read_output_data(output_fn):
    """
    Read an output file from the micromorphic filter

    :param str output_fn: The output filename
    """
    of = open(output_fn, 'r')
    
    filter_results = {}
    
    fnum = -1
        
    line = of.readline()
        
    while line != '':
        
        if (len(line)==0):
            continue
        
        line = line.strip()
        sline = line.split(',')
        
        if "*INPUT_FILE" in sline[0]:
            print("Reading data from: ", sline[1])
        
        elif "*FILTER_CONFIGURATION" in sline[0]:
            print("Filter defined in: ", sline[1])
        
        elif "*ELEMENT" in sline[0]:
            do_nothing=0
            
        elif "*NODES" in sline[0]:
            mfdata.nodes = read_values(of)
        
        elif "*TIMESTEP" in sline[0]:
            
            if fnum != -1:
                mfdata.gauss_point_info = gpinfo
                filter_results[time].update({fnum:mfdata})
            
            sline = line.split(',')
            time = float(sline[1])
            filter_results.update({time:{}})
                
        elif "*MICROMORPHIC FILTER" in line:
            if fnum != -1:
                mfdata.gauss_point_info = gpinfo
                filter_results[time].update({fnum:mfdata})
            fnum = int(sline[1])
            mfdata = MicromorphicFilterData()
            
        elif "*DOF VALUES" in sline[0]:
            mfdata.dof_values = read_values(of)
            
        elif "*GAUSS POINT INFORMATION" in sline[0]:
            gpinfo = GaussPointInformation()
            
        elif "*VOLUME" in sline[0]:
            volume = float(sline[1])
            gpinfo.volume.append(volume)
            
        elif "*DENSITY" in sline[0]:
            density = float(sline[1])
            gpinfo.density.append(density)
            
        elif "*LOCAL MASS CENTER" in sline[0]:
            lmc = [float(sl.strip()) for sl in sline[1:] if len(sl.strip())>0]
            gpinfo.local_mass_center.append(lmc)
            
        elif "*GLOBAL MASS CENTER" in sline[0]:
            gmc = [float(sl.strip()) for sl in sline[1:] if len(sl.strip())>0]
            gpinfo.global_mass_center.append(gmc)

        elif "*SURFACE AREAS" in sline[0]:
            surface_areas = read_values(of)
            faces = (surface_areas[:, 0]+0.5).flatten().astype(int)
            values = surface_areas[:, 1:].flatten()
            gpinfo.surface_area.append(dict([(f, v) for f, v in zip(faces, values)]))
            
        elif "*" in sline[0]:
            print("Warning: unknown keyword detected")
            print("  ", line)
            
        else:
            print(sline)
#            data_array.append([float(v) for v in sline])
            
        line = of.readline()
        
    if fnum != -1:
        mfdata.gauss_point_info = gpinfo
        filter_results[time].update({fnum:mfdata})

    return filter_results

def read_filter_configuration(config_fn):
    """
    Read in a filter configuration file

    :param str config_fn: The filename of the configuration file.
    """
    fc = open(config_fn, 'r')
    
    begin_data = False
    
    nodes = {}
    conn = {}
    
    for line in fc.readlines():
        
        if begin_data:
        
            data = [s.replace(',', '') for s in line.split()]
            
            if(len(data)>0):
            
                if (data[0]=='N'):
                    nid = int(data[1])
                    pos = np.array([float(data[i]) for i in range(2, len(data))])
                    nodes.update({nid:pos})
                    
                elif (data[0]=='E'):
                    eid = int(data[2])
                    nids = [int(data[i]) for i in range(3, len(data))]
                    conn.update({eid:nids})
            
        if (line.strip() == "BEGIN DATA"):
            begin_data = True
    
    fc.close()
    return nodes, conn

def read_dns_data(dns_fn):
    """
    Read data in from a DNS file

    :param str dns_fn: The filename of the DNS
    """

    fed = open(dns_fn, 'r')
    
    begin_data = False
    dns_data = {}
    for line in fed.readlines():
        
        if begin_data:
            
            if "t = " in line:
                tc = float(line[3:])
                dns_data.update({tc:np.empty((0, 3))})
            else:
                data = [s.replace(',', '') for s in line.split()]
                pos = np.array([float(data[i]) for i in range(2, 5)])
                dns_data[tc] = np.vstack([dns_data[tc], pos])
        
        if (line.strip() == "BEGIN DATA"):
            begin_data = True
        
    fed.close()
    return dns_data

def plot_hex(eid, nodes, conn, ax, u=None):
    """
    Plot a hexahedral element

    :param int eid: The element id
    :param dict nodes: The global nodes dictionary
    :param dict conn: The connectivity dictionary
    :param matplotlib.axes._subplots.Axes3DSubplot ax: The axis to plot to
    :param dict u: The displacement of the nodes.
    """

    if (u is None):
        u = dict([(key, np.zeros(nodes[key].shape)) for key in nodes.keys()])
    
    #Plot bottom square
    x = np.array([nodes[n]+u[n] for n in conn[eid][:4]] + [nodes[conn[eid][0]]+u[conn[eid][0]]])
    #print(x)
    ax.plot(*zip(*x), color='k', marker='o')
    #Plot top square
    x = np.array([nodes[n]+u[n] for n in conn[eid][4:]] + [nodes[conn[eid][4]]+u[conn[eid][4]]])
    #print(x)
    ax.plot(*zip(*x), color='k', marker='o')
    #Plot sides
    for i in range(4):
        x = [nodes[conn[eid][i]]+u[conn[eid][i]], nodes[conn[eid][i+4]]+u[conn[eid][i+4]]]
        #print(x)
        ax.plot(*zip(*x), color='k', marker='o')

def plot_hex_filter(eid, t, filter_data, ax):
    """
    Plot a hexahedral element as defined by the filter's nodes

    :param int eid: The element (filter) id
    :param float t: The time
    :param dict filter_data: The output from read_filter_data
    :param matplotlib.axes._subplots.Axes3DSubplot ax: The axis to plot to
    """

    f = filter_data[t][eid]

    #Plot bottom square
    x = np.vstack([f.nodes[:4], f.nodes[0]])
    ax.plot(*zip(*x), color='k', marker='o')

    #Plot the top square
    x = np.vstack([f.nodes[4:], f.nodes[4]])
    ax.plot(*zip(*x), color='k', marker='o')

    #Plot the sides
    for i in range(4):
        x = [f.nodes[i], f.nodes[i+4]]
        ax.plot(*zip(*x), color='k', marker='o')

def plot_filter_cgs(eid, t, filter_data, ax):
    """
    Plot the centers of gravity of the gauss domains for the given filter.

    :param int eid: The element (filter) id
    :param float t: The time
    :param dict filter_data: The output from read_filter_data
    :param matplotlib.axes._subplots.Axes3DSubplot ax: The axis to plot to
    """

    ax.scatter(*zip(*filter_data[t][eid].gauss_point_info.global_mass_center), color='r', marker='x')

def to_voigt(A):
    indices = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1), (2, 1), (2, 0), (1, 0)]
    return [A[ind] for ind in indices]

def from_voigt(Av):
    indices = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1), (2, 1), (2, 0), (1, 0)]
    
    A = np.zeros((3, 3))
    
    for i, ind in enumerate(indices):
        A[ind[0], ind[1]] = Av[i]
        
    return A

def matsqrt(A, mode=1, maxiter=20, maxls=5, tola=1e-9, tolr=1e-9):
    """
    Compute X given
        mode = 1: A = X*X
        mode = 2: A = X.T*X
        mode = 3: A = X * X.T
    """
    
    def jacobian(X, mode, indices):
        
        dXdX = np.zeros((X.shape[0], X.shape[1], X.shape[0], X.shape[1]))
        dXdXm = np.zeros((indices.shape[0], indices.shape[0]))
        
        for I, indx_1 in enumerate(indices):
            i, j = indx_1
            for J, indx_2 in enumerate(indices):
                k, l = indx_2
                if ((i==k) and (j==l)):
                    dXdX[i, j, k, l] = 1
                    
                    if ((mode==2) or (mode==3)):
                        dXdX[j, i, k, l] = dXdX[i, j, l, k] = dXdX[j, i, l, k] = 1
                        #dXdX[j, i, l, k] = 1
                    
                dXdXm[I, J] = dXdX[i, j, k, l]
#        print(dXdXm)
        
                
        
        if mode==1:
            Jac = -(np.einsum('ijlm, jk->iklm', dXdX, X) + np.einsum('ij, jklm->iklm', X, dXdX))
            
        elif (mode==2):
            #errstr = "not working"
            #raise ValueError(errstr)
            Jac = -(np.einsum('jilm, jk->iklm', dXdX, X) + np.einsum('ji, jklm->iklm', X, dXdX))
            
        elif (mode==3):
            Jac = -(np.einsum('ijlm, kj->iklm', dXdX, X) + np.einsum('ij, kjlm->iklm', X, dXdX))
            
        else:
            errstr = "mode not recognized"
            raise ValueError(errstr)
            
        Jm = np.zeros((indices.shape[0], indices.shape[0]))
        for I, indx_1 in enumerate(indices):
            i, k = indx_1
            #print('I, i, k: ',I, i, k)
            
            for J, indx_2 in enumerate(indices):
                l, m = indx_2
                #print("J, l, m: ", J, l, m)
                
                Jm[I, J] = Jac[i, k, l, m]
                
        return Jm
                    
    def fd_jacobian(Af, Xf, mode, shape, indices):
        
        eps = 1e-6
        Jest = np.zeros((Af.size, Af.size))
        R = residual(Af, Xf, mode, shape, indices, include_jacobian=False)
        for i in range(Af.size):
            delta = np.ones(Xf.size)
            delta[i] += eps
        
            Rpi = residual(Af, Xf*delta, mode, shape, indices, include_jacobian=False)
            Jest[:, i] = (Rpi-R)/np.linalg.norm(Xf - Xf*delta)
            
        return Jest
    
    def residual(Af, Xf, mode, shape, indices, include_jacobian=True):
        
        X = np.zeros(shape)
        
        for i, indx in enumerate(indices):
            X[indx[0], indx[1]] = Xf[i]
        
        
        if mode==1:
            As = X.dot(X)
            
        elif (mode==2):
            d = np.diag(X)
            X = X + X.T - np.diag(d)
            As = (X.T).dot(X)
            
        elif (mode==3):
            d = np.diag(X)
            X = X + X.T - np.diag(d)
            As = X.dot(X.T)
            
        else:
            errstr = "mode not recognized"
            raise ValueError(errstr)
            
        As = np.array([As[v[0], v[1]] for v in indices])
#        print('As: ',As)
            
        if include_jacobian:
            #if mode==1:
                
            return Af - As, jacobian(X, mode, indices)
            
            #else:
            #    return Af - As, fd_jacobian(Af, Xf, mode, shape, indices)
        else:
            return Af - As
        
        
    if ((mode==2) or (mode==3)):
        Askew = A - A.T
        if (not (np.allclose(Askew, 0))):
            errstr = "for mode 2 and 3 A must be symmetric"
            raise ValueError(errstr)
            
        indices = np.vstack(np.triu_indices(A.shape[0], k=0, m=A.shape[1])).T
    else:
        indices = [None for _ in range(A.shape[0]*A.shape[1])]
        indx = 0
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                indices[indx] = (i, j)
                indx += 1
        indices = np.array(indices)
        
#    print('indices: ', indices)
        
    #Flatten A into a vector
    Af = np.array([A[v[0], v[1]] for v in indices])
    shape = A.shape
    
    #Initialize the Newton-Raphson solution
    niter = 0
    Xf = np.eye(Af.shape[0])#copy(Af)
    Xf = np.array([Xf[v[0], v[1]] for v in indices])
    
    R, J = residual(Af, Xf, mode, shape, indices, include_jacobian=True)
    
    eps = 1e-7
    Jest = np.zeros((R.size, R.size))
    for i in range(R.size):
        delta = np.ones(Xf.size)
        delta[i] += eps
        
        Rpi = residual(Af, Xf*delta, mode, shape, indices, include_jacobian=False)
        Jest[:, i] = (Rpi-R)/np.linalg.norm(Xf - Xf*delta)
    
    Rnorm = Rp = np.linalg.norm(R)
    tol = Rnorm*tolr + tola
    
    #Begin the loop
    while ((Rnorm>tol) and (niter<maxiter)):
        dXf = -np.linalg.solve(J, R)
        
        Xf += dXf
        
        R, J = residual(Af, Xf, mode, shape, indices, include_jacobian=True)
        Rnorm = np.linalg.norm(R)
        
        #line search
        lamb = 1
        nls = 0
        while ((Rnorm>Rp) and (nls<maxls)):
            Xf -= dXf
            lamb /= 2
            dXf *= lamb
            Xf += dXf
            R, J = residual(Af, Xf, mode, shape, indices, include_jacobian=True)
            Rnorm = np.linalg.norm(R)
            nls += 1
            
        if ((nls>=maxls) and (Rnorm>Rp)):
            errstr = "failure in line search"
            raise ValueError(errstr)
        else:
            Rp = Rnorm
            
        niter += 1
        
    if (Rnorm>tol):
        errstr = "failure in Newton-Raphson"
        raise ValueError(errstr)
    
    X = np.zeros(shape)
    for i, index in enumerate(indices):
        X[index[0], index[1]] = Xf[i]
        
    if ((mode==2) or (mode==3)):
        X += X.T - np.diag(np.diag(X))
    
    return X
