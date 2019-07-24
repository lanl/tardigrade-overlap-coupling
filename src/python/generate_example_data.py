import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from tqdm import tqdm

import scipy.stats as stats
np.random.seed(123)
npts = 1000

##Random points
#x = np.array([2*(np.random.rand(3)-0.1) for _ in range(npts)])

#Structured points
x_ = np.linspace(0, 1, 10);
y_ = np.linspace(0, 1, 10);
z_ = np.linspace(0, 1, 10);
X, Y, Z = np.meshgrid(x_, y_, z_);

x = np.vstack([X.flatten(), Y.flatten(), Z.flatten()]).T

def density(x):
    return 2700.

def rotation_matrix(random=True, angles=None):
    if (random):
        thx = stats.norm(loc=np.pi/3, scale=0.05).rvs()
        thy = stats.norm(loc=-np.pi/5, scale=0.1).rvs()
        thz = stats.norm(loc=0, scale=0.2).rvs()
    else:
        thx, thy, thz = angles
    
    Rx = np.array([[1,           0,            0],\
                   [0, np.cos(thx), -np.sin(thx)],\
                   [0, np.sin(thx),  np.cos(thx)]])
    
    Ry = np.array([[ np.cos(thy), 0, np.sin(thy)],\
                   [          0, 1,           0],\
                   [-np.sin(thy), 0, np.cos(thy)]])
    
    Rz = np.array([[np.cos(thz), -np.sin(thz), 0],\
                   [np.sin(thz),  np.cos(thz), 0],\
                   [          0,            0, 1]])
    
    return Rz.dot(Ry.dot(Rx))

def transformation(x, t):
    C = t*np.array([[ 1.00, 0.10,-0.10],\
                    [ 0.10,-0.25, 0.07],\
                    [-0.10, 0.07,-0.25]]) + np.eye(3)
    b = t*np.array([1.3, -.3, 2]).reshape((-1, 1)) 
    Q = rotation_matrix(random=False, angles=t*np.array([np.pi/3, np.pi/2, np.pi/8]))
    D = np.matmul(Q, C)

    xp = (D.dot(x.reshape((-1, 1))) + b*t).flatten()
    return xp.flatten()

def stress(x, t):
    principal = np.diag([2.0*x[0], 0.0, 0.0])*t
    Q = rotation_matrix()
    return Q.dot(principal).dot(Q.T)

def voigt(A):
    order = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1), (2, 1), (2, 0), (1, 0)]
    return np.array([A[i,j] for i,j in order])

#for t in np.linspace(0, 1, 10):
#    transformation(np.zeros((3,)), t)
#raise()
    
import pickle
fn = "example_data.txt"
f = open(fn, 'w')

header  = "This file contains fake data that is not intended to be representative\n"
header += "of any true material. It exists merely to demonstrate the filter's\n"
header += "capabilities.\n"
header += "BEGIN DATA\n"
f.write(header)

time = np.linspace(0, 1, 11)
for t in time:
    f.write("t = {0}\n".format(t))
    for i, point in enumerate(tqdm(x)):
        tpoint = transformation(point, t)
        line  = "MP, {0}, {1}, {2}, {3}".format(i, tpoint[0], tpoint[1], tpoint[2])
        line += ", {0}".format(density(point))
        line += "".join([", {0}".format(s) for s in voigt(stress(point, t))])
        line += "\n"
        f.write(line)
f.close()

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(*zip(*x), label='reference')
for t in np.linspace(0, 1, 5):
    xp = np.array([transformation(xi, t) for xi in x])
    ax.scatter(*zip(*xp), label='t = {0}'.format(t))
ax.legend()
plt.savefig('transformation.pdf')
plt.show()
