import numpy as np
import matplotlib.pyplot as plt
# to plot a collection of lines
# from matplotlib.collections import LineCollection

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

from pathos.multiprocessing import ProcessingPool as Pool
from itertools import combinations

from genmesh import genmesh

delta = 0.5

plot_bonds = 0

Mesh = genmesh(P_bdry=None, meshsize=None, msh_file ='3d_sphere_unit.msh' , do_plot = False, dotsize = 10, mesh_optimize=True )

##check if nodes are on the surface
# for i in range(len(Mesh.pos)):
    # p = Mesh.pos[i]
    # print(np.sqrt(np.sum(p**2)))

## check the sum of 'volume': surface area
# print('total surface area', np.sum(Mesh.vol, axis = 0))
# print('4 pi=', 4 * np.pi * 1)


total_nodes = len(Mesh.pos)
pairs = combinations(range(total_nodes),2)

def single_bond(mesh, ij_pair):
    i = ij_pair[0]
    j = ij_pair[1]
    p_i = mesh.pos[i]
    p_j = mesh.pos[j]
    xi = p_j - p_i
    d = np.sqrt(np.sum(np.square(xi))) #norm
    if (d <= delta):
                return [i,j, xi[0], xi[1], xi[2], d]

def oneparam_f_bond(ij):
    return single_bond(Mesh,ij)

## parallel attempt
a_pool = Pool()
all_bonds = a_pool.map(oneparam_f_bond, pairs) 
# remove all None
clean = np.array([i for i in all_bonds if i is not None])

Mesh.Conn = clean[:, 0:2]
Mesh.Conn_xi = clean[:, 2:5]
Mesh.Conn_xi_norm = clean[:, 5]


# save the mesh
# np.save('data/Mesh.npy', Mesh, allow_pickle=True)

if plot_bonds:
    linewidth =1
    print('Plotting bonds')

    fig = plt.figure()
    ax = Axes3D(fig)

    Pos = Mesh.pos

    V1 = Conn[:,0]
    V2 = Conn[:,1]

    P1 = Pos[V1]
    P2 = Pos[V2]
    ls =  [ [p1, p2] for p1, p2 in zip(P1,P2)] 
    lc = Line3DCollection(ls, linewidths=linewidth, colors='b')
    ax.add_collection(lc)

    # fix the axes
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, direc in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(direc))(ctr - r, ctr + r)
    plt.show()


# save
import pickle 
with open('data/Meshdump.pkl', 'wb') as config_dictionary_file:
  pickle.dump(Mesh, config_dictionary_file)
