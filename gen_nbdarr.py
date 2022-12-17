import numpy as np
import matplotlib.pyplot as plt
# to plot a collection of lines
# from matplotlib.collections import LineCollection

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

from pathos.multiprocessing import ProcessingPool as Pool
from itertools import combinations

from genmesh import genmesh

# msh_file='mesh/3d_sphere_unit.msh'
# msh_file='mesh/3d_sphere_unit_big.msh'
## specify ngores here that was used in mesh/<filename>.geo
ngores = 30
# rad = 150/np.pi
rad = 48
# msh_file='mesh/3d_sphere_forloop.msh'
msh_file='mesh/3d_sphere_forloop_big.msh'

neighbor_type = 'nearest_neighbor'
# neighbor_type = 'peridynamic'

# delta = 0.7
h = 0.3
delta = h * rad
plot_bonds = True

Mesh =  genmesh(P_bdry=None, meshsize=None, pygmsh_geom=None, msh_file=msh_file, dimension=3, mesh_optimize=True)

# store data in mesh
Mesh.delta = delta
Mesh.ngores = ngores
Mesh.rad = rad

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
                ## curved distance: distance along a great circle on which the points lie
                d = 2* mesh.rad *  np.arcsin(d/ (2*mesh.rad))
                return [i,j, d]

def oneparam_f_bond(ij):
    return single_bond(Mesh,ij)

# gores and nodes on tendon
Mesh.ngores = ngores
Mesh.gen_nodes_on_tendon()

print('Generating pairwise reference distance and connectivity.')

if neighbor_type == 'peridynamic':
    ## parallel attempt
    a_pool = Pool()
    all_bonds = a_pool.map(oneparam_f_bond, pairs) 
    # remove all None
    clean = np.array([i for i in all_bonds if i is not None])

    Mesh.Conn = clean[:, 0:2].astype(int)
    Mesh.Conn_xi_norm = clean[:, 2]

elif neighbor_type == 'nearest_neighbor':
    Mesh.Conn = Mesh.edges 
    Mesh.Conn_xi_norm = np.sqrt(np.sum((Mesh.pos[Mesh.Conn[:,1]] - Mesh.pos[Mesh.Conn[:,0]])**2, axis=1, keepdims=False))

print('Conn', Mesh.Conn)
print('Conn_xi_norm', Mesh.Conn_xi_norm)


# save the mesh
# np.save('data/Mesh.npy', Mesh, allow_pickle=True)

if plot_bonds:
    linewidth = 1
    print('Plotting bonds')

    fig = plt.figure()
    ax = Axes3D(fig)

    Pos = Mesh.pos


    ax.scatter(Pos[:,0], Pos[:,1], Pos[:,2])

    V1 = Mesh.Conn[:,0]
    V2 = Mesh.Conn[:,1]

    P1 = Pos[V1]
    P2 = Pos[V2]
    ls =  [ [p1, p2] for p1, p2 in zip(P1,P2)] 
    lc = Line3DCollection(ls, linewidths=linewidth, colors='b', alpha=0.5)
    ax.add_collection(lc)

    mx = np.amax(np.abs(Pos))
    XYZlim = [-mx, mx]
    ax.set_xlim3d(XYZlim)
    ax.set_ylim3d(XYZlim)
    ax.set_zlim3d(XYZlim)
    ax.set_box_aspect((1, 1, 1))

    plt.show()



# save
filename = 'data/Meshdump.pkl'
print('Saving mesh and connectivity to:', filename)
Mesh.save_state(filename)
