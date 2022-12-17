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

# delta = 0.7
h = 0.3
delta = h * rad
plot_bonds = True

Mesh =  genmesh(P_bdry=None, meshsize=None, pygmsh_geom=None, msh_file=msh_file, dimension=3, mesh_optimize=True)

fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
ax = Axes3D(fig)
Pos = Mesh.pos
ax.scatter(Pos[:,0], Pos[:,1], Pos[:,2])
plt.show()

