from genmesh import genmesh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# msh_file = 'mesh/3d_sphere_unit.msh' 
msh_file = 'mesh/3d_sphere_forloop_big.msh'
mesh = genmesh(P_bdry=None, meshsize=None, msh_file =msh_file, dimension=3, mesh_optimize=True)

# mesh.plot(plot_node_text=False, highlight_bdry_nodes=False, dotsize=20)

fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
ax = Axes3D(fig)
Pos = mesh.pos
plt.scatter(Pos[:,0], Pos[:,1], Pos[:,2])
plt.show()
