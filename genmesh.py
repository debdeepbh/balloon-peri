import pygmsh
import numpy as np
import matplotlib.pyplot as plt
import pickle

# import optimesh

# 3d plot
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
# from mpl_toolkits import mplot3d

import meshio

def get_edges(T):
    """returns all the edges connecting various nodes
    : T: triangles 
    """
    E1 = np.c_[T[:,0], T[:,1]]
    E2 = np.c_[T[:,1], T[:,2]]
    E3 = np.c_[T[:,2], T[:,0]]

    E = np.r_[E1, E2, E3]
    # sort
    E.sort(axis = 1)
    # select the unique rows
    E = np.unique(E, axis = 0)
    return E


class ReturnValue(object):
    """Returns the position of nodes and nodal volume"""
    def __init__(self, pos, vol, bdry_nodes, T, edges):

        self.zero_l_nodes = None

        self.pos = pos
        self.vol = vol
        self.T = T
        self.edges = edges

        self.bdry_nodes = bdry_nodes

        self.Conn = None
        # self.Conn_xi = None
        self.Conn_xi_norm = None

        self.NArr = None
        # self.xi = None
        self.xi_norm = None

        # for tendons
        self.rad = None
        self.ngores = None
        self.nodes_tendon = []
        self.len_t = None
        # self.nodes_tendon_id = None
        self.tendon_id = None
        # self.Conn_tendon = None
        # self.Conn_xi_norm_tendon = None
        self.NArr_tendon = None
        self.xi_norm_tendon = None
        
        self.disp = None
        self.vel = None
        self.acc = None
        self.extforce = None

        self.force = None
        self.CurrPos = None

        # material properties
        self.delta = None
        self.rho = None
        self.cnot = None
        self.LinDenT = None

        # pressure properties
        self.pnot = None
        self.b = None

        # other
        self.clamped_nodes = []
        self.top_node = None
        self.bottom_node = None

        # plot related
        self.plotcounter = 1

    def gen_nodes_on_tendon(self):
        """Generates a list of nodes which are on the tendon
        """

        print('Generating nodes on tendon')
        test_dir = []
        for i in range(self.ngores):
            xval = np.cos(2*np.pi/self.ngores*i)
            yval = np.sin(2*np.pi/self.ngores*i)
            test_dir.append(np.array([xval, yval]))

        self.tendon_id = []
        for i in range(len(self.pos)):
            self.tendon_id.append([])

        for i in range(len(self.pos)):
            vec = self.pos[i][0:2]
            pos_norm = np.sqrt(np.sum(vec**2))
            if pos_norm:
                node_unit = vec / pos_norm
                for j in range(self.ngores):
                    test_unit = test_dir[j]
                    diff = np.sqrt(np.sum((node_unit  - test_unit)**2, axis = 0))
                    if diff < 1e-8:
                        self.nodes_tendon.append(i)
                        # need to treat to and bottom node separately
                        self.tendon_id[i].append(j)
            else:
                ## fix for top and bottom node
                # if pos_norm is zero, must be on the z-axis, hence top or bottom node
                self.nodes_tendon.append(i)
                # need to treat to and bottom node separately
                self.tendon_id[i].append(-1)

        self.nodes_tendon = list(set(self.nodes_tendon))
        print('Total nodes on tendon', len(self.nodes_tendon))
        # print('tendon_id', self.tendon_id)

        ####
        # generate length elements for nodes associated with the tendons #
        ####
        print('# Generating length elements for nodes associated with the tendons #')
        t_set = set(self.nodes_tendon)
        Pos = self.pos
        self.len_t = np.zeros((len(Pos), 1))
        for i in range(len(self.edges)):
            edge = self.edges[i]
            if (set(edge)).issubset(t_set):
                l = Pos[edge[1]] - Pos[edge[0]]
                cp = np.sqrt(np.sum(l**2))

                # distribute length evenly over the endpoints
                j = edge[0]
                self.len_t[j] += cp/2
                j = edge[1]
                self.len_t[j] += cp/2
        # print(self.len_t)

        


        

    def save_state(self, filename):
        """Save Mesh with current values to disk using pickle
        :filename: TODO
        :returns: TODO
        """
        with open(filename, 'wb') as file_h:
          pickle.dump(self, file_h)

    def plot(self, dotsize=10, plot_node_text=True, highlight_bdry_nodes=True):
        """plot the mesh
        :returns: TODO

        """
        # print(Pos)

        Pos = self.pos

        dimension = len(Pos[0])


        if dimension==1:
            # Plot the mesh
            plt.scatter(Pos[:,0], Pos[:,1], s = dotsize, marker = '.', linewidth = 0, cmap='viridis')
        elif dimension==2:
            # Plot the mesh
            plt.scatter(Pos[:,0], Pos[:,1], s = dotsize, marker = '.', linewidth = 0, cmap='viridis')

            # highlight boundary nodes
            if highlight_bdry_nodes:
                plt.scatter(Pos[self.bdry_nodes,0], Pos[self.bdry_nodes,1], s = dotsize*1.1,  linewidth = 0 )

            # # text
            if plot_node_text:
                for i in range(0,len(Pos)):
                    plt.annotate(str(i), (Pos[i,0], Pos[i,1]))


            plt.axis('scaled')

        elif dimension==3:

            fig = plt.figure()
            # ax = fig.add_subplot(111, projection='3d')
            ax = Axes3D(fig)

            # Plot the mesh
            ax.scatter(Pos[:,0], Pos[:,1], Pos[:,2], s = dotsize, marker = '.', linewidth = 0, cmap='viridis')

            # highlight boundary nodes
            if highlight_bdry_nodes:
                ax.scatter(Pos[self.bdry_nodes,0], Pos[self.bdry_nodes,1], Pos[self.bdry_nodes,2], s = dotsize*1.1,  linewidth = 0 )

            # # text
            if plot_node_text:
                for i in range(0,len(Pos)):
                    ax.text(Pos[i,0], Pos[i,1], Pos[i,2], str(i))


            # ax.set_box_aspect((1,1,1))
            # ax.set_aspect('equal')
            # Create cubic bounding box to simulate equal aspect ratio
            # ax.pbaspect = [1.0, 1.0, 0.25]

            ## Fix aspect ratio
            ## def axisEqual3D(ax):
            # extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
            # sz = extents[:,1] - extents[:,0]
            # centers = np.mean(extents, axis=1)
            # maxsize = max(abs(sz))
            # r = maxsize/2
            # for ctr, dim in zip(centers, 'xyz'):
                # getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

            # print(np.amax(np.abs(Pos)))
            mx = np.amax(np.abs(Pos))

            XYZlim = [-mx, mx]
            ax.set_xlim3d(XYZlim)
            ax.set_ylim3d(XYZlim)
            ax.set_zlim3d(XYZlim)
            # ax.set_aspect('equal')
            ax.set_box_aspect((1, 1, 1))

            ## Plot properties
            ax.grid(False)
            # plt.axis('off')

            # plt.savefig('mesh.png', dpi=200, bbox_inches='tight')


        plt.show()

    def saveplot(self, count, show=False):
        """TODO: Docstring for plot_full.

        :arg1: TODO
        :returns: TODO

        """

        # box_L = 1.5
        box_L = 50 * 1.5
        # print('dim, cam angle', dim,  camera_angle)
        # camera_angle = [0, 0]
        camera_angle = [90, 0]

        plot_nodes = 0
        plot_nodes_on_tendon = 0
        dotsize = 0.5
        plot_mesh = 1
        linewidth=0.1

        print(self.plotcounter)
        filename = ('img/tc_%05d.png' % self.plotcounter)

        # Creating figure
        # ax = plt.axes(projection ="3d")
        fig = plt.figure()
        ax = Axes3D(fig)
        # Creating plot
        if plot_nodes:
            ax.scatter3D(self.CurrPos[:,0],self.CurrPos[:,1],self.CurrPos[:,2], color='green', s=dotsize)
        if plot_nodes_on_tendon:
            ax.scatter3D(self.CurrPos[self.nodes_tendon,0],self.CurrPos[self.nodes_tendon,1],self.CurrPos[self.nodes_tendon,2], color='blue', s=dotsize)

        if plot_mesh:
            # debug
            edges = self.edges

            V1 = edges[:,0]
            V2 = edges[:,1]

            P1 = self.CurrPos[V1]
            P2 = self.CurrPos[V2]
            ls =  [ [p1, p2] for p1, p2 in zip(P1,P2)] 
            lc = Line3DCollection(ls, linewidths=linewidth, colors='b')
            ax.add_collection(lc)

        # f_norm = np.sqrt(np.sum(self.force**2, axis=0))
        # ax.scatter3D(self.CurrPos[:,0],self.CurrPos[:,1],self.CurrPos[:,2], c=f_norm)
        # ax.axis('equal')
        ax.view_init(elev = camera_angle[0], azim = camera_angle[1])
        ax.set_xlim3d(-box_L, box_L)
        ax.set_ylim3d(-box_L, box_L)
        ax.set_zlim3d(-box_L, box_L)
        ax.set_box_aspect((1, 1, 1))
        # plt.title("simple 3D scatter plot")
        plt.savefig(filename, dpi=200, bbox_inches='tight')
        if show:
            plt.show()
        plt.close()

        
# def genmesh(P_bdry, meshsize, pygmsh_geom=None, msh_file = None, do_plot = True, dimension = 2, dotsize = 10, mesh_optimize=True):
def genmesh(P_bdry, meshsize, pygmsh_geom=None, msh_file = None, dimension = 2, mesh_optimize=True):
    """Generate a mesh from given polygon
    :P_bdry: an array of boundary points
    :meshsize: 
    :returns: pos and vol
    """

    # mesh
    if msh_file is None:
        if pygmsh_geom is None:
            # print('Generating mesh from poly')
            # with pygmsh.geo.Geometry() as geom:
            with pygmsh.occ.Geometry() as geom:
                polygon1 = geom.add_polygon(
                    P_bdry,
                    mesh_size= meshsize,
                )
                geom.add_physical(polygon1.surface, 'surface1')
                mesh = geom.generate_mesh()
        else:
            # print(pygmsh_geom)
            print('Generating from provided pygmsh_geom object.')
            mesh = pygmsh_geom.generate_mesh()
            print(mesh)

    else:
        print('Loading mesh from file: ', msh_file)
        mesh = meshio.read(msh_file)

    if mesh_optimize:
        # print('Optimizing mesh')
        # mesh = pygmsh.optimize(mesh, method="")
        # mesh = pygmsh.optimize(mesh, method="Netgen")
        # mesh = optimesh.optimize(mesh, "CVT (block-diagonal)", 1.0e-5, 100)
        pass

    Pos = mesh.points
    total_nodes = len(Pos)
    # print('Mesh nodes: ', total_nodes)

    # elements
    if dimension==1:
        T = mesh.get_cells_type("line")
    elif dimension==2:
        T = mesh.get_cells_type("3riangle")
    elif dimension == 3:
        # T = mesh.get_cells_type("tetra")
        T = mesh.get_cells_type("triangle")
        # all lines
        # edges = mesh.get_cells_type("line")
        edges = get_edges(T)
        temp = edges.flatten()
        bdry_nodes = list(set(temp))

    # print('Total T: ', len(T))

    area = np.zeros((total_nodes, 1))
    length = np.zeros((total_nodes, 1))


    # 2D

    if dimension==1:
        Pos = Pos[:,0:2]
    elif dimension==2:
        Pos = Pos[:,0:2]
    # otherwise gmsh produces 3d anyway

    # Generate nodal volume
    if dimension==1:
        for i in range(len(T)):
            l = np.sqrt(np.sum((Pos[T[i,1]] - Pos[T[i,0]])**2))

            # distribute area evenly over the vertices
            j = T[i,0]
            area[j] += l/2
            j = T[i,1]
            area[j] += l/2

    # elif dimension==2:
    elif dimension==3:
        for i in range(len(T)):
            u_x = Pos[T[i,1], 0] - Pos[T[i,0], 0]
            u_y = Pos[T[i,1], 1] - Pos[T[i,0], 1]

            v_x = Pos[T[i,2], 0] - Pos[T[i,0], 0]
            v_y = Pos[T[i,2], 1] - Pos[T[i,0], 1]

            # cross product, half
            cp = 0.5 * abs(u_x * v_y - u_y * v_x)

            # distribute area evenly over the vertices
            j = T[i,0]
            area[j] += cp/3
            j = T[i,1]
            area[j] += cp/3
            j = T[i,2]
            area[j] += cp/3

        ## get the length elements
        # for i in range(len(edges)):
            # l = Pos[edges[i,1]] - Pos[edges[i,0]]
            # cp = np.sqrt(np.sum(l**2))

            # # distribute length evenly over the endpoints
            # j = edges[i,0]
            # length[j] += cp/2
            # j = edges[i,1]
            # length[j] += cp/2

    ## removing nodes with zero volume
    nodelist_zero_vol = np.where(area == 0)[0]
    if len(nodelist_zero_vol) > 0:
        print('Caution: there are nodes with zero volume: ', nodelist_zero_vol)
        print('Bad node that participates in generating elements (and cannot be safely removed) are ', set(range(1, len(Pos))).difference(set(T.flatten())))

        #u delete the nodes
        # Pos[nodelist_zero_vol] = []
        # area[nodelist_zero_vol] = []
        print(len(Pos))
        print(len(area))
        # Pos = np.delete(Pos, nodelist_zero_vol, axis=0)
        # area = np.delete(area, nodelist_zero_vol, axis=0)
        for jj in range(len(nodelist_zero_vol)):
            badnode = nodelist_zero_vol[jj]
            Pos = np.delete(Pos, [badnode], axis=0)
            area = np.delete(area, [badnode], axis=0)
            # shift the index of every node bigger than that by 1
            for aa in range(len(T)):
                for bb in range(len(T[aa])):
                    if T[aa][bb] > badnode:
                        T[aa][bb] -= 1
            # do the same for edges
            for aa in range(len(edges)):
                for bb in range(len(edges[aa])):
                    if edges[aa][bb] > badnode:
                        edges[aa][bb] -= 1
        ## rename an
        print(len(Pos))
        print(len(area))

    # nodelist_zero_length = np.where(length == 0)[0]
    # if len(nodelist_zero_length) > 0:
        # print('Caution: there are nodes with zero length: ', nodelist_zero_length)
        # print('Bad node that participates in generating elements (and cannot be safely removed) are ', set(range(1, len(Pos))).difference(set(edges.flatten())))

    return ReturnValue(Pos, area, bdry_nodes, T, edges)

