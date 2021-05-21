import  numpy as np
import pickle
import matplotlib.pyplot as plt
# from pathos.multiprocessing import ProcessingPool as Pool

dt = 2e-2
timesteps = 50
# timesteps = 1


# load Mesh
Mesh = pickle.load( open( "data/Meshdump.pkl", "rb" ) )
total_nodes = len(Mesh.pos)
Mesh.disp = Mesh.vel = Mesh.acc = Mesh.extforce =  np.zeros((total_nodes,3))
Mesh.CurrPos = Mesh.force =  np.zeros((total_nodes,3))

bottom_node = np.argmin(Mesh.pos[:,2])  
top_node = np.argmax(Mesh.pos[:,2])  
print('bottom node', bottom_node)

# nodes to clamp
clamped_nodes = []
# clamped_nodes = [bottom_node]

# Material properties

E = 72e20
nu = 1/3
cnot = 6*E/( np.pi * (Mesh.delta**3) * (1 - nu))

Mesh.cnot = cnot
print('cnot', cnot)
Mesh.pnot = 100
Mesh.rho = 10
Mesh.b = 50

# plot properties
box_L = 1.5
modulo = 5

## Connectivity to NbdArr
Mesh.NArr = []
# Mesh.xi = []
Mesh.xi_norm = []
for i in range(len(Mesh.pos)):
    Mesh.NArr.append([])
    Mesh.xi_norm.append([])

for i in range(len(Mesh.Conn)):
    v = Mesh.Conn[i]
    p = v[0]
    q = v[1]
    d = Mesh.Conn_xi_norm[i]
    Mesh.NArr[p].append(q)
    Mesh.NArr[q].append(p)

    Mesh.xi_norm[p].append(d)
    Mesh.xi_norm[q].append(d)

# print(Mesh.NArr)
# print(Mesh.xi_norm)

## Initial data
# Mesh.disp += [0, 0, 0]
# Mesh.disp += [0, 0, 0.5]
Mesh.disp[top_node] += [0, 0, 0.5]
# Mesh.vel += [0, 0, 1]
# Mesh.vel[top_node] += [0, 0, 0.1]
# Mesh.acc += [0, 0, 0]
# Mesh.acc[top_node] += [0, 0, 2e20]
# Mesh.extforce += [0, 0, 0]
# Mesh.extforce[top_node] += [0, 0, 2e20]

def get_peridynamic_force(Mesh):
    """Compute the peridynamic force
    :Mesh: TODO
    :returns: TODO
    """
    force = np.zeros((total_nodes, 3))

    for i in range(total_nodes):
    # def one_row(i):
        nbrs = Mesh.NArr[i]

        n_etapxi = Mesh.CurrPos[nbrs] - Mesh.CurrPos[i]
        n_etapxi_norm =  np.sqrt(np.sum(n_etapxi**2, axis=1))
        n_strain = np.array([(n_etapxi_norm - Mesh.xi_norm[i]) / Mesh.xi_norm[i]]).transpose()

        # replace negative values by zero
        n_strain[n_strain < 0] = 0

        # check if dividing by zero
        # n_unit_dir = n_etapxi / n_etapxi_norm
        n_unit_dir = np.array([p/np for p,np in zip(n_etapxi , n_etapxi_norm)])
        n_vol = Mesh.vol[nbrs]

        # print(n_unit_dir)
        # print(n_strain)
        # print('etapxi', n_CurrPos_norm)
        # print('xi', Mesh.xi_norm[i])
        # print(n_CurrPos)
        # print(n_vol)
        # print(n_strain)
        # print(n_unit_dir)

        nsum_force = np.sum(Mesh.cnot * n_strain * n_unit_dir * n_vol, axis=0)

        # return nsum_force
        force[i,:] = nsum_force

    ## parallel attempt: slow for small number of nodes
    # a_pool = Pool()
    # force = a_pool.map(one_row, range(total_nodes)) 
    # force = np.array(force)
    # print(force)

    return force
    
class PressureQ(object):
    """docstring for PressureQ"""
    def __init__(self, pforce, CurrArea, CurrNormal):
        self.pforce = pforce
        self.CurrArea = CurrArea
        self.CurrNormal = CurrNormal
        

def get_pressure(Mesh):
    """TODO: Docstring for get_pressure.

    :Mesh: TODO
    :returns: TODO

    """
    area = np.zeros((total_nodes, 1))
    unormal = pforce = np.zeros((total_nodes, 3))

    T = Mesh.T
    Pos = Mesh.CurrPos

    for i in range(len(T)):
        u_x = Pos[T[i,1], 0] - Pos[T[i,0], 0]
        u_y = Pos[T[i,1], 1] - Pos[T[i,0], 1]
        u_z = Pos[T[i,1], 2] - Pos[T[i,0], 2]

        v_x = Pos[T[i,2], 0] - Pos[T[i,0], 0]
        v_y = Pos[T[i,2], 1] - Pos[T[i,0], 1]
        v_z = Pos[T[i,2], 2] - Pos[T[i,0], 2]

        # z-value of the centroid of the triangle
        z_cent = (Pos[T[i,0],2] + Pos[T[i,1],2] + Pos[T[i,2],2])/3

        # cross product
        crossp = np.cross( [u_x, u_y, u_z], [v_x, v_y, v_z])
        crossp_norm = np.sqrt(np.sum(crossp **2))
        udir = crossp / crossp_norm

        # full area
        cp = 0.5 * crossp_norm

        # full force due to pressure: (P_0 + bz).Area.unit_normal
        pf = (Mesh.pnot + Mesh.b*z_cent) * cp  * udir

        # distribute area evenly over the vertices
        j = T[i,0]
        area[j] += cp/3
        j = T[i,1]
        area[j] += cp/3
        j = T[i,2]
        area[j] += cp/3
        
        #distribute unit normal evenly over vertices
        j = T[i,0]
        unormal[j] += udir/3
        j = T[i,1]
        unormal[j] += udir/3
        j = T[i,2]
        unormal[j] += udir/3

        #distribute pressure-force evenly over vertices
        j = T[i,0]
        pforce[j] += pf/3
        j = T[i,1]
        pforce[j] += pf/3
        j = T[i,2]
        pforce[j] += pf/3

    return PressureQ(pforce, CurrArea=area, CurrNormal=unormal)

top_nbr = Mesh.NArr[top_node]
print('top node neighbors', top_nbr)
plotcounter = 1
for t in range(timesteps):
    print('t', t)

    # initial update
    Mesh.disp += dt * Mesh.vel + (dt * dt * 0.5) * Mesh.acc
    Mesh.CurrPos = Mesh.pos + Mesh.disp

    ## compute force
    force = get_peridynamic_force(Mesh)

    print(force[top_node,:])
    print(force[top_nbr,:])
    # force += get_pressure(Mesh)
    # P = get_pressure(Mesh)
    # force += P.pforce
    # print(P.CurrNormal)

    # final update
    Mesh.acc = (1 / Mesh.rho) * (Mesh.force + Mesh.extforce)
    #	# velocity
    #	u0dot_univ{i} = uolddot_univ{i} + dt * 0.5 * uolddotdot_univ{i}
    #+ dt * 0.5 * u0dotdot_univ{i}; Mesh.vel = Mesh.vel  + (dt * 0.5)
    #*
    # Mesh.acc + (dt * 0.5) * Mesh.acc_old;
    temp_acc = Mesh.acc
    temp_acc += Mesh.acc   # Now temp_acc = (acc + acc_old)
    temp_acc *= (0.5 * dt) # now temp_acc = (dt*0.5) *(acc + acc_old)
    Mesh.vel += temp_acc

    # clamped node
    for i in range(len(clamped_nodes)):
        cnode = clamped_nodes[i]
        Mesh.disp[cnode] = [0, 0, 0]
        Mesh.vel[cnode] = [0, 0, 0]
        Mesh.acc[cnode] = [0, 0, 0]

    #plot
    if t % modulo:
        # Creating figure
        ax = plt.axes(projection ="3d")
        # Creating plot
        ax.scatter3D(Mesh.CurrPos[:,0],Mesh.CurrPos[:,1],Mesh.CurrPos[:,2], color='green')
        # ax.axis('equal')
        ax.set_xlim3d(-box_L, box_L)
        ax.set_ylim3d(-box_L, box_L)
        ax.set_zlim3d(-box_L, box_L)
        plt.title("simple 3D scatter plot")
        plt.savefig('img/tc_%05d.png' % t)
        # plt.show()
        plt.close()
        plotcounter += 1
