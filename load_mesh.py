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

# Material properties
Mesh.rho = 10
Mesh.cnot = 10e3

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
Mesh.disp += [0, 0, 0]
Mesh.vel += [1, 0, 0]
Mesh.acc += [0, 0, 0]
Mesh.extforce += [0, 0, 0]

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
    

for t in range(timesteps):
    print('t', t)

    # initial update
    Mesh.disp += dt * Mesh.vel + (dt * dt * 0.5) * Mesh.acc
    Mesh.CurrPos = Mesh.pos + Mesh.disp

    ## compute force
    force = get_peridynamic_force(Mesh)

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

    #plot
    # Creating figure
    ax = plt.axes(projection ="3d")
     
    # Creating plot
    ax.scatter3D(Mesh.CurrPos[:,0],Mesh.CurrPos[:,1],Mesh.CurrPos[:,2], color='green')

    # ax.axis('equal')
    plt.title("simple 3D scatter plot")
    plt.savefig('img/tc_%05d.png' % t)

    # plt.show()
    plt.close()
