import  numpy as np
import pickle
import copy
import matplotlib.pyplot as plt
# from pathos.multiprocessing import ProcessingPool as Pool
# from matplotlib.collections import LineCollection

# dt = 2e-3
# dt = 1e-1
dt = 1e-4
timesteps = 50000
# timesteps = 1

# resume = True
resume = False

allow_damping = 0

# plot properties
# modulo = 50
modulo = 100


if resume:
    # load Mesh
    Mesh = pickle.load( open( "savedata/Mesh_saved.pkl", "rb" ) )
    total_nodes = len(Mesh.pos)
else:
    # load Mesh
    Mesh = pickle.load( open( "data/Meshdump.pkl", "rb" ) )
    total_nodes = len(Mesh.pos)
    Mesh.disp =  np.zeros((total_nodes,3))
    Mesh.vel =  np.zeros((total_nodes,3))
    Mesh.acc =  np.zeros((total_nodes,3))
    Mesh.extforce =  np.zeros((total_nodes,3))
    Mesh.CurrPos = np.zeros((total_nodes,3))
    Mesh.force =  np.zeros((total_nodes,3))

    Mesh.bottom_node = np.argmin(Mesh.pos[:,2])  
    Mesh.top_node = np.argmax(Mesh.pos[:,2])  
    #print('bottom node', Mesh.bottom_node)

    # z-value of the bottom of the balloon
    z_0 = Mesh.pos[Mesh.bottom_node,2]
    # print('z_0', z_0)

    # nodes to clamp
    # clamped_nodes = []
    Mesh.clamped_nodes = [Mesh.bottom_node]

    # Material properties

    # E = 72e9
    E = 0.3e9 # Polyethylene (low density) LDPE 
    LinDenT =0.015205 # kg/m	# linear density of the tendon
    ETape = 761560
    nu = 1/3
    ethickness = 38.1e-06
    # BalloonFilmWeightArealDensity = ethickness * 920 * gravity	# weight of the film per unit area
    BalloonFilmMassArealDensity = ethickness * 920 

    # 2d constant: from my granular paper
    cnot = 6*E/( np.pi * (Mesh.delta**3) * (1 - nu))

    damping_coeff = 10

    # division between 1e8:too stiff and 1e10:too loose
    Mesh.cnot = cnot
    # Mesh.cnot = cnot /1e3
    print('cnot', cnot)


    # Mesh.rho = 920 # LDPE 920 kg/m^3
    ## for 2d peridynamics, we want rho to be Mass/Area
    Mesh.rho = BalloonFilmMassArealDensity

    # Is this the right unit to use?
    # Mesh.tendon_modulus = LinDenT
    Mesh.LinDenT = LinDenT
    Mesh.tendon_modulus = ETape
    # Mesh.cnot_tendon = 3 * Mesh.tendon_modulus / (Mesh.delta**3)
    Mesh.cnot_tendon = 2 * Mesh.tendon_modulus / (Mesh.delta**2)

    print('cnot_tendon', Mesh.cnot_tendon)

    # mass of the nodes
    Mesh.mass = Mesh.rho * Mesh.vol + Mesh.LinDenT * Mesh.len_t

    Mesh.allow_damping = allow_damping
    Mesh.damping_coeff = damping_coeff

    ## pressure properties
    # gravity
    g_val = -10
    # Mesh.pnot = 100
    # Mesh.b = 500
    # From Frank's code: tauVol = V_d/TargetVolume ;% 0.0125; b_d = 0.084887 % N/m^3 buoyancy = b_d * tauVol  ;
    Mesh.b = 0.104314581941182
    Mesh.pnot = 1

    ## Connectivity to NbdArr
    print('Converting connectivity to NbdArr')
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

    # NbdArr for tendons
    print('Generating tendon connectivity')
    Mesh.NArr_tendon = []
    Mesh.xi_norm_tendon = []

    print('top nodes', Mesh.top_node)
    print('bottom nodes', Mesh.bottom_node)

    for i in range(len(Mesh.pos)):
        Mesh.NArr_tendon.append([])
        Mesh.xi_norm_tendon.append([])

    for i in range(len(Mesh.Conn)):
        v = Mesh.Conn[i]
        p = v[0]
        q = v[1]
        d = Mesh.Conn_xi_norm[i]

        p_id = Mesh.tendon_id[p]
        q_id = Mesh.tendon_id[q]

        if len(p_id) and len(q_id):
            # if tendon id matches
            if (p_id[0] == q_id[0]):
                # print('tendon matches')

                # need to treat top and bottom nodes separately
                Mesh.NArr_tendon[p].append(q)
                Mesh.NArr_tendon[q].append(p)

                Mesh.xi_norm_tendon[p].append(d)
                Mesh.xi_norm_tendon[q].append(d)

            # top and bottom nodes have tendon_id = -1
            # any nonzero tendon id is a neighbor, if close enough
            if (p_id[0] == (-1)) or (q_id[0] == (-1)) :
                # print('top or bottom node', p, q)
                Mesh.NArr_tendon[p].append(q)
                Mesh.NArr_tendon[q].append(p)
                Mesh.xi_norm_tendon[p].append(d)
                Mesh.xi_norm_tendon[q].append(d)

        else:
            pass
            # print('is [] for', p,' and', q)
    print('Done')

    ## Initial data
    # Mesh.disp += [0, 0, 0]
    # Mesh.disp += [0, 0, 0.5]
    # Mesh.disp[Mesh.top_node] += [0, 0, 0.5]
    # Mesh.vel += [0, 0, 1]
    # Mesh.vel[Mesh.top_node] += [0, 0, 1e2]
    # Mesh.acc += [0, 0, 0]
    # Mesh.acc[Mesh.top_node] += [0, 0, 1e5]
    # Mesh.extforce += [0, 0, 0]
    ## gravity (not density anymore)
    Mesh.extforce = np.c_[
            np.zeros(total_nodes),
            np.zeros(total_nodes),
            g_val * Mesh.mass
            ]
    # print(Mesh.extforce)


def get_peridynamic_force_density(Mesh):
    """Compute the peridynamic force density
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
        nsum_force = np.sum(Mesh.cnot * n_strain * n_unit_dir * n_vol, axis=0)

        # return nsum_force
        force[i,:] = nsum_force

    ## parallel attempt: slow for small number of nodes
    # a_pool = Pool()
    # force = a_pool.map(one_row, range(total_nodes)) 
    # force = np.array(force)
    # print(force)

    return force

def get_peridynamic_force_density_tendon(Mesh):
    """Compute the peridynamic force density
    :Mesh: TODO
    :returns: TODO
    """
    force = np.zeros((total_nodes, 3))

    for i in range(total_nodes):
    # def one_row(i):
        nbrs = Mesh.NArr_tendon[i]
        nsum_force = 0

        if len(nbrs):
            n_etapxi = Mesh.CurrPos[nbrs] - Mesh.CurrPos[i]
            n_etapxi_norm =  np.sqrt(np.sum(n_etapxi**2, axis=1))
            n_strain = np.array([(n_etapxi_norm - Mesh.xi_norm_tendon[i]) / Mesh.xi_norm_tendon[i]]).transpose()

            # replace negative values by zero
            n_strain[n_strain < 0] = 0

            # check if dividing by zero
            # n_unit_dir = n_etapxi / n_etapxi_norm
            n_unit_dir = np.array([p/np for p,np in zip(n_etapxi , n_etapxi_norm)])
            n_vol = Mesh.len_t[nbrs]
            nsum_force = np.sum(Mesh.cnot_tendon * n_strain * n_unit_dir * n_vol, axis=0)

        # return nsum_force
        force[i,:] = nsum_force

    return force
    
class PressureQ(object):
    """docstring for PressureQ"""
    def __init__(self, pforce, CurrArea, CurrNormal):
        self.pforce = pforce
        self.CurrArea = CurrArea
        self.CurrNormal = CurrNormal
        

def get_pressure(Mesh):
    """ Returns the force due to pressure
    :Mesh: TODO
    :returns: TODO
    """
    area = np.zeros((total_nodes, 1))
    unormal = np.zeros((total_nodes, 3))
    pforce = np.zeros((total_nodes, 3))

    T = Mesh.T
    Pos = Mesh.CurrPos

    for i in range(len(T)):
        # print('i', i)
        # print('T[i,0]', T[i,0])
        # print('T[i,1]', T[i,1])
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
        # pf = (Mesh.pnot + Mesh.b*z_cent) * cp  * udir
        # editing to start counting height from the bottom of the balloon
        pf = (Mesh.pnot + Mesh.b*(z_cent - z_0) ) * cp  * udir

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

top_nbr = Mesh.NArr[Mesh.top_node]
# print('top node neighbors', top_nbr)
bottom_nbr = Mesh.NArr[Mesh.bottom_node]
# print('bottom node neighbors', bottom_nbr)

# print('Initial mean disp', np.mean(Mesh.disp, axis = 0))

for t in range(timesteps):
    # initial update
    Mesh.disp += dt * Mesh.vel + (dt * dt * 0.5) * Mesh.acc
    Mesh.CurrPos = Mesh.pos + Mesh.disp

    ### [Abandoned] compute force density
    #Mesh.force = get_peridynamic_force_density(Mesh)
    #Mesh.force += get_peridynamic_force_density_tendon(Mesh)
    #P = get_pressure(Mesh)
    ## convert force to force density
    #Mesh.force += ( P.pforce / Mesh.vol)

    ## [Full force, not the density] compute force instead of force density
    Mesh.force = get_peridynamic_force_density(Mesh) * Mesh.vol
    Mesh.force += (get_peridynamic_force_density_tendon(Mesh) * Mesh.len_t)
    Mesh.force += get_pressure(Mesh).pforce
    Mesh.force += Mesh.extforce

    if Mesh.allow_damping:
        Mesh.force -= Mesh.damping_coeff * Mesh.vel

    ## acceleration from force density
    # Mesh.acc = (1 / Mesh.rho) * (Mesh.force + Mesh.extforce)
    ## acceleration from force
    Mesh.acc = (1 / Mesh.mass) * Mesh.force
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
    for i in range(len(Mesh.clamped_nodes)):
        cnode = Mesh.clamped_nodes[i]
        Mesh.disp[cnode] = [0, 0, 0]
        Mesh.vel[cnode] = [0, 0, 0]
        Mesh.acc[cnode] = [0, 0, 0]

    #plot
    if (t % modulo)==0:
        print('c', Mesh.plotcounter)
        filename = ('output/mesh_%05d.pkl' % Mesh.plotcounter)
        Mesh.save_state(filename)

        with open('data/last_counter', 'w') as f:
            f.write(str(Mesh.plotcounter))

        # save occasionally for resuming
        if (Mesh.plotcounter % 5)==0:
            # save last 
            Mesh.save_state('savedata/Mesh_saved.pkl')

        Mesh.plotcounter += 1
