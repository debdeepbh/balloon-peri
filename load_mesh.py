import  numpy as np
import pickle

dt = 2e-3
# timesteps = 50
timesteps = 1

# load Mesh
Mesh = pickle.load( open( "data/Meshdump.pkl", "rb" ) )
total_nodes = len(Mesh.pos)
Mesh.disp = Mesh.vel = Mesh.acc = Mesh.extforce =  np.zeros((total_nodes,3))
Mesh.CurrPos = Mesh.force =  np.zeros((total_nodes,3))

## Connectivity to NbdArr
Mesh.NArr = []
Mesh.xi = []
Mesh.xi_norm = []
for i in range(len(Mesh.pos)):
    Mesh.NArr.append([])

for i in range(len(Mesh.Conn)):
    v = Mesh.Conn[i]
    p = v[0]
    q = v[1]
    Mesh.NArr[p].append(q)
    Mesh.NArr[q].append(p)
# print(Mesh.NArr)

# Mesh.xi = Mesh.pos[ arr]
for i in range(len(Mesh.pos)):
    Mesh.NArr.append([])
# Mesh.xi = Mesh.pos[ arr]
print(Mesh.xi)
# Mesh.xi_norm

# Material properties
Mesh.rho = 10

## Initial data
Mesh.disp += [0, 0, 0]
Mesh.vel += [0, 0, 0]
Mesh.acc += [0, 0, 0]
Mesh.extforce += [0, 0, 0]

def get_peridynamic_force(Mesh):
    """Compute the peridynamic force
    :Mesh: TODO
    :returns: TODO
    """

    force = np.zeros((total_nodes, 3))

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
