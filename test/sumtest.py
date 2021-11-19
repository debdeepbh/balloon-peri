from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print('rank', rank, flush=True)

total_nodes = 1000
# force = np.zeros((total_nodes, 3))
force = np.zeros(total_nodes, dtype='d')

for i in range(rank, total_nodes, size):
    # force[i] += [(i+1)**2, (i+1)**2, (i+1)**2]
    force[i] += (i+1)**2


if rank == 0:
    data = np.zeros(total_nodes, dtype='d')
else:
    data = None

print('rank before', rank, force)
comm.Allreduce(MPI.IN_PLACE, force, op=MPI.SUM)
# comm.Reduce([force, MPI.DOUBLE], [data, MPI.DOUBLE], op=MPI.SUM, root=0)
print('rank after', rank, force)
# print('after', rank, data)
