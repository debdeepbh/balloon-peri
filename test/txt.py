from mpi4py import MPI

comm = MPI.COMM_WORLD
nprocs = comm.Get_size()
rank   = comm.Get_rank()

if rank == 0:
   data = 'Hello!'
   comm.send(data, dest=nprocs-1, tag=1)
elif rank == nprocs-1:
   data = comm.recv(source=0, tag=1)
   print ('Rank ', rank, ' received ', data)
