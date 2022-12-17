import pickle
from multiprocessing import Pool
# from pathos.multiprocessing import ProcessingPool as Pool
from sys import argv

fc = 1
with open('data/last_counter', 'r') as f:
    lc = int(f.read())

if len(argv) == 3:
    fc = int(argv[1])
    lc = int(argv[2])


# box_L = None
# box_L = 300
box_L = 100
center_at_mean = True
# center_at_mean = False

def savefig(count):
    filename = ("output/mesh_%05d.pkl" % count)
    Mesh = pickle.load( open( filename, "rb" ) )
    Mesh.saveplot(count, box_L=box_L, center_at_mean=center_at_mean)
    
# parallel formulation
a_pool = Pool()
a_pool.map(savefig, range(fc, lc+1))
a_pool.close()

