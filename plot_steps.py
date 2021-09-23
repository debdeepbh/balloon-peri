import pickle
from multiprocessing import Pool
# from pathos.multiprocessing import ProcessingPool as Pool

fc = 1

with open('data/last_counter', 'r') as f:
    lc = int(f.read())

def savefig(count):
    filename = ("output/mesh_%05d.pkl" % count)
    Mesh = pickle.load( open( filename, "rb" ) )
    Mesh.saveplot(count)
    
# parallel formulation
a_pool = Pool()
a_pool.map(savefig, range(fc, lc))
a_pool.close()

