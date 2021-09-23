import pickle
from multiprocessing import Pool
# from pathos.multiprocessing import ProcessingPool as Pool
from sys import argv

fc = int(argv[1])

filename = ("output/mesh_%05d.pkl" % fc)
Mesh = pickle.load( open( filename, "rb" ) )
Mesh.saveplot(fc, show=True)
    

