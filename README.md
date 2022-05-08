# Balloon with peridynamics

Dynamics of a weather balloon launch using a state-based peridynamic model for thin sheet

# Demonstation

- State-based model with poisson ratio 0.4
![vid](demo/state-based-nu-0.4.mp4)
- State-based model with poisson ratio 0.46 and velocity-dependent damping
![vid](demo/state-based-nu-0.46-damping-1.mp4)


# Installation
```
# gmsh package
sudo apt-get install gmsh

# create a virtual environment
python3 -m venv env
source env/bin/activate

# install python packages
pip3 install -r requirements.txt
```

# Running

- Generate the surface mesh (set meshsize etc within the file) with

Full-size (150 m as half-diameter):
```
gmsh mesh/3d_sphere_forloop_big.geo -2 
```

**Unit sphere** with tendons on vertical planes 
```
gmsh mesh/3d_sphere_forloop.geo -2
```

Old (don't use):
```
gmsh mesh/3d_sphere_unit.geo -2
```
(`mesh/3d_sphere_forloop.geo` has nodes on the tendons)

- Generate neighborhood array (set `msh` file name, `ngores`, and peridynamic horizon size `delta` here)
```
python3 gen_nbdarr.py
```
* Load the mesh and simulate (and set parameters) with 
```
python3 load_mesh.py
```
- Resuming is possible by setting `resume = True` in `load_mesh.py` and running `python3 load_mesh.py` again.

- Convert saved files into images using
```
python3 plot_steps.py
```

# MPI

- Install 

```
pip3 install mpi4py
```

- Run

```
mpiexec -np 4 python3 mpi_compute.py
```
