# Balloon with peridynamics

Launch dynamics of a weather balloon using a state-based peridynamic model for thin sheet

# Demonstation

- State-based model with poisson ratio 0.4
![img](demo/nodamping.gif)
- State-based model with poisson ratio 0.46 and velocity-dependent damping
![img](demo/damping.gif)


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
gmsh mesh/3d_sphere_forloop_big.geo -3 
```

**Unit sphere** with tendons on vertical planes 

```
gmsh mesh/3d_sphere_forloop.geo -3
```

Old (don't use):

```
gmsh mesh/3d_sphere_unit.geo -3
```

(`mesh/3d_sphere_forloop.geo` has nodes on the tendons)

- Generate neighborhood array (set `msh` file name, `ngores`, and peridynamic horizon size `delta` here. Alternatively, set `neighbor_type = 'nearest_neighbor'` to use mesh edge information as the neighborhood information.)

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

# Parallel implementation with MPI

Using parallel time integration steps:

- Install 

```
pip3 install mpi4py
```

- Run

```
mpiexec -np 4 python3 mpi_compute.py
```
