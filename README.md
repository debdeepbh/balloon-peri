# balloon with peridynamics

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
mesh/3d_sphere_forloop_big.geo -2 
```

Unit sphere with tendons on vertical planes
```
mesh/3d_sphere_forloop.geo -2
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

# Todo

- [ ] Incorrect cnot in 2D; (correct expression for 2d constant `cnot` and `dt=1e-4` works)
- [ ] `cnot_1d`: the constant and integrand do not match

- [ ] Nodes on tendon
	- [ ] What is the correct modulus to use for 1d peridynamics?
	- [x] (Need to treat top and bottom nodes separately) Generate neighborhood for nodes on tendon
	- [x] compute 1d peridynamic constant: $3 \lambda/ \delta^3$ where $\lambda$ is the Lame coefficient.
	- [x] Implement 1-d peridynamics on tendon


- [x] Take reference length to be either the flat distance or the distance along the sphere.
- [x] Spherical shell, generate nodes mesh in 3d
- [x] Generate area elements from the mesh
- [x] Add air pressure based on the triangle area
- [x] clamp at one point
- [x] simulate using newton's law (Cauchy momentum with body force as pressure) 
