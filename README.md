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

New:
```
mesh/3d_sphere_forloop.geo
```

Old (don't use):
```
gmsh mesh/3d_sphere_unit.geo -2
```
(`mesh/3d_sphere_forloop.geo` has nodes on the tendons)

- Generate neighborhood array (set `ngores` and peridynamic neighborhood `delta` here)
```
python3 gen_nbdarr.py
```
* Load the mesh and simulate (and set parameters) with 
```
python3 load_mesh.py
```
- Resuming is possible by setting `resume = True` in `load_mesh.py` and running `python3 load_mesh.py` again.


# Baby case
- [ ] Nodes on tendon
	- [-] (Need to treat top and bottom nodes separately) Generate neighborhood for nodes on tendon
	- [ ] Implement 1-d peridynamics on those
- [ ] take spring constant to be 2d plastic sheet
- [ ] Take reference length to be either the flat distance or the distance along the sphere.
- [x] Spherical shell, generate nodes mesh in 3d
- [x] Generate area elements from the mesh
- [x] Add air pressure based on the triangle area
- [x] clamp at one point
- [x] simulate using newton's law (Cauchy momentum with body force as pressure) 
