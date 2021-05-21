# balloon with peridynamics

# Installation
```
# python packages
python3 -m venv env
source env/bin/activate
pip3 install -r requirements.txt
```

# Running

- Generate the mesh with

```
gmsh 3d_sphere_unit.geo -2
```
- Generate neighborhood array
```
python3 gen_nbdarr.py
```
* Load the mesh with 
```
python3 load_mesh.py
```


# Baby case
[x] Spherical shell, generate nodes mesh in 3d
[x] Generate area elements from the mesh
[ ] Add air pressure based on the triangle area
[ ] clamp at one point
[ ] take spring constant to be 2d plastic sheet
[ ] simulate using newton's law (Cauchy momentum with body force as pressure) 
