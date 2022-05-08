
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
