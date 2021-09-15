//+
radius = 150/Pi;
lc = 5;

Point(1) = {0, -0, 0, lc};
//+
SetFactory("OpenCASCADE");
// mesh size factor = 1 => 171 nodes
//Mesh.CharacteristicLengthFactor = 1;
Mesh.CharacteristicLengthFactor = 0.75;
//Mesh.CharacteristicLengthFactor = 0.5;
//Sphere(1) = {0, 0, 0, radius, -Pi/2, Pi/2, 2*Pi};
Sphere(1) = {0, 0, 0, radius};
Surface Loop(4) = {1};
Physical Surface(1) = {1};
