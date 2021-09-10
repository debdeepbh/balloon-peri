rad = 1;
lc = 0.1;
ngores = 8;

//SetFactory("OpenCASCADE");
//Mesh.CharacteristicLengthFactor = 0.1;

//+
Point(1) = {0, 0, 0, lc};
Point(2) = {0, 0, rad, lc};
Point(3) = {0, 0, -rad, lc};
Point(4) = {rad, 0, 0, lc};
Point(5) = {rad*Cos(2*Pi/ngores), rad*Sin(2*Pi/ngores), 0, lc};
//+
Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 3};
//+
Circle(3) = {2, 1, 5};
Circle(4) = {5, 1, 3};
//+
Circle(5) = {4, 1, 5};
//+
Curve Loop(1) = {3, 4, -2, -1};
//+
Surface(1) = {1};
//+
//+
//+
//Extrude {{0, 0, 1}, {0, 0, 0}, Pi/4} {
//  Curve{3}; Curve{4}; Layers{2}; Recombine;
//}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, -Pi/2} {
  Curve{1}; Curve{2}; Recombine;
}
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Curve{3}; Curve{4}; Recombine;
}
//+
Curve Loop(2) = {12, 15, -9, -6};
//+
Surface(18) = {2};
