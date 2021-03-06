//+
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, 0, 1, 1.0};
Point(3) = {0, 0, -1, 1.0};
Point(4) = {0, 1, 0, 1.0};
Point(5) = {0, -1, 0, 1.0};
Point(6) = {1, 0, 0, 1.0};
Point(7) = {-1, 0, 0, 1.0};
//+
Circle(1) = {2, 1, 4};
Circle(2) = {4, 1, 3};
//Circle(3) = {3, 1, 5};
//Circle(4) = {5, 1, 2};
//Circle(5) = {2, 1, 7};
//Circle(6) = {7, 1, 3};
//
//Circle(7) = {4, 1, 7};
//+
Curve Loop(1) = {1, 2};
//+
//Curve Loop(2) = {2, -6, -5, 1};
//+
//Surface(1) = {2};
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Curve{1}; Curve{2}; Layers{10}; Recombine;
}
