rad = 1;
// Meshsize
lc = 0.1;
ngores = 30;

// Center, north and south poles
//+
Point(1) = {0, 0, 0, lc};
Point(2) = {0, 0, rad, lc};
Point(3) = {0, 0, -rad, lc};
//+

// Cannot have arcs of angle pi, so dividing into two
For t In {1:ngores}
    Point(3+t) = {rad*Cos(t*2*Pi/ngores), rad*Sin(t*2*Pi/ngores), 0, lc};
    Circle(2*(t-1)+1) = {2, 1, 3+t};
    Circle(2*(t-1)+2) = {3+t, 1, 3};
EndFor

For c In {1:(ngores-1)}
    Curve Loop(c) = {2*(c-1)+1, 2*(c-1)+2, -(2*c+2), -(2*c+1)};
    Surface(c) = {c};
EndFor
// Last surface
Curve Loop(ngores) = {2*(ngores-1)+1, 2*(ngores-1)+2, -2, -1};
Surface(ngores) = {ngores};

