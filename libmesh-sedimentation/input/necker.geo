// Gmsh project created on Fri Oct 28 15:48:31 2016
xmin = 0.0;
xmax = 18.0;
ymin = 0.0;
ymax = 2.0;
xlock = 1.0;
lc = 0.1;

Point(1) = {xmin, ymin, 0.0, lc};
Point(2) = {xlock, ymin, 0.0, lc};
Point(3) = {xlock, ymax, 0.0, lc};
Point(4) = {xmin , ymax, 0.0, lc};
Point(5) = {xmax,  ymin, 0.0, lc};
Point(6) = {xmax , ymax, 0.0, lc};//+


Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {2, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 3};
//+
Line Loop(8) = {4, 1, 2, 3};
//+
Plane Surface(9) = {8};
//+
Line Loop(10) = {2, -7, -6, -5};
//+
Plane Surface(11) = {10};
//+
Physical Line("deposition", 1) = {1, 5};
//+
Physical Line("slipx", 2) = {4};
//+
Physical Line("noslip", 3) = {6};
//+
Physical Line("noflux", 4) = {7, 3};
//+
Physical Surface("sediment", 5) = {9};
//+
Physical Surface("water", 6) = {11};
