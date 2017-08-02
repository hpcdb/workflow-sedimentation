// Gmsh project created on Fri Oct 28 15:48:31 2016
xmin  = 0.0;
xmax  = 18.0;
ymin  = 0.0;
ymax  = 2.0;
zmin  = 0.0;
zmax  = 2.0;
xlock = 0.75;
hlock = 2.0;

lc    = 0.75;

Point(1) = {xmin,  ymin, zmin, lc};
Point(2) = {xlock, ymin, zmin, lc};
Point(3) = {xlock, ymin, zmax, lc};
Point(4) = {xmin , ymin, zmax, lc};
Point(5) = {xmax,  ymin, zmin, lc};
Point(6) = {xmax , ymin, zmax, lc};
Point(7) = {xmin,  ymax, zmin, lc};
Point(8) = {xlock, ymax, zmin, lc};
Point(9) = {xlock, ymax, zmax, lc};
Point(10) = {xmin , ymax, zmax, lc};
Point(11) = {xmax,  ymax, zmin, lc};
Point(12) = {xmax , ymax, zmax, lc};



//Transfinite Volume {44};
//Transfinite Volume {46};

//Physical Surface("LEFT", 1)  = {38};
//Physical Surface("RIGHT", 2) = {42};
//Physical Surface("BOTTOM", 3) = {22,24};
//Physical Surface("TOP", 4) = {26,28};
//Physical Surface("BACK", 5) = {30,32};
//Physical Surface("FRONT", 6) = {34,36};
//Physical Volume("WATER", 7) = {44};
//Physical Volume("SEDIMENT", 8) = {46};


// LEFT
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 10};
//+
Line(3) = {10, 7};
//+
Line(4) = {7, 1};

// Lock
//+
Line(5) = {2, 3};
//+
Line(6) = {3, 9};
//+
Line(7) = {9, 8};
//+
Line(8) = {8, 2};

// RIGHT
//+
Line(9) = {11, 5};
//+
Line(10) = {5, 6};
//+
Line(11) = {6, 12};
//+
Line(12) = {12, 11};

// BOTTOM
//+
Line(13) = {1, 2};
//+
Line(14) = {8, 7};
//+
Line(15) = {2, 5};
//+
Line(16) = {11, 8};

// TOP
//+
Line(17) = {4, 3};
//+
Line(18) = {9, 10};
//+
Line(19) = {3, 6};
//+
Line(20) = {12, 9};

//SURFACES
// LEFT
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};

// LOCK
//+
Line Loop(2) = {8, 5, 6, 7};
//+
Plane Surface(2) = {2};

// RIGHT
//+
Line Loop(3) = {12, 9, 10, 11};
//+
Plane Surface(3) = {3};

// BOTTOM

//+
Line Loop(4) = {4, 13, -8, 14};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {8, 15, -9, 16};
//+
Plane Surface(5) = {5};

// TOP

//+
Line Loop(6) = {2, -18, -6, -17};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {6, -20, -11, -19};
//+
Plane Surface(7) = {7};

// BACK

//+
Line Loop(8) = {13, 5, -17, -1};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {5, 19, -10, -15};
//+
Plane Surface(9) = {9};

// FRONT

//+
Line Loop(10) = {3, -14, -7, 18};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {7, -16, -12, 20};
//+
Plane Surface(11) = {11};

// VOLUME
// SEDIMENT

//+
Surface Loop(1) = {1, 4, 8, 6, 10, 2};
//+
Volume(1) = {1};

// WATER
//+
Surface Loop(2) = {9, 7, 11, 5, 3, 2};
//+
Volume(2) = {2};

//+

Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Transfinite Surface {5};
Transfinite Surface {6};
Transfinite Surface {7};
Transfinite Surface {8};
Transfinite Surface {9};
Transfinite Surface {10};
Transfinite Surface {11};
Transfinite Volume {1};
Transfinite Volume {2};

Physical Surface("LEFT", 1)    = {1};
Physical Surface("RIGHT", 2)   = {3};
Physical Surface("BOTTOM", 3)  = {4,5};
Physical Surface("TOP", 4)     = {6,7};
Physical Surface("BACK", 5)    = {8,9};
Physical Surface("FRONT", 6)   = {10,11};
Physical Volume("WATER", 7)    = {2};
Physical Volume("SEDIMENT", 8) = {1};

Mesh.RecombineAll = 1;
Mesh.Algorithm    = 8;
