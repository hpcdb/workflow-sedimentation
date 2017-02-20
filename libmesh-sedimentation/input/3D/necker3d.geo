// Gmsh project created on Fri Oct 28 15:48:31 2016
xmin  = 0.0;
xmax  = 18.0;
ymin  = 0.0;
ymax  = 2.0;
zmin  = 0.0;
zmax  = 2.0;
xlock = 1.0;
hlock = 1.8;

lc    = 0.075;

Point(1) = {xmin,  ymin, zmin, lc};
Point(2) = {xlock, ymin, zmin, lc};
Point(3) = {xlock, ymin, zmax, lc};
Point(4) = {xmin , ymin, zmax, lc};
Point(5) = {xmax,  ymin, zmin, lc};
Point(6) = {xmax , ymin, zmax, lc};
Point(7) = {xmin,  ymin, hlock, lc};
Point(8) = {xlock, ymin, hlock, lc};
Point(9) = {xmax , ymin, hlock, lc};


Point(10) = {xmin,  ymax, zmin, lc};
Point(11) = {xlock, ymax, zmin, lc};
Point(12) = {xlock, ymax, zmax, lc};
Point(13) = {xmin , ymax, zmax, lc};
Point(14) = {xmax,  ymax, zmin, lc};
Point(15) = {xmax , ymax, zmax, lc};
Point(16) = {xmin,  ymax, hlock, lc};
Point(17) = {xlock, ymax, hlock, lc};
Point(18) = {xmax , ymax, hlock, lc};


Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 14};
Line(4) = {14, 11};
Line(5) = {11, 10};
Line(6) = {10, 1};
Line(7) = {11, 2};
Line(8) = {4, 3};
Line(9) = {3, 6};
Line(10) = {6, 15};
Line(11) = {15, 12};
Line(12) = {12, 13};
Line(13) = {13, 4};
Line(14) = {12, 3};
Line(15) = {7, 8};
Line(16) = {8, 9};
Line(17) = {9, 18};
Line(18) = {18, 17};
Line(19) = {17, 16};
Line(20) = {16, 7};
Line(21) = {17, 8};
Line(22) = {3, 8};
Line(23) = {4, 7};
Line(24) = {7, 1};
Line(25) = {8, 2};
Line(26) = {12, 17};
Line(27) = {17, 11};
Line(28) = {13, 16};
Line(29) = {16, 10};
Line(30) = {6, 9};
Line(31) = {9, 5};
Line(32) = {15, 18};
Line(33) = {18, 14};

Line Loop(34) = {1, -7, 5, 6};
Plane Surface(35) = {34};
Line Loop(36) = {7, 2, 3, 4};
Plane Surface(37) = {36};
Line Loop(38) = {24, 1, -25, -15};
Plane Surface(39) = {38};
Line Loop(40) = {24, -6, -29, 20};
Plane Surface(41) = {40};
Line Loop(42) = {5, -29, -19, 27};
Plane Surface(43) = {42};
Line Loop(44) = {20, 15, -21, 19};
Plane Surface(45) = {44};
Line Loop(46) = {25, -7, -27, 21};
Plane Surface(47) = {46};
Line Loop(48) = {25, 2, -31, -16};
Plane Surface(49) = {48};
Line Loop(50) = {21, 16, 17, 18};
Plane Surface(51) = {50};
Line Loop(52) = {22, 16, -30, -9};
Plane Surface(53) = {52};
Line Loop(54) = {15, -22, -8, 23};
Plane Surface(55) = {54};
Line Loop(56) = {20, -23, -13, 28};
Plane Surface(57) = {56};
Line Loop(58) = {19, -28, -12, 26};
Plane Surface(59) = {58};
Line Loop(60) = {21, -22, -14, 26};
Plane Surface(61) = {60};
Line Loop(62) = {18, -26, -11, 32};
Plane Surface(63) = {62};
Line Loop(64) = {17, -32, -10, 30};
Plane Surface(65) = {64};
Line Loop(66) = {14, 9, 10, 11};
Plane Surface(67) = {66};
Line Loop(68) = {3, -33, -17, 31};
Plane Surface(69) = {68};
Line Loop(70) = {4, -27, -18, 33};
Plane Surface(71) = {70};
Line Loop(72) = {13, 8, -14, 12};
Plane Surface(73) = {72};

Transfinite Surface {35};
Transfinite Surface {37};
Transfinite Surface {39};
Transfinite Surface {41};
Transfinite Surface {43};
Transfinite Surface {45};
Transfinite Surface {47};
Transfinite Surface {49};
Transfinite Surface {51};
Transfinite Surface {53};
Transfinite Surface {55};
Transfinite Surface {57};
Transfinite Surface {59};
Transfinite Surface {61};
Transfinite Surface {63};
Transfinite Surface {65};
Transfinite Surface {67};
Transfinite Surface {69};
Transfinite Surface {71};
Transfinite Surface {73};


Surface Loop(74) = {49, 37, 69, 71, 51, 47};
Volume(75) = {74};
Surface Loop(76) = {67, 53, 65, 63, 51, 61};
Volume(77) = {76};
Surface Loop(78) = {41, 39, 35, 43, 47, 45};
Volume(79) = {78};
Surface Loop(80) = {59, 57, 55, 73, 61, 45};
Volume(81) = {80};

Transfinite Volume {75};
Transfinite Volume {77};
Transfinite Volume {79};
Transfinite Volume {81};


//Physical Surface("deposition", 1) = {35, 37};
//Physical Surface("slipx", 2)      = {41, 57};
//Physical Surface("slipz", 3)      = {73, 67};
//Physical Surface("noslip", 4)     = {65,69};
//Physical Surface("slipy", 5)      = {43, 59, 63, 71, 39, 55, 49, 53};
//Physical Surface("pnull", 6)      = {73, 67};
//Physical Volume("sediment", 7)    = {79};
//Physical Volume("water", 8)       = {75, 77, 81};



//+

Physical Surface("LEFT", 1) = {41, 57};
//+
Physical Surface("RIGHT", 2) = {69, 65};
//+
Physical Surface("BOTTOM", 3) = {35, 37};
//+
Physical Surface("TOP", 4) = {73, 67};
//+
Physical Surface("FRONT", 5) = {49, 39, 53, 55};
//+
Physical Surface("BACK", 6) = {71, 63, 43, 59};
//+
Physical Volume("WATER", 7) = {75, 77, 81};
//+
Physical Volume("SEDIMENT", 8) = {79};

Mesh.RecombineAll = 1;
Mesh.Algorithm3D  = 6;
