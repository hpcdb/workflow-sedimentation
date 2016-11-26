// Gmsh project created on Fri Oct 28 15:48:31 2016
xmin = 0.0;
xmax = 18.0;
ymin = 0.0;
ymax = 2.0;
zmin = 0.0;
zmax = 2.0;
xlock = 1.00;
hlock = 1.85;
lc = 0.25;

Point(1) = {xmin, ymin, zmin, lc};
Point(2) = {xlock,ymin, zmin, lc};
Point(3) = {xlock,ymin, zmax, lc};
Point(4) = {xmin ,ymin, zmax, lc};
Point(5) = {xmax, ymin, zmin, lc};
Point(6) = {xmax ,ymin, zmax, lc}; 

Point(7) = {xmin, ymax, zmin, lc};
Point(8) = {xlock,ymax, zmin, lc};
Point(9) = {xlock,ymax, zmax, lc};
Point(10) = {xmin ,ymax, zmax, lc};
Point(11) = {xmax, ymax, zmin, lc};
Point(12) = {xmax ,ymax, zmax, lc};  

Point(13) = {xmin ,ymin, hlock, lc}; 
Point(14) = {xmin ,ymax, hlock, lc};

Point(15) = {xlock ,ymin, hlock, lc}; 
Point(16) = {xlock ,ymax, hlock, lc};


Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 11};
Line(4) = {11, 8};
Line(5) = {8, 7};
Line(6) = {7, 1};
Line(7) = {8, 2};
Line(8) = {10, 4};
Line(9) = {4, 3};
Line(10) = {3, 6};
Line(11) = {6, 5};
Line(12) = {6, 12};
Line(13) = {12, 9};
Line(14) = {9, 10};
Line(15) = {9, 3};
Line(16) = {10, 14};
Line(17) = {4, 13};
Line(18) = {9, 16};
Line(19) = {3, 15};
Line(20) = {15, 2};
Line(21) = {16, 8};
Line(22) = {13, 1};
Line(23) = {14, 7};
Line(24) = {14, 13};
Line(25) = {13, 15};
Line(26) = {15, 16};
Line(27) = {16, 14};
Line(28) = {12, 11};
Line Loop(29) = {2, 3, 4, 7};
Plane Surface(30) = {29};
Line Loop(31) = {1, -7, 5, 6};
Plane Surface(32) = {31};
Line Loop(33) = {22, 1, -20, -25};
Plane Surface(34) = {33};
Line Loop(35) = {7, -20, 26, 21};
Plane Surface(36) = {35};
Line Loop(37) = {26, -18, 15, 19};
Plane Surface(38) = {37};
Line Loop(39) = {25, 26, 27, 24};
Plane Surface(40) = {39};
Line Loop(41) = {25, -19, -9, 17};
Plane Surface(42) = {41};
Line Loop(43) = {8, 9, -15, 14};
Plane Surface(44) = {43};
Line Loop(45) = {27, -16, -14, 18};
Plane Surface(46) = {45};
Line Loop(47) = {23, -5, -21, 27};
Plane Surface(48) = {47};
Line Loop(49) = {23, 6, -22, -24};
Plane Surface(50) = {49};
Line Loop(51) = {24, -17, -8, 16};
Plane Surface(52) = {51};
Line Loop(53) = {15, 10, 12, 13};
Plane Surface(54) = {53};
Line Loop(55) = {4, -21, -18, -13, 28};
Plane Surface(56) = {55};
Line Loop(57) = {20, 2, -11, -10, 19};
Plane Surface(58) = {57};
Line Loop(59) = {28, -3, -11, 12};
Plane Surface(60) = {59};
Surface Loop(61) = {34, 50, 48, 32, 40, 36};
Volume(62) = {61};
Surface Loop(63) = {58, 30, 60, 56, 54, 38, 36};
Volume(64) = {63};
Surface Loop(65) = {44, 52, 42, 46, 38, 40};
Volume(66) = {65};


Physical Surface("DEPOSITION") = {32, 30};
Physical Surface("SLIPX") = {50, 52};
Physical Surface("SLIPY") = {58, 34, 42, 48, 46, 56};
Physical Surface("NOSLIP") = {60};
Physical Surface("PNULL") = {54, 44};
Physical Volume("SEDIMENT") = {62};
Physical Volume("WATER") = {66, 64};

