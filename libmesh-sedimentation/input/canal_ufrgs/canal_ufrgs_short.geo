// Gmsh project created on Fri Sep 21 09:05:19 2012

Mesh.RandomFactor = 1e-6;

//tanque 5K
//scal  =  4.8;
//Ref1  =  0.005*scal;
//Ref2  =  0.04*scal;

scal  =  4.8;
Ref1  =  0.005*scal;
Ref2  =  0.04*scal;
Ref3  =  0.04*scal;
lc    =   0.1;


xc    =  13.5/lc;
yc    =   0.4/lc;
zc    =   0.5/lc;
hc    =  0.04/lc;

xt    =   1.0/lc;
yt    =   0.4/lc;
zt    =   1.5/lc;

Ref1 = 1.0*hc;
Ref2 = 1.0*hc;
Ref3 = 1.0*hc;

//channel
Point(1) = {       0,     yc/2,     0,     Ref1};
Point(2) = {       0,    -yc/2,     0,     Ref1};
Point(3) = {      xc,     yc/2,     0,     Ref1};
Point(4) = {      xc,    -yc/2,     0,     Ref1};
Point(5) = {       0,     yc/2,    zc,     Ref2};
Point(6) = {       0,    -yc/2,    zc,     Ref2};
Point(7) = {      xc,     yc/2,    zc,     Ref3};
Point(8) = {      xc,    -yc/2,    zc,     Ref3};

//inflow
Point(9)  = {     0,     yc/2,     hc,     Ref1};
Point(10) = {     0,    -yc/2,     hc,     Ref1};
Point(11) = {     xc,     yc/2,    hc,     Ref1};
Point(12) = {     xc,    -yc/2,    hc,     Ref1};


//+
Line(1) = {1, 2};
//+
Line(2) = {2, 10};
//+
Line(3) = {10, 9};
//+
Line(4) = {9, 1};
//+
Line(5) = {10, 6};
//+
Line(6) = {6, 5};
//+
Line(7) = {5, 9};
//+
Line(8) = {6, 8};
//+
Line(9) = {8, 7};
//+
Line(10) = {7, 5};
//+
Line(11) = {10, 12};
//+
Line(12) = {12, 11};
//+
Line(13) = {11, 9};
//+
Line(14) = {2, 4};
//+
Line(15) = {4, 3};
//+
Line(16) = {3, 1};
//+
Line(17) = {4, 12};
//+
Line(18) = {12, 8};
//+
Line(19) = {3, 11};
//+
Line(20) = {11, 7};


Mesh.Binary = 0;
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {3, -7, -6, -5};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {19, -12, -17, 15};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {20, -9, -18, 12};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {4, -16, 19, 13};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {13, -7, -10, -20};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {1, 14, 15, 16};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {3, -13, -12, -11};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {6, -10, -9, -8};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {5, 8, -18, -11};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {2, 11, -17, -14};
//+
Plane Surface(11) = {11};
//+
Surface Loop(1) = {2, 6, 9, 4, 10, 8};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {11, 1, 7, 3, 5, 8};
//+
Volume(2) = {2};

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

//+
Physical Surface("NOSLIP") = {11, 5, 2, 10, 6, 7};
//+
Physical Surface("INLET") = {1};
//+
Physical Surface("OUTFLOW") = {4, 3};
//+
Physical Surface("PNULL") = {9};
//+
Physical Surface("DEPOSITION") = {7};
//+
Physical Volume("FLUID") = {1, 2};
