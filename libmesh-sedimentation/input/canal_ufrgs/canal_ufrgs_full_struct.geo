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

//dissipation tank
Point(13) = {  xc+xt,     yt/2,     zc,   Ref3};
Point(14) = {  xc+xt,    -yt/2,     zc,   Ref3};
Point(15) = {  xc+xt,     yt/2,  zc-zt,   Ref2};
Point(16) = {  xc+xt,    -yt/2,  zc-zt,   Ref2};

Point(17) = {  xc+xt,     yt/2,     0,   Ref3};
Point(18) = {  xc+xt,    -yt/2,     0,   Ref3};
Point(19) = {  xc+xt,     yt/2,     hc,   Ref2};
Point(20) = {  xc+xt,    -yt/2,     hc,   Ref2};

Point(21) = {     xc,     yt/2,  zc-zt,   Ref1};
Point(22) = {     xc,    -yt/2,  zc-zt,   Ref1};
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
//+
Line(21) = {22, 4};
//+
Line(22) = {22, 21};
//+
Line(23) = {21, 3};
//+
Line(24) = {22, 16};
//+
Line(25) = {16, 15};
//+
Line(26) = {15, 21};
//+
Line(27) = {16, 18};
//+
Line(28) = {15, 17};
//+
Line(29) = {18, 17};
//+
Line(30) = {18, 4};
//+
Line(31) = {17, 3};
//+
Line(32) = {18, 20};
//+
Line(33) = {20, 12};
//+
Line(34) = {20, 19};
//+
Line(35) = {19, 11};
//+
Line(36) = {19, 17};
//+
Line(37) = {20, 14};
//+
Line(38) = {14, 8};
//+
Line(39) = {14, 13};
//+
Line(40) = {13, 7};
//+
Line(41) = {19, 13};
//+
Line Loop(1) = {5, 6, 7, -3};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {1, 2, 3, 4};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {9, -20, -12, 18};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {12, -19, -15, 17};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {23, -15, -21, 22};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {28, -29, -27, 25};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {29, -36, -34, -32};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {41, -39, -37, 34};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {20, -40, -41, 35};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {35, -19, -31, -36};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {28, 31, -23, -26};
//+
Plane Surface(11) = {11};
//+
Line Loop(12) = {9, -40, -39, 38};
//+
Plane Surface(12) = {12};
//+
Line Loop(13) = {12, -35, -34, 33};
//+
Plane Surface(13) = {13};
//+
Line Loop(14) = {15, -31, -29, 30};
//+
Plane Surface(14) = {14};
//+
Line Loop(15) = {22, -26, -25, -24};
//+
Plane Surface(15) = {15};
//+
Line Loop(16) = {18, -38, -37, 33};
//+
Plane Surface(16) = {16};
//+
Line Loop(17) = {21, -30, -27, -24};
//+
Plane Surface(17) = {17};
//+
Line Loop(18) = {17, -33, -32, 30};
//+
Plane Surface(18) = {18};
//+
Line Loop(19) = {10, 7, -13, 20};
//+
Plane Surface(19) = {19};
//+
Line Loop(20) = {13, 4, -16, 19};
//+
Plane Surface(20) = {20};
//+
Line Loop(21) = {10, -6, 8, 9};
//+
Plane Surface(21) = {21};
//+
Line Loop(22) = {13, -3, 11, 12};
//+
Plane Surface(22) = {22};
//+
Line Loop(23) = {16, 1, 14, 15};
//+
Plane Surface(23) = {23};
//+
Line Loop(24) = {8, -18, -11, 5};
//+
Plane Surface(24) = {24};
//+
Line Loop(25) = {11, -17, -14, 2};
//+
Plane Surface(25) = {25};
//+
Surface Loop(1) = {1, 24, 21, 19, 22, 3};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {25, 23, 20, 2, 22, 4};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {16, 12, 9, 8, 3, 13};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {18, 7, 10, 13, 14, 4};
//+
Volume(4) = {4};
//+
Surface Loop(5) = {17, 5, 11, 6, 15, 14};
//+
Volume(5) = {5};
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
Transfinite Surface {12};
Transfinite Surface {13};
Transfinite Surface {14};
Transfinite Surface {15};
Transfinite Surface {16};
Transfinite Surface {17};
Transfinite Surface {18};
Transfinite Surface {19};
Transfinite Surface {20};
Transfinite Surface {21};
Transfinite Surface {22};
Transfinite Surface {23};
Transfinite Surface {24};
Transfinite Surface {25};

Transfinite Volume {1};
Transfinite Volume {2};
Transfinite Volume {3};
Transfinite Volume {4};
Transfinite Volume {5};


//INLET
Physical Surface("INLET" , 1) = {2};
//+
Physical Surface("NOSLIP", 2) = {19, 20, 23, 24, 1, 25, 9, 8, 16, 10, 7, 18, 5, 11, 17, 15, 6};
//+
Physical Surface("SLIPZ") = {21, 12};
//+
Physical Surface("PNULL") = {12};
//+
Physical Surface("DEPOSITION") = {23};
//+
Physical Volume("FLUID") = {1, 2, 5, 4, 3};


Mesh.Binary = 0;
