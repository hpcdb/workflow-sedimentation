Mesh.RandomFactor = 1e-6;

lc    =   0.1;

xc    =  13.5/lc;
yc    =   0.5/lc;
hc    =  0.04/lc;

xt    =   1.0/lc;
yt    =   1.5/lc;

Ref1 = 1.0*hc;
Ref2 = 1.0*hc;
Ref3 = 1.0*hc;

//channel
Point(1) = {       0,     0,    0,   Ref1};
Point(2) = {      xc,     0,    0,   Ref1};
Point(3) = {       0,    yc,    0,   Ref2};
Point(4) = {      xc,    yc,    0,   Ref2};

//inflow
Point(5)  = {      0,    hc,    0,   Ref1};
Point(6)  = {     xc,    hc,    0,   Ref1};
Point(7)  = {  xc+xt,    hc,    0,   Ref1};
Point(8)  = {  xc+xt,     0,    0,   Ref1};

Point(11) = {  xc+xt,    yc,    0,   Ref3};
Point(13) = {    xc , yc-yt,    0,   Ref1};
Point(15) = {  xc+xt, yc-yt,    0,   Ref2};

//+
Line(1) = {1, 5};
//+
Line(2) = {5, 3};
//+
Line(3) = {1, 2};
//+
Line(4) = {5, 6};
//+
Line(5) = {3, 4};
//+
Line(6) = {4, 6};
//+
Line(7) = {6, 2};
//+
Line(8) = {4, 11};
//+
Line(9) = {11, 7};
//+
Line(10) = {7, 6};
//+
Line(11) = {2, 8};
//+
Line(12) = {8, 7};
//+
Line(13) = {2, 13};
//+
Line(14) = {13, 15};
//+
Line(15) = {15, 8};
//+
Line Loop(1) = {2, 5, 6, -4};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {4, 7, -3, 1};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {6, -10, -9, -8};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {10, 7, 11, 12};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {13, 14, 15, -11};
//+
Plane Surface(5) = {5};


Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Transfinite Surface {5};


//+
Physical Line("INLET",1) = {1};
//+
Physical Line("NOSLIP",2) = {2, 9, 13, 15, 12};
//+
Physical Line("DEPOSITION",3) = {3};
//+
Physical Line("SLIPY",4) = {5, 8};
//+
Physical Line("OUTFLOW",5) = {14};
//+
Physical Surface("FLUID",6) = {1, 2, 3, 4, 5};

Mesh.Binary = 0;
Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
