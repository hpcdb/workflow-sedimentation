Mesh.RandomFactor = 1e-6;

xc    =  10.;
yc    =   .5;
hc    =  0.008;


Ref = 10.0*hc;
Refc = 1.0*hc;

//channel
Point(1) = {       0,     0,    0,   Ref};
Point(2) = {      xc,     0,    0,   Ref};
Point(3) = {       0,    yc,    0,   Ref};
Point(4) = {      xc,    yc,    0,   Ref};

//inflow
Point(5)  = {      0,    hc,    0,   Ref};
Point(6)  = {     xc,    hc,    0,   Ref};

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
Line Loop(1) = {2, 5, 6, -4};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {4, 7, -3, 1};
//+
Plane Surface(2) = {2};


Transfinite Surface {1};
Transfinite Surface {2};


Physical Line("INLET",1) = {1};
//+
Physical Line("NOSLIP",2) = {2};
//+
Physical Line("DEPOSITION",3) = {3};
//+
Physical Line("SLIPY",4) = {5};
//+
Physical Line("OUTFLOW",5) = {6, 7};
//+
Physical Surface("FLUID",6) = {1, 2};

Mesh.Binary = 0;
Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
