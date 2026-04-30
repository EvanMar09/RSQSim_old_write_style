//Gmsh project created on Thu Oct 23 20:11:41 2025


//A script to generate a fault step over with gmesh ####
//Evan Marschall ###

//fault will be in the x-z plane#

//characteristic lengths wanted in domain
lc = 2e3;
lc_fault = 200;


//Fault Params ##
length_fault = 30e3;
fault_width = 15e3;
overlap = 5e3;
step_over_width = 3e3;
R_nucl = 1.5e3;

//Domain params 

Xmax = 100e3;
Xmin = -Xmax;
Ymax = 60e3;
Ymin = -Ymax;
Zmin = -fault_width*2;


//Volume ##

Point(1) = {Xmin, Ymin, 0, lc};
Point(2) = {Xmin, Ymax, 0, lc};
Point(3) = {Xmax, Ymax, 0, lc};
Point(4) = {Xmax, Ymin, 0, lc};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


Curve Loop(5) = {1,2,3,4};
Plane Surface(1) = {5};

Extrude {0,0, Zmin} { Surface{1}; }


//The 2 faults ##

//fault uno ###
Point(100) = {-.5 * overlap, 0.5 * step_over_width, 0, lc_fault};
Point(101) = {-.5 * overlap + length_fault, 0.5 * step_over_width, 0, lc_fault};
Point(102) = {-.5 * overlap + length_fault, 0.5 * step_over_width, -1 * fault_width, lc_fault};
Point(103) = {-.5 * overlap, 0.5 * step_over_width, -1 * fault_width, lc_fault};

Line(100) = {100, 101};
Line{100} In Surface{1};
Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 100};

Line Loop(105) = {100,101,102,103};
Plane Surface(100) = {105};

//Fault ni 

Point(200) = {.5 * overlap, -0.5 * step_over_width, 0, lc_fault};
Point(201) = {.5 * overlap - length_fault, -0.5 * step_over_width, 0, lc_fault};
Point(202) = {.5 * overlap - length_fault, -0.5 * step_over_width, -1*fault_width, lc_fault};
Point(203) = {.5 * overlap, -0.5 * step_over_width, -1*fault_width, lc_fault};

Line(200) = {201, 200};
Line{200} In Surface{1};
Line(201) = {200, 203};
Line(202) = {203, 202};
Line(203) = {202, 201};

Line Loop(204) = {200,201,202,203};
Plane Surface(200) = {204};


//Surface(101) = {105};
//Surface(201) = {204};

Surface{100,200} In Volume{1};

// Managing coarsening away from the fault
// Attractor field returns the distance to the curve (actually, the
// distance to 100 equidistant points on the curve)
Field[1] = Distance;
Field[1].FacesList = {100,200};


// Matheval field returns "distance squared + lc/20"
Field[2] = MathEval;
Field[2].F = Sprintf("0.05*F1 +(F1/2e3)^2 + %g", lc_fault);

Field[3] = Threshold;
Field[3].IField = 1;
Field[3].LcMin = lc_fault;
Field[3].LcMax = lc;
Field[3].DistMin = 2*lc_fault;
Field[3].DistMax = 2*lc_fault+0.001;


Field[4] = Restrict;
Field[4].IField = 3;
Field[4].FacesList = {100,200};



Field[5] = Min;
Field[5].FieldsList = {2,3,4};

Background Field = 5;

Physical Surface(101) = {1};
Physical Surface(102) = {100};
Physical Surface(103) = {200};
Physical Surface(105) = {14,18,22,26,27};


Physical Volume(1) = {1};
Mesh.MshFileVersion = 2.2;






