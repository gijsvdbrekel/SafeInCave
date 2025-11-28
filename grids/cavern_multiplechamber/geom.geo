//////////////////////////////////////////////////////////////
// MULTI-CHAMBER CAVERN - Three connected spherical chambers
// CORRECTED: Only define arcs actually used in quarter symmetry
//////////////////////////////////////////////////////////////

size_coarse = 65;
size_fine   = 4.5;

Lz = 660;
Ly = 450;

// Geometry parameters
h  = 205.189718;       // bottom of cavern
H1 = 70.0;             // height to middle chamber
H2 = 70.0;             // height to top chamber
R1 = 48.0;             // bottom chamber radius
R2 = 42.0;             // middle chamber radius  
R3 = 50.0;             // top chamber radius

//////////////////////////////
// OUTER BOX
//////////////////////////////

Point(1)  = {0,   0, 0,   size_coarse};
Point(2)  = {Ly,  0, 0,   size_coarse};
Point(3)  = {Ly,  0, Lz,  size_coarse};
Point(4)  = {0,   0, Lz,  size_coarse};

Point(13) = {0,   Ly, 0,   size_coarse};
Point(14) = {0,   Ly, Lz,  size_coarse};
Point(15) = {Ly,  Ly, 0,   size_coarse};
Point(16) = {Ly,  Ly, Lz,  size_coarse};

//////////////////////////////
// CAVERN GEOMETRY
//////////////////////////////

// Bottom chamber
Point(5)  = {0,   0, h,        size_fine};        // deepest point
Point(6)  = {0,   0, h+R1,     size_coarse};      // chamber 1 center
Point(9)  = {R1,  0, h+R1,     size_fine};        // x direction
Point(12) = {0,   R1, h+R1,    size_fine};        // y direction

// Middle chamber - only first quadrant
Point(17) = {0,   0, h+R1+H1,      size_coarse};  // chamber 2 center
Point(18) = {R2,  0, h+R1+H1,      size_fine};    // x direction
Point(19) = {0,   R2, h+R1+H1,     size_fine};    // y direction

// Top chamber
Point(7)  = {0,   0, h+R1+H1+H2,       size_coarse};   // chamber 3 center
Point(8)  = {0,   0, h+R1+H1+H2+R3,    size_fine};     // top point
Point(10) = {R3,  0, h+R1+H1+H2,       size_fine};     // x direction
Point(11) = {0,   R3, h+R1+H1+H2,      size_fine};     // y direction

//////////////////////////////
// OUTER BOX LINES
//////////////////////////////

Line(1)  = {1,2};
Line(2)  = {2,15};
Line(3)  = {15,13};
Line(4)  = {13,1};
Line(5)  = {1,5};
Line(6)  = {8,4};
Line(7)  = {4,3};
Line(8)  = {3,16};
Line(9)  = {16,14};
Line(10) = {14,4};
Line(11) = {13,14};
Line(12) = {15,16};
Line(13) = {2,3};

//////////////////////////////
// CAVERN CONNECTION LINES
//////////////////////////////

Line(14) = {9,18};      // chamber 1 to 2 (x)
Line(15) = {12,19};     // chamber 1 to 2 (y)
Line(16) = {18,10};     // chamber 2 to 3 (x)
Line(17) = {19,11};     // chamber 2 to 3 (y)

//////////////////////////////
// CIRCLE ARCS
//////////////////////////////

// Chamber 1 (bottom)
Circle(18) = {5, 6, 9};
Circle(19) = {5, 6, 12};
Circle(20) = {12, 6, 9};

// Chamber 2 (middle) - only one arc in first quadrant
Circle(21) = {18, 17, 19};

// Chamber 3 (top)
Circle(22) = {11, 7, 10};
Circle(23) = {10, 7, 8};
Circle(24) = {8, 7, 11};

//////////////////////////////
// CAVERN SURFACES
//////////////////////////////

Curve Loop(1) = {18, -20, -19};
Surface(1) = {1};     // bottom chamber

Curve Loop(2) = {14, 21, -15, 20};
Surface(2) = {2};     // transition 1-2

Curve Loop(3) = {16, -22, -17, -21};
Surface(3) = {3};     // transition 2-3

Curve Loop(4) = {22, 23, 24};
Surface(4) = {4};     // top chamber

//////////////////////////////
// OUTER BOX SURFACES
//////////////////////////////

Curve Loop(5) = {5, 18, 14, 16, 23, 6, 7, -13, -1};
Plane Surface(5) = {-5};

Curve Loop(6) = {5, 19, 15, 17, -24, 6, -10, -11, 4};
Plane Surface(6) = {6};

Curve Loop(7) = {4, 1, 2, 3};
Plane Surface(7) = {-7};

Curve Loop(8) = {2, 12, -8, -13};
Plane Surface(8) = {8};

Curve Loop(9) = {12, 9, -11, -3};
Plane Surface(9) = {-9};

Curve Loop(10) = {9, 10, 7, 8};
Plane Surface(10) = {10};

//////////////////////////////
// VOLUME
//////////////////////////////

Surface Loop(1) = {10, 9, 7, 8, 5, 6, 1, 2, 3, 4};
Volume(1) = {1};

//////////////////////////////
// PHYSICAL GROUPS
//////////////////////////////

Physical Surface("Top",    22) = {10};
Physical Surface("South",  23) = {5};
Physical Surface("North",  24) = {9};
Physical Surface("East",   25) = {8};
Physical Surface("West",   26) = {6};
Physical Surface("Bottom", 27) = {7};

Physical Surface("Cavern", 29) = {1, 2, 3, 4};
Physical Volume("Salt",    28) = {1};

Physical Curve("Wall_profile", 30) = {23, 16, 14, 18};
