//////////////////////////////////////////////////////////////
// TEARDROP CAVERN - Wide top tapering to narrow pointed bottom
// SCALED to 1.4M m³ volume
//////////////////////////////////////////////////////////////

size_coarse = 65;
size_fine   = 4.5;

Lz = 660;
Ly = 450;

// Geometry parameters (SCALED by 0.9816)
h  = 220.0;            // bottom of cavern (slightly higher)
H  = 176.7;            // total height (was 180.0)
Rt = 63.8;             // top bulb radius (wide) (was 65.0)
Rb = 7.9;              // bottom tip radius (very narrow) (was 8.0)

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

// Bottom tip (small)
Point(5)  = {0,   0, h,        size_fine};        // deepest point
Point(6)  = {0,   0, h+Rb,     size_coarse};      // bottom tip center
Point(9)  = {Rb,  0, h+Rb,     size_fine};        // x direction
Point(12) = {0,   Rb, h+Rb,    size_fine};        // y direction

// Top bulb (large)
Point(7)  = {0,   0, h+Rb+H,       size_coarse};   // top center
Point(8)  = {0,   0, h+Rb+H+Rt,    size_fine};     // top point
Point(10) = {Rt,  0, h+Rb+H,       size_fine};     // x direction
Point(11) = {0,   Rt, h+Rb+H,      size_fine};     // y direction

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

Line(14) = {9,10};      // bottom to top (long taper)
Line(15) = {12,11};     // bottom to top (y direction)

//////////////////////////////
// CIRCLE ARCS - BOTTOM TIP
//////////////////////////////

Circle(16) = {5, 6, 9};
Circle(17) = {5, 6, 12};
Circle(18) = {12, 6, 9};

//////////////////////////////
// CIRCLE ARCS - TOP BULB
//////////////////////////////

Circle(19) = {11, 7, 10};
Circle(20) = {10, 7, 8};
Circle(21) = {8, 7, 11};

//////////////////////////////
// CAVERN SURFACES
//////////////////////////////

Curve Loop(1) = {16, -18, -17};
Surface(1) = {1};     // bottom tip

Curve Loop(2) = {19, 20, 21};
Surface(2) = {2};     // top bulb

Curve Loop(3) = {14, -19, -15, 18};
Surface(3) = {3};     // tapered side wall

//////////////////////////////
// OUTER BOX SURFACES
//////////////////////////////

Curve Loop(4) = {5, 16, 14, 20, 6, 7, -13, -1};
Plane Surface(4) = {-4};

Curve Loop(5) = {5, 17, 15, -21, 6, -10, -11, 4};
Plane Surface(5) = {5};

Curve Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {-6};

Curve Loop(7) = {2, 12, -8, -13};
Plane Surface(7) = {7};

Curve Loop(8) = {12, 9, -11, -3};
Plane Surface(8) = {-8};

Curve Loop(9) = {9, 10, 7, 8};
Plane Surface(9) = {9};

//////////////////////////////
// VOLUME
//////////////////////////////

Surface Loop(31) = {9, 8, 6, 7, 4, 5, 1, 3, 2};
Volume(1) = {31};

//////////////////////////////
// PHYSICAL GROUPS
//////////////////////////////

Physical Surface("Top",    22) = {9};
Physical Surface("South",  23) = {4};
Physical Surface("North",  24) = {8};
Physical Surface("East",   25) = {7};
Physical Surface("West",   26) = {5};
Physical Surface("Bottom", 27) = {6};

Physical Surface("Cavern", 29) = {1, 2, 3};
Physical Volume("Salt",    28) = {1};

Physical Curve("Wall_profile", 30) = {20, 14, 16};
