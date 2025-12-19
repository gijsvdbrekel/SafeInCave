//////////////////////////////////////////////////////////////
// ASYMMETRIC CAVERN - Elliptical chambers with axis rotation
// SCALED to 0.6M m³ volume
// Bottom chamber: elongated in x-direction
// Top chamber: elongated in y-direction
//////////////////////////////////////////////////////////////

size_coarse = 65;
size_fine   = 4.5;

Lz = 660;
Ly = 450;

// Geometry parameters (SCALED by 0.708107)
h  = 205.189718;       // bottom of cavern
H = 113.318300;  // was 160.000 (original), now scaled to 600k m³
Rx_bottom = 42.522636;  // was 60.000 (original), now scaled to 600k m³
Ry_bottom = 28.348424;  // was 40.000 (original), now scaled to 600k m³
Rx_top = 28.348424;  // was 40.000 (original), now scaled to 600k m³
Ry_top = 42.522636;  // was 60.000 (original), now scaled to 600k m³

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

// Bottom elliptical chamber (wide in x, narrow in y)
Point(5)  = {0,   0, h,             size_fine};        // deepest point
Point(6)  = {0,   0, h+Ry_bottom,   size_coarse};      // bottom center
Point(9)  = {Rx_bottom,  0, h+Ry_bottom,  size_fine};  // x direction (wide)
Point(12) = {0,   Ry_bottom, h+Ry_bottom, size_fine};  // y direction (narrow)

// Top elliptical chamber (narrow in x, wide in y)
Point(7)  = {0,   0, h+Ry_bottom+H,        size_coarse};   // top center
Point(8)  = {0,   0, h+Ry_bottom+H+Ry_top, size_fine};     // top point
Point(10) = {Rx_top,  0, h+Ry_bottom+H,    size_fine};     // x direction (narrow)
Point(11) = {0,   Ry_top, h+Ry_bottom+H,   size_fine};     // y direction (wide)

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

Line(14) = {9,10};      // x-side connection
Line(15) = {12,11};     // y-side connection

//////////////////////////////
// CIRCLE ARCS - BOTTOM (elliptical)
//////////////////////////////

Ellipse(16) = {5, 6, 9, 9};      // arc in xz plane
Ellipse(17) = {5, 6, 12, 12};    // arc in yz plane
Ellipse(18) = {12, 6, 9, 9};     // connecting arc

//////////////////////////////
// CIRCLE ARCS - TOP (elliptical)
//////////////////////////////

Ellipse(19) = {11, 7, 10, 10};   // arc in yz plane
Ellipse(20) = {10, 7, 8, 8};     // arc in xz plane
Ellipse(21) = {8, 7, 11, 11};    // connecting arc

//////////////////////////////
// CAVERN SURFACES
//////////////////////////////

Curve Loop(1) = {16, -18, -17};
Surface(1) = {1};     // bottom chamber

Curve Loop(2) = {19, 20, 21};
Surface(2) = {2};     // top chamber

Curve Loop(3) = {14, -19, -15, 18};
Surface(3) = {3};     // side wall between chambers

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
