// ============================================================================
// ASYMMETRIC CAVERN - SCALED TO 1.2 MILLION m³
// ============================================================================
// Original volume: ~0.6 Mm³
// Target volume:   1.2 Mm³
// Scale factor:    1.247493
// 
// All geometric parameters (R, H, h_bottom, etc.) have been scaled
// Domain dimensions (Lx, Ly, Lz) remain 450×450×660 m
// ============================================================================

//////////////////////////////////////////////////////////////
// ASYMMETRIC CAVERN - FULL 3D (No Symmetry)
// SCALED to 1.4M m³ volume
// Bottom chamber: elongated in x-direction
// Top chamber: elongated in y-direction
// Centered in 3D domain
//////////////////////////////////////////////////////////////

size_coarse = 65;
size_fine   = 4.5;

Lz = 660;
Ly = 450;
Lx = 450;

// Geometry parameters (SCALED)
h_bottom = 192.989990;  // was 154.702263, scaled by 1.247493 for 1.2M m³
H_cavern = 141.363786;  // was 113.318300, scaled by 1.247493 for 1.2M m³
Rx_bottom = 53.046691;  // was 42.522636, scaled by 1.247493 for 1.2M m³
Ry_bottom = 35.364461;  // was 28.348424, scaled by 1.247493 for 1.2M m³
Rx_top = 35.364461;  // was 28.348424, scaled by 1.247493 for 1.2M m³
Ry_top = 53.046691;  // was 42.522636, scaled by 1.247493 for 1.2M m³

// Cavern center positions
x_center = Lx/2;
y_center = Ly/2;

//////////////////////////////
// OUTER BOX
//////////////////////////////

Point(1) = {0,  0,  0,  size_coarse};
Point(2) = {Lx, 0,  0,  size_coarse};
Point(3) = {Lx, Ly, 0,  size_coarse};
Point(4) = {0,  Ly, 0,  size_coarse};

Point(5) = {0,  0,  Lz, size_coarse};
Point(6) = {Lx, 0,  Lz, size_coarse};
Point(7) = {Lx, Ly, Lz, size_coarse};
Point(8) = {0,  Ly, Lz, size_coarse};

//////////////////////////////
// CAVERN GEOMETRY - BOTTOM ELLIPTICAL CHAMBER
//////////////////////////////

Point(100) = {x_center, y_center, h_bottom, size_fine};                          // deepest point
Point(101) = {x_center, y_center, h_bottom+Ry_bottom, size_coarse};             // bottom center

Point(102) = {x_center+Rx_bottom, y_center, h_bottom+Ry_bottom, size_fine};    // +x (wide)
Point(103) = {x_center-Rx_bottom, y_center, h_bottom+Ry_bottom, size_fine};    // -x (wide)
Point(104) = {x_center, y_center+Ry_bottom, h_bottom+Ry_bottom, size_fine};    // +y (narrow)
Point(105) = {x_center, y_center-Ry_bottom, h_bottom+Ry_bottom, size_fine};    // -y (narrow)

//////////////////////////////
// CAVERN GEOMETRY - TOP ELLIPTICAL CHAMBER
//////////////////////////////

Point(110) = {x_center, y_center, h_bottom+Ry_bottom+H_cavern+Ry_top, size_fine};      // top tip
Point(111) = {x_center, y_center, h_bottom+Ry_bottom+H_cavern, size_coarse};           // top center

Point(112) = {x_center+Rx_top, y_center, h_bottom+Ry_bottom+H_cavern, size_fine};      // +x (narrow)
Point(113) = {x_center-Rx_top, y_center, h_bottom+Ry_bottom+H_cavern, size_fine};      // -x (narrow)
Point(114) = {x_center, y_center+Ry_top, h_bottom+Ry_bottom+H_cavern, size_fine};      // +y (wide)
Point(115) = {x_center, y_center-Ry_top, h_bottom+Ry_bottom+H_cavern, size_fine};      // -y (wide)

//////////////////////////////
// OUTER BOX EDGES
//////////////////////////////

Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,1};
Line(5)  = {5,6};
Line(6)  = {6,7};
Line(7)  = {7,8};
Line(8)  = {8,5};
Line(9)  = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

//////////////////////////////
// CAVERN CONNECTION LINES
//////////////////////////////

Line(20) = {102, 112};  // +x side
Line(21) = {103, 113};  // -x side
Line(22) = {104, 114};  // +y side
Line(23) = {105, 115};  // -y side

//////////////////////////////
// ELLIPSE ARCS - BOTTOM CHAMBER
//////////////////////////////

Ellipse(30) = {100, 101, 102, 102};  // bottom to +x
Ellipse(31) = {100, 101, 104, 104};  // bottom to +y
Ellipse(32) = {100, 101, 103, 103};  // bottom to -x
Ellipse(33) = {100, 101, 105, 105};  // bottom to -y

Ellipse(34) = {102, 101, 104, 104};  // +x to +y
Ellipse(35) = {104, 101, 103, 103};  // +y to -x
Ellipse(36) = {103, 101, 105, 105};  // -x to -y
Ellipse(37) = {105, 101, 102, 102};  // -y to +x

//////////////////////////////
// ELLIPSE ARCS - TOP CHAMBER
//////////////////////////////

Ellipse(40) = {110, 111, 112, 112};  // top to +x
Ellipse(41) = {110, 111, 114, 114};  // top to +y
Ellipse(42) = {110, 111, 113, 113};  // top to -x
Ellipse(43) = {110, 111, 115, 115};  // top to -y

Ellipse(44) = {112, 111, 114, 114};  // +x to +y
Ellipse(45) = {114, 111, 113, 113};  // +y to -x
Ellipse(46) = {113, 111, 115, 115};  // -x to -y
Ellipse(47) = {115, 111, 112, 112};  // -y to +x

//////////////////////////////
// CAVERN SURFACES
//////////////////////////////

// Bottom chamber surfaces
Curve Loop(101) = {30, 34, -31};
Surface(101) = {101};

Curve Loop(102) = {31, 35, -32};
Surface(102) = {102};

Curve Loop(103) = {32, 36, -33};
Surface(103) = {103};

Curve Loop(104) = {33, 37, -30};
Surface(104) = {104};

// Side wall surfaces (connecting bottom to top)
Curve Loop(105) = {34, 22, -44, -20};
Surface(105) = {105};

Curve Loop(106) = {35, 21, -45, -22};
Surface(106) = {106};

Curve Loop(107) = {36, 23, -46, -21};
Surface(107) = {107};

Curve Loop(108) = {37, 20, -47, -23};
Surface(108) = {108};

// Top chamber surfaces
Curve Loop(109) = {40, 44, -41};
Surface(109) = {109};

Curve Loop(110) = {41, 45, -42};
Surface(110) = {110};

Curve Loop(111) = {42, 46, -43};
Surface(111) = {111};

Curve Loop(112) = {43, 47, -40};
Surface(112) = {112};

//////////////////////////////
// OUTER BOX SURFACES
//////////////////////////////

// Bottom (z=0)
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Top (z=Lz)
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

// South (y=0)
Curve Loop(3) = {1, 10, -5, -9};
Plane Surface(3) = {3};

// North (y=Ly)
Curve Loop(4) = {3, 12, -7, -11};
Plane Surface(4) = {4};

// West (x=0)
Curve Loop(5) = {4, 9, -8, -12};
Plane Surface(5) = {5};

// East (x=Lx)
Curve Loop(6) = {2, 11, -6, -10};
Plane Surface(6) = {6};

//////////////////////////////
// VOLUME
//////////////////////////////

// Cavern surface loop
Surface Loop(200) = {101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112};

// Salt volume (box minus cavern)
Surface Loop(201) = {1, 2, 3, 4, 5, 6};
Volume(1) = {201, 200};

//////////////////////////////
// PHYSICAL GROUPS
//////////////////////////////

Physical Surface("Bottom", 27) = {1};
Physical Surface("Top",    22) = {2};
Physical Surface("South",  23) = {3};
Physical Surface("North",  24) = {4};
Physical Surface("West",   26) = {5};
Physical Surface("East",   25) = {6};

Physical Surface("Cavern", 29) = {101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112};
Physical Volume("Salt",    28) = {1};

Physical Curve("Wall_profile", 30) = {40, 20, 30};
