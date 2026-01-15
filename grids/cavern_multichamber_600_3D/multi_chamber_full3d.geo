//////////////////////////////////////////////////////////////
// MULTI-CHAMBER CAVERN - FULL 3D (No Symmetry)
// Three connected spherical chambers
// SCALED to 600k m³ volume
// Centered in 3D domain
//////////////////////////////////////////////////////////////

size_coarse = 65;
size_fine   = 4.5;

Lz = 660;
Ly = 450;
Lx = 450;

// Geometry parameters (SCALED)
h_bottom = 154.702263;              // was 205.189718, scaled by 0.753947
H1 = 52.776321;                     // was 70.0, scaled by 0.753947 (chamber 1-2 spacing)
H2 = 52.776321;                     // was 70.0, scaled by 0.753947 (chamber 2-3 spacing)
R1 = 36.189477;                     // was 48.0, scaled by 0.753947 (bottom chamber)
R2 = 31.665793;                     // was 42.0, scaled by 0.753947 (middle chamber)
R3 = 37.697372;                     // was 50.0, scaled by 0.753947 (top chamber)

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
// CAVERN GEOMETRY - BOTTOM CHAMBER
//////////////////////////////

Point(100) = {x_center, y_center, h_bottom, size_fine};                   // bottom tip
Point(101) = {x_center, y_center, h_bottom+R1, size_coarse};             // chamber 1 center

Point(102) = {x_center+R1, y_center, h_bottom+R1, size_fine};            // +x
Point(103) = {x_center-R1, y_center, h_bottom+R1, size_fine};            // -x
Point(104) = {x_center, y_center+R1, h_bottom+R1, size_fine};            // +y
Point(105) = {x_center, y_center-R1, h_bottom+R1, size_fine};            // -y

//////////////////////////////
// CAVERN GEOMETRY - MIDDLE CHAMBER
//////////////////////////////

Point(110) = {x_center, y_center, h_bottom+R1+H1, size_coarse};          // chamber 2 center

Point(112) = {x_center+R2, y_center, h_bottom+R1+H1, size_fine};         // +x
Point(113) = {x_center-R2, y_center, h_bottom+R1+H1, size_fine};         // -x
Point(114) = {x_center, y_center+R2, h_bottom+R1+H1, size_fine};         // +y
Point(115) = {x_center, y_center-R2, h_bottom+R1+H1, size_fine};         // -y

//////////////////////////////
// CAVERN GEOMETRY - TOP CHAMBER
//////////////////////////////

Point(120) = {x_center, y_center, h_bottom+R1+H1+H2+R3, size_fine};      // top tip
Point(121) = {x_center, y_center, h_bottom+R1+H1+H2, size_coarse};       // chamber 3 center

Point(122) = {x_center+R3, y_center, h_bottom+R1+H1+H2, size_fine};      // +x
Point(123) = {x_center-R3, y_center, h_bottom+R1+H1+H2, size_fine};      // -x
Point(124) = {x_center, y_center+R3, h_bottom+R1+H1+H2, size_fine};      // +y
Point(125) = {x_center, y_center-R3, h_bottom+R1+H1+H2, size_fine};      // -y

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

// Chamber 1 to Chamber 2
Line(20) = {102, 112};  // +x
Line(21) = {103, 113};  // -x
Line(22) = {104, 114};  // +y
Line(23) = {105, 115};  // -y

// Chamber 2 to Chamber 3
Line(24) = {112, 122};  // +x
Line(25) = {113, 123};  // -x
Line(26) = {114, 124};  // +y
Line(27) = {115, 125};  // -y

//////////////////////////////
// CIRCLE ARCS - BOTTOM CHAMBER
//////////////////////////////

Circle(30) = {100, 101, 102};
Circle(31) = {100, 101, 104};
Circle(32) = {100, 101, 103};
Circle(33) = {100, 101, 105};

Circle(34) = {102, 101, 104};
Circle(35) = {104, 101, 103};
Circle(36) = {103, 101, 105};
Circle(37) = {105, 101, 102};

//////////////////////////////
// CIRCLE ARCS - MIDDLE CHAMBER
//////////////////////////////

Circle(40) = {112, 110, 114};
Circle(41) = {114, 110, 113};
Circle(42) = {113, 110, 115};
Circle(43) = {115, 110, 112};

//////////////////////////////
// CIRCLE ARCS - TOP CHAMBER
//////////////////////////////

Circle(50) = {120, 121, 122};
Circle(51) = {120, 121, 124};
Circle(52) = {120, 121, 123};
Circle(53) = {120, 121, 125};

Circle(54) = {122, 121, 124};
Circle(55) = {124, 121, 123};
Circle(56) = {123, 121, 125};
Circle(57) = {125, 121, 122};

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

// Transition surfaces chamber 1 to 2
Curve Loop(105) = {34, 22, -40, -20};
Surface(105) = {105};

Curve Loop(106) = {35, 21, -41, -22};
Surface(106) = {106};

Curve Loop(107) = {36, 23, -42, -21};
Surface(107) = {107};

Curve Loop(108) = {37, 20, -43, -23};
Surface(108) = {108};

// Transition surfaces chamber 2 to 3
Curve Loop(109) = {40, 26, -54, -24};
Surface(109) = {109};

Curve Loop(110) = {41, 25, -55, -26};
Surface(110) = {110};

Curve Loop(111) = {42, 27, -56, -25};
Surface(111) = {111};

Curve Loop(112) = {43, 24, -57, -27};
Surface(112) = {112};

// Top chamber surfaces
Curve Loop(113) = {50, 54, -51};
Surface(113) = {113};

Curve Loop(114) = {51, 55, -52};
Surface(114) = {114};

Curve Loop(115) = {52, 56, -53};
Surface(115) = {115};

Curve Loop(116) = {53, 57, -50};
Surface(116) = {116};

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
Surface Loop(200) = {101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};

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

Physical Surface("Cavern", 29) = {101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116};
Physical Volume("Salt",    28) = {1};

Physical Curve("Wall_profile", 30) = {50, 24, 20, 30};
