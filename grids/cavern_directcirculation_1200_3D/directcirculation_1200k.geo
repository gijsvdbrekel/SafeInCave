// ============================================================================
// IRREGULAR CAVERN - SCALED TO 1.2 MILLION m³
// ============================================================================
// Target volume:   1.2 Mm³
// Scale factor:    1.257457
// 
// All coordinates scaled (z-axis and radial distances from center)
// Domain dimensions (Lx, Ly, Lz) remain 450×450×660 m
// ============================================================================

//////////////////////////////////////////////////////////////
// IRREGULAR CAVERN - FULL 3D COMPLETE (No Symmetry)
// SCALED to 620k m³ volume (scale factor 0.940774)
// Centered in 3D domain
// FIXED: Both tips use Lines instead of Circles
//////////////////////////////////////////////////////////////

coarse_size = 65;
fine_size = 4.5;

Lz = 660.0;
Ly = 450.0;
Lx = 450.0;

x_center = Lx/2;  // 225.0
y_center = Ly/2;  // 225.0

//////////////////////////////
// OUTER BOX
//////////////////////////////

Point(1) = {0, 0, 0, coarse_size};
Point(2) = {Lx, 0, 0, coarse_size};
Point(3) = {Lx, Ly, 0, coarse_size};
Point(4) = {0, Ly, 0, coarse_size};

Point(5) = {0, 0, Lz, coarse_size};
Point(6) = {Lx, 0, Lz, coarse_size};
Point(7) = {Lx, Ly, Lz, coarse_size};
Point(8) = {0, Ly, Lz, coarse_size};

//////////////////////////////
// CAVERN GEOMETRY - 6 LEVELS
//////////////////////////////

// Level 1 - Bottom (z=180, R=16.55)
Point(100) = {225.000000, 225.000000, 212.936803, fine_size};           // tip  // scaled
Point(101) = {225.000000, 225.000000, 218.851714, coarse_size};         // center  // scaled
Point(102) = {244.578357, 225.000000, 218.851714, fine_size};   // +x  // scaled
Point(103) = {205.421643, 225.000000, 218.851714, fine_size};   // -x  // scaled
Point(104) = {225.000000, 244.578357, 218.851714, fine_size};   // +y  // scaled
Point(105) = {225.000000, 205.421643, 218.851714, fine_size};   // -y  // scaled

// Level 2 - Lower bulge (z=205, R=51.72)
Point(110) = {225.000000, 225.000000, 242.511359, coarse_size};  // scaled
Point(112) = {286.183841, 225.000000, 242.511359, fine_size};  // scaled
Point(113) = {163.816159, 225.000000, 242.511359, fine_size};  // scaled
Point(114) = {225.000000, 286.183841, 242.511359, fine_size};  // scaled
Point(115) = {225.000000, 163.816159, 242.511359, fine_size};  // scaled

// Level 3 - Maximum (z=230, R=56.0)
Point(120) = {225.000000, 225.000000, 272.085915, coarse_size};  // scaled
Point(122) = {291.247006, 225.000000, 272.085915, fine_size};  // scaled
Point(123) = {158.752994, 225.000000, 272.085915, fine_size};  // scaled
Point(124) = {225.000000, 291.247006, 272.085915, fine_size};  // scaled
Point(125) = {225.000000, 158.752994, 272.085915, fine_size};  // scaled

// Level 4 - Constriction (z=270, R=36.0)
Point(130) = {225.000000, 225.000000, 319.405204, coarse_size};  // scaled
Point(132) = {267.587360, 225.000000, 319.405204, fine_size};  // scaled
Point(133) = {182.412640, 225.000000, 319.405204, fine_size};  // scaled
Point(134) = {225.000000, 267.587360, 319.405204, fine_size};  // scaled
Point(135) = {225.000000, 182.412640, 319.405204, fine_size};  // scaled

// Level 5 - Upper (z=310, R=26.0)
Point(140) = {225.000000, 225.000000, 366.724494, coarse_size};  // scaled
Point(142) = {255.757538, 225.000000, 366.724494, fine_size};  // scaled
Point(143) = {194.242462, 225.000000, 366.724494, fine_size};  // scaled
Point(144) = {225.000000, 255.757538, 366.724494, fine_size};  // scaled
Point(145) = {225.000000, 194.242462, 366.724494, fine_size};  // scaled

// Level 6 - Top (z=335, R=8.0)
Point(150) = {225.000000, 225.000000, 402.213961, fine_size};           // tip  // scaled
Point(151) = {225.000000, 225.000000, 396.299049, coarse_size};         // center  // scaled
Point(152) = {234.463858, 225.000000, 396.299049, fine_size};  // scaled
Point(153) = {215.536142, 225.000000, 396.299049, fine_size};  // scaled
Point(154) = {225.000000, 234.463858, 396.299049, fine_size};  // scaled
Point(155) = {225.000000, 215.536142, 396.299049, fine_size};  // scaled

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

// Vertical connections (+x)
Line(20) = {102, 112};
Line(21) = {112, 122};
Line(22) = {122, 132};
Line(23) = {132, 142};
Line(24) = {142, 152};

// Vertical connections (-x)
Line(30) = {103, 113};
Line(31) = {113, 123};
Line(32) = {123, 133};
Line(33) = {133, 143};
Line(34) = {143, 153};

// Vertical connections (+y)
Line(40) = {104, 114};
Line(41) = {114, 124};
Line(42) = {124, 134};
Line(43) = {134, 144};
Line(44) = {144, 154};

// Vertical connections (-y)
Line(50) = {105, 115};
Line(51) = {115, 125};
Line(52) = {125, 135};
Line(53) = {135, 145};
Line(54) = {145, 155};

//////////////////////////////
// CIRCLE ARCS AT EACH LEVEL
//////////////////////////////

// Level 1 - Bottom (LINES from tip to ring)
Line(60) = {100, 102};
Line(61) = {100, 104};
Line(62) = {100, 103};
Line(63) = {100, 105};

Circle(64) = {102, 101, 104};
Circle(65) = {104, 101, 103};
Circle(66) = {103, 101, 105};
Circle(67) = {105, 101, 102};

// Level 2
Circle(70) = {112, 110, 114};
Circle(71) = {114, 110, 113};
Circle(72) = {113, 110, 115};
Circle(73) = {115, 110, 112};

// Level 3
Circle(80) = {122, 120, 124};
Circle(81) = {124, 120, 123};
Circle(82) = {123, 120, 125};
Circle(83) = {125, 120, 122};

// Level 4
Circle(90) = {132, 130, 134};
Circle(91) = {134, 130, 133};
Circle(92) = {133, 130, 135};
Circle(93) = {135, 130, 132};

// Level 5
Circle(100) = {142, 140, 144};
Circle(101) = {144, 140, 143};
Circle(102) = {143, 140, 145};
Circle(103) = {145, 140, 142};

// Level 6 - Top (LINES from tip to ring - FIXED!)
Line(110) = {150, 152};
Line(111) = {150, 154};
Line(112) = {150, 153};
Line(113) = {150, 155};

Circle(114) = {152, 151, 154};
Circle(115) = {154, 151, 153};
Circle(116) = {153, 151, 155};
Circle(117) = {155, 151, 152};

//////////////////////////////
// CAVERN SURFACES
//////////////////////////////

// Bottom cap (4 surfaces)
Curve Loop(201) = {60, 64, -61};
Surface(201) = {201};

Curve Loop(202) = {61, 65, -62};
Surface(202) = {202};

Curve Loop(203) = {62, 66, -63};
Surface(203) = {203};

Curve Loop(204) = {63, 67, -60};
Surface(204) = {204};

// Side walls Level 1-2 (4 surfaces)
Curve Loop(211) = {64, 40, -70, -20};
Surface(211) = {211};

Curve Loop(212) = {65, 30, -71, -40};
Surface(212) = {212};

Curve Loop(213) = {66, 50, -72, -30};
Surface(213) = {213};

Curve Loop(214) = {67, 20, -73, -50};
Surface(214) = {214};

// Side walls Level 2-3 (4 surfaces)
Curve Loop(221) = {70, 41, -80, -21};
Surface(221) = {221};

Curve Loop(222) = {71, 31, -81, -41};
Surface(222) = {222};

Curve Loop(223) = {72, 51, -82, -31};
Surface(223) = {223};

Curve Loop(224) = {73, 21, -83, -51};
Surface(224) = {224};

// Side walls Level 3-4 (4 surfaces)
Curve Loop(231) = {80, 42, -90, -22};
Surface(231) = {231};

Curve Loop(232) = {81, 32, -91, -42};
Surface(232) = {232};

Curve Loop(233) = {82, 52, -92, -32};
Surface(233) = {233};

Curve Loop(234) = {83, 22, -93, -52};
Surface(234) = {234};

// Side walls Level 4-5 (4 surfaces)
Curve Loop(241) = {90, 43, -100, -23};
Surface(241) = {241};

Curve Loop(242) = {91, 33, -101, -43};
Surface(242) = {242};

Curve Loop(243) = {92, 53, -102, -33};
Surface(243) = {243};

Curve Loop(244) = {93, 23, -103, -53};
Surface(244) = {244};

// Side walls Level 5-6 (4 surfaces)
Curve Loop(251) = {100, 44, -114, -24};
Surface(251) = {251};

Curve Loop(252) = {101, 34, -115, -44};
Surface(252) = {252};

Curve Loop(253) = {102, 54, -116, -34};
Surface(253) = {253};

Curve Loop(254) = {103, 24, -117, -54};
Surface(254) = {254};

// Top cap (4 surfaces - FIXED!)
Curve Loop(261) = {110, 114, -111};
Surface(261) = {261};

Curve Loop(262) = {111, 115, -112};
Surface(262) = {262};

Curve Loop(263) = {112, 116, -113};
Surface(263) = {263};

Curve Loop(264) = {113, 117, -110};
Surface(264) = {264};

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

// Cavern surface loop (all 32 surfaces)
Surface Loop(200) = {201, 202, 203, 204,
                     211, 212, 213, 214,
                     221, 222, 223, 224,
                     231, 232, 233, 234,
                     241, 242, 243, 244,
                     251, 252, 253, 254,
                     261, 262, 263, 264};

// Outer box surface loop
Surface Loop(201) = {1, 2, 3, 4, 5, 6};

// Salt volume (box minus cavern)
Volume(1) = {201, 200};

//////////////////////////////
// PHYSICAL GROUPS
//////////////////////////////

Physical Surface("Bottom", 27) = {1};
Physical Surface("Top", 22) = {2};
Physical Surface("South", 23) = {3};
Physical Surface("North", 24) = {4};
Physical Surface("West", 26) = {5};
Physical Surface("East", 25) = {6};

Physical Surface("Cavern", 29) = {201, 202, 203, 204,
                                  211, 212, 213, 214,
                                  221, 222, 223, 224,
                                  231, 232, 233, 234,
                                  241, 242, 243, 244,
                                  251, 252, 253, 254,
                                  261, 262, 263, 264};

Physical Volume("Salt", 28) = {1};

Physical Curve("Wall_profile", 30) = {110, 24, 23, 22, 21, 20, 60};
