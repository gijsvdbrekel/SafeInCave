//////////////////////////////////////////////////////////////
// CAVERN WITH INTERLAYER-SHAPED LEDGES - NO INTERLAYER VOLUMES
// Target volume: 600,000 m³ (scaled by 0.9577)
// Same cavern shape as cavern_interlayer_600_3D but single Salt volume
//
// This version is for comparison: same geometry but no material
// heterogeneity. Single homogeneous Salt volume.
//////////////////////////////////////////////////////////////

coarse_size = 65;
fine_size = 4.5;

Lz = 660.0;
Ly = 450.0;
Lx = 450.0;

x_center = Lx/2;  // 225.0
y_center = Ly/2;  // 225.0

//////////////////////////////////////////////////////////////
// CAVERN Z-POSITIONS (same as interlayer version)
//////////////////////////////////////////////////////////////

z_inter1_bot = 185.0;
z_inter1_top = 200.0;
z_inter2_bot = 290.0;
z_inter2_top = 305.0;

//////////////////////////////////////////////////////////////
// CAVERN GEOMETRY PARAMETERS (identical to interlayer version)
//////////////////////////////////////////////////////////////

R_bottom = 40.22;      // Radius of bottom bulge
R_ledge1 = 28.73;      // Radius at ledge 1
R_middle = 45.97;      // Radius of middle bulge
R_ledge2 = 26.82;      // Radius at ledge 2
R_top = 30.65;         // Radius of top bulge

z_bot_tip = 145.0;     // Bottom tip of cavern
z_bot_max = 170.0;     // Maximum bulge in bottom zone
z_mid_max = 240.0;     // Maximum bulge in middle zone
z_top_max = 320.0;     // Maximum bulge in top zone
z_top_tip = 345.0;     // Top tip of cavern

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
// CAVERN POINTS
//////////////////////////////

// Level 0 - Bottom tip (z = 145)
Point(100) = {x_center, y_center, z_bot_tip, fine_size};

// Level 1 - Bottom zone maximum bulge (z = 170, R = 40.22)
Point(101) = {x_center, y_center, z_bot_max, coarse_size};
Point(102) = {x_center + R_bottom, y_center, z_bot_max, fine_size};
Point(103) = {x_center - R_bottom, y_center, z_bot_max, fine_size};
Point(104) = {x_center, y_center + R_bottom, z_bot_max, fine_size};
Point(105) = {x_center, y_center - R_bottom, z_bot_max, fine_size};

// Level 2 - Bottom of ledge 1 (z = 185, R = 28.73)
Point(110) = {x_center, y_center, z_inter1_bot, coarse_size};
Point(112) = {x_center + R_ledge1, y_center, z_inter1_bot, fine_size};
Point(113) = {x_center - R_ledge1, y_center, z_inter1_bot, fine_size};
Point(114) = {x_center, y_center + R_ledge1, z_inter1_bot, fine_size};
Point(115) = {x_center, y_center - R_ledge1, z_inter1_bot, fine_size};

// Level 3 - Top of ledge 1 (z = 200, R = 28.73)
Point(120) = {x_center, y_center, z_inter1_top, coarse_size};
Point(122) = {x_center + R_ledge1, y_center, z_inter1_top, fine_size};
Point(123) = {x_center - R_ledge1, y_center, z_inter1_top, fine_size};
Point(124) = {x_center, y_center + R_ledge1, z_inter1_top, fine_size};
Point(125) = {x_center, y_center - R_ledge1, z_inter1_top, fine_size};

// Level 4 - Middle zone maximum bulge (z = 240, R = 45.97)
Point(130) = {x_center, y_center, z_mid_max, coarse_size};
Point(132) = {x_center + R_middle, y_center, z_mid_max, fine_size};
Point(133) = {x_center - R_middle, y_center, z_mid_max, fine_size};
Point(134) = {x_center, y_center + R_middle, z_mid_max, fine_size};
Point(135) = {x_center, y_center - R_middle, z_mid_max, fine_size};

// Level 5 - Bottom of ledge 2 (z = 290, R = 26.82)
Point(140) = {x_center, y_center, z_inter2_bot, coarse_size};
Point(142) = {x_center + R_ledge2, y_center, z_inter2_bot, fine_size};
Point(143) = {x_center - R_ledge2, y_center, z_inter2_bot, fine_size};
Point(144) = {x_center, y_center + R_ledge2, z_inter2_bot, fine_size};
Point(145) = {x_center, y_center - R_ledge2, z_inter2_bot, fine_size};

// Level 6 - Top of ledge 2 (z = 305, R = 26.82)
Point(150) = {x_center, y_center, z_inter2_top, coarse_size};
Point(152) = {x_center + R_ledge2, y_center, z_inter2_top, fine_size};
Point(153) = {x_center - R_ledge2, y_center, z_inter2_top, fine_size};
Point(154) = {x_center, y_center + R_ledge2, z_inter2_top, fine_size};
Point(155) = {x_center, y_center - R_ledge2, z_inter2_top, fine_size};

// Level 7 - Top zone maximum bulge (z = 320, R = 30.65)
Point(160) = {x_center, y_center, z_top_max, coarse_size};
Point(162) = {x_center + R_top, y_center, z_top_max, fine_size};
Point(163) = {x_center - R_top, y_center, z_top_max, fine_size};
Point(164) = {x_center, y_center + R_top, z_top_max, fine_size};
Point(165) = {x_center, y_center - R_top, z_top_max, fine_size};

// Level 8 - Top tip (z = 345)
Point(170) = {x_center, y_center, z_top_tip, fine_size};

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
// CAVERN VERTICAL CONNECTIONS
//////////////////////////////

// Bottom tip to level 1
Line(20) = {100, 102};
Line(21) = {100, 103};
Line(22) = {100, 104};
Line(23) = {100, 105};

// Level 1 to level 2
Line(30) = {102, 112};
Line(31) = {103, 113};
Line(32) = {104, 114};
Line(33) = {105, 115};

// Level 2 to level 3 - VERTICAL (ledge 1)
Line(40) = {112, 122};
Line(41) = {113, 123};
Line(42) = {114, 124};
Line(43) = {115, 125};

// Level 3 to level 4
Line(50) = {122, 132};
Line(51) = {123, 133};
Line(52) = {124, 134};
Line(53) = {125, 135};

// Level 4 to level 5
Line(60) = {132, 142};
Line(61) = {133, 143};
Line(62) = {134, 144};
Line(63) = {135, 145};

// Level 5 to level 6 - VERTICAL (ledge 2)
Line(70) = {142, 152};
Line(71) = {143, 153};
Line(72) = {144, 154};
Line(73) = {145, 155};

// Level 6 to level 7
Line(80) = {152, 162};
Line(81) = {153, 163};
Line(82) = {154, 164};
Line(83) = {155, 165};

// Level 7 to top tip
Line(90) = {162, 170};
Line(91) = {163, 170};
Line(92) = {164, 170};
Line(93) = {165, 170};

//////////////////////////////
// CIRCLE ARCS AT EACH LEVEL
//////////////////////////////

// Level 1 - Bottom bulge (R = 40.22)
Circle(100) = {102, 101, 104};
Circle(101) = {104, 101, 103};
Circle(102) = {103, 101, 105};
Circle(103) = {105, 101, 102};

// Level 2 - Bottom of ledge 1 (R = 28.73)
Circle(110) = {112, 110, 114};
Circle(111) = {114, 110, 113};
Circle(112) = {113, 110, 115};
Circle(113) = {115, 110, 112};

// Level 3 - Top of ledge 1 (R = 28.73)
Circle(120) = {122, 120, 124};
Circle(121) = {124, 120, 123};
Circle(122) = {123, 120, 125};
Circle(123) = {125, 120, 122};

// Level 4 - Middle bulge (R = 45.97)
Circle(130) = {132, 130, 134};
Circle(131) = {134, 130, 133};
Circle(132) = {133, 130, 135};
Circle(133) = {135, 130, 132};

// Level 5 - Bottom of ledge 2 (R = 26.82)
Circle(140) = {142, 140, 144};
Circle(141) = {144, 140, 143};
Circle(142) = {143, 140, 145};
Circle(143) = {145, 140, 142};

// Level 6 - Top of ledge 2 (R = 26.82)
Circle(150) = {152, 150, 154};
Circle(151) = {154, 150, 153};
Circle(152) = {153, 150, 155};
Circle(153) = {155, 150, 152};

// Level 7 - Top bulge (R = 30.65)
Circle(160) = {162, 160, 164};
Circle(161) = {164, 160, 163};
Circle(162) = {163, 160, 165};
Circle(163) = {165, 160, 162};

//////////////////////////////
// CAVERN SURFACES
//////////////////////////////

// Bottom cap (4 triangular surfaces)
Curve Loop(201) = {20, 100, -22};
Surface(201) = {201};
Curve Loop(202) = {22, 101, -21};
Surface(202) = {202};
Curve Loop(203) = {21, 102, -23};
Surface(203) = {203};
Curve Loop(204) = {23, 103, -20};
Surface(204) = {204};

// Side walls Level 1-2 (4 surfaces)
Curve Loop(211) = {100, 32, -110, -30};
Surface(211) = {211};
Curve Loop(212) = {101, 31, -111, -32};
Surface(212) = {212};
Curve Loop(213) = {102, 33, -112, -31};
Surface(213) = {213};
Curve Loop(214) = {103, 30, -113, -33};
Surface(214) = {214};

// Side walls Level 2-3 - VERTICAL ledge 1 (4 surfaces)
Curve Loop(221) = {110, 42, -120, -40};
Surface(221) = {221};
Curve Loop(222) = {111, 41, -121, -42};
Surface(222) = {222};
Curve Loop(223) = {112, 43, -122, -41};
Surface(223) = {223};
Curve Loop(224) = {113, 40, -123, -43};
Surface(224) = {224};

// Side walls Level 3-4 (4 surfaces)
Curve Loop(231) = {120, 52, -130, -50};
Surface(231) = {231};
Curve Loop(232) = {121, 51, -131, -52};
Surface(232) = {232};
Curve Loop(233) = {122, 53, -132, -51};
Surface(233) = {233};
Curve Loop(234) = {123, 50, -133, -53};
Surface(234) = {234};

// Side walls Level 4-5 (4 surfaces)
Curve Loop(241) = {130, 62, -140, -60};
Surface(241) = {241};
Curve Loop(242) = {131, 61, -141, -62};
Surface(242) = {242};
Curve Loop(243) = {132, 63, -142, -61};
Surface(243) = {243};
Curve Loop(244) = {133, 60, -143, -63};
Surface(244) = {244};

// Side walls Level 5-6 - VERTICAL ledge 2 (4 surfaces)
Curve Loop(251) = {140, 72, -150, -70};
Surface(251) = {251};
Curve Loop(252) = {141, 71, -151, -72};
Surface(252) = {252};
Curve Loop(253) = {142, 73, -152, -71};
Surface(253) = {253};
Curve Loop(254) = {143, 70, -153, -73};
Surface(254) = {254};

// Side walls Level 6-7 (4 surfaces)
Curve Loop(261) = {150, 82, -160, -80};
Surface(261) = {261};
Curve Loop(262) = {151, 81, -161, -82};
Surface(262) = {262};
Curve Loop(263) = {152, 83, -162, -81};
Surface(263) = {263};
Curve Loop(264) = {153, 80, -163, -83};
Surface(264) = {264};

// Top cap (4 triangular surfaces)
Curve Loop(271) = {160, 92, -90};
Surface(271) = {271};
Curve Loop(272) = {161, 91, -92};
Surface(272) = {272};
Curve Loop(273) = {162, 93, -91};
Surface(273) = {273};
Curve Loop(274) = {163, 90, -93};
Surface(274) = {274};

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
// VOLUME - SINGLE SALT BLOCK
//////////////////////////////

// Cavern surface loop (all 32 cavern surfaces)
Surface Loop(200) = {201, 202, 203, 204,
                     211, 212, 213, 214,
                     221, 222, 223, 224,
                     231, 232, 233, 234,
                     241, 242, 243, 244,
                     251, 252, 253, 254,
                     261, 262, 263, 264,
                     271, 272, 273, 274};

// Outer box surface loop
Surface Loop(201) = {1, 2, 3, 4, 5, 6};

// Salt volume (box minus cavern) - SINGLE HOMOGENEOUS VOLUME
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
                                  261, 262, 263, 264,
                                  271, 272, 273, 274};

Physical Volume("Salt", 28) = {1};

// Wall profile for visualization
Physical Curve("Wall_profile", 30) = {20, 30, 40, 50, 60, 70, 80, 90};
