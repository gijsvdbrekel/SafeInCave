// ============================================================================
// FAST-LEACHED CAVERN - SCALED TO 1.2 MILLION m³
// ============================================================================
// Rough barrel shape with oscillating radii simulating fast/uncontrolled leaching
// Target volume:   1.2 Mm³
// Analytical volume: 1200000 m³
// Domain dimensions: 450×450×660 m
// ============================================================================

size_coarse = 65;
size_fine   = 4.5;

Lz = 660.0;
Ly = 450.0;
Lx = 450.0;

x_center = Lx/2;  // 225.0
y_center = Ly/2;  // 225.0

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
// CAVERN GEOMETRY - 13 LEVELS
//////////////////////////////

// Level 1 - bottom tip (z=187.92)
Point(100) = {225.000000, 225.000000, 187.924261, size_fine};  // tip
Point(101) = {225.000000, 225.000000, 192.004098, size_coarse};  // center
Point(102) = {233.975639, 225.000000, 192.004098, size_fine};  // +x
Point(103) = {216.024361, 225.000000, 192.004098, size_fine};  // -x
Point(104) = {225.000000, 233.975639, 192.004098, size_fine};  // +y
Point(105) = {225.000000, 216.024361, 192.004098, size_fine};  // -y

// Level 2 (z=196.08, R=22.44)
Point(110) = {225.000000, 225.000000, 196.083934, size_coarse};  // center
Point(112) = {247.439099, 225.000000, 196.083934, size_fine};  // +x
Point(113) = {202.560901, 225.000000, 196.083934, size_fine};  // -x
Point(114) = {225.000000, 247.439099, 196.083934, size_fine};  // +y
Point(115) = {225.000000, 202.560901, 196.083934, size_fine};  // -y

// Level 3 (z=216.48, R=53.04)
Point(120) = {225.000000, 225.000000, 216.483114, size_coarse};  // center
Point(122) = {278.037869, 225.000000, 216.483114, size_fine};  // +x
Point(123) = {171.962131, 225.000000, 216.483114, size_fine};  // -x
Point(124) = {225.000000, 278.037869, 216.483114, size_fine};  // +y
Point(125) = {225.000000, 171.962131, 216.483114, size_fine};  // -y

// Level 4 (z=236.88, R=36.72)
Point(130) = {225.000000, 225.000000, 236.882295, size_coarse};  // center
Point(132) = {261.718525, 225.000000, 236.882295, size_fine};  // +x
Point(133) = {188.281475, 225.000000, 236.882295, size_fine};  // -x
Point(134) = {225.000000, 261.718525, 236.882295, size_fine};  // +y
Point(135) = {225.000000, 188.281475, 236.882295, size_fine};  // -y

// Level 5 (z=257.28, R=57.12)
Point(140) = {225.000000, 225.000000, 257.281475, size_coarse};  // center
Point(142) = {282.117705, 225.000000, 257.281475, size_fine};  // +x
Point(143) = {167.882295, 225.000000, 257.281475, size_fine};  // -x
Point(144) = {225.000000, 282.117705, 257.281475, size_fine};  // +y
Point(145) = {225.000000, 167.882295, 257.281475, size_fine};  // -y

// Level 6 (z=277.68, R=38.76)
Point(150) = {225.000000, 225.000000, 277.680656, size_coarse};  // center
Point(152) = {263.758443, 225.000000, 277.680656, size_fine};  // +x
Point(153) = {186.241557, 225.000000, 277.680656, size_fine};  // -x
Point(154) = {225.000000, 263.758443, 277.680656, size_fine};  // +y
Point(155) = {225.000000, 186.241557, 277.680656, size_fine};  // -y

// Level 7 (z=298.08, R=55.08)
Point(160) = {225.000000, 225.000000, 298.079836, size_coarse};  // center
Point(162) = {280.077787, 225.000000, 298.079836, size_fine};  // +x
Point(163) = {169.922213, 225.000000, 298.079836, size_fine};  // -x
Point(164) = {225.000000, 280.077787, 298.079836, size_fine};  // +y
Point(165) = {225.000000, 169.922213, 298.079836, size_fine};  // -y

// Level 8 (z=318.48, R=40.80)
Point(170) = {225.000000, 225.000000, 318.479017, size_coarse};  // center
Point(172) = {265.798361, 225.000000, 318.479017, size_fine};  // +x
Point(173) = {184.201639, 225.000000, 318.479017, size_fine};  // -x
Point(174) = {225.000000, 265.798361, 318.479017, size_fine};  // +y
Point(175) = {225.000000, 184.201639, 318.479017, size_fine};  // -y

// Level 9 (z=338.88, R=51.00)
Point(180) = {225.000000, 225.000000, 338.878197, size_coarse};  // center
Point(182) = {275.997951, 225.000000, 338.878197, size_fine};  // +x
Point(183) = {174.002049, 225.000000, 338.878197, size_fine};  // -x
Point(184) = {225.000000, 275.997951, 338.878197, size_fine};  // +y
Point(185) = {225.000000, 174.002049, 338.878197, size_fine};  // -y

// Level 10 (z=359.28, R=34.68)
Point(190) = {225.000000, 225.000000, 359.277378, size_coarse};  // center
Point(192) = {259.678607, 225.000000, 359.277378, size_fine};  // +x
Point(193) = {190.321393, 225.000000, 359.277378, size_fine};  // -x
Point(194) = {225.000000, 259.678607, 359.277378, size_fine};  // +y
Point(195) = {225.000000, 190.321393, 359.277378, size_fine};  // -y

// Level 11 (z=376.62, R=46.92)
Point(200) = {225.000000, 225.000000, 376.616681, size_coarse};  // center
Point(202) = {271.918115, 225.000000, 376.616681, size_fine};  // +x
Point(203) = {178.081885, 225.000000, 376.616681, size_fine};  // -x
Point(204) = {225.000000, 271.918115, 376.616681, size_fine};  // +y
Point(205) = {225.000000, 178.081885, 376.616681, size_fine};  // -y

// Level 12 (z=389.88, R=10.20)
Point(210) = {225.000000, 225.000000, 389.876148, size_coarse};  // center
Point(212) = {235.199590, 225.000000, 389.876148, size_fine};  // +x
Point(213) = {214.800410, 225.000000, 389.876148, size_fine};  // -x
Point(214) = {225.000000, 235.199590, 389.876148, size_fine};  // +y
Point(215) = {225.000000, 214.800410, 389.876148, size_fine};  // -y

// Level 13 - top tip (z=400.08)
Point(220) = {225.000000, 225.000000, 400.075739, size_fine};  // tip
Point(221) = {225.000000, 225.000000, 394.975943, size_coarse};  // center
Point(222) = {229.079836, 225.000000, 394.975943, size_fine};  // +x
Point(223) = {220.920164, 225.000000, 394.975943, size_fine};  // -x
Point(224) = {225.000000, 229.079836, 394.975943, size_fine};  // +y
Point(225) = {225.000000, 220.920164, 394.975943, size_fine};  // -y

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
// CAVERN LINES AND ARCS
//////////////////////////////

// Bottom tip to ring
Line(20) = {100, 102};
Line(21) = {100, 104};
Line(22) = {100, 103};
Line(23) = {100, 105};

// Bottom ring circles
Circle(24) = {102, 101, 104};
Circle(25) = {104, 101, 103};
Circle(26) = {103, 101, 105};
Circle(27) = {105, 101, 102};

// Vertical connections between levels
Line(28) = {102, 112};  // +x level 1->2
Line(29) = {103, 113};  // -x level 1->2
Line(30) = {104, 114};  // +y level 1->2
Line(31) = {105, 115};  // -y level 1->2
Line(32) = {112, 122};  // +x level 2->3
Line(33) = {113, 123};  // -x level 2->3
Line(34) = {114, 124};  // +y level 2->3
Line(35) = {115, 125};  // -y level 2->3
Line(36) = {122, 132};  // +x level 3->4
Line(37) = {123, 133};  // -x level 3->4
Line(38) = {124, 134};  // +y level 3->4
Line(39) = {125, 135};  // -y level 3->4
Line(40) = {132, 142};  // +x level 4->5
Line(41) = {133, 143};  // -x level 4->5
Line(42) = {134, 144};  // +y level 4->5
Line(43) = {135, 145};  // -y level 4->5
Line(44) = {142, 152};  // +x level 5->6
Line(45) = {143, 153};  // -x level 5->6
Line(46) = {144, 154};  // +y level 5->6
Line(47) = {145, 155};  // -y level 5->6
Line(48) = {152, 162};  // +x level 6->7
Line(49) = {153, 163};  // -x level 6->7
Line(50) = {154, 164};  // +y level 6->7
Line(51) = {155, 165};  // -y level 6->7
Line(52) = {162, 172};  // +x level 7->8
Line(53) = {163, 173};  // -x level 7->8
Line(54) = {164, 174};  // +y level 7->8
Line(55) = {165, 175};  // -y level 7->8
Line(56) = {172, 182};  // +x level 8->9
Line(57) = {173, 183};  // -x level 8->9
Line(58) = {174, 184};  // +y level 8->9
Line(59) = {175, 185};  // -y level 8->9
Line(60) = {182, 192};  // +x level 9->10
Line(61) = {183, 193};  // -x level 9->10
Line(62) = {184, 194};  // +y level 9->10
Line(63) = {185, 195};  // -y level 9->10
Line(64) = {192, 202};  // +x level 10->11
Line(65) = {193, 203};  // -x level 10->11
Line(66) = {194, 204};  // +y level 10->11
Line(67) = {195, 205};  // -y level 10->11
Line(68) = {202, 212};  // +x level 11->12
Line(69) = {203, 213};  // -x level 11->12
Line(70) = {204, 214};  // +y level 11->12
Line(71) = {205, 215};  // -y level 11->12
Line(72) = {212, 222};  // +x level 12->13
Line(73) = {213, 223};  // -x level 12->13
Line(74) = {214, 224};  // +y level 12->13
Line(75) = {215, 225};  // -y level 12->13

// Circle arcs at intermediate levels
Circle(76) = {112, 110, 114};  // level 2
Circle(77) = {114, 110, 113};
Circle(78) = {113, 110, 115};
Circle(79) = {115, 110, 112};
Circle(80) = {122, 120, 124};  // level 3
Circle(81) = {124, 120, 123};
Circle(82) = {123, 120, 125};
Circle(83) = {125, 120, 122};
Circle(84) = {132, 130, 134};  // level 4
Circle(85) = {134, 130, 133};
Circle(86) = {133, 130, 135};
Circle(87) = {135, 130, 132};
Circle(88) = {142, 140, 144};  // level 5
Circle(89) = {144, 140, 143};
Circle(90) = {143, 140, 145};
Circle(91) = {145, 140, 142};
Circle(92) = {152, 150, 154};  // level 6
Circle(93) = {154, 150, 153};
Circle(94) = {153, 150, 155};
Circle(95) = {155, 150, 152};
Circle(96) = {162, 160, 164};  // level 7
Circle(97) = {164, 160, 163};
Circle(98) = {163, 160, 165};
Circle(99) = {165, 160, 162};
Circle(100) = {172, 170, 174};  // level 8
Circle(101) = {174, 170, 173};
Circle(102) = {173, 170, 175};
Circle(103) = {175, 170, 172};
Circle(104) = {182, 180, 184};  // level 9
Circle(105) = {184, 180, 183};
Circle(106) = {183, 180, 185};
Circle(107) = {185, 180, 182};
Circle(108) = {192, 190, 194};  // level 10
Circle(109) = {194, 190, 193};
Circle(110) = {193, 190, 195};
Circle(111) = {195, 190, 192};
Circle(112) = {202, 200, 204};  // level 11
Circle(113) = {204, 200, 203};
Circle(114) = {203, 200, 205};
Circle(115) = {205, 200, 202};
Circle(116) = {212, 210, 214};  // level 12
Circle(117) = {214, 210, 213};
Circle(118) = {213, 210, 215};
Circle(119) = {215, 210, 212};

// Top tip to ring
Line(120) = {220, 222};
Line(121) = {220, 224};
Line(122) = {220, 223};
Line(123) = {220, 225};

// Top ring circles
Circle(124) = {222, 221, 224};
Circle(125) = {224, 221, 223};
Circle(126) = {223, 221, 225};
Circle(127) = {225, 221, 222};

//////////////////////////////
// CAVERN SURFACES
//////////////////////////////

// Bottom cap
Curve Loop(201) = {20, 24, -21};
Surface(201) = {201};
Curve Loop(202) = {21, 25, -22};
Surface(202) = {202};
Curve Loop(203) = {22, 26, -23};
Surface(203) = {203};
Curve Loop(204) = {23, 27, -20};
Surface(204) = {204};

// Side walls: bottom ring to level 2
Curve Loop(205) = {24, 30, -76, -28};
Surface(205) = {205};
Curve Loop(206) = {25, 29, -77, -30};
Surface(206) = {206};
Curve Loop(207) = {26, 31, -78, -29};
Surface(207) = {207};
Curve Loop(208) = {27, 28, -79, -31};
Surface(208) = {208};

// Side walls: level 2 to level 3
Curve Loop(209) = {76, 34, -80, -32};
Surface(209) = {209};
Curve Loop(210) = {77, 33, -81, -34};
Surface(210) = {210};
Curve Loop(211) = {78, 35, -82, -33};
Surface(211) = {211};
Curve Loop(212) = {79, 32, -83, -35};
Surface(212) = {212};

// Side walls: level 3 to level 4
Curve Loop(213) = {80, 38, -84, -36};
Surface(213) = {213};
Curve Loop(214) = {81, 37, -85, -38};
Surface(214) = {214};
Curve Loop(215) = {82, 39, -86, -37};
Surface(215) = {215};
Curve Loop(216) = {83, 36, -87, -39};
Surface(216) = {216};

// Side walls: level 4 to level 5
Curve Loop(217) = {84, 42, -88, -40};
Surface(217) = {217};
Curve Loop(218) = {85, 41, -89, -42};
Surface(218) = {218};
Curve Loop(219) = {86, 43, -90, -41};
Surface(219) = {219};
Curve Loop(220) = {87, 40, -91, -43};
Surface(220) = {220};

// Side walls: level 5 to level 6
Curve Loop(221) = {88, 46, -92, -44};
Surface(221) = {221};
Curve Loop(222) = {89, 45, -93, -46};
Surface(222) = {222};
Curve Loop(223) = {90, 47, -94, -45};
Surface(223) = {223};
Curve Loop(224) = {91, 44, -95, -47};
Surface(224) = {224};

// Side walls: level 6 to level 7
Curve Loop(225) = {92, 50, -96, -48};
Surface(225) = {225};
Curve Loop(226) = {93, 49, -97, -50};
Surface(226) = {226};
Curve Loop(227) = {94, 51, -98, -49};
Surface(227) = {227};
Curve Loop(228) = {95, 48, -99, -51};
Surface(228) = {228};

// Side walls: level 7 to level 8
Curve Loop(229) = {96, 54, -100, -52};
Surface(229) = {229};
Curve Loop(230) = {97, 53, -101, -54};
Surface(230) = {230};
Curve Loop(231) = {98, 55, -102, -53};
Surface(231) = {231};
Curve Loop(232) = {99, 52, -103, -55};
Surface(232) = {232};

// Side walls: level 8 to level 9
Curve Loop(233) = {100, 58, -104, -56};
Surface(233) = {233};
Curve Loop(234) = {101, 57, -105, -58};
Surface(234) = {234};
Curve Loop(235) = {102, 59, -106, -57};
Surface(235) = {235};
Curve Loop(236) = {103, 56, -107, -59};
Surface(236) = {236};

// Side walls: level 9 to level 10
Curve Loop(237) = {104, 62, -108, -60};
Surface(237) = {237};
Curve Loop(238) = {105, 61, -109, -62};
Surface(238) = {238};
Curve Loop(239) = {106, 63, -110, -61};
Surface(239) = {239};
Curve Loop(240) = {107, 60, -111, -63};
Surface(240) = {240};

// Side walls: level 10 to level 11
Curve Loop(241) = {108, 66, -112, -64};
Surface(241) = {241};
Curve Loop(242) = {109, 65, -113, -66};
Surface(242) = {242};
Curve Loop(243) = {110, 67, -114, -65};
Surface(243) = {243};
Curve Loop(244) = {111, 64, -115, -67};
Surface(244) = {244};

// Side walls: level 11 to level 12
Curve Loop(245) = {112, 70, -116, -68};
Surface(245) = {245};
Curve Loop(246) = {113, 69, -117, -70};
Surface(246) = {246};
Curve Loop(247) = {114, 71, -118, -69};
Surface(247) = {247};
Curve Loop(248) = {115, 68, -119, -71};
Surface(248) = {248};

// Side walls: level 12 to level 13
Curve Loop(249) = {116, 74, -124, -72};
Surface(249) = {249};
Curve Loop(250) = {117, 73, -125, -74};
Surface(250) = {250};
Curve Loop(251) = {118, 75, -126, -73};
Surface(251) = {251};
Curve Loop(252) = {119, 72, -127, -75};
Surface(252) = {252};

// Top cap
Curve Loop(253) = {120, 124, -121};
Surface(253) = {253};
Curve Loop(254) = {121, 125, -122};
Surface(254) = {254};
Curve Loop(255) = {122, 126, -123};
Surface(255) = {255};
Curve Loop(256) = {123, 127, -120};
Surface(256) = {256};

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

// Cavern surface loop (56 surfaces)
Surface Loop(200) = {201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256};

// Outer box surface loop
Surface Loop(201) = {1, 2, 3, 4, 5, 6};

// Salt volume (box minus cavern)
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

Physical Surface("Cavern", 29) = {201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256};
Physical Volume("Salt",    28) = {1};

Physical Curve("Wall_profile", 30) = {120, 72, 68, 64, 60, 56, 52, 48, 44, 40, 36, 32, 28, 20};
