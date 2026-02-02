//////////////////////////////////////////////////////////////
// ASYMMETRIC SHELF CAVERN - FULL 3D HETEROGENEOUS
// Target volume: 600,000 m³
//
// Design: Asymmetric cavern with offset bulges and shelf-like
// ledges created by 2 thin horizontal interlayers. The cavern
// center shifts laterally at each level, creating an irregular
// profile typical of solution mining in tilted heterogeneous layers.
//
// The interlayers create sharp shelves where leaching stops.
// Different radii in +x/-x and +y/-y directions create asymmetry.
//
// HETEROGENEOUS: 3 salt volumes separated by 2 thin interlayers
// - Salt_bottom (z=0 to 190m)
// - Interlayer_1 (z=190m to 195m) - 5m thick anhydrite
// - Salt_middle (z=195m to 280m)
// - Interlayer_2 (z=280m to 285m) - 5m thick mudstone
// - Salt_top (z=285m to 660m)
//////////////////////////////////////////////////////////////

coarse_size = 65;
fine_size = 4.5;

Lz = 660.0;
Ly = 450.0;
Lx = 450.0;

x_center = Lx/2;
y_center = Ly/2;

//////////////////////////////////////////////////////////////
// INTERLAYER Z-POSITIONS
//////////////////////////////////////////////////////////////

z_inter1_bot = 190.0;
z_inter1_top = 195.0;

z_inter2_bot = 280.0;
z_inter2_top = 285.0;

//////////////////////////////////////////////////////////////
// CAVERN GEOMETRY - Asymmetric profile with offset centers
// The cavern center shifts creating a tilted/leaning shape
//////////////////////////////////////////////////////////////

// Bottom zone - centered
z_bot_tip = 150.0;
z_bot_max = 168.0;
x_off_bot = 0.0;    // No offset at bottom
y_off_bot = 0.0;

// Asymmetric radii at bottom bulge
R_bot_px = 42.0;    // +x direction
R_bot_mx = 38.0;    // -x direction
R_bot_py = 45.0;    // +y direction
R_bot_my = 35.0;    // -y direction

// At interlayer 1 - shelf with offset
R_shelf1_px = 28.0;
R_shelf1_mx = 24.0;
R_shelf1_py = 30.0;
R_shelf1_my = 22.0;

// Middle zone - offset toward +x, +y
z_mid_max = 238.0;
x_off_mid = 8.0;    // Offset center
y_off_mid = 6.0;

R_mid_px = 50.0;    // Largest bulge
R_mid_mx = 42.0;
R_mid_py = 48.0;
R_mid_my = 40.0;

// At interlayer 2 - shelf with offset
R_shelf2_px = 26.0;
R_shelf2_mx = 22.0;
R_shelf2_py = 24.0;
R_shelf2_my = 20.0;

// Top zone - offset back
z_top_max = 320.0;
z_top_tip = 360.0;
x_off_top = -5.0;
y_off_top = 3.0;

R_top_px = 32.0;
R_top_mx = 28.0;
R_top_py = 30.0;
R_top_my = 26.0;

//////////////////////////////
// OUTER BOX CORNERS
//////////////////////////////

Point(1) = {0, 0, 0, coarse_size};
Point(2) = {Lx, 0, 0, coarse_size};
Point(3) = {Lx, Ly, 0, coarse_size};
Point(4) = {0, Ly, 0, coarse_size};

Point(11) = {0, 0, z_inter1_bot, coarse_size};
Point(12) = {Lx, 0, z_inter1_bot, coarse_size};
Point(13) = {Lx, Ly, z_inter1_bot, coarse_size};
Point(14) = {0, Ly, z_inter1_bot, coarse_size};

Point(21) = {0, 0, z_inter1_top, coarse_size};
Point(22) = {Lx, 0, z_inter1_top, coarse_size};
Point(23) = {Lx, Ly, z_inter1_top, coarse_size};
Point(24) = {0, Ly, z_inter1_top, coarse_size};

Point(31) = {0, 0, z_inter2_bot, coarse_size};
Point(32) = {Lx, 0, z_inter2_bot, coarse_size};
Point(33) = {Lx, Ly, z_inter2_bot, coarse_size};
Point(34) = {0, Ly, z_inter2_bot, coarse_size};

Point(41) = {0, 0, z_inter2_top, coarse_size};
Point(42) = {Lx, 0, z_inter2_top, coarse_size};
Point(43) = {Lx, Ly, z_inter2_top, coarse_size};
Point(44) = {0, Ly, z_inter2_top, coarse_size};

Point(51) = {0, 0, Lz, coarse_size};
Point(52) = {Lx, 0, Lz, coarse_size};
Point(53) = {Lx, Ly, Lz, coarse_size};
Point(54) = {0, Ly, Lz, coarse_size};

//////////////////////////////
// CAVERN POINTS - Asymmetric
//////////////////////////////

// Level 0 - Bottom tip (centered)
Point(100) = {x_center + x_off_bot, y_center + y_off_bot, z_bot_tip, fine_size};

// Level 1 - Bottom max bulge (asymmetric radii)
Point(101) = {x_center + x_off_bot, y_center + y_off_bot, z_bot_max, coarse_size};
Point(102) = {x_center + x_off_bot + R_bot_px, y_center + y_off_bot, z_bot_max, fine_size};
Point(103) = {x_center + x_off_bot - R_bot_mx, y_center + y_off_bot, z_bot_max, fine_size};
Point(104) = {x_center + x_off_bot, y_center + y_off_bot + R_bot_py, z_bot_max, fine_size};
Point(105) = {x_center + x_off_bot, y_center + y_off_bot - R_bot_my, z_bot_max, fine_size};

// Level 2 - Bottom of interlayer 1 (shelf)
Point(110) = {x_center + x_off_bot, y_center + y_off_bot, z_inter1_bot, coarse_size};
Point(112) = {x_center + x_off_bot + R_shelf1_px, y_center + y_off_bot, z_inter1_bot, fine_size};
Point(113) = {x_center + x_off_bot - R_shelf1_mx, y_center + y_off_bot, z_inter1_bot, fine_size};
Point(114) = {x_center + x_off_bot, y_center + y_off_bot + R_shelf1_py, z_inter1_bot, fine_size};
Point(115) = {x_center + x_off_bot, y_center + y_off_bot - R_shelf1_my, z_inter1_bot, fine_size};

// Level 3 - Top of interlayer 1 (same radius for vertical wall)
Point(120) = {x_center + x_off_bot, y_center + y_off_bot, z_inter1_top, coarse_size};
Point(122) = {x_center + x_off_bot + R_shelf1_px, y_center + y_off_bot, z_inter1_top, fine_size};
Point(123) = {x_center + x_off_bot - R_shelf1_mx, y_center + y_off_bot, z_inter1_top, fine_size};
Point(124) = {x_center + x_off_bot, y_center + y_off_bot + R_shelf1_py, z_inter1_top, fine_size};
Point(125) = {x_center + x_off_bot, y_center + y_off_bot - R_shelf1_my, z_inter1_top, fine_size};

// Level 4 - Middle max bulge (offset center, largest asymmetric bulge)
Point(130) = {x_center + x_off_mid, y_center + y_off_mid, z_mid_max, coarse_size};
Point(132) = {x_center + x_off_mid + R_mid_px, y_center + y_off_mid, z_mid_max, fine_size};
Point(133) = {x_center + x_off_mid - R_mid_mx, y_center + y_off_mid, z_mid_max, fine_size};
Point(134) = {x_center + x_off_mid, y_center + y_off_mid + R_mid_py, z_mid_max, fine_size};
Point(135) = {x_center + x_off_mid, y_center + y_off_mid - R_mid_my, z_mid_max, fine_size};

// Level 5 - Bottom of interlayer 2 (shelf)
Point(140) = {x_center + x_off_mid, y_center + y_off_mid, z_inter2_bot, coarse_size};
Point(142) = {x_center + x_off_mid + R_shelf2_px, y_center + y_off_mid, z_inter2_bot, fine_size};
Point(143) = {x_center + x_off_mid - R_shelf2_mx, y_center + y_off_mid, z_inter2_bot, fine_size};
Point(144) = {x_center + x_off_mid, y_center + y_off_mid + R_shelf2_py, z_inter2_bot, fine_size};
Point(145) = {x_center + x_off_mid, y_center + y_off_mid - R_shelf2_my, z_inter2_bot, fine_size};

// Level 6 - Top of interlayer 2
Point(150) = {x_center + x_off_mid, y_center + y_off_mid, z_inter2_top, coarse_size};
Point(152) = {x_center + x_off_mid + R_shelf2_px, y_center + y_off_mid, z_inter2_top, fine_size};
Point(153) = {x_center + x_off_mid - R_shelf2_mx, y_center + y_off_mid, z_inter2_top, fine_size};
Point(154) = {x_center + x_off_mid, y_center + y_off_mid + R_shelf2_py, z_inter2_top, fine_size};
Point(155) = {x_center + x_off_mid, y_center + y_off_mid - R_shelf2_my, z_inter2_top, fine_size};

// Level 7 - Top max bulge (offset back)
Point(160) = {x_center + x_off_top, y_center + y_off_top, z_top_max, coarse_size};
Point(162) = {x_center + x_off_top + R_top_px, y_center + y_off_top, z_top_max, fine_size};
Point(163) = {x_center + x_off_top - R_top_mx, y_center + y_off_top, z_top_max, fine_size};
Point(164) = {x_center + x_off_top, y_center + y_off_top + R_top_py, z_top_max, fine_size};
Point(165) = {x_center + x_off_top, y_center + y_off_top - R_top_my, z_top_max, fine_size};

// Level 8 - Top tip
Point(170) = {x_center + x_off_top, y_center + y_off_top, z_top_tip, fine_size};

//////////////////////////////
// HORIZONTAL EDGES
//////////////////////////////

Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};
Line(11) = {11, 12}; Line(12) = {12, 13}; Line(13) = {13, 14}; Line(14) = {14, 11};
Line(21) = {21, 22}; Line(22) = {22, 23}; Line(23) = {23, 24}; Line(24) = {24, 21};
Line(31) = {31, 32}; Line(32) = {32, 33}; Line(33) = {33, 34}; Line(34) = {34, 31};
Line(41) = {41, 42}; Line(42) = {42, 43}; Line(43) = {43, 44}; Line(44) = {44, 41};
Line(51) = {51, 52}; Line(52) = {52, 53}; Line(53) = {53, 54}; Line(54) = {54, 51};

//////////////////////////////
// VERTICAL EDGES (DOMAIN BOUNDARIES)
//////////////////////////////

Line(61) = {1, 11}; Line(62) = {11, 21}; Line(63) = {21, 31}; Line(64) = {31, 41}; Line(65) = {41, 51};
Line(71) = {2, 12}; Line(72) = {12, 22}; Line(73) = {22, 32}; Line(74) = {32, 42}; Line(75) = {42, 52};
Line(81) = {3, 13}; Line(82) = {13, 23}; Line(83) = {23, 33}; Line(84) = {33, 43}; Line(85) = {43, 53};
Line(91) = {4, 14}; Line(92) = {14, 24}; Line(93) = {24, 34}; Line(94) = {34, 44}; Line(95) = {44, 54};

//////////////////////////////
// CAVERN VERTICAL CONNECTIONS
//////////////////////////////

Line(200) = {100, 102}; Line(201) = {100, 103}; Line(202) = {100, 104}; Line(203) = {100, 105};
Line(210) = {102, 112}; Line(211) = {103, 113}; Line(212) = {104, 114}; Line(213) = {105, 115};
Line(220) = {112, 122}; Line(221) = {113, 123}; Line(222) = {114, 124}; Line(223) = {115, 125};
Line(230) = {122, 132}; Line(231) = {123, 133}; Line(232) = {124, 134}; Line(233) = {125, 135};
Line(240) = {132, 142}; Line(241) = {133, 143}; Line(242) = {134, 144}; Line(243) = {135, 145};
Line(250) = {142, 152}; Line(251) = {143, 153}; Line(252) = {144, 154}; Line(253) = {145, 155};
Line(260) = {152, 162}; Line(261) = {153, 163}; Line(262) = {154, 164}; Line(263) = {155, 165};
Line(270) = {162, 170}; Line(271) = {163, 170}; Line(272) = {164, 170}; Line(273) = {165, 170};

//////////////////////////////
// ELLIPSE ARCS AT EACH CAVERN LEVEL (asymmetric)
//////////////////////////////

// Level 1 - Bottom bulge (use Splines for asymmetric shape)
Ellipse(300) = {102, 101, 102, 104};
Ellipse(301) = {104, 101, 104, 103};
Ellipse(302) = {103, 101, 103, 105};
Ellipse(303) = {105, 101, 105, 102};

// Level 2 - Shelf 1 bottom
Ellipse(310) = {112, 110, 112, 114};
Ellipse(311) = {114, 110, 114, 113};
Ellipse(312) = {113, 110, 113, 115};
Ellipse(313) = {115, 110, 115, 112};

// Level 3 - Shelf 1 top
Ellipse(320) = {122, 120, 122, 124};
Ellipse(321) = {124, 120, 124, 123};
Ellipse(322) = {123, 120, 123, 125};
Ellipse(323) = {125, 120, 125, 122};

// Level 4 - Middle bulge
Ellipse(330) = {132, 130, 132, 134};
Ellipse(331) = {134, 130, 134, 133};
Ellipse(332) = {133, 130, 133, 135};
Ellipse(333) = {135, 130, 135, 132};

// Level 5 - Shelf 2 bottom
Ellipse(340) = {142, 140, 142, 144};
Ellipse(341) = {144, 140, 144, 143};
Ellipse(342) = {143, 140, 143, 145};
Ellipse(343) = {145, 140, 145, 142};

// Level 6 - Shelf 2 top
Ellipse(350) = {152, 150, 152, 154};
Ellipse(351) = {154, 150, 154, 153};
Ellipse(352) = {153, 150, 153, 155};
Ellipse(353) = {155, 150, 155, 152};

// Level 7 - Top bulge
Ellipse(360) = {162, 160, 162, 164};
Ellipse(361) = {164, 160, 164, 163};
Ellipse(362) = {163, 160, 163, 165};
Ellipse(363) = {165, 160, 165, 162};

//////////////////////////////
// CAVERN SURFACES
//////////////////////////////

// Bottom cap
Curve Loop(401) = {200, 300, -202}; Surface(401) = {401};
Curve Loop(402) = {202, 301, -201}; Surface(402) = {402};
Curve Loop(403) = {201, 302, -203}; Surface(403) = {403};
Curve Loop(404) = {203, 303, -200}; Surface(404) = {404};

// Walls L1-L2
Curve Loop(411) = {300, 212, -310, -210}; Surface(411) = {411};
Curve Loop(412) = {301, 211, -311, -212}; Surface(412) = {412};
Curve Loop(413) = {302, 213, -312, -211}; Surface(413) = {413};
Curve Loop(414) = {303, 210, -313, -213}; Surface(414) = {414};

// Walls L2-L3 (interlayer 1)
Curve Loop(421) = {310, 222, -320, -220}; Surface(421) = {421};
Curve Loop(422) = {311, 221, -321, -222}; Surface(422) = {422};
Curve Loop(423) = {312, 223, -322, -221}; Surface(423) = {423};
Curve Loop(424) = {313, 220, -323, -223}; Surface(424) = {424};

// Walls L3-L4
Curve Loop(431) = {320, 232, -330, -230}; Surface(431) = {431};
Curve Loop(432) = {321, 231, -331, -232}; Surface(432) = {432};
Curve Loop(433) = {322, 233, -332, -231}; Surface(433) = {433};
Curve Loop(434) = {323, 230, -333, -233}; Surface(434) = {434};

// Walls L4-L5
Curve Loop(441) = {330, 242, -340, -240}; Surface(441) = {441};
Curve Loop(442) = {331, 241, -341, -242}; Surface(442) = {442};
Curve Loop(443) = {332, 243, -342, -241}; Surface(443) = {443};
Curve Loop(444) = {333, 240, -343, -243}; Surface(444) = {444};

// Walls L5-L6 (interlayer 2)
Curve Loop(451) = {340, 252, -350, -250}; Surface(451) = {451};
Curve Loop(452) = {341, 251, -351, -252}; Surface(452) = {452};
Curve Loop(453) = {342, 253, -352, -251}; Surface(453) = {453};
Curve Loop(454) = {343, 250, -353, -253}; Surface(454) = {454};

// Walls L6-L7
Curve Loop(461) = {350, 262, -360, -260}; Surface(461) = {461};
Curve Loop(462) = {351, 261, -361, -262}; Surface(462) = {462};
Curve Loop(463) = {352, 263, -362, -261}; Surface(463) = {463};
Curve Loop(464) = {353, 260, -363, -263}; Surface(464) = {464};

// Top cap
Curve Loop(471) = {360, 272, -270}; Surface(471) = {471};
Curve Loop(472) = {361, 271, -272}; Surface(472) = {472};
Curve Loop(473) = {362, 273, -271}; Surface(473) = {473};
Curve Loop(474) = {363, 270, -273}; Surface(474) = {474};

//////////////////////////////
// DOMAIN BOUNDARY SURFACES
//////////////////////////////

Curve Loop(501) = {1, 2, 3, 4}; Plane Surface(501) = {501};
Curve Loop(502) = {51, 52, 53, 54}; Plane Surface(502) = {502};

Curve Loop(511) = {1, 71, -11, -61}; Plane Surface(511) = {511};
Curve Loop(512) = {11, 72, -21, -62}; Plane Surface(512) = {512};
Curve Loop(513) = {21, 73, -31, -63}; Plane Surface(513) = {513};
Curve Loop(514) = {31, 74, -41, -64}; Plane Surface(514) = {514};
Curve Loop(515) = {41, 75, -51, -65}; Plane Surface(515) = {515};

Curve Loop(521) = {-3, 81, 13, -91}; Plane Surface(521) = {521};
Curve Loop(522) = {-13, 82, 23, -92}; Plane Surface(522) = {522};
Curve Loop(523) = {-23, 83, 33, -93}; Plane Surface(523) = {523};
Curve Loop(524) = {-33, 84, 43, -94}; Plane Surface(524) = {524};
Curve Loop(525) = {-43, 85, 53, -95}; Plane Surface(525) = {525};

Curve Loop(531) = {-4, 91, 14, -61}; Plane Surface(531) = {531};
Curve Loop(532) = {-14, 92, 24, -62}; Plane Surface(532) = {532};
Curve Loop(533) = {-24, 93, 34, -63}; Plane Surface(533) = {533};
Curve Loop(534) = {-34, 94, 44, -64}; Plane Surface(534) = {534};
Curve Loop(535) = {-44, 95, 54, -65}; Plane Surface(535) = {535};

Curve Loop(541) = {2, 81, -12, -71}; Plane Surface(541) = {541};
Curve Loop(542) = {12, 82, -22, -72}; Plane Surface(542) = {542};
Curve Loop(543) = {22, 83, -32, -73}; Plane Surface(543) = {543};
Curve Loop(544) = {32, 84, -42, -74}; Plane Surface(544) = {544};
Curve Loop(545) = {42, 85, -52, -75}; Plane Surface(545) = {545};

//////////////////////////////
// INTERNAL HORIZONTAL SURFACES
//////////////////////////////

Curve Loop(551) = {11, 12, 13, 14};
Curve Loop(552) = {310, 311, 312, 313};
Plane Surface(551) = {551, 552};

Curve Loop(561) = {21, 22, 23, 24};
Curve Loop(562) = {320, 321, 322, 323};
Plane Surface(561) = {561, 562};

Curve Loop(571) = {31, 32, 33, 34};
Curve Loop(572) = {340, 341, 342, 343};
Plane Surface(571) = {571, 572};

Curve Loop(581) = {41, 42, 43, 44};
Curve Loop(582) = {350, 351, 352, 353};
Plane Surface(581) = {581, 582};

//////////////////////////////
// VOLUMES
//////////////////////////////

Surface Loop(601) = {501, 511, 521, 531, 541, 551,
                     401, 402, 403, 404, 411, 412, 413, 414};
Volume(1) = {601};

Surface Loop(602) = {551, 512, 522, 532, 542, 561,
                     421, 422, 423, 424};
Volume(2) = {602};

Surface Loop(603) = {561, 513, 523, 533, 543, 571,
                     431, 432, 433, 434, 441, 442, 443, 444};
Volume(3) = {603};

Surface Loop(604) = {571, 514, 524, 534, 544, 581,
                     451, 452, 453, 454};
Volume(4) = {604};

Surface Loop(605) = {581, 515, 525, 535, 545, 502,
                     461, 462, 463, 464, 471, 472, 473, 474};
Volume(5) = {605};

//////////////////////////////
// PHYSICAL GROUPS
//////////////////////////////

Physical Surface("Bottom", 27) = {501};
Physical Surface("Top", 22) = {502};
Physical Surface("South", 23) = {511, 512, 513, 514, 515};
Physical Surface("North", 24) = {521, 522, 523, 524, 525};
Physical Surface("West", 26) = {531, 532, 533, 534, 535};
Physical Surface("East", 25) = {541, 542, 543, 544, 545};

Physical Surface("Cavern", 29) = {401, 402, 403, 404,
                                  411, 412, 413, 414,
                                  421, 422, 423, 424,
                                  431, 432, 433, 434,
                                  441, 442, 443, 444,
                                  451, 452, 453, 454,
                                  461, 462, 463, 464,
                                  471, 472, 473, 474};

Physical Volume("Salt_bottom", 31) = {1};
Physical Volume("Interlayer_1", 32) = {2};
Physical Volume("Salt_middle", 33) = {3};
Physical Volume("Interlayer_2", 34) = {4};
Physical Volume("Salt_top", 35) = {5};

Physical Curve("Wall_profile", 30) = {200, 210, 220, 230, 240, 250, 260, 270};
