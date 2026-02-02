//////////////////////////////////////////////////////////////
// VERTICAL INTRUSION CAVERN - FULL 3D HETEROGENEOUS
// Target volume: 600,000 m³
//
// Design: Cavern with 2 horizontal interlayers that create
// ledges, PLUS an asymmetric shape in the upper zone that
// represents the effect of a vertical interlayer (anhydrite vein)
// intersecting the leaching zone from above.
//
// The vertical interlayer effect is represented by:
// - Reduced radius in the +x direction in the upper zone
// - Sharp notch where the vertical layer stopped leaching
//
// This creates a realistic irregular cavern shape where
// heterogeneous layers (both horizontal and vertical) have
// affected the solution mining process.
//
// HETEROGENEOUS: 3 salt volumes separated by 2 thin interlayers
// - Salt_bottom (z=0 to 185m)
// - Interlayer_1 (z=185m to 190m) - 5m thick mudstone
// - Salt_middle (z=190m to 285m)
// - Interlayer_2 (z=285m to 290m) - 5m thick anhydrite
// - Salt_top (z=290m to 660m)
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

z_inter1_bot = 185.0;
z_inter1_top = 190.0;

z_inter2_bot = 285.0;
z_inter2_top = 290.0;

//////////////////////////////////////////////////////////////
// CAVERN GEOMETRY
// Upper zone has asymmetric shape due to vertical interlayer effect
//////////////////////////////////////////////////////////////

z_bot_tip = 150.0;
z_bot_max = 168.0;
R_bottom = 40.0;

R_ledge1 = 26.0;

z_mid_max = 238.0;
R_middle = 46.0;

R_ledge2 = 24.0;

// Upper zone - asymmetric due to vertical interlayer
// +x direction has reduced radius (notch from vertical layer)
z_upper_max = 335.0;
R_upper_px = 18.0;    // Reduced - notch from vertical interlayer
R_upper_mx = 32.0;    // Normal
R_upper_py = 30.0;    // Normal
R_upper_my = 28.0;    // Normal

z_top_tip = 375.0;

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
// CAVERN POINTS
//////////////////////////////

// Level 0 - Bottom tip
Point(100) = {x_center, y_center, z_bot_tip, fine_size};

// Level 1 - Bottom max bulge
Point(101) = {x_center, y_center, z_bot_max, coarse_size};
Point(102) = {x_center + R_bottom, y_center, z_bot_max, fine_size};
Point(103) = {x_center - R_bottom, y_center, z_bot_max, fine_size};
Point(104) = {x_center, y_center + R_bottom, z_bot_max, fine_size};
Point(105) = {x_center, y_center - R_bottom, z_bot_max, fine_size};

// Level 2 - Horizontal interlayer 1 bottom (ledge)
Point(110) = {x_center, y_center, z_inter1_bot, coarse_size};
Point(112) = {x_center + R_ledge1, y_center, z_inter1_bot, fine_size};
Point(113) = {x_center - R_ledge1, y_center, z_inter1_bot, fine_size};
Point(114) = {x_center, y_center + R_ledge1, z_inter1_bot, fine_size};
Point(115) = {x_center, y_center - R_ledge1, z_inter1_bot, fine_size};

// Level 3 - Horizontal interlayer 1 top
Point(120) = {x_center, y_center, z_inter1_top, coarse_size};
Point(122) = {x_center + R_ledge1, y_center, z_inter1_top, fine_size};
Point(123) = {x_center - R_ledge1, y_center, z_inter1_top, fine_size};
Point(124) = {x_center, y_center + R_ledge1, z_inter1_top, fine_size};
Point(125) = {x_center, y_center - R_ledge1, z_inter1_top, fine_size};

// Level 4 - Middle max bulge
Point(130) = {x_center, y_center, z_mid_max, coarse_size};
Point(132) = {x_center + R_middle, y_center, z_mid_max, fine_size};
Point(133) = {x_center - R_middle, y_center, z_mid_max, fine_size};
Point(134) = {x_center, y_center + R_middle, z_mid_max, fine_size};
Point(135) = {x_center, y_center - R_middle, z_mid_max, fine_size};

// Level 5 - Horizontal interlayer 2 bottom (ledge)
Point(140) = {x_center, y_center, z_inter2_bot, coarse_size};
Point(142) = {x_center + R_ledge2, y_center, z_inter2_bot, fine_size};
Point(143) = {x_center - R_ledge2, y_center, z_inter2_bot, fine_size};
Point(144) = {x_center, y_center + R_ledge2, z_inter2_bot, fine_size};
Point(145) = {x_center, y_center - R_ledge2, z_inter2_bot, fine_size};

// Level 6 - Horizontal interlayer 2 top
Point(150) = {x_center, y_center, z_inter2_top, coarse_size};
Point(152) = {x_center + R_ledge2, y_center, z_inter2_top, fine_size};
Point(153) = {x_center - R_ledge2, y_center, z_inter2_top, fine_size};
Point(154) = {x_center, y_center + R_ledge2, z_inter2_top, fine_size};
Point(155) = {x_center, y_center - R_ledge2, z_inter2_top, fine_size};

// Level 7 - Upper max bulge (ASYMMETRIC - vertical interlayer effect)
Point(160) = {x_center, y_center, z_upper_max, coarse_size};
Point(162) = {x_center + R_upper_px, y_center, z_upper_max, fine_size};  // Notch!
Point(163) = {x_center - R_upper_mx, y_center, z_upper_max, fine_size};
Point(164) = {x_center, y_center + R_upper_py, z_upper_max, fine_size};
Point(165) = {x_center, y_center - R_upper_my, z_upper_max, fine_size};

// Level 8 - Top tip
Point(170) = {x_center, y_center, z_top_tip, fine_size};

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
// VERTICAL EDGES
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
// CIRCLE/ELLIPSE ARCS
//////////////////////////////

// Level 1 - Bottom bulge (symmetric)
Circle(300) = {102, 101, 104}; Circle(301) = {104, 101, 103};
Circle(302) = {103, 101, 105}; Circle(303) = {105, 101, 102};

// Level 2 - Ledge 1 bottom
Circle(310) = {112, 110, 114}; Circle(311) = {114, 110, 113};
Circle(312) = {113, 110, 115}; Circle(313) = {115, 110, 112};

// Level 3 - Ledge 1 top
Circle(320) = {122, 120, 124}; Circle(321) = {124, 120, 123};
Circle(322) = {123, 120, 125}; Circle(323) = {125, 120, 122};

// Level 4 - Middle bulge (symmetric)
Circle(330) = {132, 130, 134}; Circle(331) = {134, 130, 133};
Circle(332) = {133, 130, 135}; Circle(333) = {135, 130, 132};

// Level 5 - Ledge 2 bottom
Circle(340) = {142, 140, 144}; Circle(341) = {144, 140, 143};
Circle(342) = {143, 140, 145}; Circle(343) = {145, 140, 142};

// Level 6 - Ledge 2 top
Circle(350) = {152, 150, 154}; Circle(351) = {154, 150, 153};
Circle(352) = {153, 150, 155}; Circle(353) = {155, 150, 152};

// Level 7 - Upper bulge (ASYMMETRIC - vertical interlayer notch)
Ellipse(360) = {162, 160, 162, 164};
Ellipse(361) = {164, 160, 164, 163};
Ellipse(362) = {163, 160, 163, 165};
Ellipse(363) = {165, 160, 165, 162};

//////////////////////////////
// CAVERN SURFACES
//////////////////////////////

Curve Loop(401) = {200, 300, -202}; Surface(401) = {401};
Curve Loop(402) = {202, 301, -201}; Surface(402) = {402};
Curve Loop(403) = {201, 302, -203}; Surface(403) = {403};
Curve Loop(404) = {203, 303, -200}; Surface(404) = {404};

Curve Loop(411) = {300, 212, -310, -210}; Surface(411) = {411};
Curve Loop(412) = {301, 211, -311, -212}; Surface(412) = {412};
Curve Loop(413) = {302, 213, -312, -211}; Surface(413) = {413};
Curve Loop(414) = {303, 210, -313, -213}; Surface(414) = {414};

Curve Loop(421) = {310, 222, -320, -220}; Surface(421) = {421};
Curve Loop(422) = {311, 221, -321, -222}; Surface(422) = {422};
Curve Loop(423) = {312, 223, -322, -221}; Surface(423) = {423};
Curve Loop(424) = {313, 220, -323, -223}; Surface(424) = {424};

Curve Loop(431) = {320, 232, -330, -230}; Surface(431) = {431};
Curve Loop(432) = {321, 231, -331, -232}; Surface(432) = {432};
Curve Loop(433) = {322, 233, -332, -231}; Surface(433) = {433};
Curve Loop(434) = {323, 230, -333, -233}; Surface(434) = {434};

Curve Loop(441) = {330, 242, -340, -240}; Surface(441) = {441};
Curve Loop(442) = {331, 241, -341, -242}; Surface(442) = {442};
Curve Loop(443) = {332, 243, -342, -241}; Surface(443) = {443};
Curve Loop(444) = {333, 240, -343, -243}; Surface(444) = {444};

Curve Loop(451) = {340, 252, -350, -250}; Surface(451) = {451};
Curve Loop(452) = {341, 251, -351, -252}; Surface(452) = {452};
Curve Loop(453) = {342, 253, -352, -251}; Surface(453) = {453};
Curve Loop(454) = {343, 250, -353, -253}; Surface(454) = {454};

Curve Loop(461) = {350, 262, -360, -260}; Surface(461) = {461};
Curve Loop(462) = {351, 261, -361, -262}; Surface(462) = {462};
Curve Loop(463) = {352, 263, -362, -261}; Surface(463) = {463};
Curve Loop(464) = {353, 260, -363, -263}; Surface(464) = {464};

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
