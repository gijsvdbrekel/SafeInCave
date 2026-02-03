//////////////////////////////////////////////////////////////
// CAVERN A5 ZUIDWENDING - FULL 3D WITH INTERLAYERS
// Based on sonar survey October 2023
// Target volume: ~965,531 m³
// Roof: 1140 m depth, Bottom: 1510 m depth, Height: 370 m
//
// HETEROGENEOUS: 2 interlayers (0.2m thick each) at constriction points
// creating ledges from non-leachable layers in heterogeneous salt
//
// Structure (5 volumes):
// - Salt_bottom (z=0 to 269.9m)
// - Interlayer_1 (z=269.9m to 270.1m) - 0.2m thick
// - Salt_middle (z=270.1m to 369.9m)
// - Interlayer_2 (z=369.9m to 370.1m) - 0.2m thick
// - Salt_top (z=370.1m to 660m)
//////////////////////////////////////////////////////////////

coarse_size = 65;
fine_size = 4.5;

Lz = 660.0;
Ly = 450.0;
Lx = 450.0;

x_center = Lx/2;  // 225.0
y_center = Ly/2;  // 225.0

//////////////////////////////////////////////////////////////
// INTERLAYER Z-POSITIONS (0.2m thick each)
//////////////////////////////////////////////////////////////

z_inter1_bot = 269.9;   // Bottom of interlayer 1
z_inter1_top = 270.1;   // Top of interlayer 1 (0.2m thick)
z_inter2_bot = 369.9;   // Bottom of interlayer 2
z_inter2_top = 370.1;   // Top of interlayer 2 (0.2m thick)

//////////////////////////////////////////////////////////////
// CAVERN GEOMETRY PARAMETERS - A5 PROFILE
//////////////////////////////////////////////////////////////

// Vertical positions (model coordinates)
z_bot_tip = 145.0;     // Bottom tip of cavern
z_level1 = 175.0;      // Lower expansion
z_level2 = 220.0;      // Maximum lower bulge
// Interlayer 1 at z=270 (constriction point)
z_level4 = 320.0;      // Middle bulge
// Interlayer 2 at z=370 (constriction point)
z_level6 = 420.0;      // Upper body
z_level7 = 460.0;      // Neck transition
z_level8 = 490.0;      // Chimney
z_top_tip = 515.0;     // Top tip of cavern

// Radii at each level
R_level1 = 35.0;       // Lower expansion
R_level2 = 45.0;       // Maximum lower bulge
R_inter1 = 30.0;       // At interlayer 1 (constriction - vertical walls)
R_level4 = 38.0;       // Middle bulge
R_inter2 = 28.0;       // At interlayer 2 (constriction - vertical walls)
R_level6 = 32.0;       // Upper body
R_level7 = 20.0;       // Neck transition
R_level8 = 12.0;       // Chimney (narrow)

//////////////////////////////
// OUTER BOX CORNERS AT ALL Z-LEVELS
//////////////////////////////

// z=0 plane (bottom)
Point(1) = {0, 0, 0, coarse_size};
Point(2) = {Lx, 0, 0, coarse_size};
Point(3) = {Lx, Ly, 0, coarse_size};
Point(4) = {0, Ly, 0, coarse_size};

// z=z_inter1_bot plane
Point(11) = {0, 0, z_inter1_bot, coarse_size};
Point(12) = {Lx, 0, z_inter1_bot, coarse_size};
Point(13) = {Lx, Ly, z_inter1_bot, coarse_size};
Point(14) = {0, Ly, z_inter1_bot, coarse_size};

// z=z_inter1_top plane
Point(21) = {0, 0, z_inter1_top, coarse_size};
Point(22) = {Lx, 0, z_inter1_top, coarse_size};
Point(23) = {Lx, Ly, z_inter1_top, coarse_size};
Point(24) = {0, Ly, z_inter1_top, coarse_size};

// z=z_inter2_bot plane
Point(31) = {0, 0, z_inter2_bot, coarse_size};
Point(32) = {Lx, 0, z_inter2_bot, coarse_size};
Point(33) = {Lx, Ly, z_inter2_bot, coarse_size};
Point(34) = {0, Ly, z_inter2_bot, coarse_size};

// z=z_inter2_top plane
Point(41) = {0, 0, z_inter2_top, coarse_size};
Point(42) = {Lx, 0, z_inter2_top, coarse_size};
Point(43) = {Lx, Ly, z_inter2_top, coarse_size};
Point(44) = {0, Ly, z_inter2_top, coarse_size};

// z=Lz plane (top)
Point(51) = {0, 0, Lz, coarse_size};
Point(52) = {Lx, 0, Lz, coarse_size};
Point(53) = {Lx, Ly, Lz, coarse_size};
Point(54) = {0, Ly, Lz, coarse_size};

//////////////////////////////
// CAVERN POINTS
//////////////////////////////

// Level 0 - Bottom tip
Point(100) = {x_center, y_center, z_bot_tip, fine_size};

// Level 1 - Lower expansion (R = 35)
Point(110) = {x_center, y_center, z_level1, coarse_size};
Point(111) = {x_center + R_level1, y_center, z_level1, fine_size};
Point(112) = {x_center - R_level1, y_center, z_level1, fine_size};
Point(113) = {x_center, y_center + R_level1, z_level1, fine_size};
Point(114) = {x_center, y_center - R_level1, z_level1, fine_size};

// Level 2 - Maximum lower bulge (R = 45)
Point(120) = {x_center, y_center, z_level2, coarse_size};
Point(121) = {x_center + R_level2, y_center, z_level2, fine_size};
Point(122) = {x_center - R_level2, y_center, z_level2, fine_size};
Point(123) = {x_center, y_center + R_level2, z_level2, fine_size};
Point(124) = {x_center, y_center - R_level2, z_level2, fine_size};

// Level 3a - Bottom of interlayer 1 (R = 30)
Point(130) = {x_center, y_center, z_inter1_bot, coarse_size};
Point(131) = {x_center + R_inter1, y_center, z_inter1_bot, fine_size};
Point(132) = {x_center - R_inter1, y_center, z_inter1_bot, fine_size};
Point(133) = {x_center, y_center + R_inter1, z_inter1_bot, fine_size};
Point(134) = {x_center, y_center - R_inter1, z_inter1_bot, fine_size};

// Level 3b - Top of interlayer 1 (same R = 30, vertical wall)
Point(135) = {x_center, y_center, z_inter1_top, coarse_size};
Point(136) = {x_center + R_inter1, y_center, z_inter1_top, fine_size};
Point(137) = {x_center - R_inter1, y_center, z_inter1_top, fine_size};
Point(138) = {x_center, y_center + R_inter1, z_inter1_top, fine_size};
Point(139) = {x_center, y_center - R_inter1, z_inter1_top, fine_size};

// Level 4 - Middle bulge (R = 38)
Point(140) = {x_center, y_center, z_level4, coarse_size};
Point(141) = {x_center + R_level4, y_center, z_level4, fine_size};
Point(142) = {x_center - R_level4, y_center, z_level4, fine_size};
Point(143) = {x_center, y_center + R_level4, z_level4, fine_size};
Point(144) = {x_center, y_center - R_level4, z_level4, fine_size};

// Level 5a - Bottom of interlayer 2 (R = 28)
Point(150) = {x_center, y_center, z_inter2_bot, coarse_size};
Point(151) = {x_center + R_inter2, y_center, z_inter2_bot, fine_size};
Point(152) = {x_center - R_inter2, y_center, z_inter2_bot, fine_size};
Point(153) = {x_center, y_center + R_inter2, z_inter2_bot, fine_size};
Point(154) = {x_center, y_center - R_inter2, z_inter2_bot, fine_size};

// Level 5b - Top of interlayer 2 (same R = 28, vertical wall)
Point(155) = {x_center, y_center, z_inter2_top, coarse_size};
Point(156) = {x_center + R_inter2, y_center, z_inter2_top, fine_size};
Point(157) = {x_center - R_inter2, y_center, z_inter2_top, fine_size};
Point(158) = {x_center, y_center + R_inter2, z_inter2_top, fine_size};
Point(159) = {x_center, y_center - R_inter2, z_inter2_top, fine_size};

// Level 6 - Upper body (R = 32)
Point(160) = {x_center, y_center, z_level6, coarse_size};
Point(161) = {x_center + R_level6, y_center, z_level6, fine_size};
Point(162) = {x_center - R_level6, y_center, z_level6, fine_size};
Point(163) = {x_center, y_center + R_level6, z_level6, fine_size};
Point(164) = {x_center, y_center - R_level6, z_level6, fine_size};

// Level 7 - Neck transition (R = 20)
Point(170) = {x_center, y_center, z_level7, coarse_size};
Point(171) = {x_center + R_level7, y_center, z_level7, fine_size};
Point(172) = {x_center - R_level7, y_center, z_level7, fine_size};
Point(173) = {x_center, y_center + R_level7, z_level7, fine_size};
Point(174) = {x_center, y_center - R_level7, z_level7, fine_size};

// Level 8 - Chimney (R = 12)
Point(180) = {x_center, y_center, z_level8, coarse_size};
Point(181) = {x_center + R_level8, y_center, z_level8, fine_size};
Point(182) = {x_center - R_level8, y_center, z_level8, fine_size};
Point(183) = {x_center, y_center + R_level8, z_level8, fine_size};
Point(184) = {x_center, y_center - R_level8, z_level8, fine_size};

// Level 9 - Top tip
Point(190) = {x_center, y_center, z_top_tip, fine_size};

//////////////////////////////
// HORIZONTAL EDGES AT EACH Z-LEVEL
//////////////////////////////

// z=0 (bottom)
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// z=z_inter1_bot
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 11};

// z=z_inter1_top
Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 24};
Line(24) = {24, 21};

// z=z_inter2_bot
Line(31) = {31, 32};
Line(32) = {32, 33};
Line(33) = {33, 34};
Line(34) = {34, 31};

// z=z_inter2_top
Line(41) = {41, 42};
Line(42) = {42, 43};
Line(43) = {43, 44};
Line(44) = {44, 41};

// z=Lz (top)
Line(51) = {51, 52};
Line(52) = {52, 53};
Line(53) = {53, 54};
Line(54) = {54, 51};

//////////////////////////////
// VERTICAL EDGES (DOMAIN BOUNDARIES)
//////////////////////////////

// x=0, y=0 (SW corner)
Line(61) = {1, 11};
Line(62) = {11, 21};
Line(63) = {21, 31};
Line(64) = {31, 41};
Line(65) = {41, 51};

// x=Lx, y=0 (SE corner)
Line(71) = {2, 12};
Line(72) = {12, 22};
Line(73) = {22, 32};
Line(74) = {32, 42};
Line(75) = {42, 52};

// x=Lx, y=Ly (NE corner)
Line(81) = {3, 13};
Line(82) = {13, 23};
Line(83) = {23, 33};
Line(84) = {33, 43};
Line(85) = {43, 53};

// x=0, y=Ly (NW corner)
Line(91) = {4, 14};
Line(92) = {14, 24};
Line(93) = {24, 34};
Line(94) = {34, 44};
Line(95) = {44, 54};

//////////////////////////////
// CAVERN VERTICAL CONNECTIONS
//////////////////////////////

// Bottom tip to level 1
Line(200) = {100, 111};
Line(201) = {100, 112};
Line(202) = {100, 113};
Line(203) = {100, 114};

// Level 1 to level 2
Line(210) = {111, 121};
Line(211) = {112, 122};
Line(212) = {113, 123};
Line(213) = {114, 124};

// Level 2 to interlayer 1 bottom
Line(220) = {121, 131};
Line(221) = {122, 132};
Line(222) = {123, 133};
Line(223) = {124, 134};

// Through interlayer 1 (VERTICAL - same radius)
Line(225) = {131, 136};
Line(226) = {132, 137};
Line(227) = {133, 138};
Line(228) = {134, 139};

// Interlayer 1 top to level 4
Line(230) = {136, 141};
Line(231) = {137, 142};
Line(232) = {138, 143};
Line(233) = {139, 144};

// Level 4 to interlayer 2 bottom
Line(240) = {141, 151};
Line(241) = {142, 152};
Line(242) = {143, 153};
Line(243) = {144, 154};

// Through interlayer 2 (VERTICAL - same radius)
Line(245) = {151, 156};
Line(246) = {152, 157};
Line(247) = {153, 158};
Line(248) = {154, 159};

// Interlayer 2 top to level 6
Line(250) = {156, 161};
Line(251) = {157, 162};
Line(252) = {158, 163};
Line(253) = {159, 164};

// Level 6 to level 7
Line(260) = {161, 171};
Line(261) = {162, 172};
Line(262) = {163, 173};
Line(263) = {164, 174};

// Level 7 to level 8
Line(270) = {171, 181};
Line(271) = {172, 182};
Line(272) = {173, 183};
Line(273) = {174, 184};

// Level 8 to top tip
Line(280) = {181, 190};
Line(281) = {182, 190};
Line(282) = {183, 190};
Line(283) = {184, 190};

//////////////////////////////
// CIRCLE ARCS AT EACH CAVERN LEVEL
//////////////////////////////

// Level 1 (R = 35)
Circle(300) = {111, 110, 113};
Circle(301) = {113, 110, 112};
Circle(302) = {112, 110, 114};
Circle(303) = {114, 110, 111};

// Level 2 (R = 45)
Circle(310) = {121, 120, 123};
Circle(311) = {123, 120, 122};
Circle(312) = {122, 120, 124};
Circle(313) = {124, 120, 121};

// Interlayer 1 bottom (R = 30)
Circle(320) = {131, 130, 133};
Circle(321) = {133, 130, 132};
Circle(322) = {132, 130, 134};
Circle(323) = {134, 130, 131};

// Interlayer 1 top (R = 30)
Circle(325) = {136, 135, 138};
Circle(326) = {138, 135, 137};
Circle(327) = {137, 135, 139};
Circle(328) = {139, 135, 136};

// Level 4 (R = 38)
Circle(330) = {141, 140, 143};
Circle(331) = {143, 140, 142};
Circle(332) = {142, 140, 144};
Circle(333) = {144, 140, 141};

// Interlayer 2 bottom (R = 28)
Circle(340) = {151, 150, 153};
Circle(341) = {153, 150, 152};
Circle(342) = {152, 150, 154};
Circle(343) = {154, 150, 151};

// Interlayer 2 top (R = 28)
Circle(345) = {156, 155, 158};
Circle(346) = {158, 155, 157};
Circle(347) = {157, 155, 159};
Circle(348) = {159, 155, 156};

// Level 6 (R = 32)
Circle(350) = {161, 160, 163};
Circle(351) = {163, 160, 162};
Circle(352) = {162, 160, 164};
Circle(353) = {164, 160, 161};

// Level 7 (R = 20)
Circle(360) = {171, 170, 173};
Circle(361) = {173, 170, 172};
Circle(362) = {172, 170, 174};
Circle(363) = {174, 170, 171};

// Level 8 (R = 12)
Circle(370) = {181, 180, 183};
Circle(371) = {183, 180, 182};
Circle(372) = {182, 180, 184};
Circle(373) = {184, 180, 181};

//////////////////////////////
// CAVERN SURFACES
//////////////////////////////

// Bottom cap (4 triangular surfaces)
Curve Loop(401) = {200, 300, -202};
Surface(401) = {401};
Curve Loop(402) = {202, 301, -201};
Surface(402) = {402};
Curve Loop(403) = {201, 302, -203};
Surface(403) = {403};
Curve Loop(404) = {203, 303, -200};
Surface(404) = {404};

// Side walls Level 1-2 (4 surfaces)
Curve Loop(411) = {300, 212, -310, -210};
Surface(411) = {411};
Curve Loop(412) = {301, 211, -311, -212};
Surface(412) = {412};
Curve Loop(413) = {302, 213, -312, -211};
Surface(413) = {413};
Curve Loop(414) = {303, 210, -313, -213};
Surface(414) = {414};

// Side walls Level 2 to interlayer 1 bottom (4 surfaces)
Curve Loop(421) = {310, 222, -320, -220};
Surface(421) = {421};
Curve Loop(422) = {311, 221, -321, -222};
Surface(422) = {422};
Curve Loop(423) = {312, 223, -322, -221};
Surface(423) = {423};
Curve Loop(424) = {313, 220, -323, -223};
Surface(424) = {424};

// Side walls through interlayer 1 - VERTICAL (4 surfaces)
Curve Loop(425) = {320, 227, -325, -225};
Surface(425) = {425};
Curve Loop(426) = {321, 226, -326, -227};
Surface(426) = {426};
Curve Loop(427) = {322, 228, -327, -226};
Surface(427) = {427};
Curve Loop(428) = {323, 225, -328, -228};
Surface(428) = {428};

// Side walls interlayer 1 top to level 4 (4 surfaces)
Curve Loop(431) = {325, 232, -330, -230};
Surface(431) = {431};
Curve Loop(432) = {326, 231, -331, -232};
Surface(432) = {432};
Curve Loop(433) = {327, 233, -332, -231};
Surface(433) = {433};
Curve Loop(434) = {328, 230, -333, -233};
Surface(434) = {434};

// Side walls Level 4 to interlayer 2 bottom (4 surfaces)
Curve Loop(441) = {330, 242, -340, -240};
Surface(441) = {441};
Curve Loop(442) = {331, 241, -341, -242};
Surface(442) = {442};
Curve Loop(443) = {332, 243, -342, -241};
Surface(443) = {443};
Curve Loop(444) = {333, 240, -343, -243};
Surface(444) = {444};

// Side walls through interlayer 2 - VERTICAL (4 surfaces)
Curve Loop(445) = {340, 247, -345, -245};
Surface(445) = {445};
Curve Loop(446) = {341, 246, -346, -247};
Surface(446) = {446};
Curve Loop(447) = {342, 248, -347, -246};
Surface(447) = {447};
Curve Loop(448) = {343, 245, -348, -248};
Surface(448) = {448};

// Side walls interlayer 2 top to level 6 (4 surfaces)
Curve Loop(451) = {345, 252, -350, -250};
Surface(451) = {451};
Curve Loop(452) = {346, 251, -351, -252};
Surface(452) = {452};
Curve Loop(453) = {347, 253, -352, -251};
Surface(453) = {453};
Curve Loop(454) = {348, 250, -353, -253};
Surface(454) = {454};

// Side walls Level 6-7 (4 surfaces)
Curve Loop(461) = {350, 262, -360, -260};
Surface(461) = {461};
Curve Loop(462) = {351, 261, -361, -262};
Surface(462) = {462};
Curve Loop(463) = {352, 263, -362, -261};
Surface(463) = {463};
Curve Loop(464) = {353, 260, -363, -263};
Surface(464) = {464};

// Side walls Level 7-8 (4 surfaces)
Curve Loop(471) = {360, 272, -370, -270};
Surface(471) = {471};
Curve Loop(472) = {361, 271, -371, -272};
Surface(472) = {472};
Curve Loop(473) = {362, 273, -372, -271};
Surface(473) = {473};
Curve Loop(474) = {363, 270, -373, -273};
Surface(474) = {474};

// Top cap (4 triangular surfaces)
Curve Loop(481) = {370, 282, -280};
Surface(481) = {481};
Curve Loop(482) = {371, 281, -282};
Surface(482) = {482};
Curve Loop(483) = {372, 283, -281};
Surface(483) = {483};
Curve Loop(484) = {373, 280, -283};
Surface(484) = {484};

//////////////////////////////
// DOMAIN BOUNDARY SURFACES (EXTERNAL)
//////////////////////////////

// Bottom (z=0)
Curve Loop(501) = {1, 2, 3, 4};
Plane Surface(501) = {501};

// Top (z=Lz)
Curve Loop(502) = {51, 52, 53, 54};
Plane Surface(502) = {502};

// South faces (y=0) - 5 segments
Curve Loop(511) = {1, 71, -11, -61};
Plane Surface(511) = {511};
Curve Loop(512) = {11, 72, -21, -62};
Plane Surface(512) = {512};
Curve Loop(513) = {21, 73, -31, -63};
Plane Surface(513) = {513};
Curve Loop(514) = {31, 74, -41, -64};
Plane Surface(514) = {514};
Curve Loop(515) = {41, 75, -51, -65};
Plane Surface(515) = {515};

// North faces (y=Ly) - 5 segments
Curve Loop(521) = {-3, 81, 13, -91};
Plane Surface(521) = {521};
Curve Loop(522) = {-13, 82, 23, -92};
Plane Surface(522) = {522};
Curve Loop(523) = {-23, 83, 33, -93};
Plane Surface(523) = {523};
Curve Loop(524) = {-33, 84, 43, -94};
Plane Surface(524) = {524};
Curve Loop(525) = {-43, 85, 53, -95};
Plane Surface(525) = {525};

// West faces (x=0) - 5 segments
Curve Loop(531) = {-4, 91, 14, -61};
Plane Surface(531) = {531};
Curve Loop(532) = {-14, 92, 24, -62};
Plane Surface(532) = {532};
Curve Loop(533) = {-24, 93, 34, -63};
Plane Surface(533) = {533};
Curve Loop(534) = {-34, 94, 44, -64};
Plane Surface(534) = {534};
Curve Loop(535) = {-44, 95, 54, -65};
Plane Surface(535) = {535};

// East faces (x=Lx) - 5 segments
Curve Loop(541) = {2, 81, -12, -71};
Plane Surface(541) = {541};
Curve Loop(542) = {12, 82, -22, -72};
Plane Surface(542) = {542};
Curve Loop(543) = {22, 83, -32, -73};
Plane Surface(543) = {543};
Curve Loop(544) = {32, 84, -42, -74};
Plane Surface(544) = {544};
Curve Loop(545) = {42, 85, -52, -75};
Plane Surface(545) = {545};

//////////////////////////////
// INTERNAL HORIZONTAL SURFACES (with cavern holes)
//////////////////////////////

// z=z_inter1_bot (bottom of interlayer 1) - with hole for cavern
Curve Loop(551) = {11, 12, 13, 14};
Curve Loop(552) = {320, 321, 322, 323};
Plane Surface(551) = {551, 552};

// z=z_inter1_top (top of interlayer 1) - with hole for cavern
Curve Loop(561) = {21, 22, 23, 24};
Curve Loop(562) = {325, 326, 327, 328};
Plane Surface(561) = {561, 562};

// z=z_inter2_bot (bottom of interlayer 2) - with hole for cavern
Curve Loop(571) = {31, 32, 33, 34};
Curve Loop(572) = {340, 341, 342, 343};
Plane Surface(571) = {571, 572};

// z=z_inter2_top (top of interlayer 2) - with hole for cavern
Curve Loop(581) = {41, 42, 43, 44};
Curve Loop(582) = {345, 346, 347, 348};
Plane Surface(581) = {581, 582};

//////////////////////////////
// VOLUMES
//////////////////////////////

// Salt_bottom (z=0 to z_inter1_bot) - contains bottom of cavern
Surface Loop(601) = {501, 511, 521, 531, 541, 551,
                     401, 402, 403, 404, 411, 412, 413, 414, 421, 422, 423, 424};
Volume(1) = {601};

// Interlayer_1 (z_inter1_bot to z_inter1_top) - cavern passes through
Surface Loop(602) = {551, 512, 522, 532, 542, 561,
                     425, 426, 427, 428};
Volume(2) = {602};

// Salt_middle (z_inter1_top to z_inter2_bot) - contains middle bulge
Surface Loop(603) = {561, 513, 523, 533, 543, 571,
                     431, 432, 433, 434, 441, 442, 443, 444};
Volume(3) = {603};

// Interlayer_2 (z_inter2_bot to z_inter2_top) - cavern passes through
Surface Loop(604) = {571, 514, 524, 534, 544, 581,
                     445, 446, 447, 448};
Volume(4) = {604};

// Salt_top (z_inter2_top to Lz) - contains top of cavern
Surface Loop(605) = {581, 515, 525, 535, 545, 502,
                     451, 452, 453, 454, 461, 462, 463, 464, 471, 472, 473, 474, 481, 482, 483, 484};
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
                                  425, 426, 427, 428,
                                  431, 432, 433, 434,
                                  441, 442, 443, 444,
                                  445, 446, 447, 448,
                                  451, 452, 453, 454,
                                  461, 462, 463, 464,
                                  471, 472, 473, 474,
                                  481, 482, 483, 484};

Physical Volume("Salt_bottom", 31) = {1};
Physical Volume("Interlayer_1", 32) = {2};
Physical Volume("Salt_middle", 33) = {3};
Physical Volume("Interlayer_2", 34) = {4};
Physical Volume("Salt_top", 35) = {5};

// Wall profile for visualization
Physical Curve("Wall_profile", 30) = {200, 210, 220, 225, 230, 240, 245, 250, 260, 270, 280};

//////////////////////////////
// MESH SIZE FIELDS FOR INTERLAYER REFINEMENT
// Fine mesh near cavern wall in interlayers, transitioning to coarse
//////////////////////////////

// Mesh sizes for interlayer refinement
interlayer_fine = 0.05;    // Very fine at interlayer cavern wall (0.05m = 5cm)
interlayer_trans = 2.0;    // Transition size
dist_min = 0;              // Start fine mesh at wall
dist_max = 50;             // Transition complete at 50m from wall

// Field 1: Distance from interlayer 1 cavern wall surfaces (425, 426, 427, 428)
Field[1] = Distance;
Field[1].SurfacesList = {425, 426, 427, 428};

// Field 2: Distance from interlayer 2 cavern wall surfaces (445, 446, 447, 448)
Field[2] = Distance;
Field[2].SurfacesList = {445, 446, 447, 448};

// Field 3: Threshold for interlayer 1 - fine near wall, coarse far away
Field[3] = Threshold;
Field[3].InField = 1;
Field[3].SizeMin = interlayer_fine;
Field[3].SizeMax = coarse_size;
Field[3].DistMin = dist_min;
Field[3].DistMax = dist_max;

// Field 4: Threshold for interlayer 2 - fine near wall, coarse far away
Field[4] = Threshold;
Field[4].InField = 2;
Field[4].SizeMin = interlayer_fine;
Field[4].SizeMax = coarse_size;
Field[4].DistMin = dist_min;
Field[4].DistMax = dist_max;

// Field 5: Distance from entire cavern surface (for general refinement)
Field[5] = Distance;
Field[5].SurfacesList = {401, 402, 403, 404,
                         411, 412, 413, 414,
                         421, 422, 423, 424,
                         425, 426, 427, 428,
                         431, 432, 433, 434,
                         441, 442, 443, 444,
                         445, 446, 447, 448,
                         451, 452, 453, 454,
                         461, 462, 463, 464,
                         471, 472, 473, 474,
                         481, 482, 483, 484};

// Field 6: General cavern refinement (standard fine_size near wall)
Field[6] = Threshold;
Field[6].InField = 5;
Field[6].SizeMin = fine_size;
Field[6].SizeMax = coarse_size;
Field[6].DistMin = 0;
Field[6].DistMax = 100;

// Field 7: Combine all fields - take minimum (finest mesh where needed)
Field[7] = Min;
Field[7].FieldsList = {3, 4, 6};

// Apply the combined field as background mesh
Background Field = 7;
