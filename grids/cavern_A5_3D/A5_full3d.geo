//////////////////////////////////////////////////////////////
// CAVERN A5 ZUIDWENDING - FULL 3D
// Based on sonar survey October 2023
// Target volume: 965,531 m³
// Roof: 1140 m depth, Bottom: 1510 m depth, Height: 370 m
// Shape: Elongated irregular with narrow chimney at top
//////////////////////////////////////////////////////////////

coarse_size = 65;
fine_size = 4.5;

Lz = 660.0;
Ly = 450.0;
Lx = 450.0;

x_center = Lx/2;  // 225.0
y_center = Ly/2;  // 225.0

//////////////////////////////////////////////////////////////
// CAVERN GEOMETRY PARAMETERS - A5 PROFILE
// Derived from sonar image showing irregular elongated shape
//////////////////////////////////////////////////////////////

// Vertical positions (model coordinates)
z_bot_tip = 145.0;     // Bottom tip of cavern
z_level1 = 175.0;      // Lower expansion
z_level2 = 220.0;      // Maximum lower bulge
z_level3 = 270.0;      // First constriction
z_level4 = 320.0;      // Middle bulge
z_level5 = 370.0;      // Upper constriction
z_level6 = 420.0;      // Upper body
z_level7 = 460.0;      // Neck transition
z_level8 = 490.0;      // Chimney
z_top_tip = 515.0;     // Top tip of cavern

// Radii at each level (scaled for ~965k m³ volume)
R_level1 = 35.0;       // Lower expansion
R_level2 = 45.0;       // Maximum lower bulge
R_level3 = 30.0;       // First constriction
R_level4 = 38.0;       // Middle bulge
R_level5 = 28.0;       // Upper constriction
R_level6 = 32.0;       // Upper body
R_level7 = 20.0;       // Neck transition
R_level8 = 12.0;       // Chimney (narrow)

//////////////////////////////
// OUTER BOX CORNERS
//////////////////////////////

// z=0 plane (bottom)
Point(1) = {0, 0, 0, coarse_size};
Point(2) = {Lx, 0, 0, coarse_size};
Point(3) = {Lx, Ly, 0, coarse_size};
Point(4) = {0, Ly, 0, coarse_size};

// z=Lz plane (top)
Point(5) = {0, 0, Lz, coarse_size};
Point(6) = {Lx, 0, Lz, coarse_size};
Point(7) = {Lx, Ly, Lz, coarse_size};
Point(8) = {0, Ly, Lz, coarse_size};

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

// Level 3 - First constriction (R = 30)
Point(130) = {x_center, y_center, z_level3, coarse_size};
Point(131) = {x_center + R_level3, y_center, z_level3, fine_size};
Point(132) = {x_center - R_level3, y_center, z_level3, fine_size};
Point(133) = {x_center, y_center + R_level3, z_level3, fine_size};
Point(134) = {x_center, y_center - R_level3, z_level3, fine_size};

// Level 4 - Middle bulge (R = 38)
Point(140) = {x_center, y_center, z_level4, coarse_size};
Point(141) = {x_center + R_level4, y_center, z_level4, fine_size};
Point(142) = {x_center - R_level4, y_center, z_level4, fine_size};
Point(143) = {x_center, y_center + R_level4, z_level4, fine_size};
Point(144) = {x_center, y_center - R_level4, z_level4, fine_size};

// Level 5 - Upper constriction (R = 28)
Point(150) = {x_center, y_center, z_level5, coarse_size};
Point(151) = {x_center + R_level5, y_center, z_level5, fine_size};
Point(152) = {x_center - R_level5, y_center, z_level5, fine_size};
Point(153) = {x_center, y_center + R_level5, z_level5, fine_size};
Point(154) = {x_center, y_center - R_level5, z_level5, fine_size};

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
// OUTER BOX EDGES
//////////////////////////////

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

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

// Level 2 to level 3
Line(220) = {121, 131};
Line(221) = {122, 132};
Line(222) = {123, 133};
Line(223) = {124, 134};

// Level 3 to level 4
Line(230) = {131, 141};
Line(231) = {132, 142};
Line(232) = {133, 143};
Line(233) = {134, 144};

// Level 4 to level 5
Line(240) = {141, 151};
Line(241) = {142, 152};
Line(242) = {143, 153};
Line(243) = {144, 154};

// Level 5 to level 6
Line(250) = {151, 161};
Line(251) = {152, 162};
Line(252) = {153, 163};
Line(253) = {154, 164};

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

// Level 3 (R = 30)
Circle(320) = {131, 130, 133};
Circle(321) = {133, 130, 132};
Circle(322) = {132, 130, 134};
Circle(323) = {134, 130, 131};

// Level 4 (R = 38)
Circle(330) = {141, 140, 143};
Circle(331) = {143, 140, 142};
Circle(332) = {142, 140, 144};
Circle(333) = {144, 140, 141};

// Level 5 (R = 28)
Circle(340) = {151, 150, 153};
Circle(341) = {153, 150, 152};
Circle(342) = {152, 150, 154};
Circle(343) = {154, 150, 151};

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

// Side walls Level 2-3 (4 surfaces)
Curve Loop(421) = {310, 222, -320, -220};
Surface(421) = {421};
Curve Loop(422) = {311, 221, -321, -222};
Surface(422) = {422};
Curve Loop(423) = {312, 223, -322, -221};
Surface(423) = {423};
Curve Loop(424) = {313, 220, -323, -223};
Surface(424) = {424};

// Side walls Level 3-4 (4 surfaces)
Curve Loop(431) = {320, 232, -330, -230};
Surface(431) = {431};
Curve Loop(432) = {321, 231, -331, -232};
Surface(432) = {432};
Curve Loop(433) = {322, 233, -332, -231};
Surface(433) = {433};
Curve Loop(434) = {323, 230, -333, -233};
Surface(434) = {434};

// Side walls Level 4-5 (4 surfaces)
Curve Loop(441) = {330, 242, -340, -240};
Surface(441) = {441};
Curve Loop(442) = {331, 241, -341, -242};
Surface(442) = {442};
Curve Loop(443) = {332, 243, -342, -241};
Surface(443) = {443};
Curve Loop(444) = {333, 240, -343, -243};
Surface(444) = {444};

// Side walls Level 5-6 (4 surfaces)
Curve Loop(451) = {340, 252, -350, -250};
Surface(451) = {451};
Curve Loop(452) = {341, 251, -351, -252};
Surface(452) = {452};
Curve Loop(453) = {342, 253, -352, -251};
Surface(453) = {453};
Curve Loop(454) = {343, 250, -353, -253};
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
// OUTER BOX SURFACES
//////////////////////////////

// Bottom (z=0)
Curve Loop(501) = {1, 2, 3, 4};
Plane Surface(501) = {501};

// Top (z=Lz)
Curve Loop(502) = {5, 6, 7, 8};
Plane Surface(502) = {502};

// South (y=0)
Curve Loop(503) = {1, 10, -5, -9};
Plane Surface(503) = {503};

// North (y=Ly)
Curve Loop(504) = {3, 12, -7, -11};
Plane Surface(504) = {504};

// West (x=0)
Curve Loop(505) = {4, 9, -8, -12};
Plane Surface(505) = {505};

// East (x=Lx)
Curve Loop(506) = {2, 11, -6, -10};
Plane Surface(506) = {506};

//////////////////////////////
// VOLUME
//////////////////////////////

// Cavern surface loop (all cavern surfaces)
Surface Loop(600) = {401, 402, 403, 404,
                     411, 412, 413, 414,
                     421, 422, 423, 424,
                     431, 432, 433, 434,
                     441, 442, 443, 444,
                     451, 452, 453, 454,
                     461, 462, 463, 464,
                     471, 472, 473, 474,
                     481, 482, 483, 484};

// Salt volume (box minus cavern)
Surface Loop(601) = {501, 502, 503, 504, 505, 506};
Volume(1) = {601, 600};

//////////////////////////////
// PHYSICAL GROUPS
//////////////////////////////

Physical Surface("Bottom", 27) = {501};
Physical Surface("Top", 22) = {502};
Physical Surface("South", 23) = {503};
Physical Surface("North", 24) = {504};
Physical Surface("West", 26) = {505};
Physical Surface("East", 25) = {506};

Physical Surface("Cavern", 29) = {401, 402, 403, 404,
                                  411, 412, 413, 414,
                                  421, 422, 423, 424,
                                  431, 432, 433, 434,
                                  441, 442, 443, 444,
                                  451, 452, 453, 454,
                                  461, 462, 463, 464,
                                  471, 472, 473, 474,
                                  481, 482, 483, 484};

Physical Volume("Salt", 28) = {1};

// Wall profile for visualization
Physical Curve("Wall_profile", 30) = {200, 210, 220, 230, 240, 250, 260, 270, 280};
