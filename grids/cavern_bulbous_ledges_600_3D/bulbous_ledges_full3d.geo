//////////////////////////////////////////////////////////////
// BULBOUS LEDGES CAVERN - FULL 3D HETEROGENEOUS
// Target volume: 600,000 m³
//
// Design: Multiple stacked bulges with 3 thin horizontal interlayers
// that create sharp inward ledges (non-leachable mudstone/anhydrite)
// The cavern has an irregular profile with varying bulge sizes
// resembling real solution-mined caverns in heterogeneous salt
//
// HETEROGENEOUS: 4 salt volumes separated by 3 thin interlayers
// - Salt_bottom (z=0 to 175m)
// - Interlayer_1 (z=175m to 180m) - 5m thick mudstone
// - Salt_lower_mid (z=180m to 235m)
// - Interlayer_2 (z=235m to 240m) - 5m thick anhydrite
// - Salt_upper_mid (z=240m to 300m)
// - Interlayer_3 (z=300m to 305m) - 5m thick mudstone
// - Salt_top (z=305m to 660m)
//////////////////////////////////////////////////////////////

coarse_size = 65;
fine_size = 4.5;

Lz = 660.0;
Ly = 450.0;
Lx = 450.0;

x_center = Lx/2;  // 225.0
y_center = Ly/2;  // 225.0

//////////////////////////////////////////////////////////////
// INTERLAYER Z-POSITIONS (thin layers - 5m each)
//////////////////////////////////////////////////////////////

z_inter1_bot = 177.0;
z_inter1_top = 178.0;   // 1m thick

z_inter2_bot = 237.0;
z_inter2_top = 238.0;   // 1m thick

z_inter3_bot = 302.0;
z_inter3_top = 303.0;   // 1m thick

//////////////////////////////////////////////////////////////
// CAVERN GEOMETRY - Irregular bulbous profile
// 4 bulging zones with varying radii
//////////////////////////////////////////////////////////////

// Bottom bulge zone (below interlayer 1)
z_bot_tip = 145.0;
z_bot_max = 160.0;      // Maximum radius point
R_bottom = 38.5;        // Large bottom bulge

// At interlayer 1 - sharp inward ledge
R_ledge1 = 22.0;        // Smaller radius creates ledge

// Lower-middle bulge zone (between interlayers 1 and 2)
z_lmid_max = 208.0;
R_lower_mid = 48.0;     // Largest bulge

// At interlayer 2 - another sharp ledge
R_ledge2 = 25.0;

// Upper-middle bulge zone (between interlayers 2 and 3)
z_umid_max = 270.0;
R_upper_mid = 42.0;

// At interlayer 3 - top ledge
R_ledge3 = 20.0;

// Top bulge zone (above interlayer 3)
z_top_max = 325.0;
R_top = 28.0;
z_top_tip = 355.0;

//////////////////////////////
// OUTER BOX CORNERS AT EACH Z-LEVEL
//////////////////////////////

// z=0 plane (bottom)
Point(1) = {0, 0, 0, coarse_size};
Point(2) = {Lx, 0, 0, coarse_size};
Point(3) = {Lx, Ly, 0, coarse_size};
Point(4) = {0, Ly, 0, coarse_size};

// z=z_inter1_bot
Point(11) = {0, 0, z_inter1_bot, coarse_size};
Point(12) = {Lx, 0, z_inter1_bot, coarse_size};
Point(13) = {Lx, Ly, z_inter1_bot, coarse_size};
Point(14) = {0, Ly, z_inter1_bot, coarse_size};

// z=z_inter1_top
Point(21) = {0, 0, z_inter1_top, coarse_size};
Point(22) = {Lx, 0, z_inter1_top, coarse_size};
Point(23) = {Lx, Ly, z_inter1_top, coarse_size};
Point(24) = {0, Ly, z_inter1_top, coarse_size};

// z=z_inter2_bot
Point(31) = {0, 0, z_inter2_bot, coarse_size};
Point(32) = {Lx, 0, z_inter2_bot, coarse_size};
Point(33) = {Lx, Ly, z_inter2_bot, coarse_size};
Point(34) = {0, Ly, z_inter2_bot, coarse_size};

// z=z_inter2_top
Point(41) = {0, 0, z_inter2_top, coarse_size};
Point(42) = {Lx, 0, z_inter2_top, coarse_size};
Point(43) = {Lx, Ly, z_inter2_top, coarse_size};
Point(44) = {0, Ly, z_inter2_top, coarse_size};

// z=z_inter3_bot
Point(51) = {0, 0, z_inter3_bot, coarse_size};
Point(52) = {Lx, 0, z_inter3_bot, coarse_size};
Point(53) = {Lx, Ly, z_inter3_bot, coarse_size};
Point(54) = {0, Ly, z_inter3_bot, coarse_size};

// z=z_inter3_top
Point(61) = {0, 0, z_inter3_top, coarse_size};
Point(62) = {Lx, 0, z_inter3_top, coarse_size};
Point(63) = {Lx, Ly, z_inter3_top, coarse_size};
Point(64) = {0, Ly, z_inter3_top, coarse_size};

// z=Lz plane (top)
Point(71) = {0, 0, Lz, coarse_size};
Point(72) = {Lx, 0, Lz, coarse_size};
Point(73) = {Lx, Ly, Lz, coarse_size};
Point(74) = {0, Ly, Lz, coarse_size};

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

// Level 2 - Bottom of interlayer 1 (ledge)
Point(110) = {x_center, y_center, z_inter1_bot, coarse_size};
Point(112) = {x_center + R_ledge1, y_center, z_inter1_bot, fine_size};
Point(113) = {x_center - R_ledge1, y_center, z_inter1_bot, fine_size};
Point(114) = {x_center, y_center + R_ledge1, z_inter1_bot, fine_size};
Point(115) = {x_center, y_center - R_ledge1, z_inter1_bot, fine_size};

// Level 3 - Top of interlayer 1 (same radius - vertical wall)
Point(120) = {x_center, y_center, z_inter1_top, coarse_size};
Point(122) = {x_center + R_ledge1, y_center, z_inter1_top, fine_size};
Point(123) = {x_center - R_ledge1, y_center, z_inter1_top, fine_size};
Point(124) = {x_center, y_center + R_ledge1, z_inter1_top, fine_size};
Point(125) = {x_center, y_center - R_ledge1, z_inter1_top, fine_size};

// Level 4 - Lower-middle max bulge
Point(130) = {x_center, y_center, z_lmid_max, coarse_size};
Point(132) = {x_center + R_lower_mid, y_center, z_lmid_max, fine_size};
Point(133) = {x_center - R_lower_mid, y_center, z_lmid_max, fine_size};
Point(134) = {x_center, y_center + R_lower_mid, z_lmid_max, fine_size};
Point(135) = {x_center, y_center - R_lower_mid, z_lmid_max, fine_size};

// Level 5 - Bottom of interlayer 2 (ledge)
Point(140) = {x_center, y_center, z_inter2_bot, coarse_size};
Point(142) = {x_center + R_ledge2, y_center, z_inter2_bot, fine_size};
Point(143) = {x_center - R_ledge2, y_center, z_inter2_bot, fine_size};
Point(144) = {x_center, y_center + R_ledge2, z_inter2_bot, fine_size};
Point(145) = {x_center, y_center - R_ledge2, z_inter2_bot, fine_size};

// Level 6 - Top of interlayer 2
Point(150) = {x_center, y_center, z_inter2_top, coarse_size};
Point(152) = {x_center + R_ledge2, y_center, z_inter2_top, fine_size};
Point(153) = {x_center - R_ledge2, y_center, z_inter2_top, fine_size};
Point(154) = {x_center, y_center + R_ledge2, z_inter2_top, fine_size};
Point(155) = {x_center, y_center - R_ledge2, z_inter2_top, fine_size};

// Level 7 - Upper-middle max bulge
Point(160) = {x_center, y_center, z_umid_max, coarse_size};
Point(162) = {x_center + R_upper_mid, y_center, z_umid_max, fine_size};
Point(163) = {x_center - R_upper_mid, y_center, z_umid_max, fine_size};
Point(164) = {x_center, y_center + R_upper_mid, z_umid_max, fine_size};
Point(165) = {x_center, y_center - R_upper_mid, z_umid_max, fine_size};

// Level 8 - Bottom of interlayer 3 (ledge)
Point(170) = {x_center, y_center, z_inter3_bot, coarse_size};
Point(172) = {x_center + R_ledge3, y_center, z_inter3_bot, fine_size};
Point(173) = {x_center - R_ledge3, y_center, z_inter3_bot, fine_size};
Point(174) = {x_center, y_center + R_ledge3, z_inter3_bot, fine_size};
Point(175) = {x_center, y_center - R_ledge3, z_inter3_bot, fine_size};

// Level 9 - Top of interlayer 3
Point(180) = {x_center, y_center, z_inter3_top, coarse_size};
Point(182) = {x_center + R_ledge3, y_center, z_inter3_top, fine_size};
Point(183) = {x_center - R_ledge3, y_center, z_inter3_top, fine_size};
Point(184) = {x_center, y_center + R_ledge3, z_inter3_top, fine_size};
Point(185) = {x_center, y_center - R_ledge3, z_inter3_top, fine_size};

// Level 10 - Top max bulge
Point(190) = {x_center, y_center, z_top_max, coarse_size};
Point(192) = {x_center + R_top, y_center, z_top_max, fine_size};
Point(193) = {x_center - R_top, y_center, z_top_max, fine_size};
Point(194) = {x_center, y_center + R_top, z_top_max, fine_size};
Point(195) = {x_center, y_center - R_top, z_top_max, fine_size};

// Level 11 - Top tip
Point(200) = {x_center, y_center, z_top_tip, fine_size};

//////////////////////////////
// HORIZONTAL EDGES AT EACH Z-LEVEL
//////////////////////////////

// z=0
Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};

// z=z_inter1_bot
Line(11) = {11, 12}; Line(12) = {12, 13}; Line(13) = {13, 14}; Line(14) = {14, 11};

// z=z_inter1_top
Line(21) = {21, 22}; Line(22) = {22, 23}; Line(23) = {23, 24}; Line(24) = {24, 21};

// z=z_inter2_bot
Line(31) = {31, 32}; Line(32) = {32, 33}; Line(33) = {33, 34}; Line(34) = {34, 31};

// z=z_inter2_top
Line(41) = {41, 42}; Line(42) = {42, 43}; Line(43) = {43, 44}; Line(44) = {44, 41};

// z=z_inter3_bot
Line(51) = {51, 52}; Line(52) = {52, 53}; Line(53) = {53, 54}; Line(54) = {54, 51};

// z=z_inter3_top
Line(61) = {61, 62}; Line(62) = {62, 63}; Line(63) = {63, 64}; Line(64) = {64, 61};

// z=Lz
Line(71) = {71, 72}; Line(72) = {72, 73}; Line(73) = {73, 74}; Line(74) = {74, 71};

//////////////////////////////
// VERTICAL EDGES (DOMAIN BOUNDARIES)
//////////////////////////////

// SW corner (x=0, y=0)
Line(81) = {1, 11}; Line(82) = {11, 21}; Line(83) = {21, 31};
Line(84) = {31, 41}; Line(85) = {41, 51}; Line(86) = {51, 61}; Line(87) = {61, 71};

// SE corner (x=Lx, y=0)
Line(91) = {2, 12}; Line(92) = {12, 22}; Line(93) = {22, 32};
Line(94) = {32, 42}; Line(95) = {42, 52}; Line(96) = {52, 62}; Line(97) = {62, 72};

// NE corner (x=Lx, y=Ly)
Line(101) = {3, 13}; Line(102) = {13, 23}; Line(103) = {23, 33};
Line(104) = {33, 43}; Line(105) = {43, 53}; Line(106) = {53, 63}; Line(107) = {63, 73};

// NW corner (x=0, y=Ly)
Line(111) = {4, 14}; Line(112) = {14, 24}; Line(113) = {24, 34};
Line(114) = {34, 44}; Line(115) = {44, 54}; Line(116) = {54, 64}; Line(117) = {64, 74};

//////////////////////////////
// CAVERN VERTICAL CONNECTIONS
//////////////////////////////

// Bottom tip to level 1
Line(300) = {100, 102}; Line(301) = {100, 103}; Line(302) = {100, 104}; Line(303) = {100, 105};

// Level 1 to level 2
Line(310) = {102, 112}; Line(311) = {103, 113}; Line(312) = {104, 114}; Line(313) = {105, 115};

// Level 2 to level 3 (through interlayer 1)
Line(320) = {112, 122}; Line(321) = {113, 123}; Line(322) = {114, 124}; Line(323) = {115, 125};

// Level 3 to level 4
Line(330) = {122, 132}; Line(331) = {123, 133}; Line(332) = {124, 134}; Line(333) = {125, 135};

// Level 4 to level 5
Line(340) = {132, 142}; Line(341) = {133, 143}; Line(342) = {134, 144}; Line(343) = {135, 145};

// Level 5 to level 6 (through interlayer 2)
Line(350) = {142, 152}; Line(351) = {143, 153}; Line(352) = {144, 154}; Line(353) = {145, 155};

// Level 6 to level 7
Line(360) = {152, 162}; Line(361) = {153, 163}; Line(362) = {154, 164}; Line(363) = {155, 165};

// Level 7 to level 8
Line(370) = {162, 172}; Line(371) = {163, 173}; Line(372) = {164, 174}; Line(373) = {165, 175};

// Level 8 to level 9 (through interlayer 3)
Line(380) = {172, 182}; Line(381) = {173, 183}; Line(382) = {174, 184}; Line(383) = {175, 185};

// Level 9 to level 10
Line(390) = {182, 192}; Line(391) = {183, 193}; Line(392) = {184, 194}; Line(393) = {185, 195};

// Level 10 to top tip
Line(400) = {192, 200}; Line(401) = {193, 200}; Line(402) = {194, 200}; Line(403) = {195, 200};

//////////////////////////////
// CIRCLE ARCS AT EACH CAVERN LEVEL
//////////////////////////////

// Level 1 - Bottom bulge
Circle(410) = {102, 101, 104}; Circle(411) = {104, 101, 103};
Circle(412) = {103, 101, 105}; Circle(413) = {105, 101, 102};

// Level 2 - Ledge 1 bottom
Circle(420) = {112, 110, 114}; Circle(421) = {114, 110, 113};
Circle(422) = {113, 110, 115}; Circle(423) = {115, 110, 112};

// Level 3 - Ledge 1 top
Circle(430) = {122, 120, 124}; Circle(431) = {124, 120, 123};
Circle(432) = {123, 120, 125}; Circle(433) = {125, 120, 122};

// Level 4 - Lower-mid bulge
Circle(440) = {132, 130, 134}; Circle(441) = {134, 130, 133};
Circle(442) = {133, 130, 135}; Circle(443) = {135, 130, 132};

// Level 5 - Ledge 2 bottom
Circle(450) = {142, 140, 144}; Circle(451) = {144, 140, 143};
Circle(452) = {143, 140, 145}; Circle(453) = {145, 140, 142};

// Level 6 - Ledge 2 top
Circle(460) = {152, 150, 154}; Circle(461) = {154, 150, 153};
Circle(462) = {153, 150, 155}; Circle(463) = {155, 150, 152};

// Level 7 - Upper-mid bulge
Circle(470) = {162, 160, 164}; Circle(471) = {164, 160, 163};
Circle(472) = {163, 160, 165}; Circle(473) = {165, 160, 162};

// Level 8 - Ledge 3 bottom
Circle(480) = {172, 170, 174}; Circle(481) = {174, 170, 173};
Circle(482) = {173, 170, 175}; Circle(483) = {175, 170, 172};

// Level 9 - Ledge 3 top
Circle(490) = {182, 180, 184}; Circle(491) = {184, 180, 183};
Circle(492) = {183, 180, 185}; Circle(493) = {185, 180, 182};

// Level 10 - Top bulge
Circle(500) = {192, 190, 194}; Circle(501) = {194, 190, 193};
Circle(502) = {193, 190, 195}; Circle(503) = {195, 190, 192};

//////////////////////////////
// CAVERN SURFACES
//////////////////////////////

// Bottom cap
Curve Loop(601) = {300, 410, -302}; Surface(601) = {601};
Curve Loop(602) = {302, 411, -301}; Surface(602) = {602};
Curve Loop(603) = {301, 412, -303}; Surface(603) = {603};
Curve Loop(604) = {303, 413, -300}; Surface(604) = {604};

// Side walls L1-L2
Curve Loop(611) = {410, 312, -420, -310}; Surface(611) = {611};
Curve Loop(612) = {411, 311, -421, -312}; Surface(612) = {612};
Curve Loop(613) = {412, 313, -422, -311}; Surface(613) = {613};
Curve Loop(614) = {413, 310, -423, -313}; Surface(614) = {614};

// Side walls L2-L3 (interlayer 1)
Curve Loop(621) = {420, 322, -430, -320}; Surface(621) = {621};
Curve Loop(622) = {421, 321, -431, -322}; Surface(622) = {622};
Curve Loop(623) = {422, 323, -432, -321}; Surface(623) = {623};
Curve Loop(624) = {423, 320, -433, -323}; Surface(624) = {624};

// Side walls L3-L4
Curve Loop(631) = {430, 332, -440, -330}; Surface(631) = {631};
Curve Loop(632) = {431, 331, -441, -332}; Surface(632) = {632};
Curve Loop(633) = {432, 333, -442, -331}; Surface(633) = {633};
Curve Loop(634) = {433, 330, -443, -333}; Surface(634) = {634};

// Side walls L4-L5
Curve Loop(641) = {440, 342, -450, -340}; Surface(641) = {641};
Curve Loop(642) = {441, 341, -451, -342}; Surface(642) = {642};
Curve Loop(643) = {442, 343, -452, -341}; Surface(643) = {643};
Curve Loop(644) = {443, 340, -453, -343}; Surface(644) = {644};

// Side walls L5-L6 (interlayer 2)
Curve Loop(651) = {450, 352, -460, -350}; Surface(651) = {651};
Curve Loop(652) = {451, 351, -461, -352}; Surface(652) = {652};
Curve Loop(653) = {452, 353, -462, -351}; Surface(653) = {653};
Curve Loop(654) = {453, 350, -463, -353}; Surface(654) = {654};

// Side walls L6-L7
Curve Loop(661) = {460, 362, -470, -360}; Surface(661) = {661};
Curve Loop(662) = {461, 361, -471, -362}; Surface(662) = {662};
Curve Loop(663) = {462, 363, -472, -361}; Surface(663) = {663};
Curve Loop(664) = {463, 360, -473, -363}; Surface(664) = {664};

// Side walls L7-L8
Curve Loop(671) = {470, 372, -480, -370}; Surface(671) = {671};
Curve Loop(672) = {471, 371, -481, -372}; Surface(672) = {672};
Curve Loop(673) = {472, 373, -482, -371}; Surface(673) = {673};
Curve Loop(674) = {473, 370, -483, -373}; Surface(674) = {674};

// Side walls L8-L9 (interlayer 3)
Curve Loop(681) = {480, 382, -490, -380}; Surface(681) = {681};
Curve Loop(682) = {481, 381, -491, -382}; Surface(682) = {682};
Curve Loop(683) = {482, 383, -492, -381}; Surface(683) = {683};
Curve Loop(684) = {483, 380, -493, -383}; Surface(684) = {684};

// Side walls L9-L10
Curve Loop(691) = {490, 392, -500, -390}; Surface(691) = {691};
Curve Loop(692) = {491, 391, -501, -392}; Surface(692) = {692};
Curve Loop(693) = {492, 393, -502, -391}; Surface(693) = {693};
Curve Loop(694) = {493, 390, -503, -393}; Surface(694) = {694};

// Top cap
Curve Loop(701) = {500, 402, -400}; Surface(701) = {701};
Curve Loop(702) = {501, 401, -402}; Surface(702) = {702};
Curve Loop(703) = {502, 403, -401}; Surface(703) = {703};
Curve Loop(704) = {503, 400, -403}; Surface(704) = {704};

//////////////////////////////
// DOMAIN BOUNDARY SURFACES
//////////////////////////////

// Bottom (z=0)
Curve Loop(801) = {1, 2, 3, 4}; Plane Surface(801) = {801};

// Top (z=Lz)
Curve Loop(802) = {71, 72, 73, 74}; Plane Surface(802) = {802};

// South faces (y=0) - 7 segments
Curve Loop(811) = {1, 91, -11, -81}; Plane Surface(811) = {811};
Curve Loop(812) = {11, 92, -21, -82}; Plane Surface(812) = {812};
Curve Loop(813) = {21, 93, -31, -83}; Plane Surface(813) = {813};
Curve Loop(814) = {31, 94, -41, -84}; Plane Surface(814) = {814};
Curve Loop(815) = {41, 95, -51, -85}; Plane Surface(815) = {815};
Curve Loop(816) = {51, 96, -61, -86}; Plane Surface(816) = {816};
Curve Loop(817) = {61, 97, -71, -87}; Plane Surface(817) = {817};

// North faces (y=Ly)
Curve Loop(821) = {-3, 101, 13, -111}; Plane Surface(821) = {821};
Curve Loop(822) = {-13, 102, 23, -112}; Plane Surface(822) = {822};
Curve Loop(823) = {-23, 103, 33, -113}; Plane Surface(823) = {823};
Curve Loop(824) = {-33, 104, 43, -114}; Plane Surface(824) = {824};
Curve Loop(825) = {-43, 105, 53, -115}; Plane Surface(825) = {825};
Curve Loop(826) = {-53, 106, 63, -116}; Plane Surface(826) = {826};
Curve Loop(827) = {-63, 107, 73, -117}; Plane Surface(827) = {827};

// West faces (x=0)
Curve Loop(831) = {-4, 111, 14, -81}; Plane Surface(831) = {831};
Curve Loop(832) = {-14, 112, 24, -82}; Plane Surface(832) = {832};
Curve Loop(833) = {-24, 113, 34, -83}; Plane Surface(833) = {833};
Curve Loop(834) = {-34, 114, 44, -84}; Plane Surface(834) = {834};
Curve Loop(835) = {-44, 115, 54, -85}; Plane Surface(835) = {835};
Curve Loop(836) = {-54, 116, 64, -86}; Plane Surface(836) = {836};
Curve Loop(837) = {-64, 117, 74, -87}; Plane Surface(837) = {837};

// East faces (x=Lx)
Curve Loop(841) = {2, 101, -12, -91}; Plane Surface(841) = {841};
Curve Loop(842) = {12, 102, -22, -92}; Plane Surface(842) = {842};
Curve Loop(843) = {22, 103, -32, -93}; Plane Surface(843) = {843};
Curve Loop(844) = {32, 104, -42, -94}; Plane Surface(844) = {844};
Curve Loop(845) = {42, 105, -52, -95}; Plane Surface(845) = {845};
Curve Loop(846) = {52, 106, -62, -96}; Plane Surface(846) = {846};
Curve Loop(847) = {62, 107, -72, -97}; Plane Surface(847) = {847};

//////////////////////////////
// INTERNAL HORIZONTAL SURFACES (with cavern holes)
//////////////////////////////

// z=z_inter1_bot
Curve Loop(851) = {11, 12, 13, 14};
Curve Loop(852) = {420, 421, 422, 423};
Plane Surface(851) = {851, 852};

// z=z_inter1_top
Curve Loop(861) = {21, 22, 23, 24};
Curve Loop(862) = {430, 431, 432, 433};
Plane Surface(861) = {861, 862};

// z=z_inter2_bot
Curve Loop(871) = {31, 32, 33, 34};
Curve Loop(872) = {450, 451, 452, 453};
Plane Surface(871) = {871, 872};

// z=z_inter2_top
Curve Loop(881) = {41, 42, 43, 44};
Curve Loop(882) = {460, 461, 462, 463};
Plane Surface(881) = {881, 882};

// z=z_inter3_bot
Curve Loop(891) = {51, 52, 53, 54};
Curve Loop(892) = {480, 481, 482, 483};
Plane Surface(891) = {891, 892};

// z=z_inter3_top
Curve Loop(901) = {61, 62, 63, 64};
Curve Loop(902) = {490, 491, 492, 493};
Plane Surface(901) = {901, 902};

//////////////////////////////
// VOLUMES
//////////////////////////////

// Salt_bottom (z=0 to z_inter1_bot)
Surface Loop(1001) = {801, 811, 821, 831, 841, 851,
                      601, 602, 603, 604, 611, 612, 613, 614};
Volume(1) = {1001};

// Interlayer_1 (z_inter1_bot to z_inter1_top)
Surface Loop(1002) = {851, 812, 822, 832, 842, 861,
                      621, 622, 623, 624};
Volume(2) = {1002};

// Salt_lower_mid (z_inter1_top to z_inter2_bot)
Surface Loop(1003) = {861, 813, 823, 833, 843, 871,
                      631, 632, 633, 634, 641, 642, 643, 644};
Volume(3) = {1003};

// Interlayer_2 (z_inter2_bot to z_inter2_top)
Surface Loop(1004) = {871, 814, 824, 834, 844, 881,
                      651, 652, 653, 654};
Volume(4) = {1004};

// Salt_upper_mid (z_inter2_top to z_inter3_bot)
Surface Loop(1005) = {881, 815, 825, 835, 845, 891,
                      661, 662, 663, 664, 671, 672, 673, 674};
Volume(5) = {1005};

// Interlayer_3 (z_inter3_bot to z_inter3_top)
Surface Loop(1006) = {891, 816, 826, 836, 846, 901,
                      681, 682, 683, 684};
Volume(6) = {1006};

// Salt_top (z_inter3_top to Lz)
Surface Loop(1007) = {901, 817, 827, 837, 847, 802,
                      691, 692, 693, 694, 701, 702, 703, 704};
Volume(7) = {1007};

//////////////////////////////
// PHYSICAL GROUPS
//////////////////////////////

Physical Surface("Bottom", 27) = {801};
Physical Surface("Top", 22) = {802};
Physical Surface("South", 23) = {811, 812, 813, 814, 815, 816, 817};
Physical Surface("North", 24) = {821, 822, 823, 824, 825, 826, 827};
Physical Surface("West", 26) = {831, 832, 833, 834, 835, 836, 837};
Physical Surface("East", 25) = {841, 842, 843, 844, 845, 846, 847};

Physical Surface("Cavern", 29) = {601, 602, 603, 604,
                                  611, 612, 613, 614,
                                  621, 622, 623, 624,
                                  631, 632, 633, 634,
                                  641, 642, 643, 644,
                                  651, 652, 653, 654,
                                  661, 662, 663, 664,
                                  671, 672, 673, 674,
                                  681, 682, 683, 684,
                                  691, 692, 693, 694,
                                  701, 702, 703, 704};

Physical Volume("Salt_bottom", 31) = {1};
Physical Volume("Interlayer_1", 32) = {2};
Physical Volume("Salt_lower_mid", 33) = {3};
Physical Volume("Interlayer_2", 34) = {4};
Physical Volume("Salt_upper_mid", 35) = {5};
Physical Volume("Interlayer_3", 36) = {6};
Physical Volume("Salt_top", 37) = {7};

Physical Curve("Wall_profile", 30) = {300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400};
