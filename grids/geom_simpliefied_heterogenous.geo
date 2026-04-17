// Note: cavern height is 224.82 meters.
 
coarse_size_1 = 200;
coarse_size_2 = 35;
fine_size = 8.5;
 
Lz = 660.0;
Ly = 2000.0;
Lx = Ly;
 
ovb_thickness = 400;
salt_thickness = 1000;
hanging_wall = 400;
 
depth = ovb_thickness + hanging_wall;
depth_aux = 430 + depth;
 
// Cavern profile
Point(0) = {0.000000, 0.000000, 430.000000 - depth_aux, fine_size};
Point(1) = {0.000000, 0.000000, 205.189718 - depth_aux, fine_size};
Point(2) = {0.000000, 0.000000, 221.346389 - depth_aux, coarse_size_1};
Point(3) = {0.000000, 0.000000, 424.920441 - depth_aux, coarse_size_1};
Point(4) = {0.000000, 0.000000, 235.887393 - depth_aux, coarse_size_1};
Point(5) = {0.000000, 0.000000, 393.414933 - depth_aux, coarse_size_1};
Point(6) = {0.000000, 0.000000, 311.823745 - depth_aux, coarse_size_1};
Point(7) = {0.000000, 0.000000, 301.321909 - depth_aux, coarse_size_1};
Point(8) = {0.000000, 0.000000, 380.489596 - depth_aux, coarse_size_1};
Point(9) = {0.000000, 0.000000, 344.944920 - depth_aux, coarse_size_1};
Point(10) = {0.000000, 0.000000, 412.802938 - depth_aux, coarse_size_1};
Point(11) = {0.000000, 0.000000, 227.809058 - depth_aux, coarse_size_1};
Point(12) = {0.000000, 0.000000, 402.301102 - depth_aux, coarse_size_1};
Point(13) = {0.000000, 0.000000, 331.211750 - depth_aux, coarse_size_1};
Point(14) = {0.000000, 0.000000, 251.236230 - depth_aux, coarse_size_1};
Point(15) = {0.000000, 0.000000, 356.254590 - depth_aux, coarse_size_1};
Point(16) = {0.000000, 0.000000, 285.165239 - depth_aux, coarse_size_1};
Point(17) = {0.000000, 0.000000, 267.392901 - depth_aux, coarse_size_1};
Point(28) = {45.000000, 0.000000, 344.944920 - depth_aux, fine_size};
Point(29) = {68.597561, 0.000000, 251.236230 - depth_aux, fine_size};
Point(30) = {47.195122, 0.000000, 331.211750 - depth_aux, fine_size};
Point(31) = {72.987805, 0.000000, 285.165239 - depth_aux, fine_size};
Point(32) = {47.743902, 0.000000, 235.887393 - depth_aux, fine_size};
Point(33) = {47.743902, 0.000000, 380.489596 - depth_aux, fine_size};
Point(34) = {0.000000, 48.292683, 356.254590 - depth_aux, fine_size};
Point(35) = {48.292683, 0.000000, 356.254590 - depth_aux, fine_size};
Point(36) = {0.000000, 42.804878, 393.414933 - depth_aux, fine_size};
Point(37) = {0.000000, 68.597561, 251.236230 - depth_aux, fine_size};
Point(38) = {0.000000, 56.524390, 311.823745 - depth_aux, fine_size};
Point(39) = {0.000000, 74.634146, 267.392901 - depth_aux, fine_size};
Point(40) = {0.000000, 45.000000, 344.944920 - depth_aux, fine_size};
Point(41) = {0.000000, 72.987805, 285.165239 - depth_aux, fine_size};
Point(42) = {56.524390, 0.000000, 311.823745 - depth_aux, fine_size};
Point(43) = {7.134146, 0.000000, 221.346389 - depth_aux, fine_size};
Point(44) = {34.573171, 0.000000, 402.301102 - depth_aux, fine_size};
Point(45) = {0.000000, 19.756098, 412.802938 - depth_aux, fine_size};
Point(46) = {0.000000, 47.195122, 331.211750 - depth_aux, fine_size};
Point(47) = {57.621951, 0.000000, 301.321909 - depth_aux, fine_size};
Point(48) = {74.634146, 0.000000, 267.392901 - depth_aux, fine_size};
Point(49) = {0.000000, 7.134146, 221.346389 - depth_aux, fine_size};
Point(50) = {0.000000, 10.426829, 424.920441 - depth_aux, fine_size};
Point(51) = {19.756098, 0.000000, 412.802938 - depth_aux, fine_size};
Point(52) = {0.000000, 21.951220, 227.809058 - depth_aux, fine_size};
Point(53) = {0.000000, 47.743902, 235.887393 - depth_aux, fine_size};
Point(54) = {0.000000, 47.743902, 380.489596 - depth_aux, fine_size};
Point(55) = {10.426829, 0.000000, 424.920441 - depth_aux, fine_size};
Point(56) = {42.804878, 0.000000, 393.414933 - depth_aux, fine_size};
Point(57) = {21.951220, 0.000000, 227.809058 - depth_aux, fine_size};
Point(58) = {0.000000, 57.621951, 301.321909 - depth_aux, fine_size};
Point(59) = {0.000000, 34.573171, 402.301102 - depth_aux, fine_size};
 
// Overburden layer
Point(60) = {0.0, 0.0, 0.0, coarse_size_2};
Point(61) = {Lx, 0.0, 0.0, coarse_size_1};
Point(62) = {Lx, Ly, 0.0, coarse_size_1};
Point(63) = {0.0, Ly, 0.0, coarse_size_1};
Point(64) = {0.0, 0.0, -ovb_thickness, coarse_size_2};
Point(65) = {Lx, 0.0, -ovb_thickness, coarse_size_1};
Point(66) = {Lx, Ly, -ovb_thickness, coarse_size_1};
Point(67) = {0.0, Ly, -ovb_thickness, coarse_size_1};
 
 
// Salt layer
Point(68) = {0.0, 0.0, -ovb_thickness-salt_thickness, coarse_size_2};
Point(69) = {Lx, 0.0, -ovb_thickness-salt_thickness, coarse_size_1};
Point(70) = {Lx, Ly, -ovb_thickness-salt_thickness, coarse_size_1};
Point(71) = {0.0, Ly, -ovb_thickness-salt_thickness, coarse_size_1};
Point(72) = {Lx, 0.0, 301.321909 - depth_aux, coarse_size_1};
Point(73) = {Lx, 0.0, 311.823745 - depth_aux, coarse_size_1};
Point(74) = {Lx, Ly, 301.321909 - depth_aux, coarse_size_1};
Point(75) = {Lx, Ly, 311.823745 - depth_aux, coarse_size_1};
Point(76) = {0.0, Ly, 301.321909 - depth_aux, coarse_size_1};
Point(77) = {0.0, Ly, 311.823745 - depth_aux, coarse_size_1};
 
 
Line(1) = {1, 43};
Line(2) = {43, 57};
Line(3) = {57, 32};
Line(4) = {32, 29};
Line(5) = {29, 48};
Line(6) = {48, 31};
Line(7) = {31, 47};
Line(8) = {47, 42};
Line(9) = {42, 30};
Line(10) = {30, 28};
Line(11) = {28, 35};
Line(12) = {35, 33};
Line(13) = {33, 56};
Line(14) = {56, 44};
Line(15) = {44, 51};
Line(16) = {51, 55};
Line(17) = {55, 0};
Line(18) = {0, 50};
Line(19) = {50, 45};
Line(20) = {45, 59};
Line(21) = {59, 36};
Line(22) = {36, 54};
Line(23) = {54, 34};
Line(24) = {34, 40};
Line(25) = {40, 46};
Line(26) = {46, 38};
Line(27) = {38, 58};
Line(28) = {58, 41};
Line(29) = {41, 39};
Line(30) = {39, 37};
Line(31) = {37, 53};
Line(32) = {53, 52};
Line(33) = {52, 49};
Line(34) = {49, 1};
Line(35) = {1, 68};
Line(36) = {68, 69};
Line(37) = {69, 70};
Line(38) = {70, 71};
Line(39) = {67, 66};
Line(40) = {66, 65};
Line(41) = {65, 64};
Line(42) = {64, 67};
Line(43) = {67, 63};
Line(44) = {63, 60};
Line(45) = {60, 61};
Line(46) = {61, 62};
Line(47) = {62, 66};
Line(48) = {62, 63};
Line(49) = {60, 64};
Line(50) = {64, 0};
Line(51) = {61, 65};
Line(52) = {68, 71};
Line(53) = {69, 72};
Line(54) = {73, 72};
Line(55) = {65, 73};
Line(56) = {72, 74};
Line(57) = {75, 73};
Line(58) = {75, 74};
Line(59) = {70, 74};
Line(60) = {75, 66};
Line(61) = {74, 76};
Line(62) = {75, 77};
Line(63) = {71, 76};
Line(64) = {76, 77};
Line(65) = {77, 67};
Line(66) = {42, 73};
Line(67) = {72, 47};
Line(68) = {38, 77};
Line(69) = {76, 58};
//+
Circle(70) = {49, 2, 43};
//+
Circle(71) = {52, 11, 57};
//+
Circle(72) = {53, 4, 32};
//+
Circle(73) = {37, 14, 29};
//+
Circle(74) = {39, 17, 48};
//+
Circle(75) = {41, 16, 31};
//+
Circle(76) = {58, 7, 47};
//+
Circle(77) = {38, 6, 42};
//+
Circle(78) = {46, 13, 30};
//+
Circle(79) = {40, 9, 28};
//+
Circle(80) = {34, 15, 35};
//+
Circle(81) = {54, 8, 33};
//+
Circle(82) = {36, 5, 56};
//+
Circle(83) = {59, 12, 44};
//+
Circle(84) = {45, 10, 51};
//+
Circle(85) = {50, 3, 55};
//+
Curve Loop(1) = {52, -38, -37, -36};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {44, 45, 46, 48};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {43, -48, 47, -39};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {60, -39, -65, -62};
//+
Plane Surface(4) = {-4};
//+
Curve Loop(5) = {62, -64, -61, -58};
//+
Plane Surface(5) = {-5};
//+
Curve Loop(6) = {63, -61, -59, 38};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {53, 56, -59, -37};
//+
Plane Surface(7) = {-7};
//+
Curve Loop(8) = {54, 56, -58, 57};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {60, 40, 55, -57};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {51, -40, -47, -46};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {41, -49, 45, 51};
//+
Plane Surface(11) = {-11};
//+
Curve Loop(12) = {44, 49, 42, 43};
//+
Plane Surface(12) = {-12};
//+
Curve Loop(13) = {52, 63, 69, 28, 29, 30, 31, 32, 33, 34, 35};
//+
Plane Surface(13) = {-13};
//+
Curve Loop(14) = {69, -27, 68, -64};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {68, 65, -42, 50, 18, 19, 20, 21, 22, 23, 24, 25, 26};
//+
Plane Surface(15) = {-15};
//+
Curve Loop(16) = {35, 36, 53, 67, -7, -6, -5, -4, -3, -2, -1};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {67, 8, 66, 54};
//+
Plane Surface(17) = {-17};
//+
Curve Loop(18) = {66, -55, 41, 50, -17, -16, -15, -14, -13, -12, -11, -10, -9};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {42, 39, 40, 41};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {68, -62, 57, -66, -77};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {69, 76, -67, 56, 61};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {34, 1, -70};
//+
Surface(22) = {22};
//+
Curve Loop(23) = {33, 70, 2, -71};
//+
Surface(23) = {23};
//+
Curve Loop(24) = {32, 71, 3, -72};
//+
Surface(24) = {24};
//+
Curve Loop(25) = {73, -4, -72, -31};
//+
Surface(25) = {-25};
//+
Curve Loop(26) = {30, 73, 5, -74};
//+
Surface(26) = {26};
//+
Curve Loop(27) = {6, -75, 29, 74};
//+
Surface(27) = {27};
//+
Curve Loop(28) = {7, -76, 28, 75};
//+
Surface(28) = {28};
//+
Curve Loop(29) = {8, -77, 27, 76};
//+
Surface(29) = {29};
//+
Curve Loop(30) = {77, 9, -78, 26};
//+
Surface(30) = {30};
//+
Curve Loop(31) = {10, -79, 25, 78};
//+
Surface(31) = {31};
//+
Curve Loop(32) = {80, -11, -79, -24};
//+
Surface(32) = {-32};
//+
Curve Loop(33) = {12, -81, 23, 80};
//+
Surface(33) = {33};
//+
Curve Loop(34) = {13, -82, 22, 81};
//+
Surface(34) = {34};
//+
Curve Loop(35) = {14, -83, 21, 82};
//+
Surface(35) = {35};
//+
Curve Loop(36) = {15, -84, 20, 83};
//+
Surface(36) = {36};
//+
Curve Loop(37) = {16, -85, 19, 84};
//+
Surface(37) = {37};
//+
Curve Loop(38) = {17, 18, 85};
//+
Surface(38) = {38};
//+
Surface Loop(1) = {2, 12, 3, 10, 11, 19};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {18, 9, 4, 15, 30, 31, 32, 33, 34, 35, 36, 37, 38, 20, 19};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {14, 29, 17, 8, 5, 20, 21};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {1, 13, 16, 7, 6, 28, 27, 26, 25, 24, 23, 22, 21};
//+
Volume(4) = {4};
//+
Physical Surface("Top", 86) = {2};
//+
Physical Surface("Bottom", 87) = {1};
//+
Physical Surface("West", 88) = {15, 12, 13, 14};
//+
Physical Surface("South", 89) = {11, 18, 16, 17};
//+
Physical Surface("East", 90) = {10, 9, 8, 7};
//+
Physical Surface("North", 91) = {3, 4, 6, 5};
//+
Physical Surface("Cavern", 92) = {22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38};
//+
Physical Volume("Overburden", 93) = {1};
//+
Physical Volume("Salt_top", 94) = {2};
//+
Physical Volume("Interlayer", 95) = {3};
//+
Physical Volume("Salt_bottom", 96) = {4};
//+
Physical Curve("Wall_profile", 97) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
