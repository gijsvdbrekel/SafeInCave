SetFactory("OpenCASCADE");

// ------------------------------------------------------------
// Mesh sizes
// ------------------------------------------------------------
size_coarse = 65;
size_fine   = 4.5;

// ------------------------------------------------------------
// Domain box
// ------------------------------------------------------------
Lx = 450;
Ly = 450;
Lz = 660;

// salt block
Box(1) = {0, 0, 0, Lx, Ly, Lz};

// ------------------------------------------------------------
// Cylindrical cavern with a small angle
// ------------------------------------------------------------

// Cavern parameters
R_cav  = 54.42210251952307; // radius
L_cav  = 230;               // length of cavern axis (3D)
alpha  = 5*Pi/180;          // tilt angle from vertical (~5°)

// Start point of cavern axis (in the middle of the block)
xc0 = Lx/2;
yc0 = Ly/2;
zc0 = 380;      // depth of top of cavern axis

// Direction vector of the axis: a bit in +x, mostly downward in -z
dx =  L_cav*Sin(alpha);
dy =  0;
dz = -L_cav*Cos(alpha);

// Cylinder along tilted axis
// Cylinder(tag) = {x0, y0, z0, dx, dy, dz, R, angle};
Cylinder(2) = {xc0, yc0, zc0, dx, dy, dz, R_cav, 2*Pi};

// ------------------------------------------------------------
// Boolean difference: salt MINUS cavern
// ------------------------------------------------------------
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

// The result volume ID is 3
SaltVol = 3;

// ------------------------------------------------------------
// Helper line: cavern axis (wall profile)
// ------------------------------------------------------------
Point(100) = {xc0,        yc0,        zc0,        size_fine};
Point(101) = {xc0 + dx,   yc0 + dy,   zc0 + dz,   size_fine};
Line(200)  = {100, 101};

// ------------------------------------------------------------
// Surface Selection
// ------------------------------------------------------------

// Get all surfaces of the resulting volume
allSurfaces[] = Surface{:};

// Initialize arrays
tops[] = {};
bottoms[] = {};
souths[] = {};
norths[] = {};
wests[] = {};
easts[] = {};
cavs[] = {};

// Loop through all surfaces and classify by center-of-mass
For i In {0 : #allSurfaces[]-1}
    surf = allSurfaces[i];
    bbox[] = BoundingBox Surface{surf};
    
    // bbox format: {xmin, ymin, zmin, xmax, ymax, zmax}
    xcenter = (bbox[0] + bbox[3]) / 2;
    ycenter = (bbox[1] + bbox[4]) / 2;
    zcenter = (bbox[2] + bbox[5]) / 2;
    
    tol = 1.0;
    
    // Check if it's on one of the box faces
    If (Fabs(zcenter - Lz) < tol)
        tops[] += {surf};
    EndIf
    If (Fabs(zcenter - 0) < tol)
        bottoms[] += {surf};
    EndIf
    If (Fabs(ycenter - 0) < tol)
        souths[] += {surf};
    EndIf
    If (Fabs(ycenter - Ly) < tol)
        norths[] += {surf};
    EndIf
    If (Fabs(xcenter - 0) < tol)
        wests[] += {surf};
    EndIf
    If (Fabs(xcenter - Lx) < tol)
        easts[] += {surf};
    EndIf
    
    // If not on any box face, it's probably the cavern
    If (Fabs(zcenter - Lz) > tol && Fabs(zcenter - 0) > tol &&
        Fabs(ycenter - 0) > tol && Fabs(ycenter - Ly) > tol &&
        Fabs(xcenter - 0) > tol && Fabs(xcenter - Lx) > tol)
        cavs[] += {surf};
    EndIf
EndFor

// ------------------------------------------------------------
// Physical groups
// ------------------------------------------------------------
Physical Volume("Salt", 28) = {SaltVol};

If (#tops[] > 0)
    Physical Surface("Top", 22) = {tops[]};
EndIf
If (#bottoms[] > 0)
    Physical Surface("Bottom", 27) = {bottoms[]};
EndIf
If (#souths[] > 0)
    Physical Surface("South", 23) = {souths[]};
EndIf
If (#norths[] > 0)
    Physical Surface("North", 24) = {norths[]};
EndIf
If (#wests[] > 0)
    Physical Surface("West", 26) = {wests[]};
EndIf
If (#easts[] > 0)
    Physical Surface("East", 25) = {easts[]};
EndIf
If (#cavs[] > 0)
    Physical Surface("Cavern", 29) = {cavs[]};
EndIf

Physical Curve("Wall_profile", 30) = {200};

// ------------------------------------------------------------
// Mesh refinement ONLY near cavern wall (like original)
// ------------------------------------------------------------

// Refine near cavern axis (=> near cavern wall)
Field[1] = Distance;
Field[1].CurvesList = {200};   // axis line
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField  = 1;
Field[2].SizeMin  = size_fine;     // fine near cavern (4.5 m)
Field[2].SizeMax  = size_coarse;   // coarse far away (65 m)
Field[2].DistMin  = R_cav;         // up to radius ~54 m: fine
Field[2].DistMax  = 3*R_cav;       // gradual transition up to ~163 m

// Set this as the background field
Background Field = 2;

// Override mesh size constraints
Mesh.CharacteristicLengthMin = size_fine;
Mesh.CharacteristicLengthMax = size_coarse;
Mesh.CharacteristicLengthFromPoints     = 0;
Mesh.CharacteristicLengthFromCurvature  = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;