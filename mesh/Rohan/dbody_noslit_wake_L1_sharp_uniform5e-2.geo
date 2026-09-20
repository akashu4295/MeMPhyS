SetFactory("OpenCASCADE");

// =====================================================================
//  POINTS
// =====================================================================

// --- domain corners --------------------------------------------------
Point(1) = {-27.5, -10, 0, 1.0};
Point(2) = {-27.5, 10, 0, 1.0};
Point(3) = {27.5, -10, 0, 1.0};
Point(4) = {27.5, 10, 0, 1.0};

// --- body --------------------------------------------------------------
Point(5) = {-22.5, 0, 0, 1.0};       // nose arc CENTRE (construction only)
Point(6) = {-22.5, 0.5, 0, 1.0};     // top of semicircle
Point(7) = {-20.5, 0.5, 0, 1.0};     // TOP TRAILING CORNER (sharp)
Point(8) = {-20.5, -0.5, 0, 1.0};    // BOTTOM TRAILING CORNER (sharp)
Point(9) = {-22.5, -0.5, 0, 1.0};    // bottom of semicircle
Point(10) = {-23, 0, 0, 1.0};        // nose tip

// =====================================================================
//  CURVES
// =====================================================================

// --- domain ----------------------------------------------------------
Line(1) = {1, 3};        // bottom wall
Line(2) = {3, 4};        // outlet (right edge)
Line(3) = {4, 2};        // top wall
Line(4) = {2, 1};        // inlet (left edge)

// --- body ------------------------------------------------------------
Line(5) = {6, 7};        // top flat side
Line(6) = {7, 8};        // blunt trailing edge (solid base, no slits)
Line(7) = {8, 9};        // bottom flat side
Circle(8) = {9, 5, 10};  // nose, lower quarter
Circle(9) = {10, 5, 6};  // nose, upper quarter

// =====================================================================
//  SURFACE  (one only: domain with the body as a hole)
// =====================================================================
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8, 9};
Plane Surface(1) = {1, 2};

// =====================================================================
//  MESH SIZE & FIELD REFINEMENT
// =====================================================================

// 1. Calculate distance from the body boundary curves
Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8, 9};
Field[1].Sampling = 100;

// 2. Define a threshold to transition size based on distance
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.02;   // Fine element size right at the body surface
Field[2].SizeMax = 1.5;    // Coarse element size far away at domain boundaries
Field[2].DistMin = 0.2;    // Distance up to which SizeMin is maintained
Field[2].DistMax = 12.0;   // Distance at which element size reaches SizeMax

// Set the threshold field as the active background field
Background Field = 2;

// Gmsh global configuration settings
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Smoothing = 20;
Mesh.MshFileVersion = 2.2;

// =====================================================================
//  PHYSICAL GROUPS
// =====================================================================
Physical Curve("inlet") = {4};
Physical Curve("outlet") = {2};
Physical Curve("top_wall") = {3};
Physical Curve("bottom_wall") = {1};
Physical Curve("wall") = {5, 6, 7, 8, 9};
Physical Surface("interior") = {1};