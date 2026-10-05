// Gmsh project - D-shaped bluff body, NO SLITS, SHARP corners  [L = 1]
// Narrow size range (max/min = 0.105/0.025 = 4.2, i.e. < 5) and gentle growth.
// REFINED 2x vs dbody_noslit_sharp_L1_ratio4_fine.geo: every cell size halved;
// distances, thicknesses and box bounds unchanged.
//
// MODIFIED: nose tip is now 12 D from the inlet (D = body height = 1),
// i.e. body moved +7.5 vs the original.  Second (far-wake) box REMOVED;
// the single wake box follows the body.
//
//   gmsh dbody_noslit_sharp_L1_ratio4_fine2x_12D.geo -2 -o dbody_noslit_sharp_L1_ratio4_fine2x_12D.msh
//
// Geometry after Ananthu, Nair & Narayanan, Phys. Fluids 36, 105157 (2024)
// All coordinates in units of L, with L = 1.  Use nu = 1e-4 for Re_L = 1e4.
SetFactory("OpenCASCADE");

// =====================================================================
//  SIZE PARAMETERS  (edit these to refine / coarsen)
// =====================================================================
s_wall  = 0.025;    // at the body wall             (1x)
s_wake1 = 0.035;    // wake box                     (~91 cells / wavelength 3.2 L)
s_far   = 0.105;    // free stream                  (4.2x)

// =====================================================================
//  POINTS
// =====================================================================
//+
Point(1) = {-27.5, -10, 0, 1.0};
//+
Point(2) = {-27.5, 10, 0, 1.0};
//+
Point(3) = {27.5, -10, 0, 1.0};
//+
Point(4) = {27.5, 10, 0, 1.0};

//+
Point(5) = {-15, 0, 0, 1.0};           // nose arc CENTRE (construction only)
//+
Point(6) = {-15, 0.5, 0, 1.0};       // top of semicircle
//+
Point(7) = {-13, 0.5, 0, 1.0};       // top trailing corner (sharp)
//+
Point(8) = {-13, -0.5, 0, 1.0};      // bottom trailing corner (sharp)
//+
Point(9) = {-15, -0.5, 0, 1.0};      // bottom of semicircle
//+
Point(10) = {-15.5, 0, 0, 1.0};        // nose tip  (12 D from inlet)

// =====================================================================
//  CURVES
// =====================================================================
//+
Line(1) = {1, 3};        // bottom wall
//+
Line(2) = {3, 4};        // outlet (right edge)
//+
Line(3) = {4, 2};        // top wall
//+
Line(4) = {2, 1};        // inlet (left edge)

//+
Line(5) = {6, 7};        // top flat side
//+
Line(6) = {7, 8};        // blunt trailing edge
//+
Line(7) = {8, 9};        // bottom flat side
//+
Circle(8) = {9, 5, 10};  // nose, lower quarter
//+
Circle(9) = {10, 5, 6};  // nose, upper quarter

// =====================================================================
//  SURFACE
// =====================================================================
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5, 6, 7, 8, 9};
//+
Plane Surface(1) = {1, 2};

// =====================================================================
//  MESH SIZE FIELDS
// =====================================================================
// --- body wall band --------------------------------------------------
// 0.025 -> 0.105 over 1.7 L  =>  gradient 0.047 (~5 % growth per cell)
//+
Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8, 9};
Field[1].Sampling = 1200;
//+
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = s_wall;
Field[2].SizeMax = s_far;
Field[2].DistMin = 0.07;
Field[2].DistMax = 1.77;

// --- wake box ----------------------------------------------------------
// 0.035 -> 0.105 over 1.8 L  =>  gradient 0.039
//+
Field[3] = Box;                 // wake: -13.5 <= x/L <= 27.5 (to outlet), |y/L| <= 3
Field[3].VIn = s_wake1;
Field[3].VOut = s_far;
Field[3].XMin = -13.5;  Field[3].XMax = 27.5;
Field[3].YMin = -3.0;   Field[3].YMax = 3.0;
Field[3].ZMin = -1;     Field[3].ZMax = 1;
Field[3].Thickness = 1.8;
//+
Field[5] = Min;
Field[5].FieldsList = {2, 3};
//+
Background Field = 5;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;      // the "1.0" on each Point is ignored
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeMax = s_far;
Mesh.Smoothing = 20;
Mesh.MshFileVersion = 2.2;        // MSH 2.2 ASCII

// =====================================================================
//  PHYSICAL GROUPS
// =====================================================================
//+
Physical Curve("inlet") = {4};
//+
Physical Curve("outlet") = {2};
//+
Physical Curve("top_wall") = {3};
//+
Physical Curve("bottom_wall") = {1};
//+
Physical Curve("wall") = {5, 6, 7, 8, 9};
//+
Physical Surface("interior") = {1};
