// ---------------------------------------------------------------------
// 1. PARAMETERS
// ---------------------------------------------------------------------
DefineConstant[
  D       = {1.0,   Name "Geometry/1Cylinder diameter"},
  Lup     = {5.0,   Name "Geometry/2Upstream length (in D)"},
  Ldown   = {12.0,  Name "Geometry/3Downstream length (in D)"},
  Lside   = {5.0,   Name "Geometry/4Half-height (in D)"},

  hwall   = {5e-3,  Name "Mesh/1Cell size at the cylinder wall"},
  hfar    = {0.08,  Name "Mesh/2Farfield cell size (matches the uniform case)"},
  growth  = {1.10,  Name "Mesh/3Cell-to-cell growth ratio between wall and farfield"}
];

R = D/2.0;
xmin = -Lup*D;   xmax =  Ldown*D;
ymin = -Lside*D; ymax =  Lside*D;

// ---------------------------------------------------------------------
// 2. GEOMETRY (OpenCASCADE Engine)
// ---------------------------------------------------------------------
SetFactory("OpenCASCADE");

// Outer domain points and boundary lines
Point(11) = {xmin, ymin, 0.0};
Point(12) = {xmax, ymin, 0.0};
Point(13) = {xmax, ymax, 0.0};
Point(14) = {xmin, ymax, 0.0};

Line(11) = {11, 12};   // bottom
Line(12) = {12, 13};   // outlet
Line(13) = {13, 14};   // top
Line(14) = {14, 11};   // inlet

// Cylinder curve defined directly by {Center_X, Center_Y, Center_Z, Radius}
Circle(1) = {0.0, 0.0, 0.0, R};

Curve Loop(1) = {11, 12, 13, 14};
Curve Loop(2) = {1};

// Reversing inner loop orientation (-2) for proper hole definition
Plane Surface(1) = {1, -2};

// ---------------------------------------------------------------------
// 3. MESH SIZE FIELD
// ---------------------------------------------------------------------
Field[1] = Distance;
Field[1].CurvesList = {1};
Field[1].Sampling   = 800;

Field[2] = MathEval;
Field[2].F = Sprintf("min(%.10g, %.10g + %.10g*F1)", hfar, hwall, growth - 1.0);

Background Field = 2;

// ---------------------------------------------------------------------
// 4. MESH OPTIONS
// ---------------------------------------------------------------------
Mesh.MeshSizeExtendFromBoundary = 0;   // size driven exclusively by background field
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;
Mesh.MeshSizeFactor             = 1.0;

Mesh.Algorithm    = 6;    // Frontal-Delaunay
Mesh.Optimize     = 1;
Mesh.Smoothing    = 20;
Mesh.ElementOrder = 1;

Mesh.MshFileVersion = 2.2;
Mesh.Binary         = 0;
Mesh.SaveAll        = 0;

// ---------------------------------------------------------------------
// 5. PHYSICAL GROUPS
// ---------------------------------------------------------------------
Physical Curve("cylinder") = {1};
Physical Curve("inlet")    = {14};
Physical Curve("outlet")   = {12};
Physical Curve("bottom")   = {11};
Physical Curve("top")      = {13};
Physical Surface("fluid")  = {1};