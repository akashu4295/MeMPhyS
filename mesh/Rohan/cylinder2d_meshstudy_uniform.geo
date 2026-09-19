// 1. PARAMETERS
// ---------------------------------------------------------------------
DefineConstant[
  D      = {1.0,  Name "Geometry/1Cylinder diameter"},
  Lup    = {5.0,  Name "Geometry/2Upstream length (in D)"},
  Ldown  = {12.0, Name "Geometry/3Downstream length (in D)"},
  Lside  = {5.0,  Name "Geometry/4Half-height (in D)"},

  h      = {0.08, Name "Mesh/1Uniform cell size everywhere (farfield == wall)"}
];

R = D/2.0;
xmin = -Lup*D;   xmax =  Ldown*D;
ymin = -Lside*D; ymax =  Lside*D;

// ---------------------------------------------------------------------
// 2. GEOMETRY
// ---------------------------------------------------------------------
Point(1) = { 0.0, 0.0, 0.0};       // centre
Point(2) = {  R,  0.0, 0.0};
Point(3) = { 0.0,  R,  0.0};
Point(4) = { -R,  0.0, 0.0};
Point(5) = { 0.0, -R,  0.0};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Point(11) = {xmin, ymin, 0.0};
Point(12) = {xmax, ymin, 0.0};
Point(13) = {xmax, ymax, 0.0};
Point(14) = {xmin, ymax, 0.0};

Line(11) = {11, 12};   // bottom
Line(12) = {12, 13};   // outlet
Line(13) = {13, 14};   // top
Line(14) = {14, 11};   // inlet

Curve Loop(1) = {11, 12, 13, 14};
Curve Loop(2) = {1, 2, 3, 4};
Plane Surface(1) = {1, 2};

// ---------------------------------------------------------------------
// 3. SIZE FIELD  --  one constant value everywhere, no distance/growth
// ---------------------------------------------------------------------
Field[1] = MathEval;
Field[1].F = Sprintf("%.10g", h);

Background Field = 1;

// ---------------------------------------------------------------------
// 4. MESH OPTIONS
// ---------------------------------------------------------------------
Mesh.MeshSizeExtendFromBoundary = 0;   // nothing but the field sets sizes
Mesh.MeshSizeFromPoints         = 0;
Mesh.MeshSizeFromCurvature      = 0;
Mesh.MeshSizeFactor             = 1.0;

Mesh.Algorithm    = 6;    // Frontal-Delaunay: most isotropic triangles
Mesh.Optimize     = 1;
Mesh.Smoothing    = 20;   // Laplacian smoothing -> evens out residual size variation
Mesh.ElementOrder = 1;

Mesh.MshFileVersion = 2.2;
Mesh.Binary         = 0;
Mesh.SaveAll        = 0;

// ---------------------------------------------------------------------
// 5. PHYSICAL GROUPS
// ---------------------------------------------------------------------
Physical Curve("cylinder") = {1, 2, 3, 4};
Physical Curve("inlet")    = {14};
Physical Curve("outlet")   = {12};
Physical Curve("bottom")   = {11};
Physical Curve("top")      = {13};
Physical Surface("fluid")  = {1};
