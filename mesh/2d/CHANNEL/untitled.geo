//+
SetFactory("OpenCASCADE");
// Geometry
Point(1) = {0, 0, 0, 1};
Point(2) = {10, 0, 0, 1};
Point(3) = {10, 1, 0, 1};
Point(4) = {0, 1, 0, 1};

Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // top
Line(4) = {4, 1}; // left

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Physical groups
Physical Line("bottom") = {1};
Physical Line("outlet") = {2};
Physical Line("top")    = {3};
Physical Line("inlet")  = {4};

Physical Surface("fluid") = {1};

Mesh.MeshSizeMax = 0.5;