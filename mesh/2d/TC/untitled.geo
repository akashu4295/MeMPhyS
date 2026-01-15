//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
//+
Circle(2) = {0, 0, 0, 2, 0, 2*Pi};
//+
Physical Curve("inner_wall", 3) = {1};
//+
Physical Curve("outer_wall", 4) = {2};
//+
Curve Loop(1) = {2};
//+
Curve Loop(2) = {1};
//+
Plane Surface(1) = {1, 2};
//+
Mesh.MeshSizeMax = 0.06;