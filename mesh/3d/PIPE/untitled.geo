SetFactory("OpenCASCADE");

// Geometry parameters
R  = 1.0;
L  = 10.0;
lc = 0.3;

// Create cylinder (volume)
Cylinder(1) = {0, 0, 0,   0, 0, L,   R};

// Physical groups (optional but recommended)
Physical Volume("Fluid") = {1};

// Identify boundary surfaces automatically
Physical Surface("Wall")   = {1};
Physical Surface("Outlet") = {2};
Physical Surface("Inlet")  = {3};

// --- Mesh controls ---
Mesh.Algorithm  = 6;   // 2D: Delaunay
Mesh.Algorithm3D = 1;  // 3D: TetGen (robust, no Netgen)
Mesh.Optimize = 1;

Mesh.MeshSizeMax = 0.06;