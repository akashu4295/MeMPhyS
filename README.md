# Multigrid_meshless_NS

IN-WORK
1. A multigrid (multi-level) accelerated Navier Stokes solver to achieve steady state using TIMPLE a meshless framework
2. A compressible flow solver
   
Current Status: 
1. The meshless (first and second derivative matrices) Dx, Dy, Dz, and Laplacian matrices for a Gmsh ASCII (version 2) .msh file implemented.
2. Poly Harmonic Spline Radial Basis Function with appended polynomials implemented. (  MULTIQUADRICS AND GAUSSIAN IN-WORK )
3. Local Interpolation with point clouds identified through a kd-tree algorithm
4. Heat conduction problem added with multigrid implementation and for Dirchlet boundary conditions
5. SOR solver used for Pressure Poisson and Heat conduction equations
6. Heat conduction equations solved for Neumann boundary condition
7. Navier Stokes Solver implemented with Fractional Step
8. Navier Stokes Solver implemented with Time Implicit solver
9. Multigrided Poisson solver

To run the code use a terminal and compile with gcc compiler:
command: gcc <code_file_name.c> -lm

Make necessary changes in the grid_filenames.csv and flow_parameters.csv
grid_filenames.csv has the mesh filenames in the order from finest grid to coarse grid
flow_parameters.csv has the details like polynomial degree, phs degree, cloud_size_multiplier and others
