# Multigrid_meshless_NS
(IN-WORK) A multigrid (multi-level) accelerated Navier Stokes solver within a meshless framework

At Present: 
1. The meshless (first and second derivative matrices) Dx, Dy, Dz, Dxx, Dyy, Dzz, and Laplacian matrices for a Gmsh ASCII (version 2) .msh file implemented.
2. Poly Harmonic Spline Radial Basis Function with appended polynomials implemented. (  MULTIQUADRICS AND GAUSSIAN IN-WORK )
3. Local Interpolation with point clouds identified through a kd-tree algorithm
4. Heat conduction problem added with multigrid implementation and for Dirchlet boundary conditions
5. SOR solver used for Pressure Poisson and Heat conduction equations

In work
1. Heat conduction equations to be solved for Neumann boundary condition
2. Navier Stokes Solver to be implemented


To run the code use a terminal and compile with gcc compiler:
command: gcc multigrid_heat_conduction.c -lm

Make necessary changes in the grid_filenames.csv and flow_parameters.csv
grid_filenames.csv has the mesh filenames in the order from finest grid to coarse grid
flow_parameters.csv has the details like polynomial degree, phs degree, cloud_size_multiplier and others
