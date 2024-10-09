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
