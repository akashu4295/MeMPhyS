# Multigrid_meshless_NS
(IN-WORK) A multigrid (multi-level) accelerated Navier Stokes solver within a meshless framework

At Present: 
1. The meshless (first and second derivative matrices) Dx, Dy, Dz, Dxx, Dyy, Dzz, and Laplacian matrices for a Gmsh ASCII (version 2) .msh file implemented.
2. Poly Harmonic Spline Radial Basis Function with appended polynomials implemented. (  MULTIQUADRICS AND GAUSSIAN IN-WORK )
3. Local Interpolation with point clouds identified through a kd-tree algorithm
4. Navier Stokes equations are solved numerically using a fractional step method in 2D, lid-driven cavity problem added
5. Heat conduction problem added
6. SOR solver used for Pressure Poisson and Heat conduction equations

BUGS:
1. TESTED UP TO POLYNOMIAL DEGREE 4; HIGHER POLYNOMIAL DEGREES CREATE ERRONEOUS DATA
2. SOLUTION CONVERGES WELL, HOWEVER ACCURACY IN COMPARISON WITH ANALYTICAL DATA IS NOT GREAT
