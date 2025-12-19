# Meshless Method Navier-Stokes Solver (C Code)

This repository contains a C-based implementation of a **Meshless Method** for solving the Navier-Stokes equations. The solver utilizes a radial basis function (RBF) based meshless approach with optional multicore CPU acceleration using OpenACC.

---

### Quick Notes on status of the code

### IN-WORK/ TODO
1. Implementing various boundary conditions
2. Extending OpenACC pragmas to implement GPU acceleration.
3. A multigrid (multi-level) accelerated Navier Stokes solver to achieve steady state using TIMPLE in meshless framework
4. A compressible flow solver
   
### Current Status: 
1. The meshless (first and second derivative matrices) Dx, Dy, Dz, and Laplacian matrices for a Gmsh ASCII (version 2) .msh file implemented.
2. Poly Harmonic Spline Radial Basis Function with appended polynomials implemented. (  MULTIQUADRICS AND GAUSSIAN IN-WORK )
3. Local Interpolation with point clouds identified through a kd-tree algorithm
4. Heat conduction problem added with multigrid implementation and for Dirchlet boundary conditions
5. SOR solver used for Pressure Poisson and Heat conduction equations
6. Heat conduction equations solved for Neumann boundary condition
7. Navier Stokes Solver implemented with Fractional Step
8. Navier Stokes Solver implemented with Time Implicit solver
9. Multigrided Poisson solver

**Make necessary changes in the grid_filenames.csv and flow_parameters.csv
grid_filenames.csv has the mesh filenames in the order from finest grid to coarse grid
flow_parameters.csv has the details like polynomial degree, phs degree, cloud_size_multiplier and others**

---

### Compilation

To compile the code for different hardware architectures:

**Single-core CPU:**
```bash
gcc @sources.txt -lm
```

**Multi-core CPU:**

```bash
nvc -acc -ta=multicore @sources.txt
```

**GPU (NVIDIA):**

```bash
nvc -acc -gpu=managed @sources.txt
```

---

### Running the Code

After compiling:

```bash
./a.out
```

---

## Installation of GUI wrapper with python
Since we are in the developing stage, we are using a conda environment. If you don't have conda installed in the system, please do the following:

### Linux
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
```bash
bash Miniconda3-latest-Linux-x86_64.sh
```
```bash
source ~/.bashrc
```
```bash
rm ./Miniconda3-latest-Linux-x86_64.sh
```

### Windows (Powershell commands)
```bash
curl -o Miniconda3-latest-Windows-x86_64.exe `
     https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe
```
```bash
Start-Process -Wait .\Miniconda3-latest-Windows-x86_64.exe `
  -ArgumentList "/InstallationType=JustMe",
                "/AddToPath=0",
                "/RegisterPython=0",
                "/S",
                "/D=$env:USERPROFILE\Miniconda3"
```
```bash
$env:USERPROFILE\Miniconda3\Scripts\conda.exe init powershell
```
```bash
Remove-Item Miniconda3-latest-Windows-x86_64.exe 
```
```PS
conda install -c conda-forge m2w64-gcc
```

After installing conda environment run this (for both Linux and Windows system):
```bash
conda env create -f environment.yml
conda activate memphys_gui
```bash
To start the GUI:
```bash
python gui.py
```
---

### Folder Structure

* `header_files/`: Contains all the custom C header files used in the solver.
* `init/`: Holds sample initialization scripts.
* `flow_parameters.csv`: Contains simulation parameters. Can be edited by the user.
* `grid_filenames.csv`: Lists available mesh files for use in the simulation.

---

### `flow_parameters.csv`

This file defines the core physical and numerical parameters for the simulation. Each line is of the format `parameter_name,value`. Here are some key parameters:

* `domain_dimensions`: Dimensionality of the domain (2 or 3).
* `poly_deg`: Degree of the appended polynomial (recommended: \[2–15]).
* `phs_deg`: Degree of the polyharmonic spline (odd integers: 3, 5, 7, 9...).
* `cloud_size_multiplier`: Controls the local stencil size; typically between 1.5 and 2.5 for 2d cases and between 3 to 4 for 3d cases.
* `courant_number`: CFL condition; controls the time step size.
* `steady_tolerance`: Tolerance for steady-state convergence.
* `poisson_solver_tolerance`: Tolerance for pressure Poisson solver.
* `num_vcycles`: Number of V-cycles in multigrid.
* `num_relax`: Number of relaxation steps.
* `num_time_steps`: Total simulation time steps.
* `write_interval`: How often to write output.
* `Re`: Reynolds number.
* `restart_filename`: If `restart` is 1, this file will be loaded.

You can customize these settings to suit your problem configuration.

---

### `grid_filenames.csv`

* First line: `num_levels` – Number of mesh refinement levels or hierarchies.
* Following lines: Paths to `.msh` mesh files generated using Gmsh from fine grid to coarse grid.

#### Using Custom Mesh Files

To use your own mesh:

1. Create a mesh in Gmsh.
2. Export it in **ASCII format (version 2)**.
3. Add its path to `grid_filenames.csv`.
4. Ensure the mesh is consistent with expected domain dimensions.

---

### Notes

* Modify `flow_parameters.csv` as needed before compilation.

---

### Sample Cases

Several 2d and 3d mesh cases are included under the paths:

* `mesh/2d/TC/`
* `mesh/2d/SQ/`
* `mesh/3d/CUBE/`
* `mesh/3d/SP_in_SP/`

---

### Author

Dr. Akash Unnikrishnan developed this code as part of his PhD work with the guidance of Prof. Surya Pratap Vanka from University of Illinois at Urbana Champaign and Prof. Vinod Narayanan at Indian Institute of Technology Gandhinagar. Special mention has to go for Dr. Shantanu Shahane, who developed the first version of memphys which is available here: github.com/shahaneshantanu/memphys
