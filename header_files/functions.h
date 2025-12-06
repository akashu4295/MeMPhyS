// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <time.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "structures.h"

// Mesh Function declarations
void read_PointStructure(PointStructure* myPointStruct);
void read_flow_parameters(char *filename);
void calculate_parameters(PointStructure* myPointStruct);
void correct_normal_directions(PointStructure* myPointStruct);
void read_grid_filenames(PointStructure** myPointStruct, char* filename, short* num_levels);
void read_complete_mesh_data(PointStructure* myPointStruct, short num_levels);   
void create_restriction_matrix(PointStructure* myPointStruct, PointStructure* myPointStruct1);
void create_prolongation_matrix(PointStructure* myPointStruct, PointStructure* myPointStruct1);
void rcm_reordering(PointStructure* myPointStruct);
double calculate_dt(PointStructure* myPointStruct);

// RBF Function declarations
double calculate_phs_rbf(double *x, double *c, int phs, int dimension);
void create_A_matrix_from_cloud_indices_vectorised(PointStructure* myPointStruct, double* A, int cloud_index);
void gradx_matrix_vectorised(PointStructure* myPointStruct, double* A, int cloud_index);
void grady_matrix_vectorised(PointStructure* myPointStruct, double* A, int cloud_index);
void gradz_matrix_vectorised(PointStructure* myPointStruct, double* A, int cloud_index);
void laplacian_matrix_vectorised(PointStructure* myPointStruct, double* A, int cloud_index);
void create_full_gradx_matrix_vectorised(PointStructure* myPointStruct);
void create_full_grady_matrix_vectorised(PointStructure* myPointStruct);
void create_full_gradz_matrix_vectorised(PointStructure* myPointStruct);
void create_full_laplacian_matrix_vectorised(PointStructure* myPointStruct);
void create_derivative_matrices_vectorised(PointStructure* myPointStruct);

// Functions to test and calculate errors with a manufactured problem
void calculate_errors_2d(PointStructure* myPointStruct, double* f, double* fx, double* fy, double* lapf, double* F, double* Fx, double* Fy, double* lapF);
void calculate_errors_3d(PointStructure* myPointStruct, double* f, double* fx, double* fy, double* fz, double* lapf, double* F, double* Fx, double* Fy, double* Fz, double* lapF);
void calculate_test_vectors_2d(PointStructure* myPointStruct, double** f, double** fx, double** fy, double** lapf);
void calculate_test_vectors_3d(PointStructure* myPointStruct, double** f, double** fx, double** fy, double** fz, double** lapf);
void set_manufactured_solution_2d(PointStructure* myPointStruct, double** f, double** fx, double** fy, double** lapf, int k);
void set_manufactured_solution_3d(PointStructure* myPointStruct, double** f, double** fx, double** fy, double** fz, double** lapf, int k);
void test_derivatives(PointStructure* myPointStruct, short num_levels, short domain_dimension);
void free_test_vectors(double** f, double** fx, double** fy, double** fz, double** lapf, double*** Dx, double*** Dy, double*** Dz, double*** lap, int num_nodes);
double l2_norm_gen(PointStructure* myPointStruct, double *A, double *B, int n_rows_A);
void test_derivatives(PointStructure* myPointStruct, short num_levels, short domain_dimension);
double calculate_torque_vectorised(PointStructure* myPointStruct, FieldVariables* field);

// Write files Function Declarations
void write_normals(PointStructure* myPointStruct, char* filename);
void write_boundary_tags(PointStructure* myPointStruct, char* filename);
void write_corner_tags(PointStructure* myPointStruct, char* filename);
void write_coordinates(PointStructure* myPointStruct, char* filename);
void write_cloud_index(PointStructure* myPointStruct, char* filename);
void write_prolongation_and_restriction_points(PointStructure* myPointStruct, char* filename);
void write_test_files(double* f, double* fx, double* fy, double* fz, double* lapf, double* fxx, double* fyy, double* fzz, int num_nodes, char* folder1);
void write_processed_grid_data(PointStructure* myPointStruct, int num_levels);
void make_directory(const char* name);
int write_vtk(char *gmsh_filename, FieldVariables *field, PointStructure* myPS);


////////////////////////////////////////////////////////////////////////////////////
// Function Declarations
//////////////////////////////////////////////////////////////////////////////////

// Time Implicit Solver Modules

double time_implicit_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void calculate_intermediate_velocity_implicit_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void calculate_mass_residual_implicit_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void multigrid_Poisson_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void update_velocity_implicit_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void relaxation_vectorised(PointStructure* mypointstruct, FieldVariables* field);
void calculate_residuals_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void restrict_residuals_vectorised(PointStructure* myPointStruct_f,PointStructure* myPointStruct, FieldVariables* field_f, FieldVariables* field);
void prolongate_corrections_vectorised(PointStructure* myPointStruct_f,PointStructure* myPointStruct, FieldVariables* field_f, FieldVariables* field);
void update_boundary_pressure_vectorised(PointStructure* mypointstruct, FieldVariables* field);
void update_boundary_pprime_vectorised(PointStructure* mypointstruct, FieldVariables* field);

double time_implicit_solver_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void calculate_intermediate_velocity_implicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void calculate_mass_residual_implicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void multigrid_Poisson_solver_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void update_velocity_implicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void relaxation_vectorised_2d(PointStructure* mypointstruct, FieldVariables* field);
void calculate_residuals_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void restrict_residuals_vectorised_2d(PointStructure* myPointStruct_f,PointStructure* myPointStruct, FieldVariables* field_f, FieldVariables* field);
void prolongate_corrections_vectorised_2d(PointStructure* myPointStruct_f,PointStructure* myPointStruct, FieldVariables* field_f, FieldVariables* field);
void update_boundary_pressure_vectorised_2d(PointStructure* mypointstruct, FieldVariables* field);
void update_boundary_pprime_vectorised_2d(PointStructure* mypointstruct, FieldVariables* field);



// Fractional Step Explicit Solver Modules

double fractional_step_explicit_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_intermediate_velocity_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_mass_residual_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_boundary_dpdn_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_multigrid_Poisson_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_relaxation_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_residuals_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_restrict_residuals_vectorised(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, FieldVariables* field_f, FieldVariables* field_c);
void FS_prolongate_corrections_vectorised(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, FieldVariables* field_f, FieldVariables* field_c);
void FS_update_velocity_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_update_boundary_pressure_vectorised(PointStructure* myPointStruct, FieldVariables* field);

double fractional_step_explicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_intermediate_velocity_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_mass_residual_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_boundary_dpdn_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_multigrid_Poisson_solver_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_relaxation_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_residuals_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_restrict_residuals_vectorised_2d(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, FieldVariables* field_f, FieldVariables* field_c);
void FS_prolongate_corrections_vectorised_2d(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, FieldVariables* field_f, FieldVariables* field_c);
void FS_update_velocity_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_update_boundary_pressure_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);

#endif