// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef FUNCTIONS_SUPPLEMENTARY_H
#define FUNCTIONS_SUPPLEMENTARY_H

#include "structures.h"

/////////////////////////////////////////////////////////////////////////////////
// Function Declarations

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
int write_vtk(char* filename, FieldVariables* field, PointStructure* myPointStruct);

#endif
