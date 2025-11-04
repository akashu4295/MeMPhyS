// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef FUNCTIONS_H
#define FUNCTIONS_H
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

#endif