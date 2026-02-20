// Author :  Akash Unnikrishnan
// Affiliation : Indian Institute of Technology Gandhinagar

#ifndef MATH_LIBRARY_H
#define MATH_LIBRARY_H
#include "structures.h"

////////////////////////////////////////////////////////////////////////
// Function Declarations
////////////////////////////////////////////////////////////////////////

void create_matrix(double ***A, int n_rows, int n_cols);
void create_matrix_int(int ***A, int n_rows, int n_cols);
double** create_matrix1(int n_rows, int n_cols);
double* create_matrix_vectorised(int n_rows, int n_cols);
double* create_vector(int n_rows);
void free_matrix(double **A, int n_rows);
void free_matrix_int(int **A, int n_rows);
void free_vector(double *A);
void print_matrix(double **A, int n_rows, int n_cols);
void multiply_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A, int n_cols_B);
void multiply_matrices_vectorised(double *A, double *B, double *C, int n_rows_A, int n_cols_A, int n_cols_B);
void multiply_matrix_vector(double **A, double *B, double *C, int n_rows_A, int n_cols_A);
void multiply_vector_matrix(double *B, double **A, double **C, int n_rows_A, int n_cols_A);
void multiply_scalar_matrix(double scalar, double **A, double **B, int n_rows_A, int n_cols_A);
void multiply_scalar_vector(double scalar, double *A, double *B, int n_rows_A);
void multiply_vector_matrix_columnwise(double *B, double **A, double *C, int n_rows_A, int n_cols_A);
void multiply_vector_matrix_columnwise_vectorised(double *B, double *A, double *C, int n_rows_A, int n_cols_A);
void add_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A);
void add_matrices_to_first(double **A, double **B, int n_rows_A, int n_cols_A);
void add_vectors(double *A, double *B, double *C, int n_rows_A);
void subtract_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A);
void subtract_vectors(double *A, double *B, double *C, int n_rows_A);
void transpose_matrix(double **A, double **B, int n_rows_A, int n_cols_A);
void copy_matrix(double **A, double **B, int n_rows_A, int n_cols_A);
void copy_vector(double *A, double *B, int n_rows_A);
void swap_vectors(double *A, double *B, int n_rows_A);
void vector_outer_product(double *A, double *B, double **C, int n_rows_A, int n_rows_B);
double vector_inner_product(double *A, double *B, int n_rows_A);
double vector_norm(double *A, int n_rows_A);
double** identity_matrix(int n_rows);
void matrixInverse_Gauss_Jordan(double** matrix1, double** matrix, int order);
void matrixInverse_Gauss_Jordan_vectorised(double* matrix1, double* matrix, int order);
void luDecomposition(double **A, int n, int *P);
void forwardSubstitution(double **L, double *y, double *b, int n, int *P);
void backwardSubstitution(double **U, double *x, double *y, int n);
void matrixInverse_LU(double **A, double **A_inv, int n);
void printMatrix(double **A, int n);
void matrixInverse_Gauss_Jordan2(double** matrix1, double** inverse, int order);
void write_matrix_to_file(double **A, int n_rows, int n_cols, char *filename);
double l2_norm(double *A, double *B, int n_rows_A);
void linear_system_solver(double** A, double* x, double* b, int n);
double min3(double a, double b, double c);

void multiply_sparse_matrix_vector(double** D_coeff, double* f, double* dfdx, int** cloud, int n_rows_D, int n_cols_D);
void multiply_sparse_vector_matrix(double* f, double** D_coeff, double** ftimesD, int n_rows_D, int n_cols_D);
void multiply_sparse_matrix_vector_gpu(double** D_coeff, double* f, double* dfdx, int** cloud, int n_rows_D, int n_cols_D);
void multiply_sparse_matrix_vector_vectorised(double* D_coeff, double* f, double* dfdx, int* cloud, int n_rows_D, int n_cols_D);
void multiply_sparse_matrix_vector_vectorised_gpu(double* D_coeff, double* f, double* dfdx, int* cloud, int n_rows_D, int n_cols_D);
void multiply_sparse_matrix_vector_vectorised_gpu_async(double* D_coeff, double* f, double* dfdx, int* cloud, int n_rows_D, int n_cols_D, int async_queue);


void SpMV_Laplace_2D(PointStructure* ps, const double* x, double* y);
int BiCGStab_Solve(PointStructure* ps, const double* b, double* x, int max_iter, double tol);

#endif