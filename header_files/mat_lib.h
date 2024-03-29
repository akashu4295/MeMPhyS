#include <stdio.h>
#include <stdlib.h>
#include <math.h>  

double** create_matrix(int n_rows, int n_cols);
int** create_matrix_int(int n_rows, int n_cols);
double* create_vector(int n_rows);
void free_matrix(double **A, int n_rows);
void free_matrix_int(int **A, int n_rows);
void free_vector(double *A);
void print_matrix(double **A, int n_rows, int n_cols);
void multiply_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A, int n_cols_B);
void multiply_matrix_vector(double **A, double *B, double *C, int n_rows_A, int n_cols_A);
void multiply_scalar_matrix(double scalar, double **A, double **B, int n_rows_A, int n_cols_A);
void multiply_scalar_vector(double scalar, double *A, double *B, int n_rows_A);
void add_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A);
void add_vectors(double *A, double *B, double *C, int n_rows_A);
void subtract_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A);
void subtract_vectors(double *A, double *B, double *C, int n_rows_A);
void transpose_matrix(double **A, double **B, int n_rows_A, int n_cols_A);
void copy_matrix(double **A, double **B, int n_rows_A, int n_cols_A);
void copy_vector(double *A, double *B, int n_rows_A);
void vector_outer_product(double *A, double *B, double **C, int n_rows_A, int n_rows_B);
double vector_inner_product(double *A, double *B, int n_rows_A);
double vector_norm(double *A, int n_rows_A);
double** identity_matrix(int n_rows);
// double** InverseOfMatrix(double** matrix, int order);
double** inverseOfMatrix(double** matrix, int order);


double** create_matrix(int n_rows, int n_cols)
{
    int i;
    double **A;
    A = (double **)malloc(n_rows * sizeof(double *));
    for (i = 0; i < n_rows; i++)
    {
        A[i] = (double *)malloc(n_cols * sizeof(double));
    }
    return A;
}

int** create_matrix_int(int n_rows, int n_cols)
{
    int i;
    int **A;
    A = (int **)malloc(n_rows * sizeof(int *));
    for (i = 0; i < n_rows; i++)
    {
        A[i] = (int *)malloc(n_cols * sizeof(int));
    }
    return A;
}

double* create_vector(int n_rows)
{
    double *A;
    A = (double *)malloc(n_rows * sizeof(double));
    return A;
}

void free_matrix(double **A, int n_rows)
{
    int i;
    for (i = 0; i < n_rows; i++)
    {
        free(A[i]);
    }
    free(A);
}

void free_matrix_int(int **A, int n_rows)
{
    int i;
    for (i = 0; i < n_rows; i++)
    {
        free(A[i]);
    }
    free(A);
}

void free_vector(double *A)
{
    free(A);
}

void print_matrix(double **A, int n_rows, int n_cols)
{
    int i, j;
    for (i = 0; i < n_rows; i++)
    {
        for (j = 0; j < n_cols; j++)
        {
            printf("%g ", A[i][j]);
        }
        printf("\n");
    }
}

void multiply_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A, int n_cols_B)
{
    int i, j, k;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_B; j++)
        {
            C[i][j] = 0;
            for (k = 0; k < n_cols_A; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void multiply_matrix_vector(double **A, double *B, double *C, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        C[i] = 0;
        for (j = 0; j < n_cols_A; j++)
        {
            C[i] += A[i][j] * B[j];
        }
    }
}

void multiply_scalar_matrix(double scalar, double **A, double **B, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_A; j++)
        {
            B[i][j] = scalar * A[i][j];
        }
    }
}

void multiply_scalar_vector(double scalar, double *A, double *B, int n_rows_A)
{
    int i;
    for (i = 0; i < n_rows_A; i++)
    {
        B[i] = scalar * A[i];
    }
}

void add_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_A; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

void add_vectors(double *A, double *B, double *C, int n_rows_A)
{
    int i;
    for (i = 0; i < n_rows_A; i++)
    {
        C[i] = A[i] + B[i];
    }
}

void subtract_matrices(double **A, double **B, double **C, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_A; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

void subtract_vectors(double *A, double *B, double *C, int n_rows_A)
{
    int i;
    for (i = 0; i < n_rows_A; i++)
    {
        C[i] = A[i] - B[i];
    }
}

void transpose_matrix(double **A, double **B, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_A; j++)
        {
            B[j][i] = A[i][j];
        }
    }
}

void copy_matrix(double **A, double **B, int n_rows_A, int n_cols_A)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_cols_A; j++)
        {
            B[i][j] = A[i][j];
        }
    }
}

void copy_vector(double *A, double *B, int n_rows_A)
{
    int i;
    for (i = 0; i < n_rows_A; i++)
    {
        B[i] = A[i];
    }
}

void vector_outer_product(double *A, double *B, double **C, int n_rows_A, int n_rows_B)
{
    int i, j;
    for (i = 0; i < n_rows_A; i++)
    {
        for (j = 0; j < n_rows_B; j++)
        {
            C[i][j] = A[i] * B[j];
        }
    }
}

double vector_inner_product(double *A, double *B, int n_rows_A)
{
    int i;
    double result = 0;
    for (i = 0; i < n_rows_A; i++)
    {
        result += A[i] * B[i];
    }
    return result;
}

double vector_norm(double *A, int n_rows_A)
{
    int i;
    double result = 0;
    for (i = 0; i < n_rows_A; i++)
    {
        result += A[i] * A[i];
    }
    return sqrt(result);
}

double** identity_matrix(int n_rows)
{
    int i, j;
    double **A;
    A = (double **)malloc(n_rows * sizeof(double *));
    for (i = 0; i < n_rows; i++)
    {
        A[i] = (double *)malloc(n_rows * sizeof(double));
    }
    for (i = 0; i < n_rows; i++)
    {
        for (j = 0; j < n_rows; j++)
        {
            if (i == j)
            {
                A[i][j] = 1;
            }
            else
            {
                A[i][j] = 0;
            }
        }
    }
    return A;
}

// // Function to perform the inverse operation on the matrix.
// double** InverseOfMatrix(double** matrix, int order)
// {
//     double* temp; // for swapping rows
//     double** augment = create_matrix(order, 2*order);
//     // Create the augmented matrix
//     for (int i = 0; i < order; i++) {
//         for (int j = 0; j < 2 * order; j++) {
//             if (j == (i + order))
//                 augment[i][j] = 1;
//             else if (j < order)
//                 augment[i][j] = matrix[i][j];
//             else
//                 augment[i][j] = 0;
//         }
//     }
 
//     // Interchange the row of matrix,
//     // interchanging of row will start from the last row
//     for (int i = order - 1; i > 0; i--) {
//         if (augment[i - 1][0] < augment[i][0]) {
//             double* temp = augment[i];
//             augment[i] = augment[i - 1];
//             augment[i - 1] = temp;
//         }
//     }
//     free(temp);

//     double temp1;
//     // Replace a row by sum of itself and a
//     // constant multiple of another row of the matrix
//     for (int i = 0; i < order; i++) {
//         for (int j = 0; j < order; j++) {
//             if (j != i) {
//                 temp1 = augment[j][i] / augment[i][i];
//                 for (int k = 0; k < 2 * order; k++) {
//                     augment[j][k] -= augment[i][k] * temp1;
//                 }
//             }
//         }
//     }
 
//     // Multiply each row by a nonzero integer.
//     // Divide row element by the diagonal element
//     for (int i = 0; i < order; i++) {
//         temp1 = augment[i][i];
//         for (int j = 0; j < 2 * order; j++) {
//             augment[i][j] = augment[i][j] / temp1;
//         }
//     }
//     double** inverse = create_matrix(order, order);
//     for (int i = 0; i < order; i++) {
//         for (int j = 0; j < order; j++) {
//             inverse[i][j] = augment[i][j + order];
//         }
//     }
//     free_matrix(augment, order);
//     return inverse;
// }

double** inverseOfMatrix(double** matrix1, int order)
{
    double temp; 
    double** matrix = create_matrix(order, 2*order);
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            matrix[i][j] = matrix1[i][j];
        }
    }
    // Create the augmented matrix
    for (int i = 0; i < order; i++) {
        for (int j = order; j < 2 * order; j++) {
            // Add '1' at the diagonal places of
            // the matrix to create a identity matrix
            if (j == (i + order))
                matrix[i][j] = 1;
            else
                matrix[i][j] = 0;
        }
    }
 
    
    print_matrix(matrix, order,2* order);
    // Interchange the row of matrix,
    // interchanging of row will start from the last row
    for (int i = order - 1; i > 0; i--) {
 
        // Swapping each and every element of the two rows
        // if (matrix[i - 1][0] < matrix[i][0])
        // for (int j = 0; j < 2 * order; j++) {
        //
        //        // Swapping of the row, if above
        //        // condition satisfied.
        // temp = matrix[i][j];
        // matrix[i][j] = matrix[i - 1][j];
        // matrix[i - 1][j] = temp;
        //    }
 
        // Directly swapping the rows using pointers saves
        // time
 
        if (matrix[i - 1][0] < matrix[i][0]) {
            double* temp = matrix[i];
            matrix[i] = matrix[i - 1];
            matrix[i - 1] = temp;
        }
    }
 
    // Print matrix after interchange operations.
    printf("\n=== Augmented Matrix ===\n");
    print_matrix(matrix, order, order * 2);
 
    // Replace a row by sum of itself and a
    // constant multiple of another row of the matrix
    for (int i = 0; i < order; i++) {
 
        for (int j = 0; j < order; j++) {
 
            if (j != i) {
 
                temp = matrix[j][i] / matrix[i][i];
                for (int k = 0; k < 2 * order; k++) {
 
                    matrix[j][k] -= matrix[i][k] * temp;
                }
            }
        }
    }
 
    // Multiply each row by a nonzero integer.
    // Divide row element by the diagonal element
    for (int i = 0; i < order; i++) {
 
        temp = matrix[i][i];
        for (int j = 0; j < 2 * order; j++) {
 
            matrix[i][j] = matrix[i][j] / temp;
        }
    }
    double** inverse = create_matrix(order, order);
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            inverse[i][j] = matrix[i][j + order];
        }
    }
    free_matrix(matrix, order);
    return inverse;
}