// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef RBF_H
#define RBF_H

#include "structures.h"
#include "general_functions.h"
#include "mesh_functions.h"
#include "kdtree_functions.h"
#include "mat_lib.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Function declarations

double calculate_phs_rbf(double *x, double *c, int phs, int dimension);
void create_A_matrix_from_cloud_indices(PointStructure* myPointStruct, double** A, int* cloud);
void gradx_matrix(PointStructure* myPointStruct, double** gradx, int* cloud);
void grady_matrix(PointStructure* myPointStruct, double** grady, int* cloud);
void gradz_matrix(PointStructure* myPointStruct, double** gradz, int* cloud);
void laplacian_matrix(PointStructure* myPointStruct, double** lap, int* cloud);
void create_full_gradx_matrix(PointStructure* myPointStruct);
void create_full_grady_matrix(PointStructure* myPointStruct);
void create_full_gradz_matrix(PointStructure* myPointStruct);
void create_full_laplacian_matrix(PointStructure* myPointStruct);
void create_derivative_matrices(PointStructure* myPointStruct);

// Function definitions
double calculate_phs_rbf(double *pt1, double *pt2, int phs, int dimension) {
    double sum = 0;
    for (int i = 0; i < dimension; i++) {
        sum += pow(pt1[i] - pt2[i], 2);
    }
    if (sum == 0)
        return 0;
    return pow(sum, phs*0.5);
}

void create_A_matrix_from_cloud_indices(PointStructure* myPointStruct, double** A, int* cloud) {
    int m = myPointStruct->num_cloud_points;
    int n = myPointStruct->num_poly_terms;
    double pt1[3], pt2[3];
    double seed_pt[3] = {myPointStruct->x[cloud[0]], myPointStruct->y[cloud[0]], myPointStruct->z[cloud[0]]};

    // Parallelize the outer loop and the inner loop for the phi matrix with PHS RBF function
    for (int i = 0; i < m; i++) {
        A[i][i] = 0;
        pt1[0] = myPointStruct->x[cloud[i]];
        pt1[1] = myPointStruct->y[cloud[i]];
        pt1[2] = myPointStruct->z[cloud[i]];
        for (int j = i + 1; j < m; j++) {
            pt2[0] = myPointStruct->x[cloud[j]];
            pt2[1] = myPointStruct->y[cloud[j]];
            pt2[2] = myPointStruct->z[cloud[j]];
            A[i][j] = calculate_phs_rbf(pt1, pt2, parameters.phs_degree, parameters.dimension);
            A[j][i] = A[i][j];
        }
    }

    // Parallelize the loop for adding polynomial terms
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j + m] = 1;
            if (myPointStruct->pow_x[j] != 0)
                A[i][j + m] *= pow(myPointStruct->x[cloud[i]] - seed_pt[0], myPointStruct->pow_x[j]);
            if (myPointStruct->pow_y[j] != 0)
                A[i][j + m] *= pow(myPointStruct->y[cloud[i]] - seed_pt[1], myPointStruct->pow_y[j]);
            if (myPointStruct->pow_z[j] != 0)
                A[i][j + m] *= pow(myPointStruct->z[cloud[i]] - seed_pt[2], myPointStruct->pow_z[j]);

            A[j + m][i] = A[i][j + m];
        }
    }
}

void gradx_matrix(PointStructure* myPointStruct, double** grad, int* cloud) {
    int m = myPointStruct->num_cloud_points;
    int n = myPointStruct->num_poly_terms;
    double pt1[3], pt2[3];
    double seed_pt[3] = {myPointStruct->x[cloud[0]], myPointStruct->y[cloud[0]], myPointStruct->z[cloud[0]]};

    // Parallelizing the outer and inner loops for PHS RBF function
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]];
        pt1[1] = myPointStruct->y[cloud[i]];
        pt1[2] = myPointStruct->z[cloud[i]];
        for (int j = 0; j < m; j++) {
            pt2[0] = myPointStruct->x[cloud[j]];
            pt2[1] = myPointStruct->y[cloud[j]];
            pt2[2] = myPointStruct->z[cloud[j]];
            grad[i][j] = parameters.phs_degree * (pt1[0] - pt2[0]) * calculate_phs_rbf(pt1, pt2, parameters.phs_degree - 2, parameters.dimension);
        }
    }

    // Parallelizing the polynomial terms loop
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]];
        pt1[1] = myPointStruct->y[cloud[i]];
        pt1[2] = myPointStruct->z[cloud[i]];
        grad[i][m] = 0;
        for (int j = 1; j < n; j++) {
            if (myPointStruct->pow_x[j] == 0) {
                grad[i][j + m] = 0;
            } else {
                grad[i][j + m] = myPointStruct->pow_x[j] * pow(pt1[0] - seed_pt[0], myPointStruct->pow_x[j] - 1);
                grad[i][j + m] *= pow(pt1[1] - seed_pt[1], myPointStruct->pow_y[j]);
                grad[i][j + m] *= pow(pt1[2] - seed_pt[2], myPointStruct->pow_z[j]);
            }
            grad[j + m][i] = grad[i][j + m];
        }
    }
}

void grady_matrix(PointStructure* myPointStruct, double** grad, int* cloud) {
    int m = myPointStruct->num_cloud_points;
    int n = myPointStruct->num_poly_terms;
    double pt1[3], pt2[3];
    double seed_pt[3] = {myPointStruct->x[cloud[0]],myPointStruct->y[cloud[0]],myPointStruct->z[cloud[0]]};

    // Parallelizing the outer and inner loops for PHS RBF function
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]]; 
        pt1[1] = myPointStruct->y[cloud[i]]; 
        pt1[2] = myPointStruct->z[cloud[i]];
        for (int j = 0; j < m; j++) { 
            pt2[0] = myPointStruct->x[cloud[j]];
            pt2[1] = myPointStruct->y[cloud[j]];
            pt2[2] = myPointStruct->z[cloud[j]];
            grad[i][j] = parameters.phs_degree*(pt1[1]-pt2[1]) * calculate_phs_rbf(pt1, pt2, parameters.phs_degree-2, parameters.dimension);
        }
    }

    // add polynomial terms
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]]; 
        pt1[1] = myPointStruct->y[cloud[i]]; 
        pt1[2] = myPointStruct->z[cloud[i]]; 
        grad[i][0] = 0;
        for (int j = 1; j < n; j++) {
            if (myPointStruct->pow_y[j] == 0)
                grad[i][j+m] = 0;
            else{
                grad[i][j+m] = myPointStruct->pow_y[j]*pow(pt1[1]-seed_pt[1], myPointStruct->pow_y[j]-1);
                grad[i][j+m] = grad[i][j+m] * pow(pt1[0]-seed_pt[0], myPointStruct->pow_x[j]);
                grad[i][j+m] = grad[i][j+m] * pow(pt1[2]-seed_pt[2], myPointStruct->pow_z[j]);
            }
            grad[j+m][i] = grad[i][j+m];
            // printf("grad[%d][%d] = %lf\n", i, j+m, grad[i][j+m]);
        }
    }
}

void gradz_matrix(PointStructure* myPointStruct, double** grad, int* cloud) {
    int m = myPointStruct->num_cloud_points;
    int n = myPointStruct->num_poly_terms;
    double pt1[3], pt2[3];
    double seed_pt[3] = {myPointStruct->x[cloud[0]],myPointStruct->y[cloud[0]],myPointStruct->z[cloud[0]]};

    // Parallelizing the outer and inner loops for PHS RBF function
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]]; 
        pt1[1] = myPointStruct->y[cloud[i]]; 
        pt1[2] = myPointStruct->z[cloud[i]];
        for (int j = 0; j < m; j++){ 
            pt2[0] = myPointStruct->x[cloud[j]];
            pt2[1] = myPointStruct->y[cloud[j]];
            pt2[2] = myPointStruct->z[cloud[j]];
            grad[i][j] = parameters.phs_degree*(pt1[2]-pt2[2]) * calculate_phs_rbf(pt1, pt2, parameters.phs_degree-2, parameters.dimension);
        }
    }

    // add polynomial terms
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]]; 
        pt1[1] = myPointStruct->y[cloud[i]]; 
        pt1[2] = myPointStruct->z[cloud[i]];
        grad[i][0] = 0;
        for (int j = 1; j < n; j++) {
            if (myPointStruct->pow_z[j] == 0)
                grad[i][j+m] = 0;
            else{
                grad[i][j+m] = myPointStruct->pow_z[j]*pow(pt1[2]-seed_pt[2], myPointStruct->pow_z[j]-1);
                grad[i][j+m] = grad[i][j+m] * pow(pt1[1]-seed_pt[1], myPointStruct->pow_y[j]);
                grad[i][j+m] = grad[i][j+m] * pow(pt1[0]-seed_pt[0], myPointStruct->pow_x[j]);
            }
            grad[j+m][i] = grad[i][j+m];
        }
    }
}

void laplacian_matrix(PointStructure* myPointStruct, double** lap, int* cloud) {
    int m = myPointStruct->num_cloud_points;
    int n = myPointStruct->num_poly_terms;
    double pt1[3], pt2[3];
    double seed_pt[3] = {myPointStruct->x[cloud[0]], myPointStruct->y[cloud[0]], myPointStruct->z[cloud[0]]};

    // Parallelizing the outer loop for calculating the Laplacian matrix
    for (int i = 0; i < m; i++) {
        pt1[0] = myPointStruct->x[cloud[i]]; 
        pt1[1] = myPointStruct->y[cloud[i]]; 
        pt1[2] = myPointStruct->z[cloud[i]];
        for (int j = 0; j < m; j++) {
            pt2[0] = myPointStruct->x[cloud[j]];
            pt2[1] = myPointStruct->y[cloud[j]];
            pt2[2] = myPointStruct->z[cloud[j]];
            lap[i][j] = parameters.phs_degree * parameters.phs_degree *
                         calculate_phs_rbf(pt1, pt2, parameters.phs_degree - 2, parameters.dimension);
        }
    }

    // Adding polynomial terms
    double dtempx, dtempy, dtempz;
    
    if (parameters.dimension == 3) {
        for (int i = 0; i < m; i++) {
            pt1[0] = myPointStruct->x[cloud[i]] - seed_pt[0]; 
            pt1[1] = myPointStruct->y[cloud[i]] - seed_pt[1]; 
            pt1[2] = myPointStruct->z[cloud[i]] - seed_pt[2];
            for (int j = 0; j < n; j++) {      
                dtempx = 0; dtempy = 0; dtempz = 0;
                if (myPointStruct->pow_x[j] > 1) {  
                    dtempx =  myPointStruct->pow_x[j] * (myPointStruct->pow_x[j] - 1) *
                              pow(pt1[0], myPointStruct->pow_x[j] - 2);
                    dtempx *= pow(pt1[1], myPointStruct->pow_y[j]);
                    dtempx *= pow(pt1[2], myPointStruct->pow_z[j]);
                }
                else if (myPointStruct->pow_y[j] > 1) {  
                    dtempy =  myPointStruct->pow_y[j] * (myPointStruct->pow_y[j] - 1) *
                              pow(pt1[1], myPointStruct->pow_y[j] - 2);
                    dtempy *= pow(pt1[0], myPointStruct->pow_x[j]);
                    dtempy *= pow(pt1[2], myPointStruct->pow_z[j]);
                }
                else if (myPointStruct->pow_z[j] > 1) {  
                    dtempz =  myPointStruct->pow_z[j] * (myPointStruct->pow_z[j] - 1) *
                              pow(pt1[2], myPointStruct->pow_z[j] - 2);
                    dtempz *= pow(pt1[0], myPointStruct->pow_x[j]);
                    dtempz *= pow(pt1[1], myPointStruct->pow_y[j]);
                }
                lap[i][j + m] = dtempx + dtempy + dtempz;
                lap[j + m][i] = lap[i][j + m];
            }
        }
    }

    if (parameters.dimension == 2) {
        for (int i = 0; i < m; i++) {
            pt1[0] = myPointStruct->x[cloud[i]] - seed_pt[0]; 
            pt1[1] = myPointStruct->y[cloud[i]] - seed_pt[1]; 
            pt1[2] = myPointStruct->z[cloud[i]] - seed_pt[2];
            for (int j = 0; j < n; j++) {        
                dtempx = 0; dtempy = 0;
                if (myPointStruct->pow_x[j] > 1) {  
                    dtempx =  myPointStruct->pow_x[j] * (myPointStruct->pow_x[j] - 1) *
                              pow(pt1[0], myPointStruct->pow_x[j] - 2);
                    dtempx *= pow(pt1[1], myPointStruct->pow_y[j]);
                }
                else if (myPointStruct->pow_y[j] > 1) {  
                    dtempy =  myPointStruct->pow_y[j] * (myPointStruct->pow_y[j] - 1) *
                              pow(pt1[1], myPointStruct->pow_y[j] - 2);
                    dtempy *= pow(pt1[0], myPointStruct->pow_x[j]);
                }
                lap[i][j + m] = dtempx + dtempy;
                lap[j + m][i] = lap[i][j + m];
            }
        }
    }
}

// Following function creates the full derivative matrices for the mesh

void create_full_gradx_matrix(PointStructure* myPointStruct) {
    double **grad, **A, **A_inv, **B1;
    int m = myPointStruct->num_cloud_points + myPointStruct->num_poly_terms;
    
    // Allocate memory for matrices
    create_matrix(&grad, m, m);
    create_matrix(&A, m, m);
    create_matrix(&A_inv, m, m);
    create_matrix(&B1, m, m);

    // Parallelize the outer loop with OpenACC
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        if (myPointStruct->corner_tag[i] == 0) {
            // Sequential parts remain within the loop
            create_A_matrix_from_cloud_indices(myPointStruct, A, myPointStruct->cloud_index[i]);
            gradx_matrix(myPointStruct, grad, myPointStruct->cloud_index[i]);
            matrixInverse_Gauss_Jordan(A, A_inv, m);
            multiply_matrices(grad, A_inv, B1, m, m, m);
            
            // Parallelize the inner loop (assignment) if it's large enough
            for (int j = 0; j < myPointStruct->num_cloud_points; j++) {
                myPointStruct->Dx[i][j] = B1[0][j];
            }
        }
    }

    // Free matrices
    free_matrix(A, m);
    free_matrix(grad, m);
    free_matrix(A_inv, m);
    free_matrix(B1, m);
}

void create_full_grady_matrix(PointStructure* myPointStruct) {
    double **grad, **A, **A_inv, **B1;
    create_matrix(&grad, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    create_matrix(&A, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    create_matrix(&A_inv, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    create_matrix(&B1, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);

    // Parallelize the outer loop with OpenACC
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        if (myPointStruct->corner_tag[i] == 0){
            create_A_matrix_from_cloud_indices(myPointStruct, A, myPointStruct->cloud_index[i]);
            grady_matrix(myPointStruct, grad, myPointStruct->cloud_index[i]);
            matrixInverse_Gauss_Jordan(A, A_inv, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
            multiply_matrices(grad, A_inv, B1, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
            for (int j = 0; j < myPointStruct->num_cloud_points; j++) {
                myPointStruct->Dy[i][j] = B1[0][j];
            }
        }
    }
    free_matrix(A, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    free_matrix(grad, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    free_matrix(A_inv, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    free_matrix(B1, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
}

void create_full_gradz_matrix(PointStructure* myPointStruct) {
    double **grad, **A, **A_inv, **B1;
    create_matrix(&grad, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    create_matrix(&A, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    create_matrix(&A_inv, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    create_matrix(&B1, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);

    // Parallelize the outer loop with OpenACC
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        if (myPointStruct->corner_tag[i] == 0){
            create_A_matrix_from_cloud_indices(myPointStruct, A, myPointStruct->cloud_index[i]);
            gradz_matrix(myPointStruct, grad, myPointStruct->cloud_index[i]);
            matrixInverse_Gauss_Jordan(A, A_inv, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
            multiply_matrices(grad, A_inv, B1, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
            for (int j = 0; j < myPointStruct->num_cloud_points; j++) {
                myPointStruct->Dz[i][j] = B1[0][j];
            }
        }
    }
    free_matrix(A, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    free_matrix(grad, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    free_matrix(A_inv, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    free_matrix(B1, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
}

void create_full_laplacian_matrix(PointStructure* myPointStruct) {
    double **lap, **A, **A_inv, **B1;
    create_matrix(&lap, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    create_matrix(&A, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    create_matrix(&A_inv, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    create_matrix(&B1, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    
    // Parallelize the outer loop with OpenACC
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        if (myPointStruct->corner_tag[i] == 0){
            create_A_matrix_from_cloud_indices(myPointStruct, A, myPointStruct->cloud_index[i]);
            laplacian_matrix(myPointStruct, lap, myPointStruct->cloud_index[i]);
            matrixInverse_Gauss_Jordan(A, A_inv, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
            multiply_matrices(lap, A_inv, B1, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
            for (int j = 0; j < myPointStruct->num_cloud_points; j++) {
                myPointStruct->lap[i][j] = B1[0][j];
            }
        }
    }
    free_matrix(A, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    free_matrix(lap, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    free_matrix(A_inv, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
    free_matrix(B1, myPointStruct->num_cloud_points+myPointStruct->num_poly_terms);
}

void create_laplacian_for_Poisson(PointStructure* myPointStruct) {
    // Parallelize the outer loop with OpenACC
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        if (myPointStruct->corner_tag[i] == 0){
            for (int j = 0; j < myPointStruct->num_cloud_points; j++) {
                myPointStruct->lap_Poison[i][j] = myPointStruct->lap[i][j];
            }
        }
        if (myPointStruct->boundary_tag[i] == 1){
            for (int j = 0; j < myPointStruct->num_cloud_points; j++) {
                myPointStruct->lap_Poison[i][j] = myPointStruct->Dx[i][j] * myPointStruct->x_normal[i] + myPointStruct->Dy[i][j] * myPointStruct->y_normal[i];
                if (parameters.dimension == 3){
                    myPointStruct->lap_Poison[i][j] += myPointStruct->Dz[i][j] * myPointStruct->z_normal[i];
                }
            }
        }
    }
    // for (int i = 0; i < myPointStruct->num_nodes; i++) {
    //     for (int j = 0; j < myPointStruct->num_cloud_points; j++) {
    //         myPointStruct->lap_Poison[i][j] = myPointStruct->lap_Poison[i][j]/myPointStruct->lap_Poison[i][0];
    //     }
    // }
}

void create_derivative_matrices(PointStructure* myPointStruct){
    create_matrix(&(myPointStruct->Dx), myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    create_matrix(&(myPointStruct->Dy), myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    create_matrix(&(myPointStruct->lap), myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    if (parameters.dimension == 3){
        create_matrix(&(myPointStruct->Dz), myPointStruct->num_nodes, myPointStruct->num_cloud_points);
        create_full_gradx_matrix(myPointStruct);
        create_full_grady_matrix(myPointStruct);
        create_full_gradz_matrix(myPointStruct);
        create_full_laplacian_matrix(myPointStruct);
    }
    if (parameters.dimension == 2){
        create_full_gradx_matrix(myPointStruct);
        create_full_grady_matrix(myPointStruct);
        create_full_laplacian_matrix(myPointStruct);
    }
 //   if (!parameters.compressible){
        create_matrix(&(myPointStruct->lap_Poison), myPointStruct->num_nodes, myPointStruct->num_cloud_points);
        create_laplacian_for_Poisson(myPointStruct);
 //   }
}
#endif
