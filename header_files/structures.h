// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <stdbool.h>
#include <stdlib.h>

///////////////////////////////////////////////////////////////////////////////
// Structures
//////////////////////////////////////////////////////////////////////////////

// Structure to represent the common parameters for all grids
struct parameters
    {
    short dimension; //dimension of the problem, default is 3
    short poly_degree; //degree of the polynomial basis functions
    short phs_degree; //degree of the PHS basis functions
    short cloud_size_multiplier; //multiplier for the number of cloud points
    short test;
    }parameters; // child created only once and used globally

// Structure to represent the mesh data
typedef struct PointStructure {   
    char mesh_filename[50]; // name of the mesh file
    int num_nodes; // number of nodes
    int num_corners; // number of corner nodes
    int num_boundary_nodes; // number of boundary nodes
    int num_elem; // number of elements
    double d_avg; // average distance between nodes
    short num_poly_terms; //number of polynomial terms  //////NUM_POLY_TERMS
    short num_cloud_points; //number of cloud points in the domain
    short poly_degree; //degree of the polynomial basis functions
    double* x;  // x coordinates of the nodes
    double* y;  // y coordinates of the nodes
    double* z;  // z coordinates of the nodes
    int* point_index; // index of the point in the original dataset
    double* x_normal; // x component of the normal vector
    double* y_normal; // y component of the normal vector
    double* z_normal; // z component of the normal vector
    bool* boundary_tag; // Boolean tag for boundary nodes
    bool* corner_tag; // Boolean tag for corner nodes
    short* pow_x; //power of x in the polynomial basis functions
    short* pow_y; //power of y in the polynomial basis functions
    short* pow_z; //power of z in the polynomial basis functions
    int** cloud_index; // indices of the INTERPOLATION cloud points 
    int* prolongation_points; // indices of the PROLONGATION cloud points
    int* restriction_points; // indices of the RESTRICTION cloud points
    double** Dx; // x-derivative matrix
    double** Dy; // y-derivative matrix
    double** Dz; // z-derivative matrix
    double** lap; // laplacian matrix
    double** restr_mat;// restriction matrix
    double** prol_mat;// prolongation matrix
}PointStructure;

typedef struct FieldVariables {
    double* u; // x-component of the velocity field
    double* v; // y-component of the velocity field
    double* w; // z-component of the velocity field
    double* p; // pressure field
    double* T; // temperature field
    double rho; // density field
    double mu; // dynamic viscosity field
    double Re; // Reynolds number field
} FieldVariables;

// Structure to represent a point in 3-dimensional space
typedef struct Point {
    double coords[3];
    int index; // Index of the point in the original dataset
} Point;

void AllocateMemoryPointStructure(PointStructure* myPointStruct, int nodes) {
    myPointStruct->x_normal = (double*)malloc(nodes * sizeof(double));
    myPointStruct->y_normal = (double*)malloc(nodes * sizeof(double));
    myPointStruct->z_normal = (double*)malloc(nodes * sizeof(double));
    myPointStruct->x = (double*)malloc(nodes * sizeof(double));
    myPointStruct->y = (double*)malloc(nodes * sizeof(double));
    myPointStruct->z = (double*)malloc(nodes * sizeof(double));
    myPointStruct->point_index = (int*)malloc(nodes * sizeof(int));
    myPointStruct->boundary_tag = (bool*)malloc(nodes * sizeof(bool));
    myPointStruct->corner_tag = (bool*)malloc(nodes * sizeof(bool));
    myPointStruct->cloud_index = (int**)malloc(nodes * sizeof(int*)); // Allocate memory for the cloud index
    for (int i = 0; i < nodes; i++) 
        myPointStruct->cloud_index[i] = (int*)malloc(myPointStruct->num_cloud_points * sizeof(int));
    myPointStruct->num_nodes = nodes;
    myPointStruct->num_elem = 0;
}

void free_PointStructure(PointStructure* myPointStruct) {
    free(myPointStruct->x);
    free(myPointStruct->y);
    free(myPointStruct->z);
    free(myPointStruct->point_index);
    free(myPointStruct->x_normal);
    free(myPointStruct->y_normal);
    free(myPointStruct->z_normal);
    free(myPointStruct->boundary_tag);
    free(myPointStruct->corner_tag);
    free(myPointStruct);
    free(myPointStruct->pow_x);
    free(myPointStruct->pow_y);
    free(myPointStruct->pow_z);
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        free(myPointStruct->cloud_index[i]);
        free(myPointStruct->Dx[i]);
        free(myPointStruct->Dy[i]);
        free(myPointStruct->Dz[i]);
        free(myPointStruct->lap[i]);
        free(myPointStruct->restr_mat[i]);
        free(myPointStruct->prol_mat[i]);
    }
    free(myPointStruct->Dx);
    free(myPointStruct->Dy);
    free(myPointStruct->Dz);
    free(myPointStruct->lap);
    free(myPointStruct->cloud_index);
    free(myPointStruct->prolongation_points);
    free(myPointStruct->restriction_points);
    free(myPointStruct->restr_mat);
    free(myPointStruct->prol_mat);
}

void AllocateMemoryFieldVariables(FieldVariables** field, PointStructure* myPointStruct, int num_levels) {
    *field = (FieldVariables*)malloc(num_levels * sizeof(FieldVariables));
    for (int ii = 0; ii < num_levels; ii++) {
        (*field)[ii].u = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].v = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].w = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].p = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].T = (double*)malloc(myPointStruct[ii].num_nodes * sizeof(double));
    }
}

#endif