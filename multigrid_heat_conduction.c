// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign
// compile command: nvcc read_mesh_data.c
// Run command: ./a.out

///////////////////////////////////////////////////////////////////////////////
////////////// Header files

#include <time.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "header_files/structures.h"
#include "header_files/general_functions.h"
#include "header_files/mesh_functions.h"
#include "header_files/kdtree_functions.h"
#include "header_files/mat_lib.h"
#include "header_files/rbf.h"
#include "header_files/solvers.h"

///////////////////////////////////////////////////////////////////////////////
////////////// Main Program

int main()
{
    clock_t clock_start, clock_program_begin;
    clock_program_begin = clock();    // Start the clock
    short num_levels;
    FILE *file;
    PointStructure* myPointStruct;

///// Read Input File Name and Parameters
    read_flow_parameters("flow_parameters.csv");     // Read the parameters from the file
    read_grid_filenames(&myPointStruct, "grid_filenames.csv", &num_levels);     // Read the parameters from the file
    myPointStruct[0].poly_degree = parameters.poly_degree; // setting polynomial degree for all grids
	for (short ii = 1; ii<num_levels ; ii = ii +1)
        myPointStruct[ii].poly_degree = 3; // setting polynomial degree for all grids

///// Read Mesh data on all levels and compute derivative matrices
    clock_start = clock();    // Start the clock
    read_complete_mesh_data(myPointStruct, num_levels);
    printf("Time taken to read the mesh data: %lf\n", (double)(clock()-clock_start)/CLOCKS_PER_SEC);

    clock_start = clock();    // Start the clock
    for (int ii = 0; ii<num_levels ; ii = ii +1)
        create_derivative_matrices(&myPointStruct[ii]);
    printf("Time taken to create derivative matrices: %lf\n", (double)(clock()-clock_start)/CLOCKS_PER_SEC);

    clock_start = clock();    // Start the clock
    if(parameters.test>0){
        test_derivatives (myPointStruct, num_levels, parameters.dimension);
        printf("Time taken to test the derivatives: %lf\n", (double)(clock()-clock_start)/CLOCKS_PER_SEC);
    }

///// Multilevel solution starts here
    clock_start = clock();    // Start the clock
    double **T, **source, **res, **zeros;
    int k = 1; // wave number
    T = (double**)malloc(num_levels*sizeof(double*));
    source = (double**)malloc(num_levels*sizeof(double*));
    res = (double**)malloc(num_levels*sizeof(double*));
    zeros = (double**)malloc(num_levels*sizeof(double*));
    double* temp = (double*)malloc(myPointStruct[0].num_nodes*sizeof(double)); // For final l2 norm check

    for (short ii = 0; ii<num_levels ; ii = ii +1){
        T[ii] = (double*)malloc(myPointStruct[ii].num_nodes*sizeof(double));
        source[ii] = (double*)malloc(myPointStruct[ii].num_nodes*sizeof(double));
        res[ii] = (double*)malloc(myPointStruct[ii].num_nodes*sizeof(double));
        zeros[ii] = (double*)malloc(myPointStruct[ii].num_nodes*sizeof(double));
    }

    // Set source term for the finest grid
    for (int i = 0; i<num_levels; i++){
        for (int j = 0;j<myPointStruct[i].num_nodes;j++){
            T[i][j] = 0.0; source[i][j] = 0.0; res[i][j] = 0.0; zeros[i][j] = 0.0;
        }
    }
    set_source_term(&myPointStruct[0], source[0], k);

    for (int i = 0; i < myPointStruct[0].num_nodes; i++)
    {
        temp[i] = T[0][i];
        if (!myPointStruct[0].boundary_tag[i])
            T[0][i] = 1;
    }

    clock_start = clock();    // Start the clock

    int num_vcycle = 1;
    int num_iter = 5000;
    double omega = 1.0;
    clock_t start_vcycle, end_vcycle;
    start_vcycle = clock();    // Start the clock
    for (int icycle = 0; icycle < num_vcycle; icycle++){
        clock_start = clock();    // Start the clock
        for (int ilev = 0; ilev < num_levels; ilev++){
            relax_vcycle(&myPointStruct[ilev], T[ilev], source[ilev], num_iter, omega);
            calculate_residual(&myPointStruct[ilev], res[ilev], T[ilev], source[ilev]);
            printf("Residual norm at level %d: %e\n", ilev, l2_norm(res[ilev], zeros[ilev], myPointStruct[ilev].num_nodes));
            if (ilev != num_levels-1){
                restrict_residuals(&myPointStruct[ilev], &myPointStruct[ilev+1], res[ilev], source[ilev+1]);
                }
            }
        printf("\n");

        for (int ilev =num_levels-1; ilev > 0; ilev--){
            prolongate_correction (&myPointStruct[ilev-1], &myPointStruct[ilev], T[ilev-1], T[ilev]);
            //printf("Correction norm at level %d: %e\n", ilev, l2_norm(T[ilev], zeros[ilev], myPointStruct[ilev].num_nodes));
            if (ilev != 0) 
                relax_vcycle(&myPointStruct[ilev], T[ilev], source[ilev], num_iter, omega);
            }
    }
    end_vcycle = clock();  // End the clock
    printf("Time taken for %d vcycles: %lf\n", num_vcycle, (double)(end_vcycle-start_vcycle)/CLOCKS_PER_SEC);
    printf("Time for execution (total): %lf\n", (double)(clock()-clock_program_begin)/CLOCKS_PER_SEC);

// Write the solution to a file
    file = fopen("grid_1.csv", "w");
    for (int ii = 0; ii<myPointStruct[0].num_nodes; ii++)
        fprintf(file, "%lf, %lf, %lf, %lf\n", myPointStruct[0].x[ii], myPointStruct[0].y[ii], myPointStruct[0].z[ii], T[0][ii]);
    file = fopen("grid_2.csv", "w");
    for (int ii = 0; ii<myPointStruct[1].num_nodes; ii++)
        fprintf(file, "%lf, %lf, %lf, %lf\n", myPointStruct[1].x[ii], myPointStruct[1].y[ii], myPointStruct[1].z[ii], T[1][ii]);
    file = fopen("grid_3.csv", "w");
    for (int ii = 0; ii<myPointStruct[2].num_nodes; ii++)
        fprintf(file, "%lf, %lf, %lf, %lf\n", myPointStruct[2].x[ii], myPointStruct[2].y[ii], myPointStruct[2].z[ii], T[2][ii]);  
    fclose (file);


// Free all the memory
    for (int ii = 0; ii<num_levels ; ii = ii +1){
        free(T[ii]);
        free(source[ii]);
        free(res[ii]);
        free(zeros[ii]);
        free_PointStructure(myPointStruct);
    }
    free(T);
    free(source);
    free(myPointStruct);
    free(temp);
    free(res);
    free(zeros);
    
    return 0;
}