// Author :  Akash Unnikrishnan
// Affiliation : Indian Institute of Technology Gandhinagar
// compile command: gcc main.c -lm
// Run command: ./a.out

#include "header_files/read_mesh_and_parameters.c"
#include "header_files/rbf.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

// extern int errno;

struct stat st;

int main()
{   
    char mesh_filename[50]; // = "mesh/Square_n_60_unstruc.msh";    //Mesh file name
    char folder[50], folder1[50], temp[50];
    FILE *file;
    int e;
    for (int ii = 10; ii<110 ; ii = ii +10){
        sprintf(mesh_filename, "mesh/Square_n_%d_unstruc.msh", ii);
        printf("Mesh file name: %s\n", mesh_filename);
        read_parameters("parameters.csv", mesh_filename, &parameters);     // Read the parameters from the file
        sprintf(folder, "mesh_%d", ii);
        sprintf(folder1, "%s/poly_%d", folder,parameters.poly_degree);

        e = stat(folder1, &st);
        if (errno == ENOENT) {
            e = stat(folder, &st);
            if (errno == ENOENT) 
                mkdir(folder);
            mkdir(folder1);
            }
        
        
        mesh* m1; // Create and read mesh data
        readmesh(&m1, mesh_filename);    

        // Manufactured RHS function
        int k = 1; // wave number
        double* f = (double*)malloc(*m1->n_count*sizeof(double));
        for (int i = 0; i < *m1->n_count; i++) { 
            f[i] = sin(2*3.14*k*m1->points[i].coords[0]) + sin(2*3.14*k*m1->points[i].coords[1]);
        }

        printf("RHS function created\n");

        // create gradient matrices
        double **Dx = create_full_gradx_matrix(m1, &parameters);
        double **Dy = create_full_grady_matrix(m1, &parameters);
        double **Dxx = create_full_gradx_matrix(m1, &parameters);
        double **Dyy = create_full_grady_matrix(m1, &parameters);
        double **lap = create_full_laplacian_matrix(m1, &parameters);
        printf("Full gradient matrices created\n");

        double *fx, *fy, *lapf, *fxx, *fyy;
        fx = create_vector(*m1->n_count);
        fy = create_vector(*m1->n_count);
        lapf = create_vector(*m1->n_count);
        fxx = create_vector(*m1->n_count);
        fyy = create_vector(*m1->n_count);
        multiply_matrix_vector(Dx, f, fx, *m1->n_count, *m1->n_count);
        multiply_matrix_vector(Dy, f, fy, *m1->n_count, *m1->n_count);
        multiply_matrix_vector(Dxx, f, fxx, *m1->n_count, *m1->n_count);
        multiply_matrix_vector(Dyy, f, fyy, *m1->n_count, *m1->n_count);
        multiply_matrix_vector(lap, f, lapf, *m1->n_count, *m1->n_count);

        strcat(folder1,"/");
        strcpy(temp,folder1);
        file = fopen(strcat(temp,"f.csv"), "w");
        for (int i = 0; i < *m1->n_count; i++) {
            fprintf(file, "%f\n", f[i]);
        }
        fclose(file);
        
        strcpy(temp,folder1);
        file = fopen(strcat(temp,"fx.csv"), "w");
        for (int i = 0; i < *m1->n_count; i++) {
            fprintf(file, "%f\n", fx[i]);
        }
        fclose(file);
        
        strcpy(temp,folder1);
        file = fopen(strcat(temp,"fy.csv"), "w");
        for (int i = 0; i < *m1->n_count; i++) {
            fprintf(file, "%f\n", fy[i]);
        }
        fclose(file);
        
        strcpy(temp,folder1);
        file = fopen(strcat(temp,"fxx.csv"), "w");
        for (int i = 0; i < *m1->n_count; i++) {
            fprintf(file, "%f\n", fxx[i]);
        }
        fclose(file);
        
        strcpy(temp,folder1);
        file = fopen(strcat(temp,"fyy.csv"), "w");
        for (int i = 0; i < *m1->n_count; i++) {
            fprintf(file, "%f\n", fyy[i]);
        }
        fclose(file);
        
        strcpy(temp,folder1);
        file = fopen(strcat(temp,"lapf.csv"), "w");
        for (int i = 0; i < *m1->n_count; i++) {
            fprintf(file, "%f\n", lapf[i]);
        }
        fclose(file);

        // write_corner_and_boundary_tags(m1);   // Write the corner and boundary tags
        // write_normals(m1);  // Write the normals at boundary nodes

        free(fx);
        free(fy);
        free(f);
        free(fxx);
        free(fyy);
        free(lapf);
        free_matrix(Dx, *m1->n_count);
        free_matrix(Dy, *m1->n_count);
        free_matrix(Dxx, *m1->n_count);
        free_matrix(Dyy, *m1->n_count);
        free_matrix(lap, *m1->n_count);
        free_mesh(m1);
    }
    return(0);
}