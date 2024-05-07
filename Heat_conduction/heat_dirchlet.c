// Author :  Akash Unnikrishnan
// Affiliation : Indian Institute of Technology Gandhinagar
// compile command: gcc main.c -lm
// Run command: ./a.out

#include "../header_files/read_mesh_and_parameters.c"
#include "../header_files/rbf.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <time.h>

#define M_PI 3.14159265358979323846
// extern int errno;

struct stat st;

int main()
{   
    char mesh_filename[50] = "mesh/Square_n_20_unstruc.msh";    //Mesh file name
    char folder[50], folder1[50], temp[50];
    FILE *file;
    int e;

    printf("Mesh file name: %s\n", mesh_filename);
    read_parameters("parameters.csv", mesh_filename, &parameters);     // Read the parameters from the file
    sprintf(folder, "solution");
    sprintf(folder1, "%s/poly_%d", folder,parameters.poly_degree);

    e = stat(folder1, &st);
    // printf("Error number: %d\n", errno);
    if (errno == ENOENT) {
        e = stat(folder, &st);
        if (errno == ENOENT) 
            mkdir(folder);
        mkdir(folder1);
        }
    
    mesh* m1; // Create and read mesh data
    readmesh(&m1, mesh_filename);    

    printf("RHS function created\n");

    // create gradient matrices
    double **Dxx = create_full_dxx_matrix(m1, &parameters);
    double **Dyy = create_full_dyy_matrix(m1, &parameters);
    double **lap = create_full_laplacian_matrix(m1, &parameters);

    double *f, *fxx, *fyy, *lapf;
    fxx = create_vector(*m1->n_count);
    fyy = create_vector(*m1->n_count);
    lapf = create_vector(*m1->n_count);
    f = create_vector(*m1->n_count);

    // Iterative solver for heat conduction equation
    double *T = create_vector(*m1->n_count);
    double *T_new = create_vector(*m1->n_count);
    double *Source = create_vector(*m1->n_count);
    
    double steady_state_err = 0;
    double alpha = 0.01;
    double tol = 1e-6;
    int max_iter = 10000;
    int k = 2; // wave number
    double omega = 1.5;
    // Initialisation
    for (int i = 0; i < *m1->n_count; i++) {
        f[i] = sin(M_PI*k*m1->points[i].coords[0]) + sin(M_PI*k*m1->points[i].coords[1]);
        T[i] =  0;//sin(M_PI*k*m1->points[i].coords[0]) + sin(M_PI*k*m1->points[i].coords[1]);
        Source[i] = 2*M_PI*M_PI*k*k*alpha*(sin(M_PI*k*m1->points[i].coords[0])*sin(M_PI*k*m1->points[i].coords[1]));
        T_new[i] = T[i];
        if (m1->boundary_tag[i] == 1) {
            T[i] = 0;
            T_new[i] = 0;
            Source[i] = 0;
        }
    }

    // Boundary conditions enforcing on the derivative matrix
    for (int i = 0; i < *m1->n_count; i++) {
        if (m1->boundary_tag[i] == 1) {
            for (int j = 0; j < *m1->n_count; j++) {
                Dxx[i][j] = 0;
                Dyy[i][j] = 0;
                lap[i][j] = 0;
            }
            Dxx[i][i] = 1;
            Dyy[i][i] = 1;
            lap[i][i] = 1;
        }
    }
    int iter = 0;
    // Time stepping
    for (iter = 0; iter < max_iter; iter++) {
        // for (int j = 0; j < *m1->n_count; j++) {
        //     if (m1->boundary_tag[j] == 1) {
        //         Source[j] = 0;
        //     }
        // }
        // Calculate the derivatives
        // multiply_matrix_vector(Dx, T, fx, *m1->n_count, *m1->n_count);
        // multiply_matrix_vector(Dy, T, fy, *m1->n_count, *m1->n_count);
        // multiply_matrix_vector(Dxx, T, fxx, *m1->n_count, *m1->n_count);
        // multiply_matrix_vector(Dyy, T, fyy, *m1->n_count, *m1->n_count);
        // multiply_matrix_vector(lap, T, lapf, *m1->n_count, *m1->n_count);

        // Update the temperature
        // for (int j = 0; j < *m1->n_count; j++) {
        //     T_new[j] = T[j] + omega*dt*(alpha*(lapf[j]) + Source[j]);
        // }

        // SOR method
            for (int j = 0; j < *m1->n_count; j++) {
                double sum = 0;
                if (m1->boundary_tag[j] == 1) {
                    T_new[j] = 0;
                    T[j] = 0;
                }
                else {
                    for (int jj = 1; jj<parameters.n_cloud_points; jj++) {
                        sum += lap[j][m1->cloud_index[j][jj]]*T_new[m1->cloud_index[j][jj]];
                    }
                T_new[j] = (1-omega)*T[j] + (omega/lap[j][j])*(Source[j]/alpha - sum);
                }
            }

        // Calculate the residual
        steady_state_err = 0;
        for (int j = 0; j < *m1->n_count; j++) {
            steady_state_err += fabs(T_new[j]-T[j]);
        }

        printf("Iteration: %d, Error: %lf\n", iter, steady_state_err);

        // Update the temperature
        for (int j = 0; j < *m1->n_count; j++) {
            T[j] = T_new[j];
        }
        if (fabs(steady_state_err) < tol) {
            break;
        }   

    }
    double solution_error = 0;
    for (int i = 0; i < *m1->n_count; i++) 
        solution_error += fabs(T[i] - f[i]);

    printf ("No of iterations: %d\n", iter);
    printf("Max iterations: %d\n", max_iter);
    printf("Steady state error: %lf\n", steady_state_err);
    printf("Tolerance: %lf\n", tol);
    printf("No. of grid points: %d\n", *m1->n_count);
    printf("SOR parameter: %lf\n", omega);
    printf("Solution error: %lf\n", sqrt(solution_error));

    // Write the solution to a file
    file = fopen(strcat(folder1,"/T.csv"), "w");
    for (int i = 0; i < *m1->n_count; i++) {
        fprintf(file, "%f\n", T[i]);
    }
    fclose(file);

    free(fxx);
    free(fyy);
    free(T);
    free(T_new);
    free(Source);
    free(lapf);
    free_matrix(Dxx, *m1->n_count);
    free_matrix(Dyy, *m1->n_count);
    free_mesh(m1);
    free(f);
    return 0;
}