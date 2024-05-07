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
    double **Dx = create_full_gradx_matrix(m1, &parameters);
    double **Dy = create_full_grady_matrix(m1, &parameters);
    double **lap = create_full_laplacian_matrix(m1, &parameters);
    double **inv;
    create_matrix(&inv, *m1->n_count, *m1->n_count);

    double *u, *v, *p, *u_star, *v_star, *u_new, *v_new, *p_new, *rhs_u, *rhs_v, *rhs_p, *conv, *adv1, *adv2;
    u = create_vector(*m1->n_count);
    v = create_vector(*m1->n_count);
    p = create_vector(*m1->n_count);
    u_star = create_vector(*m1->n_count);
    v_star = create_vector(*m1->n_count);
    u_new = create_vector(*m1->n_count);
    v_new = create_vector(*m1->n_count);
    p_new = create_vector(*m1->n_count);
    rhs_u = create_vector(*m1->n_count);
    rhs_v = create_vector(*m1->n_count);
    rhs_p = create_vector(*m1->n_count);
    conv = create_vector(*m1->n_count);
    adv1 = create_vector(*m1->n_count);
    adv2 = create_vector(*m1->n_count);

    double steady_state_err = 0;
    double dt = 0.01;
    double Re = 100;
    double tol = 1e-6;
    double p_tol = 1e-6;
    int max_iter = 1000;

    // Boundary conditions enforcing on the derivative matrix
    for (int i = 0; i < *m1->n_count; i++) {
        if (m1->boundary_tag[i] == 1) {
        }
    }
    int iter = 0;
    // Time stepping
    for (iter = 0; iter < max_iter; iter++) {
        // Enforcing boundary conditions
        for (int i = 0; i < *m1->n_count; i++) {
            if (m1->boundary_tag[i] == 1) {
                u[i] = 0; u_star[i] = 0; u_new[i] = 0; rhs_u[i] = 0;
                v[i] = 0; v_star[i] = 0; v_new[i] = 0; rhs_v[i] = 0;
                p[i] = 0; rhs_p[i] = 0;
                if (fabs(m1->points[i].coords[1] - 1) < 1e-5){
                    u[i] = 1; u_star[i] = 1; u_new[i] = 1; rhs_u[i] = 1;
                }
            }
        }
    // u_star = u + dt*((1/Re)*conv(u) - adv(u))
        multiply_matrix_vector(lap, u, conv, *m1->n_count, *m1->n_count);
        multiply_matrix_vector(Dx, u, adv1, *m1->n_count, *m1->n_count);
        multiply_matrix_vector(Dy, u, adv2, *m1->n_count, *m1->n_count);
        for (int i = 0; i < *m1->n_count; i++) {
            u_star[i] = u[i] + dt*((1/Re)*conv[i] - u[i]*adv1[i] - v[i]*adv2[i]);
        }
    // v_star = v + dt*((1/Re)*conv(v) - adv(v))
        multiply_matrix_vector(lap, v, conv, *m1->n_count, *m1->n_count);
        multiply_matrix_vector(Dx, v, adv1, *m1->n_count, *m1->n_count);
        multiply_matrix_vector(Dy, v, adv2, *m1->n_count, *m1->n_count);
        for (int i = 0; i < *m1->n_count; i++) {
            v_star[i] = v[i] + dt*((1/Re)*conv[i] - u[i]*adv1[i] - v[i]*adv2[i]);
        }

    // Pressure Poison
        multiply_matrix_vector(Dx, u_star, adv1, *m1->n_count, *m1->n_count);
        multiply_matrix_vector(Dy, v_star, adv2, *m1->n_count, *m1->n_count);
        for (int i = 0; i < *m1->n_count; i++) {
            rhs_p[i] = (1/dt)*(adv1[i] + adv2[i]);
        }

    // SOR method
    for (int p_iter = 0; p_iter<100; p_iter++) {
        double omega = 1.4;
        double sum = 0;
        for (int i = 0; i < *m1->n_count; i++) {
            sum = 0;
            if (m1->boundary_tag[i] == 1) {
                p[i] = 0; p_new[i] = 0;
                continue;
            }
            for (int j = 1; j < parameters.n_cloud_points; j++) {
                sum += lap[i][m1->cloud_index[i][j]]*p_new[m1->cloud_index[i][j]];
            }
            p_new[i] = (1-omega)*p[i] + (omega/lap[i][i])*(rhs_p[i] - sum);
        }
        double p_err = 0;
        for (int i = 0; i < *m1->n_count; i++) {
            p_err += fabs(p_new[i] - p[i]);
        }
        for (int i = 0; i < *m1->n_count; i++) {
            p[i] = p_new[i];
        }
        if (p_err < p_tol) {
            break;
        }
    }

    // Update the velocity
        multiply_matrix_vector(Dx, p, adv1, *m1->n_count, *m1->n_count);
        multiply_matrix_vector(Dy, p, adv2, *m1->n_count, *m1->n_count);
        for (int i = 0; i < *m1->n_count; i++) {
            u_new[i] = u_star[i] - dt*adv1[i];
            v_new[i] = v_star[i] - dt*adv2[i];
        }
    
    // Calculate the residual
        steady_state_err = 0;
        for (int j = 0; j < *m1->n_count; j++) {
            steady_state_err += fabs(u_new[j]-u[j]);
            steady_state_err += fabs(v_new[j]-v[j]);
        }

        printf("Iteration: %d, Error: %lf\n", iter, steady_state_err);

        // Update the solution
        for (int i = 0; i < *m1->n_count; i++) {
            u[i] = u_new[i];
            v[i] = v_new[i];
        }

        if (fabs(steady_state_err) < tol) {
            break;
        }   

    }

    printf ("No of iterations: %d\n", iter);
    printf("Max iterations: %d\n", max_iter);
    printf("Steady state error: %lf\n", steady_state_err);
    printf("Tolerance: %lf\n", tol);
    printf("No. of grid points: %d\n", *m1->n_count);

    // Write the solution to a file
    file = fopen(strcat(folder1,"/velocity.csv"), "w");
    for (int i = 0; i < *m1->n_count; i++) {
        fprintf(file, "%f, %f\n", u[i], v[i]);
    }
    fclose(file);

    free_matrix(Dx, *m1->n_count);
    free_matrix(Dy, *m1->n_count);
    free_matrix(lap, *m1->n_count);
    free_matrix(inv, *m1->n_count);
    free(u); free(v); free(p); free(u_star); free(v_star); free(u_new); free(v_new); free(rhs_u); free(rhs_v); free(rhs_p); free(conv); free(adv1); free(adv2);
    free_mesh(m1);
    return 0;
}