// Author :  Akash Unnikrishnan
// Affiliation : Indian Institute of Technology Gandhinagar
// compile command: gcc main.c -lm
// Run command: ./a.out

#include "header_files/read_mesh_and_parameters.c"
#include "header_files/rbf.h"

int main()
{   
    char mesh_filename[50] = "mesh/Square_n_100_unstruc.msh";    //Mesh file name

    read_parameters("parameters.csv", mesh_filename, &parameters);     // Read the parameters from the file
    
    mesh* m1; // Create and read mesh data
    readmesh(&m1, mesh_filename);    

    // Testing the code
    int test_index = 4;
    printf("Test index: %d\n", test_index);
    printf("Cloud points: ");
    for (int i = 0; i < parameters.n_cloud_points; i++) {
        printf("%d ", m1->cloud_index[test_index][i]);
    }
    printf("\n");
    // Manufactured RHS function
    int k = 1; // wave number
    double* f = (double*)malloc(*m1->n_count*sizeof(double));
    for (int i = 0; i < *m1->n_count; i++) { 
        f[i] = sin(2*3.14*k*m1->points[i].coords[0]) + sin(2*3.14*k*m1->points[i].coords[1]);
    }

    // create A matrix and rhs vector
    double **A = create_A_matrix_from_cloud_indices(m1, m1->cloud_index[test_index], parameters.n_cloud_points, &parameters);
    double **B = gradx_matrix(m1, m1->cloud_index[test_index], parameters.n_cloud_points, &parameters);
    double *rhs = create_rhs_vector_from_cloud_indices(m1, f, m1->cloud_index[test_index], parameters.n_cloud_points, &parameters);
    write_matrix_to_file(A, parameters.n_cloud_points+parameters.n_poly_terms, parameters.n_cloud_points+parameters.n_poly_terms, "A_matrix.csv");
    write_matrix_to_file(B, parameters.n_cloud_points+parameters.n_poly_terms, parameters.n_cloud_points+parameters.n_poly_terms, "B_matrix.csv");
    
    // Inverse of A matrix
    double **A_inv;
    create_matrix(&A_inv, parameters.n_cloud_points+parameters.n_poly_terms, parameters.n_cloud_points+parameters.n_poly_terms);
    matrixInverse_Gauss_Jordan(A, A_inv, parameters.n_cloud_points+parameters.n_poly_terms);
    //matrixInverse_LU(A, A_inv, parameters.n_cloud_points+parameters.n_poly_terms);
    // write_matrix_to_file(A_inv, parameters.n_cloud_points+parameters.n_poly_terms, parameters.n_cloud_points+parameters.n_poly_terms, "A_inv_matrix.csv");

    double **B1;
    create_matrix(&B1, parameters.n_cloud_points+parameters.n_poly_terms, parameters.n_cloud_points+parameters.n_poly_terms);
    multiply_matrices(B, A_inv, B1, parameters.n_cloud_points+parameters.n_poly_terms, parameters.n_cloud_points+parameters.n_poly_terms, parameters.n_cloud_points+parameters.n_poly_terms); 
    // write_matrix_to_file(B1, parameters.n_cloud_points+parameters.n_poly_terms, parameters.n_cloud_points+parameters.n_poly_terms, "B1_matrix.csv");


    // Lambdas and gammas in sol vector
    double* lambda= create_vector(parameters.n_cloud_points+parameters.n_poly_terms);
    multiply_matrix_vector(A_inv, rhs, lambda, parameters.n_cloud_points+parameters.n_poly_terms, parameters.n_cloud_points+parameters.n_poly_terms);

    double pt[3] = {0,0,0}; // central point for testing
    for (int i = 0; i < parameters.n_cloud_points; i++) {
        for (int j = 0; j < 3; j++) {
            pt[j] += m1->points[m1->cloud_index[test_index][i]].coords[j];
        }
    }
    for (int j = 0; j < 3; j++) {
        pt[j] /= parameters.n_cloud_points;
    }

    // Numerical solution
    double* sol = create_vector(parameters.n_cloud_points+parameters.n_poly_terms);
    for (int i = 0; i < parameters.n_cloud_points; i++) {
        sol[i] = calculate_phs_rbf(&m1->points[m1->cloud_index[test_index][i]], &pt[0], parameters.phs_degree, parameters.dimension);
    }
    for (int i = 0; i < parameters.n_poly_terms; i++) {
        sol[i+parameters.n_cloud_points] = pow(pt[0], parameters.pow_x[i])
                        *pow(pt[1], parameters.pow_y[i])
                        *pow(pt[2], parameters.pow_z[i]);
    }
    double sol_num = vector_inner_product(sol, lambda, parameters.n_cloud_points+parameters.n_poly_terms);
    
    printf("Numerical Solution: %lf\n", sol_num);
    printf("Exact Solution: %lf\n", sin(2*3.14*k*pt[0]) + sin(2*3.14*k*pt[1]));
    free_matrix(A, parameters.n_cloud_points+parameters.n_poly_terms);
    free_matrix(A_inv, parameters.n_cloud_points+parameters.n_poly_terms);

    // write_normals(m1, "normals.txt");  // Write the normals at boundary nodes
    // write_corner_tags(m1, "corner_tags.txt");   // Write the corner and boundary tags
    // write_boundary_tags(m1, "boundary_tags.txt"); // Write the boundary tags
    // write_cloud_index(m1, "cloud_index.txt");  // Write the cloud index

    printf("Average edge length = %lf\n", m1->d_avg); // Average edge length of the mesh
    free_mesh(m1); // free the mesh data 
    free(f);
    free(rhs);
    free(lambda);
    free(sol);
    free_matrix(B, parameters.n_cloud_points);
    free_matrix(B1, parameters.n_cloud_points);

    return(0);
}