// Author :  Akash Unnikrishnan
// Affiliation : Indian Institute of Technology Gandhinagar

double calculate_phs_rbf(Point *x, double *c, int phs, int dim) {
    double sum = 0;
    for (int i = 0; i < dim; i++) {
        sum += pow(x->coords[i] - c[i], 2);
    }
    return pow(sum, phs*0.5);
}

double* create_rhs_vector_from_cloud_indices(mesh* m1, double* f, int* cloud, int n, struct parameters *parameters)
{
    double* rhs = (double*)malloc((parameters->n_cloud_points+parameters->n_poly_terms) * sizeof(double));
    for (int j = 0; j < parameters->n_cloud_points; j++)
        rhs[j] = f[cloud[j]];
    for (int j = 0; j < parameters->n_poly_terms; j++)
        rhs[j+parameters->n_cloud_points] = 0;
    return rhs;
}

double** create_A_matrix_from_cloud_indices(mesh* m1, int* cloud, int n, struct parameters *parameters){
    int m = parameters->n_cloud_points;
    if (n<m)
        m = n;
    double **A;
    create_matrix(&A, m+parameters->n_poly_terms, m+parameters->n_poly_terms);
    // phi matrix with phs rbf function
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            A[i][j] = calculate_phs_rbf(&m1->points[cloud[i]], &m1->points[cloud[j]].coords[0], parameters->phs_degree, parameters->dimension);
            // if (i == 0)
            //     printf("%d\n", cloud[j]);
        }
    }
    // add polynomial terms
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < parameters->n_poly_terms; j++) {
            A[i][j+m] = 1;
            // printf("%d\n",i);
            // printf("%lf, %lf, %lf\n", fabs(m1->points[cloud[i]].coords[0] - 0.0), fabs(m1->points[cloud[i]].coords[1] - 0.0), fabs(m1->points[cloud[i]].coords[2] - 0.0));
            if (parameters->pow_x[j] == 0){
                if (fabs(m1->points[cloud[i]].coords[0] - 0.0)>1e-5)
                    A[i][j+m] *= pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]);
                else
                    A[i][j+m] *= 1;
            }
            else
                A[i][j+m] *= pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]);
            
            if (parameters->pow_y[j] == 0){    
                if (fabs(m1->points[cloud[i]].coords[1] - 0.0)>1e-5)   
                    A[i][j+m] *= pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]);
                else
                    A[i][j+m] *= 1;
            }
            else
                A[i][j+m] *= pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]);
            
            if (parameters->pow_z[j] == 0){
                if (fabs(m1->points[cloud[i]].coords[2] - 0.0)>1e-5)  
                    A[i][j+m] *= pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]);
                else
                    A[i][j+m] *= 1;
            }
            else
                A[i][j+m] *= pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]);

            A[j+m][i] = A[i][j+m];
        }
    }
    return A;
}

double** gradx_matrix(mesh* m1, int* cloud, int n, struct parameters *parameters) {
    int m = parameters->n_cloud_points;
    if (n<m)
        m = n;
    double **grad;
    create_matrix(&grad, m+parameters->n_poly_terms, m+parameters->n_poly_terms);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++)
            grad[i][j] = parameters->phs_degree*(m1->points[cloud[i]].coords[0] - m1->points[cloud[j]].coords[0]) * calculate_phs_rbf(&m1->points[cloud[i]], &m1->points[cloud[j]].coords[0], parameters->phs_degree-2, parameters->dimension);
    }
    // add polynomial terms
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < parameters->n_poly_terms; j++) {
            if (parameters->pow_x[j]-1<0){
                grad[i][j+m] = 0; 
                continue;
            }
            else if (parameters->pow_x[j]-1 == 0)
                grad[i][j+m] = 1;
            else
                grad[i][j+m] = parameters->pow_x[j]*pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]-1);
            if (parameters->pow_y[j] != 0)
                grad[i][j+m] = grad[i][j+m] * pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]);
            if (parameters->pow_z[j] != 0)
                grad[i][j+m] = grad[i][j+m] * pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]);
            grad[j+m][i] = grad[i][j+m];
        }
    }
    return grad;
}

double** grady_matrix(mesh* m1, int* cloud, int n, struct parameters *parameters) {
    int m = parameters->n_cloud_points;
    if (n<m)
        m = n;
    double **grad;
    create_matrix(&grad, m+parameters->n_poly_terms, m+parameters->n_poly_terms);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) 
            grad[i][j] = parameters->phs_degree*(m1->points[cloud[i]].coords[1] - m1->points[cloud[j]].coords[1]) * calculate_phs_rbf(&m1->points[cloud[i]], &m1->points[cloud[j]].coords[0], parameters->phs_degree-2, parameters->dimension);
    }
    // add polynomial terms
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < parameters->n_poly_terms; j++) {
            if (parameters->pow_y[j]-1<0){
                grad[i][j+m] = 0; 
                continue;
            }
            else if (parameters->pow_y[j]-1 == 0)
                grad[i][j+m] = 1;
            else
                grad[i][j+m] = parameters->pow_y[j]*pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]-1);
            if (parameters->pow_x[j] != 0)
                grad[i][j+m] = grad[i][j+m] * pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]);
            if (parameters->pow_z[j] != 0)
                grad[i][j+m] = grad[i][j+m] * pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]);
            grad[j+m][i] = grad[i][j+m];
        }
    }
    return grad;
}

double** gradz_matrix(mesh* m1, int* cloud, int n, struct parameters *parameters) {
    int m = parameters->n_cloud_points;
    if (n<m)
        m = n;
    double **grad;
    create_matrix(&grad, m+parameters->n_poly_terms, m+parameters->n_poly_terms);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) 
            grad[i][j] = parameters->phs_degree*(m1->points[cloud[i]].coords[2] - m1->points[cloud[j]].coords[2]) * calculate_phs_rbf(&m1->points[cloud[i]], &m1->points[cloud[j]].coords[0], parameters->phs_degree-2, parameters->dimension);
    }
    // add polynomial terms
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < parameters->n_poly_terms; j++) {
            if (parameters->pow_z[j]-1<0){
                grad[i][j+m] = 0; 
                continue;
            }
            else if (parameters->pow_z[j]-1 == 0)
                grad[i][j+m] = 1;
            else
                grad[i][j+m] = parameters->pow_z[j]*pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]-1);
            if (parameters->pow_y[j] != 0)
                grad[i][j+m] = grad[i][j+m] * pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]);
            if (parameters->pow_x[j] != 0)
                grad[i][j+m] = grad[i][j+m] * pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]);
            grad[j+m][i] = grad[i][j+m];
        }
    }
    return grad;
}

// double** laplacian_matrix(mesh* m1, int* cloud, int n, struct parameters *parameters) {
//     int m = parameters->n_cloud_points;
//     if (n<m)
//         m = n;
//     double **lap;
//     create_matrix(&lap, m+parameters->n_poly_terms, m+parameters->n_poly_terms);
//     for (int i = 0; i < m; i++) {
//         for (int j = 0; j < m; j++) 
//             lap[i][j] = parameters->phs_degree*parameters->phs_degree*calculate_phs_rbf(&m1->points[cloud[i]], &m1->points[cloud[j]].coords[0], parameters->phs_degree-2, parameters->dimension);
//     }
//     // add polynomial terms
//     for (int i = 0; i < m; i++) {
//         for (int j = 0; j < parameters->n_poly_terms; j++) {
//             if (parameters->dimension == 3)
//                 lap[i][j+m] = parameters->pow_x[j]*(parameters->pow_x[j]-1)*pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]-2)*pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j])*pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j])
//                         + parameters->pow_y[j]*(parameters->pow_y[j]-1)*pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j])*pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]-2)*pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j])
//                         + parameters->pow_z[j]*(parameters->pow_z[j]-1)*pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j])*pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j])*pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]-2);
//             else
//                 lap[i][j+m] = parameters->pow_x[j]*(parameters->pow_x[j]-1)*pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]-2)*pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j])
//                         + parameters->pow_y[j]*(parameters->pow_y[j]-1)*pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j])*pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]-2);
//             if (isnan(lap[i][j+m]))
//                 lap[i][j+m] = 0;
//             lap[j+m][i] = lap[i][j+m];
//         }
//     }
//     return lap;
// }

double** laplacian_matrix(mesh* m1, int* cloud, int n, struct parameters *parameters) {
    int m = parameters->n_cloud_points;
    if (n<m)
        m = n;
    double **lap;
    double dtemp;
    create_matrix(&lap, m+parameters->n_poly_terms, m+parameters->n_poly_terms);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) 
            lap[i][j] = parameters->phs_degree*parameters->phs_degree*calculate_phs_rbf(&m1->points[cloud[i]], &m1->points[cloud[j]].coords[0], parameters->phs_degree-2, parameters->dimension);
    }
    // add polynomial terms
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < parameters->n_poly_terms; j++) {
            if (parameters->dimension == 3){
                if (parameters->pow_x[j]-2<0 || parameters->pow_x[j]-1<0)
                    lap[i][j+m] = 0; 
                else{
                    lap[i][j+m] = parameters->pow_x[j]*(parameters->pow_x[j]-1)*pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]-2);
                    if (parameters->pow_y[j] != 0)
                        lap[i][j+m] *= pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]);
                    if (parameters->pow_z[j] != 0)
                        lap[i][j+m] *= pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]);
                }
                if (parameters->pow_y[j]-2<0 || parameters->pow_y[j]-1<0)
                    lap[i][j+m] += 0; 
                else {
                    dtemp = parameters->pow_y[j]*(parameters->pow_y[j]-1)*pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]-2);
                    if (parameters->pow_x[j] != 0)
                        dtemp *= pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]);
                    if (parameters->pow_z[j] != 0)
                        dtemp *= pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]);
                    lap[i][j+m] += dtemp;
                }
                if (parameters->pow_z[j]-2<0 || parameters->pow_z[j]-1<0)
                    lap[i][j+m] += 0; 
                else {
                    dtemp = parameters->pow_z[j]*(parameters->pow_z[j]-1)*pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]-2);
                    if (parameters->pow_x[j] != 0)
                        dtemp *= pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]);
                    if (parameters->pow_y[j] != 0)
                        dtemp *= pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]);
                    lap[i][j+m] += dtemp;
                }
            }
            else{
                if (parameters->pow_x[j]-2<0 || parameters->pow_x[j]-1<0)
                    lap[i][j+m] = 0; 
                else {
                    lap[i][j+m] = parameters->pow_x[j]*(parameters->pow_x[j]-1)*pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]-2);
                    if (parameters->pow_y[j] != 0)
                        lap[i][j+m] *= pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]);
                }
                if (parameters->pow_y[j]-2<0 || parameters->pow_y[j]-1<0)
                    lap[i][j+m] += 0; 
                else {
                    dtemp = parameters->pow_y[j]*(parameters->pow_y[j]-1)*pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]-2);
                    if (parameters->pow_x[j] != 0)
                        dtemp *= pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]);
                    lap[i][j+m] += dtemp;
                }
            }
            // if (isnan(lap[i][j+m]))
            //     lap[i][j+m] = 0;
            lap[j+m][i] = lap[i][j+m];
        }
    }
    return lap;
}

double** dxx_matrix(mesh* m1, int* cloud, int n, struct parameters *parameters) {
    int m = parameters->n_cloud_points;
    if (n<m)
        m = n;
    double **grad;
    create_matrix(&grad, m+parameters->n_poly_terms, m+parameters->n_poly_terms);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++)
            grad[i][j] = parameters->phs_degree*(parameters->phs_degree-2)*pow((m1->points[cloud[i]].coords[0] - m1->points[cloud[j]].coords[0]),2) * calculate_phs_rbf(&m1->points[cloud[i]], &m1->points[cloud[j]].coords[0], parameters->phs_degree-4, parameters->dimension) + parameters->phs_degree * calculate_phs_rbf(&m1->points[cloud[i]], &m1->points[cloud[j]].coords[0], parameters->phs_degree-2, parameters->dimension);
    }
    // add polynomial terms
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < parameters->n_poly_terms; j++) {
            if (fabs(m1->points[cloud[i]].coords[0] - 0)>1e-5)
                grad[i][j+m] = parameters->pow_x[j]*(parameters->pow_x[j]-1)*pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]-2);
            if (fabs(m1->points[cloud[i]].coords[1] - 0)>1e-5)
                grad[i][j+m] = grad[i][j+m] * pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]);
            if (fabs(m1->points[cloud[i]].coords[2] - 0)>1e-5)
                {printf("%f\n", pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]));
                grad[i][j+m] = grad[i][j+m] * pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]);}
            grad[j+m][i] = grad[i][j+m];
        }
    }
    return grad;
}

double** dyy_matrix(mesh* m1, int* cloud, int n, struct parameters *parameters) {
    int m = parameters->n_cloud_points;
    if (n<m)
        m = n;
    double **grad;
    create_matrix(&grad, m+parameters->n_poly_terms, m+parameters->n_poly_terms);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++)
            grad[i][j] = parameters->phs_degree*(parameters->phs_degree-2)*pow((m1->points[cloud[i]].coords[1] - m1->points[cloud[j]].coords[1]),2) * calculate_phs_rbf(&m1->points[cloud[i]], &m1->points[cloud[j]].coords[0], parameters->phs_degree-4, parameters->dimension) + parameters->phs_degree * calculate_phs_rbf(&m1->points[cloud[i]], &m1->points[cloud[j]].coords[0], parameters->phs_degree-2, parameters->dimension);
    }
    // add polynomial terms
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < parameters->n_poly_terms; j++) {
            if (fabs(m1->points[cloud[i]].coords[0] - 0)>1e-5)
                grad[i][j+m] = pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]);
            if (fabs(m1->points[cloud[i]].coords[1] - 0)>1e-5)
                grad[i][j+m] = grad[i][j+m] * parameters->pow_y[j]*(parameters->pow_y[j]-1)*pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]-2);
            if (fabs(m1->points[cloud[i]].coords[2] - 0)>1e-5)
                {printf("%f\n", pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]));
                grad[i][j+m] = grad[i][j+m] * pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]);}
            grad[j+m][i];
        }
    }
    return grad;
}

double** dzz_matrix(mesh* m1, int* cloud, int n, struct parameters *parameters) {
    int m = parameters->n_cloud_points;
    if (n<m)
        m = n;
    double **grad;
    create_matrix(&grad, m+parameters->n_poly_terms, m+parameters->n_poly_terms);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++)
            grad[i][j] = parameters->phs_degree*(parameters->phs_degree-2)*pow((m1->points[cloud[i]].coords[2] - m1->points[cloud[j]].coords[2]),2) * calculate_phs_rbf(&m1->points[cloud[i]], &m1->points[cloud[j]].coords[0], parameters->phs_degree-4, parameters->dimension) + parameters->phs_degree * calculate_phs_rbf(&m1->points[cloud[i]], &m1->points[cloud[j]].coords[0], parameters->phs_degree-2, parameters->dimension);
    }
    // add polynomial terms
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < parameters->n_poly_terms; j++) {
            if (fabs(m1->points[cloud[i]].coords[0] - 0)>1e-5)
                grad[i][j+m] = pow(m1->points[cloud[i]].coords[0], parameters->pow_x[j]);
            if (fabs(m1->points[cloud[i]].coords[1] - 0)>1e-5)
                grad[i][j+m] = grad[i][j+m] * pow(m1->points[cloud[i]].coords[1], parameters->pow_y[j]);
            if (fabs(m1->points[cloud[i]].coords[2] - 0)>1e-5)
                {printf("%f\n", pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]));
                grad[i][j+m] = grad[i][j+m] * parameters->pow_z[j]*(parameters->pow_z[j]-1)*pow(m1->points[cloud[i]].coords[2], parameters->pow_z[j]-2);}
            grad[j+m][i] = grad[i][j+m];
        }
    }
    return grad;
}



// Following function creates the full derivative matrices for the mesh

double** create_full_gradx_matrix(mesh* m1, struct parameters *parameters) {
    double **grad, **gradx, **A, **A_inv, **B1;
    create_matrix(&grad, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&A, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&A_inv, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&B1, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&gradx, *m1->n_count, *m1->n_count);

    for (int i = 0; i < *m1->n_count; i++) {
        if (m1->corner_tag[i] == 0){
            A = create_A_matrix_from_cloud_indices(m1, m1->cloud_index[i], parameters->n_cloud_points, parameters);
            grad = gradx_matrix(m1, m1->cloud_index[i], parameters->n_cloud_points, parameters);
            matrixInverse_Gauss_Jordan(A, A_inv, parameters->n_cloud_points+parameters->n_poly_terms);
            multiply_matrices(grad, A_inv, B1, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
            for (int j = 0; j < parameters->n_cloud_points; j++) {
                gradx[i][m1->cloud_index[i][j]] = B1[0][j];
            }
        }
    }
    free_matrix(A, parameters->n_cloud_points+parameters->n_poly_terms);
    free_matrix(grad, parameters->n_cloud_points+parameters->n_poly_terms);
    free_matrix(A_inv, parameters->n_cloud_points+parameters->n_poly_terms);
    free_matrix(B1, parameters->n_cloud_points+parameters->n_poly_terms);
    return gradx;
}

double** create_full_grady_matrix(mesh* m1, struct parameters *parameters) {
    double **grad, **grady, **A, **A_inv, **B1;
    create_matrix(&grad, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&A, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&A_inv, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&B1, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&grady, *m1->n_count, *m1->n_count);
    for (int i = 0; i < *m1->n_count; i++) {
        if (m1->corner_tag[i] == 0){
            A = create_A_matrix_from_cloud_indices(m1, m1->cloud_index[i], parameters->n_cloud_points, parameters);
            grad = grady_matrix(m1, m1->cloud_index[i], parameters->n_cloud_points, parameters);
            matrixInverse_Gauss_Jordan(A, A_inv, parameters->n_cloud_points+parameters->n_poly_terms);
            multiply_matrices(grad, A_inv, B1, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
            for (int j = 0; j < parameters->n_cloud_points; j++) {
                grady[i][m1->cloud_index[i][j]] = B1[0][j];
            }
        }
    }
    free_matrix(A, parameters->n_cloud_points+parameters->n_poly_terms);
    free_matrix(grad, parameters->n_cloud_points+parameters->n_poly_terms);
    free_matrix(A_inv, parameters->n_cloud_points+parameters->n_poly_terms);
    free_matrix(B1, parameters->n_cloud_points+parameters->n_poly_terms);
    return grady;
}

double** create_full_gradz_matrix(mesh* m1, struct parameters *parameters) {
    double **grad, **gradz, **A, **A_inv, **B1;
    create_matrix(&grad, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&A, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&A_inv, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&B1, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&gradz, *m1->n_count, *m1->n_count);
    for (int i = 0; i < *m1->n_count; i++) {
        if (m1->corner_tag[i] == 0){
            A = create_A_matrix_from_cloud_indices(m1, m1->cloud_index[i], parameters->n_cloud_points, parameters);
            grad = gradz_matrix(m1, m1->cloud_index[i], parameters->n_cloud_points, parameters);
            matrixInverse_Gauss_Jordan(A, A_inv, parameters->n_cloud_points+parameters->n_poly_terms);
            multiply_matrices(grad, A_inv, B1, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
            for (int j = 0; j < parameters->n_cloud_points; j++) {
                gradz[i][m1->cloud_index[i][j]] = B1[0][j];
            }
        }
    }
    free_matrix(A, parameters->n_cloud_points+parameters->n_poly_terms);
    free_matrix(grad, parameters->n_cloud_points+parameters->n_poly_terms);
    free_matrix(A_inv, parameters->n_cloud_points+parameters->n_poly_terms);
    free_matrix(B1, parameters->n_cloud_points+parameters->n_poly_terms);
    return gradz;
}

double** create_full_laplacian_matrix(mesh* m1, struct parameters *parameters) {
    double **lap, **laplacian, **A, **A_inv, **B1;
    create_matrix(&lap, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&A, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&A_inv, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&B1, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
    create_matrix(&laplacian, *m1->n_count, *m1->n_count);
    for (int i = 0; i < *m1->n_count; i++) {
        if (m1->corner_tag[i] == 0){
            A = create_A_matrix_from_cloud_indices(m1, m1->cloud_index[i], parameters->n_cloud_points, parameters);
            lap = laplacian_matrix(m1, m1->cloud_index[i], parameters->n_cloud_points, parameters);
            matrixInverse_Gauss_Jordan(A, A_inv, parameters->n_cloud_points+parameters->n_poly_terms);
            multiply_matrices(lap, A_inv, B1, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms, parameters->n_cloud_points+parameters->n_poly_terms);
            for (int j = 0; j < parameters->n_cloud_points; j++) {
                laplacian[i][m1->cloud_index[i][j]] = B1[0][j];
            }
        }
    }
    free_matrix(A, parameters->n_cloud_points+parameters->n_poly_terms);
    free_matrix(lap, parameters->n_cloud_points+parameters->n_poly_terms);
    free_matrix(A_inv, parameters->n_cloud_points+parameters->n_poly_terms);
    free_matrix(B1, parameters->n_cloud_points+parameters->n_poly_terms);
    return laplacian;
}
