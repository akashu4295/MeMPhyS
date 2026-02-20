// void FS_relaxation_vectorised_BiCGStab(PointStructure* ps, const double* b, double* x, int max_iter, double tol)
// {
//     int N = ps->num_nodes;
//     int n = ps->num_cloud_points;

//     // Extract diagonal for preconditioning
//     double *diag_inv = malloc(N * sizeof(double));
//     for (int i = 0; i < N; i++) {
//         double diag = ps->lap_Poison[i*n + 0];  // Assuming diagonal is at j=0
//         if (fabs(diag) < 1e-14) diag = 1.0;     // Avoid division by zero
//         diag_inv[i] = 1.0 / diag;
//     }

//     // Allocate temporary vectors
//     double *r = malloc(N * sizeof(double));
//     double *r0 = malloc(N * sizeof(double));
//     double *p = malloc(N * sizeof(double));
//     double *v = malloc(N * sizeof(double));
//     double *s = malloc(N * sizeof(double));
//     double *t = malloc(N * sizeof(double));
//     double *z = malloc(N * sizeof(double));  // For preconditioning
//     double *y = malloc(N * sizeof(double));  // For preconditioning

//     // Compute initial residual: r = b - A*x
//     for (int i = 0; i < N; i++) {
//         double sum = 0.0;
//         for (int j = 0; j < n; j++) {
//             int col = ps->cloud_index[i*n + j];
//             sum += ps->lap_Poison[i*n + j] * x[col];
//         }
//         r[i] = b[i] - sum;
//         r0[i] = r[i];
//         p[i] = 0.0;
//         v[i] = 0.0;
//     }

//     // Compute initial residual norm
//     double norm_r0 = 0.0;
//     for (int i = 0; i < N; i++) {
//         norm_r0 += r[i] * r[i];
//     }
//     norm_r0 = sqrt(norm_r0);
//     // printf("BiCGStab initial residual = %e\n", norm_r0);

//     double rho = 1.0, alpha = 1.0, omega = 1.0;
//     double rho_new;

//     for (int iter = 0; iter < max_iter; iter++) {

//         // rho_new = r0^T * r
//         rho_new = 0.0;
//         for (int i = 0; i < N; i++) {
//             rho_new += r0[i] * r[i];
//         }

//         if (fabs(rho_new) < 1e-30) {
//             printf("BiCGStab: rho too small, breaking at iter %d\n", iter);
//             break;
//         }

//         // beta = (rho_new/rho) * (alpha/omega)
//         double beta = (rho_new / rho) * (alpha / omega);
        
//         // p = r + beta*(p - omega*v)
//         for (int i = 0; i < N; i++) {
//             p[i] = r[i] + beta * (p[i] - omega * v[i]);
//         }

//         // Precondition: z = M^-1 * p (Jacobi: z = diag^-1 * p)
//         for (int i = 0; i < N; i++) {
//             z[i] = diag_inv[i] * p[i];
//         }

//         // v = A*z
//         for (int i = 0; i < N; i++) {
//             double sum = 0.0;
//             for (int j = 0; j < n; j++) {
//                 int col = ps->cloud_index[i*n + j];
//                 sum += ps->lap_Poison[i*n + j] * z[col];
//             }
//             v[i] = sum;
//         }

//         // alpha = rho_new / (r0^T * v)
//         double r0v = 0.0;
//         for (int i = 0; i < N; i++) {
//             r0v += r0[i] * v[i];
//         }

//         if (fabs(r0v) < 1e-30) {
//             printf("BiCGStab: r0v too small, breaking at iter %d\n", iter);
//             break;
//         }

//         alpha = rho_new / r0v;

//         // s = r - alpha*v
//         for (int i = 0; i < N; i++) {
//             s[i] = r[i] - alpha * v[i];
//         }

//         // Precondition: y = M^-1 * s
//         for (int i = 0; i < N; i++) {
//             y[i] = diag_inv[i] * s[i];
//         }

//         // t = A*y
//         for (int i = 0; i < N; i++) {
//             double sum = 0.0;
//             for (int j = 0; j < n; j++) {
//                 int col = ps->cloud_index[i*n + j];
//                 sum += ps->lap_Poison[i*n + j] * y[col];
//             }
//             t[i] = sum;
//         }

//         // omega = (t^T * s) / (t^T * t)
//         double ts = 0.0, tt = 0.0;
//         for (int i = 0; i < N; i++) {
//             ts += t[i] * s[i];
//             tt += t[i] * t[i];
//         }

//         if (fabs(tt) < 1e-30) {
//             printf("BiCGStab: tt too small, breaking at iter %d\n", iter);
//             break;
//         }

//         omega = ts / tt;

//         // x = x + alpha*z + omega*y
//         for (int i = 0; i < N; i++) {
//             x[i] += alpha * z[i] + omega * y[i];
//         }

//         // r = s - omega*t
//         for (int i = 0; i < N; i++) {
//             r[i] = s[i] - omega * t[i];
//         }

//         // Check convergence
//         double norm_r = 0.0;
//         for (int i = 0; i < N; i++) {
//             norm_r += r[i] * r[i];
//         }
//         norm_r = sqrt(norm_r);

//         if (norm_r / norm_r0 < tol) {
//             printf("BiCGStab converged at iter %d, relative residual = %e\n", 
//                    iter, norm_r/norm_r0);
//             break;
//         }

//         rho = rho_new;
        
//         // Print progress
//         // if (iter % (max_iter-1) == 0) {
//         //     printf("BiCGStab iter %d, relative residual = %e\n", iter, norm_r/norm_r0);
//         // }
//     }

//     free(r); free(r0); free(p); free(v); free(s); free(t);
//     free(z); free(y); free(diag_inv);
// }

// int BiCGStab_Solve(PointStructure* ps, const double* b, double* x, int max_iter, double tol)
// {
//     int N = ps->num_nodes;
//     int n = ps->num_cloud_points;

//     // Allocate temporary vectors
//     double *r = malloc(N * sizeof(double));
//     double *r0 = malloc(N * sizeof(double));
//     double *p = malloc(N * sizeof(double));
//     double *v = malloc(N * sizeof(double));
//     double *s = malloc(N * sizeof(double));
//     double *t = malloc(N * sizeof(double));

//     // Compute initial residual: r = b - A*x
//     for (int i = 0; i < N; i++) {
//         double sum = 0.0;
//         for (int j = 0; j < n; j++) {
//             int col = ps->cloud_index[i*n + j];
//             sum += ps->lap_Poison[i*n + j] * x[col];
//         }
//         r[i] = b[i] - sum;
//         r0[i] = r[i];
//         p[i] = 0.0;
//         v[i] = 0.0;
//     }

//     double rho = 1.0, alpha = 1.0, omega = 1.0;
//     double rho_new;

//     for (int iter = 0; iter < max_iter; iter++) {

//         // rho_new = r0^T * r
//         rho_new = 0.0;
//         for (int i = 0; i < N; i++) {
//             rho_new += r0[i] * r[i];
//         }

//         if (fabs(rho_new) < 1e-30) {
//             printf("BiCGStab: rho too small, breaking at iter %d\n", iter);
//             break;
//         }

//         // beta = (rho_new/rho) * (alpha/omega)
//         double beta = (rho_new / rho) * (alpha / omega);
        
//         // p = r + beta*(p - omega*v)
//         for (int i = 0; i < N; i++) {
//             p[i] = r[i] + beta * (p[i] - omega * v[i]);
//         }

//         // v = A*p
//         for (int i = 0; i < N; i++) {
//             double sum = 0.0;
//             for (int j = 0; j < n; j++) {
//                 int col = ps->cloud_index[i*n + j];
//                 sum += ps->lap_Poison[i*n + j] * p[col];
//             }
//             v[i] = sum;
//         }

//         // alpha = rho_new / (r0^T * v)
//         double r0v = 0.0;
//         for (int i = 0; i < N; i++) {
//             r0v += r0[i] * v[i];
//         }

//         if (fabs(r0v) < 1e-30) {
//             printf("BiCGStab: r0v too small, breaking at iter %d\n", iter);
//             break;
//         }

//         alpha = rho_new / r0v;

//         // s = r - alpha*v
//         for (int i = 0; i < N; i++) {
//             s[i] = r[i] - alpha * v[i];
//         }

//         // Check if |s| small enough → converged early
//         double norm_s = 0.0;
//         for (int i = 0; i < N; i++) {
//             norm_s += s[i] * s[i];
//         }

//         if (sqrt(norm_s) < tol) {
//             for (int i = 0; i < N; i++) {
//                 x[i] += alpha * p[i];
//             }
//             printf("BiCGStab converged early at iter %d, norm_s = %e\n", iter, sqrt(norm_s));
//             break;
//         }

//         // t = A*s
//         for (int i = 0; i < N; i++) {
//             double sum = 0.0;
//             for (int j = 0; j < n; j++) {
//                 int col = ps->cloud_index[i*n + j];
//                 sum += ps->lap_Poison[i*n + j] * s[col];
//             }
//             t[i] = sum;
//         }

//         // omega = (t^T * s) / (t^T * t)
//         double ts = 0.0, tt = 0.0;
//         for (int i = 0; i < N; i++) {
//             ts += t[i] * s[i];
//             tt += t[i] * t[i];
//         }

//         if (fabs(tt) < 1e-30) {
//             printf("BiCGStab: tt too small, breaking at iter %d\n", iter);
//             break;
//         }

//         omega = ts / tt;

//         // x = x + alpha*p + omega*s
//         for (int i = 0; i < N; i++) {
//             x[i] += alpha * p[i] + omega * s[i];
//         }

//         // r = s - omega*t
//         for (int i = 0; i < N; i++) {
//             r[i] = s[i] - omega * t[i];
//         }

//         // Check convergence
//         double norm_r = 0.0;
//         for (int i = 0; i < N; i++) {
//             norm_r += r[i] * r[i];
//         }

//         if (sqrt(norm_r) < tol) {
//             printf("BiCGStab converged at iter %d, residual = %e\n", iter, sqrt(norm_r));
//             break;
//         }

//         rho = rho_new;
        
//         // Optional: print progress
//         if (iter % 100 == 0) {
//             printf("BiCGStab iter %d, residual = %e\n", iter, sqrt(norm_r));
//         }
//     }

//     free(r);
//     free(r0);
//     free(p);
//     free(v);
//     free(s);
//     free(t);

//     return 0; // success
// }