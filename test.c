#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Wall-clock timer
static double now() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// Build 2D Poisson matrix (5-point stencil) and RHS using analytical solution
void build_poisson_problem(int N, double **A, double **b) {
    int n = N * N;
    double h = 1.0 / (N + 1);

    *A = (double*)calloc(n * n, sizeof(double));
    *b = (double*)calloc(n, sizeof(double));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int id = i * N + j;
            (*A)[id*n + id] = 4.0;

            if (i > 0) (*A)[id*n + (id - N)] = -1.0;
            if (i < N-1) (*A)[id*n + (id + N)] = -1.0;
            if (j > 0) (*A)[id*n + (id - 1)] = -1.0;
            if (j < N-1) (*A)[id*n + (id + 1)] = -1.0;

            double x = (i + 1) * h;
            double y = (j + 1) * h;
            (*b)[id] = 2.0 * M_PI * M_PI * sin(M_PI*x) * sin(M_PI*y) * h * h; // scaled by h^2
        }
    }
}

// BiCGSTAB solver
int bicgstab(int n, double *A, double *b, double *x, int max_it, double tol) {
    double *r = calloc(n, sizeof(double));
    double *r0 = calloc(n, sizeof(double));
    double *p = calloc(n, sizeof(double));
    double *v = calloc(n, sizeof(double));
    double *s = calloc(n, sizeof(double));
    double *t = calloc(n, sizeof(double));

    // r = b - A x
    for (int i = 0; i < n; i++) {
        double Ax = 0.0;
        for (int j = 0; j < n; j++) Ax += A[i*n + j] * x[j];
        r[i] = b[i] - Ax;
        r0[i] = r[i];
    }

    double rho_old = 1, alpha = 1, omega = 1;
    for (int i = 0; i < n; i++) p[i] = 0.0;

    for (int it = 0; it < max_it; it++) {
        double rho_new = 0;
        for (int i = 0; i < n; i++) rho_new += r0[i] * r[i];
        if (fabs(rho_new) < 1e-30) break;

        double beta = (rho_new / rho_old) * (alpha / omega);
        for (int i = 0; i < n; i++) p[i] = r[i] + beta * (p[i] - omega * v[i]);

        // v = A p
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) sum += A[i*n + j] * p[j];
            v[i] = sum;
        }

        double denom = 0;
        for (int i = 0; i < n; i++) denom += r0[i] * v[i];
        alpha = rho_new / denom;

        for (int i = 0; i < n; i++) s[i] = r[i] - alpha * v[i];

        double norm_s = 0;
        for (int i = 0; i < n; i++) norm_s += s[i] * s[i];
        norm_s = sqrt(norm_s);
        if (norm_s < tol) {
            for (int i = 0; i < n; i++) x[i] += alpha * p[i];
            break;
        }

        // t = A s
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) sum += A[i*n + j] * s[j];
            t[i] = sum;
        }

        double t_s = 0, t_t = 0;
        for (int i = 0; i < n; i++) {
            t_s += t[i] * s[i];
            t_t += t[i] * t[i];
        }
        omega = t_s / t_t;

        for (int i = 0; i < n; i++) x[i] += alpha * p[i] + omega * s[i];
        for (int i = 0; i < n; i++) r[i] = s[i] - omega * t[i];

        double norm_r = 0;
        for (int i = 0; i < n; i++) norm_r += r[i] * r[i];
        norm_r = sqrt(norm_r);
        if (norm_r < tol) break;

        rho_old = rho_new;
    }

    free(r); free(r0); free(p); free(v); free(s); free(t);
    return 0;
}

// Compute L2 error against analytical solution
double l2_error(double* x, int N) {
    double h = 1.0 / (N+1);
    double error = 0.0;
    int n = N*N;
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            int id = i*N + j;
            double xi = (i+1)*h;
            double yj = (j+1)*h;
            double u_exact = sin(M_PI*xi)*sin(M_PI*yj);
            double diff = x[id] - u_exact;
            error += diff*diff;
        }
    }
    return sqrt(error / n);
}

// Placeholder AMG
void AMG_solve_placeholder() {
    printf("[AMG placeholder] AMG solver not linked. Add library to enable.\n");
}

// Main
int main() {
    int N = 40; // grid size
    int n = N*N;
    double *A, *b, *x;

    build_poisson_problem(N, &A, &b);
    x = calloc(n, sizeof(double)); // initial guess = 0

    printf("Matrix size = %d x %d (unknowns = %d)\n", n, n, n);

    double t0 = now();
    bicgstab(n, A, b, x, 2000, 1e-8);
    double t1 = now();

    double err = l2_error(x, N);
    printf("L2 error = %.6e\n", err);
    printf("BiCGSTAB time: %.6f sec\n", t1 - t0);

    printf("Running AMG...\n");
    AMG_solve_placeholder();

    free(A); free(b); free(x);
    return 0;
}
