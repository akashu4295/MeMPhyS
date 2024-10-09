// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef SOLVERS_H
#define SOLVERS_H

#include "structures.h"
#include "general_functions.h"
#include "mesh_functions.h"
#include "mat_lib.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/////////////////////////////////////////////////////////////////////////////
// Function Declarations
////////////////////////////////////////////////////////////////////////////

// Fractional Step Explicit Solver Modules
void calculate_intermediate_velocity(PointStructure* myPointStruct, double* u_intermediate, double* v_intermediate, double* w_intermediate, double* u_n, double* v_n, double* w_n, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian);
void calculate_mass_residual(PointStructure* myPointStruct, double* u_intermediate, double* v_intermediate, double* w_intermediate, double** Dx, double** Dy, double** Dz, double rho, double dt, double* mass_residual);
void calculate_pressure(PointStructure* myPointStruct, double* p, double* u_intermediate, double* v_intermediate, double* w_intermediate, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian, double* mass_residual);
void update_velocity(PointStructure* myPointStruct, double* u_n, double* v_n, double* w_n, double* u_intermediate, double* v_intermediate, double* w_intermediate, double* p, double rho, double dt, double** Dx, double** Dy, double** Dz);
void fractional_step_explicit(PointStructure* myPointStruct, FieldVariables* myFieldVariables_n, double dt);
void calculate_intermediate_velocity_update(PointStructure* myPointStruct, double* u_intermediate, double* v_intermediate, double* w_intermediate,double* u_n, double* v_n, double* w_n, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian);

// Time-Implicit Solver Modules
void calculate_intermediate_velocity_implicit(PointStructure* myPointStruct, double* u_new, double* v_new, double* w_new, double* u_old, double* v_old, double* w_old, double* u_n, double* v_n, double* w_n, double* p_n, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian, int iter_momentum);
void calculate_mass_residual_implicit(PointStructure* myPointStruct, double* u_new, double* v_new, double* w_new, double** Dx, double** Dy, double** Dz, double rho, double dt, double* mass_residual);
void calculate_pressure_correction_implicit(PointStructure* myPointStruct, double* p, double* u_new, double* v_new, double* w_new, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian, double* mass_residual, int iter_pressure);
void update_velocity_implicit(PointStructure* myPointStruct, double* u_old, double* v_old, double* w_old, double* p, double* u_new, double* v_new, double* w_new, double* pprime, double rho, double dt, double** Dx, double** Dy, double** Dz);
void time_implicit_solver(PointStructure* myPointStruct, FieldVariables* myFieldVariables_n, double dt, int iter_momentum, int iter_pressure, int iter_timple);

// Poisson solver
void Poisson_solver(double* p, double* rhs, double** lhs, int** cloud_index, int num_nodes, int num_cloud_points, int num_boundary_points, int num_iter);

// Heat Conduction Solver
void restrict_residuals(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, double* res_f, double* res_c);
void prolongate_correction(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, double* corr_f, double* corr_c);
void calculate_residual(PointStructure* myPointStruct, double* T, double* source, double* res);
void relax_vcycle(PointStructure* myPointStruct, double* T, double* source, int num_iter, double omega);

/////////////////////////////////////////////////////////////////////////////
// Function Definitions
////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// Fractional Step Explicit Solver Modules
/////////////////////////////////////////////////////////////////////////////
void calculate_intermediate_velocity_update(PointStructure* myPointStruct, double* u_intermediate, double* v_intermediate, double* w_intermediate,double* u_n, double* v_n, double* w_n, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian)
{
    // Calculate the intermediate velocity
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    int num_boundary_points = myPointStruct->num_boundary_nodes;
    double* diffusion = create_vector(num_nodes);
    double* advection = create_vector(num_nodes);
    double *temp1 = create_vector(num_nodes);
    double *temp2 = create_vector(num_nodes);
    double *temp3 = create_vector(num_nodes);
    double nu = mu/rho;

    multiply_sparse_matrix_vector(Dx, u_n, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(Dy, u_n, temp2, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector(Dz, u_n, temp3, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    }
    for (int i = num_boundary_points; i < num_nodes; i++){
        if (!myPointStruct->boundary_tag[i]) 
            u_intermediate[i] = u_n[i] - dt * (u_n[i] * temp1[i] + v_n[i] * temp2[i] + w_n[i] * temp3[i]);
    }

    multiply_sparse_matrix_vector(laplacian, u_n, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    for (int i = num_boundary_points; i < num_nodes; i++){
        if (!myPointStruct->boundary_tag[i])
            u_intermediate[i] = u_intermediate[i] + dt*nu* temp1[i];
    }

//        v velocity

    multiply_sparse_matrix_vector(Dx, v_n, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(Dy, v_n, temp2, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector(Dz, v_n, temp3, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    }
    for (int i = num_boundary_points; i < num_nodes; i++){
        if (!myPointStruct->boundary_tag[i]) 
            v_intermediate[i] = v_n[i] - dt * (u_n[i] * temp1[i] + v_n[i] * temp2[i] + w_n[i] * temp3[i]);
    }

    multiply_sparse_matrix_vector(laplacian, v_n, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes,num_cloud_points);
    for (int i = num_boundary_points; i < num_nodes; i++){
        if (!myPointStruct->boundary_tag[i])
            v_intermediate[i] = v_intermediate[i] + dt * nu* temp1[i];
    }

//           w-velocity

    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector(Dx, w_n, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector(Dy, w_n, temp2, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector(Dz, w_n, temp3, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    }
    for (int i = num_boundary_points; i < num_nodes; i++){
        if (!myPointStruct->boundary_tag[i]) 
            w_intermediate[i] = w_n[i] - dt * (u_n[i] * temp1[i] + v_n[i] * temp2[i] + w_n[i] * temp3[i]);
    }

    multiply_sparse_matrix_vector(laplacian, w_n, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes,num_cloud_points);
    for (int i = num_boundary_points; i < num_nodes; i++){
        if (!myPointStruct->boundary_tag[i])
            w_intermediate[i] = w_intermediate[i] + dt * nu* temp1[i];
    }
 
    free_vector(temp1);
    free_vector(temp2);
    free_vector(temp3);
    // for (int i = 0; i < num_nodes; i++){
    //     printf("%lf %lf %lf %lf\n", u_intermediate[i], v_intermediate[i], u_n[i], v_n[i]);
    // }
}

void calculate_intermediate_velocity(PointStructure* myPointStruct, double* u_intermediate, double* v_intermediate, double* w_intermediate, double* u_n, double* v_n, double* w_n, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian)
{
    // Calculate the intermediate velocity
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    int num_boundary_points = myPointStruct->num_boundary_nodes;
    double** advection_operator = create_matrix1(num_nodes, myPointStruct->num_cloud_points);
    double* diffusion = create_vector(num_nodes);
    double** temp1 = create_matrix1(num_nodes, myPointStruct->num_cloud_points);
    double* advection = create_vector(num_nodes);

    multiply_sparse_vector_matrix(u_n, Dx, temp1, num_nodes, num_cloud_points);
    add_matrices_to_first(advection_operator, temp1, num_nodes, num_cloud_points);
    multiply_sparse_vector_matrix(v_n, Dy, temp1, num_nodes, num_cloud_points);
    add_matrices_to_first(advection_operator, temp1, num_nodes, num_cloud_points);
    if (parameters.dimension == 3){
        multiply_sparse_vector_matrix(w_n, Dz, temp1, num_nodes, num_cloud_points);
        add_matrices_to_first(advection_operator, temp1, num_nodes, num_cloud_points);
    }
    double nu = mu/rho;
    
    multiply_sparse_matrix_vector(advection_operator, u_n, advection, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(laplacian, u_n, diffusion, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    for (int i = num_boundary_points; i < num_nodes; i++){
        u_intermediate[i] = u_n[i] + dt * (-advection[i] + nu* diffusion[i]);
    }

    multiply_sparse_matrix_vector(advection_operator, v_n, advection, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(laplacian, v_n, diffusion, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    for (int i = num_boundary_points; i < num_nodes; i++){
        v_intermediate[i] = v_n[i] + dt * (-advection[i] + nu* diffusion[i]);
    }

    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector(advection_operator, w_n, advection, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector(laplacian, w_n, diffusion, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
        for (int i = num_boundary_points; i < num_nodes; i++){
            w_intermediate[i] = w_n[i] + dt * (-advection[i] + nu* diffusion[i]);
        }
    }
    
    free_matrix(advection_operator, num_nodes);
    free_vector(diffusion);
    free_matrix(temp1, num_nodes);
    free_vector(advection);

    // for (int i = 0; i < num_nodes; i++){
    //     printf("%lf %lf %lf %lf\n", u_intermediate[i], v_intermediate[i], u_n[i], v_n[i]);
    // }
}

void calculate_mass_residual(PointStructure* myPointStruct, double* u_intermediate, double*v_intermediate, double* w_intermediate, double** Dx, double** Dy, double** Dz, double rho, double dt, double* mass_residual)
{
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    int num_boundary_points = myPointStruct->num_boundary_nodes;
    double* temp1 = create_vector(num_nodes);
    double* temp2 = create_vector(num_nodes);
    multiply_sparse_matrix_vector(Dx, u_intermediate, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(Dy, v_intermediate, temp2, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    add_vectors(temp1, temp2, mass_residual, num_nodes);
    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector(Dz, w_intermediate, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
        add_vectors(mass_residual, temp1, mass_residual, num_nodes);
    }
    for (int i = 0; i < num_boundary_points; i++){
        mass_residual[i] = 0;
    }
    for (int i = num_boundary_points; i < num_nodes; i++){
        mass_residual[i] = rho*mass_residual[i]/dt;
    }
    free_vector(temp1);
    free_vector(temp2);
}

void calculate_pressure(PointStructure* myPointStruct, double* p, double* u_intermediate, double* v_intermediate, double* w_intermediate, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian, double* mass_residual)
{
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    int num_boundary_points = myPointStruct->num_boundary_nodes;

    double* ptemp = create_vector(num_nodes);
    double* temp1 = create_vector(num_nodes);
    double* temp2 = create_vector(num_nodes);
    
    // Solve the system of equations
    double sum = 0, Ap = 0;
    calculate_mass_residual(myPointStruct, u_intermediate, v_intermediate, w_intermediate, Dx, Dy, Dz, rho, dt, mass_residual);

    for (int iter = 0; iter < 100; iter++)
    {       
        for (int i = 0; i < num_nodes; i++){
            ptemp[i] = p[i];
        }
    // Loop over interior nodes
        for (int i = num_boundary_points; i < num_nodes; i++)
        {
            sum = 0;
            for (int j = 1; j < num_cloud_points; j++){
                sum += laplacian[i][j]*p[myPointStruct->cloud_index[i][j]];
            }
            p[i] = 1.0*((mass_residual[i]-sum)/laplacian[i][0]) - 0.0*p[i];
        }
        printf("L2_norm = %E\n", l2_norm(ptemp, p, num_nodes));
    
        // Loop over boundary nodes, set the Neumann Boundary conditions: dpdn = 0
        for (int i = myPointStruct->num_corners; i < num_boundary_points; i++)
        {
            sum = 0; Ap = 0;
            for (int j = 1; j < num_cloud_points; j++){
                sum += Dx[i][j]*p[myPointStruct->cloud_index[i][j]]*myPointStruct->x_normal[i];
                sum += Dy[i][j]*p[myPointStruct->cloud_index[i][j]]*myPointStruct->y_normal[i];
                if (parameters.dimension == 3){
                    sum += Dz[i][j]*p[myPointStruct->cloud_index[i][j]]*myPointStruct->z_normal[i];
                }
            }
            Ap += Dx[i][0]*myPointStruct->x_normal[i];
            Ap += Dy[i][0]*myPointStruct->y_normal[i];
            if (parameters.dimension == 3){
                Ap += Dz[i][0]*myPointStruct->z_normal[i];
            }
            p[i] = -sum/Ap;
            // p[i] = p[i];
            // printf("%lf %lf %lf\n", sum, Ap, ptemp[i]);
        
        // Update the pressure
        // printf("L2_norm = %E\n", l2_norm(ptemp, p, num_nodes));
        // for (int i = 0; i < num_nodes; i++){
        //     p[i] = ptemp[i];
        }
    }
    
    free(ptemp);
    free(temp1);
    free(temp2);
}

void update_velocity(PointStructure* myPointStruct, double* u_n, double* v_n, double* w_n, double* u_intermediate, double* v_intermediate, double* w_intermediate, double* p, double rho, double dt, double** Dx, double** Dy, double** Dz)
{
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    int num_boundary_points = myPointStruct->num_boundary_nodes;
    double* dpdx = create_vector(num_nodes);
    double* dpdy = create_vector(num_nodes);
    double* dpdz = create_vector(num_nodes);
    multiply_sparse_matrix_vector(Dx, p, dpdx, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(Dy, p, dpdy, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector(Dz, p, dpdz, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    }
    // Update Interior nodes
    for (int i = num_boundary_points; i < num_nodes; i++){
        u_n[i] = u_intermediate[i] - dt * dpdx[i]/rho;
        v_n[i] = v_intermediate[i] - dt * dpdy[i]/rho;
        if (parameters.dimension == 3){
            w_n[i] = w_intermediate[i] - dt * dpdz[i]/rho;
        }
    }
    // Update Boundary nodes.. Needs to determine the boundary conditions depending on the problem

}

void fractional_step_explicit(PointStructure* myPointStruct, FieldVariables* myFieldVariables_n, double dt)
{   
    double* u_intermediate = create_vector(myPointStruct->num_nodes);
    double* v_intermediate = create_vector(myPointStruct->num_nodes);
    double* w_intermediate = create_vector(myPointStruct->num_nodes);
    double* mass_residual = create_vector(myPointStruct->num_nodes);
    calculate_intermediate_velocity(myPointStruct, u_intermediate, v_intermediate, w_intermediate, myFieldVariables_n->u, myFieldVariables_n->v, myFieldVariables_n->w, myFieldVariables_n->rho, myFieldVariables_n->mu, dt, myPointStruct->Dx, myPointStruct->Dy, myPointStruct->Dz, myPointStruct->lap);

    calculate_pressure(myPointStruct, myFieldVariables_n->p, u_intermediate, v_intermediate, w_intermediate, myFieldVariables_n->rho, myFieldVariables_n->mu, dt, myPointStruct->Dx, myPointStruct->Dy, myPointStruct->Dz, myPointStruct->lap, mass_residual);

    update_velocity(myPointStruct, myFieldVariables_n->u, myFieldVariables_n->v, myFieldVariables_n->w, u_intermediate, v_intermediate, w_intermediate, myFieldVariables_n->p, myFieldVariables_n->rho, dt, myPointStruct->Dx, myPointStruct->Dy, myPointStruct->Dz);
    free(u_intermediate);
    free(v_intermediate);
    free(w_intermediate);
    free(mass_residual);
}

////////////////////////////////////////////////////////////////////////////////////
// Time-Implicit Solver Modules
////////////////////////////////////////////////////////////////////////////////////
void calculate_intermediate_velocity_implicit(PointStructure* myPointStruct, 
double* u_new, double* v_new, double* w_new, double* u_old, double* v_old, 
double* w_old, double *u_n, double *v_n, double *w_n,double* p_n, 
double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian, int iter_momentum){
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    int num_boundary_points = myPointStruct->num_boundary_nodes;
    
    double nu = mu/rho;
    double denom = 0;
    double temp1, temp2, temp3, temp4;
    double *dpdx, *dpdy, *dpdz;
    dpdx = create_vector(num_nodes);
    dpdy = create_vector(num_nodes);

    multiply_sparse_matrix_vector(Dx, p_n, dpdx, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(Dy, p_n, dpdy, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    if (parameters.dimension == 3){
        dpdz = create_vector(num_nodes);
        multiply_sparse_matrix_vector(Dz, p_n, dpdz, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    }
    for (int iter = 0; iter<iter_momentum; iter++){
        for (int i = num_boundary_points; i < num_nodes; i++){
        // u-velocity
            temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0;
            
            denom  = rho* (u_old[i] * Dx[i][0] + v_old[i] * Dy[i][0]) 
                                - mu * laplacian [i][0] + rho / dt ;
            for (int j = 1; j < num_cloud_points; j++){
                temp1 += Dx[i][j] * u_new[myPointStruct->cloud_index[i][j]];
                temp2 += Dy[i][j] * u_new[myPointStruct->cloud_index[i][j]];
                temp4 += laplacian[i][j] * u_new[myPointStruct->cloud_index[i][j]];
            }

            if (parameters.dimension == 3){
                denom += rho * w_old[i] * Dz[i][0];
                for (int j = 1; j < num_cloud_points; j++){
                    temp3 += Dz[i][j] * u_new[myPointStruct->cloud_index[i][j]];
                    temp4 += laplacian[i][j] * u_new[myPointStruct->cloud_index[i][j]];
                }
            }
            
            u_new[i] = (rho * (u_n[i]/dt -  (u_old[i] * temp1 + v_old[i] * temp2 
                                        + w_old[i] * temp3)) + mu*temp4 - dpdx[i]) / denom;

        // v-velocity
            temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0;
            
            denom = rho* (u_old[i] * Dx[i][0] + v_old[i] * Dy[i][0]) - mu * laplacian [i][0] + rho / dt ;
            for (int j = 1; j < num_cloud_points; j++){
                temp1 += Dx[i][j] * v_new[myPointStruct->cloud_index[i][j]];
                temp2 += Dy[i][j] * v_new[myPointStruct->cloud_index[i][j]];
                temp4 += laplacian[i][j] * v_new[myPointStruct->cloud_index[i][j]];
            }

            if (parameters.dimension == 3){
                denom += rho * w_old[i] * Dz[i][0];
                for (int j = 1; j < num_cloud_points; j++){
                    temp3 += Dz[i][j] * v_new[myPointStruct->cloud_index[i][j]];
                    temp4 += laplacian[i][j] * v_new[myPointStruct->cloud_index[i][j]];
                }
            }

            v_new[i] = (rho * (v_n[i]/dt -  (u_old[i] * temp1 + v_old[i] * temp2 
                                        + w_old[i] * temp3)) + mu*temp4 - dpdy[i]) / denom;

        // w-velocity
            if (parameters.dimension == 3){
                temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0;
                
                denom  = rho* (u_old[i] * Dx[i][0] + v_old[i] * Dy[i][0] + w_old[i] * Dz[i][0]) - mu * laplacian [i][0] + rho / dt;
                for (int j = 1; j < num_cloud_points; j++){
                    temp1 += Dx[i][j] * w_new[myPointStruct->cloud_index[i][j]];
                    temp2 += Dy[i][j] * w_new[myPointStruct->cloud_index[i][j]];
                    temp3 += Dz[i][j] * w_new[myPointStruct->cloud_index[i][j]];
                    temp4 += laplacian[i][j] * w_new[myPointStruct->cloud_index[i][j]];
                }

                w_new[i] = (rho * (w_n[i]/dt -  (u_old[i] * temp1 + v_old[i] * temp2 
                                            + w_old[i] * temp3)) + mu*temp4 - dpdz[i]) / denom;
            }
        }
    }
    free(dpdx);
    free(dpdy);
    if (parameters.dimension == 3){
        free(dpdz);
    }
}


void calculate_mass_residual_implicit(PointStructure* myPointStruct, 
            double* u_new, double* v_new, double* w_new, double** Dx, double** Dy, double** Dz,
            double rho, double dt, double* mass_residual){
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    int num_boundary_points = myPointStruct->num_boundary_nodes;
    double* temp1 = create_vector(num_nodes);
    double* temp2 = create_vector(num_nodes);
    multiply_sparse_matrix_vector(Dx, u_new, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(Dy, v_new, temp2, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    add_vectors(temp1, temp2, mass_residual, num_nodes);
    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector(Dz, w_new, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
        add_vectors(mass_residual, temp1, mass_residual, num_nodes);
    }
    for (int i = 0; i < num_boundary_points; i++){
        mass_residual[i] = 0;
    }
    for (int i = num_boundary_points; i < num_nodes; i++){
        mass_residual[i] = rho*mass_residual[i]/dt;
    }
    free(temp1);
    free(temp2);
}

void calculate_pressure_correction_implicit(PointStructure* myPointStruct, double* pprime, double* u_new, double* v_new, double* w_new, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian, double* mass_residual, int iter_pressure){
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    int num_boundary_points = myPointStruct->num_boundary_nodes;

    double* ptemp = create_vector(num_nodes);
    double* temp1 = create_vector(num_nodes);
    double* temp2 = create_vector(num_nodes);
    
    // Solve the system of equations
    double sum = 0, Ap = 0;
    
    calculate_mass_residual_implicit(myPointStruct, u_new, v_new, w_new, Dx, Dy, Dz, rho, dt, mass_residual);

    for (int iter = 0; iter < iter_pressure; iter++)
    {       
        for (int i = 0; i < num_nodes; i++){
            ptemp[i] = pprime[i]; // for L2 norm calculation
        }
    // Loop over interior nodes
        for (int i = num_boundary_points; i < num_nodes; i++)
        {
            sum = 0;
            for (int j = 1; j < num_cloud_points; j++){
                sum += laplacian[i][j]*pprime[myPointStruct->cloud_index[i][j]];
            }
            pprime[i] = 1.0*((mass_residual[i]-sum)/laplacian[i][0]) - 0.0*pprime[i];
        }
        printf("L2_norm = %E\n", l2_norm(ptemp, pprime, num_nodes));
    
        // Loop over boundary nodes, set the Neumann Boundary conditions: dpdn = 0
        for (int i = myPointStruct->num_corners; i < num_boundary_points; i++)
        {
            sum = 0; Ap = 0;
            for (int j = 1; j < num_cloud_points; j++){
                sum += Dx[i][j]*pprime[myPointStruct->cloud_index[i][j]]*myPointStruct->x_normal[i];
                sum += Dy[i][j]*pprime[myPointStruct->cloud_index[i][j]]*myPointStruct->y_normal[i];
                if (parameters.dimension == 3){
                    sum += Dz[i][j]*pprime[myPointStruct->cloud_index[i][j]]*myPointStruct->z_normal[i];
                }
            }
            Ap += Dx[i][0]*myPointStruct->x_normal[i];
            Ap += Dy[i][0]*myPointStruct->y_normal[i];
            if (parameters.dimension == 3){
                Ap += Dz[i][0]*myPointStruct->z_normal[i];
            }
            pprime[i] = -sum/Ap;
        }
    }

    free(ptemp);
    free(temp1);
    free(temp2);
}

void update_velocity_implicit(PointStructure* myPointStruct,
 double* u_old, double* v_old, double* w_old, double* p, 
 double* u_new, double* v_new, double* w_new, double* pprime, 
 double rho, double dt, double** Dx, double** Dy, double** Dz){
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    int num_boundary_points = myPointStruct->num_boundary_nodes;
    double* dpdx = create_vector(num_nodes);
    double* dpdy = create_vector(num_nodes);
    double* dpdz = create_vector(num_nodes);
    multiply_sparse_matrix_vector(Dx, pprime, dpdx, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(Dy, pprime, dpdy, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector(Dz, pprime, dpdz, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    }
    // Update Interior nodes
    for (int i = num_boundary_points; i < num_nodes; i++){
        u_old[i] = u_new[i] - dt * dpdx[i]/rho;
        v_old[i] = v_new[i] - dt * dpdy[i]/rho;
        if (parameters.dimension == 3){
            w_old[i] = w_new[i] - dt * dpdz[i]/rho;
        }
        p[i] = p[i] + pprime[i];
        u_new[i] = u_old[i];
        v_new[i] = v_old[i];
        if (parameters.dimension == 3){
            w_new[i] = w_old[i];
        }
    }
    for (int i = 0; i < num_boundary_points; i++){
        p[i] = p[i] + pprime[i];
    }
}


void time_implicit_solver(PointStructure* myPointStruct, FieldVariables* myFieldVariables_n, double dt, int iter_momentum, int iter_pressure, int iter_timple){
    double* u_new = create_vector(myPointStruct->num_nodes);
    double* v_new = create_vector(myPointStruct->num_nodes);
    double* w_new = create_vector(myPointStruct->num_nodes);
    double* u_old = create_vector(myPointStruct->num_nodes);
    double* v_old = create_vector(myPointStruct->num_nodes);
    double* w_old = create_vector(myPointStruct->num_nodes);
    double* pprime = create_vector(myPointStruct->num_nodes);
    double* mass_residual = create_vector(myPointStruct->num_nodes);
    
    for (int iter = 0; iter<iter_timple; iter++){
        calculate_intermediate_velocity_implicit(myPointStruct, u_new, v_new, w_new, 
                        u_old, v_old, w_old, myFieldVariables_n->u, 
                        myFieldVariables_n->v, myFieldVariables_n->w, 
                        myFieldVariables_n->p, myFieldVariables_n->rho, 
                        myFieldVariables_n->mu, dt, myPointStruct->Dx, 
                        myPointStruct->Dy, myPointStruct->Dz, myPointStruct->lap, 
                        iter_momentum);
        calculate_pressure_correction_implicit(myPointStruct, pprime, u_new, v_new, w_new, myFieldVariables_n->rho, myFieldVariables_n->mu, dt, myPointStruct->Dx, myPointStruct->Dy, myPointStruct->Dz, myPointStruct->lap, mass_residual, iter_pressure);
        update_velocity_implicit(myPointStruct, u_old, v_old, w_old, myFieldVariables_n->p, u_new, v_new, w_new, pprime, myFieldVariables_n->rho, dt, myPointStruct->Dx, myPointStruct->Dy, myPointStruct->Dz);
    }
    copy_vector(u_old, myFieldVariables_n->u, myPointStruct->num_nodes);
    copy_vector(v_old, myFieldVariables_n->v, myPointStruct->num_nodes);
    copy_vector(w_old, myFieldVariables_n->w, myPointStruct->num_nodes);

    free(u_new);
    free(v_new);
    free(w_new);
    free(mass_residual);
    free(u_old);
    free(v_old);
    free(w_old);
    free(pprime);
}

///////////////////////////////////////////////////////////////////////////////

void Poisson_solver(double* p, double* rhs, double** lhs, int** cloud_index, int num_nodes, int num_cloud_points, int num_boundary_points, int num_iter){
    double* ptemp = create_vector(num_nodes);
    double* temp = create_vector(num_nodes);
    double sum = 0;
    for (int iter = 0; iter < num_iter; iter++){
        for (int i = 0; i < num_nodes; i++){
            ptemp[i] = p[i];
        }
        for (int i = num_boundary_points; i < num_nodes; i++){
            sum = 0;
            for (int j = 1; j < num_cloud_points; j++){
                sum += lhs[i][j]*p[cloud_index[i][j]];
            }
            p[i] = 1.6*((rhs[i]-sum)/lhs[i][0]) - 0.6*p[i];
        }
    }
    free(ptemp);
    free(temp);
}

void restrict_residuals(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, 
                                                    double *res_f, double* res_c){
    double results = 0.0;
    int i_restr;
    for (int i = 0; i < mypointStruct_c->num_nodes; i++){
        if (mypointStruct_c->boundary_tag[i]==false){
            i_restr = mypointStruct_c->restriction_points[i];
            for (int j = 0; j < mypointStruct_f->num_cloud_points; j++){
                results += mypointStruct_c->restr_mat[i][j] * res_f[mypointStruct_f->cloud_index[i_restr][j]];
            }
        }
        res_c[i] = results;
        results = 0.0;
    }
}

void prolongate_correction(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, 
                                                    double *cor_f, double* cor_c){
    double results = 0.0;
    int i_prol;
    for (int i = 0; i < mypointStruct_f->num_nodes; i++){
        if (mypointStruct_f->boundary_tag[i]==false){
            i_prol = mypointStruct_f->prolongation_points[i];
            for (int j = 0; j < mypointStruct_c->num_cloud_points; j++){
                results += mypointStruct_f->prol_mat[i][j] * cor_c[mypointStruct_c->cloud_index[i_prol][j]];
            }
        }
        cor_f[i] = cor_f[i] + results;
        results = 0.0;
    }
    for (int i = 0; i<mypointStruct_c->num_nodes; i++){
        if (mypointStruct_c->boundary_tag[i]==false){
            cor_c[i] = 0.0;
        }
    }
}

void calculate_residual(PointStructure* mypointStruct, double* res, double* p, double* rhs){
    int num_nodes = mypointStruct->num_nodes;
    int num_cloud_points = mypointStruct->num_cloud_points;
    int num_boundary_points = mypointStruct->num_boundary_nodes;
    double sum = 0;
    for (int i = num_boundary_points; i < num_nodes; i++){
        sum = 0;
        for (int j = 0; j < num_cloud_points; j++){
            sum += mypointStruct->lap[i][j]*p[mypointStruct->cloud_index[i][j]];
        }
        res[i] = rhs[i] - sum;
    }
}

void relax_vcycle(PointStructure* mypointstruct, double* T, double* rhs, int num_iter, double omega){
    int num_nodes = mypointstruct->num_nodes;
    int num_cloud_points = mypointstruct->num_cloud_points;
    int num_boundary_points = mypointstruct->num_boundary_nodes;
    double* Ttemp = create_vector(num_nodes);
    double sum = 0;

    for (int iter = 0; iter < num_iter; iter++){
        for (int i = 0; i < num_nodes; i++){
            Ttemp[i] = T[i];
        }
        for (int i = num_boundary_points; i < num_nodes; i++){
            sum = 0;
            for (int j = 1; j < num_cloud_points; j++){
                sum += mypointstruct->lap[i][j]*T[mypointstruct->cloud_index[i][j]];
            }
            T[i] = omega*((rhs[i]-sum)/mypointstruct->lap[i][0]) + (1-omega)*T[i];
        }
    }
    free(Ttemp);
}
#endif