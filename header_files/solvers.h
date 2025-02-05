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
#include "../init.c"

/////////////////////////////////////////////////////////////////////////////
// Function Declarations
////////////////////////////////////////////////////////////////////////////

// Fractional Step Explicit Solver Modules
void update_velocity(PointStructure* myPointStruct, double* u_n, double* v_n, double* w_n, double* u_intermediate, double* v_intermediate, double* w_intermediate, double* p, double rho, double dt, double** Dx, double** Dy, double** Dz);
void fractional_step_explicit(PointStructure* myPointStruct, FieldVariables* myFieldVariables_n, double dt, int iter_pressure, double omega, bool neumann_flag_boundary);
void calculate_intermediate_velocity(PointStructure* myPointStruct, double* u_intermediate, double* v_intermediate, double* w_intermediate,double* u_n, double* v_n, double* w_n, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian);

// Poisson solver
void solve_poisson_equation(PointStructure* myPointStruct, double* P, double* source, double* dpdn, int iter_pressure, double omega, bool neumann_flag_boundary);
void dpdn_from_momentum_equation(double* dpdn, PointStructure* myPointStruct, double* u_n, double* v_n, double* w_n, double* u_intermediate, double* v_intermediate, double* w_intermediate, double rho, double dt);
void dpdn_from_complete_momentum_equation(double* dpdn, PointStructure* myPointStruct, FieldVariables* field, double* u_intermediate, double* v_intermediate, double* w_intermediate, double dt);
void calculate_mass_residual(PointStructure* myPointStruct, double* u_new, double* v_new, double* w_new, double rho, double dt, double* mass_residual);

// Heat Conduction Solver
void restrict_residuals(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, double* res_f, double* res_c);
void prolongate_correction(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, double* corr_f, double* corr_c);
void calculate_residual(PointStructure* myPointStruct, double* T, double* source, double* res);
void relax_vcycle(PointStructure* myPointStruct, double* T, double* source, int num_iter, double omega, bool neumann_flag_boundary);

// Time-Implicit Solver Modules
void calculate_intermediate_velocity_implicit(PointStructure* myPointStruct, double* u_new, double* v_new, double* w_new, double* u_old, double* v_old, double* w_old, double* u_n, double* v_n, double* w_n, double* p_n, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian, int iter_momentum);
void update_velocity_implicit(PointStructure* myPointStruct, double* u_old, double* v_old, double* w_old, double* p, double* u_new, double* v_new, double* w_new, double* pprime, double rho, double dt, double** Dx, double** Dy, double** Dz);
void time_implicit_solver(PointStructure* myPointStruct, FieldVariables* myFieldVariables_n, double dt, int iter_momentum, int iter_pressure, int iter_timple, double omega, bool neumann_flag_boundary);

// Temporary
// void calculate_intermediate_velocity(PointStructure* myPointStruct, double* u_intermediate, double* v_intermediate, double* w_intermediate, double* u_n, double* v_n, double* w_n, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian);
// void calculate_mass_residual(PointStructure* myPointStruct, double* u_intermediate, double* v_intermediate, double* w_intermediate, double** Dx, double** Dy, double** Dz, double rho, double dt, double* mass_residual);
// void calculate_pressure(PointStructure* myPointStruct, double* p, double* u_intermediate, double* v_intermediate, double* w_intermediate, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian, double* mass_residual);
// void calculate_pressure_correction_implicit(PointStructure* myPointStruct, double* p, double* u_new, double* v_new, double* w_new, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian, double* mass_residual, int iter_pressure, double omega);
// void Poisson_solver(double* p, double* rhs, double** lhs, int** cloud_index, int num_nodes, int num_cloud_points, int num_boundary_points, int num_iter);

/////////////////////////////////////////////////////////////////////////////
// Function Definitions
////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
// Fractional Step Explicit Solver Modules
/////////////////////////////////////////////////////////////////////////////
void fractional_step_explicit(PointStructure* myPointStruct, FieldVariables* field, double dt, int iter_pressure, double omega, bool neumann_flag_boundary)
{   
    double* u_intermediate = create_vector(myPointStruct->num_nodes);
    double* v_intermediate = create_vector(myPointStruct->num_nodes);
    double* w_intermediate = create_vector(myPointStruct->num_nodes);
    double* mass_residual = create_vector(myPointStruct->num_nodes);
    double* dpdn = create_vector(myPointStruct->num_boundary_nodes);

    calculate_intermediate_velocity(myPointStruct, u_intermediate, v_intermediate, w_intermediate, field->u, field->v, field->w, field->rho, field->mu, dt, myPointStruct->Dx, myPointStruct->Dy, myPointStruct->Dz, myPointStruct->lap);
    calculate_mass_residual(myPointStruct, u_intermediate, v_intermediate, w_intermediate, field->rho, dt, mass_residual);
    dpdn_from_momentum_equation(dpdn, myPointStruct, field->u, field->v, field->w, u_intermediate, v_intermediate, w_intermediate, field->rho, dt);
    // dpdn_from_complete_momentum_equation(dpdn, myPointStruct, field, u_intermediate, v_intermediate, w_intermediate, dt);
    solve_poisson_equation(myPointStruct, field->p, mass_residual, dpdn, iter_pressure, omega, neumann_flag_boundary);
    update_velocity(myPointStruct, field->u, field->v, field->w, u_intermediate, v_intermediate, w_intermediate, field->p, field->rho, dt, myPointStruct->Dx, myPointStruct->Dy, myPointStruct->Dz);
    
    free(u_intermediate);
    free(v_intermediate);
    free(w_intermediate);
    free(mass_residual);
    free(dpdn);
}

void calculate_intermediate_velocity(PointStructure* myPointStruct, double* u_intermediate, double* v_intermediate, double* w_intermediate,double* u_n, double* v_n, double* w_n, double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, double** laplacian)
{
    // Calculate the intermediate velocity
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    int num_boundary_points = myPointStruct->num_boundary_nodes;
    double *temp1 = create_vector(num_nodes);
    double *temp2 = create_vector(num_nodes);
    double *temp3 = create_vector(num_nodes);
    double *temp4 = create_vector(num_nodes);
    double nu = mu/rho;

// x-momentum
    multiply_sparse_matrix_vector(Dx, u_n, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(Dy, u_n, temp2, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    if (parameters.dimension == 3)
        multiply_sparse_matrix_vector(Dz, u_n, temp3, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(laplacian, u_n, temp4, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    for (int i = 0; i < num_nodes; i++)
        u_intermediate[i] = u_n[i] - dt * (u_n[i] * temp1[i] + v_n[i] * temp2[i] + w_n[i] * temp3[i] - nu *temp4[i]);

// y-momentum
    multiply_sparse_matrix_vector(Dx, v_n, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(Dy, v_n, temp2, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    if (parameters.dimension == 3)
        multiply_sparse_matrix_vector(Dz, v_n, temp3, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(laplacian, v_n, temp4, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes,num_cloud_points);
    for (int i = 0; i < num_nodes; i++)
        v_intermediate[i] = v_n[i] - dt * (u_n[i] * temp1[i] + v_n[i] * temp2[i] + w_n[i] * temp3[i] - nu * temp4[i]);

//  z-momentum
    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector(Dx, w_n, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector(Dy, w_n, temp2, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector(Dz, w_n, temp3, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector(laplacian, w_n, temp4, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes,num_cloud_points);
        for (int i = 0; i < num_nodes; i++)
            w_intermediate[i] = w_n[i] - dt * (u_n[i] * temp1[i] + v_n[i] * temp2[i] + w_n[i] * temp3[i] -nu * temp4[i]);
    }
 
    free_vector(temp1);
    free_vector(temp2);
    free_vector(temp3);
    free_vector(temp4);
}

void calculate_mass_residual(PointStructure* myPointStruct, double* u_new, 
            double* v_new, double* w_new, double rho, double dt, double* mass_residual){
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    int num_boundary_points = myPointStruct->num_boundary_nodes;
    double* temp1 = create_vector(num_nodes);
    double* temp2 = create_vector(num_nodes);
    double* temp3 = create_vector(num_nodes);
    multiply_sparse_matrix_vector(myPointStruct->Dx, u_new, temp1, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector(myPointStruct->Dy, v_new, temp2, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector(myPointStruct->Dz, w_new, temp3, myPointStruct->cloud_index, myPointStruct->boundary_tag, num_nodes, num_cloud_points);
    }
    for (int i = num_boundary_points; i < num_nodes; i++){
        mass_residual[i] = rho*(temp1[i]+temp2[i]+temp3[i])/dt;
    }
    free(temp1);
    free(temp2);
    free(temp3);
}

void dpdn_from_momentum_equation(double* dpdn, PointStructure* myPointStruct, double* u_n, double* v_n, double* w_n, double* u_intermediate, double* v_intermediate, double* w_intermediate, double rho, double dt){
    double dpdx, dpdy, dpdz;
    for(int i = 0; i < myPointStruct->num_boundary_nodes; i++){
        dpdx = (u_intermediate[i] - u_n[i]) * rho/dt; 
        dpdy = (v_intermediate[i] - v_n[i]) * rho/dt;
        dpdn[i] = 0; 
        if (parameters.dimension == 3){
            dpdz = (w_intermediate[i] - w_n[i]) * rho/dt; 
            dpdn[i] += dpdz*myPointStruct->z_normal[i];
        }
        dpdn[i] += dpdx*myPointStruct->x_normal[i] + dpdy*myPointStruct->y_normal[i];
    }
}

void dpdn_from_complete_momentum_equation(double* dpdn, PointStructure* myPointStruct, FieldVariables* field, double* u_intermediate, double* v_intermediate, double* w_intermediate, double dt){
    // double rho = rho;
    double temp1, temp2, temp3, temp4;
    double nu = field->mu / field->rho;
    for(int i = 0; i < myPointStruct->num_boundary_nodes; i++){
        temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0;
        for(int j = 0; j<myPointStruct->num_cloud_points; j++){
            temp1 = nu * myPointStruct->lap[i][j] * field->u[myPointStruct->cloud_index[i][j]];
            temp1 -= field->u[i]* field->u[myPointStruct->cloud_index[i][j]] * myPointStruct->Dx[i][j];
            temp1 -= field->v[i]* field->u[myPointStruct->cloud_index[i][j]] * myPointStruct->Dy[i][j];
            if (parameters.dimension == 3){
                temp1 -= field->w[i]* field->u[myPointStruct->cloud_index[i][j]] * myPointStruct->Dz[i][j];
            }
            
            temp2 = nu* myPointStruct->lap[i][j] * field->v[myPointStruct->cloud_index[i][j]];
            temp2 -= field->u[i]* field->v[myPointStruct->cloud_index[i][j]] * myPointStruct->Dx[i][j];
            temp2 -= field->v[i]* field->v[myPointStruct->cloud_index[i][j]] * myPointStruct->Dy[i][j];
            if (parameters.dimension == 3){
                temp2 -= field->w[i]* field->v[myPointStruct->cloud_index[i][j]] * myPointStruct->Dz[i][j];
            }
            
            if (parameters.dimension == 3){
                temp3 = nu* myPointStruct->lap[i][j] * field->w[myPointStruct->cloud_index[i][j]];
                temp3 -= field->u[i]* field->w[myPointStruct->cloud_index[i][j]] * myPointStruct->Dx[i][j];
                temp3 -= field->v[i]* field->w[myPointStruct->cloud_index[i][j]] * myPointStruct->Dy[i][j];
                temp3 -= field->w[i]* field->w[myPointStruct->cloud_index[i][j]] * myPointStruct->Dz[i][j];
            }
        }
        dpdn[i] = field->rho * (temp1*myPointStruct->x_normal[i] + temp2*myPointStruct->y_normal[i] + temp3*myPointStruct->z_normal[i]);
        // printf("%e, %e, %e, %e\n", temp1, temp2, temp3, dpdn[i]);
    }
}

void solve_poisson_equation(PointStructure* myPointStruct, double* P, double* source, double* dpdn, int iter_pressure, double omega, bool neumann_flag_boundary){
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    int num_boundary_points = myPointStruct->num_boundary_nodes;

    double* ptemp = create_vector(num_nodes);

    // Solve the system of equations
    double sum, sumx, sumy, sumz, Ap, pref;// dpdn = 0;
    
    for (int iter = 0; iter < iter_pressure; iter++)
    {       
        for (int i = 0; i < num_nodes; i++){
            ptemp[i] = P[i]; // for L2 norm calculation
        }
    // Loop over interior nodes
        for (int i = num_boundary_points; i < num_nodes; i++){
            sum = 0;
            for (int j = 1; j < num_cloud_points; j++){
                sum += myPointStruct->lap[i][j]*P[myPointStruct->cloud_index[i][j]];
            }
            P[i] = omega*((source[i]-sum)/myPointStruct->lap[i][0]) - (1-omega)*P[i];
        }
    
        // Loop over boundary nodes, set the Neumann Boundary conditions: dpdn = 0
        if (neumann_flag_boundary){
            for (int i = myPointStruct->num_corners; i < num_boundary_points; i++){
                sumx = 0; sumy = 0; sumz = 0; Ap = 0;
                for (int j = 1; j < num_cloud_points; j++){
                    sumx += myPointStruct->Dx[i][j]*P[myPointStruct->cloud_index[i][j]];
                    sumy += myPointStruct->Dy[i][j]*P[myPointStruct->cloud_index[i][j]];
                    if (parameters.dimension == 3){
                        sumz += myPointStruct->Dz[i][j]*P[myPointStruct->cloud_index[i][j]];
                    }
                }
                Ap += myPointStruct->Dx[i][0]*myPointStruct->x_normal[i];
                Ap += myPointStruct->Dy[i][0]*myPointStruct->y_normal[i];
                if (parameters.dimension == 3){
                    Ap += myPointStruct->Dz[i][0]*myPointStruct->z_normal[i];
                }
                P[i] = (dpdn[i]-sumx*myPointStruct->x_normal[i] -sumy*myPointStruct->y_normal[i] -sumz*myPointStruct->z_normal[i])/Ap;
            }
            pref = P[0];
            for (int i = 0; i < num_nodes; i++){
                P[i] = P[i] - pref;
            }
        }
    }
    free(ptemp);
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
    
    free(dpdx);
    free(dpdy);
    free(dpdz);
}

////////////////////////////////////////////////////////////////////////////////////
// Time-Implicit Solver Modules
////////////////////////////////////////////////////////////////////////////////////
void time_implicit_solver(PointStructure* myPointStruct, FieldVariables* field, double dt, int iter_momentum, int iter_pressure, int iter_timple, double omega, bool neumann_flag_boundary){
    double* u_new = create_vector(myPointStruct->num_nodes);
    double* v_new = create_vector(myPointStruct->num_nodes);
    double* w_new = create_vector(myPointStruct->num_nodes);
    double* u_old = create_vector(myPointStruct->num_nodes);
    double* v_old = create_vector(myPointStruct->num_nodes);
    double* w_old = create_vector(myPointStruct->num_nodes);
    double* pprime = create_vector(myPointStruct->num_nodes);
    double* dpdn = create_vector(myPointStruct->num_boundary_nodes);
    double* mass_residual = create_vector(myPointStruct->num_nodes);

    for (int i = 0; i < myPointStruct->num_nodes; i++){
        u_old[i] = field->u[i]; u_new[i] = field->u[i];
        v_old[i] = field->v[i]; v_new[i] = field->v[i];
        w_old[i] = field->w[i]; w_new[i] = field->w[i];
    }
    
    for (int iter = 0; iter<iter_timple; iter++){
        calculate_intermediate_velocity_implicit(myPointStruct, u_new, v_new, w_new, 
                        u_old, v_old, w_old, field->u, field->v, field->w, field->p, 
                        field->rho, field->mu, dt, myPointStruct->Dx, myPointStruct->Dy, 
                        myPointStruct->Dz, myPointStruct->lap, iter_momentum);
        calculate_mass_residual(myPointStruct, u_new, v_new, w_new, field->rho, dt, mass_residual);
        // dpdn_from_momentum_equation(dpdn, myPointStruct, u_old, v_old, w_old, u_new, v_new, w_new, field->rho, dt);
        solve_poisson_equation(myPointStruct, pprime, mass_residual, dpdn, iter_pressure, omega, neumann_flag_boundary);
        update_velocity_implicit(myPointStruct, u_old, v_old, w_old, field->p, u_new, v_new, w_new, pprime, field->rho, dt, myPointStruct->Dx, myPointStruct->Dy, myPointStruct->Dz);
    }
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        field->u[i] = u_new[i];
        field->v[i] = v_new[i];
        field->w[i] = w_new[i];
    }

    free(u_new);
    free(v_new);
    free(w_new);
    free(mass_residual);
    free(u_old);
    free(v_old);
    free(w_old);
    free(pprime);
    free(dpdn);
}

void calculate_intermediate_velocity_implicit(PointStructure* myPointStruct, 
            double* u_new, double* v_new, double* w_new, double* u_old, double* v_old, 
            double* w_old, double *u_n, double *v_n, double *w_n,double* p_n, 
            double rho, double mu, double dt, double** Dx, double** Dy, double** Dz, 
            double** laplacian, int iter_momentum){
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
                    // temp4 += laplacian[i][j] * u_new[myPointStruct->cloud_index[i][j]];
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
                    // temp4 += laplacian[i][j] * v_new[myPointStruct->cloud_index[i][j]];
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
        u_new[i] = u_new[i] - dt * dpdx[i]/rho;
        v_new[i] = v_new[i] - dt * dpdy[i]/rho;
        if (parameters.dimension == 3){
            w_new[i] = w_new[i] - dt * dpdz[i]/rho;
        }
        p[i] = p[i] + pprime[i];
        u_old[i] = u_new[i];
        v_old[i] = v_new[i];
        if (parameters.dimension == 3){
            w_old[i] = w_new[i];
        }
    }
    for (int i = 0; i < num_boundary_points; i++){
        p[i] = p[i] + pprime[i];
    }
}


/////////////////////////////////////////////////////////////////////////////
// Test Heat conduction V-Cycle Solver Modules
/////////////////////////////////////////////////////////////////////////////

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

void relax_vcycle(PointStructure* mypointstruct, double* T, double* rhs, int num_iter, double omega, bool neumann_flag_boundary){
    int num_nodes = mypointstruct->num_nodes;
    int num_cloud_points = mypointstruct->num_cloud_points;
    int num_boundary_points = mypointstruct->num_boundary_nodes;
    double* Ttemp = create_vector(num_nodes);
    double sum = 0, Ap = 0, dpdn = 0;
    double k = 1.0;

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
        if (neumann_flag_boundary){// && (iter%5 == 0)){
            for (int i = 0; i < num_boundary_points; i++)
            {
                sum = 0; Ap = 0;
                for (int j = 1; j < num_cloud_points; j++){
                    sum += mypointstruct->Dx[i][j]*T[mypointstruct->cloud_index[i][j]]*mypointstruct->x_normal[i];
                    sum += mypointstruct->Dy[i][j]*T[mypointstruct->cloud_index[i][j]]*mypointstruct->y_normal[i];
                    if (parameters.dimension == 3){
                        sum += mypointstruct->Dz[i][j]*T[mypointstruct->cloud_index[i][j]]*mypointstruct->z_normal[i];
                    }
                }
                Ap += mypointstruct->Dx[i][0]*mypointstruct->x_normal[i];
                Ap += mypointstruct->Dy[i][0]*mypointstruct->y_normal[i];
                if (parameters.dimension == 3){
                    Ap += mypointstruct->Dz[i][0]*mypointstruct->z_normal[i];
                }
                dpdn = 2*3.14*k*(cos(2*3.14*k*mypointstruct->x[i])*sin(2*3.14*k*mypointstruct->y[i])*mypointstruct->x_normal[i]
                        +cos(2*3.14*k*mypointstruct->y[i])*sin(2*3.14*k*mypointstruct->x[i])*mypointstruct->y_normal[i]);
                T[i] = (dpdn-sum)/Ap;
            }
        }
    }
    free(Ttemp);
}



#endif
