// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#include <time.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "structures.h"
#include "mat_lib.h"
#include "solvers.h"

/////////////////////////////////////////////////////////////////////////////
// Fractional Step Explicit Solver Modules
/////////////////////////////////////////////////////////////////////////////
double calculate_dt(PointStructure* myPointStruct){
    double dt = 1e10;
    double d = myPointStruct->d_avg;
    double termd =    2*parameters.nu/(d*d);
    double termc =    1/d;
    dt = 0.5/(termd + termc);
    return (dt*parameters.courant_number);
}

double fractional_step_explicit_vectorised(PointStructure* myPointStruct, FieldVariables* field){   
    double steady_state_error = 0.0;
    # pragma acc parallel loop present(field, myPointStruct)
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        field[0].u_old[i] = field[0].u[i];
        field[0].v_old[i] = field[0].v[i];
        field[0].w_old[i] = field[0].w[i];
    }

    FS_calculate_intermediate_velocity_vectorised(&myPointStruct[0], &field[0]);
    FS_calculate_mass_residual_vectorised(&myPointStruct[0], &field[0]);
    FS_calculate_boundary_dpdn_vectorised(&myPointStruct[0], &field[0]);
    FS_multigrid_Poisson_solver_vectorised(myPointStruct, field);
    FS_update_velocity_vectorised(&myPointStruct[0], &field[0]);

    # pragma acc parallel loop present(field, myPointStruct) reduction(+:steady_state_error)
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        steady_state_error += pow((field[0].u[i] - field[0].u_old[i]), 2) + pow((field[0].v[i] - field[0].v_old[i]), 2) + pow((field[0].w[i] - field[0].w_old[i]), 2);
    }

    steady_state_error = sqrt(steady_state_error/myPointStruct[0].num_nodes);
    return steady_state_error;
}

double fractional_step_explicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field)
{   
    double steady_state_error = 0.0;

    #pragma acc parallel loop present(field, myPointStruct)
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        field[0].u_old[i] = field[0].u[i];
        field[0].v_old[i] = field[0].v[i];
    }

    #pragma acc data present(field, myPointStruct, parameters)
    {
        FS_calculate_intermediate_velocity_vectorised_2d(&myPointStruct[0], &field[0]);
        FS_calculate_mass_residual_vectorised_2d(&myPointStruct[0], &field[0]);
        FS_calculate_boundary_dpdn_vectorised_2d(&myPointStruct[0], &field[0]);
        FS_multigrid_Poisson_solver_vectorised_2d(myPointStruct, field);
        FS_update_velocity_vectorised_2d(&myPointStruct[0], &field[0]);
    }

    #pragma acc parallel loop present(field, myPointStruct) reduction(+:steady_state_error)
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        double du = field[0].u[i] - field[0].u_old[i];
        double dv = field[0].v[i] - field[0].v_old[i];
        steady_state_error += du*du + dv*dv;
        // printf("Node: %d, u: %lf, u_old: %lf, v: %lf, v_old: %lf, du: %lf, dv: %lf\n", i, field[0].u[i], field[0].u_old[i], field[0].v[i], field[0].v_old[i], du, dv);
    }

    steady_state_error = sqrt(steady_state_error / myPointStruct[0].num_nodes);
    return steady_state_error;
}


void FS_calculate_intermediate_velocity_vectorised(PointStructure* myPointStruct, FieldVariables* field)
{
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
// x-momentum
    #pragma acc data present(field, myPointStruct, parameters)
    {
        multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->u, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->u, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->u, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector_vectorised(myPointStruct->lap, field->u, field->dpdn, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    }    
    # pragma acc parallel loop gang vector present(field, myPointStruct, parameters)
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        field->u_new[i] = field->u[i] - parameters.dt * (field->u[i] * field->dpdx[i] + field->v[i] * field->dpdy[i] + field->w[i] * field->dpdz[i] - parameters.nu *field->dpdn[i]);

// y-momentum
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->v, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->v, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->v, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->lap, field->v, field->dpdn, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    
    # pragma acc parallel loop gang vector present(field, myPointStruct, parameters)
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        field->v_new[i] = field->v[i] - parameters.dt * (field->u[i] * field->dpdx[i] + field->v[i] * field->dpdy[i] + field->w[i] * field->dpdz[i] - parameters.nu * field->dpdn[i]);

//  z-momentum
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->w, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->w, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->w, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->lap, field->w, field->dpdn, myPointStruct->cloud_index, myPointStruct->num_nodes,myPointStruct->num_cloud_points);
    
    # pragma acc parallel loop gang vector present(field, myPointStruct, parameters)
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        field->w_new[i] = field->w[i] - parameters.dt * (field->u[i] * field->dpdx[i] + field->v[i] * field->dpdy[i] + field->w[i] * field->dpdz[i] - parameters.nu * field->dpdn[i]);

    // Remove this line unless you need to debug - it causes expensive GPU-to-CPU transfers
    // # pragma acc update host(field[0].u_new[:num_nodes], field[0].v_new[:num_nodes], field[0].w_new[:num_nodes])
}

void FS_calculate_intermediate_velocity_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field)
{
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
// x-momentum
    multiply_sparse_matrix_vector_vectorised_gpu(myPointStruct->Dx, field->u, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised_gpu(myPointStruct->Dy, field->u, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised_gpu(myPointStruct->lap, field->u, field->dpdn, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    
    # pragma acc parallel loop gang vector present(field, myPointStruct, parameters)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        field->u_new[i] = field->u[i] - parameters.dt * (field->u[i] * field->dpdx[i] + field->v[i] * field->dpdy[i] - parameters.nu *field->dpdn[i]);
    }
// y-momentum
    multiply_sparse_matrix_vector_vectorised_gpu(myPointStruct->Dx, field->v, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised_gpu(myPointStruct->Dy, field->v, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised_gpu(myPointStruct->lap, field->v, field->dpdn, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    
    # pragma acc parallel loop gang vector present(field, myPointStruct, parameters)
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        field->v_new[i] = field->v[i] - parameters.dt * (field->u[i] * field->dpdx[i] + field->v[i] * field->dpdy[i] - parameters.nu * field->dpdn[i]);

    // # pragma acc update host(field[0].u_new[0:num_nodes], field[0].v_new[0:num_nodes])
    // for (int i = 0; i < num_nodes; i++){
    //     printf("DEBUG inside FS u_new[%d]: %lf, v_new[%d]: %lf\n", i, field[0].u_new[i], i, field[0].v_new[i]);   
    // }
}

void FS_calculate_mass_residual_vectorised(PointStructure* myPointStruct, FieldVariables* field)
{
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
    multiply_sparse_matrix_vector_vectorised_gpu(myPointStruct->Dx, field->u_new, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised_gpu(myPointStruct->Dy, field->v_new, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised_gpu(myPointStruct->Dz, field->w_new, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    
    # pragma acc parallel loop gang vector present(field, myPointStruct, parameters)
    for (int i = 0; i < num_nodes; i++)
        if (!myPointStruct->boundary_tag[i])
            field->source[i] = parameters.rho*(field->dpdx[i]+field->dpdy[i]+field->dpdz[i])/parameters.dt;
}

void FS_calculate_mass_residual_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field)
{
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    # pragma acc data present(field, myPointStruct, parameters)
    {
        multiply_sparse_matrix_vector_vectorised_gpu(myPointStruct->Dx, field->u_new, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector_vectorised_gpu(myPointStruct->Dy, field->v_new, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    }

    # pragma acc parallel loop present(field, myPointStruct, parameters)
    for (int i = 0; i < num_nodes; i++)
        if (!myPointStruct->boundary_tag[i]){
            field->source[i] = parameters.rho*(field->dpdx[i]+field->dpdy[i])/parameters.dt;
        }
    // #pragma acc update host(field->source[0:num_nodes])
    // for (int i = 0; i < num_nodes; i++){
    //     printf("DEBUG source[%d]: %lf\n", i, field->source[i]);   
    // }
}

void FS_calculate_boundary_dpdn_vectorised(PointStructure* myPointStruct, FieldVariables* field)
{
    double dpdx, dpdy, dpdz;
    # pragma acc parallel loop gang vector present(field, myPointStruct, parameters)
    for(int i = 0; i < myPointStruct->num_nodes; i++){
        if (myPointStruct->boundary_tag[i]){
            dpdx = (field->u_new[i] - field->u[i]) * parameters.rho/parameters.dt; 
            dpdy = (field->v_new[i] - field->v[i]) * parameters.rho/parameters.dt;
            dpdz = (field->w_new[i] - field->w[i]) * parameters.rho/parameters.dt; 
            field->dpdn[i] = dpdx*myPointStruct->x_normal[i] + dpdy*myPointStruct->y_normal[i] + dpdz*myPointStruct->z_normal[i];
        }
    }
}

void FS_calculate_boundary_dpdn_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field)
{
    double dpdx, dpdy;
    // Copy only boundary normals to gpu
    # pragma acc parallel loop gang vector present(field, myPointStruct, parameters)
    for(int i = 0; i < myPointStruct->num_nodes; i++){
        if (myPointStruct->boundary_tag[i]){
            dpdx = (field->u_new[i] - field->u[i]) * parameters.rho/parameters.dt; 
            dpdy = (field->v_new[i] - field->v[i]) * parameters.rho/parameters.dt;
            field->dpdn[i] = dpdx*myPointStruct->x_normal[i] + dpdy*myPointStruct->y_normal[i];
        }
    }
}

void FS_multigrid_Poisson_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field)
{
    for (int icycle = 0; icycle < parameters.num_vcycles; icycle++){
        for (int ilev = 0; ilev < parameters.num_levels; ilev++){
            FS_relaxation_vectorised(&myPointStruct[ilev], &field[ilev]);
            FS_calculate_residuals_vectorised(&myPointStruct[ilev], &field[ilev]);
            if (ilev != parameters.num_levels-1){
                FS_restrict_residuals_vectorised(&myPointStruct[ilev], &myPointStruct[ilev+1], &field[ilev], &field[ilev+1]);
            }
        }
        for (int ilev = parameters.num_levels-1; ilev > 0; ilev--){
            FS_prolongate_corrections_vectorised(&myPointStruct[ilev-1], &myPointStruct[ilev], &field[ilev-1], &field[ilev]);
            if (ilev != 1) 
                FS_relaxation_vectorised(&myPointStruct[ilev-1], &field[ilev-1]);
        }   
        // update_boundary_pressure(&myPointStruct[0], &field[0]);
    } 
}

void FS_multigrid_Poisson_solver_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field)
{
    for (int icycle = 0; icycle < parameters.num_vcycles; icycle++){
        for (int ilev = 0; ilev < parameters.num_levels; ilev++){
            FS_relaxation_vectorised_2d(&myPointStruct[ilev], &field[ilev]);
            FS_calculate_residuals_vectorised_2d(&myPointStruct[ilev], &field[ilev]);
            if (ilev != parameters.num_levels-1){
                FS_restrict_residuals_vectorised_2d(&myPointStruct[ilev], &myPointStruct[ilev+1], &field[ilev], &field[ilev+1]);
                }
            }
        for (int ilev = parameters.num_levels-1; ilev > 0; ilev--){
            FS_prolongate_corrections_vectorised_2d(&myPointStruct[ilev-1], &myPointStruct[ilev], &field[ilev-1], &field[ilev]);
            if (ilev != 1) 
                FS_relaxation_vectorised_2d(&myPointStruct[ilev-1], &field[ilev-1]);
        }   
        // update_boundary_pressure(&myPointStruct[0], &field[0]);
    } 
}


void FS_relaxation_vectorised(PointStructure* mypointstruct, FieldVariables* field)
{
    int n = mypointstruct->num_cloud_points;
    for (int iter = 0; iter < parameters.num_relax; iter++){
        # pragma acc parallel loop gang vector present(field, mypointstruct, parameters)
        for (int i = 0; i < mypointstruct->num_nodes; i++){
            double sum = 0.0;
            # pragma acc loop reduction(+:sum)
            for (int j = 1; j < n; j++){
                sum += mypointstruct->lap_Poison[i*n + j]*field->p[mypointstruct->cloud_index[i*n + j]];
            }
            field->p[i] = parameters.omega*((field->source[i]-sum)/mypointstruct->lap_Poison[i*n]) + (1-parameters.omega)*field->p[i]; 
        }
        
        double pref = field->p[0];
        # pragma acc parallel loop gang vector present(field, mypointstruct)
        for (int i = 0; i < mypointstruct->num_nodes; i++){
            field->p[i] = field->p[i] - pref;
        }
        
        # pragma acc parallel loop gang vector present(field, mypointstruct, parameters)
        for (int i = 0; i< mypointstruct->num_nodes; i++){
            double res_val = field->source[i];
            # pragma acc loop reduction(-:res_val)
            for (int j = 0; j < n; j++){
                res_val -= mypointstruct->lap[i*n + j]*field->p[mypointstruct->cloud_index[i*n + j]];
            }
            field->res[i] = res_val;
        }
        // printf("Pressure Residual: %e\n", l2_norm_gen(mypointstruct, field->res, zeros, mypointstruct->num_nodes));
    }
}

void FS_relaxation_vectorised_2d(PointStructure* mypointstruct, FieldVariables* field)
{
    int n = mypointstruct->num_cloud_points;
    // if (!acc_is_present(mypointstruct->lap_Poison, mypointstruct->num_nodes * mypointstruct->num_cloud_points * sizeof(double)))
    //     printf("lap_Poison not present on GPU!\n");

    // if (!acc_is_present(mypointstruct->cloud_index, mypointstruct->num_nodes * mypointstruct->num_cloud_points * sizeof(int)))
    //     printf("cloud_index not present on GPU!\n");

    // if (!acc_is_present(field->p, mypointstruct->num_nodes * sizeof(double)))
    //     printf("field->p not present on GPU!\n");

    for (int iter = 0; iter < parameters.num_relax; iter++){
        # pragma acc parallel loop gang vector present(field, mypointstruct, parameters)
        for (int i = 0; i < mypointstruct->num_nodes; i++){
            double sum = 0.0;
            int k = i*n;
            # pragma acc loop reduction(+:sum)
            for (int j = 1; j < n; j++){
                k += 1;
                sum += mypointstruct->lap_Poison[k]*field->p[mypointstruct->cloud_index[k]];
            }
            field->p[i] = parameters.omega*((field->source[i]-sum)/mypointstruct->lap_Poison[i*n]) + (1-parameters.omega)*field->p[i]; 
        }
        double pref = field->p[0];
        # pragma acc parallel loop gang vector present(field, mypointstruct)
        for (int i = 0; i < mypointstruct->num_nodes; i++){
            field->p[i] = field->p[i] - pref;
        }
        
        # pragma acc parallel loop gang vector present(field, mypointstruct, parameters)
        for (int i = 0; i< mypointstruct->num_nodes; i++){
            field->res[i]=field->source[i];
            for (int j = 0; j < n; j++){
                int k = i*n + j;
                field->res[i] = field->res[i] - mypointstruct->lap[k]*field->p[mypointstruct->cloud_index[k]];
            }
        }
        // printf("Pressure Residual: %e\n", l2_norm_gen(mypointstruct, field->res, zeros, mypointstruct->num_nodes));
    
        // printf("Pressure Residual: %e\n", l2_norm_gen(mypointstruct, field->p, zeros, mypointstruct->num_nodes));
    }
}

void FS_restrict_residuals_vectorised(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, FieldVariables* field_f, FieldVariables* field_c)
{
    int n = mypointStruct_f->num_cloud_points;
    # pragma acc parallel loop gang vector present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_c->num_nodes; i++){
        if (mypointStruct_c->boundary_tag[i]==false){
            double results = 0.0;
            int i_restr = mypointStruct_c->restriction_points[i];
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < n; j++){
                results += mypointStruct_c->restr_mat[i*n + j] * field_f->res[mypointStruct_f->cloud_index[i_restr*n +j]];
            }
            field_c->source[i] = results;
        }
    }
}

void FS_restrict_residuals_vectorised_2d(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, FieldVariables* field_f, FieldVariables* field_c)
{
    int n = mypointStruct_f->num_cloud_points;
    # pragma acc parallel loop gang vector present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_c->num_nodes; i++){
        if (mypointStruct_c->boundary_tag[i]==false){
            double results = 0.0;
            int i_restr = mypointStruct_c->restriction_points[i];
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < n; j++){
                results += mypointStruct_c->restr_mat[i*n + j] * field_f->res[mypointStruct_f->cloud_index[i_restr*n +j]];
            }
            field_c->source[i] = results;
        }
    }
}

void FS_prolongate_corrections_vectorised(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, FieldVariables* field_f, FieldVariables* field_c)
{
    int n = mypointStruct_c->num_cloud_points;
    # pragma acc parallel loop gang vector present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_f->num_nodes; i++){
        if (mypointStruct_f->boundary_tag[i]==false){
            int i_prol = mypointStruct_f->prolongation_points[i];
            double results = 0.0;
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < n; j++){
                results += mypointStruct_f->prol_mat[i*n + j] * field_c->p[mypointStruct_c->cloud_index[i_prol*n +j]];
            }
            field_f->p[i] = field_f->p[i] + results;
        }
    }
    # pragma acc parallel loop gang vector present(field_c, mypointStruct_c)
    for (int i = 0; i<mypointStruct_c->num_nodes; i++){
        if (mypointStruct_c->boundary_tag[i]==false){
            field_c->p[i] = 0.0;
        }
    }
}
void FS_prolongate_corrections_vectorised_2d(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, FieldVariables* field_f, FieldVariables* field_c)
{
    int n = mypointStruct_c->num_cloud_points;
    # pragma acc parallel loop gang vector present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_f->num_nodes; i++){
        if (mypointStruct_f->boundary_tag[i]==false){
            int i_prol = mypointStruct_f->prolongation_points[i];
            double results = 0.0;
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < n; j++){
                results += mypointStruct_f->prol_mat[i*n + j] * field_c->p[mypointStruct_c->cloud_index[i_prol*n +j]];
            }
        field_f->p[i] = field_f->p[i] + results;
        }
    }
    # pragma acc parallel loop gang vector present(field_c, mypointStruct_c)
    for (int i = 0; i<mypointStruct_c->num_nodes; i++){
        if (mypointStruct_c->boundary_tag[i]==false){
            field_c->p[i] = 0.0;
        }
    }
}


void FS_calculate_residuals_vectorised(PointStructure* mypointStruct, FieldVariables* field)
{
    int n = mypointStruct->num_cloud_points;
    # pragma acc parallel loop gang vector present(field, mypointStruct)
    for (int i = 0; i < mypointStruct->num_nodes; i++){
        if (!mypointStruct->boundary_tag[i]){
            double sum = 0;
            # pragma acc loop reduction(+:sum)
            for (int j = 0; j < n; j++){
                sum += mypointStruct->lap[i*n + j]*field->p[mypointStruct->cloud_index[i*n + j]];
            }
            field->res[i] = field->source[i] - sum;
        }
    }
}
void FS_calculate_residuals_vectorised_2d(PointStructure* mypointStruct, FieldVariables* field){
    int n = mypointStruct->num_cloud_points;
    # pragma acc parallel loop gang vector present(field, mypointStruct)
    for (int i = 0; i < mypointStruct->num_nodes; i++){
        if (!mypointStruct->boundary_tag[i]){
            double sum = 0;
            # pragma acc loop reduction(+:sum)
            for (int j = 0; j < n; j++){
                sum += mypointStruct->lap[i*n+j]*field->p[mypointStruct->cloud_index[i*n+j]];
            }
            field->res[i] = field->source[i] - sum;
        }
    }
}

void FS_update_velocity_vectorised(PointStructure* myPointStruct, FieldVariables* field)
{
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->p, field->dpdx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->p, field->dpdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->p, field->dpdz, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    
    // Update Interior nodes
    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (!myPointStruct->boundary_tag[i]){
            field->u[i] = field->u_new[i] - parameters.dt * field->dpdx[i]/parameters.rho;
            field->v[i] = field->v_new[i] - parameters.dt * field->dpdy[i]/parameters.rho;
            field->w[i] = field->w_new[i] - parameters.dt * field->dpdz[i]/parameters.rho;
        }
    }
    
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->u, field->dpdx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->v, field->dpdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->w, field->dpdz, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    
    double sum = 0.0;
    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct) reduction(+:sum)
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        sum += parameters.rho*fabs(field->dpdx[i]+field->dpdy[i]+field->dpdz[i]);

    printf("Mass residual: %e\n", (sum)/myPointStruct->num_nodes);
}

void FS_update_velocity_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field)
{
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->p, field->dpdx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->p, field->dpdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    
    // Update Interior nodes
    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (!myPointStruct->boundary_tag[i]){
            field->u[i] = field->u_new[i] - parameters.dt * field->dpdx[i]/parameters.rho;
            field->v[i] = field->v_new[i] - parameters.dt * field->dpdy[i]/parameters.rho;
        }
    }
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->u, field->dpdx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->v, field->dpdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    
    double sum;
    sum = 0.0;
    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct) reduction(+:sum)
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        sum += parameters.rho*fabs(field->dpdx[i]+field->dpdy[i]);

    printf("Mass residual: %e\n", (sum)/myPointStruct->num_nodes);
}


void FS_update_boundary_pressure_vectorised(PointStructure* mypointstruct, FieldVariables* field)
{
    double sumx, sumy, sumz, Ap, term;
    int n = mypointstruct->num_cloud_points;
    if (parameters.neumann_flag_boundary){
        # pragma acc parallel loop gang vector present(field, parameters, mypointstruct)
        for (int i = 0; i < mypointstruct->num_boundary_nodes; i++)
        {
            sumx = 0; sumy = 0; sumz = 0; Ap = 0;
            # pragma acc loop reduction(+:sumx, sumy, sumz, Ap)
            for (int j = 1; j < n; j++){
                sumx += mypointstruct->Dx[i*n + j]*field->p[mypointstruct->cloud_index[i*n + j]];
                sumy += mypointstruct->Dy[i*n + j]*field->p[mypointstruct->cloud_index[i*n + j]];
                sumz += mypointstruct->Dz[i*n + j]*field->p[mypointstruct->cloud_index[i*n + j]];
            }
            Ap += mypointstruct->Dx[i*n]*mypointstruct->x_normal[i];
            Ap += mypointstruct->Dy[i*n]*mypointstruct->y_normal[i];
            Ap += mypointstruct->Dz[i*n]*mypointstruct->z_normal[i];

            term = (field->dpdn[i]-sumx*mypointstruct->x_normal[i] -sumy*mypointstruct->y_normal[i] -sumz*mypointstruct->z_normal[i])/Ap;
            field->p[i] = field->p[i] * 0.5 + 0.5 * term;
        }
    }
}

void FS_update_boundary_pressure_vectorised_2d(PointStructure* mypointstruct, FieldVariables* field){
    double sumx, sumy, Ap, term;
    int n = mypointstruct->num_cloud_points;
    if (parameters.neumann_flag_boundary){
        # pragma acc parallel loop gang vector present(field, parameters, mypointstruct)
        for (int i = 0; i < mypointstruct->num_boundary_nodes; i++)
        {
            sumx = 0; sumy = 0; Ap = 0;
            # pragma acc loop reduction(+:sumx, sumy, Ap)
            for (int j = 1; j < n; j++){
                int k = i*n + j;
                sumx += mypointstruct->Dx[k]*field->p[mypointstruct->cloud_index[k]];
                sumy += mypointstruct->Dy[k]*field->p[mypointstruct->cloud_index[k]];
            }
            Ap += mypointstruct->Dx[i*n]*mypointstruct->x_normal[i];
            Ap += mypointstruct->Dy[i*n]*mypointstruct->y_normal[i];
            term = (field->dpdn[i]-sumx*mypointstruct->x_normal[i] -sumy*mypointstruct->y_normal[i])/Ap;
            field->p[i] = field->p[i] * 0.5 + 0.5 * term;
        }
    }
}