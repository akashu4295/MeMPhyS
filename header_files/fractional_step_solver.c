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
#include "functions.h"


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
    # pragma acc parallel loop present(field[0], myPointStruct[0])
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

    # pragma acc parallel loop present(field[0], myPointStruct[0]) reduction(+:steady_state_error)
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        steady_state_error += pow((field[0].u[i] - field[0].u_old[i]), 2) + pow((field[0].v[i] - field[0].v_old[i]), 2) + pow((field[0].w[i] - field[0].w_old[i]), 2);
    }
    printf("Steady state error (before sqrt): %e\n", steady_state_error);
    steady_state_error = sqrt(steady_state_error/myPointStruct[0].num_nodes);
    return steady_state_error;
}


double fractional_step_explicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field)
{   
    double steady_state_error = 0.0;

    #pragma acc parallel loop present(field[0], myPointStruct[0])
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        field[0].u_old[i] = field[0].u[i];
        field[0].v_old[i] = field[0].v[i];
        // printf("u, v, u_old, v_old: %lf, %lf, %lf, %lf\n", field[0].u[i], field[0].v[i], field[0].u_old[i], field[0].v_old[i]);
    }

    #pragma acc data present(field[:parameters.num_levels], myPointStruct[:parameters.num_levels], parameters)
    {
        FS_calculate_intermediate_velocity_vectorised_2d(&myPointStruct[0], &field[0]);
        FS_calculate_mass_residual_vectorised_2d(&myPointStruct[0], &field[0]);
        FS_calculate_boundary_dpdn_vectorised_2d(&myPointStruct[0], &field[0]);
        FS_multigrid_Poisson_solver_vectorised(myPointStruct, field);
        FS_update_velocity_vectorised_2d(&myPointStruct[0], &field[0]);
    }

    #pragma acc parallel loop present(field[0], myPointStruct[0]) reduction(+:steady_state_error)
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        double du = field[0].u[i] - field[0].u_old[i];
        double dv = field[0].v[i] - field[0].v_old[i];
        steady_state_error += du*du + dv*dv;
    }
    // printf("Steady state error (before sqrt): %e\n", steady_state_error);
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

}

// # pragma acc routine
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

    /* ---- ENFORCE BOUNDARY CONDITIONS ON u* ---- */
    apply_velocity_bc_intermediate(myPointStruct, field);
}

void apply_velocity_bc_intermediate(PointStructure* myPointStruct, FieldVariables* field)
{
    int n = myPointStruct->num_nodes;
    # pragma acc parallel loop gang vector present(field, myPointStruct)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (myPointStruct->boundary_tag[i] && !myPointStruct->corner_tag[i]){
            if (myPointStruct->node_bc[i].type == BC_VELOCITY_INLET 
                    || myPointStruct->node_bc[i].type == BC_WALL 
                    || myPointStruct->node_bc[i].type == BC_VELOCITY_OUTLET){
                field->u_new[i] = myPointStruct->node_bc[i].u;
                field->v_new[i] = myPointStruct->node_bc[i].v;
            }
        }
    }
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
        // if (!myPointStruct->boundary_tag[i])
        field->source[i] = parameters.rho*(field->dpdx[i]+field->dpdy[i]+field->dpdz[i])/parameters.dt;
}

// # pragma acc routine
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
    for (int i = 0; i < num_nodes; i++){
        if (myPointStruct->corner_tag[i]) // skip corners
            continue;
        else if (!myPointStruct->boundary_tag[i]) // interior nodes
            field->source[i] = parameters.rho*(field->dpdx[i]+field->dpdy[i])/parameters.dt;
        else if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET) // pressure outlet boundary nodes
            field->source[i] = myPointStruct->node_bc[i].p;
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
        if (myPointStruct->boundary_tag[i] && !myPointStruct->corner_tag[i]){
            dpdx = (field->u_new[i] - field->u[i]) * parameters.rho/parameters.dt; 
            dpdy = (field->v_new[i] - field->v[i]) * parameters.rho/parameters.dt;
            dpdz = (field->w_new[i] - field->w[i]) * parameters.rho/parameters.dt; 
            field->dpdn[i] = dpdx*myPointStruct->x_normal[i] + dpdy*myPointStruct->y_normal[i] + dpdz*myPointStruct->z_normal[i];
        }
    }
}

// # pragma acc routine
void FS_calculate_boundary_dpdn_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field)
{
    double dpdx, dpdy;
    // Copy only boundary normals to gpu
    # pragma acc parallel loop gang vector present(field, myPointStruct, parameters)
    for(int i = 0; i < myPointStruct->num_nodes; i++){
        if (myPointStruct->corner_tag[i]) // skip corners
            continue;
        if ((myPointStruct->boundary_tag[i]) && (myPointStruct->node_bc[i].type != BC_PRESSURE_OUTLET)){
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

void FS_relaxation_vectorised(PointStructure* mypointstruct, FieldVariables* field)
{
    int n = mypointstruct->num_cloud_points;
    int N = mypointstruct->num_nodes;

    #pragma acc data present(field, mypointstruct, parameters)
    for (int iter = 0; iter < parameters.num_relax; iter++) {
        // Relaxation step
        #pragma acc parallel loop gang vector_length(128) present(field, mypointstruct, parameters)
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            #pragma acc loop vector reduction(+:sum)
            for (int j = 1; j < n; j++) {
                sum += mypointstruct->lap_Poison[i*n + j] *
                       field->p[mypointstruct->cloud_index[i*n + j]];
            }
            field->p[i] = parameters.omega *
                ((field->source[i] - sum) / mypointstruct->lap_Poison[i*n])
                + (1 - parameters.omega) * field->p[i];
        }

        if (!mypointstruct->flag_outlets) {
            double pref = field->p[0];
            // Normalization + residual in one pass
            #pragma acc parallel loop gang vector_length(128) present(field, mypointstruct)
            for (int i = 0; i < N; i++) {
                double val = field->p[i] - pref;
                field->p[i] = val;
                double res = field->source[i];
                #pragma acc loop vector reduction(+:res)
                for (int j = 0; j < n; j++) {
                    res -= mypointstruct->lap[i*n + j] *
                        field->p[mypointstruct->cloud_index[i*n + j]];
                }
                field->res[i] = res;
            }
        }
        else {
            // Just residual calculation
            #pragma acc parallel loop gang vector_length(128) present(field, mypointstruct)
            for (int i = 0; i < N; i++) {
                double res = field->source[i];
                #pragma acc loop vector reduction(+:res)
                for (int j = 0; j < n; j++) {
                    res -= mypointstruct->lap[i*n + j] *
                        field->p[mypointstruct->cloud_index[i*n + j]];
                }
                field->res[i] = res;
            }
        }
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

void FS_update_velocity_vectorised(PointStructure* myPointStruct, FieldVariables* field)
{
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->p, field->dpdx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->p, field->dpdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->p, field->dpdz, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    
    // Update Interior nodes
    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (!myPointStruct->boundary_tag[i] || myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET){
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
        if (!myPointStruct->boundary_tag[i] || myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET){
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



// Obsolete functions

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
                int k = i*n + j;
                sumx += mypointstruct->Dx[k]*field->p[mypointstruct->cloud_index[k]];
                sumy += mypointstruct->Dy[k]*field->p[mypointstruct->cloud_index[k]];
                sumz += mypointstruct->Dz[k]*field->p[mypointstruct->cloud_index[k]];
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

int BiCGStab_Solve(PointStructure* ps,
                   const double* b,   // RHS = field->source
                   double* x,         // solution = field->p
                   int max_iter,
                   double tol)
{
    int N = ps->num_nodes;
    int n = ps->num_cloud_points;

    // Allocate temporary vectors
    double *r = malloc(N * sizeof(double));
    double *r0 = malloc(N * sizeof(double));
    double *p = malloc(N * sizeof(double));
    double *v = malloc(N * sizeof(double));
    double *t = malloc(N * sizeof(double));
    double *Ax = malloc(N * sizeof(double));

    // Compute initial residual: r = b - A*x    
    for (int i = 0; i < N; i++) {
        double Ax = 0.0;
        for (int j = 0; j < n; j++) {
            int col = ps->cloud_index[i*n + j];
            Ax += ps->lap_Poison[i*n + j] * x[col];
        }
        r[i] = b[i] - Ax;
        r0[i] = r[i];
    }

    for (int i = 0; i < N; i++) {
        r[i] = b[i] - Ax[i];
        r0[i] = r[i];
        p[i] = 0.0;
        v[i] = 0.0;
    }

    double rho = 1, alpha = 1, omega = 1;
    double rho_new;

    for (int iter = 0; iter < max_iter; iter++) {

        // rho = r0^T r
        rho_new = 0.0;
        for (int i = 0; i < N; i++)
            rho_new += r0[i] * r[i];

        if (fabs(rho_new) < 1e-30) break;

        // p = r + (rho_new/rho)*(alpha/omega)*(p - omega*v)
        double beta = (rho_new/rho)*(alpha/omega);
        for (int i = 0; i < N; i++)
            p[i] = r[i] + beta*(p[i] - omega*v[i]);

        // v = A*p
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                int col = ps->cloud_index[i*n + j];
                sum += ps->lap_Poison[i*n + j] * p[col];
            }
            v[i] = sum;
        }   

        // alpha = rho_new / (r0^T v)
        double r0v = 0.0;
        for (int i = 0; i < N; i++)
            r0v += r0[i] * v[i];

        alpha = rho_new / r0v;

        // s = r - alpha*v
        double *s = t; // reuse t buffer for s
        for (int i = 0; i < N; i++)
            s[i] = r[i] - alpha*v[i];

        // Check if |s| small enough → converged early
        double norm_s = 0.0;
        for (int i = 0; i < N; i++)
            norm_s += s[i]*s[i];

        if (sqrt(norm_s) < tol) {
            for (int i = 0; i < N; i++)
                x[i] += alpha*p[i];
            break;
        }

        // t = A*s
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                int col = ps->cloud_index[i*n + j];
                sum += ps->lap_Poison[i*n + j] * s[col];
            }
            t[i] = sum;
        }

        // omega = (t^T s) / (t^T t)
        double ts = 0.0, tt = 0.0;
        for (int i = 0; i < N; i++) {
            ts += t[i] * s[i];
            tt += t[i] * t[i];
        }
        omega = ts / tt;

        // x = x + alpha*p + omega*s
        for (int i = 0; i < N; i++)
            x[i] += alpha*p[i] + omega*s[i];

        // r = s - omega*t
        for (int i = 0; i < N; i++)
            r[i] = s[i] - omega*t[i];

        // Check convergence
        double norm_r = 0.0;
        for (int i = 0; i < N; i++)
            norm_r += r[i]*r[i];

        if (sqrt(norm_r) < tol) break;

        rho = rho_new;
    }

    free(r); free(r0); free(p);
    free(v); free(t); free(Ax);

    return 0; // success
}