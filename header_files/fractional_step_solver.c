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
        // if (parameters.poisson_solver_type ==1)
            field[0].p_old[i] = field[0].p[i];
    }

    FS_calculate_intermediate_velocity_vectorised(&myPointStruct[0], &field[0]);
    FS_calculate_mass_residual_vectorised(&myPointStruct[0], &field[0]);
    FS_multigrid_Poisson_solver_vectorised(myPointStruct, field);
    FS_update_velocity_vectorised(&myPointStruct[0], &field[0]);

    # pragma acc parallel loop present(field[0], myPointStruct[0]) reduction(+:steady_state_error)
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        double du = field[0].u[i] - field[0].u_old[i];
        double dv = field[0].v[i] - field[0].v_old[i];
        double dw = field[0].w[i] - field[0].w_old[i];
        steady_state_error += du*du + dv*dv + dw*dw;
    }
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
        // if (parameters.poisson_solver_type == 1)
            field[0].p_old[i] = field[0].p[i];
    }

    #pragma acc data present(field[:parameters.num_levels], myPointStruct[:parameters.num_levels], parameters)
    {
        FS_calculate_intermediate_velocity_vectorised_2d(&myPointStruct[0], &field[0]);
        FS_calculate_mass_residual_vectorised_2d(&myPointStruct[0], &field[0]);
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

    /* ---- ENFORCE BOUNDARY CONDITIONS ON u* ---- */
    # pragma acc parallel loop gang vector present(field, myPointStruct)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (myPointStruct->boundary_tag[i] && !myPointStruct->corner_tag[i]){
            if (myPointStruct->node_bc[i].type == BC_VELOCITY_INLET 
                    || myPointStruct->node_bc[i].type == BC_WALL 
                    || myPointStruct->node_bc[i].type == BC_VELOCITY_OUTLET){
                field->u_new[i] = myPointStruct->node_bc[i].u;
                field->v_new[i] = myPointStruct->node_bc[i].v;
                field->w_new[i] = myPointStruct->node_bc[i].w;
            }
        }
    }
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
    for (int i = 0; i < num_nodes; i++){
        if (myPointStruct->corner_tag[i]) // skip corners
            continue;
        if (!myPointStruct->boundary_tag[i])
            field->source[i] = parameters.rho*(field->dpdx[i]+field->dpdy[i]+field->dpdz[i])/parameters.dt;
        else if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET) // pressure outlet boundary nodes
            field->source[i] = myPointStruct->node_bc[i].p;
        else
            field->source[i] = parameters.rho*((field->u_new[i] - field->u[i])*myPointStruct->x_normal[i] + (field->v_new[i] - field->v[i])*myPointStruct->y_normal[i] + (field->w_new[i] - field->w[i])*myPointStruct->z_normal[i])/parameters.dt;
    }
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
        if (!myPointStruct->boundary_tag[i])
            field->source[i] = parameters.rho*(field->dpdx[i]+field->dpdy[i])/parameters.dt;
        else if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET) // pressure outlet boundary nodes
            field->source[i] = myPointStruct->node_bc[i].p;
        else
            field->source[i] = parameters.rho*((field->u_new[i] - field->u[i])*myPointStruct->x_normal[i] + (field->v_new[i] - field->v[i])*myPointStruct->y_normal[i])/parameters.dt;
    }
}

void FS_multigrid_Poisson_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field)
{
    for (int icycle = 0; icycle < parameters.num_vcycles; icycle++){
        for (int ilev = 0; ilev < parameters.num_levels; ilev++){
            if (parameters.poisson_solver_type == 1)    
                FS_relaxation_vectorised_Jacobi(&myPointStruct[ilev], &field[ilev]);
            else if (parameters.poisson_solver_type == 2)
                FS_relaxation_vectorised_Gauss_Seidel(&myPointStruct[ilev], &field[ilev]);
            else if (parameters.poisson_solver_type == 3)
                FS_relaxation_vectorised_bicgstab(&myPointStruct[ilev], &field[ilev]);
            FS_calculate_residuals_vectorised(&myPointStruct[ilev], &field[ilev]);
            if (ilev != parameters.num_levels-1){
                FS_restrict_residuals_vectorised(&myPointStruct[ilev], &myPointStruct[ilev+1], &field[ilev], &field[ilev+1]);
            }
        }
        for (int ilev = parameters.num_levels-1; ilev > 0; ilev--){
            FS_prolongate_corrections_vectorised(&myPointStruct[ilev-1], &myPointStruct[ilev], &field[ilev-1], &field[ilev]);
            if (ilev != 1) {
                if (parameters.poisson_solver_type == 1)    
                    FS_relaxation_vectorised_Jacobi(&myPointStruct[ilev-1], &field[ilev-1]);
                else if (parameters.poisson_solver_type == 2)
                    FS_relaxation_vectorised_Gauss_Seidel(&myPointStruct[ilev-1], &field[ilev-1]);
                else if (parameters.poisson_solver_type == 3)
                    FS_relaxation_vectorised_bicgstab(&myPointStruct[ilev-1], &field[ilev-1]);   
            } 
        }
    } 
}

void FS_relaxation_vectorised_Jacobi(PointStructure* mypointstruct, FieldVariables* field)
{
    int n = mypointstruct->num_cloud_points;
    int N = mypointstruct->num_nodes;
    #pragma acc parallel loop present(field)
    for (int i = 0; i < N; i++) {
        field->p_old[i] = field->p[i];
    }
    for (int iter = 0; iter < parameters.num_relax; iter++) {
        #pragma acc parallel loop gang vector present(field, mypointstruct, parameters)
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 1; j < n; j++) {
                int neighbor_idx = mypointstruct->cloud_index[i*n + j];
                sum += mypointstruct->lap_Poison[i*n + j] * field->p_old[neighbor_idx];
            }
            field->p[i] = parameters.omega * (field->source[i] - sum) / mypointstruct->lap_Poison[i*n] 
                            + (1.0 - parameters.omega) * field->p_old[i];
        }
        #pragma acc parallel loop present(field)
        for (int i = 0; i < N; i++) {
            field->p_old[i] = field->p[i];
        }
    }
}

void FS_relaxation_vectorised_Gauss_Seidel(PointStructure* mypointstruct, FieldVariables* field)
{
    int n = mypointstruct->num_cloud_points;
    int N = mypointstruct->num_nodes;
    #pragma acc data present(field, mypointstruct, parameters)
    for (int iter = 0; iter < parameters.num_relax; iter++) {
        #pragma acc parallel loop gang vector_length(128) present(field, mypointstruct, parameters)
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            #pragma acc loop vector reduction(+:sum)
            for (int j = 1; j < n; j++) {
                sum += mypointstruct->lap_Poison[i*n + j] *
                       field->p[mypointstruct->cloud_index[i*n + j]];
            }
            field->p[i] = parameters.omega * ((field->source[i] - sum) / mypointstruct->lap_Poison[i*n])
                + (1 - parameters.omega) * field->p[i];
        }
    }
}

void FS_relaxation_vectorised_bicgstab(PointStructure* mypointstruct, FieldVariables* field)
{
    int n = mypointstruct->num_cloud_points;
    int N = mypointstruct->num_nodes;

    int temp = BiCGStab_Solve_Preconditioned(mypointstruct, field->source, field->p, parameters.num_relax, parameters.poisson_solver_tolerance);
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
    double sum_res = 0.0;
    for (int i = 0; i < mypointStruct->num_nodes; i++){
        if (!mypointStruct->boundary_tag[i]){
            double sum = 0;
            # pragma acc loop reduction(+:sum)
            for (int j = 0; j < n; j++){
                sum += mypointStruct->lap_Poison[i*n + j]*field->p[mypointStruct->cloud_index[i*n + j]];
            }
            field->res[i] = field->source[i] - sum;
            sum_res += fabs(field->res[i]);
        }
    }
    printf("Poisson residual: %e\n", sum_res);
}

void FS_update_velocity_vectorised(PointStructure* myPointStruct, FieldVariables* field)
{
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->p, field->dpdx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->p, field->dpdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->p, field->dpdz, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    
    // Update Interior nodes
    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET){
            // normal velocity derivatve set to 0
            double sumux = 0.0, sumuy = 0.0, sumvx = 0.0, sumvy = 0.0, sumwx = 0.0, sumwy = 0.0;
            for (int j = 1; j < myPointStruct->num_cloud_points; j++){
                int k = i*myPointStruct->num_cloud_points + j;
                sumux -= myPointStruct->Dx[k] * field->u[myPointStruct->cloud_index[k]];
                sumuy -= myPointStruct->Dy[k] * field->u[myPointStruct->cloud_index[k]];
                sumvx -= myPointStruct->Dx[k] * field->v[myPointStruct->cloud_index[k]];
                sumvy -= myPointStruct->Dy[k] * field->v[myPointStruct->cloud_index[k]];
                sumwx -= myPointStruct->Dx[k] * field->w[myPointStruct->cloud_index[k]];
                sumwy -= myPointStruct->Dy[k] * field->w[myPointStruct->cloud_index[k]];
                }
            double Ap = myPointStruct->Dx[i*myPointStruct->num_cloud_points] * myPointStruct->x_normal[i]
                        + myPointStruct->Dy[i*myPointStruct->num_cloud_points] * myPointStruct->y_normal[i]
                        + myPointStruct->Dz[i*myPointStruct->num_cloud_points] * myPointStruct->z_normal[i];
            field->u[i] = (sumux * myPointStruct->x_normal[i] + sumuy * myPointStruct->y_normal[i]) / Ap;
            field->v[i] = (sumvx * myPointStruct->x_normal[i] + sumvy * myPointStruct->y_normal[i]) / Ap;
            field->w[i] = (sumwx * myPointStruct->x_normal[i] + sumwy * myPointStruct->y_normal[i]) / Ap;
        }
        else if (!myPointStruct->boundary_tag[i]){
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
        if (!myPointStruct->corner_tag[i])
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
        if (myPointStruct->node_bc[i].type == BC_PRESSURE_OUTLET){
            // normal velocity derivatve set to 0
            double sumux = 0.0, sumuy = 0.0, sumvx = 0.0, sumvy = 0.0;
            for (int j = 1; j < myPointStruct->num_cloud_points; j++){
                int k = i*myPointStruct->num_cloud_points + j;
                sumux -= myPointStruct->Dx[k] * field->u[myPointStruct->cloud_index[k]];
                sumuy -= myPointStruct->Dy[k] * field->u[myPointStruct->cloud_index[k]];
                sumvx -= myPointStruct->Dx[k] * field->v[myPointStruct->cloud_index[k]];
                sumvy -= myPointStruct->Dy[k] * field->v[myPointStruct->cloud_index[k]];
                }
            double Ap = myPointStruct->Dx[i*myPointStruct->num_cloud_points] * myPointStruct->x_normal[i]
                        + myPointStruct->Dy[i*myPointStruct->num_cloud_points] * myPointStruct->y_normal[i];
            field->u[i] = (sumux * myPointStruct->x_normal[i] + sumuy * myPointStruct->y_normal[i]) / Ap;
            field->v[i] = (sumvx * myPointStruct->x_normal[i] + sumvy * myPointStruct->y_normal[i]) / Ap;
        }
        else if (!myPointStruct->boundary_tag[i]){
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
        if (!myPointStruct->corner_tag[i])
            sum += parameters.rho*fabs(field->dpdx[i]+field->dpdy[i]);

    printf("Mass residual: %e\n", (sum)/myPointStruct->num_nodes);
}

int BiCGStab_Solve(PointStructure* ps, const double* b, double* x, int max_iter, double tol)
{
    int N = ps->num_nodes;
    int n = ps->num_cloud_points;

    // Allocate temporary vectors
    double *r = malloc(N * sizeof(double));
    double *r0 = malloc(N * sizeof(double));
    double *p = malloc(N * sizeof(double));
    double *v = malloc(N * sizeof(double));
    double *s = malloc(N * sizeof(double));
    double *t = malloc(N * sizeof(double));

    // Compute initial residual: r = b - A*x
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            int col = ps->cloud_index[i*n + j];
            sum += ps->lap_Poison[i*n + j] * x[col];
        }
        r[i] = b[i] - sum;
        r0[i] = r[i];
        p[i] = 0.0;
        v[i] = 0.0;
    }

    double rho = 1.0, alpha = 1.0, omega = 1.0;
    double rho_new;

    for (int iter = 0; iter < max_iter; iter++) {

        // rho_new = r0^T * r
        rho_new = 0.0;
        for (int i = 0; i < N; i++) {
            rho_new += r0[i] * r[i];
        }

        if (fabs(rho_new) < 1e-30) {
            printf("BiCGStab: rho too small, breaking at iter %d\n", iter);
            break;
        }

        // beta = (rho_new/rho) * (alpha/omega)
        double beta = (rho_new / rho) * (alpha / omega);
        
        // p = r + beta*(p - omega*v)
        for (int i = 0; i < N; i++) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        // v = A*p
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                int col = ps->cloud_index[i*n + j];
                sum += ps->lap_Poison[i*n + j] * p[col];
            }
            v[i] = sum;
        }

        // alpha = rho_new / (r0^T * v)
        double r0v = 0.0;
        for (int i = 0; i < N; i++) {
            r0v += r0[i] * v[i];
        }

        if (fabs(r0v) < 1e-30) {
            printf("BiCGStab: r0v too small, breaking at iter %d\n", iter);
            break;
        }

        alpha = rho_new / r0v;

        // s = r - alpha*v
        for (int i = 0; i < N; i++) {
            s[i] = r[i] - alpha * v[i];
        }

        // Check if |s| small enough → converged early
        double norm_s = 0.0;
        for (int i = 0; i < N; i++) {
            norm_s += s[i] * s[i];
        }

        if (sqrt(norm_s) < tol) {
            for (int i = 0; i < N; i++) {
                x[i] += alpha * p[i];
            }
            printf("BiCGStab converged early at iter %d, norm_s = %e\n", iter, sqrt(norm_s));
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

        // omega = (t^T * s) / (t^T * t)
        double ts = 0.0, tt = 0.0;
        for (int i = 0; i < N; i++) {
            ts += t[i] * s[i];
            tt += t[i] * t[i];
        }

        if (fabs(tt) < 1e-30) {
            printf("BiCGStab: tt too small, breaking at iter %d\n", iter);
            break;
        }

        omega = ts / tt;

        // x = x + alpha*p + omega*s
        for (int i = 0; i < N; i++) {
            x[i] += alpha * p[i] + omega * s[i];
        }

        // r = s - omega*t
        for (int i = 0; i < N; i++) {
            r[i] = s[i] - omega * t[i];
        }

        // Check convergence
        double norm_r = 0.0;
        for (int i = 0; i < N; i++) {
            norm_r += r[i] * r[i];
        }

        if (sqrt(norm_r) < tol) {
            printf("BiCGStab converged at iter %d, residual = %e\n", iter, sqrt(norm_r));
            break;
        }

        rho = rho_new;
        
        // Optional: print progress
        if (iter % 100 == 0) {
            printf("BiCGStab iter %d, residual = %e\n", iter, sqrt(norm_r));
        }
    }

    free(r);
    free(r0);
    free(p);
    free(v);
    free(s);
    free(t);

    return 0; // success
}

int BiCGStab_Solve_Preconditioned(PointStructure* ps, const double* b, double* x, int max_iter, double tol)
{
    int N = ps->num_nodes;
    int n = ps->num_cloud_points;

    // Extract diagonal for preconditioning
    double *diag_inv = malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) {
        double diag = ps->lap_Poison[i*n + 0];  // Assuming diagonal is at j=0
        if (fabs(diag) < 1e-14) diag = 1.0;     // Avoid division by zero
        diag_inv[i] = 1.0 / diag;
    }

    // Allocate temporary vectors
    double *r = malloc(N * sizeof(double));
    double *r0 = malloc(N * sizeof(double));
    double *p = malloc(N * sizeof(double));
    double *v = malloc(N * sizeof(double));
    double *s = malloc(N * sizeof(double));
    double *t = malloc(N * sizeof(double));
    double *z = malloc(N * sizeof(double));  // For preconditioning
    double *y = malloc(N * sizeof(double));  // For preconditioning

    // Compute initial residual: r = b - A*x
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            int col = ps->cloud_index[i*n + j];
            sum += ps->lap_Poison[i*n + j] * x[col];
        }
        r[i] = b[i] - sum;
        r0[i] = r[i];
        p[i] = 0.0;
        v[i] = 0.0;
    }

    // Compute initial residual norm
    double norm_r0 = 0.0;
    for (int i = 0; i < N; i++) {
        norm_r0 += r[i] * r[i];
    }
    norm_r0 = sqrt(norm_r0);
    // printf("BiCGStab initial residual = %e\n", norm_r0);

    double rho = 1.0, alpha = 1.0, omega = 1.0;
    double rho_new;

    for (int iter = 0; iter < max_iter; iter++) {

        // rho_new = r0^T * r
        rho_new = 0.0;
        for (int i = 0; i < N; i++) {
            rho_new += r0[i] * r[i];
        }

        if (fabs(rho_new) < 1e-30) {
            printf("BiCGStab: rho too small, breaking at iter %d\n", iter);
            break;
        }

        // beta = (rho_new/rho) * (alpha/omega)
        double beta = (rho_new / rho) * (alpha / omega);
        
        // p = r + beta*(p - omega*v)
        for (int i = 0; i < N; i++) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        // Precondition: z = M^-1 * p (Jacobi: z = diag^-1 * p)
        for (int i = 0; i < N; i++) {
            z[i] = diag_inv[i] * p[i];
        }

        // v = A*z
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                int col = ps->cloud_index[i*n + j];
                sum += ps->lap_Poison[i*n + j] * z[col];
            }
            v[i] = sum;
        }

        // alpha = rho_new / (r0^T * v)
        double r0v = 0.0;
        for (int i = 0; i < N; i++) {
            r0v += r0[i] * v[i];
        }

        if (fabs(r0v) < 1e-30) {
            printf("BiCGStab: r0v too small, breaking at iter %d\n", iter);
            break;
        }

        alpha = rho_new / r0v;

        // s = r - alpha*v
        for (int i = 0; i < N; i++) {
            s[i] = r[i] - alpha * v[i];
        }

        // Precondition: y = M^-1 * s
        for (int i = 0; i < N; i++) {
            y[i] = diag_inv[i] * s[i];
        }

        // t = A*y
        for (int i = 0; i < N; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                int col = ps->cloud_index[i*n + j];
                sum += ps->lap_Poison[i*n + j] * y[col];
            }
            t[i] = sum;
        }

        // omega = (t^T * s) / (t^T * t)
        double ts = 0.0, tt = 0.0;
        for (int i = 0; i < N; i++) {
            ts += t[i] * s[i];
            tt += t[i] * t[i];
        }

        if (fabs(tt) < 1e-30) {
            printf("BiCGStab: tt too small, breaking at iter %d\n", iter);
            break;
        }

        omega = ts / tt;

        // x = x + alpha*z + omega*y
        for (int i = 0; i < N; i++) {
            x[i] += alpha * z[i] + omega * y[i];
        }

        // r = s - omega*t
        for (int i = 0; i < N; i++) {
            r[i] = s[i] - omega * t[i];
        }

        // Check convergence
        double norm_r = 0.0;
        for (int i = 0; i < N; i++) {
            norm_r += r[i] * r[i];
        }
        norm_r = sqrt(norm_r);

        if (norm_r / norm_r0 < tol) {
            printf("BiCGStab converged at iter %d, relative residual = %e\n", 
                   iter, norm_r/norm_r0);
            break;
        }

        rho = rho_new;
        
        // Print progress
        // if (iter % (max_iter-1) == 0) {
        //     printf("BiCGStab iter %d, relative residual = %e\n", iter, norm_r/norm_r0);
        // }
    }

    free(r); free(r0); free(p); free(v); free(s); free(t);
    free(z); free(y); free(diag_inv);

    return 0;
}