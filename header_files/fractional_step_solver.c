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
        if (parameters.dimension == 3)
            field[0].w_old[i] = field[0].w[i];
    }

    FS_calculate_intermediate_velocity_vectorised(&myPointStruct[0], &field[0]);
    FS_calculate_mass_residual_vectorised(&myPointStruct[0], &field[0]);
    FS_calculate_boundary_dpdn_vectorised(&myPointStruct[0], &field[0]);
    FS_multigrid_Poisson_solver_vectorised(myPointStruct, field);
    FS_update_velocity_vectorised(&myPointStruct[0], &field[0]);

    # pragma acc parallel loop present(field, myPointStruct) reduction(+:steady_state_error)
    for (int i = 0; i < myPointStruct[0].num_nodes; i++){
        steady_state_error += (field[0].u[i] - field[0].u_old[i])*(field[0].u[i] - field[0].u_old[i]) + (field[0].v[i] - field[0].v_old[i])*(field[0].v[i] - field[0].v_old[i]);
        if (parameters.dimension == 3)
            steady_state_error += (field[0].w[i] - field[0].w_old[i])*(field[0].w[i] - field[0].w_old[i]);
    }
    steady_state_error = sqrt(steady_state_error/myPointStruct[0].num_nodes);
    return steady_state_error;
}

void FS_calculate_intermediate_velocity_vectorised(PointStructure* myPointStruct, FieldVariables* field){
    
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
// x-momentum
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->u, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->u, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    if (parameters.dimension == 3)
        multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->u, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->lap, field->u, field->dpdn, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    
    # pragma acc parallel loop present(field, parameters, myPointStruct)
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        field->u_new[i] = field->u[i] - parameters.dt * (field->u[i] * field->dpdx[i] + field->v[i] * field->dpdy[i] + field->w[i] * field->dpdz[i] - parameters.nu *field->dpdn[i]);

// y-momentum
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->v, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->v, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    if (parameters.dimension == 3)
        multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->v, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->lap, field->v, field->dpdn, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    
    # pragma acc parallel loop present(field, parameters, myPointStruct)
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        field->v_new[i] = field->v[i] - parameters.dt * (field->u[i] * field->dpdx[i] + field->v[i] * field->dpdy[i] + field->w[i] * field->dpdz[i] - parameters.nu * field->dpdn[i]);

//  z-momentum
    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->w, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->w, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->w, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
        multiply_sparse_matrix_vector_vectorised(myPointStruct->lap, field->w, field->dpdn, myPointStruct->cloud_index, myPointStruct->num_nodes,myPointStruct->num_cloud_points);
        
        # pragma acc parallel loop present(field, parameters, myPointStruct)
        for (int i = 0; i < myPointStruct->num_nodes; i++)
            field->w_new[i] = field->w[i] - parameters.dt * (field->u[i] * field->dpdx[i] + field->v[i] * field->dpdy[i] + field->w[i] * field->dpdz[i] - parameters.nu * field->dpdn[i]);
    }
    
    # pragma acc update host(field[0].u_new[:num_nodes], field[0].v_new[:num_nodes])
}

void FS_calculate_mass_residual_vectorised(PointStructure* myPointStruct, FieldVariables* field){
    
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;
    
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->u_new, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->v_new, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    if (parameters.dimension == 3)
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->w_new, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    
    # pragma acc parallel loop present(field, parameters, myPointStruct, temp1, temp2, temp3)
    for (int i = 0; i < num_nodes; i++)
        if (!myPointStruct->boundary_tag[i])
            field->source[i] = parameters.rho*(field->dpdx[i]+field->dpdy[i]+field->dpdz[i])/parameters.dt;
}

void FS_calculate_boundary_dpdn_vectorised(PointStructure* myPointStruct, FieldVariables* field){
    double dpdx, dpdy, dpdz;
    // Copy only boundary normals to gpu
    # pragma acc parallel loop present(field, parameters, myPointStruct, zeros[:myPointStruct->num_nodes])
    for(int i = 0; i < myPointStruct->num_nodes; i++){
        if (myPointStruct->boundary_tag[i]){
            dpdx = (field->u_new[i] - field->u[i]) * parameters.rho/parameters.dt; 
            dpdy = (field->v_new[i] - field->v[i]) * parameters.rho/parameters.dt;
            field->dpdn[i] = dpdx*myPointStruct->x_normal[i] + dpdy*myPointStruct->y_normal[i];
            if (parameters.dimension == 3){
                dpdz = (field->w_new[i] - field->w[i]) * parameters.rho/parameters.dt; 
                field->dpdn[i] += dpdz*myPointStruct->z_normal[i];
            }
        }
    }
}

void FS_multigrid_Poisson_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field){
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

void FS_relaxation_vectorised(PointStructure* mypointstruct, FieldVariables* field){
    # pragma acc loop
    int n = mypointstruct->num_cloud_points;
    for (int iter = 0; iter < parameters.num_relax; iter++){
        # pragma acc parallel loop present(field, parameters, mypointstruct)
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
        for (int i = 0; i < mypointstruct->num_nodes; i++){
            field->p[i] = field->p[i] - pref;
        }
        
        # pragma acc parallel loop present(field, mypointstruct)
        for (int i = 0; i< mypointstruct->num_nodes; i++){
            field->res[i]=field->source[i];
            int k = i*n;
            for (int j = 0; j < n; j++){
                field->res[i] = field->res[i] - mypointstruct->lap[k]*field->p[mypointstruct->cloud_index[k]];
                k += 1;
            }
        }
        // printf("Pressure Residual: %e\n", l2_norm_gen(mypointstruct, field->res, zeros, mypointstruct->num_nodes));
    
        // printf("Pressure Residual: %e\n", l2_norm_gen(mypointstruct, field->p, zeros, mypointstruct->num_nodes));
    }
}

void FS_restrict_residuals_vectorised(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, FieldVariables* field_f, FieldVariables* field_c){
    int n = mypointStruct_f->num_cloud_points;
    # pragma acc parallel loop present(field_f, field_c, mypointStruct_f, mypointStruct_c)
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

void FS_prolongate_corrections_vectorised(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, FieldVariables* field_f, FieldVariables* field_c){
    int n = mypointStruct_c->num_cloud_points;
    # pragma acc parallel loop present(field_f, field_c, mypointStruct_f, mypointStruct_c)
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
    # pragma acc parallel loop present(field_c, mypointStruct_c)
    for (int i = 0; i<mypointStruct_c->num_nodes; i++){
        if (mypointStruct_c->boundary_tag[i]==false){
            field_c->p[i] = 0.0;
        }
    }
}

void FS_calculate_residuals_vectorised(PointStructure* mypointStruct, FieldVariables* field){
    int n = mypointStruct->num_cloud_points;
    for (int i = 0; i < mypointStruct->num_nodes; i++){
        if (!mypointStruct->boundary_tag[i]){
            double sum = 0;
            int k = i*n;
            for (int j = 0; j < n; j++){
                sum += mypointStruct->lap[k]*field->p[mypointStruct->cloud_index[k]];
                k += 1;
            }
            field->res[i] = field->source[i] - sum;
        }
    }
}

void FS_update_velocity_vectorised(PointStructure* myPointStruct, FieldVariables* field)
{
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->p, field->dpdx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->p, field->dpdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    if (parameters.dimension == 3){
        multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->p, field->dpdz, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    }
    // Update Interior nodes
    # pragma acc parallel loop present(field, parameters, myPointStruct, dpdx, dpdy, dpdz)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (!myPointStruct->boundary_tag[i]){
            field->u[i] = field->u_new[i] - parameters.dt * field->dpdx[i]/parameters.rho;
            field->v[i] = field->v_new[i] - parameters.dt * field->dpdy[i]/parameters.rho;
            if (parameters.dimension == 3){
                field->w[i] = field->w_new[i] - parameters.dt * field->dpdz[i]/parameters.rho;
            }
        }
    }
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->u, field->dpdx, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->v, field->dpdy, myPointStruct->cloud_index, myPointStruct->num_nodes, myPointStruct->num_cloud_points);
    
    double sum;
    sum = 0.0;
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        sum += parameters.rho*fabs(field->dpdx[i]+field->dpdy[i]);

    printf("Mass residual: %e\n", (sum)/myPointStruct->num_nodes);
}


void FS_update_boundary_pressure_vectorised(PointStructure* mypointstruct, FieldVariables* field){
    double sumx, sumy, sumz, Ap, term;
    int n = mypointstruct->num_cloud_points;
    if (parameters.neumann_flag_boundary){
    # pragma acc parallel loop present(field, parameters, mypointstruct)
    for (int i = 0; i < mypointstruct->num_boundary_nodes; i++)
    {
        sumx = 0; sumy = 0; sumz = 0; Ap = 0;
        # pragma acc loop reduction(+:sumx, sumy, sumz, Ap)
        int k = i*n;
        for (int j = 1; j < n; j++){
            k += 1;
            sumx += mypointstruct->Dx[k]*field->p[mypointstruct->cloud_index[k]];
            sumy += mypointstruct->Dy[k]*field->p[mypointstruct->cloud_index[k]];
            if (parameters.dimension == 3){
                sumz += mypointstruct->Dz[k]*field->p[mypointstruct->cloud_index[k]];
            }
        }
        k = i*n;
        Ap += mypointstruct->Dx[k]*mypointstruct->x_normal[i];
        Ap += mypointstruct->Dy[k]*mypointstruct->y_normal[i];
        if (parameters.dimension == 3){
            Ap += mypointstruct->Dz[k]*mypointstruct->z_normal[i];
        }
        term = (field->dpdn[i]-sumx*mypointstruct->x_normal[i] -sumy*mypointstruct->y_normal[i] -sumz*mypointstruct->z_normal[i])/Ap;
        field->p[i] = field->p[i] * 0.5 + 0.5 * term;
    }
    }
}