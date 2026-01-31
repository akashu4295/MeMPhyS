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
#include "functions.h"
#include "mat_lib.h"

////////////////////////////////////////////////////////////////////////////////////
// Time-Implicit Solver Modules
////////////////////////////////////////////////////////////////////////////////////

double time_implicit_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field){
    double steady_state_error = 0.0;

    # pragma acc parallel loop present(field, myPointStruct)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        field[0].u_old[i] = field[0].u[i]; field[0].u_new[i] = field[0].u[i];
        field[0].v_old[i] = field[0].v[i]; field[0].v_new[i] = field[0].v[i];
        field[0].w_old[i] = field[0].w[i]; field[0].w_new[i] = field[0].w[i];
        field[0].pprime[i] = 0;
    }
    for (int iter = 0; iter<parameters.iter_timple; iter++){ 
        # pragma acc parallel loop present(field, myPointStruct)
        for (int i = 0; i < myPointStruct->num_nodes; i++){
            field[0].u_new[i] = field[0].u[i];
            field[0].v_new[i] = field[0].v[i];
            field[0].w_new[i] = field[0].w[i];
            field[0].pprime[i] = 0; 
        }
        calculate_intermediate_velocity_implicit_vectorised(myPointStruct, field);
        calculate_mass_residual_implicit_vectorised(myPointStruct, field);
        multigrid_Poisson_solver_vectorised(myPointStruct, field);
        update_velocity_implicit_vectorised(myPointStruct, field);
        update_boundary_pressure_vectorised(myPointStruct, field);
    }

    # pragma acc parallel loop present(field, myPointStruct) reduction(+:steady_state_error)
    for (int i=0; i<myPointStruct[0].num_nodes; i++){
        steady_state_error += pow(field->u[i]-field->u_old[i],2) + pow(field->v[i]-field->v_old[i],2) + pow(field->w[i]-field->w_old[i],2);
    }
    steady_state_error = sqrt(steady_state_error/myPointStruct[0].num_nodes);
    return steady_state_error;
}

double time_implicit_solver_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field){
    double steady_state_error = 0.0;

    # pragma acc parallel loop present(field, myPointStruct)
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        field[0].u_old[i] = field[0].u[i]; 
        field[0].u_new[i] = field[0].u[i];
        field[0].v_old[i] = field[0].v[i]; 
        field[0].v_new[i] = field[0].v[i];
        field[0].pprime[i] = 0;
    }
    
    for (int iter = 0; iter<parameters.iter_timple; iter++){ 
        # pragma acc parallel loop present(field, myPointStruct)
        for (int i = 0; i < myPointStruct->num_nodes; i++){
            field[0].u_new[i] = field[0].u[i];
            field[0].v_new[i] = field[0].v[i];
            field[0].pprime[i] = 0; 
        }
        
        calculate_intermediate_velocity_implicit_vectorised_2d(myPointStruct, field);
        calculate_mass_residual_implicit_vectorised_2d(myPointStruct, field);
        multigrid_Poisson_solver_vectorised(myPointStruct, field);
        update_velocity_implicit_vectorised_2d(myPointStruct, field);
        update_boundary_pressure_vectorised_2d(myPointStruct, field);
    }

    # pragma acc parallel loop present(field, myPointStruct) reduction(+:steady_state_error)
    for (int i=0; i<myPointStruct[0].num_nodes; i++){
        steady_state_error += pow(field->u[i]-field->u_old[i],2) + pow(field->v[i]-field->v_old[i],2);
    }
    steady_state_error = sqrt(steady_state_error/myPointStruct[0].num_nodes);
    return steady_state_error;
}


void calculate_intermediate_velocity_implicit_vectorised(PointStructure* myPointStruct, FieldVariables* field){
    
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    double denom = 0;
    double temp1, temp2, temp3, temp4;
    double advection_x, advection_y, advection_z, advection, diffusion, unsteady;

    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->p, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->p, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->p, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);
  
    
    for (int iter = 0; iter<parameters.iter_momentum; iter++){
        # pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
        for (int i = 0; i < num_nodes; i++){
            int k = i*num_cloud_points;
            if(!myPointStruct->boundary_tag[i]){
        // u-velocity
                temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0;			
                advection_x = field->u_new[i] * myPointStruct->Dx[k];                      
                advection_y = field->v_new[i] * myPointStruct->Dy[k]; 
                advection_z = field->w_new[i] * myPointStruct->Dz[k];
                diffusion   = - parameters.facRe * parameters.mu * myPointStruct->lap[k];
                unsteady    = 1./(parameters.dt*parameters.facdt);			
                denom       = parameters.rho* (advection_x + advection_y + advection_z + unsteady) + diffusion; 
                # pragma acc loop reduction(+:temp1, temp2, temp3, temp4)
                for (int j = 1; j < num_cloud_points; j++){
                    k += 1;
                    temp1  += myPointStruct->Dx[k] * field->u[myPointStruct->cloud_index[k]];
                    temp2  += myPointStruct->Dy[k] * field->u[myPointStruct->cloud_index[k]];
                    temp3  += myPointStruct->Dz[k] * field->u[myPointStruct->cloud_index[k]];
                    temp4  += myPointStruct->lap[k] * field->u[myPointStruct->cloud_index[k]];
                }
                k = i*num_cloud_points;
                advection = field->u_new[i] * temp1 + field->v_new[i] * temp2 + field->w_new[i] * temp3;
                unsteady  = field->u_old[i]/(parameters.dt*parameters.facdt);
                field->u[i] = (parameters.rho * (unsteady - advection) - (parameters.facRe-1) *   parameters.mu* myPointStruct->lap[i*num_cloud_points] * field->u[i] + 
                               parameters.mu*temp4 - field->dpdx[i]) / denom;
        // v-velocity
                temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0;			
	            advection_x = field->u_new[i] * myPointStruct->Dx[k];                      
	            advection_y = field->v_new[i] * myPointStruct->Dy[k]; 
                advection_z = field->w_new[i] * myPointStruct->Dz[k];  
                diffusion   = - parameters.facRe * parameters.mu * myPointStruct->lap[k];
                unsteady    = 1./(parameters.dt*parameters.facdt);			
                denom       = parameters.rho* (advection_x + advection_y + advection_z + unsteady) + diffusion; 
                # pragma acc loop
                int k = i *num_cloud_points;
                # pragma acc loop reduction(+:temp1, temp2, temp3, temp4)
                for (int j = 1; j < num_cloud_points; j++){
                    k +=1;
                    temp1 += myPointStruct->Dx[k] * field->v[myPointStruct->cloud_index[k]];
                    temp2 += myPointStruct->Dy[k] * field->v[myPointStruct->cloud_index[k]];
                    temp3 += myPointStruct->Dz[k] * field->v[myPointStruct->cloud_index[k]];
                    temp4 += myPointStruct->lap[k] * field->v[myPointStruct->cloud_index[k]];
                }
                advection = field->u_new[i] * temp1 + field->v_new[i] * temp2 + field->w_new[i] * temp3;
                unsteady  = field->v_old[i]/(parameters.dt*parameters.facdt);
                field->v[i] = (parameters.rho * (unsteady - advection) - (parameters.facRe-1) * parameters.mu* myPointStruct->lap[i*num_cloud_points] * field->v[i] +
                                    parameters.mu*temp4 - field->dpdy[i]) / denom;

            // w-velocity
                temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0;	
                k = i*num_cloud_points;		
                advection_x = field->u_new[i] * myPointStruct->Dx[k];                      
                advection_y = field->v_new[i] * myPointStruct->Dy[k];
                advection_z = field->w_new[i] * myPointStruct->Dz[k];  
                diffusion   = - parameters.facRe * parameters.mu * myPointStruct->lap[k];
                unsteady    = 1./(parameters.dt*parameters.facdt);			
                denom       = parameters.rho* (advection_x + advection_y + advection_z + unsteady) + diffusion; 
                # pragma acc loop reduction(+:temp1, temp2, temp3, temp4)
                for (int j = 1; j < num_cloud_points; j++){
                    k += 1;
                    temp1 += myPointStruct->Dx[k] * field->w[myPointStruct->cloud_index[k]];
                    temp2 += myPointStruct->Dy[k] * field->w[myPointStruct->cloud_index[k]];
                    temp3 += myPointStruct->Dz[k] * field->w[myPointStruct->cloud_index[k]];
                    temp4 += myPointStruct->lap[k]* field->w[myPointStruct->cloud_index[k]];
                }
                advection = field->u_new[i] * temp1 + field->v_new[i] * temp2 + field->w_new[i] * temp3;
                unsteady  = field->w_old[i]/(parameters.dt*parameters.facdt);
                field->w[i] = (parameters.rho * (unsteady - advection)- (parameters.facRe-1) * parameters.mu* myPointStruct->lap[i*num_cloud_points] * field->w[i] + 
                                parameters.mu*temp4 - field->dpdz[i]) / denom;
            }
        }
    }
}

void calculate_intermediate_velocity_implicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field){
    
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->p, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->p, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    
    for (int iter = 0; iter<parameters.iter_momentum; iter++){
        # pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
        for (int i = 0; i < num_nodes; i++){
            if(!myPointStruct->boundary_tag[i]){
                double temp1 = 0, temp2 = 0, temp4 = 0;
                double advection_x, advection_y, advection, diffusion, unsteady, denom;
                
        // u-velocity
                advection_x = field->u_new[i] * myPointStruct->Dx[i*num_cloud_points];                      
                advection_y = field->v_new[i] * myPointStruct->Dy[i*num_cloud_points]; 
                diffusion   = - parameters.facRe * parameters.mu * myPointStruct->lap[i*num_cloud_points];
                unsteady    = 1./(parameters.dt*parameters.facdt);			
                denom       = parameters.rho* (advection_x + advection_y + unsteady) + diffusion; 
                
                # pragma acc loop reduction(+:temp1, temp2, temp4)
                for (int j = 1; j < num_cloud_points; j++){
                    int k = i*num_cloud_points + j;
                    temp1  += myPointStruct->Dx[k] * field->u[myPointStruct->cloud_index[k]];
                    temp2  += myPointStruct->Dy[k] * field->u[myPointStruct->cloud_index[k]];
                    temp4  += myPointStruct->lap[k]* field->u[myPointStruct->cloud_index[k]];
                }
                
                advection = field->u_new[i] * temp1 + field->v_new[i] * temp2;
                unsteady  = field->u_old[i]/(parameters.dt*parameters.facdt);
                field->u[i] = (parameters.rho * (unsteady - advection) - (parameters.facRe-1) * parameters.mu* myPointStruct->lap[i*num_cloud_points] * field->u[i] + 
                               parameters.mu*temp4 - field->dpdx[i]) / denom;
                
        // v-velocity
                temp1 = 0; temp2 = 0; temp4 = 0;
                advection_x = field->u_new[i] * myPointStruct->Dx[i*num_cloud_points];                      
                advection_y = field->v_new[i] * myPointStruct->Dy[i*num_cloud_points]; 
                diffusion   = - parameters.facRe * parameters.mu * myPointStruct->lap[i*num_cloud_points];
                unsteady    = 1./(parameters.dt*parameters.facdt);			
                denom       = parameters.rho* (advection_x + advection_y + unsteady) + diffusion; 
                
                # pragma acc loop reduction(+:temp1, temp2, temp4)
                for (int j = 1; j < num_cloud_points; j++){
                    int k = i*num_cloud_points + j;
                    temp1 += myPointStruct->Dx[k] * field->v[myPointStruct->cloud_index[k]];
                    temp2 += myPointStruct->Dy[k] * field->v[myPointStruct->cloud_index[k]];
                    temp4 += myPointStruct->lap[k]* field->v[myPointStruct->cloud_index[k]];
                }
                
                advection = field->u_new[i] * temp1 + field->v_new[i] * temp2;
                unsteady  = field->v_old[i]/(parameters.dt*parameters.facdt);
                field->v[i] = (parameters.rho * (unsteady - advection) - (parameters.facRe-1) * parameters.mu* myPointStruct->lap[i*num_cloud_points] * field->v[i] +
                                    parameters.mu*temp4 - field->dpdy[i]) / denom;                  
            }
        }
    }
}


void calculate_mass_residual_implicit_vectorised(PointStructure* myPointStruct, FieldVariables* field){

    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    double sum = 0.0;
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->u, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->v, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->w, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);

    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct) reduction(+:sum)
    for (int i = 0; i < num_nodes; i++){
        if(!myPointStruct->boundary_tag[i]){
            field->source[i] = parameters.rho*(field->dpdx[i]+field->dpdy[i]+field->dpdz[i])/parameters.dt;
            sum += (field->source[i]*parameters.dt);
      	}
      	else{
      	    field->source[i] = 0;
    	}
    }
    printf("Mass residual: %e\n", fabs(sum)/num_nodes);         
    # pragma acc update host(field[0].source[:num_nodes])
}

void calculate_mass_residual_implicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field){

    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    double sum = 0.0;
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->u, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->v, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);

    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct) reduction(+:sum)
    for (int i = 0; i < num_nodes; i++){
        if(!myPointStruct->boundary_tag[i]){
            field->source[i] = parameters.rho*(field->dpdx[i]+field->dpdy[i])/parameters.dt;
            sum += (field->source[i]*parameters.dt);
        }
        else{
            field->source[i] = 0;
        }
    }
    printf("Mass residual: %e\n", fabs(sum)/num_nodes);         
    // Remove this line unless debugging - causes expensive GPU-to-CPU transfer
    // # pragma acc update host(field[0].source[:num_nodes])
}

void multigrid_Poisson_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field){
    for (int icycle = 0; icycle < parameters.num_vcycles; icycle++){
        for (int ilev = 0; ilev < parameters.num_levels; ilev++){
            relaxation_vectorised(&myPointStruct[ilev], &field[ilev]);  
            calculate_residuals_vectorised(&myPointStruct[ilev], &field[ilev]);
            if (ilev != parameters.num_levels-1){
                restrict_residuals_vectorised(&myPointStruct[ilev], &myPointStruct[ilev+1], &field[ilev], &field[ilev+1]);
            }
        }
        for (int ilev = parameters.num_levels-1; ilev > 0; ilev--){
            prolongate_corrections_vectorised(&myPointStruct[ilev-1], &myPointStruct[ilev], &field[ilev-1], &field[ilev]);
            if (ilev != 1)
                relaxation_vectorised(&myPointStruct[ilev-1], &field[ilev-1]);
        }
    }        
}

void update_boundary_pprime_vectorised(PointStructure* mypointstruct, FieldVariables* field){
    double sumx, sumy, sumz, Ap;
    int n = mypointstruct->num_cloud_points;
    if (parameters.neumann_flag_boundary){
        # pragma acc parallel loop gang vector present(field, parameters, mypointstruct)
        for (int i = 0; i < mypointstruct->num_boundary_nodes; i++){
            sumx = 0; sumy = 0; sumz = 0; Ap = 0;
            # pragma acc loop reduction(+:sumx, sumy, sumz)
            for (int j = 1; j < mypointstruct->num_cloud_points; j++){
                int k = i*n + j;
                sumx += mypointstruct->Dx[k]*field->pprime[mypointstruct->cloud_index[k]];
                sumy += mypointstruct->Dy[k]*field->pprime[mypointstruct->cloud_index[k]];
                sumz += mypointstruct->Dz[k]*field->pprime[mypointstruct->cloud_index[k]];
            }
            int k = i*n;
            Ap += mypointstruct->Dx[k]*mypointstruct->x_normal[i];
            Ap += mypointstruct->Dy[k]*mypointstruct->y_normal[i];
            Ap += mypointstruct->Dz[k]*mypointstruct->z_normal[i];
            field->pprime[i] = (-sumx*mypointstruct->x_normal[i] -sumy*mypointstruct->y_normal[i] -sumz*mypointstruct->z_normal[i])/Ap;
        }
    }
}

void update_boundary_pprime_vectorised_2d(PointStructure* mypointstruct, FieldVariables* field){
    int n = mypointstruct->num_cloud_points;
    if (parameters.neumann_flag_boundary){
        # pragma acc parallel loop gang vector present(field, parameters, mypointstruct)
        for (int i = 0; i < mypointstruct->num_boundary_nodes; i++){
            double sumx = 0, sumy = 0, Ap = 0;
            
            # pragma acc loop reduction(+:sumx, sumy)
            for (int j = 1; j < mypointstruct->num_cloud_points; j++){
                int k = i*n + j;
                sumx += mypointstruct->Dx[k]*field->pprime[mypointstruct->cloud_index[k]];
                sumy += mypointstruct->Dy[k]*field->pprime[mypointstruct->cloud_index[k]];
            }
            
            Ap += mypointstruct->Dx[i*n]*mypointstruct->x_normal[i];
            Ap += mypointstruct->Dy[i*n]*mypointstruct->y_normal[i];
            field->pprime[i] = (-sumx*mypointstruct->x_normal[i] -sumy*mypointstruct->y_normal[i])/Ap;
        }
    }
}

void relaxation_vectorised(PointStructure* mypointstruct, FieldVariables* field){
     //double* zeros=create_vector(mypointstruct->num_nodes);
    int n = mypointstruct->num_cloud_points;
    //  # pragma acc loop
    for (int iter = 0; iter < parameters.num_relax; iter++){
        # pragma acc parallel loop gang vector present(field, parameters, mypointstruct)
        for (int i = 0; i < mypointstruct->num_nodes; i++){
            double sum = 0.0;
            int k = i*n;
            # pragma acc loop reduction(+:sum)
            for (int j = 1; j < n; j++){
                k += 1;
                sum += mypointstruct->lap_Poison[k]*field->pprime[mypointstruct->cloud_index[k]];
            }
            field->pprime[i] = parameters.omega*((field->source[i]-sum)/mypointstruct->lap_Poison[i*n]) + (1-parameters.omega)*field->pprime[i]; 
        }
        double pref = field->pprime[0];
        # pragma acc parallel loop gang vector present(field, mypointstruct)
        for (int i = 0; i < mypointstruct->num_nodes; i++){
            field->pprime[i] = field->pprime[i] - pref;
        }
        
        # pragma acc parallel loop gang vector present(field, mypointstruct)
        for (int i = 0; i< mypointstruct->num_nodes; i++){
            if(!mypointstruct->boundary_tag[i]){
                double res_val = field->source[i];
                # pragma acc loop reduction(-:res_val)
                for (int j = 0; j < mypointstruct->num_cloud_points; j++){
                    res_val -= mypointstruct->lap_Poison[i*n + j]*field->pprime[mypointstruct->cloud_index[i*n + j]];
                }
                field->res[i] = res_val;
            }
            else{
                field->res[i] = 0;
            }
        }
           //printf("Iteration: %d, Pressure Residual: %e\n", iter, l2_norm_gen(mypointstruct, field->res, zeros, mypointstruct->num_nodes));    
    }
           // printf(" Pressure Residual: %e\n", l2_norm_gen(mypointstruct, field->res, zeros, mypointstruct->num_nodes));    
           //printf("Pressure corrections: %e\n", l2_norm_gen(mypointstruct, field->pprime, zeros, mypointstruct->num_nodes));
    //free(zeros);
}

void calculate_residuals_vectorised(PointStructure* mypointStruct, FieldVariables* field){
    // double* zeros=create_vector(mypointStruct->num_nodes);
    int n = mypointStruct->num_cloud_points;
    # pragma acc parallel loop gang vector present(field, mypointStruct)
    for (int i = 0; i < mypointStruct->num_nodes; i++){
           field->res[i] = 0.0;
        if (!mypointStruct->boundary_tag[i]){
            double sum = 0;
            # pragma acc loop reduction(+:sum)
            for (int j = 0; j < n; j++){
                sum += mypointStruct->lap_Poison[i*n +j]*field->pprime[mypointStruct->cloud_index[i*n +j]];
            }
            field->res[i] = field->source[i] - sum;
        }
    }
     //printf("Pressure Residual: %e\n", l2_norm_gen(mypointStruct, field->res, zeros, mypointStruct->num_nodes));
    //  free(zeros);
}

void restrict_residuals_vectorised(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, 
                                                    FieldVariables* field_f, FieldVariables* field_c){
    int n = mypointStruct_f->num_cloud_points;
    # pragma acc parallel loop gang vector present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_c->num_nodes; i++){
        double results = 0.0;
        if (mypointStruct_c->boundary_tag[i]==false){
            int k = i*n;
            int i_restr_n = mypointStruct_c->restriction_points[i]*n;
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < mypointStruct_f->num_cloud_points; j++){
                results += mypointStruct_c->restr_mat[k] * field_f->res[mypointStruct_f->cloud_index[i_restr_n + j]];
                k += 1;
            }
        }
                field_c->source[i] = results;
    }
}

void prolongate_corrections_vectorised(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, FieldVariables* field_f, FieldVariables* field_c){
    int n = mypointStruct_c->num_cloud_points;
    # pragma acc parallel loop gang vector present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_f->num_nodes; i++){
        if (mypointStruct_f->boundary_tag[i]==false){
            int i_prol_n = mypointStruct_f->prolongation_points[i]*n;
            int k = i*n;
            double results = 0.0;
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < n; j++){
                results += mypointStruct_f->prol_mat[k] * field_c->pprime[mypointStruct_c->cloud_index[i_prol_n +j]];
                k += 1;
            }
            field_f->pprime[i] = field_f->pprime[i] + results;
        }
    }
    # pragma acc parallel loop gang vector present(field_c, mypointStruct_c)
    for (int i = 0; i<mypointStruct_c->num_nodes; i++){
       // if (mypointStruct_c->boundary_tag[i]==false){
            field_c->pprime[i] = 0.0;
       // }
    }
}

void update_velocity_implicit_vectorised(PointStructure* myPointStruct, FieldVariables* field){
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->pprime, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->pprime, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dz, field->pprime, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points);

    // Update Interior nodes
    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
    for (int i = 0; i < num_nodes; i++){
        if(!myPointStruct->boundary_tag[i]){
            field->u[i] = field->u[i] - parameters.dt * field->dpdx[i]/parameters.rho;
            field->v[i] = field->v[i] - parameters.dt * field->dpdy[i]/parameters.rho;
            field->w[i] = field->w[i] - parameters.dt * field->dpdz[i]/parameters.rho;
	        field->p[i] = field->p[i] + 1.0*field->pprime[i];
        }
    }
    double pref = field->p[0];
    # pragma acc parallel loop gang vector present(field, myPointStruct)
    for (int i = 0; i < num_nodes; i++){
        field->p[i] = field->p[i] - pref;
    }
}

void update_velocity_implicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field){
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dx, field->pprime, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points);
    multiply_sparse_matrix_vector_vectorised(myPointStruct->Dy, field->pprime, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points);

    // Update Interior nodes
    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
    for (int i = 0; i < num_nodes; i++){
        if(!myPointStruct->boundary_tag[i]){
            field->u[i] = field->u[i] - parameters.dt * field->dpdx[i]/parameters.rho;
            field->v[i] = field->v[i] - parameters.dt * field->dpdy[i]/parameters.rho;
            field->p[i] = field->p[i] + 1.0*field->pprime[i];
        }
    }
    
    double pref = field->p[0];
    # pragma acc parallel loop gang vector present(field, myPointStruct)
    for (int i = 0; i < num_nodes; i++){
        field->p[i] = field->p[i] - pref;
    }
}

void update_boundary_pressure_vectorised(PointStructure* myPointStruct, FieldVariables* field){
    double temp1, temp2, temp3;
    double sumx, sumy, sumz, Ap;
    double nu = parameters.mu / parameters.rho;
    int n = myPointStruct->num_cloud_points;

    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
    for(int i = 0; i < myPointStruct->num_nodes; i++){
        if(myPointStruct->boundary_tag[i]){
            int k = i*n;
            temp1 = 0; temp2 = 0; temp3 = 0;
            # pragma acc loop reduction(+:temp1, temp2, temp3)
            for(int j = 0; j<myPointStruct->num_cloud_points; j++){
                temp1 += nu * myPointStruct->lap[k] * field->u_new[myPointStruct->cloud_index[k]];
                temp1 -= field->u_new[i]* field->u[myPointStruct->cloud_index[k]] * myPointStruct->Dx[k];
                temp1 -= field->v_new[i]* field->u[myPointStruct->cloud_index[k]] * myPointStruct->Dy[k];
                temp1 -= field->w_new[i]* field->u[myPointStruct->cloud_index[k]] * myPointStruct->Dz[k];
                
                temp2 += nu* myPointStruct->lap[k] * field->v[myPointStruct->cloud_index[k]];
                temp2 -= field->u_new[i]* field->v[myPointStruct->cloud_index[k]] * myPointStruct->Dx[k];
                temp2 -= field->v_new[i]* field->v[myPointStruct->cloud_index[k]] * myPointStruct->Dy[k];     
                temp2 -= field->w_new[i]* field->v[myPointStruct->cloud_index[k]] * myPointStruct->Dz[k];
                
                temp3 += nu* myPointStruct->lap[k] * field->w[myPointStruct->cloud_index[k]];
                temp3 -= field->u_new[i]* field->w[myPointStruct->cloud_index[k]] * myPointStruct->Dx[k];
                temp3 -= field->v_new[i]* field->w[myPointStruct->cloud_index[k]] * myPointStruct->Dy[k];
                temp3 -= field->w_new[i]* field->w[myPointStruct->cloud_index[k]] * myPointStruct->Dz[k];
                k +=1;
            }
        field->dpdn[i] = parameters.rho * (temp1*myPointStruct->x_normal[i] + temp2*myPointStruct->y_normal[i] + temp3*myPointStruct->z_normal[i]);
		
        sumx = 0; sumy = 0; sumz = 0; Ap = 0.0;
        # pragma acc loop reduction(+:sumx, sumy, sumz)
        for (int j = 1; j < myPointStruct->num_cloud_points; j++){
            k = i*n + j;
            sumx += myPointStruct->Dx[k]*field->p[myPointStruct->cloud_index[k]];
            sumy += myPointStruct->Dy[k]*field->p[myPointStruct->cloud_index[k]];
            sumz += myPointStruct->Dz[k]*field->p[myPointStruct->cloud_index[k]];
        }
        k = i*n;    
        Ap += myPointStruct->Dx[k]*myPointStruct->x_normal[i];
        Ap += myPointStruct->Dy[k]*myPointStruct->y_normal[i];
        Ap += myPointStruct->Dz[k]*myPointStruct->z_normal[i];
        field->p[i] =(field->dpdn[i] - sumx*myPointStruct->x_normal[i] -sumy*myPointStruct->y_normal[i] - sumz*myPointStruct->z_normal[i])/Ap;
    	}
    }
}

void update_boundary_pressure_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field){
    double nu = parameters.mu / parameters.rho;
    int n = myPointStruct->num_cloud_points;

    # pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
    for(int i = 0; i < myPointStruct->num_nodes; i++){
        if(myPointStruct->boundary_tag[i]){
            double temp1 = 0, temp2 = 0;
            
            # pragma acc loop reduction(+:temp1, temp2)
            for(int j = 0; j<myPointStruct->num_cloud_points; j++){
                int k = i*n + j;
                double t1_contrib = nu * myPointStruct->lap[k] * field->u_new[myPointStruct->cloud_index[k]];
                t1_contrib -= field->u_new[i]* field->u[myPointStruct->cloud_index[k]] * myPointStruct->Dx[k];
                t1_contrib -= field->v_new[i]* field->u[myPointStruct->cloud_index[k]] * myPointStruct->Dy[k];
                
                double t2_contrib = nu* myPointStruct->lap[k] * field->v[myPointStruct->cloud_index[k]];
                t2_contrib -= field->u_new[i]* field->v[myPointStruct->cloud_index[k]] * myPointStruct->Dx[k];
                t2_contrib -= field->v_new[i]* field->v[myPointStruct->cloud_index[k]] * myPointStruct->Dy[k];     
                
                temp1 += t1_contrib;
                temp2 += t2_contrib;
            }
            field->dpdn[i] = parameters.rho * (temp1*myPointStruct->x_normal[i] + temp2*myPointStruct->y_normal[i]);

            double sumx = 0, sumy = 0, Ap = 0.0;
            # pragma acc loop reduction(+:sumx, sumy)
            for (int j = 1; j < myPointStruct->num_cloud_points; j++){
                int k = i*n + j;
                sumx += myPointStruct->Dx[k]*field->p[myPointStruct->cloud_index[k]];
                sumy += myPointStruct->Dy[k]*field->p[myPointStruct->cloud_index[k]];
            }
            
            Ap += myPointStruct->Dx[i*n]*myPointStruct->x_normal[i];
            Ap += myPointStruct->Dy[i*n]*myPointStruct->y_normal[i];
            field->p[i] =(field->dpdn[i] - sumx*myPointStruct->x_normal[i] -sumy*myPointStruct->y_normal[i])/Ap;
        }
    }
}