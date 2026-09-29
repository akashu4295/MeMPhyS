// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#include "functions.h"

void solve_momentum(PointStructure* pStruct, FieldVariables* field);
void solve_pprime(PointStructure* pStruct, FieldVariables* field);
void apply_smoother(PointStructure* pStruct, FieldVariables* field);

void restrict_fas(PointStructure* fineStruct, PointStructure* coarseStruct, FieldVariables* field_f, FieldVariables* field_c);
void calculate_nonlinear_residual(PointStructure* myPointStruct, FieldVariables* field);
void restrict_solution_fas(PointStructure* fineStruct, PointStructure* coarseStruct, FieldVariables* field_f, FieldVariables* field_c);
void restrict_residual_fas(PointStructure* fineStruct, PointStructure* coarseStruct, FieldVariables* field_f, FieldVariables* field_c);
void calculate_coarse_rhs_fas(PointStructure* coarseStruct, FieldVariables* field_c);
void prolongate_corrections_fas(PointStructure* fineStruct, PointStructure* coarseStruct, FieldVariables* field_f, FieldVariables* field_c);

// Main Time-Implicit Solver (FAS Multigrid)
// ============================================================================

double multigrid_time_implicit_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field){
    int num_levels = parameters.num_levels; 

    if (parameters.dimension == 3) {
        #pragma acc parallel loop present(field[0], myPointStruct[0])
        for (int i = 0; i < myPointStruct[0].num_nodes; i++) {
            field[0].u_old[i] = field[0].u[i]; 
            field[0].v_old[i] = field[0].v[i]; 
            field[0].w_old[i] = field[0].w[i]; 
        }
    } else {
        #pragma acc parallel loop present(field[0], myPointStruct[0])
        for (int i = 0; i < myPointStruct[0].num_nodes; i++) {
            field[0].u_old[i] = field[0].u[i]; 
            field[0].v_old[i] = field[0].v[i]; 
        }
    }
    
    #pragma acc data present(field[:num_levels], myPointStruct[:num_levels], parameters)
    {
        for (int iter = 0; iter < parameters.iter_timple; iter++){ 
            
            for (int ilev = 0; ilev < num_levels - 1; ilev++) {
                apply_smoother(&myPointStruct[ilev], &field[ilev]);
                restrict_fas(&myPointStruct[ilev], &myPointStruct[ilev+1], &field[ilev], &field[ilev+1]);
            }

            int coarsest = num_levels - 1;
            apply_smoother(&myPointStruct[coarsest], &field[coarsest]);

            for (int ilev = coarsest; ilev > 0; ilev--) {
                prolongate_corrections_fas(&myPointStruct[ilev-1], &myPointStruct[ilev], &field[ilev-1], &field[ilev]);
                apply_smoother(&myPointStruct[ilev-1], &field[ilev-1]);
            }
        }
    }

    double steady_state_error_par = 0.0;    
    
    if (parameters.dimension == 3) {
        #pragma acc parallel loop present(field[0], myPointStruct[0]) reduction(+:steady_state_error_par) 
        for (int i = 0; i < myPointStruct[0].num_nodes; i++){
            double du = field[0].u[i] - field[0].u_old[i];
            double dv = field[0].v[i] - field[0].v_old[i];
            double dw = field[0].w[i] - field[0].w_old[i];
            steady_state_error_par += fabs(du) + fabs(dv) + fabs(dw);
        }
    } else {
        #pragma acc parallel loop present(field[0], myPointStruct[0]) reduction(+:steady_state_error_par) 
        for (int i = 0; i < myPointStruct[0].num_nodes; i++){
            double du = field[0].u[i] - field[0].u_old[i];
            double dv = field[0].v[i] - field[0].v_old[i];
            steady_state_error_par += fabs(du) + fabs(dv);
        }
    }
    
    return (steady_state_error_par/(parameters.dimension * myPointStruct[0].num_nodes * parameters.dt));
}


void solve_momentum(PointStructure* pStruct, FieldVariables* field) {
    if (parameters.dimension == 3) {
        calculate_intermediate_velocity_implicit_vectorised(pStruct, field);
    } else {
        calculate_intermediate_velocity_implicit_vectorised_2d(pStruct, field);
    }
}

void solve_pprime(PointStructure* pStruct, FieldVariables* field) {
    if (parameters.dimension == 3) {
        calculate_mass_residual_implicit_vectorised(pStruct, field);
        single_grid_Poisson_solver_vectorised(pStruct, field);
        update_velocity_implicit_vectorised(pStruct, field);
        update_boundary_pressure_vectorised(pStruct, field);
    } else {
        calculate_mass_residual_implicit_vectorised_2d(pStruct, field);
        single_grid_Poisson_solver_vectorised(pStruct, field);
        update_velocity_implicit_vectorised_2d(pStruct, field);
        update_boundary_pressure_vectorised_2d(pStruct, field);
    }
}

void apply_smoother(PointStructure* pStruct, FieldVariables* field) {
    solve_momentum(pStruct, field);
    solve_pprime(pStruct, field);
}

void single_grid_Poisson_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field) {
    if (parameters.poisson_solver_type == 1)    
        relaxation_vectorised_Jacobi(myPointStruct, field->source, field->pprime, field->p_old);
    else if (parameters.poisson_solver_type == 3)
        relaxation_vectorised_BiCGStab(myPointStruct, field->source, field->pprime, parameters.num_relax, parameters.poisson_solver_tolerance);
    else
        relaxation_vectorised_GaussSeidel(myPointStruct, field->source, field->pprime);  
}

// FAS Multigrid Subroutines
// ============================================================================

void restrict_fas(PointStructure* fineStruct, PointStructure* coarseStruct, FieldVariables* field_f, FieldVariables* field_c) {
    calculate_nonlinear_residual(fineStruct, field_f);
    restrict_solution_fas(fineStruct, coarseStruct, field_f, field_c);
    restrict_residual_fas(fineStruct, coarseStruct, field_f, field_c);
    calculate_coarse_rhs_fas(coarseStruct, field_c);
}

void restrict_solution_fas(PointStructure* fineStruct, PointStructure* coarseStruct, FieldVariables* field_f, FieldVariables* field_c) {
    int n_f = fineStruct->num_cloud_points;
    int num_nodes_c = coarseStruct->num_nodes;

    #pragma acc parallel loop gang vector present(field_f, field_c, fineStruct, coarseStruct)
    for (int i = 0; i < num_nodes_c; i++) {
        if (!coarseStruct->boundary_tag[i]) {
            int base_c = i * n_f; 
            int i_restr_node = coarseStruct->restriction_points[i];
            int base_f = i_restr_node * n_f;
            
            double u_res = 0.0, v_res = 0.0, w_res = 0.0, p_res = 0.0;

            #pragma acc loop seq
            for (int j = 0; j < n_f; j++) {
                double weight = coarseStruct->restr_mat[base_c + j];
                int idx = fineStruct->cloud_index[base_f + j];
                u_res += weight * field_f->u[idx];
                v_res += weight * field_f->v[idx];
                p_res += weight * field_f->p[idx];
                if (parameters.dimension == 3) {
                    w_res += weight * field_f->w[idx];
                }
            }
            // Save original restricted solution to compute \delta\Phi^H later
            field_c->u_restricted[i] = u_res; field_c->u[i] = u_res;
            field_c->v_restricted[i] = v_res; field_c->v[i] = v_res;
            if (parameters.dimension == 3) {
                field_c->w_restricted[i] = w_res; field_c->w[i] = w_res;
            }
            field_c->p[i] = p_res;
        }
    }
}

void prolongate_corrections_fas(PointStructure* fineStruct, PointStructure* coarseStruct, FieldVariables* field_f, FieldVariables* field_c) {
    int n_c = coarseStruct->num_cloud_points;
    int num_nodes_f = fineStruct->num_nodes;

    #pragma acc parallel loop gang vector present(field_f, field_c, fineStruct, coarseStruct)
    for (int i = 0; i < num_nodes_f; i++) {
        if (!fineStruct->boundary_tag[i]) {
            int i_prol_node = fineStruct->prolongation_points[i];
            int base_f = i * n_c;
            int base_c = i_prol_node * n_c;
            
            double du = 0.0, dv = 0.0, dw = 0.0;

            #pragma acc loop seq
            for (int j = 0; j < n_c; j++) {
                double weight = fineStruct->prol_mat[base_f + j];
                int idx = coarseStruct->cloud_index[base_c + j];
                // \delta\Phi^H = \Phi^H_{new} - I_h^H\Phi^h
                du += weight * (field_c->u[idx] - field_c->u_restricted[idx]);
                dv += weight * (field_c->v[idx] - field_c->v_restricted[idx]);
                if (parameters.dimension == 3) {
                    dw += weight * (field_c->w[idx] - field_c->w_restricted[idx]);
                }
            }
            
            field_f->u[i] += du;
            field_f->v[i] += dv;
            if (parameters.dimension == 3) {
                field_f->w[i] += dw;
            }
        }
    }
}

void calculate_nonlinear_residual(PointStructure* myPointStruct, FieldVariables* field) {
    int num_nodes = myPointStruct->num_nodes;
    int num_cloud_points = myPointStruct->num_cloud_points;

    #pragma acc data present(field, parameters, myPointStruct)
    {
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->p, field->dpdx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 1);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->p, field->dpdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 2);
        if (parameters.dimension == 3) {
            multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dz, field->p, field->dpdz, myPointStruct->cloud_index, num_nodes, num_cloud_points, 3);
        }
        
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dx, field->u, field->dudx, myPointStruct->cloud_index, num_nodes, num_cloud_points, 4);
        multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dy, field->v, field->dvdy, myPointStruct->cloud_index, num_nodes, num_cloud_points, 5);
        if (parameters.dimension == 3) {
            multiply_sparse_matrix_vector_vectorised_gpu_async(myPointStruct->Dz, field->w, field->dwdz, myPointStruct->cloud_index, num_nodes, num_cloud_points, 6);
        }
    }
    #pragma acc wait(1,2,3,4,5,6)

    #pragma acc parallel loop gang vector present(field, parameters, myPointStruct)
    for (int i = 0; i < num_nodes; i++) {
        if (!myPointStruct->boundary_tag[i] && !myPointStruct->corner_tag[i]) {
            double t1=0, t2=0, t3=0, t4=0;  // u advection/diffusion terms
            double t5=0, t6=0, t7=0, t8=0;  // v advection/diffusion terms
            double t9=0, t10=0, t11=0, t12=0; // w advection/diffusion terms
            int base = i * num_cloud_points;

            #pragma acc loop seq
            for (int j = 0; j < num_cloud_points; j++) {
                int k = base + j;
                int idx = myPointStruct->cloud_index[k];
                
                // U terms
                t1 += myPointStruct->Dx[k] * field->u[idx];
                t2 += myPointStruct->Dy[k] * field->u[idx];
                t4 += myPointStruct->lap[k] * field->u[idx];
                
                // V terms
                t5 += myPointStruct->Dx[k] * field->v[idx];
                t6 += myPointStruct->Dy[k] * field->v[idx];
                t8 += myPointStruct->lap[k] * field->v[idx];
                
                // W terms
                if (parameters.dimension == 3) {
                    t3 += myPointStruct->Dz[k] * field->u[idx];
                    t7 += myPointStruct->Dz[k] * field->v[idx];
                    t9 += myPointStruct->Dx[k] * field->w[idx];
                    t10 += myPointStruct->Dy[k] * field->w[idx];
                    t11 += myPointStruct->Dz[k] * field->w[idx];
                    t12 += myPointStruct->lap[k] * field->w[idx];
                }
            }

            // Construct Nonlinear Operator L^h(u,v,w,p)
            
            // X-Momentum operator
            double unst_u = parameters.rho * (field->u[i] - field->u_old[i]) / parameters.dt;
            double adv_u  = parameters.rho * (field->u[i] * t1 + field->v[i] * t2);
            adv_u += parameters.rho * (parameters.dimension == 3 ? field->w[i] * t3 : 0.0);
            double diff_u = parameters.mu * t4;
            double L_u    = unst_u + adv_u - diff_u + field->dpdx[i];
            
            // Y-Momentum operator
            double unst_v = parameters.rho * (field->v[i] - field->v_old[i]) / parameters.dt;
            double adv_v  = parameters.rho * (field->u[i] * t5 + field->v[i] * t6);
            adv_v += parameters.rho * (parameters.dimension == 3 ? field->w[i] * t7 : 0.0);
            double diff_v = parameters.mu * t8;
            double L_v    = unst_v + adv_v - diff_v + field->dpdy[i];

            // Continuity operator 
            double L_p = parameters.rho * (field->dudx[i] + field->dvdy[i]) / parameters.dt;
            
            field->res_u[i] = field->source_u[i] - L_u;
            field->res_v[i] = field->source_v[i] - L_v;

            if (parameters.dimension == 3) {
                double unst_w = parameters.rho * (field->w[i] - field->w_old[i]) / parameters.dt;
                double adv_w  = parameters.rho * (field->u[i] * t9 + field->v[i] * t10 + field->w[i] * t11);
                double diff_w = parameters.mu * t12;
                double L_w    = unst_w + adv_w - diff_w + field->dpdz[i];
                
                L_p += parameters.rho * field->dwdz[i] / parameters.dt;
                field->res_w[i] = field->source_w[i] - L_w;
            }

            field->res_p[i] = field->source_p[i] - L_p; 
            
        } else {
            field->res_u[i] = 0.0;
            field->res_v[i] = 0.0;
            if (parameters.dimension == 3) {
                field->res_w[i] = 0.0;
            }
            field->res_p[i] = 0.0;
        }
    }
}

void restrict_residual_fas(PointStructure* fineStruct, PointStructure* coarseStruct, FieldVariables* field_f, FieldVariables* field_c) {
    int n_f = fineStruct->num_cloud_points;
    int num_nodes_c = coarseStruct->num_nodes;

    #pragma acc parallel loop gang vector present(field_f, field_c, fineStruct, coarseStruct)
    for (int i = 0; i < num_nodes_c; i++) {
        if (!coarseStruct->boundary_tag[i] && !coarseStruct->corner_tag[i]) {
            int base_c = i * n_f; 
            int i_restr_node = coarseStruct->restriction_points[i];
            int base_f = i_restr_node * n_f;
            
            double ru = 0.0, rv = 0.0, rw = 0.0, rp = 0.0;

            #pragma acc loop seq
            for (int j = 0; j < n_f; j++) {
                double weight = coarseStruct->restr_mat[base_c + j];
                int idx = fineStruct->cloud_index[base_f + j];
                
                ru += weight * field_f->res_u[idx];
                rv += weight * field_f->res_v[idx];
                if (parameters.dimension == 3) {
                    rw += weight * field_f->res_w[idx];
                }
                rp += weight * field_f->res_p[idx];
            }
            
            // Store the restricted fine residual I_h^H R^h
            field_c->res_u_restricted[i] = ru;
            field_c->res_v_restricted[i] = rv;
            if (parameters.dimension == 3) {
                field_c->res_w_restricted[i] = rw;
            }
            field_c->res_p_restricted[i] = rp;
        }
    }
}

void calculate_coarse_rhs_fas(PointStructure* coarseStruct, FieldVariables* field_c) {
    int num_nodes = coarseStruct->num_nodes;
    
    calculate_nonlinear_residual(coarseStruct, field_c);
    
    #pragma acc parallel loop gang vector present(field_c, coarseStruct)
    for (int i = 0; i < num_nodes; i++) {
        if (!coarseStruct->boundary_tag[i] && !coarseStruct->corner_tag[i]) {
            // New Source = Old Source + Restricted Fine Residual - Coarse Residual
            field_c->source_u[i] = field_c->source_u[i] + field_c->res_u_restricted[i] - field_c->res_u[i];
            field_c->source_v[i] = field_c->source_v[i] + field_c->res_v_restricted[i] - field_c->res_v[i];
            if (parameters.dimension == 3) {
                field_c->source_w[i] = field_c->source_w[i] + field_c->res_w_restricted[i] - field_c->res_w[i];
            }
            field_c->source_p[i] = field_c->source_p[i] + field_c->res_p_restricted[i] - field_c->res_p[i];
            
            // Update the mass source explicitly if your TIMPLE solver reads from 'field->source'
            field_c->source[i] = field_c->source_p[i]; 
        }
    }
}