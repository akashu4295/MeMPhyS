#include "structures.h"
#include <stdlib.h>
#include <stdbool.h>

void AllocateMemoryPointStructure(PointStructure* myPointStruct, int nodes) {
    myPointStruct->x_normal = (double*)malloc(nodes * sizeof(double));
    myPointStruct->y_normal = (double*)malloc(nodes * sizeof(double));
    myPointStruct->z_normal = (double*)malloc(nodes * sizeof(double));
    myPointStruct->x = (double*)malloc(nodes * sizeof(double));
    myPointStruct->y = (double*)malloc(nodes * sizeof(double));
    myPointStruct->z = (double*)malloc(nodes * sizeof(double));
    myPointStruct->point_index = (int*)malloc(nodes * sizeof(int));
    myPointStruct->boundary_tag = (bool*)malloc(nodes * sizeof(bool));
    myPointStruct->corner_tag = (bool*)malloc(nodes * sizeof(bool));
    myPointStruct->node_bc = (BCValue*)malloc(nodes * sizeof(BCValue));
    myPointStruct->pow_x = (short*)malloc(myPointStruct->num_poly_terms * sizeof(short));
    myPointStruct->pow_y = (short*)malloc(myPointStruct->num_poly_terms * sizeof(short));
    myPointStruct->pow_z = (short*)malloc(myPointStruct->num_poly_terms * sizeof(short));
    myPointStruct->num_nodes = nodes;
    myPointStruct->num_elem = 0;
    myPointStruct->rcm_order = (int*)malloc(nodes*sizeof(int));
}

void AllocateMemoryFieldVariables(FieldVariables** field, PointStructure* myPointStruct, int num_levels) {
    *field = (FieldVariables*) malloc(num_levels * sizeof(FieldVariables));
    
    for (int ii = 0; ii < num_levels; ii++) {
        int N = myPointStruct[ii].num_nodes;
        
        // ========== INCOMPRESSIBLE FLOW VARIABLES (ALWAYS ALLOCATED) ==========
        (*field)[ii].u = (double*) malloc(N * sizeof(double));
        (*field)[ii].v = (double*) malloc(N * sizeof(double));
        (*field)[ii].u_new = (double*) malloc(N * sizeof(double));
        (*field)[ii].v_new = (double*) malloc(N * sizeof(double));
        (*field)[ii].u_old = (double*) malloc(N * sizeof(double));
        (*field)[ii].v_old = (double*) malloc(N * sizeof(double));
        (*field)[ii].p = (double*) malloc(N * sizeof(double));
        (*field)[ii].p_old = (double*) malloc(N * sizeof(double));
        (*field)[ii].pprime = (double*) malloc(N * sizeof(double));
        (*field)[ii].res = (double*) malloc(N * sizeof(double));
        (*field)[ii].source = (double*) malloc(N * sizeof(double));
        (*field)[ii].dpdn = (double*) malloc(N * sizeof(double));
        (*field)[ii].dpdx = (double*) malloc(N * sizeof(double));
        (*field)[ii].dpdy = (double*) malloc(N * sizeof(double));
        (*field)[ii].dudx = (double*) malloc(N * sizeof(double));
        (*field)[ii].dudy = (double*) malloc(N * sizeof(double));
        (*field)[ii].dvdx = (double*) malloc(N * sizeof(double));
        (*field)[ii].dvdy = (double*) malloc(N * sizeof(double));
        (*field)[ii].lapu = (double*) malloc(N * sizeof(double));
        (*field)[ii].lapv = (double*) malloc(N * sizeof(double));
        
        if (parameters.dimension == 3) {
            (*field)[ii].w = (double*) malloc(N * sizeof(double));
            (*field)[ii].w_new = (double*) malloc(N * sizeof(double));
            (*field)[ii].w_old = (double*) malloc(N * sizeof(double));
            (*field)[ii].dpdz = (double*) malloc(N * sizeof(double));
            (*field)[ii].dudz = (double*) malloc(N * sizeof(double));
            (*field)[ii].dvdz = (double*) malloc(N * sizeof(double));
            (*field)[ii].dwdx = (double*) malloc(N * sizeof(double));
            (*field)[ii].dwdy = (double*) malloc(N * sizeof(double));
            (*field)[ii].dwdz = (double*) malloc(N * sizeof(double));
            (*field)[ii].lapw = (double*) malloc(N * sizeof(double));
        }
        else {
            (*field)[ii].w = NULL;
            (*field)[ii].w_new = NULL;
            (*field)[ii].w_old = NULL;
            (*field)[ii].dpdz = NULL;
            (*field)[ii].dudz = NULL;
            (*field)[ii].dvdz = NULL;
            (*field)[ii].dwdx = NULL;
            (*field)[ii].dwdy = NULL;
            (*field)[ii].dwdz = NULL;
            (*field)[ii].lapw = NULL;
        }
        
        // ========== COMPRESSIBLE FLOW VARIABLES (CONDITIONAL) ==========
        if (parameters.compressible_flow) {
            // Primary thermodynamic variables
            (*field)[ii].rho = (double*) malloc(N * sizeof(double));
            (*field)[ii].rho_old = (double*) malloc(N * sizeof(double));
            (*field)[ii].rho_new = (double*) malloc(N * sizeof(double));
            (*field)[ii].T = (double*) malloc(N * sizeof(double));
            (*field)[ii].T_old = (double*) malloc(N * sizeof(double));
            (*field)[ii].T_new = (double*) malloc(N * sizeof(double));
            (*field)[ii].e = (double*) malloc(N * sizeof(double));
            (*field)[ii].e_old = (double*) malloc(N * sizeof(double));
            
            // Gradient variables
            (*field)[ii].drhodx = (double*) malloc(N * sizeof(double));
            (*field)[ii].drhody = (double*) malloc(N * sizeof(double));
            (*field)[ii].drhodz = (double*) malloc(N * sizeof(double));
            (*field)[ii].dTdx = (double*) malloc(N * sizeof(double));
            (*field)[ii].dTdy = (double*) malloc(N * sizeof(double));
            (*field)[ii].dTdz = (double*) malloc(N * sizeof(double));
            (*field)[ii].dedx = (double*) malloc(N * sizeof(double));
            (*field)[ii].dedy = (double*) malloc(N * sizeof(double));
            (*field)[ii].dedz = (double*) malloc(N * sizeof(double));
            
            // Transport properties
            (*field)[ii].mu = (double*) malloc(N * sizeof(double));
            (*field)[ii].kappa = (double*) malloc(N * sizeof(double));
            
            // Stress tensor components
            (*field)[ii].tau_xx = (double*) malloc(N * sizeof(double));
            (*field)[ii].tau_yy = (double*) malloc(N * sizeof(double));
            (*field)[ii].tau_zz = (double*) malloc(N * sizeof(double));
            (*field)[ii].tau_xy = (double*) malloc(N * sizeof(double));
            (*field)[ii].tau_xz = (double*) malloc(N * sizeof(double));
            (*field)[ii].tau_yz = (double*) malloc(N * sizeof(double));
            (*field)[ii].div_tau_x = (double*) malloc(N * sizeof(double));
            (*field)[ii].div_tau_y = (double*) malloc(N * sizeof(double));
            (*field)[ii].div_tau_z = (double*) malloc(N * sizeof(double));
            
            // Source terms
            (*field)[ii].Q_visc = (double*) malloc(N * sizeof(double));
            (*field)[ii].Q_source = (double*) malloc(N * sizeof(double));
        } else {
            // Set pointers to NULL for incompressible flow
            (*field)[ii].rho = NULL;
            (*field)[ii].rho_old = NULL;
            (*field)[ii].rho_new = NULL;
            (*field)[ii].T_old = NULL;
            (*field)[ii].T_new = NULL;
            (*field)[ii].e = NULL;
            (*field)[ii].e_old = NULL;
            (*field)[ii].drhodx = NULL;
            (*field)[ii].drhody = NULL;
            (*field)[ii].drhodz = NULL;
            (*field)[ii].dTdx = NULL;
            (*field)[ii].dTdy = NULL;
            (*field)[ii].dTdz = NULL;
            (*field)[ii].dedx = NULL;
            (*field)[ii].dedy = NULL;
            (*field)[ii].dedz = NULL;
            (*field)[ii].mu = NULL;
            (*field)[ii].kappa = NULL;
            (*field)[ii].tau_xx = NULL;
            (*field)[ii].tau_yy = NULL;
            (*field)[ii].tau_zz = NULL;
            (*field)[ii].tau_xy = NULL;
            (*field)[ii].tau_xz = NULL;
            (*field)[ii].tau_yz = NULL;
            (*field)[ii].div_tau_x = NULL;
            (*field)[ii].div_tau_y = NULL;
            (*field)[ii].div_tau_z = NULL;
            (*field)[ii].Q_visc = NULL;
            (*field)[ii].Q_source = NULL;
        }
        
        // ========== INITIALIZE ALL ALLOCATED VARIABLES ==========
        for (int i = 0; i < N; i++) {
            // Incompressible variables
            (*field)[ii].u[i] = 0.0;
            (*field)[ii].v[i] = 0.0;
            (*field)[ii].u_new[i] = 0.0;
            (*field)[ii].v_new[i] = 0.0;
            (*field)[ii].u_old[i] = 0.0;
            (*field)[ii].v_old[i] = 0.0;
            (*field)[ii].p[i] = 0.0;
            (*field)[ii].p_old[i] = 0.0;
            (*field)[ii].pprime[i] = 0.0;
            (*field)[ii].res[i] = 0.0;
            (*field)[ii].source[i] = 0.0;
            (*field)[ii].dpdn[i] = 0.0;
            (*field)[ii].dpdx[i] = 0.0;
            (*field)[ii].dpdy[i] = 0.0;
            (*field)[ii].lapu[i] = 0.0;
            (*field)[ii].lapv[i] = 0.0;
            (*field)[ii].dudx[i] = 0.0;
            (*field)[ii].dudy[i] = 0.0;
            (*field)[ii].dvdx[i] = 0.0;
            (*field)[ii].dvdy[i] = 0.0;

            if (parameters.dimension == 3) {
                (*field)[ii].w[i] = 0.0;
                (*field)[ii].w_new[i] = 0.0;
                (*field)[ii].w_old[i] = 0.0;
                (*field)[ii].dpdz[i] = 0.0;
                (*field)[ii].dudz[i] = 0.0;
                (*field)[ii].dvdz[i] = 0.0;
                (*field)[ii].dwdx[i] = 0.0;
                (*field)[ii].dwdy[i] = 0.0;
                (*field)[ii].dwdz[i] = 0.0;
                (*field)[ii].lapw[i] = 0.0;
            }
            
            // Compressible variables (if allocated)
            if (parameters.compressible_flow) {
                (*field)[ii].rho[i] = parameters.rho_ref;
                (*field)[ii].rho_old[i] = parameters.rho_ref;
                (*field)[ii].rho_new[i] = parameters.rho_ref;
                (*field)[ii].T[i] = parameters.T_ref;  // Initialize to reference temperature
                (*field)[ii].T_old[i] = parameters.T_ref;
                (*field)[ii].T_new[i] = parameters.T_ref;
                (*field)[ii].e[i] = parameters.cv * parameters.T_ref;
                (*field)[ii].e_old[i] = parameters.cv * parameters.T_ref;
                
                (*field)[ii].drhodx[i] = 0.0;
                (*field)[ii].drhody[i] = 0.0;
                (*field)[ii].drhodz[i] = 0.0;
                (*field)[ii].dTdx[i] = 0.0;
                (*field)[ii].dTdy[i] = 0.0;
                (*field)[ii].dTdz[i] = 0.0;
                (*field)[ii].dedx[i] = 0.0;
                (*field)[ii].dedy[i] = 0.0;
                (*field)[ii].dedz[i] = 0.0;
                
                (*field)[ii].mu[i] = parameters.mu_ref;
                (*field)[ii].kappa[i] = parameters.mu_ref * parameters.cp / parameters.Pr;
                
                (*field)[ii].tau_xx[i] = 0.0;
                (*field)[ii].tau_yy[i] = 0.0;
                (*field)[ii].tau_zz[i] = 0.0;
                (*field)[ii].tau_xy[i] = 0.0;
                (*field)[ii].tau_xz[i] = 0.0;
                (*field)[ii].tau_yz[i] = 0.0;
                (*field)[ii].div_tau_x[i] = 0.0;
                (*field)[ii].div_tau_y[i] = 0.0;
                (*field)[ii].div_tau_z[i] = 0.0;
                
                (*field)[ii].Q_visc[i] = 0.0;
                (*field)[ii].Q_source[i] = 0.0;
            }
        }
    }
}

void free_PointStructure(PointStructure* myPointStruct, int num_levels) {
    for (int i = 0; i < num_levels; i++) {
        free(myPointStruct[i].x);
        free(myPointStruct[i].y);
        free(myPointStruct[i].z);
        free(myPointStruct[i].point_index);
        free(myPointStruct[i].x_normal);
        free(myPointStruct[i].y_normal);
        free(myPointStruct[i].z_normal);
        free(myPointStruct[i].boundary_tag);
        free(myPointStruct[i].node_bc);
        free(myPointStruct[i].corner_tag);
        free(myPointStruct[i].pow_x);
        free(myPointStruct[i].pow_y);
        free(myPointStruct[i].pow_z);
        free(myPointStruct[i].cloud_index);
        free(myPointStruct[i].Dx);
        free(myPointStruct[i].Dy);
        if (parameters.dimension == 3)
            free(myPointStruct[i].Dz);
        free(myPointStruct[i].lap);
        free(myPointStruct[i].lap_Poison);
        
        if (i != num_levels-1) {
            free(myPointStruct[i].prolongation_points);
            free(myPointStruct[i].prol_mat);
        }
        if (i != 0) {
            free(myPointStruct[i].restriction_points);
            free(myPointStruct[i].restr_mat);
        }
    }
    free(myPointStruct);
}

void free_field(FieldVariables* field, int num_levels) {
    for (int i = 0; i < num_levels; i++) {
        // Free incompressible variables (always allocated)
        free(field[i].u);
        free(field[i].v);
        free(field[i].u_new);
        free(field[i].v_new);
        free(field[i].u_old);
        free(field[i].v_old);
        free(field[i].p);
        free(field[i].p_old);
        free(field[i].pprime);
        free(field[i].res);
        free(field[i].source);
        free(field[i].dpdn);
        free(field[i].dpdx);
        free(field[i].dpdy);
        free(field[i].dudx);
        free(field[i].dudy);
        free(field[i].dvdx);
        free(field[i].dvdy);
        free(field[i].lapu);
        free(field[i].lapv);

        if (parameters.dimension == 3) {
            free(field[i].w);
            free(field[i].w_new);
            free(field[i].w_old);
            free(field[i].dpdz);
            free(field[i].dudz);
            free(field[i].dvdz);
            free(field[i].dwdx);
            free(field[i].dwdy);
            free(field[i].dwdz);
            free(field[i].lapw);
        }

        // Free compressible variables (if allocated)
        if (parameters.compressible_flow) {
            free(field[i].T);
            free(field[i].rho);
            free(field[i].rho_old);
            free(field[i].rho_new);
            free(field[i].T_old);
            free(field[i].T_new);
            free(field[i].e);
            free(field[i].e_old);
            free(field[i].drhodx);
            free(field[i].drhody);
            free(field[i].drhodz);
            free(field[i].dTdx);
            free(field[i].dTdy);
            free(field[i].dTdz);
            free(field[i].dedx);
            free(field[i].dedy);
            free(field[i].dedz);
            free(field[i].mu);
            free(field[i].kappa);
            free(field[i].tau_xx);
            free(field[i].tau_yy);
            free(field[i].tau_zz);
            free(field[i].tau_xy);
            free(field[i].tau_xz);
            free(field[i].tau_yz);
            free(field[i].div_tau_x);
            free(field[i].div_tau_y);
            free(field[i].div_tau_z);
            free(field[i].Q_visc);
            free(field[i].Q_source);
        }
    }
    free(field);
}
