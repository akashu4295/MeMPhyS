// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#include "functions.h"

// Function definitions
void copyin_parameters_to_gpu(){
    # pragma acc enter data copyin(parameters)
}

void copyin_pointstructure_to_gpu(PointStructure *myPointStruct)
{
    int L = parameters.num_levels;

    #pragma acc enter data copyin(myPointStruct[:L])
    for (int l = 0; l < L; l++) {
        int N  = myPointStruct[l].num_nodes;
        int Nc = myPointStruct[l].num_cloud_points;

        #pragma acc enter data copyin( \
            myPointStruct[l].x_normal[:N], \
            myPointStruct[l].y_normal[:N], \
            myPointStruct[l].corner_tag[:N], \
            myPointStruct[l].boundary_tag[:N], \
            myPointStruct[l].node_bc[:N], \
            myPointStruct[l].cloud_index[:N * Nc], \
            myPointStruct[l].Dx[:N * Nc], \
            myPointStruct[l].Dy[:N * Nc], \
            myPointStruct[l].lap[:N * Nc], \
            myPointStruct[l].lap_Poison[:N * Nc], \
            myPointStruct[l].restriction_points[:N], \
            myPointStruct[l].prolongation_points[:N] )

        if (parameters.dimension == 3) {
            #pragma acc enter data copyin(\
                myPointStruct[l].z_normal[:N], \
                myPointStruct[l].Dz[:N * Nc])
        }
        
        if (myPointStruct[l].prol_mat != NULL && l + 1 < L)  {
            size_t prolongation_size = (size_t)N * myPointStruct[l + 1].num_cloud_points;
            #pragma acc enter data copyin(myPointStruct[l].prol_mat[:prolongation_size])
        }
        if (myPointStruct[l].restr_mat != NULL && l > 0)  {
            size_t restriction_size = (size_t)N * myPointStruct[l - 1].num_cloud_points;
            #pragma acc enter data copyin(myPointStruct[l].restr_mat[:restriction_size])
        }
    }
}


void copyin_field_to_gpu(FieldVariables *field, PointStructure *myPointStruct)
{
    int L = parameters.num_levels;

    #pragma acc enter data copyin(field[:L])
    for (int l = 0; l < L; l++) {
        int N = myPointStruct[l].num_nodes;

        #pragma acc enter data copyin( \
            field[l].u[:N], \
            field[l].v[:N], \
            field[l].u_old[:N], \
            field[l].v_old[:N], \
            field[l].u_new[:N], \
            field[l].v_new[:N], \
            field[l].dudx[:N], \
            field[l].dudy[:N], \
            field[l].dvdx[:N], \
            field[l].dvdy[:N], \
            field[l].lapu[:N], \
            field[l].lapv[:N], \
            field[l].p[:N], \
            field[l].p_old[:N], \
            field[l].pprime[:N], \
            field[l].dpdn[:N], \
            field[l].dpdx[:N], \
            field[l].dpdy[:N], \
            field[l].hyper_u[:N], \
            field[l].hyper_v[:N], \
            field[l].lapu_old[:N], \
            field[l].lapv_old[:N], \
            field[l].res[:N], \
            field[l].source[:N] )

        if (parameters.dimension == 3) {
            #pragma acc enter data copyin( \
                field[l].w[:N], \
                field[l].w_old[:N], \
                field[l].w_new[:N], \
                field[l].dwdx[:N], \
                field[l].dwdy[:N], \
                field[l].dwdz[:N], \
                field[l].dvdz[:N], \
                field[l].dudz[:N], \
                field[l].lapw[:N], \
                field[l].dpdz[:N], \
                field[l].hyper_w[:N], \
                field[l].lapw_old[:N] )
        }

        if (parameters.poisson_solver_type == 2) {
            #pragma acc enter data copyin(\
                field[l].color_offsets[0:field[l].num_colors + 1], \
                field[l].color_node_list[:N])
        }
        /* ---------------- compressible variables ---------------- */
        if (parameters.compressible_flow) {

            #pragma acc enter data copyin( \
                field[l].rho[:N], field[l].rho_old[:N], field[l].rho_new[:N], \
                field[l].T_new[:N], field[l].T[:N], field[l].T_old[:N], \
                field[l].e[:N], field[l].e_old[:N], \
                field[l].drhodx[:N], field[l].drhody[:N], field[l].drhodz[:N], \
                field[l].dTdx[:N], field[l].dTdy[:N], field[l].dTdz[:N], \
                field[l].dedx[:N], field[l].dedy[:N], field[l].dedz[:N], \
                field[l].mu[:N], field[l].kappa[:N], \
                field[l].tau_xx[:N], field[l].tau_yy[:N], field[l].tau_zz[:N], \
                field[l].tau_xy[:N], field[l].tau_xz[:N], field[l].tau_yz[:N], \
                field[l].div_tau_x[:N], field[l].div_tau_y[:N], field[l].div_tau_z[:N], \
                field[l].Q_visc[:N], field[l].Q_source[:N] )
        }
    }
}

void free_pointstructure_from_gpu(PointStructure* myPointStruct){
    int L = parameters.num_levels;
    #pragma acc exit data delete(myPointStruct[:L])
}

void free_field_from_gpu(FieldVariables* field, PointStructure* myPointStruct) {
    int L = parameters.num_levels;
    #pragma acc exit data delete(field[:L])
}

void copy_all_data_to_gpu(PointStructure* myPointStruct, FieldVariables* field){
    copyin_parameters_to_gpu();
    copyin_pointstructure_to_gpu(myPointStruct);
    copyin_field_to_gpu(field, myPointStruct);
}

void free_all_data_from_gpu(PointStructure* myPointStruct, FieldVariables* field){
    free_pointstructure_from_gpu(myPointStruct);
    free_field_from_gpu(field, myPointStruct);
}