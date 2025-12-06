// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef OPENACC_FUNCTIONS_H
#define OPENACC_FUNCTIONS_H

#include "structures.h"
#include "functions.h"
#include "mat_lib.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Function declarations
//void copyin_pointstructure_to_gpu(PointStructure* myPointStruct);
//void copyin_field_to_gpu(FieldVariables* field, PointStructure* myPointStruct);
//void copyin_parameters_to_gpu();
//void copypout_pointstructure_from_gpu(PointStructure* myPointStruct);
//void copypout_field_from_gpu(FieldVariables* field, PointStructure* myPointStruct);


// Function definitions
void copyin_parameters_to_gpu(){
    # pragma acc enter data copyin(parameters)
}

void copyin_pointstructure_to_gpu(PointStructure* myPointStruct) {
    #pragma acc enter data copyin(myPointStruct[:parameters.num_levels])

    for (int i = 0; i < parameters.num_levels; i++) {

        #pragma acc enter data copyin( myPointStruct[i].x[:myPointStruct[i].num_nodes],  \
                                       myPointStruct[i].y[:myPointStruct[i].num_nodes],  \
                                       myPointStruct[i].z[:myPointStruct[i].num_nodes],  \
                                       myPointStruct[i].x_normal[:myPointStruct[i].num_nodes], \
                                       myPointStruct[i].y_normal[:myPointStruct[i].num_nodes], \
                                       myPointStruct[i].z_normal[:myPointStruct[i].num_nodes], \
                                       myPointStruct[i].corner_tag[:myPointStruct[i].num_nodes], \
                                       myPointStruct[i].boundary_tag[:myPointStruct[i].num_nodes], \
                                       myPointStruct[i].point_index[:myPointStruct[i].num_nodes], \
                                       myPointStruct[i].cloud_index[:myPointStruct[i].num_nodes * myPointStruct[i].num_cloud_points], \
                                       myPointStruct[i].Dx[:myPointStruct[i].num_nodes * myPointStruct[i].num_cloud_points], \
                                       myPointStruct[i].Dy[:myPointStruct[i].num_nodes * myPointStruct[i].num_cloud_points], \
                                       myPointStruct[i].lap[:myPointStruct[i].num_nodes * myPointStruct[i].num_cloud_points], \
                                       myPointStruct[i].lap_Poison[:myPointStruct[i].num_nodes * myPointStruct[i].num_cloud_points], \
                                       myPointStruct[i].prol_mat[:myPointStruct[i].num_nodes * myPointStruct[i].num_cloud_points], \
                                       myPointStruct[i].restr_mat[:myPointStruct[i].num_nodes * myPointStruct[i].num_cloud_points] )

        if (parameters.dimension == 3) {
            #pragma acc enter data copyin(myPointStruct[i].Dz[:myPointStruct[i].num_nodes * myPointStruct[i].num_cloud_points])
        }

        #pragma acc enter data attach( myPointStruct[i].x, myPointStruct[i].y, myPointStruct[i].z, \
                                       myPointStruct[i].x_normal, myPointStruct[i].y_normal, myPointStruct[i].z_normal, \
                                       myPointStruct[i].corner_tag, myPointStruct[i].boundary_tag, \
                                       myPointStruct[i].point_index, myPointStruct[i].cloud_index, \
                                       myPointStruct[i].Dx, myPointStruct[i].Dy, myPointStruct[i].lap, \
                                       myPointStruct[i].lap_Poison, myPointStruct[i].prol_mat, myPointStruct[i].restr_mat )

        if (parameters.dimension == 3) {
            #pragma acc enter data attach(myPointStruct[i].Dz)
        }
    }
}

void copyin_field_to_gpu(FieldVariables* field, PointStructure* myPointStruct) {
    #pragma acc enter data copyin(field[:parameters.num_levels])

    for (int i = 0; i < parameters.num_levels; i++) {

        // Copy in field arrays
        #pragma acc enter data copyin( field[i].u[:myPointStruct[i].num_nodes],  \
                                       field[i].v[:myPointStruct[i].num_nodes],  \
                                       field[i].w[:myPointStruct[i].num_nodes],  \
                                       field[i].p[:myPointStruct[i].num_nodes],  \
                                       field[i].u_new[:myPointStruct[i].num_nodes], \
                                       field[i].v_new[:myPointStruct[i].num_nodes], \
                                       field[i].w_new[:myPointStruct[i].num_nodes], \
                                       field[i].u_old[:myPointStruct[i].num_nodes], \
                                       field[i].v_old[:myPointStruct[i].num_nodes], \
                                       field[i].w_old[:myPointStruct[i].num_nodes], \
                                       field[i].p_old[:myPointStruct[i].num_nodes], \
                                       field[i].pprime[:myPointStruct[i].num_nodes], \
                                       field[i].res[:myPointStruct[i].num_nodes], \
                                       field[i].source[:myPointStruct[i].num_nodes], \
                                       field[i].T[:myPointStruct[i].num_nodes], \
                                       field[i].dpdn[:myPointStruct[i].num_nodes], \
                                       field[i].dpdx[:myPointStruct[i].num_nodes], \
                                       field[i].dpdy[:myPointStruct[i].num_nodes], \
                                       field[i].dpdz[:myPointStruct[i].num_nodes] )

        // Attach pointer members
        #pragma acc enter data attach( field[i].u, field[i].v, field[i].w, field[i].p,  \
                                       field[i].u_new, field[i].v_new, field[i].w_new,  \
                                       field[i].u_old, field[i].v_old, field[i].w_old,  \
                                       field[i].p_old, field[i].pprime, field[i].res,   \
                                       field[i].source, field[i].T, field[i].dpdn, field[i].dpdx, \
                                       field[i].dpdy, field[i].dpdz )
    }
}

// void copyin_pointstructure_to_gpu(PointStructure* myPointStruct){
//     # pragma acc enter data copyin(myPointStruct[:parameters.num_levels])
//     for (int i = 0; i<parameters.num_levels; i++){
//         # pragma acc enter data copyin(myPointStruct[i].num_nodes)
//         # pragma acc enter data copyin(myPointStruct[i].num_poly_terms)
//         # pragma acc enter data copyin(myPointStruct[i].num_cloud_points)
//         # pragma acc enter data copyin(myPointStruct[i].x[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(myPointStruct[i].y[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(myPointStruct[i].z[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(myPointStruct[i].x_normal[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(myPointStruct[i].y_normal[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(myPointStruct[i].z_normal[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(myPointStruct[i].corner_tag[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(myPointStruct[i].cloud_index[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
//         # pragma acc enter data copyin(myPointStruct[i].pow_x[:myPointStruct[i].num_poly_terms])
//         # pragma acc enter data copyin(myPointStruct[i].pow_y[:myPointStruct[i].num_poly_terms])
//         # pragma acc enter data copyin(myPointStruct[i].pow_z[:myPointStruct[i].num_poly_terms])
//         # pragma acc enter data copyin(myPointStruct[i].point_index[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(myPointStruct[i].boundary_tag[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(myPointStruct[i].prolongation_points[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(myPointStruct[i].restriction_points[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(myPointStruct[i].prol_mat[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
//         # pragma acc enter data copyin(myPointStruct[i].restr_mat[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
//         # pragma acc enter data copyin(myPointStruct[i].Dx[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
//         # pragma acc enter data copyin(myPointStruct[i].Dy[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
//         # pragma acc enter data copyin(myPointStruct[i].lap[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
//         if (parameters.dimension == 3){
//             # pragma acc enter data copyin(myPointStruct[i].Dz[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
//         }
//         # pragma acc enter data copyin(myPointStruct[i].lap_Poison[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
//     }
// }

// void copyin_field_to_gpu(FieldVariables* field, PointStructure* myPointStruct){
//     # pragma acc enter data copyin(field[:parameters.num_levels])
//     for (int i = 0; i<parameters.num_levels; i++){
//         # pragma acc enter data copyin(field[i].u[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].v[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].w[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].pprime[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].u_new[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].v_new[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].w_new[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].u_old[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].v_old[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].w_old[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].p_old[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].T[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].p[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].res[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].source[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].dpdn[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].dpdx[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].dpdy[:myPointStruct[i].num_nodes])
//         # pragma acc enter data copyin(field[i].dpdz[:myPointStruct[i].num_nodes])
//     }
// }

void copypout_pointstructure_from_gpu(PointStructure* myPointStruct){
    # pragma acc exit data copyout(myPointStruct[:parameters.num_levels])
}

void copypout_field_from_gpu(FieldVariables* field, PointStructure* myPointStruct){
    # pragma acc exit data copyout(field[:parameters.num_levels])
}

void copyin_essentials_to_gpu(PointStructure* myPointStruct){
    # pragma acc enter data copyin(myPointStruct[:parameters.num_levels])
    for (int i = 0; i<parameters.num_levels; i++){
        # pragma acc enter data copyin(myPointStruct[i].num_nodes)
        # pragma acc enter data copyin(myPointStruct[i].num_poly_terms)
        # pragma acc enter data copyin(myPointStruct[i].num_cloud_points)
        # pragma acc enter data copyin(myPointStruct[i].x_normal[:myPointStruct[i].num_nodes])
        # pragma acc enter data copyin(myPointStruct[i].y_normal[:myPointStruct[i].num_nodes])
        # pragma acc enter data copyin(myPointStruct[i].z_normal[:myPointStruct[i].num_nodes])
        # pragma acc enter data copyin(myPointStruct[i].corner_tag[:myPointStruct[i].num_nodes])
        # pragma acc enter data copyin(myPointStruct[i].cloud_index[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
        # pragma acc enter data copyin(myPointStruct[i].boundary_tag[:myPointStruct[i].num_nodes])
        # pragma acc enter data copyin(myPointStruct[i].Dx[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
        # pragma acc enter data copyin(myPointStruct[i].Dy[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
        # pragma acc enter data copyin(myPointStruct[i].lap[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
        # pragma acc enter data copyin(myPointStruct[i].lap_Poison[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
        # pragma acc enter data copyin(myPointStruct[i].prol_mat[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
        # pragma acc enter data copyin(myPointStruct[i].restr_mat[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
        # pragma acc enter data copyin(myPointStruct[i].restriction_points[:myPointStruct[i].num_nodes])
        # pragma acc enter data copyin(myPointStruct[i].prolongation_points[:myPointStruct[i].num_nodes])
        if (parameters.dimension == 3){
            # pragma acc enter data copyin(myPointStruct[i].Dz[:myPointStruct[i].num_nodes][:myPointStruct[i].num_cloud_points])
        }
    }
}

#endif
