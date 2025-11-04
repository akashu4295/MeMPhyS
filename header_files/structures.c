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
    myPointStruct->num_nodes = nodes;
    myPointStruct->num_elem = 0;
    myPointStruct->rcm_order = (int*)malloc(nodes*sizeof(int));
}

void AllocateMemoryFieldVariables(FieldVariables** field, PointStructure* myPointStruct, int num_levels) {
    *field = (FieldVariables*) malloc(num_levels * sizeof(FieldVariables));
    for (int ii = 0; ii < num_levels; ii++) {
        (*field)[ii].u = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].v = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].w = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].u_new = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].v_new = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].w_new = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].pprime = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].p = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].u_old = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].v_old = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].w_old = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].p_old = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].res = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].source = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].dpdn = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].dpdx = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].dpdy = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].dpdz = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        (*field)[ii].T = (double*) malloc(myPointStruct[ii].num_nodes * sizeof(double));
        for (int i = 0; i < myPointStruct[ii].num_nodes; i++){
            (*field)[ii].u[i] = 0.0;
            (*field)[ii].v[i] = 0.0;
            (*field)[ii].w[i] = 0.0;
            (*field)[ii].u_new[i] = 0.0;
            (*field)[ii].v_new[i] = 0.0;
            (*field)[ii].w_new[i] = 0.0;
            (*field)[ii].u_old[i] = 0.0;
            (*field)[ii].v_old[i] = 0.0;
            (*field)[ii].w_old[i] = 0.0;
            (*field)[ii].pprime[i] = 0.0;
            (*field)[ii].p_old[i] = 0.0;
            (*field)[ii].p[i] = 0.0;
            (*field)[ii].res[i] = 0.0;
            (*field)[ii].source[i] = 0.0;
            (*field)[ii].dpdn[i] = 0.0;
            (*field)[ii].dpdx[i] = 0.0;
            (*field)[ii].dpdy[i] = 0.0;
            (*field)[ii].dpdz[i] = 0.0;
            (*field)[ii].T[i] = 0.0;
        }
    }
}

void free_PointStructure(PointStructure* myPointStruct, int num_levels) {
    for (int i = 0; i<num_levels; i++){
        free(myPointStruct[i].x);
        free(myPointStruct[i].y);
        free(myPointStruct[i].z);
        free(myPointStruct[i].point_index);
        free(myPointStruct[i].x_normal);
        free(myPointStruct[i].y_normal);
        free(myPointStruct[i].z_normal);
        free(myPointStruct[i].boundary_tag);
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
        if (i == num_levels-1){}
            //free(myPointStruct[i].prolongation_points);
        else{
            free(myPointStruct[i].prolongation_points);
            free(myPointStruct[i].prol_mat);
        }
        if (i == 0){}
            // free(myPointStruct[i].restriction_points);
        else{
            free(myPointStruct[i].restriction_points);
            free(myPointStruct[i].restr_mat);
        }
    }
    free(myPointStruct);
}

void free_field(FieldVariables* field, int num_levels) {
    for (int i = 0; i < num_levels; i++){
        free(field[i].u);
        free(field[i].v);
        free(field[i].w);
        free(field[i].pprime);
        free(field[i].T);
        free(field[i].res);
        free(field[i].source);
        free(field[i].dpdn);
        free(field[i].u_new);
        free(field[i].v_new);
        free(field[i].w_new);
    }
    free(field);
}
