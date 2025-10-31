#ifndef HEAT_CONDUCTION_H
#define HEAT_CONDUCTION_H


#include <time.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "structures.h"
#include "test_functions.h"
#include "functions.h"
#include "mat_lib.h"
#include "kdtree_functions.h"
#include "rbf.h"


/////////////////////////////////////////////////////////////////////////////
// Test Heat conduction V-Cycle Solver Modules
/////////////////////////////////////////////////////////////////////////////

void relaxation(PointStructure* mypointstruct, FieldVariables* field){
    double* zeros=create_vector(mypointstruct->num_nodes);
    # pragma acc loop
    int n = mypointstruct->num_cloud_points;
    for (int iter = 0; iter < parameters.num_relax; iter++){
        for (int i = 0; i < mypointstruct->num_nodes; i++){
            int k = i*n;
            double sum = 0.0;
            for (int j = 1; j < mypointstruct->num_cloud_points; j++){
                k += 1;
                sum += mypointstruct->lap_Poison[k]*field->p[mypointstruct->cloud_index[k]];
            }
            k = i*n;
            field->p[i] = parameters.omega*((field->source[i]-sum)/mypointstruct->lap_Poison[k]) + (1-parameters.omega)*field->p[i]; 
            //if (mypointstruct->boundary_tag[i] && mypointstruct->x[i]-1.0<1e-9){
               // field->p[i] = 0.0;
            }
        
        /*# pragma acc parallel loop present(field, mypointstruct)
        for (int i = 0; i< mypointstruct->num_nodes; i++){
          if (mypointstruct->boundary_tag[i]==false){ 
            field->res[i]=field->source[i];
            for (int j = 0; j < mypointstruct->num_cloud_points; j++){
                field->res[i] = field->res[i] - mypointstruct->lap_Poison[i][j]*field->p[mypointstruct->cloud_index[i][j]];
            }
          }
            else
            field->res[i] = 0;
        }*/
        
       //printf("Iter, Pressure Residual: %d %e\n", iter,l2_norm_gen(mypointstruct, field->res, zeros, mypointstruct->num_nodes)); 
     
    }
    printf("Pressure Residual: %e\n", l2_norm_gen(mypointstruct, field->res, zeros, mypointstruct->num_nodes)); 
    free(zeros);
}

void restrict_residuals(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, 
                                                    FieldVariables* field_f, FieldVariables* field_c){
    double* zeros=create_vector(mypointStruct_c->num_nodes);
    int n = mypointStruct_f->num_cloud_points;
    # pragma acc parallel loop present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_c->num_nodes; i++){
        double results = 0.0;
        if (mypointStruct_c->boundary_tag[i]==false){   
            int k = i*n;       
            int i_restr_n = mypointStruct_c->restriction_points[i] *n ;
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < mypointStruct_f->num_cloud_points; j++){
                results += mypointStruct_c->restr_mat[k] * field_f->res[mypointStruct_f->cloud_index[i_restr_n +j]];
                k += 1;
            }           
        }
        
        field_c->source[i] = results;
        //printf("i, Restricted residual: %d,%e\n",i, l2_norm_gen(mypointStruct_c, field_c->source, zeros, mypointStruct_c->num_nodes));
    }
    //printf("Restricted residual: %e\n",l2_norm_gen(mypointStruct_c, field_c->source, zeros, mypointStruct_c->num_nodes));
    free (zeros); 
}

void prolongate_corrections(PointStructure* mypointStruct_f, PointStructure* mypointStruct_c, FieldVariables* field_f, FieldVariables* field_c){
    int n = mypointStruct_c->num_cloud_points;                                       
    # pragma acc parallel loop present(field_f, field_c, mypointStruct_f, mypointStruct_c)
    for (int i = 0; i < mypointStruct_f->num_nodes; i++){
        if (mypointStruct_f->boundary_tag[i]==false){
            int i_prol_n = mypointStruct_f->prolongation_points[i]*n;;
            double results = 0.0;
            int k = i*n;
            # pragma acc loop reduction(+:results)
            for (int j = 0; j < n; j++){
                results += mypointStruct_f->prol_mat[k] * field_c->p[mypointStruct_c->cloud_index[i_prol_n +j]];
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

void calculate_residuals(PointStructure* mypointStruct, FieldVariables* field){
    int n = mypointStruct->num_cloud_points;
    for (int i = 0; i < mypointStruct->num_nodes; i++){
        field->res[i] = 0.0;
        if (!mypointStruct->boundary_tag[i]){
            int k = i*n;
            double sum = 0;
            for (int j = 0; j < mypointStruct->num_cloud_points; j++){
                sum += mypointStruct->lap_Poison[k]*field->p[mypointStruct->cloud_index[k]];
                k += 1;
            }
            field->res[i] = field->source[i] - sum;
        }
    }
}



#endif
