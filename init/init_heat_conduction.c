// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef INIT_C
#define INIT_C

#include "../header_files/structures_vectors.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void initial_conditions(PointStructure* myPointStruct, FieldVariables* field, int numlevels, int k);
void boundary_conditions(PointStructure* myPointStruct, FieldVariables* field, int numlevels, int k);
void set_boundary_rhs(PointStructure* myPointStruct, FieldVariables* field, int numlevels, int k);


void boundary_conditions(PointStructure* myPointStruct, FieldVariables* field, int numlevels, int k){
    double x;
    for (int ilev = 0; ilev < numlevels; ilev++){
        for (int i = 0; i < myPointStruct[ilev].num_nodes; i++){
            if (myPointStruct[ilev].boundary_tag[i] == true){
                x = myPointStruct[ilev].x[i];
                if (x == 0){
                    field[ilev].T[i] = 0.0;
                    field[ilev].p[i] = 0.0;
                }
            }
        }
    }
}

void initial_conditions(PointStructure* myPointStruct, FieldVariables* field, int numlevels, int k){
    for (int ilev = 0; ilev < numlevels; ilev++){
        for (int i = 0; i < myPointStruct[ilev].num_nodes; i++){
            field[ilev].res[i] = 0.0;
            if (ilev==0){
                field[ilev].source[i] = -8 * 3.14 * 3.14 * k * k * 
                     sin(2 * 3.14 * k * myPointStruct->x[i]) * 
                     sin(2 * 3.14 * k * myPointStruct->y[i]);
            }
            else{
                field[ilev].source[i] = 0.0;
            }
            field[ilev].p[i] = 0.0;
            field[ilev].T[i] = 0.0;
        }
    }
}

void set_boundary_rhs(PointStructure* myPointStruct, FieldVariables* field, int numlevels, int k){
    int n =myPointStruct->num_cloud_points;
    for (int i = 0; i < myPointStruct->num_nodes; i++){
        if (myPointStruct->boundary_tag[i] == true){
            if (myPointStruct->x[i] == 0){
                field[0].source[i] = 1.0;
                int k = i*n;
                for (int j = 0; j< n; j++){
                    myPointStruct[0].lap_Poison[k] = 0.0;
                    k +=1;
                }
                myPointStruct[0].lap_Poison[i*n] = 1.0;
            }
            else
                field[0].source[i] = 2*3.14*k*(sin(2*3.14*k*myPointStruct->x[i])*cos(2*3.14*k*myPointStruct->y[i])*myPointStruct->y_normal[i] +
                                                cos(2*3.14*k*myPointStruct->x[i])*sin(2*3.14*k*myPointStruct->y[i])*myPointStruct->x_normal[i]);
        }
    }
}

#endif
