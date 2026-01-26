// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef INIT_C
#define INIT_C

#include "header_files/structures.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void initial_conditions(PointStructure* myPointStruct, FieldVariables* myfieldvariables, int numlevels);
void boundary_conditions(PointStructure* myPointStruct, FieldVariables* myfieldvariables, int numlevels);

void initial_conditions(PointStructure* myPointStruct, FieldVariables* myfieldvariables, int numlevels)
{
    for (int ii = 0; ii < numlevels; ii++){
        for (int i = 0; i < myPointStruct[ii].num_nodes; i++){
            myPointStruct[ii].node_bc[i].type = BC_INTERIOR;
            myPointStruct[ii].node_bc[i].u = 0;
            myPointStruct[ii].node_bc[i].v = 0;
            myPointStruct[ii].node_bc[i].w = 0;
            myPointStruct[ii].node_bc[i].p = 0;
            myfieldvariables[ii].u[i] = 0;
            myfieldvariables[ii].v[i] = 0;
            myfieldvariables[ii].w[i] = 0;
            myfieldvariables[ii].p[i] = 0;
            myfieldvariables[ii].p_old[i] = 0;
        }
    }
}

void boundary_conditions(PointStructure* myPointStruct, FieldVariables* myfieldvariables, int numlevels)
{
    double x,y;
    for (int ii = 0; ii < numlevels; ii++){
        for (int i = 0; i < myPointStruct[ii].num_nodes; i++){
            x = myPointStruct[ii].x[i];
            y = myPointStruct[ii].y[i]; 
            if (fabs(x-0.0)<1e-9){
                myPointStruct[ii].node_bc[i].type = BC_VELOCITY_INLET;
                myPointStruct[ii].node_bc[i].u = 1.0;
                myPointStruct[ii].node_bc[i].v = 0.0;
                myPointStruct[ii].node_bc[i].p = 0.0;
                if(parameters.dimension==3)
                    myPointStruct[ii].node_bc[i].w = 0;
            }
            else if (fabs(x-10.0)<1e-9){
                myPointStruct[ii].node_bc[i].type = BC_PRESSURE_OUTLET;
                myPointStruct[ii].node_bc[i].u = 0.0;
                myPointStruct[ii].node_bc[i].v = 0.0;
                myPointStruct[ii].node_bc[i].p = 0.0;
                if(parameters.dimension==3)
                    myPointStruct[ii].node_bc[i].w = 0;
            }
            else if (fabs(y-0.0)<1e-9 || fabs(y-1.0)<1e-9){
                myPointStruct[ii].node_bc[i].type = BC_WALL;
                myPointStruct[ii].node_bc[i].u = 0.0;
                myPointStruct[ii].node_bc[i].v = 0.0;
                myPointStruct[ii].node_bc[i].p = 0.0;
                if(parameters.dimension==3)
                    myPointStruct[ii].node_bc[i].w = 0;
            }
        }
    }
}



////// Taylor Couette flow
// void boundary_conditions(PointStructure* myPointStruct, FieldVariables* myfieldvariables, int numlevels){
//     double x, y, r;
//     double r_i = 1.0; // Inner radius
//     double r_o = 2.0; // Outer radius
//     for (int ii = 0; ii < numlevels; ii++){
//         for (int i = 0; i < myPointStruct[ii].num_nodes; i++){
//             x = myPointStruct[ii].x[i]; y = myPointStruct[ii].y[i];
//             r = sqrt(x*x + y*y);
//             if (fabs(r-r_i)<1e-9){
//                 myPointStruct[ii].node_bc[i].type = BC_VELOCITY_INLET;
//                 myfieldvariables[ii].u[i] = -y/r;
//                 myfieldvariables[ii].v[i] = x/r;
//                 myfieldvariables[ii].p[i] = 0.0;
//                 if(parameters.dimension==3){
//                     myfieldvariables[ii].w[i] = 0.0;
//                 }
//             }
//             else if (fabs(r-r_o)<1e-9){
//                 myPointStruct[ii].node_bc[i].type = BC_WALL;
//                 myfieldvariables[ii].u[i] = 0.0;
//                 myfieldvariables[ii].v[i] = 0.0;
//                 myfieldvariables[ii].p[i] = 0.0;
//                 if(parameters.dimension==3){
//                     myfieldvariables[ii].w[i] = 0.0;
//                 }   
//             }
//         }
//     }
// }



#endif