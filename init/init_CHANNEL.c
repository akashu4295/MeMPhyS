// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef INIT_C
#define INIT_C

#include "../header_files/structures.h"
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
                myfieldvariables[ii].u[i] = 1.0;
                myfieldvariables[ii].v[i] = 0.0;
                myfieldvariables[ii].p[i] = 0.0;
                if(parameters.dimension==3)
                    myfieldvariables[ii].w[i] = 0;
            }
            else if (fabs(x-10.0)<1e-9){
                myPointStruct[ii].node_bc[i].type = BC_PRESSURE_OUTLET;
                myfieldvariables[ii].u[i] = 1.0;
                myfieldvariables[ii].v[i] = 0.0;
                myfieldvariables[ii].p[i] = 0.0;
                if(parameters.dimension==3)
                    myfieldvariables[ii].w[i] = 0;
            }
            else if (fabs(y-0.0)<1e-9 || fabs(y-1.0)<1e-9){
                myPointStruct[ii].node_bc[i].type = BC_WALL;
                myfieldvariables[ii].u[i] = 0.0;
                myfieldvariables[ii].v[i] = 0.0;
                myfieldvariables[ii].p[i] = 0.0;
                if(parameters.dimension==3)
                    myfieldvariables[ii].w[i] = 0;
            }
        }
    }
}


#endif
