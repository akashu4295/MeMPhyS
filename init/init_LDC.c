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
    double y;
    for (int ii = 0; ii < numlevels; ii++){
        for (int i = 0; i < myPointStruct[ii].num_nodes; i++){
            y = myPointStruct[ii].y[i]; 
            if (fabs(y-1.0)<1e-9){
                myfieldvariables[ii].u[i] = 1.0;
                myfieldvariables[ii].v[i] = 0.0;
                myfieldvariables[ii].p[i] = 0.0;
                if(parameters.dimension==3)
                    myfieldvariables[ii].w[i] = 0;
            }
        }
    }
}


#endif
