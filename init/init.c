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
void analytical_solution(PointStructure myPointStruct, double* u_ana, double* v_ana, double* p_ana);

void initial_conditions(PointStructure* myPointStruct, FieldVariables* myfieldvariables, int numlevels)
{
    double x,y,z,r;
    double u_theta = -1, r_o = 2.0, r_i = 1.0, omega_i = -1.0;
    double A = r_i * omega_i * r_o * r_i / ((r_o * r_o) - (r_i * r_i)); 
    for (int ii = 0; ii < numlevels; ii++)
    {
        for (int i = 0; i < myPointStruct[ii].num_nodes; i++)
        {
            x = myPointStruct[ii].x[i]; y = myPointStruct[ii].y[i]; z = myPointStruct[ii].z[i];
            r = sqrt(x*x + y*y + z*z);
            u_theta = A * ((r_o / r) - (r / r_o));
            // myfieldvariables[ii].u[i] = (-u_theta * y / r);
            // myfieldvariables[ii].v[i] = (u_theta * x / r);
            // myfieldvariables[ii].w[i] = 0;
            // myfieldvariables[ii].p[i] = (r * r / (2 * r_o * r_o)) - (2 * log(r)) - (r_o * r_o / (2 * r * r));
            // myfieldvariables[ii].p[i] = myfieldvariables[ii].p[i] * A * A;
            myfieldvariables[ii].u[i] = 0;
            myfieldvariables[ii].v[i] = 0;
            myfieldvariables[ii].w[i] = 0;
            myfieldvariables[ii].p[i] = 0;
        }
    }
}

void boundary_conditions(PointStructure* myPointStruct, FieldVariables* myfieldvariables, int numlevels)
{
    double x, y, z, r;
    double u_theta = -1, r_o = 2, r_i = 1, omega_i = 1;
    for (int ii = 0; ii < numlevels; ii++)
    {
        for (int i = 0; i < myPointStruct[ii].num_nodes; i++)
        {
            x = myPointStruct[ii].x[i]; y = myPointStruct[ii].y[i]; z = myPointStruct[ii].z[i]; 
            // r = sqrt(x*x + y*y + z*z);    
            // double A = r_i * omega_i * r_o * r_i / ((r_o * r_o) - (r_i * r_i));
            // u_theta = A * ((r_o / r) - (r / r_o));
            // if (fabs(r-r_i)<1e-5){
            //     myfieldvariables[ii].u[i] = -y/r;
            //     myfieldvariables[ii].v[i] = x/r;
            //     myfieldvariables[ii].p[i] = A*A*(r*r/(2*r_o*r_o) - 2*log(r) - r_o*r_o/(2*r*r));
            //     if(parameters.dimension==3)
            //         myfieldvariables[ii].w[i] = 0;
            // }
            // else{
            //     myfieldvariables[ii].u[i] = -u_theta * y / r;
            //     myfieldvariables[ii].v[i] = u_theta * x / r;
            //     myfieldvariables[ii].p[i] = A*A*(r*r/(2*r_o*r_o) - 2*log(r) - r_o*r_o/(2*r*r));
            //     if(parameters.dimension==3)
            //         myfieldvariables[ii].w[i] = 0;
            // }
            if (fabs(y-1.0)<1e-5){
                myfieldvariables[ii].u[i] = 1;
                myfieldvariables[ii].v[i] = 0;
                myfieldvariables[ii].p[i] = 0;
                if(parameters.dimension==3)
                    myfieldvariables[ii].w[i] = 0;
            }
            // else{
            //     myfieldvariables[ii].u[i] = 0;
            //     myfieldvariables[ii].v[i] = 0;
            //     myfieldvariables[ii].p[i] = 0;  
            //     if(parameters.dimension==3)
            //         myfieldvariables[ii].w[i] = 0;
            // }
        }
    }
}

void analytical_solution(PointStructure myPointStruct, double* u_ana, double* v_ana, double* p_ana)
{
    double x,y,z,r;
    double u_theta = -1, r_o = 2, r_i = 1, omega_i = 1;
    double A = r_i * omega_i * r_o * r_i / ((r_o * r_o) - (r_i * r_i)); 
    for (int i = 0; i < myPointStruct.num_nodes; i++)
    {
        x = myPointStruct.x[i]; y = myPointStruct.y[i]; z = myPointStruct.z[i];
        r = sqrt(x*x + y*y + z*z);
        u_theta = A * ((r_o / r) - (r / r_o));
        u_ana[i] = (-u_theta * y / r);
        v_ana[i] = (u_theta * x / r);
        p_ana[i] = (r * r / (2 * r_o * r_o)) - (2 * log(r)) - (r_o * r_o / (2 * r * r));
        p_ana[i] = p_ana[i] * A * A;
    }
}


#endif
