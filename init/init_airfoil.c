// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef INIT_C
#define INIT_C

#include "../header_files/structures.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void set_up_boundary_tags(PointStructure* myPointStruct, int numlevels);
void initial_conditions(PointStructure* myPointStruct, Compressible_FieldVariables* field, int numlevels);
void boundary_conditions(PointStructure* myPointStruct, Compressible_FieldVariables* field, int numlevels);

void set_up_boundary_tags(PointStructure* myPointStruct, int numlevels)
{
    for (int ii = 0; ii < numlevels; ii++)
    {
        for (int i = 0; i < myPointStruct[ii].num_boundary_nodes; i++)
        {
            if (myPointStruct[ii].x[i] - 20 < 1e-9)
                myPointStruct[ii].bc_outflow[i] = 1;
            else if (myPointStruct[ii].y[i] - 4 < 1e-9)
                myPointStruct[ii].bc_inflow[i] = 1;
            else if (myPointStruct[ii].y[i] + 4 < 1e-9)
                myPointStruct[ii].bc_inflow[i] = 1;
            else if (myPointStruct[ii].x[i] + 2 < 1e-9)
                myPointStruct[ii].bc_inflow[i] = 1;
            else
                myPointStruct[ii].bc_wall[i] = 1;
        }
    }
}

void initial_conditions(PointStructure* myPointStruct, Compressible_FieldVariables* field, int numlevels)
{
    for (int ii = 0; ii < numlevels; ii++)
    {
        for (int i = 0; i < myPointStruct[ii].num_nodes; i++)
        {
            field[ii].M[i] = 0.6;
            field[ii].T[i] = 300.0;
            field[ii].rho[i] = 1.21;
            field[ii].u[i] = field[ii].M[i]*sqrt(parameters.gamma*field[ii].T[i]*parameters.R);
            field[ii].v[i] = 0.0;
            field[ii].w[i] = 0.0;
            field[ii].p[i] = field[ii].rho[i]*parameters.R*field[ii].T[i];
            field[ii].e[i] = field[ii].p[i]/((parameters.gamma-1)*field[ii].rho[i]);
            field[ii].tau_xx[i] = 0;
            field[ii].tau_yy[i] = 0;
            field[ii].tau_zz[i] = 0;
            field[ii].tau_xy[i] = 0;
            field[ii].tau_xz[i] = 0;
            field[ii].tau_yz[i] = 0;
            field[ii].div_u[i] = 0;
            field[ii].q1[i] = 1.21;
            field[ii].q2[i] = field[ii].rho[i]*field[ii].u[i];
            field[ii].q3[i] = field[ii].rho[i]*field[ii].v[i];
            field[ii].q4[i] = field[ii].rho[i]*field[ii].w[i];
            field[ii].q5[i] = field[ii].rho[i]*field[ii].e[i];
            field[ii].q1_step[i] = field[ii].q1[i];
            field[ii].q2_step[i] = field[ii].q2[i];
            field[ii].q3_step[i] = field[ii].q3[i];
            field[ii].q4_step[i] = field[ii].q4[i];
            field[ii].q5_step[i] = field[ii].q5[i];
            field[ii].q1_old[i] = field[ii].q1[i];
            field[ii].q2_old[i] = field[ii].q2[i];
            field[ii].q3_old[i] = field[ii].q3[i];
            field[ii].q4_old[i] = field[ii].q4[i];
            field[ii].q5_old[i] = field[ii].q5[i];
            field[ii].r1[i] = 0;
            field[ii].r2[i] = 0;
            field[ii].r3[i] = 0;
            field[ii].r4[i] = 0;
            field[ii].r5[i] = 0;
            field[ii].u_old[i] = field[ii].u[i];
            field[ii].v_old[i] = field[ii].v[i];
            field[ii].w_old[i] = field[ii].w[i];
            field[ii].res[i] = 0;
            field[ii].source[i] = 0;
            field[ii].dpdn[i] = 0;
            field[ii].mu[i] = 1.8e-5;
            field[ii].nu[i] = field[ii].mu[i]/field[ii].rho[i];
            field[ii].k[i] = 0.0262;
            field[ii].temp1[i] = 0;
            field[ii].temp2[i] = 0;
            field[ii].temp3[i] = 0;
            field[ii].temp4[i] = 0;
            field[ii].temp5[i] = 0;
        }
    }
}

void boundary_conditions(PointStructure* myPointStruct, Compressible_FieldVariables* field, int numlevels)
{
    double x, y, z;
    for (int ii = 0; ii < numlevels; ii++)
    {
        for (int i = 0; i < myPointStruct[ii].num_boundary_nodes; i++)
        {
            x = myPointStruct[ii].x[i]; y = myPointStruct[ii].y[i]; z = myPointStruct[ii].z[i]; 
            if (myPointStruct[ii].bc_inflow[i] == 1)
            {
                field[ii].M[i] = 0.6;
                field[ii].T[i] = 300.0;
                field[ii].rho[i] = 1.21;
                field[ii].u[i] = field[ii].M[i]*sqrt(parameters.gamma*field[ii].T[i]*parameters.R);
                field[ii].v[i] = 0.0;
                field[ii].w[i] = 0.0;
                field[ii].p[i] = field[ii].rho[i]*parameters.R*field[ii].T[i];
                field[ii].e[i] = field[ii].p[i]/((parameters.gamma-1)*field[ii].rho[i]);
            }
            else if (myPointStruct[ii].bc_wall[i] == 1)
            {
                field[ii].u[i] = 0.0;
                field[ii].v[i] = 0.0;
                field[ii].w[i] = 0.0;
                field[ii].T[i] = 273;
                field[ii].rho[i] = 1.0;

                double Ap, sumx, sumy, sumz;
                sumx = 0; sumy = 0; sumz = 0; Ap = 0;
                for (int j = 1; j < myPointStruct[ii].num_cloud_points; j++){
                    sumx += myPointStruct[ii].Dx[i][j]*field[ii].p[myPointStruct[ii].cloud_index[i][j]];
                    sumy += myPointStruct[ii].Dy[i][j]*field[ii].p[myPointStruct[ii].cloud_index[i][j]];
                    if (parameters.dimension == 3){
                        sumz += myPointStruct[ii].Dz[i][j]*field[ii].p[myPointStruct[ii].cloud_index[i][j]];
                    }
                }
                Ap += myPointStruct[ii].Dx[i][0]*myPointStruct[ii].x_normal[i];
                Ap += myPointStruct[ii].Dy[i][0]*myPointStruct[ii].y_normal[i];
                if (parameters.dimension == 3){
                    Ap += myPointStruct[ii].Dz[i][0]*myPointStruct[ii].z_normal[i];
                }
                field[ii].p[i] = (-sumx*myPointStruct[ii].x_normal[i] -sumy*myPointStruct[ii].y_normal[i] -sumz*myPointStruct[ii].z_normal[i])/Ap;
                
                field[ii].e[i] = field[ii].p[i]/((parameters.gamma-1)*field[ii].rho[i]);
            }
            else if (myPointStruct[ii].bc_outflow[i] == 1)
            {
                field[ii].u[i] = 1.0;
                field[ii].v[i] = 0.0;
                field[ii].w[i] = 0.0;
                field[ii].p[i] = 1.0;
                field[ii].T[i] = 1.0;
                field[ii].rho[i] = 1.0;
            }
            else if (myPointStruct[ii].bc_symmetry[i] == 1)
            {
                double Ap, sumx, sumy, sumz;
                sumx = 0; sumy = 0; sumz = 0; Ap = 0;
                for (int j = 1; j < myPointStruct[ii].num_cloud_points; j++){
                    sumx += myPointStruct[ii].Dx[i][j]*field[ii].u[myPointStruct[ii].cloud_index[i][j]];
                    sumy += myPointStruct[ii].Dy[i][j]*field[ii].u[myPointStruct[ii].cloud_index[i][j]];
                    if (parameters.dimension == 3){
                        sumz += myPointStruct[ii].Dz[i][j]*field[ii].u[myPointStruct[ii].cloud_index[i][j]];
                    }
                }
                Ap += myPointStruct[ii].Dx[i][0]*myPointStruct[ii].x_normal[i];
                Ap += myPointStruct[ii].Dy[i][0]*myPointStruct[ii].y_normal[i];
                if (parameters.dimension == 3){
                    Ap += myPointStruct[ii].Dz[i][0]*myPointStruct[ii].z_normal[i];
                }
                field[ii].u[i] = (-sumx*myPointStruct[ii].x_normal[i] -sumy*myPointStruct[ii].y_normal[i] -sumz*myPointStruct[ii].z_normal[i])/Ap;
                
                for (int j = 1; j < myPointStruct[ii].num_cloud_points; j++){
                    sumx += myPointStruct[ii].Dx[i][j]*field[ii].v[myPointStruct[ii].cloud_index[i][j]];
                    sumy += myPointStruct[ii].Dy[i][j]*field[ii].v[myPointStruct[ii].cloud_index[i][j]];
                    if (parameters.dimension == 3){
                        sumz += myPointStruct[ii].Dz[i][j]*field[ii].v[myPointStruct[ii].cloud_index[i][j]];
                    }
                }
                Ap += myPointStruct[ii].Dx[i][0]*myPointStruct[ii].x_normal[i];
                Ap += myPointStruct[ii].Dy[i][0]*myPointStruct[ii].y_normal[i];
                if (parameters.dimension == 3){
                    Ap += myPointStruct[ii].Dz[i][0]*myPointStruct[ii].z_normal[i];
                }
                field[ii].v[i] = (-sumx*myPointStruct[ii].x_normal[i] -sumy*myPointStruct[ii].y_normal[i] -sumz*myPointStruct[ii].z_normal[i])/Ap;

                for (int j = 1; j < myPointStruct[ii].num_cloud_points; j++){
                    sumx += myPointStruct[ii].Dx[i][j]*field[ii].rho[myPointStruct[ii].cloud_index[i][j]];
                    sumy += myPointStruct[ii].Dy[i][j]*field[ii].rho[myPointStruct[ii].cloud_index[i][j]];
                    if (parameters.dimension == 3){
                        sumz += myPointStruct[ii].Dz[i][j]*field[ii].rho[myPointStruct[ii].cloud_index[i][j]];
                    }
                }
                Ap += myPointStruct[ii].Dx[i][0]*myPointStruct[ii].x_normal[i];
                Ap += myPointStruct[ii].Dy[i][0]*myPointStruct[ii].y_normal[i];
                if (parameters.dimension == 3){
                    Ap += myPointStruct[ii].Dz[i][0]*myPointStruct[ii].z_normal[i];
                }
                field[ii].rho[i] = (-sumx*myPointStruct[ii].x_normal[i] -sumy*myPointStruct[ii].y_normal[i] -sumz*myPointStruct[ii].z_normal[i])/Ap;
            
                for (int j = 1; j < myPointStruct[ii].num_cloud_points; j++){
                    sumx += myPointStruct[ii].Dx[i][j]*field[ii].T[myPointStruct[ii].cloud_index[i][j]];
                    sumy += myPointStruct[ii].Dy[i][j]*field[ii].T[myPointStruct[ii].cloud_index[i][j]];
                    if (parameters.dimension == 3){
                        sumz += myPointStruct[ii].Dz[i][j]*field[ii].T[myPointStruct[ii].cloud_index[i][j]];
                    }
                }
                Ap += myPointStruct[ii].Dx[i][0]*myPointStruct[ii].x_normal[i];
                Ap += myPointStruct[ii].Dy[i][0]*myPointStruct[ii].y_normal[i];
                if (parameters.dimension == 3){
                    Ap += myPointStruct[ii].Dz[i][0]*myPointStruct[ii].z_normal[i];
                }
                field[ii].T[i] = (-sumx*myPointStruct[ii].x_normal[i] -sumy*myPointStruct[ii].y_normal[i] -sumz*myPointStruct[ii].z_normal[i])/Ap;
            
            }
        }
    }
}

#endif
