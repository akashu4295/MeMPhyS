// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign
// Functions used to write the output to files

#ifndef WRITE_FUNCTIONS_H
#define WRITE_FUNCTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "structures.h"

//////////////////////////////////////////////////////////////////////
// Function Declarations
//////////////////////////////////////////////////////////////////////
void make_directory(const char* name);
void write_normals(PointStructure* myPointStruct, char* filename);
void write_boundary_tags(PointStructure* myPointStruct, char* filename);
void write_corner_tags(PointStructure* myPointStruct, char* filename);
void write_coordinates(PointStructure* myPointStruct, char* filename);
void write_cloud_index(PointStructure* myPointStruct, char* filename);
void write_prolongation_and_restriction_points(PointStructure* myPointStruct, char* filename);
void write_test_files(double* f, double* fx, double* fy, double* fz, double* lapf, double* fxx, double* fyy, double* fzz, int num_nodes, char* folder1);
void write_grid_filenames(PointStructure* myPointStruct, int num_levels);


//////////////////////////////////////////////////////////////////////
// Function Definitions
//////////////////////////////////////////////////////////////////////

void write_normals(PointStructure* myPointStruct, char* filename)
{
    FILE *file;
    printf("Writing normals to file %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        fprintf(file, "%d %lf %lf %lf\n", i, myPointStruct->x_normal[i], myPointStruct->y_normal[i], myPointStruct->z_normal[i]);
    fclose(file);
}

void write_boundary_tags(PointStructure* myPointStruct, char* filename)
{
    FILE *file;
    printf("Writing boundary tags to file %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        fprintf(file, "%d %d\n", i, myPointStruct->boundary_tag[i]);
    fclose(file);
}

void write_corner_tags(PointStructure* myPointStruct, char* filename)
{
    FILE *file;
    printf("Writing corner tags to file %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        fprintf(file, "%d %d\n", i, myPointStruct->corner_tag[i]);
    fclose(file);
}

void write_coordinates(PointStructure* myPointStruct, char* filename)
{
    FILE *file;
    printf("Writing coordinates to %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < myPointStruct->num_nodes; i++)
        fprintf(file, "%lf %lf %lf\n", myPointStruct->x[i], myPointStruct->y[i], myPointStruct->z[i]);
    fclose(file);
}

void write_cloud_index(PointStructure* myPointStruct, char* filename)
{
    FILE *file;
    printf("Writing cloud index to %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        fprintf(file, "%d ", i);
        for (int j = 0; j < myPointStruct->num_cloud_points; j++)
            fprintf(file, "%d ", myPointStruct->cloud_index[i][j]);
        fprintf(file, "\n");
    }
    fclose(file);
}

void write_prolongation_and_restriction_points(PointStructure* myPointStruct, char* filename)
{
    FILE *file;
    printf("Writing prolongation and restriction points to %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        fprintf(file, "%d %d %d\n", i, myPointStruct->prolongation_points[i], myPointStruct->restriction_points[i]);
    }
    fclose(file);
}

void write_test_files(double* f, double* fx, double* fy, double* fz, double* lapf, double* fxx, double* fyy, double* fzz, int num_nodes, char* folder1)
{
    FILE *file;
    char temp[100];
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"f.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", f[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"fx.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", fx[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"fy.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", fy[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"fz.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", fz[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"lapf.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", lapf[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"fxx.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", fxx[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"fyy.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", fyy[i]);
    }
    fclose(file);
    
    strcpy(temp,folder1);
    file = fopen(strcat(temp,"fzz.csv"), "w");
    for (int i = 0; i < num_nodes; i++) {
        fprintf(file, "%f\n", fzz[i]);
    }
    fclose(file);
    printf("Files written\n");
}

void make_directory(const char* name) {
   #ifdef __linux__
       mkdir(name, 777); 
   #else
       mkdir(name);
   #endif
}

void write_grid_filenames(PointStructure* myPointStruct, int ii)
{   
    FILE *file;
    char filename[50];
    sprintf(filename, "normals_%d.csv", ii);
    write_normals(myPointStruct, filename); // Write normals of all points
    sprintf(filename, "boundary_tags_%d.csv", ii);
    write_boundary_tags(myPointStruct, filename); // Write boundary tags of all points
    sprintf(filename, "corner_tags_%d.csv", ii);
    write_corner_tags(myPointStruct, filename); // Write corner tags of all points
    sprintf(filename, "coordinates_%d.csv", ii);
    write_coordinates(myPointStruct, filename); // Write coordinates of all points
    sprintf(filename, "cloud_index_%d.csv", ii);
    write_cloud_index(myPointStruct, filename); // Write coordinates of all points
    printf("\n\n");
}

#endif