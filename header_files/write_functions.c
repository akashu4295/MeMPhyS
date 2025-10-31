// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign
// Functions used to write the output to files

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functions_supplementary.h"

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

    int n = myPointStruct->num_cloud_points;
    for (int i = 0; i < myPointStruct->num_nodes; i++) {
        fprintf(file, "%d ", i);
        for (int j = 0; j < n; j++)
            fprintf(file, "%d ", myPointStruct->cloud_index[i*n +j]);
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

void write_processed_grid_data(PointStructure* myPointStruct, int ii)
{   
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
    sprintf(filename, "prolongation_and_restriction_%d.csv", ii);
    write_prolongation_and_restriction_points(myPointStruct, filename);
    printf("\n\n");
}

void write_vtk(PointStructure* myPointStruct, FieldVariables* field)
{
    FILE *fp = fopen("solution.vtk", "w");
    int n_nodes = myPointStruct[0].num_nodes;
    int n_elems = myPointStruct[0].num_elem;
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Meshless solver output\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // --- Write points ---
    fprintf(fp, "POINTS %d float\n", n_nodes);
    for (int i = 0; i < n_nodes; i++)
        fprintf(fp, "%f %f %f\n", myPointStruct[0].x[i], myPointStruct[0].y[i], myPointStruct[0].z[i]);

    // --- Write cells ---
    int nodes_per_elem = 3; // for triangle
    fprintf(fp, "CELLS %d %d\n", n_elems, n_elems*(nodes_per_elem+1));
    for (int i = 0; i < n_elems; i++) {
        fprintf(fp, "%d", nodes_per_elem);
        for (int j = 0; j < nodes_per_elem; j++) {
            fprintf(fp, " %d", elem_conn[i*nodes_per_elem + j] - 1); // convert 1-based -> 0-based
        }
        fprintf(fp, "\n");
    }

    // --- Write cell types ---
    fprintf(fp, "CELL_TYPES %d\n", n_elems);
    for (int i = 0; i < n_elems; i++)
        fprintf(fp, "5\n");  // VTK_TRIANGLE

    // --- Write point data (e.g., pressure) ---
    fprintf(fp, "POINT_DATA %d\n", n_nodes);
    fprintf(fp, "SCALARS pressure float 1\nLOOKUP_TABLE default\n");
    for (int i = 0; i < n_nodes; i++)
        fprintf(fp, "%f\n", field[0].p[i]);

    fclose(fp);

}