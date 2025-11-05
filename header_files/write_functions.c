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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Structure to hold field data (velocities and pressure)
typedef struct {
    double *u;  // x-velocity
    double *v;  // y-velocity
    double *w;  // z-velocity
    double *p;  // pressure
} Field;

// Function to convert Gmsh mesh to VTK with field data
int write_vtk(char *gmsh_filename, FieldVariables *field, PointStructure* myPS) 
{
    FILE *fp_in, *fp_out;
    char line[256];
    char vtk_filename[256] = "Solution.vtk";
    int num_nodes = 0, num_elements = 0;
    int i, node_id;
    double x, y, z;
    
    // Open Gmsh file
    fp_in = fopen(gmsh_filename, "r");
    if (!fp_in) {
        fprintf(stderr, "Error: Cannot open file %s\n", gmsh_filename);
        return -1;
    }

    // First pass: count nodes and elements
    while (fgets(line, sizeof(line), fp_in)) {
        if (strstr(line, "$Nodes")) {
            if (fscanf(fp_in, "%d", &num_nodes) != 1) {
                fprintf(stderr, "Error reading number of nodes\n");
                fclose(fp_in);
                return -1;
            }
            break;
        }
    }
    
    // Allocate memory for node coordinates
    double *nodes_x = (double*)malloc(num_nodes * sizeof(double));
    double *nodes_y = (double*)malloc(num_nodes * sizeof(double));
    double *nodes_z = (double*)malloc(num_nodes * sizeof(double));
    
    if (!nodes_x || !nodes_y || !nodes_z) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        fclose(fp_in);
        return -1;
    }
    
    // Read node coordinates
    for (i = 0; i < num_nodes; i++) {
        if (fscanf(fp_in, "%d %lf %lf %lf", &node_id, &x, &y, &z) != 4) {
            fprintf(stderr, "Error reading node %d\n", i);
            free(nodes_x); free(nodes_y); free(nodes_z);
            fclose(fp_in);
            return -1;
        }
        nodes_x[node_id - 1] = x;
        nodes_y[node_id - 1] = y;
        nodes_z[node_id - 1] = z;
    }
    
    // Find $Elements section
    while (fgets(line, sizeof(line), fp_in)) {
        if (strstr(line, "$Elements")) {
            if (fscanf(fp_in, "%d", &num_elements) != 1) {
                fprintf(stderr, "Error reading number of elements\n");
                free(nodes_x); free(nodes_y); free(nodes_z);
                fclose(fp_in);
                return -1;
            }
            break;
        }
    }
    
    // Count triangular elements (type 2)
    int num_triangles = 0;
    int *triangle_conn = (int*)malloc(num_elements * 3 * sizeof(int));
    long elements_pos = ftell(fp_in);
    
    for (i = 0; i < num_elements; i++) {
        int elem_id, elem_type, num_tags, tag1, n1, n2, n3;
        if (fscanf(fp_in, "%d %d %d", &elem_id, &elem_type, &num_tags) != 3) {
            fprintf(stderr, "Error reading element %d\n", i);
            free(nodes_x); free(nodes_y); free(nodes_z); free(triangle_conn);
            fclose(fp_in);
            return -1;
        }
        
        // Skip tags
        for (int j = 0; j < num_tags; j++) {
            fscanf(fp_in, "%d", &tag1);
        }
        
        if (elem_type == 2) {  // Triangle
            if (fscanf(fp_in, "%d %d %d", &n1, &n2, &n3) != 3) {
                fprintf(stderr, "Error reading triangle nodes\n");
                free(nodes_x); free(nodes_y); free(nodes_z); free(triangle_conn);
                fclose(fp_in);
                return -1;
            }
            triangle_conn[num_triangles * 3 + 0] = n1 - 1;
            triangle_conn[num_triangles * 3 + 1] = n2 - 1;
            triangle_conn[num_triangles * 3 + 2] = n3 - 1;
            num_triangles++;
        } else {
            // Skip other element types
            fgets(line, sizeof(line), fp_in);
        }
    }
    
    fclose(fp_in);
    
    int num_corner = myPS->num_corners;
    // Write VTK file
    fp_out = fopen(vtk_filename, "w");
    if (!fp_out) {
        fprintf(stderr, "Error: Cannot create file %s\n", vtk_filename);
        free(nodes_x); free(nodes_y); free(nodes_z); free(triangle_conn);
        return -1;
    }
    
    // Write VTK header
    fprintf(fp_out, "# vtk DataFile Version 3.0\n");
    fprintf(fp_out, "Gmsh mesh with velocity and pressure fields\n");
    fprintf(fp_out, "ASCII\n");
    fprintf(fp_out, "DATASET UNSTRUCTURED_GRID\n");
    
    // Write points
    fprintf(fp_out, "POINTS %d double\n", num_nodes);
    for (i = 0; i < num_nodes; i++) {
        fprintf(fp_out, "%.16e %.16e %.16e\n", nodes_x[i], nodes_y[i], nodes_z[i]);
    }
    
    // Write cells (triangles)
    fprintf(fp_out, "\nCELLS %d %d\n", num_triangles, num_triangles * 4);
    for (i = 0; i < num_triangles; i++) {
        fprintf(fp_out, "3 %d %d %d\n", 
                triangle_conn[i * 3 + 0],
                triangle_conn[i * 3 + 1],
                triangle_conn[i * 3 + 2]);
    }
    
    // Write cell types (5 = triangle)
    fprintf(fp_out, "\nCELL_TYPES %d\n", num_triangles);
    for (i = 0; i < num_triangles; i++) {
        fprintf(fp_out, "5\n");
    }
    
    // Write point data
    fprintf(fp_out, "\nPOINT_DATA %d\n", num_nodes);
    
    // Write velocity vector field
    fprintf(fp_out, "VECTORS velocity double\n");
    for (i = 0; i < num_corner; i++) {
        fprintf(fp_out, "%.16e %.16e %.16e\n", 0.0, 0.0, 0.0);
    }
    for (i = num_corner; i < num_nodes; i++) {
        int k = myPS[0].rcm_order[i-num_corner];
        fprintf(fp_out, "%.16e %.16e %.16e\n", field->u[k], field->v[k], field->w[k]);
    }
    
    // Write pressure scalar field
    fprintf(fp_out, "\nSCALARS pressure double 1\n");
    fprintf(fp_out, "LOOKUP_TABLE default\n");
    for (i = 0; i < num_corner; i++) {
        fprintf(fp_out, "%.16e\n", 0.0);
    }
    for (i = num_corner; i < num_nodes; i++) {
        int k = myPS[0].rcm_order[i-num_corner];
        fprintf(fp_out, "%.16e\n", field->p[k]);
    }
    
    fclose(fp_out);
    
    // Cleanup
    free(nodes_x);
    free(nodes_y);
    free(nodes_z);
    free(triangle_conn);
    
    printf("Successfully converted %s to %s\n", gmsh_filename, vtk_filename);
    printf("Nodes: %d, Triangles: %d\n", num_nodes, num_triangles);
    
    return 0;
}
