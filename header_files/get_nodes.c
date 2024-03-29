/*! gcc -Wall -g -o get_nodes get_nodes.c kdtree.c */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>  
#include <assert.h> 
#include "kdtree.h"
#include "kdtree.c"
#include "mat_lib.h"
#include "parameters.h"

#define DIMENSION_MAX 3 // Maximum dimensionality
double dtemp; int itemp; char temp[50]; // temporary variables used multiple times in the code

static double dist_sq( double *a1, double *a2, int dims ) {
    double dist_sq = 0, diff;
    while( --dims >= 0 ) {
        diff = (a1[dims] - a2[dims]);
        dist_sq += diff*diff;
    }
    return dist_sq;
}

// Structure to represent a point in DIMENSION_MAX-dimensional space
typedef struct Point {
    double coords[DIMENSION_MAX];
    int index; // Index of the point in the original dataset
} Point;

typedef struct mesh
    {   
        char* mesh_filename;
        int* n_count;
        int* dim;
        int* e_count;
        Point* points;
        double* nx;
        double* ny;
        bool* boundary_tag;
        bool* corner_tag;
        int** cloud_index;
    }mesh;

mesh* AllocateMemoryMesh(int nodes) {
    mesh* m1 = (struct mesh*)malloc(sizeof(struct mesh));
    m1->n_count = (int*)malloc(sizeof(int*));
    m1->dim = (int*)malloc(sizeof(int*));
    m1->e_count = (int*)malloc(sizeof(int*));
    m1->mesh_filename = (char*)malloc(50 * sizeof(char*));
    m1->points = (Point *)malloc(nodes * sizeof(Point));
    m1->nx = (double*)malloc(nodes * sizeof(double*));
    m1->ny = (double*)malloc(nodes * sizeof(double*));
    m1->boundary_tag = (bool*)malloc(nodes * sizeof(bool*));
    m1->corner_tag = (bool*)malloc(nodes * sizeof(bool*));
    m1->cloud_index = create_matrix_int(nodes, parameters.n_cloud_points);
    // m1->cloud_index = (int*)malloc(parameters.n_cloud_points*nodes * sizeof(int*));
    *m1->n_count = nodes;
    return m1;
}

void free_mesh(mesh* m1) {
    free(m1->n_count);
    free(m1->dim);
    free(m1->e_count);
    free(m1->points);
    free(m1->nx);
    free(m1->ny);
    free(m1->boundary_tag);
    free(m1->corner_tag);
    free_matrix_int(m1->cloud_index, *m1->n_count);
    free(m1);
}

void calculate_normals(mesh* m1)
{
    double dx, dy, mag, dot_product;
    for (int i = 0; i < *m1->n_count; i++)
    {
        if(m1->boundary_tag[i] == true)
        {
            for (int j = 0; j < parameters.n_cloud_points; j++)
            {
                if (m1->boundary_tag[m1->cloud_index[i][j]] == false)
                {
                    dx = m1->points[m1->cloud_index[i][j]].coords[0] - m1->points[i].coords[0];
                    dy = m1->points[m1->cloud_index[i][j]].coords[1] - m1->points[i].coords[1];
                    dot_product = dx * m1->nx[i] + dy * m1->ny[i];
                    if (dot_product > 0)
                    {
                        m1->nx[i] = -m1->nx[i];
                        m1->ny[i] = -m1->ny[i];
                    }
                    break;
                }
            }
            mag = sqrt(m1->nx[i] * m1->nx[i] + m1->ny[i] * m1->ny[i]);
            m1->nx[i] = m1->nx[i] / mag;
            m1->ny[i] = m1->ny[i] / mag;
        }
    }
}

mesh* readmesh(char* filename)
{
    FILE *file;
    double sum_z = 0;
    file = fopen(filename, "r");
    double x_avg = 0, y_avg = 0, z_avg = 0;
    double x_min = 1e6, y_min = 1e6, z_min = 1e6;
    double x_max = -1e6, y_max = -1e6, z_max = -1e6;
    if (file == NULL)
    {
        printf("Error: Unable to open the file\n");
        exit(1);
    }

    // Read the mesh file nodes
    while (true)
    {
        fscanf(file, "%s ", temp);
        if (strcmp(temp, "$Nodes") == 0)
            break;
    }
    fscanf(file, "%i ", &itemp);
    mesh* m1 = AllocateMemoryMesh(itemp);
    strcpy(m1->mesh_filename, filename);

    for (int i = 0; i < *m1->n_count; i++)
    {
        fscanf(file, "%i ", &itemp); 
        fscanf(file, "%lf ", &dtemp); 
        m1->points[i].coords[0]=dtemp;
        fscanf(file, "%lf ", &dtemp); 
        m1->points[i].coords[1]=dtemp;
        fscanf(file, "%lf ", &dtemp); 
        m1->points[i].coords[2]=dtemp;
        m1->points[i].index=i; 
        m1->boundary_tag[i]=false;
        m1->nx[i]=0;
        m1->ny[i]=0;

        sum_z += abs(m1->points[i].coords[2]);

        if (m1->points[i].coords[0] < x_min)
            x_min = m1->points[i].coords[0];
        if (m1->points[i].coords[1] < y_min)
            y_min = m1->points[i].coords[1];
        if (m1->points[i].coords[2] < z_min)
            z_min = m1->points[i].coords[2];
        if (m1->points[i].coords[0] > x_max)
            x_max = m1->points[i].coords[0];
        if (m1->points[i].coords[1] > y_max)
            y_max = m1->points[i].coords[1];
        if (m1->points[i].coords[2] > z_max)
            z_max = m1->points[i].coords[2];
         
    }
    if (sum_z > 1e-4) //if all z coordinates are zero, then set dimension of the problem to 2
        {parameters.dimension = 2;
        m1->dim = &parameters.dimension;}
    else
        {parameters.dimension = 3;
        m1->dim = &parameters.dimension;}

    // Read the mesh file elements and boundary nodes
    int e_type, tag_int, e_node1, e_node2, e_node3, e_node4;
    double dx,dy;
    
    while (true)
    {
        fscanf(file, "%s ", temp);
        if (strcmp(temp, "$Elements") == 0)
            break;
    }
    fscanf(file, "%i ", &itemp);
    *m1->e_count = itemp;
    
    if (parameters.dimension==2) //if 2D problem, then read only line elements as boundary elements
    {
        for (int ie = 0; ie < *m1->e_count; ie++)
        {
            fscanf(file, "%i ", &itemp);   //element number
            fscanf(file, "%i ", &e_type); //type of the element
        
            if (e_type == 1) //reading line elements
            { 
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); // edge/boundary tag, if a square, elements on one edge will have same tag
                fscanf(file, "%i ", &e_node1);     //node number
                (*m1).boundary_tag[e_node1 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node2);     //node number
                (*m1).boundary_tag[e_node2 - 1] = true; //set boundary tag to true
                dx = (*m1).points[e_node2 - 1].coords[0] - (*m1).points[e_node1 - 1].coords[0];
                dy = (*m1).points[e_node2 - 1].coords[1] - (*m1).points[e_node1 - 1].coords[1];
                (*m1).nx[e_node1 - 1] += -dy; (*m1).ny[e_node1 - 1] += dx;
                (*m1).nx[e_node2 - 1] += -dy; (*m1).ny[e_node2 - 1] += dx;
            }
            else if (e_type == 15) //reading corner nodes
            { 
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); // edge/boundary tag, if a square elements on one edge will have same tag
                fscanf(file, "%i ", &itemp);     //node number
                (*m1).corner_tag[itemp - 1] = true; //set corner tag to true
            }
            else
            {                              //not reading any other kind of elements
                fscanf(file, "%*[^\n]\n"); //skip reading remaining row
            }
        }
    }
    else //if 3D problem, then read only surface elements as boundary elements
    {
        for (int ie = 0; ie < *m1->e_count; ie++)
        {
            fscanf(file, "%i ", &itemp);   //element number
            fscanf(file, "%i ", &e_type); //element type
            if (e_type == 2)
            { //important 3 node triangle CV on the boundary
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); //element tag
                fscanf(file, "%i ", &e_node1);     //node number
                (*m1).boundary_tag[e_node1 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node2);     //node number
                (*m1).boundary_tag[e_node2 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node3);     //node number
                (*m1).boundary_tag[e_node3 - 1] = true; //set boundary tag to true
                dx = (*m1).points[e_node2 - 1].coords[0] - (*m1).points[e_node1 - 1].coords[0];
                dy = (*m1).points[e_node2 - 1].coords[1] - (*m1).points[e_node1 - 1].coords[1];
                (*m1).nx[e_node1 - 1] += -dy; (*m1).ny[e_node1 - 1] += dx;
                (*m1).nx[e_node2 - 1] += -dy; (*m1).ny[e_node2 - 1] += dx;
                dx = (*m1).points[e_node3 - 1].coords[0] - (*m1).points[e_node2 - 1].coords[0];
                dy = (*m1).points[e_node3 - 1].coords[1] - (*m1).points[e_node2 - 1].coords[1];
                (*m1).nx[e_node2 - 1] += -dy; (*m1).ny[e_node2 - 1] += dx;
                (*m1).nx[e_node3 - 1] += -dy; (*m1).ny[e_node3 - 1] += dx;
                dx = (*m1).points[e_node1 - 1].coords[0] - (*m1).points[e_node3 - 1].coords[0];
                dy = (*m1).points[e_node1 - 1].coords[1] - (*m1).points[e_node3 - 1].coords[1];
                (*m1).nx[e_node3 - 1] += -dy; (*m1).ny[e_node3 - 1] += dx;
                (*m1).nx[e_node1 - 1] += -dy; (*m1).ny[e_node1 - 1] += dx;
            }
            else if (e_type == 3)
            { //important 4 node quad element on the boundary
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); //element tag
                fscanf(file, "%i ", &e_node1);     //node number
                (*m1).boundary_tag[e_node1 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node2);     //node number
                (*m1).boundary_tag[e_node2 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node3);     //node number
                (*m1).boundary_tag[e_node3 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node4);     //node number
                (*m1).boundary_tag[e_node4 - 1] = true; //set boundary tag to true
                dx = (*m1).points[e_node2 - 1].coords[0] - (*m1).points[e_node1 - 1].coords[0];
                dy = (*m1).points[e_node2 - 1].coords[1] - (*m1).points[e_node1 - 1].coords[1];
                (*m1).nx[e_node1 - 1] += -dy; (*m1).ny[e_node1 - 1] += dx;
                (*m1).nx[e_node2 - 1] += -dy; (*m1).ny[e_node2 - 1] += dx;
                dx = (*m1).points[e_node3 - 1].coords[0] - (*m1).points[e_node2 - 1].coords[0];
                dy = (*m1).points[e_node3 - 1].coords[1] - (*m1).points[e_node2 - 1].coords[1];
                (*m1).nx[e_node2 - 1] += -dy; (*m1).ny[e_node2 - 1] += dx;
                (*m1).nx[e_node3 - 1] += -dy; (*m1).ny[e_node3 - 1] += dx;
                dx = (*m1).points[e_node4 - 1].coords[0] - (*m1).points[e_node3 - 1].coords[0];
                dy = (*m1).points[e_node4 - 1].coords[1] - (*m1).points[e_node3 - 1].coords[1];
                (*m1).nx[e_node3 - 1] += -dy; (*m1).ny[e_node3 - 1] += dx;
                (*m1).nx[e_node4 - 1] += -dy; (*m1).ny[e_node4 - 1] += dx;
                dx = (*m1).points[e_node1 - 1].coords[0] - (*m1).points[e_node4 - 1].coords[0];
                dy = (*m1).points[e_node1 - 1].coords[1] - (*m1).points[e_node4 - 1].coords[1];
                (*m1).nx[e_node4 - 1] += -dy; (*m1).ny[e_node4 - 1] += dx;
                (*m1).nx[e_node1 - 1] += -dy; (*m1).ny[e_node1 - 1] += dx;
            }
            else if (e_type == 1)
            { //2 node line: log these 2 vertex numbers to delete them later
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); //element tag
                fscanf(file, "%i ", &itemp);   //first point
                (*m1).corner_tag[itemp - 1] = true;
                fscanf(file, "%i ", &itemp); //second point
                (*m1).corner_tag[itemp - 1] = true;
            }
            else
            {                              //not reading any other kind of element
                fscanf(file, "%*[^\n]\n"); //skip reading remaining row
            }
        }
    }
    fclose(file);
    printf("No of nodes = %d \nNo of elements = %d \n", *m1->n_count, *m1->e_count);

    int K = parameters.n_cloud_points; // Number of nearest neighbors to find
    
    printf("\nBuilding KDTree!\n");
    void* ptree;
    ptree = kd_create( 3 ); // create a k-d tree for 3-dimensional points 
    for(int i=0; i<*m1->n_count; i++ ) {
        assert(kd_insert3(ptree, m1->points[i].coords[0], m1->points[i].coords[1], m1->points[i].coords[2], &m1->points[i].index) == 0);
    }
    printf("Done!\n");
    printf("\nSearching for neighbours!\n");
    // Search for neighbours
    double pos[3], dist;
    double pt[3] = { 0.0, 0.0, 0.0 };
    double radius = 0.3; // Search radius
    double *distance = (double*)malloc(K * sizeof(double*));
    int *ind = (int*)malloc(K * sizeof(int*));

    filename = "cloud_index.txt";
    file = fopen(filename, "w");
    for(int i=0; i<*m1->n_count; i++ ) {
        struct kdres *presults;
        int* pch; // Pointer to the index of the nearest neighbor
            //initialise the distance and index arrays
        for (int i =0; i<K ; i++) {
            distance[i] = 100.0;
            ind[i] = 0;
        }
        pt[0] = m1->points[i].coords[0];
        pt[1] = m1->points[i].coords[1];
        pt[2] = m1->points[i].coords[2];
        // printf("(%g, %g, %g) -> ", pt[0], pt[1], pt[2]); 
        presults = kd_nearest_range( ptree, pt, radius ); // find points closest to the point pt and within distance radius

        while (!kd_res_end( presults ))
        {
            pch = (int*)kd_res_item( presults, pos );
            dist = sqrt( dist_sq( pt, pos, 3 ) );
            // printf("node at (%.3f, %.3f, %.3f) is %.3f away and has data=%d\n", pos[0], pos[1], pos[2], dist, *pch);
            // printf("(%g, %g, %g)\n", pt[0], pt[1], pt[2]); 
            for (int j = 0; j < K; ++j) {
                if (dist < distance[j]) {
                    for (int k = K-1; k > j; k--) {
                        distance[k] = distance[k-1];
                        ind[k] = ind[k-1];
                    }
                    distance[j] = dist;
                    ind[j] = *pch;
                    break;
                }
            }
            kd_res_next( presults );
        }

        for (int j = 0; j < K; ++j) {
            m1->cloud_index[i][j] = ind[j];
            // printf("%d,", m1->cloud_index[i*K+j]);
            fprintf(file, "%d,", m1->cloud_index[i][j]);
        }      

        fprintf(file, "\n");
        kd_res_free( presults );
    }
    fclose(file);
    printf("Done! Neighbour indices written to file %s!\n", filename);
    kd_free( ptree );
    free(distance);
    free(ind);
    
    calculate_normals(m1);
    return m1;
}
