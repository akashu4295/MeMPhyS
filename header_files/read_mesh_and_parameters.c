// Author :  Akash Unnikrishnan
// Affiliation : Indian Institute of Technology Gandhinagar
#include <time.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structures.h"
#include <assert.h> 
#include "kdtree.h"
#include "kdtree.c"
#include "mat_lib.h"

// Function declarations
int* find_neighbours(double* p, void* ptree, double radius, struct parameters *parameters);
int* find_neighbours_for_mesh(Point* p, void* ptree, double radius, struct parameters *parameters, mesh* m1);
void* create_kdtree(mesh* m1);
void free_kdtree(void* ptree);
void readmesh(mesh** m1, char* filename);
void read_parameters(char *filename, char* mesh_filename, struct parameters *parameters);
void calculate_parameters(struct parameters *parameters);
int get_no_of_nodes(char* filename);
int get_domain_dimension(char* filename);
void AllocateMemoryMesh(mesh** m1, int nodes);
void free_mesh(mesh* m1);
void calculate_normals(mesh* m1);
void write_normals(mesh* m1, char* filename);
void write_cloud_indices(mesh* m1, char* filename);
void write_boundary_tags(mesh* m1, char* filename);
void write_corner_tags(mesh* m1, char* filename);


// Function definitions
void AllocateMemoryMesh(mesh** m1, int nodes) {
    *m1 = (struct mesh*)malloc(sizeof(struct mesh));
    (*m1)->mesh_filename = (char*)malloc(50*sizeof(char));
    (*m1)->n_count = (int*)malloc(sizeof(int));
    (*m1)->e_count = (int*)malloc(sizeof(int));
    (*m1)->points = (Point*)malloc(nodes * sizeof(Point));
    (*m1)->nx = (double*)malloc(nodes * sizeof(double));
    (*m1)->ny = (double*)malloc(nodes * sizeof(double));
    (*m1)->nz = (double*)malloc(nodes * sizeof(double));
    (*m1)->boundary_tag = (bool*)malloc(nodes * sizeof(bool));
    (*m1)->corner_tag = (bool*)malloc(nodes * sizeof(bool));
    // clound index matrix initialisation
    (*m1)->cloud_index = (int**)malloc(nodes * sizeof(int*));
    for (int i = 0; i < nodes; i++)
        (*m1)->cloud_index[i] = (int*)malloc(parameters.n_cloud_points * sizeof(int));
    *(*m1)->n_count = nodes;
    *(*m1)->e_count = 0;
    }

void free_mesh(mesh* m1) {
    free(m1->n_count);
    free(m1->e_count);
    free(m1->points);
    free(m1->nx);
    free(m1->ny);
    free(m1->boundary_tag);
    free(m1->corner_tag);
    free_matrix_int(m1->cloud_index, *m1->n_count);
    free(m1);
}

void read_parameters(char *filename, char* mesh_filename, struct parameters *parameters) {
    FILE *file;
    char ctemp[100];
    int itemp;
    double dtemp;
    file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("Error: Unable to open the file\n");
        exit(1);
    }
    fscanf(file, "%[^,],%d\n", ctemp, &parameters->poly_degree);
    printf("PARAMETERS: %s = %d\n", ctemp, parameters->poly_degree);
    fscanf(file, "%[^,],%d\n", ctemp, &parameters->phs_degree);
    printf("PARAMETERS: %s = %d\n", ctemp, parameters->phs_degree);
    fscanf(file, "%[^,],%d\n", ctemp, &parameters->cloud_size_multiplier);
    printf("PARAMETERS: %s = %d\n", ctemp, parameters->cloud_size_multiplier);
    fclose(file);
    int nodes = get_no_of_nodes(mesh_filename);    // Get the number of nodes in the mesh
    parameters->dimension = get_domain_dimension(mesh_filename);   // Set the dimension of the problem
    calculate_parameters(parameters);   // Calculate n_poly_terms and n_cloud_points and powers of monomials
}

void calculate_parameters(struct parameters *parameters) {
    if (parameters->dimension == 3)
        parameters->n_poly_terms = (parameters->poly_degree + 1) * (parameters->poly_degree + 2) * (parameters->poly_degree + 3)/ 6;
    else
        parameters->n_poly_terms = (parameters->poly_degree + 1) * (parameters->poly_degree + 2)/ 2;
    parameters->n_cloud_points = parameters->cloud_size_multiplier * parameters->n_poly_terms;
    printf("PARAMETERS: n_poly_terms = %d\n", parameters->n_poly_terms);
    printf("PARAMETERS: n_cloud_points = %d\n", parameters->n_cloud_points);
    printf("PARAMETERS: dimension = %d\n", parameters->dimension);
    parameters->pow_x = (int*)malloc(parameters->n_poly_terms * sizeof(int));
    parameters->pow_y = (int*)malloc(parameters->n_poly_terms * sizeof(int));
    parameters->pow_z = (int*)malloc(parameters->n_poly_terms * sizeof(int));
    int count = 0;
    
    if (parameters->dimension==2)
        for (int i = 0; i <= parameters->poly_degree; i++)
        {
            for (int j = 0; j <= parameters->poly_degree - i; j++)
            {
                parameters->pow_x[count] = i;
                parameters->pow_y[count] = j;
                parameters->pow_z[count] = 0;
                count++;
            }
        }
    else if (parameters->dimension==3)
        for (int i = 0; i <= parameters->poly_degree; i++)
        {
            for (int j = 0; j <= parameters->poly_degree - i; j++)
            {
                for (int k = 0; k <= parameters->poly_degree - i - j; k++)
                {
                    parameters->pow_x[count] = i;
                    parameters->pow_y[count] = j;
                    parameters->pow_z[count] = k;
                    count++;
                }
            }
        }
    }

int get_no_of_nodes(char* filename)
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
    int nodes = itemp;
    fclose(file);
    return nodes;
}

int get_domain_dimension(char* filename)
{   
    FILE *file;
    double sum_z = 0;
    file = fopen(filename, "r");
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
    int nodes = itemp;
    for (int i = 0; i < nodes; i++)
    {
        fscanf(file, "%i ", &itemp); //node number
        fscanf(file, "%lf ", &dtemp);
        fscanf(file, "%lf ", &dtemp); 
        fscanf(file, "%lf ", &dtemp); 
        sum_z += abs(dtemp);    
    }
    fclose(file);
    if (sum_z < 1e-4) //if all z coordinates are zero, then set dimension of the problem to 2
        return 2;
    else
        return 3;
}   


static double dist_sq( double *a1, double *a2, int dims ) { // used in kdtree sorting
    double dist_sq = 0, diff;
    while( --dims >= 0 ) {
        diff = (a1[dims] - a2[dims]);
        dist_sq += diff*diff;
    }
    return dist_sq;
}

void calculate_normals(mesh* m1)
{
    double dx, dy, dz, mag, dot_product;
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
                    dz = m1->points[m1->cloud_index[i][j]].coords[2] - m1->points[i].coords[2];
                    dot_product = dx * m1->nx[i] + dy * m1->ny[i] + dz * m1->nz[i];
                    if (dot_product > 0)
                    {
                        m1->nx[i] = -m1->nx[i];
                        m1->ny[i] = -m1->ny[i];
                        m1->nz[i] = -m1->nz[i];
                    }
                    break;
                }
            }
            mag = sqrt(m1->nx[i] * m1->nx[i] + m1->ny[i] * m1->ny[i] + m1->nz[i] * m1->nz[i]);
            m1->nx[i] = m1->nx[i] / mag;
            m1->ny[i] = m1->ny[i] / mag;
            m1->nz[i] = m1->nz[i] / mag;
        }
    }
}

void* create_kdtree(mesh* m1)
{
    void* ptree;
    ptree = kd_create( 3 ); // create a k-d tree for 3-dimensional points 
    for(int i=0; i<*m1->n_count; i++ ) {
        if(m1->corner_tag[i] == false)
            assert(kd_insert3(ptree, m1->points[i].coords[0], m1->points[i].coords[1], m1->points[i].coords[2], &m1->points[i].index) == 0);
    }
    return ptree;
}

void free_kdtree(void* ptree)
{
    kd_free( ptree );
}

int* find_neighbours(double* p, void* ptree, double radius, struct parameters *parameters){
    double pos[3], dist;
    double pt[3] = { 0.0, 0.0, 0.0 };
    double *distance = (double*)malloc(parameters->n_cloud_points * sizeof(double*));
    int *ind = (int*)malloc(parameters->n_cloud_points * sizeof(int*));
    struct kdres *presults;
    int* pch; // Pointer to the index of the nearest neighbor
    //initialise the distance and index arrays
    for (int i =0; i<parameters->n_cloud_points ; i++) {
        distance[i] = 100.0;
        ind[i] = 0;
    }
    pt[0] = p[0];
    pt[1] = p[1];
    pt[2] = p[2];
    // printf("(%g, %g, %g) -> ", pt[0], pt[1], pt[2]); 
    presults = kd_nearest_range( ptree, pt, radius ); // find points closest to the point pt and within distance radius

    while (!kd_res_end( presults ))
    {   
        pch = (int*)kd_res_item( presults, pos );
        dist = sqrt( dist_sq( pt, pos, 3 ) );
        for (int j = 0; j < parameters->n_cloud_points; ++j) {
            if (dist < distance[j]) {
                for (int k = parameters->n_cloud_points-1; k > j; k--) {
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
    if (distance[parameters->n_cloud_points-1]==100) {
        radius = 1.2*radius;
        kd_res_free( presults );
        free(distance);
        free(ind);
        return find_neighbours(p, ptree, radius, parameters);
    }
    else {
        free(distance);
        kd_res_free( presults );
        return ind;
    }
}

int* find_neighbours_for_mesh(Point* p, void* ptree, double radius, struct parameters *parameters, mesh* m1){
    double pos[3], dist;
    double pt[3] = { 0.0, 0.0, 0.0 };
    double *distance = (double*)malloc(parameters->n_cloud_points * sizeof(double*));
    int *ind = (int*)malloc(parameters->n_cloud_points * sizeof(int*));
    struct kdres *presults;
    int* pch; // Pointer to the index of the nearest neighbor
    //initialise the distance and index arrays
    for (int i =0; i<parameters->n_cloud_points ; i++) {
        distance[i] = 100.0;
        ind[i] = 0;
    }
    pt[0] = p->coords[0];
    pt[1] = p->coords[1];
    pt[2] = p->coords[2];
    // printf("(%g, %g, %g) -> ", pt[0], pt[1], pt[2]); 
    presults = kd_nearest_range( ptree, pt, radius ); // find points closest to the point pt and within distance radius

    while (!kd_res_end( presults ))
    {   
        pch = (int*)kd_res_item( presults, pos );
        if (m1->boundary_tag[p->index] && m1->boundary_tag[*pch] && !(p->index == *pch)) {
            kd_res_next( presults );
            continue;
        }
        dist = sqrt( dist_sq( pt, pos, 3 ) );
        for (int j = 0; j < parameters->n_cloud_points; ++j) {
            if (dist < distance[j]) {
                for (int k = parameters->n_cloud_points-1; k > j; k--) {
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
    if (distance[parameters->n_cloud_points-1]==100) {
        radius = 1.2*radius;
        kd_res_free( presults );
        free(distance);
        free(ind);
        return find_neighbours_for_mesh(p, ptree, radius, parameters, m1);
    }
    else {
        free(distance);
        kd_res_free( presults );
        return ind;
    }
}

void readmesh(mesh** m1, char* filename)
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
    AllocateMemoryMesh(m1, itemp);
    strcpy((*m1)->mesh_filename, filename);
    
    printf("Reading mesh file %s\n", (*m1)->mesh_filename);
    
    for (int i = 0; i < *(*m1)->n_count; i++)
    {
        fscanf(file, "%i ", &itemp); //node number
        fscanf(file, "%lf ", &dtemp);
        (*m1)->points[i].coords[0]=dtemp;
        fscanf(file, "%lf ", &dtemp); 
        (*m1)->points[i].coords[1]=dtemp;
        fscanf(file, "%lf ", &dtemp); 
        (*m1)->points[i].coords[2]=dtemp;
        //printf("%d %lf %lf %lf\n", i, (*m1)->points[i].coords[0], (*m1)->points[i].coords[1], (*m1)->points[i].coords[2]);
        
        (*m1)->points[i].index=i; 
        (*m1)->boundary_tag[i]=false;
        (*m1)->corner_tag[i]=false;
        (*m1)->nx[i]=0;
        (*m1)->ny[i]=0;
        (*m1)->nz[i]=0;
        sum_z += abs((*m1)->points[i].coords[2]);

        if ((*m1)->points[i].coords[0] < x_min)
            x_min = (*m1)->points[i].coords[0];
        if ((*m1)->points[i].coords[1] < y_min)
            y_min = (*m1)->points[i].coords[1];
        if ((*m1)->points[i].coords[2] < z_min)
            z_min = (*m1)->points[i].coords[2];
        if ((*m1)->points[i].coords[0] > x_max)
            x_max = (*m1)->points[i].coords[0];
        if ((*m1)->points[i].coords[1] > y_max)
            y_max = (*m1)->points[i].coords[1];
        if ((*m1)->points[i].coords[2] > z_max)
            z_max = (*m1)->points[i].coords[2];     
    }
    (*m1)->x_min = x_min;
    (*m1)->y_min = y_min;
    (*m1)->z_min = z_min;
    (*m1)->x_max = x_max;
    (*m1)->y_max = y_max;
    (*m1)->z_max = z_max;

    while (true)
    {
        fscanf(file, "%s ", temp);
        if (strcmp(temp, "$Elements") == 0)
            break;
    }
    fscanf(file, "%i ", &itemp);
    *(*m1)->e_count = itemp;

    if (sum_z < 1e-4){ //if all z coordinates are zero, then set dimension of the problem to 2
        parameters.dimension = 2;
        (*m1)->d_avg = sqrt((x_max - x_min)*(y_max - y_min) / *(*m1)->e_count);
        printf("x_min = %lf, x_max = %lf\n", x_min, x_max);
        printf("y_min = %lf, y_max = %lf\n", y_min, y_max);
        // printf("e_count = %d\n", *(*m1)->e_count);
        printf("Average distance between nodes = %lf\n", (*m1)->d_avg);
        }
    else {
        parameters.dimension = 3;
        (*m1)->d_avg = sqrt((x_max - x_min)*(y_max - y_min)*(z_max - z_min) / *(*m1)->e_count);
        printf("x_min = %lf, x_max = %lf\n", x_min, x_max);
        printf("y_min = %lf, y_max = %lf\n", y_min, y_max);
        printf("z_min = %lf, z_max = %lf\n", z_min, z_max);
        printf("average distance between nodes = %lf\n", (*m1)->d_avg);
    }

    // Read the mesh file elements and boundary nodes
    int e_type, tag_int, e_node1, e_node2, e_node3, e_node4; //variables for reading elements
    double dx,dy,dz,dx1,dy1,dz1; // vectors for calculating normal to the surface
    int count_e = 0;
    (*m1)->d_avg = 0;
    if (parameters.dimension==2) //if 2D problem, then read only line elements as boundary elements
    {
        for (int ie = 0; ie < *(*m1)->e_count; ie++)
        {
            fscanf(file, "%i ", &itemp);   //element number
            fscanf(file, "%i ", &e_type); //type of the element
        
            if (e_type == 1) //reading line elements
            { 
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); // edge/boundary tag, if a square, elements on one edge will have same tag
                fscanf(file, "%i ", &e_node1);     //node number
                (*m1)->boundary_tag[e_node1 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node2);     //node number
                (*m1)->boundary_tag[e_node2 - 1] = true; //set boundary tag to true
                dx = (*m1)->points[e_node2 - 1].coords[0] - (*m1)->points[e_node1 - 1].coords[0]; //calculate normal to the line
                dy = (*m1)->points[e_node2 - 1].coords[1] - (*m1)->points[e_node1 - 1].coords[1]; 
                (*m1)->nx[e_node1 - 1] += -dy; (*m1)->ny[e_node1 - 1] += dx; //add the normals to each node
                (*m1)->nx[e_node2 - 1] += -dy; (*m1)->ny[e_node2 - 1] += dx; //add the normals to each node
                (*m1)->d_avg += sqrt(dx*dx + dy*dy); //average distance between nodes
                count_e++;
            }
            else if (e_type == 15) //reading corner nodes
            { 
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); // edge/boundary tag; if a square domain, elements on one edge will have same tag
                fscanf(file, "%i ", &itemp);   //node number
                (*m1)->corner_tag[itemp - 1] = true; //set corner tag to true
            }
            else
            {                              //not reading any other kind of elements
                fscanf(file, "%*[^\n]\n"); //skip reading remaining row
            }
        }
    }
    else //if 3D problem, then read only surface elements as boundary elements
    {
        for (int ie = 0; ie < *(*m1)->e_count; ie++)
        {
            fscanf(file, "%i ", &itemp);   //element number
            fscanf(file, "%i ", &e_type); //element type
            if (e_type == 2)
            { //important 3 node triangle CV on the boundary
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); //element tag
                fscanf(file, "%i ", &e_node1);     //node number
                (*m1)->boundary_tag[e_node1 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node2);     //node number
                (*m1)->boundary_tag[e_node2 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node3);     //node number
                (*m1)->boundary_tag[e_node3 - 1] = true; //set boundary tag to true
                //calculate normal to the surface of the triangle
                // by calculating cross product of the two vectors
                dx = (*m1)->points[e_node2 - 1].coords[0] - (*m1)->points[e_node1 - 1].coords[0]; 
                dy = (*m1)->points[e_node2 - 1].coords[1] - (*m1)->points[e_node1 - 1].coords[1];
                dz = (*m1)->points[e_node2 - 1].coords[2] - (*m1)->points[e_node1 - 1].coords[2];
                dx1 = (*m1)->points[e_node3 - 1].coords[0] - (*m1)->points[e_node1 - 1].coords[0];
                dy1 = (*m1)->points[e_node3 - 1].coords[1] - (*m1)->points[e_node1 - 1].coords[1];
                dz1 = (*m1)->points[e_node3 - 1].coords[2] - (*m1)->points[e_node1 - 1].coords[2];
                dtemp = dy*dz1 - dz*dy1;
                (*m1)->nx[e_node1 - 1] += dtemp;
                (*m1)->nx[e_node2 - 1] += dtemp;
                (*m1)->nx[e_node3 - 1] += dtemp;
                dtemp = dz*dx1 - dx*dz1;
                (*m1)->ny[e_node1 - 1] += dtemp;
                (*m1)->ny[e_node2 - 1] += dtemp;
                (*m1)->ny[e_node3 - 1] += dtemp;
                dtemp = dx*dy1 - dy*dx1;
                (*m1)->nz[e_node1 - 1] += dtemp;
                (*m1)->nz[e_node2 - 1] += dtemp;
                (*m1)->nz[e_node3 - 1] += dtemp;
                (*m1)->d_avg += sqrt(dx*dx + dy*dy + dz*dz); //average distance between nodes
                count_e++;
                (*m1)->d_avg += sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1); //average distance between nodes
                count_e++;
            }
            else if (e_type == 3)
            { //important 4 node quad element on the boundary
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); //element tag
                fscanf(file, "%i ", &e_node1);     //node number
                (*m1)->boundary_tag[e_node1 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node2);     //node number
                (*m1)->boundary_tag[e_node2 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node3);     //node number
                (*m1)->boundary_tag[e_node3 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node4);     //node number
                (*m1)->boundary_tag[e_node4 - 1] = true; //set boundary tag to true
                // calculate normal to the surface of the quad element
                // by calculating cross product of any two vectors/lines on the surface
                dx = (*m1)->points[e_node2 - 1].coords[0] - (*m1)->points[e_node1 - 1].coords[0];
                dy = (*m1)->points[e_node2 - 1].coords[1] - (*m1)->points[e_node1 - 1].coords[1];
                dz = (*m1)->points[e_node2 - 1].coords[2] - (*m1)->points[e_node1 - 1].coords[2];
                dx1 = (*m1)->points[e_node3 - 1].coords[0] - (*m1)->points[e_node1 - 1].coords[0];
                dy1 = (*m1)->points[e_node3 - 1].coords[1] - (*m1)->points[e_node1 - 1].coords[1];
                dz1 = (*m1)->points[e_node3 - 1].coords[2] - (*m1)->points[e_node1 - 1].coords[2];
                dtemp = dy*dz1 - dz*dy1;
                (*m1)->nx[e_node1 - 1] += dtemp;
                (*m1)->nx[e_node2 - 1] += dtemp;
                (*m1)->nx[e_node3 - 1] += dtemp;
                (*m1)->nx[e_node4 - 1] += dtemp;
                dtemp = dz*dx1 - dx*dz1;
                (*m1)->ny[e_node1 - 1] += dtemp;
                (*m1)->ny[e_node2 - 1] += dtemp;
                (*m1)->ny[e_node3 - 1] += dtemp;
                (*m1)->ny[e_node4 - 1] += dtemp;
                dtemp = dx*dy1 - dy*dx1;
                (*m1)->nz[e_node1 - 1] += dtemp;
                (*m1)->nz[e_node2 - 1] += dtemp;
                (*m1)->nz[e_node3 - 1] += dtemp;
                (*m1)->nz[e_node4 - 1] += dtemp;
                (*m1)->d_avg += sqrt(dx*dx + dy*dy + dz*dz); //average distance between nodes
                count_e++;
                (*m1)->d_avg += sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1); //average distance between nodes
                count_e++;
            }
            else if (e_type == 1)
            { //2 node line: log these 2 vertex numbers to delete them later
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); //element tag
                fscanf(file, "%i ", &itemp);   //first point
                (*m1)->corner_tag[itemp - 1] = true;
                fscanf(file, "%i ", &itemp); //second point
                (*m1)->corner_tag[itemp - 1] = true;
            }
            else if (e_type == 15) //reading corner nodes
            { 
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); // edge/boundary tag; if a square domain, elements on one edge will have same tag
                fscanf(file, "%i ", &itemp);   //node number
                (*m1)->corner_tag[itemp - 1] = true; //set corner tag to true
            }
            else
            {                              //not reading any other kind of element
                fscanf(file, "%*[^\n]\n"); //skip reading remaining row
            }
        }
    }
    (*m1)->d_avg = (*m1)->d_avg / count_e;
    fclose(file);
    printf("No of nodes = %d \nNo of elements = %d \n", *(*m1)->n_count, *(*m1)->e_count);
    
    // for (int i =0; i<*(*m1)->n_count; i++)
    // {
    //     (*m1)->points[i].coords[0] = (*m1)->points[i].coords[0] +1;
    //     (*m1)->points[i].coords[1] = (*m1)->points[i].coords[1] +1;
    // }

    printf("\nBuilding KDTree...\n");
    void* ptree=create_kdtree((*m1));

    printf("Searching for neighbours...\n");
    // Search for neighbours
    double radius = parameters.n_cloud_points*(*m1)->d_avg; // Search radius
    printf("Search radius = %g\n", radius);
    int *ind = (int*)malloc(parameters.n_cloud_points * sizeof(int*));

    for(int i=0; i<*(*m1)->n_count; i++ ) {
        ind = find_neighbours_for_mesh(&(*m1)->points[i], ptree, radius, &parameters, *m1);
        for (int j = 0; j < parameters.n_cloud_points; j++) {
            (*m1)->cloud_index[i][j] = ind[j];
        }
    }

    free_kdtree(ptree);
    free(ind);
    calculate_normals((*m1));  // Calculate the normals at boundary nodes
}

void write_cloud_index(mesh* m1, char* filename)
{
    FILE *file;
    printf("Writing cloud index to file %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for(int i=0; i<*m1->n_count; i++ ) {
        for (int j = 0; j < parameters.n_cloud_points; j++)
            fprintf(file, "%d,", m1->cloud_index[i][j]);
        fprintf(file, "\n");
    }
    fclose(file);
    printf("Cloud index written to file %s!\n\n", filename);
}

void write_normals(mesh* m1, char* filename)
{
    FILE *file;
    printf("Writing normals to file %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < *m1->n_count; i++)
        fprintf(file, "%d %lf %lf %lf\n", i, m1->nx[i], m1->ny[i], m1->nz[i]);
    fclose(file);
}

void write_boundary_tags(mesh* m1, char* filename)
{
    FILE *file;
    printf("Writing boundary tags to file %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < *m1->n_count; i++)
        fprintf(file, "%d %d\n", i, m1->boundary_tag[i]);
    fclose(file);
}

void write_corner_tags(mesh* m1, char* filename)
{
    FILE *file;
    printf("Writing corner tags to file %s\n", filename);
    file = fopen(filename, "w");
    if(file==NULL)
    {
        printf("Error: Unable to open the file %s\n",filename);
        exit(1);
    }
    for (int i = 0; i < *m1->n_count; i++)
        fprintf(file, "%d %d\n", i, m1->corner_tag[i]);
    fclose(file);
}
