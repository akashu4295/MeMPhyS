// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef MESH_FUNCTIONS_H
#define MESH_FUNCTIONS_H

#include "structures.h"
#include "kdtree_functions.h"
#include "general_functions.h"
#include "write_functions.h"
#include "rbf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double dtemp; int itemp; char temp[50]; // temporary variables used multiple times in the code

///////////////////////////////////////////////////////////////////////////////
// Function definitions
///////////////////////////////////////////////////////////////////////////////

void read_PointStructure(PointStructure* myPointStruct);
void read_flow_parameters(char *filename);
void calculate_parameters(PointStructure* myPointStruct);
int get_no_of_nodes(char* filename);
void correct_normal_directions(PointStructure* myPointStruct);
void read_grid_filenames(PointStructure** myPointStruct, char* filename, short* num_levels);
void read_complete_mesh_data(PointStructure* myPointStruct, short num_levels);   
void create_restriction_matrix(PointStructure* myPointStruct, PointStructure* myPointStruct1);
void create_prolongation_matrix(PointStructure* myPointStruct, PointStructure* myPointStruct1);

///////////////////////////////////////////////////////////////////////////////
// Function definitions
///////////////////////////////////////////////////////////////////////////////

void read_flow_parameters(char *filename) {
    FILE *file;
    char ctemp[100];

    file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("Error: Unable to open the file\n");
        exit(1);
    }
    fscanf(file, "%[^,],%d\n", ctemp, &parameters.dimension);
    printf("%s = %d\n", ctemp, parameters.dimension);
    fscanf(file, "%[^,],%d\n", ctemp, &parameters.poly_degree);
    printf("PARAMETERS: %s = %d\n", ctemp, parameters.poly_degree);
    fscanf(file, "%[^,],%d\n", ctemp, &parameters.phs_degree);
    printf("PARAMETERS: %s = %d\n", ctemp, parameters.phs_degree);
    fscanf(file, "%[^,],%d\n", ctemp, &parameters.cloud_size_multiplier);
    printf("PARAMETERS: %s = %d\n", ctemp, parameters.cloud_size_multiplier);
    fscanf(file, "%[^,],%d\n", ctemp, &parameters.test);
    printf("PARAMETERS: %s = %d\n", ctemp, parameters.test);
    fclose(file);
}

void calculate_parameters(PointStructure* myPointStruct) {
    if (parameters.dimension == 3)
        myPointStruct->num_poly_terms = (myPointStruct->poly_degree + 1) * (myPointStruct->poly_degree + 2) * (myPointStruct->poly_degree + 3)/ 6;
    else
        myPointStruct->num_poly_terms = (myPointStruct->poly_degree + 1) * (myPointStruct->poly_degree + 2)/ 2;
    myPointStruct->num_cloud_points = parameters.cloud_size_multiplier * myPointStruct->num_poly_terms;
    printf("PARAMETERS: num_poly_terms = %d\n", myPointStruct->num_poly_terms);
    printf("PARAMETERS: num_cloud_points = %d\n", myPointStruct->num_cloud_points);
    printf("PARAMETERS: dimension = %d\n", parameters.dimension);
    myPointStruct->pow_x = (short*)malloc(myPointStruct->num_poly_terms * sizeof(short));
    myPointStruct->pow_y = (short*)malloc(myPointStruct->num_poly_terms * sizeof(short));
    myPointStruct->pow_z = (short*)malloc(myPointStruct->num_poly_terms * sizeof(short));
    int count = 0;
    
    if (parameters.dimension==2)
        for (int i = 0; i <= myPointStruct->poly_degree; i++)
        {
            for (int j = 0; j <= myPointStruct->poly_degree - i; j++)
            {
                myPointStruct->pow_x[count] = i;
                myPointStruct->pow_y[count] = j;
                myPointStruct->pow_z[count] = 0;
                count++;
            }
        }
    else if (parameters.dimension==3)
        for (int i = 0; i <= myPointStruct->poly_degree; i++)
        {
            for (int j = 0; j <= myPointStruct->poly_degree - i; j++)
            {
                for (int k = 0; k <= myPointStruct->poly_degree - i - j; k++)
                {
                    myPointStruct->pow_x[count] = i;
                    myPointStruct->pow_y[count] = j;
                    myPointStruct->pow_z[count] = k;
                    count++;
                }
            }
        }
}

int get_no_of_nodes(char* filename) {   
    FILE *file;
    double sum_z = 0;
    file = fopen(filename, "r");
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

void read_PointStructure(PointStructure* myPointStruct)
{   
    FILE *file;
    
    char filename[50];
    sprintf(filename, myPointStruct->mesh_filename);

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

    // Allocate memory for the members of point structure
    AllocateMemoryPointStructure(myPointStruct, itemp);
    
    for (int i = 0; i < myPointStruct->num_nodes; i++)
    {
        fscanf(file, "%i ", &itemp); //node number
        fscanf(file, "%lf ", &dtemp);
        myPointStruct->x[i]=dtemp;
        fscanf(file, "%lf ", &dtemp); 
        myPointStruct->y[i]=dtemp;
        fscanf(file, "%lf ", &dtemp); 
        myPointStruct->z[i]=dtemp;
        
        myPointStruct->point_index[i] = i;
        myPointStruct->boundary_tag[i]=false;
        myPointStruct->corner_tag[i]=false;
        myPointStruct->x_normal[i]=0;
        myPointStruct->y_normal[i]=0;
        myPointStruct->z_normal[i]=0;

        if (myPointStruct->x[i] < x_min)
            x_min = myPointStruct->x[i];
        if (myPointStruct->y[i] < y_min)
            y_min = myPointStruct->y[i];
        if (myPointStruct->z[i] < z_min)
            z_min = myPointStruct->z[i];
        if (myPointStruct->x[i] > x_max)
            x_max = myPointStruct->x[i];
        if (myPointStruct->y[i] > y_max)
            y_max = myPointStruct->y[i];
        if (myPointStruct->z[i] > z_max)
            z_max = myPointStruct->z[i];
    }

    while (true)
    {
        fscanf(file, "%s ", temp);
        if (strcmp(temp, "$Elements") == 0)
            break;
    }
    fscanf(file, "%i ", &itemp);
    myPointStruct->num_elem = itemp;

    if (parameters.dimension == 2){ //if all z coordinates are zero, then set dimension of the problem to 2
        myPointStruct->d_avg = sqrt((x_max - x_min)*(y_max - y_min) / myPointStruct->num_elem);
        printf("Average distance between nodes = %lf\n", myPointStruct->d_avg);
        }
    else {
        myPointStruct->d_avg = cbrt((x_max - x_min)*(y_max - y_min)*(z_max - z_min)/myPointStruct->num_elem);
        printf("Average distance between nodes = %lf\n", myPointStruct->d_avg);
    }

    // Read the mesh file elements and boundary nodes
    int e_type, tag_int, e_node1, e_node2, e_node3, e_node4; //variables for reading elements
    double dx,dy,dz,dx1,dy1,dz1; // vectors for calculating normal to the surface

    if (parameters.dimension==2) //if 2D problem, then read only line elements as boundary elements
    {
        for (int ie = 0; ie < myPointStruct->num_elem; ie++)
        {
            fscanf(file, "%i ", &itemp);   //element number
            fscanf(file, "%i ", &e_type); //type of the element
        
            if (e_type == 1) //reading line elements
            { 
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); // edge/boundary tag, if a square, elements on one edge will have same tag
                fscanf(file, "%i ", &e_node1);     //node number
                myPointStruct->boundary_tag[e_node1 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node2);     //node number
                myPointStruct->boundary_tag[e_node2 - 1] = true; //set boundary tag to true
                dx = myPointStruct->x[e_node2 - 1] - myPointStruct->x[e_node1 - 1]; //calculate normal to the line
                dy = myPointStruct->y[e_node2 - 1] - myPointStruct->y[e_node1 - 1]; 
                myPointStruct->x_normal[e_node1 - 1] += -dy; myPointStruct->y_normal[e_node1 - 1] += dx; //add the normals to each node
                myPointStruct->x_normal[e_node2 - 1] += -dy; myPointStruct->y_normal[e_node2 - 1] += dx; //add the normals to each node
            }
            else if (e_type == 15) //reading corner nodes
            { 
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); // edge/boundary tag; if a square domain, elements on one edge will have same tag
                fscanf(file, "%i ", &itemp);   //node number
                myPointStruct->corner_tag[itemp - 1] = true; //set corner tag to true
            }
            else
            {                              //not reading any other kind of elements
                fscanf(file, "%*[^\n]\n"); //skip reading remaining row
            }
        }
    }
    else //if 3D problem, then read only surface elements as boundary elements
    {   
        for (int ie = 0; ie < myPointStruct->num_elem; ie++)
        {
            fscanf(file, "%i ", &itemp);   //element number
            fscanf(file, "%i ", &e_type); //element type
            if (e_type == 2)
            { //important 3 node triangle CV on the boundary
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); //element tag
                fscanf(file, "%i ", &e_node1);     //node number
                myPointStruct->boundary_tag[e_node1 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node2);     //node number
                myPointStruct->boundary_tag[e_node2 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node3);     //node number
                myPointStruct->boundary_tag[e_node3 - 1] = true; //set boundary tag to true
                //calculate normal to the surface of the triangle
                // by calculating cross product of the two vectors
                dx = myPointStruct->x[e_node2 - 1] - myPointStruct->x[e_node1 - 1]; 
                dy = myPointStruct->y[e_node2 - 1] - myPointStruct->y[e_node1 - 1];
                dz = myPointStruct->z[e_node2 - 1] - myPointStruct->z[e_node1 - 1];
                dx1 = myPointStruct->x[e_node3 - 1] - myPointStruct->x[e_node1 - 1];
                dy1 = myPointStruct->y[e_node3 - 1] - myPointStruct->y[e_node1 - 1];
                dz1 = myPointStruct->z[e_node3 - 1] - myPointStruct->z[e_node1 - 1];
                dtemp = dy*dz1 - dz*dy1;
                myPointStruct->x_normal[e_node1 - 1] += dtemp;
                myPointStruct->x_normal[e_node2 - 1] += dtemp;
                myPointStruct->x_normal[e_node3 - 1] += dtemp;
                dtemp = dz*dx1 - dx*dz1;
                myPointStruct->y_normal[e_node1 - 1] += dtemp;
                myPointStruct->y_normal[e_node2 - 1] += dtemp;
                myPointStruct->y_normal[e_node3 - 1] += dtemp;
                dtemp = dx*dy1 - dy*dx1;
                myPointStruct->z_normal[e_node1 - 1] += dtemp;
                myPointStruct->z_normal[e_node2 - 1] += dtemp;
                myPointStruct->z_normal[e_node3 - 1] += dtemp;
            }
            else if (e_type == 3)
            { //important 4 node quad element on the boundary
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); //element tag
                fscanf(file, "%i ", &e_node1);     //node number
                myPointStruct->boundary_tag[e_node1 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node2);     //node number
                myPointStruct->boundary_tag[e_node2 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node3);     //node number
                myPointStruct->boundary_tag[e_node3 - 1] = true; //set boundary tag to true
                fscanf(file, "%i ", &e_node4);     //node number
                myPointStruct->boundary_tag[e_node4 - 1] = true; //set boundary tag to true
                // calculate normal to the surface of the quad element
                // by calculating cross product of any two vectors/lines on the surface
                dx = myPointStruct->x[e_node2 - 1] - myPointStruct->x[e_node1 - 1];
                dy = myPointStruct->y[e_node2 - 1] - myPointStruct->y[e_node1 - 1];
                dz = myPointStruct->z[e_node2 - 1] - myPointStruct->z[e_node1 - 1];
                dx1 = myPointStruct->x[e_node3 - 1] - myPointStruct->x[e_node1 - 1];
                dy1 = myPointStruct->y[e_node3 - 1] - myPointStruct->y[e_node1 - 1];
                dz1 = myPointStruct->z[e_node3 - 1] - myPointStruct->z[e_node1 - 1];
                dtemp = dy*dz1 - dz*dy1;
                myPointStruct->x_normal[e_node1 - 1] += dtemp;
                myPointStruct->x_normal[e_node2 - 1] += dtemp;
                myPointStruct->x_normal[e_node3 - 1] += dtemp;
                myPointStruct->x_normal[e_node4 - 1] += dtemp;
                dtemp = dz*dx1 - dx*dz1;
                myPointStruct->y_normal[e_node1 - 1] += dtemp;
                myPointStruct->y_normal[e_node2 - 1] += dtemp;
                myPointStruct->y_normal[e_node3 - 1] += dtemp;
                myPointStruct->y_normal[e_node4 - 1] += dtemp;
                dtemp = dx*dy1 - dy*dx1;
                myPointStruct->z_normal[e_node1 - 1] += dtemp;
                myPointStruct->z_normal[e_node2 - 1] += dtemp;
                myPointStruct->z_normal[e_node3 - 1] += dtemp;
                myPointStruct->z_normal[e_node4 - 1] += dtemp;
            }
            else if (e_type == 1)
            { //2 node line: log these 2 vertex numbers to delete them later
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); //element tag
                fscanf(file, "%i ", &itemp);   //first point
                myPointStruct->corner_tag[itemp - 1] = true;
                fscanf(file, "%i ", &itemp); //second point
                myPointStruct->corner_tag[itemp - 1] = true;
            }
            else if (e_type == 15) //reading corner nodes
            { 
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &itemp);
                fscanf(file, "%i ", &tag_int); // edge/boundary tag; if a square domain, elements on one edge will have same tag
                fscanf(file, "%i ", &itemp);   //node number
                myPointStruct->corner_tag[itemp - 1] = true; //set corner tag to true
            }
            else
            {                              //not reading any other kind of element
                fscanf(file, "%*[^\n]\n"); //skip reading remaining row
            }
        }
    }
    fclose(file);
    myPointStruct->num_corners = 0;
    myPointStruct->num_boundary_nodes = 0;
    for (int i = 0; i < myPointStruct->num_nodes; i++)
    {
        if (myPointStruct->boundary_tag[i] == true)
        {
            myPointStruct->num_boundary_nodes++;
        }
        if (myPointStruct->corner_tag[i] == true)
        {
            myPointStruct->num_corners++;
        }
    }
    printf("No of nodes = %d \nNo of elements = %d \n", myPointStruct->num_nodes, myPointStruct->num_elem);
    printf("No of boundary nodes = %d \nNo of corner nodes = %d \n", myPointStruct->num_boundary_nodes, myPointStruct->num_corners);
}

void correct_normal_directions(PointStructure *myPointStruct)
{
    double dx, dy, dz, mag, dot_product, x_avg, y_avg, z_avg;
    short count = 0;
    for (int i = 0; i < myPointStruct->num_nodes; i++)
    {   
        if(myPointStruct->boundary_tag[i] == true)
        {   
            x_avg = 0; y_avg = 0; z_avg = 0; count = 0;
            // for (int j = 0; j < myPointStruct->num_cloud_points; j++){
            //     if ((myPointStruct->boundary_tag[myPointStruct->cloud_index[i][j]] == false)&&(myPointStruct->corner_tag[myPointStruct->cloud_index[i][j]] == false))
            //     {
            //         x_avg += myPointStruct->x[myPointStruct->cloud_index[i][j]];
            //         y_avg += myPointStruct->y[myPointStruct->cloud_index[i][j]];
            //         z_avg += myPointStruct->z[myPointStruct->cloud_index[i][j]];
            //         count++;
            //         // if (count == 10)
            //         //     break;
            //     }
            // }
            for (int j = myPointStruct->num_cloud_points-1; j > 0; j--){
                if (myPointStruct->boundary_tag[myPointStruct->cloud_index[i][j]] == false){
                    x_avg += myPointStruct->x[myPointStruct->cloud_index[i][j]];
                    y_avg += myPointStruct->y[myPointStruct->cloud_index[i][j]];
                    z_avg += myPointStruct->z[myPointStruct->cloud_index[i][j]];
                    count++;
                    if (count == 1)
                        break;
                }
            }
            x_avg = x_avg/count; y_avg = y_avg/count; z_avg = z_avg/count;

            dx = x_avg - myPointStruct->x[i];
            dy = y_avg - myPointStruct->y[i];
            dz = z_avg - myPointStruct->z[i];
            dot_product = dx * myPointStruct->x_normal[i] + dy * myPointStruct->y_normal[i] + dz * myPointStruct->z_normal[i];
            if (dot_product > 0)
            {
                myPointStruct->x_normal[i] = -myPointStruct->x_normal[i];
                myPointStruct->y_normal[i] = -myPointStruct->y_normal[i];
                myPointStruct->z_normal[i] = -myPointStruct->z_normal[i];
            }

            mag = sqrt(myPointStruct->x_normal[i] * myPointStruct->x_normal[i] + myPointStruct->y_normal[i] * myPointStruct->y_normal[i] + myPointStruct->z_normal[i] * myPointStruct->z_normal[i]);
            myPointStruct->x_normal[i] = myPointStruct->x_normal[i] / mag;
            myPointStruct->y_normal[i] = myPointStruct->y_normal[i] / mag;
            myPointStruct->z_normal[i] = myPointStruct->z_normal[i] / mag;
        }
    }
}

void read_grid_filenames(PointStructure** myPointStruct, char* filename, short* num_levels)
{   
    FILE *file;
    file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("Error: Unable to open the file\n");
        exit(1);
    }
    fscanf(file, "%[^,],%d\n", temp, num_levels);
    printf("%s = %d\n", temp, *num_levels);
    // Allocate memory for point structure for all levels and read meshfile names
    *myPointStruct = (PointStructure*)malloc((*num_levels) * sizeof(PointStructure));
    for (short ii = 0; ii<*num_levels ; ii = ii +1)
        fscanf(file, "%s\n", (*myPointStruct)[ii].mesh_filename);
    fclose(file);
}

void read_complete_mesh_data(PointStructure* myPointStruct, short num_levels)
{   
    for (short ii = 0; ii<num_levels ; ii = ii +1){
        printf("\nReading mesh data for level %d\n", ii+1);
        calculate_parameters(&myPointStruct[ii]);
        read_PointStructure(&myPointStruct[ii]);
        find_cloud_index(&myPointStruct[ii], &myPointStruct[ii]);
        correct_normal_directions(&myPointStruct[ii]);
    }
    
    printf("\nIdentifying prolongation and restriction points\n");
    for (short ii = 0; ii<num_levels ; ii = ii +1){
        printf("Level %d\n", ii+1);
        if (ii == num_levels-1)
            myPointStruct[ii].prolongation_points = find_nearest_point(&myPointStruct[ii], &myPointStruct[ii], 1);
        else{
            myPointStruct[ii].prolongation_points = find_nearest_point(&myPointStruct[ii], &myPointStruct[ii+1], 1);
            create_prolongation_matrix(&myPointStruct[ii], &myPointStruct[ii+1]);
        }
        if (ii == 0)
            myPointStruct[ii].restriction_points = find_nearest_point(&myPointStruct[ii], &myPointStruct[ii], 1);  
        else{
            myPointStruct[ii].restriction_points = find_nearest_point(&myPointStruct[ii], &myPointStruct[ii-1], 1);
            create_restriction_matrix(&myPointStruct[ii-1], &myPointStruct[ii]);
        }
    }
    printf("Prolongation and restriction points identified\n");
}

//////// Restriction matrix creation //////////////

void create_restriction_matrix(PointStructure* myPointStruct_f, 
                                    PointStructure* myPointStruct_c)
{
    // allocate memory for restriction matrix
    
    short m = myPointStruct_f->num_cloud_points;
    short n = myPointStruct_f->num_poly_terms;
    short mpn = m+n;

    myPointStruct_c->restr_mat = (double**)malloc(myPointStruct_c->num_nodes * sizeof(double*));

    for (int i = 0; i < myPointStruct_c->num_nodes; i++)
        myPointStruct_c->restr_mat[i] = (double*)malloc(m * sizeof(double));

    // Initialise to zeros
    for (int i = 0; i < myPointStruct_c->num_nodes; i++)
    {
        for (short j = 0; j < m; j++)
        {
            myPointStruct_c->restr_mat[i][j] = 0;
        }
    }
    
    double **A_inv = create_matrix1(mpn,mpn);
    double **A = create_matrix1(mpn,mpn);
    double *temp = create_vector(mpn);
    double *temp1 = create_vector(m); // coefficients of restriction matrix
    int i_restr;

    double point1[3], point2[3];
    for (int i = myPointStruct_c->num_boundary_nodes; 
                    i < myPointStruct_c->num_nodes; i++)
    {   
        i_restr = myPointStruct_c->restriction_points[i];
        create_A_matrix_from_cloud_indices(myPointStruct_f, A, myPointStruct_f->cloud_index[i_restr]);
        matrixInverse_Gauss_Jordan(A, A_inv, m+n);
        
        point1[0] = myPointStruct_c->x[i];
        point1[1] = myPointStruct_c->y[i];
        point1[2] = myPointStruct_c->z[i];
        
        for (short j = 0; j < m; j++) {
            point2[0] = myPointStruct_f->x[myPointStruct_f->cloud_index[i_restr][j]];
            point2[1] = myPointStruct_f->y[myPointStruct_f->cloud_index[i_restr][j]];
            point2[2] = myPointStruct_f->z[myPointStruct_f->cloud_index[i_restr][j]];
            temp[j] = calculate_phs_rbf(point1, point2, parameters.phs_degree, 
                                                        parameters.dimension);
        }
        double seed_pt[3] = {myPointStruct_f->x[i_restr],
                             myPointStruct_f->y[i_restr],
                             myPointStruct_f->z[i_restr]};
        for (short j = 0; j < n; j++) 
            temp[j+m] = pow(point1[0]-seed_pt[0], myPointStruct_f->pow_x[j])
                                * pow(point1[1]-seed_pt[1], myPointStruct_f->pow_y[j])
                                * pow(point1[2]-seed_pt[2], myPointStruct_f->pow_z[j]);
        
        multiply_vector_matrix_columnwise(temp, A_inv, temp1, m+n, m);
        
        for (short j = 0; j < m; j++)
        {
            myPointStruct_c->restr_mat[i][j] = temp1[j];
        }
    }
    free(temp);
    free(temp1);
    free_matrix(A, mpn);
    free_matrix(A_inv, mpn);
}

void create_prolongation_matrix(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c)
{
    // allocate memory for restriction matrix
    
    short m = myPointStruct_c->num_cloud_points;
    short n = myPointStruct_c->num_poly_terms;
    short mpn = m+n;
    
    myPointStruct_f->prol_mat = (double**)malloc(myPointStruct_f->num_nodes * sizeof(double*));
    
    for (int i = 0;
                     i < myPointStruct_f->num_nodes; i++)
        myPointStruct_f->prol_mat[i] = (double*)malloc(m * sizeof(double));

    // Initialise to zeros
    for (int i = 0;
                     i < myPointStruct_f->num_nodes; i++)
    {
        for (short j = 0; j < m; j++)
        {
            myPointStruct_f->prol_mat[i][j] = 0;
        }
    }
    
    double **A_inv = create_matrix1(mpn,mpn);
    double **A = create_matrix1(mpn,mpn);
    double *temp = create_vector(mpn);
    double *temp1 = create_vector(m);
    int i_prol;

    double point1[3], point2[3], seed_pt[3];
    for (int i = myPointStruct_f->num_boundary_nodes;
                     i < myPointStruct_f->num_nodes; i++)
    {   
        i_prol = myPointStruct_f->prolongation_points[i];
        create_A_matrix_from_cloud_indices(myPointStruct_c, A, myPointStruct_c->cloud_index[i_prol]);
        matrixInverse_Gauss_Jordan(A, A_inv, m+n);
        
        point1[0] = myPointStruct_f->x[i];
        point1[1] = myPointStruct_f->y[i];
        point1[2] = myPointStruct_f->z[i];
        
        for (short j = 0; j < m; j++) {
            point2[0] = myPointStruct_c->x[myPointStruct_c->cloud_index[i_prol][j]];
            point2[1] = myPointStruct_c->y[myPointStruct_c->cloud_index[i_prol][j]];
            point2[2] = myPointStruct_c->z[myPointStruct_c->cloud_index[i_prol][j]];
            temp[j] = calculate_phs_rbf(point1, point2, parameters.phs_degree, 
                                                        parameters.dimension);
        }
        seed_pt[0]= myPointStruct_c->x[i_prol];
        seed_pt[1]= myPointStruct_c->y[i_prol];
        seed_pt[2]= myPointStruct_c->z[i_prol];
        for (short j = 0; j < n; j++) 
            temp[j+m] = pow(point1[0]-seed_pt[0], myPointStruct_c->pow_x[j])
                                * pow(point1[1]-seed_pt[1], myPointStruct_c->pow_y[j])
                                * pow(point1[2]-seed_pt[2], myPointStruct_c->pow_z[j]);
        
        multiply_vector_matrix_columnwise(temp, A_inv, temp1, m+n, m);
        
        for (short j = 0; j < (m); j++)
        {
            myPointStruct_f->prol_mat[i][j] = temp1[j];
        }
    }
    free(temp);
    free(temp1);
    free_matrix(A, mpn);
    free_matrix(A_inv, mpn);
}


#endif