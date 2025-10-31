#include "structures.h"
#include "kdtree.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>

//////////////////////////////////////////////////////////////////////////////////
// Function definitions
//////////////////////////////////////////////////////////////////////////////////
static double dist_sq( double *a1, double *a2 ) { // used in kdtree sorting
    double dist_sq = 0, diff;
    for (int dims = 0; dims < 3; dims++) {
        diff = (a1[dims] - a2[dims]);
        dist_sq += diff*diff;
    }
    return dist_sq;
}

void* create_kdtree(PointStructure* myPointStruct){
    void* ptree;
    ptree = kd_create(3); // create a k-d tree for 3-dimensional points 
    for(int i=0; i<myPointStruct->num_nodes; i++ ) {
        if(myPointStruct->corner_tag[i] == false)
            assert(kd_insert3(ptree, myPointStruct->x[i], myPointStruct->y[i], myPointStruct->z[i], &myPointStruct->point_index[i]) == 0);
    }
    return ptree;
}

void* create_kdtree_without_boundarynodes(PointStructure* myPointStruct){
    void* ptree;
    ptree = kd_create(3); // create a k-d tree for 3-dimensional points 
    for(int i=0; i<myPointStruct->num_nodes; i++ ) {
        if(myPointStruct->boundary_tag[i] == false)
            assert(kd_insert3(ptree, myPointStruct->x[i], myPointStruct->y[i], myPointStruct->z[i], &myPointStruct->point_index[i]) == 0);
    }
    return ptree;
}

void free_kdtree(void* ptree){
    kd_free( ptree );
}

int* find_neighbours(double* p, void* ptree, double radius, int num_cloud_points){
    double pos[3], dist;    // The position of the nearest neighbor and the distance from the search point

    // Arrays to store the distances and indices of the nearest neighbors for sorting purposes and counting the number of points
    double* distance = (double*)malloc(num_cloud_points * sizeof(double*));  // Array to store the distances of the nearest neighbors
    int* ind = (int*)malloc(num_cloud_points * sizeof(int*));  // Array to store the indices of the nearest neighbors
    
    struct kdres *presults; // Kdtree structure to store results of search
    
    int* pch; // Pointer to the index of the nearest neighbor
    // Initialise the distance and index arrays
    for (int i =0; i<num_cloud_points ; i++) {
        distance[i] = 100.0;
        ind[i] = 0;
    }

    // find points closest to the point pt and within distance radius
    presults = kd_nearest_range( ptree, p, radius ); 

    // browse through the results of the search
    while (!kd_res_end( presults ))
    {   
        pch = (int*)kd_res_item( presults, pos ); // get the neighbour
        dist = sqrt( dist_sq( p, pos) ); // 3 is the dimension of the problem, we will keep it 3 irrespective of 2d or 3d problem

        for (int j = 0; j < num_cloud_points; j++) {
            if (dist < distance[j]) {
                for (int k = num_cloud_points-1; k > j; k--) {
                    distance[k] = distance[k-1];
                    ind[k] = ind[k-1];
                }
                distance[j] = dist;
                ind[j] = *pch;
                break;
            }
        }
        kd_res_next( presults ); // move to the next neighbour
    }
    if (distance[num_cloud_points-1]==100) {  // If the last distance is 100, increase the radius and search again
        radius = 1.2*radius;
        kd_res_free( presults );
        free(distance);
        free(ind);
        return find_neighbours(p, ptree, radius, num_cloud_points);
    }
    else {
        free(distance);
        kd_res_free( presults );
        return ind;
    }
}

void find_cloud_index(PointStructure* myPointStruct1){
    // Create a kdtree for the cloud points
    int n = myPointStruct1->num_cloud_points;
    
    myPointStruct1->cloud_index = (int*)malloc(myPointStruct1->num_nodes * n * sizeof(int)); // Allocate memory for the cloud index

    void* ptree = create_kdtree(myPointStruct1);
    double radius = (myPointStruct1->d_avg) * n ; // Initial radius for the search
    
    double pt[3]; // Array to store the coordinates of the search point
    for (int i = myPointStruct1->num_boundary_nodes; i < myPointStruct1->num_nodes; i++) {
        int* neighbours; // Array to store the indices of the nearest neighbours
        pt[0] = myPointStruct1->x[i]; pt[1] = myPointStruct1->y[i]; pt[2] = myPointStruct1->z[i];
        neighbours = find_neighbours(pt, ptree, radius, n);
        for (int j = 0; j < n; j++) 
            myPointStruct1->cloud_index[i*n + j] = neighbours[j];
        free(neighbours);
    }
    free_kdtree(ptree);

    void* ptree2 = create_kdtree_without_boundarynodes(myPointStruct1);
    
    for (int i = 0; i < myPointStruct1->num_boundary_nodes; i++) {
        int* neighbours; // Array to store the indices of the nearest neighbours
        pt[0] = myPointStruct1->x[i]; pt[1] = myPointStruct1->y[i]; pt[2] = myPointStruct1->z[i];
        neighbours = find_neighbours(pt, ptree2, radius, n);
        myPointStruct1->cloud_index[i*n] = i;
        for (int j = 1; j < n; j++) 
            myPointStruct1->cloud_index[i*n +j] = neighbours[j-1];
        free(neighbours);
    }
    free_kdtree(ptree2);
}

int* find_nearest_point(PointStructure* myPointStruct1, PointStructure* myPointStruct2, int num_cloud_points){
    // Create a kdtree for the cloud points
    void* ptree = create_kdtree(myPointStruct2);
  
    double radius = (myPointStruct2->d_avg) * 10; // Initial radius for the search
    int* neighbour; // Array to store the indices of the nearest neighbours
    int* temp; // Temporary array to store the indices of the nearest neighbours
    double pt[3]; // Array to store the coordinates of the search point
    neighbour = (int*)malloc(myPointStruct1->num_nodes * sizeof(int*));
    temp = (int*)malloc(10 * sizeof(int*));
    for (int i = 0; i < myPointStruct1->num_nodes; i++) {
        if (myPointStruct1->corner_tag[i] == false) {
            pt[0] = myPointStruct1->x[i]; pt[1] = myPointStruct1->y[i]; pt[2] = myPointStruct1->z[i];
            temp = find_neighbours(pt, ptree, radius, 3);
            neighbour[i] = temp[0];
        }
        else
            neighbour[i] = 0;
    }
    free_kdtree(ptree);
    free(temp);
    return neighbour;
}