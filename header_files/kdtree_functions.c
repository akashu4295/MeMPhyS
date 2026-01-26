#include "structures.h"
#include "kdtree.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <stdio.h>

//////////////////////////////////////////////////////////////////////////////////
// Function definitions
//////////////////////////////////////////////////////////////////////////////////
// static double dist_sq( double *a1, double *a2 ) { // used in kdtree sorting
//     double dist_sq = 0, diff;
//     for (int dims = 0; dims < 3; dims++) {
//         diff = (a1[dims] - a2[dims]);
//         dist_sq += diff*diff;
//     }
//     return dist_sq;
// }

// void* create_kdtree(PointStructure* myPointStruct){
//     void* ptree;
//     ptree = kd_create(3); // create a k-d tree for 3-dimensional points 
//     for(int i=0; i<myPointStruct->num_nodes; i++ ) {
//         if(myPointStruct->corner_tag[i] == false)
//             assert(kd_insert3(ptree, myPointStruct->x[i], myPointStruct->y[i], myPointStruct->z[i], &myPointStruct->point_index[i]) == 0);
//     }
//     return ptree;
// }

// void* create_kdtree_without_boundarynodes(PointStructure* myPointStruct){
//     void* ptree;
//     ptree = kd_create(3); // create a k-d tree for 3-dimensional points 
//     for(int i=0; i<myPointStruct->num_nodes; i++ ) {
//         if(myPointStruct->boundary_tag[i] == false)
//             assert(kd_insert3(ptree, myPointStruct->x[i], myPointStruct->y[i], myPointStruct->z[i], &myPointStruct->point_index[i]) == 0);
//     }
//     return ptree;
// }

// void free_kdtree(void* ptree){
//     kd_free( ptree );
// }

// int* find_neighbours(double* p, void* ptree, double radius, int num_cloud_points){
//     double pos[3], dist;    // The position of the nearest neighbor and the distance from the search point

//     // Arrays to store the distances and indices of the nearest neighbors for sorting purposes and counting the number of points
//     double* distance = (double*)malloc(num_cloud_points * sizeof(double*));  // Array to store the distances of the nearest neighbors
//     int* ind = (int*)malloc(num_cloud_points * sizeof(int*));  // Array to store the indices of the nearest neighbors
    
//     struct kdres *presults; // Kdtree structure to store results of search
    
//     int* pch; // Pointer to the index of the nearest neighbor
//     // Initialise the distance and index arrays
//     for (int i =0; i<num_cloud_points ; i++) {
//         distance[i] = 100.0;
//         ind[i] = 0;
//     }

//     // find points closest to the point pt and within distance radius
//     presults = kd_nearest_range( ptree, p, radius ); 

//     // browse through the results of the search
//     while (!kd_res_end( presults ))
//     {   
//         pch = (int*)kd_res_item( presults, pos ); // get the neighbour
//         dist = sqrt( dist_sq( p, pos) ); // 3 is the dimension of the problem, we will keep it 3 irrespective of 2d or 3d problem

//         for (int j = 0; j < num_cloud_points; j++) {
//             if (dist < distance[j]) {
//                 for (int k = num_cloud_points-1; k > j; k--) {
//                     distance[k] = distance[k-1];
//                     ind[k] = ind[k-1];
//                 }
//                 distance[j] = dist;
//                 ind[j] = *pch;
//                 break;
//             }
//         }
//         kd_res_next( presults ); // move to the next neighbour
//     }
//     if (distance[num_cloud_points-1]==100) {  // If the last distance is 100, increase the radius and search again
//         radius = 1.2*radius;
//         kd_res_free( presults );
//         free(distance);
//         free(ind);
//         return find_neighbours(p, ptree, radius, num_cloud_points);
//     }
//     else {
//         free(distance);
//         kd_res_free( presults );
//         return ind;
//     }
// }

// void find_cloud_index(PointStructure* myPointStruct1){
//     // Create a kdtree for the cloud points
//     int n = myPointStruct1->num_cloud_points;
    
//     myPointStruct1->cloud_index = (int*)malloc(myPointStruct1->num_nodes * n * sizeof(int)); // Allocate memory for the cloud index

//     void* ptree = create_kdtree(myPointStruct1);
//     double radius = (myPointStruct1->d_avg) * n ; // Initial radius for the search
    
//     double pt[3]; // Array to store the coordinates of the search point
//     for (int i = 0; i < myPointStruct1->num_nodes; i++) {
//         if (myPointStruct1->boundary_tag[i])
//             continue;
//         int* neighbours; // Array to store the indices of the nearest neighbours
//         pt[0] = myPointStruct1->x[i]; pt[1] = myPointStruct1->y[i]; pt[2] = myPointStruct1->z[i];
//         neighbours = find_neighbours(pt, ptree, radius, n);
//         for (int j = 0; j < n; j++) 
//             myPointStruct1->cloud_index[i*n + j] = neighbours[j];
//         free(neighbours);
//     }
//     free_kdtree(ptree);

//     void* ptree2 = create_kdtree_without_boundarynodes(myPointStruct1);
    
//     for (int i = 0; i < myPointStruct1->num_nodes; i++) {
//         if (!myPointStruct1->boundary_tag[i])
//             continue;
//         int* neighbours; // Array to store the indices of the nearest neighbours
//         pt[0] = myPointStruct1->x[i]; pt[1] = myPointStruct1->y[i]; pt[2] = myPointStruct1->z[i];
//         neighbours = find_neighbours(pt, ptree2, radius, n);
//         myPointStruct1->cloud_index[i*n] = i;
//         for (int j = 1; j < n; j++) 
//             myPointStruct1->cloud_index[i*n +j] = neighbours[j-1];
//         free(neighbours);
//     }
//     free_kdtree(ptree2);
// }

// int* find_nearest_point(PointStructure* myPointStruct1, PointStructure* myPointStruct2, int num_cloud_points){
//     // Create a kdtree for the cloud points
//     void* ptree = create_kdtree(myPointStruct2);
  
//     double radius = (myPointStruct2->d_avg) * 10; // Initial radius for the search
//     int* neighbour; // Array to store the indices of the nearest neighbours
//     int* temp; // Temporary array to store the indices of the nearest neighbours
//     double pt[3]; // Array to store the coordinates of the search point
//     neighbour = (int*)malloc(myPointStruct1->num_nodes * sizeof(int*));
//     temp = (int*)malloc(10 * sizeof(int*));
//     for (int i = 0; i < myPointStruct1->num_nodes; i++) {
//         if (myPointStruct1->corner_tag[i] == false) {
//             pt[0] = myPointStruct1->x[i]; pt[1] = myPointStruct1->y[i]; pt[2] = myPointStruct1->z[i];
//             temp = find_neighbours(pt, ptree, radius, 3);
//             neighbour[i] = temp[0];
//         }
//         else
//             neighbour[i] = 0;
//     }
//     free_kdtree(ptree);
//     free(temp);
//     return neighbour;
// }


// /* -------------------------------------------------- */
// /* Utility                                            */
// /* -------------------------------------------------- */

static double dist_sq(double *a1, double *a2)
{
    double d0 = a1[0] - a2[0];
    double d1 = a1[1] - a2[1];
    double d2 = a1[2] - a2[2];
    return d0*d0 + d1*d1 + d2*d2;
}

/* -------------------------------------------------- */
/* KD-tree creation                                   */
/* -------------------------------------------------- */

void* create_kdtree(PointStructure* ps)
{
    void* ptree = kd_create(3);
    for (int i = 0; i < ps->num_nodes; i++) {
        if (!ps->corner_tag[i]) {
            assert(
                kd_insert3(ptree,
                           ps->x[i], ps->y[i], ps->z[i],
                           &ps->point_index[i]) == 0
            );
        }
    }
    return ptree;
}

void* create_kdtree_without_boundarynodes(PointStructure* ps)
{
    void* ptree = kd_create(3);
    for (int i = 0; i < ps->num_nodes; i++) {
        if (!ps->boundary_tag[i]) {
            assert(
                kd_insert3(ptree,
                           ps->x[i], ps->y[i], ps->z[i],
                           &ps->point_index[i]) == 0
            );
        }
    }
    return ptree;
}

void free_kdtree(void* ptree)
{
    kd_free(ptree);
}

/* -------------------------------------------------- */
/* Neighbour search                                   */
/* -------------------------------------------------- */

int* find_neighbours(double* p, void* ptree, double radius, int num_cloud_points)
{
    const double INF = DBL_MAX;
    const int MAX_EXPANDS = 50;

    double pos[3], dist;
    struct kdres *presults;

    double* distance = malloc(num_cloud_points * sizeof(double));
    int* ind = malloc(num_cloud_points * sizeof(int));

    if (!distance || !ind) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < num_cloud_points; i++) {
        distance[i] = INF;
        ind[i] = -1;
    }

    for (int expand = 0; expand < MAX_EXPANDS; expand++) {

        presults = kd_nearest_range(ptree, p, radius);

        while (!kd_res_end(presults)) {
            int* pch = (int*)kd_res_item(presults, pos);
            dist = sqrt(dist_sq(p, pos));

            for (int j = 0; j < num_cloud_points; j++) {
                if (dist < distance[j]) {
                    for (int k = num_cloud_points - 1; k > j; k--) {
                        distance[k] = distance[k - 1];
                        ind[k] = ind[k - 1];
                    }
                    distance[j] = dist;
                    ind[j] = *pch;
                    break;
                }
            }
            kd_res_next(presults);
        }

        kd_res_free(presults);

        if (distance[num_cloud_points - 1] < INF) {
            free(distance);
            return ind;
        }

        radius *= 1.5;
    }

    fprintf(stderr, "KD-tree search failed: insufficient neighbors\n");
    free(distance);
    free(ind);
    return NULL;
}

/* -------------------------------------------------- */
/* Cloud construction                                 */
/* -------------------------------------------------- */

void* create_kdtree_no_corners(PointStructure* ps)
{
    void* ptree = kd_create(3);
    for (int i = 0; i < ps->num_nodes; i++) {
        if (!ps->corner_tag[i]) {
            kd_insert3(ptree, ps->x[i], ps->y[i], ps->z[i],
                       &ps->point_index[i]);
        }
    }
    return ptree;
}

void* create_kdtree_interior_only(PointStructure* ps)
{
    void* ptree = kd_create(3);
    for (int i = 0; i < ps->num_nodes; i++) {
        if (!ps->boundary_tag[i] && !ps->corner_tag[i]) {
            kd_insert3(ptree, ps->x[i], ps->y[i], ps->z[i],
                       &ps->point_index[i]);
        }
    }
    return ptree;
}

void find_cloud_index(PointStructure* ps)
{
    int n = ps->num_cloud_points;
    int N = ps->num_nodes;

    ps->cloud_index = malloc(N * n * sizeof(int));
    if (!ps->cloud_index) abort();

    double radius = ps->d_avg * n;
    double pt[3];

    /* -------- Interior nodes (exclude corners) -------- */
    void* ptree_all = create_kdtree_no_corners(ps);

    for (int i = 0; i < N; i++) {
        if (ps->corner_tag[i]) continue;
        if (ps->boundary_tag[i]) continue;

        pt[0] = ps->x[i];
        pt[1] = ps->y[i];
        pt[2] = ps->z[i];

        int* neigh = find_neighbours(pt, ptree_all, radius, n);

        for (int j = 0; j < n; j++)
            ps->cloud_index[i*n + j] = neigh[j];

        free(neigh);
    }
    free_kdtree(ptree_all);

    /* -------- Boundary nodes → interior-only cloud -------- */
    void* ptree_int = create_kdtree_interior_only(ps);

    for (int i = 0; i < ps->num_nodes; i++) {
        if (ps->corner_tag[i]) continue;
        if (!ps->boundary_tag[i]) continue;

        pt[0] = ps->x[i];
        pt[1] = ps->y[i];
        pt[2] = ps->z[i];

        int* neigh = find_neighbours(pt, ptree_int, radius, n - 1);

        ps->cloud_index[i*n] = i;  /* self */
        for (int j = 1; j < n; j++)
            ps->cloud_index[i*n + j] = neigh[j - 1];

        free(neigh);
    }
    free_kdtree(ptree_int);
}

/* -------------------------------------------------- */
/* Nearest-point mapping                              */
/* -------------------------------------------------- */

int* find_nearest_point(PointStructure* ps1,
                        PointStructure* ps2,
                        int num_cloud_points)
{
    int* neighbour = malloc(ps1->num_nodes * sizeof(int));
    if (!neighbour) abort();

    void* ptree = create_kdtree(ps2);
    double radius = ps2->d_avg * 10.0;
    double pt[3];

    for (int i = 0; i < ps1->num_nodes; i++) {
        if (!ps1->corner_tag[i]) {
            pt[0] = ps1->x[i];
            pt[1] = ps1->y[i];
            pt[2] = ps1->z[i];

            int* tmp = find_neighbours(pt, ptree, radius, 3);
            neighbour[i] = tmp[0];
            free(tmp);
        } else {
            neighbour[i] = -1;
        }
    }

    free_kdtree(ptree);
    return neighbour;
}