// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign
// Functions used to write the output to files

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "functions.h"

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

int write_vtk(char *gmsh_filename, FieldVariables *field, PointStructure* myPS, int step)
{
    FILE *fp_in, *fp_out;
    char line[256];
    char vtk_filename[256];
    sprintf(vtk_filename, "Solution_%06d.vtk", step);

    int num_nodes = 0, num_elements = 0;
    int i, node_id;
    double x, y, z;

    fp_in = fopen(gmsh_filename, "r");
    if (!fp_in) {
        fprintf(stderr, "Error: Cannot open %s\n", gmsh_filename);
        return -1;
    }

    /* ---------------- READ NODES ---------------- */
    while (fgets(line, sizeof(line), fp_in)) {
        if (strstr(line, "$Nodes")) {
            fscanf(fp_in, "%d", &num_nodes);
            break;
        }
    }

    double *nodes_x = malloc(num_nodes * sizeof(double));
    double *nodes_y = malloc(num_nodes * sizeof(double));
    double *nodes_z = malloc(num_nodes * sizeof(double));

    for (i = 0; i < num_nodes; i++) {
        fscanf(fp_in, "%d %lf %lf %lf", &node_id, &x, &y, &z);
        nodes_x[node_id - 1] = x;
        nodes_y[node_id - 1] = y;
        nodes_z[node_id - 1] = z;
    }

    /* ---------------- READ ELEMENTS ---------------- */
    while (fgets(line, sizeof(line), fp_in)) {
        if (strstr(line, "$Elements")) {
            fscanf(fp_in, "%d", &num_elements);
            break;
        }
    }

    int max_conn = (parameters.dimension == 2) ? 3 : 4;
    int vtk_type  = (parameters.dimension == 2) ? 5 : 10;

    int *conn = malloc(num_elements * max_conn * sizeof(int));
    int cell_count = 0;

    for (i = 0; i < num_elements; i++) {
        int elem_id, elem_type, num_tags;

        fscanf(fp_in, "%d %d %d", &elem_id, &elem_type, &num_tags);

        for (int j = 0; j < num_tags; j++) {
            int tmp;
            fscanf(fp_in, "%d", &tmp);
        }

        if (parameters.dimension == 2 && elem_type == 2) {
            int n1, n2, n3;
            fscanf(fp_in, "%d %d %d", &n1, &n2, &n3);

            conn[cell_count*3+0] = n1-1;
            conn[cell_count*3+1] = n2-1;
            conn[cell_count*3+2] = n3-1;
            cell_count++;
        }
        else if (parameters.dimension == 3 && elem_type == 4) {
            int n1, n2, n3, n4;
            fscanf(fp_in, "%d %d %d %d", &n1, &n2, &n3, &n4);

            conn[cell_count*4+0] = n1-1;
            conn[cell_count*4+1] = n2-1;
            conn[cell_count*4+2] = n3-1;
            conn[cell_count*4+3] = n4-1;
            cell_count++;
        }
        else {
            fgets(line, sizeof(line), fp_in);
        }
    }

    fclose(fp_in);

    /* ---------------- WRITE VTK ---------------- */
    fp_out = fopen(vtk_filename, "w");

    fprintf(fp_out, "# vtk DataFile Version 3.0\n");
    fprintf(fp_out, "Mesh\n");
    fprintf(fp_out, "ASCII\n");
    fprintf(fp_out, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp_out, "POINTS %d double\n", num_nodes);
    for (i = 0; i < num_nodes; i++) {
        fprintf(fp_out, "%lf %lf %lf\n",
                nodes_x[i], nodes_y[i], nodes_z[i]);
    }

    int vtk_cell_size = max_conn + 1;

    fprintf(fp_out, "\nCELLS %d %d\n",
            cell_count,
            cell_count * vtk_cell_size);

    for (i = 0; i < cell_count; i++) {
        fprintf(fp_out, "%d ", max_conn);
        for (int j = 0; j < max_conn; j++) {
            fprintf(fp_out, "%d ", conn[i*max_conn + j]);
        }
        fprintf(fp_out, "\n");
    }

    fprintf(fp_out, "\nCELL_TYPES %d\n", cell_count);
    for (i = 0; i < cell_count; i++) {
        fprintf(fp_out, "%d\n", vtk_type);
    }

    /* ---------------- POINT DATA ---------------- */
    fprintf(fp_out, "\nPOINT_DATA %d\n", num_nodes);

    /* Velocity */
    fprintf(fp_out, "VECTORS velocity double\n");

    int num_corner = myPS->num_corners;

    for (i = 0; i < num_corner; i++) {
        fprintf(fp_out, "0.0 0.0 0.0\n");
    }

    for (i = num_corner; i < num_nodes; i++) {
        int k = myPS[0].rcm_order[i - num_corner];

        if (parameters.dimension == 2) {
            fprintf(fp_out, "%.16e %.16e %.16e\n",
                    field->u[k], field->v[k], 0.0);
        } else {
            fprintf(fp_out, "%.16e %.16e %.16e\n",
                    field->u[k], field->v[k], field->w[k]);
        }
    }

    /* Pressure */
    fprintf(fp_out, "\nSCALARS pressure double 1\n");
    fprintf(fp_out, "LOOKUP_TABLE default\n");

    for (i = 0; i < num_corner; i++) {
        fprintf(fp_out, "0.0\n");
    }

    for (i = num_corner; i < num_nodes; i++) {
        int k = myPS[0].rcm_order[i - num_corner];
        fprintf(fp_out, "%.16e\n", field->p[k]);
    }

    /* Peclet number: |velocity| * (distance to nearest cloud neighbour) / nu, per node */
    fprintf(fp_out, "\nSCALARS peclet double 1\n");
    fprintf(fp_out, "LOOKUP_TABLE default\n");

    for (i = 0; i < num_corner; i++) {
        fprintf(fp_out, "0.0\n");
    }

    int ncp = myPS->num_cloud_points;
    for (i = num_corner; i < num_nodes; i++) {
        int k = myPS[0].rcm_order[i - num_corner];
        int nearest = myPS->cloud_index[k*ncp + 1];   // index 0 in the cloud is always the node itself

        double dx = myPS->x[k] - myPS->x[nearest];
        double dy = myPS->y[k] - myPS->y[nearest];
        double dz = (parameters.dimension == 3) ? (myPS->z[k] - myPS->z[nearest]) : 0.0;
        double h  = sqrt(dx*dx + dy*dy + dz*dz);
        
        // actual vmag
        double vmag = (parameters.dimension == 3)
                        ? sqrt(field->u[k]*field->u[k] + field->v[k]*field->v[k] + field->w[k]*field->w[k])
                        : sqrt(field->u[k]*field->u[k] + field->v[k]*field->v[k]);

        // manually overwrite
        // double vmag = 1;
        double pe = (parameters.nu > 0.0) ? vmag * h / parameters.nu : 0.0;
        fprintf(fp_out, "%.16e\n", pe);
    }    

    fclose(fp_out);

    free(nodes_x);
    free(nodes_y);
    free(nodes_z);
    free(conn);

    printf("VTK written: %d cells\n", cell_count);

    return 0;
}

int read_vtk_restart(char *vtk_filename, FieldVariables *field, PointStructure* myPS)
{
    FILE *fp_in;
    char line[256];
    int num_points = 0;

    fp_in = fopen(vtk_filename, "r");
    if (!fp_in) {
        fprintf(stderr, "Error: Cannot open restart file %s\n", vtk_filename);
        return -1;
    }

    while (fgets(line, sizeof(line), fp_in)) {
        if (strstr(line, "POINT_DATA")) {
            sscanf(line, "POINT_DATA %d", &num_points);
            break;
        }
    }

    int num_corner = myPS->num_corners;
    int expected_points = num_corner + myPS->num_nodes;
    if (num_points != expected_points) {
        fprintf(stderr, "Error: restart file %s has %d points, expected %d -- different mesh?\n",
                vtk_filename, num_points, expected_points);
        fclose(fp_in);
        return -1;
    }

    /* Velocity: VTK point i (after the corner points) is solver node rcm_order[i - num_corner] */
    while (fgets(line, sizeof(line), fp_in)) {
        if (strstr(line, "VECTORS velocity")) break;
    }
    for (int i = 0; i < num_points; i++) {
        double u, v, w;
        if (fscanf(fp_in, "%lf %lf %lf", &u, &v, &w) != 3) {
            fprintf(stderr, "Error reading velocity at point %d in %s\n", i, vtk_filename);
            fclose(fp_in);
            return -1;
        }
        if (i >= num_corner) {
            int k = myPS->rcm_order[i - num_corner];
            field->u[k] = u;
            field->v[k] = v;
            if (parameters.dimension == 3)
                field->w[k] = w;
        }
    }

    /* Pressure */
    while (fgets(line, sizeof(line), fp_in)) {
        if (strstr(line, "SCALARS pressure")) break;
    }
    fgets(line, sizeof(line), fp_in);   // skip the LOOKUP_TABLE line
    for (int i = 0; i < num_points; i++) {
        double p;
        if (fscanf(fp_in, "%lf", &p) != 1) {
            fprintf(stderr, "Error reading pressure at point %d in %s\n", i, vtk_filename);
            fclose(fp_in);
            return -1;
        }
        if (i >= num_corner) {
            int k = myPS->rcm_order[i - num_corner];
            field->p[k] = p;
        }
    }

    fclose(fp_in);
    printf("Restart: loaded u, v, p from %s\n", vtk_filename);
    return 0;
}

static const char* get_poisson_solver_name(int type) {
    switch (type) {
        case 1:  return "Jacobi";
        case 2:  return "Gauss-Seidel";
        case 3:  return "BiCGStab";
        default: return "Unknown";
    }
}

static const char* get_time_scheme_name(int scheme) {
    switch (scheme) {
        case 0:  return "Explicit";
        case 1:  return "Crank-Nicolson";
        default: return "Unknown";
    }
}

static const char* get_algorithm_name(int fractional_step) {
    switch (fractional_step) {
        case 0:  return "Time Implicit";
        case 1:  return "Fractional Step";
        default: return "Unknown Algorithm";
    }
}


void write_solver_data(const PointStructure* point_struct, double steady_state_error, int iteration)
{
    const char* filename = "Solver_data.txt";
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Unable to open output file '%s'\n", filename);
        exit(EXIT_FAILURE);
    }   
    // Header
    fprintf(file, "================================================================================\n");
    fprintf(file, "                          SOLVER LOG & RUN DATA SUMMARY                         \n");
    fprintf(file, "================================================================================\n\n");

    // Section 1: Mesh & Geometry
    fprintf(file, "--- MESH & GEOMETRY INFORMATION ------------------------------------------------\n");
    fprintf(file, "  %-30s : %s\n", "Mesh Filename",       point_struct->mesh_filename);
    fprintf(file, "  %-30s : %d\n", "Number of Nodes",     point_struct->num_nodes);
    fprintf(file, "  %-30s : %d\n", "Number of Corners",   point_struct->num_corners);
    fprintf(file, "  %-30s : %d\n", "Cloud Points Count",  point_struct->num_cloud_points);
    fprintf(file, "  %-30s : %.6e\n", "Avg Point Distance (d_avg)", point_struct->d_avg);
    fprintf(file, "  %-30s : %.6e\n", "Min Point Distance (d_min)", point_struct->d_min);
    fprintf(file, "  %-30s : %.6e\n", "Max Point Distance (d_max)", point_struct->d_max);
    fprintf(file, "\n");

    // Section 2: Physical & Domain Parameters
    fprintf(file, "--- PHYSICAL & DOMAIN PROPERTIES -----------------------------------------------\n");
    fprintf(file, "  %-30s : %dD\n",  "Dimension",            parameters.dimension);
    fprintf(file, "  %-30s : %.6e\n", "Kinematic Viscosity (nu)", parameters.nu);
    fprintf(file, "  %-30s : %.6f\n", "Fluid Density (rho)",      parameters.rho);
    fprintf(file, "\n");

    // Section 3: Numerical Model & Basis Function
    fprintf(file, "--- NUMERICAL MODEL PARAMETERS -------------------------------------------------\n");
    fprintf(file, "  %-30s : %s\n",  "Algorithm Type",       get_algorithm_name(parameters.fractional_step));
    fprintf(file, "  %-30s : %d\n",  "Polynomial Degree",    parameters.poly_degree);
    fprintf(file, "  %-30s : %d\n",  "PHS Degree",           parameters.phs_degree);
    fprintf(file, "  %-30s : %.6f\n", "Courant Number",       parameters.courant_number);
    fprintf(file, "  %-30s : %d\n",  "Number of Levels",     parameters.num_levels);
    fprintf(file, "  %-30s : %d\n",  "Test Parameter",       parameters.test);
    fprintf(file, "  %-30s : %d\n",  "Restart Flag",         parameters.restart);
    fprintf(file, "\n");

    // Section 4: Time Integration Scheme
    fprintf(file, "--- TIME DISCRETIZATION SCHEME ------------------------------------------------\n");
    fprintf(file, "  %-30s : %s\n",   "Time Scheme",          get_time_scheme_name(parameters.time_scheme));
    fprintf(file, "  %-30s : %.4f\n", "Theta Parameter",      parameters.theta);
    fprintf(file, "  %-30s : %.6e\n", "Time Step Size (dt)",  parameters.dt);
    fprintf(file, "  %-30s : %d\n",   "Total Time Steps",    parameters.num_time_steps);
    fprintf(file, "\n");

    // Section 5: Poisson Solver Settings
    fprintf(file, "--- POISSON SOLVER SETTINGS ----------------------------------------------------\n");
    fprintf(file, "  %-30s : %s\n",   "Poisson Solver Type",  get_poisson_solver_name(parameters.poisson_solver_type));
    fprintf(file, "  %-30s : %.6e\n", "Steady State Tolerance", parameters.steady_state_tolerance);
    fprintf(file, "  %-30s : %.6e\n", "Poisson Tolerance",    parameters.poisson_solver_tolerance);
    fprintf(file, "  %-30s : %.4f\n", "SOR Relaxation (omega)", parameters.omega);
    fprintf(file, "  %-30s : %d\n",   "Max SOR Iterations",   parameters.num_relax);
    fprintf(file, "\n");

    // Section 6: Execution Profiling
    fprintf(file, "--- EXECUTION TIMING PROFILE (seconds) -----------------------------------------\n");
    fprintf(file, "  %-30s : %.4f s\n", "Read Grids & Flow Params", timer.time_read);
    fprintf(file, "  %-30s : %.4f s\n", "Create Derivative Matrices", timer.time_derivatives);
    fprintf(file, "  %-30s : %.4f s\n", "Initialization Time",     timer.time_initialisation);
    fprintf(file, "  %-30s : %.4f s\n", "GPU Copy Time",          timer.time_copy_to_gpu);
    fprintf(file, "  %-30s : %.4f s\n", "Total Solver Loop Time",  timer.time_total);
    fprintf(file, "\n");

    // Section 7: Current Convergence Snapshot
    fprintf(file, "--- CONVERGENCE SNAPSHOT -------------------------------------------------------\n");
    fprintf(file, "  %-30s : %d\n",   "Current Iteration",    iteration);
    fprintf(file, "  %-30s : %.6e\n", "Steady State Error Residual", steady_state_error);
    fprintf(file, "================================================================================\n");

    fclose(file);
}