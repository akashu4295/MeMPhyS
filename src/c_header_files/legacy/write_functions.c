
// Helper function to find nearest neighbor and interpolate
void interpolate_field_data(double x_target, double y_target, double z_target,
                            PointStructure *myPS, FieldVariables *field,
                            double *u_out, double *v_out, double *w_out, double *p_out)
{
    // Find nearest neighbor (simple approach)
    // For better results, could use inverse distance weighting with multiple neighbors
    double min_dist = 1e30;
    int nearest_idx = 0;
    
    for (int i = 0; i < myPS->num_nodes; i++) {
        double dx = myPS->x[i] - x_target;
        double dy = myPS->y[i] - y_target;
        double dz = myPS->z[i] - z_target;
        double dist = sqrt(dx*dx + dy*dy + dz*dz);
        
        if (dist < min_dist) {
            min_dist = dist;
            nearest_idx = i;
        }
    }
    
    // Return field values at nearest point
    *u_out = field->u[nearest_idx];
    *v_out = field->v[nearest_idx];
    *w_out = field->w[nearest_idx];
    *p_out = field->p[nearest_idx];
}

// Alternative: Inverse distance weighted interpolation (more accurate)
void interpolate_field_data_idw(double x_target, double y_target, double z_target,
                                PointStructure *myPS, FieldVariables *field,
                                double *u_out, double *v_out, double *w_out, double *p_out,
                                int num_neighbors)
{
    // Find k nearest neighbors
    double *distances = (double*)malloc(myPS->num_nodes * sizeof(double));
    int *indices = (int*)malloc(myPS->num_nodes * sizeof(int));
    
    // Calculate all distances
    for (int i = 0; i < myPS->num_nodes; i++) {
        double dx = myPS->x[i] - x_target;
        double dy = myPS->y[i] - y_target;
        double dz = myPS->z[i] - z_target;
        distances[i] = sqrt(dx*dx + dy*dy + dz*dz);
        indices[i] = i;
    }
    
    // Simple selection sort for k nearest (could use heap for efficiency)
    for (int i = 0; i < num_neighbors && i < myPS->num_nodes; i++) {
        for (int j = i + 1; j < myPS->num_nodes; j++) {
            if (distances[j] < distances[i]) {
                double temp_d = distances[i];
                distances[i] = distances[j];
                distances[j] = temp_d;
                
                int temp_i = indices[i];
                indices[i] = indices[j];
                indices[j] = temp_i;
            }
        }
    }
    
    // Inverse distance weighting
    double sum_weights = 0.0;
    double u_sum = 0.0, v_sum = 0.0, w_sum = 0.0, p_sum = 0.0;
    
    int k = (num_neighbors < myPS->num_nodes) ? num_neighbors : myPS->num_nodes;
    
    for (int i = 0; i < k; i++) {
        double dist = distances[i];
        if (dist < 1e-10) {
            // Target point coincides with a data point
            *u_out = field->u[indices[i]];
            *v_out = field->v[indices[i]];
            *w_out = field->w[indices[i]];
            *p_out = field->p[indices[i]];
            free(distances);
            free(indices);
            return;
        }
        
        double weight = 1.0 / (dist * dist);
        sum_weights += weight;
        
        u_sum += weight * field->u[indices[i]];
        v_sum += weight * field->v[indices[i]];
        w_sum += weight * field->w[indices[i]];
        p_sum += weight * field->p[indices[i]];
    }
    
    *u_out = u_sum / sum_weights;
    *v_out = v_sum / sum_weights;
    *w_out = w_sum / sum_weights;
    *p_out = p_sum / sum_weights;
    
    free(distances);
    free(indices);
}


int write_vtk_test(char *gmsh_filename, FieldVariables *field, PointStructure* myPS)
{
    FILE *fp_in, *fp_out;
    char line[256];
    char vtk_filename[256] = "Solution_interpolated.vtk";

    int num_nodes = 0, num_elements = 0;
    int i, node_id;
    double x, y, z;

    fp_in = fopen(gmsh_filename, "r");
    if (!fp_in) {
        fprintf(stderr, "Error: Cannot open %s\n", gmsh_filename);
        return -1;
    }

    /* ---------------- 1. READ NODES FROM GMSH ---------------- */
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

    /* ---------------- 2. READ ELEMENTS FROM GMSH ---------------- */
    while (fgets(line, sizeof(line), fp_in)) {
        if (strstr(line, "$Elements")) {
            fscanf(fp_in, "%d", &num_elements);
            break;
        }
    }

    int max_conn = (parameters.dimension == 2) ? 3 : 4;   // Triangle vs Tet
    int vtk_type  = (parameters.dimension == 2) ? 5 : 10; // VTK_TRIANGLE vs VTK_TETRA
    int *conn = malloc(num_elements * max_conn * sizeof(int));
    int cell_count = 0;

    for (i = 0; i < num_elements; i++) {
        int elem_id, elem_type, num_tags;
        fscanf(fp_in, "%d %d %d", &elem_id, &elem_type, &num_tags);

        for (int j = 0; j < num_tags; j++) {
            int tmp; fscanf(fp_in, "%d", &tmp);
        }

        if (parameters.dimension == 2 && elem_type == 2) {
            fscanf(fp_in, "%d %d %d", &conn[cell_count*3+0], &conn[cell_count*3+1], &conn[cell_count*3+2]);
            conn[cell_count*3+0] -= 1; conn[cell_count*3+1] -= 1; conn[cell_count*3+2] -= 1;
            cell_count++;
        }
        else if (parameters.dimension == 3 && elem_type == 4) {
            fscanf(fp_in, "%d %d %d %d", &conn[cell_count*4+0], &conn[cell_count*4+1], &conn[cell_count*4+2], &conn[cell_count*4+3]);
            conn[cell_count*4+0] -= 1; conn[cell_count*4+1] -= 1; conn[cell_count*4+2] -= 1; conn[cell_count*4+3] -= 1;
            cell_count++;
        }
        else {
            fgets(line, sizeof(line), fp_in);
        }
    }
    fclose(fp_in);

    /* ---------------- 3. INTERPOLATE DATA FOR ALL NODES ---------------- */
    printf("Interpolating data for %d nodes (including corners)...\n", num_nodes);
    double *u_interp = malloc(num_nodes * sizeof(double));
    double *v_interp = malloc(num_nodes * sizeof(double));
    double *w_interp = malloc(num_nodes * sizeof(double));
    double *p_interp = malloc(num_nodes * sizeof(double));

    for (i = 0; i < num_nodes; i++) {
        // Using IDW for better quality at corners. Set neighbors to 4 or 8.
        interpolate_field_data_idw(nodes_x[i], nodes_y[i], nodes_z[i], 
                                   myPS, field, 
                                   &u_interp[i], &v_interp[i], &w_interp[i], &p_interp[i], 8);
    }

    /* ---------------- 4. WRITE CONSOLIDATED VTK ---------------- */
    fp_out = fopen(vtk_filename, "w");
    fprintf(fp_out, "# vtk DataFile Version 3.0\nInterpolated Solution\nASCII\nDATASET UNSTRUCTURED_GRID\n");

    // Write Points
    fprintf(fp_out, "POINTS %d double\n", num_nodes);
    for (i = 0; i < num_nodes; i++) {
        fprintf(fp_out, "%lf %lf %lf\n", nodes_x[i], nodes_y[i], nodes_z[i]);
    }

    // Write Cells
    fprintf(fp_out, "\nCELLS %d %d\n", cell_count, cell_count * (max_conn + 1));
    for (i = 0; i < cell_count; i++) {
        fprintf(fp_out, "%d ", max_conn);
        for (int j = 0; j < max_conn; j++) {
            fprintf(fp_out, "%d ", conn[i*max_conn + j]);
        }
        fprintf(fp_out, "\n");
    }

    fprintf(fp_out, "\nCELL_TYPES %d\n", cell_count);
    for (i = 0; i < cell_count; i++) fprintf(fp_out, "%d\n", vtk_type);

    /* ---------------- 5. WRITE INTERPOLATED POINT DATA ---------------- */
    fprintf(fp_out, "\nPOINT_DATA %d\n", num_nodes);

    // Velocity Vectors
    fprintf(fp_out, "VECTORS velocity double\n");
    for (i = 0; i < num_nodes; i++) {
        fprintf(fp_out, "%.16e %.16e %.16e\n", u_interp[i], v_interp[i], (parameters.dimension == 2) ? 0.0 : w_interp[i]);
    }

    // Pressure Scalars
    fprintf(fp_out, "\nSCALARS pressure double 1\nLOOKUP_TABLE default\n");
    for (i = 0; i < num_nodes; i++) {
        fprintf(fp_out, "%.16e\n", p_interp[i]);
    }

    fclose(fp_out);

    // Cleanup
    free(nodes_x); free(nodes_y); free(nodes_z); free(conn);
    free(u_interp); free(v_interp); free(w_interp); free(p_interp);

    printf("VTK written with interpolated data: %s\n", vtk_filename);
    return 0;
}

// // Function to convert Gmsh mesh to VTK with field data using interpolation
// int write_vtk_test(char *gmsh_filename, FieldVariables *field, PointStructure* myPS) 
// {
//     FILE *fp_in, *fp_out;
//     char line[256];
//     char vtk_filename[256] = "Solution_test.vtk";
//     int num_nodes = 0, num_elements = 0;
//     int i, node_id;
//     double x, y, z;
    
//     // Open Gmsh file
//     fp_in = fopen(gmsh_filename, "r");
//     if (!fp_in) {
//         fprintf(stderr, "Error: Cannot open file %s\n", gmsh_filename);
//         return -1;
//     }

//     // First pass: count nodes and elements
//     while (fgets(line, sizeof(line), fp_in)) {
//         if (strstr(line, "$Nodes")) {
//             if (fscanf(fp_in, "%d", &num_nodes) != 1) {
//                 fprintf(stderr, "Error reading number of nodes\n");
//                 fclose(fp_in);
//                 return -1;
//             }
//             break;
//         }
//     }
    
//     // Allocate memory for node coordinates
//     double *nodes_x = (double*)malloc(num_nodes * sizeof(double));
//     double *nodes_y = (double*)malloc(num_nodes * sizeof(double));
//     double *nodes_z = (double*)malloc(num_nodes * sizeof(double));
    
//     if (!nodes_x || !nodes_y || !nodes_z) {
//         fprintf(stderr, "Error: Memory allocation failed\n");
//         fclose(fp_in);
//         return -1;
//     }
    
//     // Read node coordinates
//     for (i = 0; i < num_nodes; i++) {
//         if (fscanf(fp_in, "%d %lf %lf %lf", &node_id, &x, &y, &z) != 4) {
//             fprintf(stderr, "Error reading node %d\n", i);
//             free(nodes_x); free(nodes_y); free(nodes_z);
//             fclose(fp_in);
//             return -1;
//         }
//         nodes_x[node_id - 1] = x;
//         nodes_y[node_id - 1] = y;
//         nodes_z[node_id - 1] = z;
//     }
    
//     // Find $Elements section
//     while (fgets(line, sizeof(line), fp_in)) {
//         if (strstr(line, "$Elements")) {
//             if (fscanf(fp_in, "%d", &num_elements) != 1) {
//                 fprintf(stderr, "Error reading number of elements\n");
//                 free(nodes_x); free(nodes_y); free(nodes_z);
//                 fclose(fp_in);
//                 return -1;
//             }
//             break;
//         }
//     }
    
//     // Count triangular elements (type 2)
//     int num_triangles = 0;
//     int *triangle_conn = (int*)malloc(num_elements * 3 * sizeof(int));
    
//     for (i = 0; i < num_elements; i++) {
//         int elem_id, elem_type, num_tags, tag1, n1, n2, n3;
//         if (fscanf(fp_in, "%d %d %d", &elem_id, &elem_type, &num_tags) != 3) {
//             fprintf(stderr, "Error reading element %d\n", i);
//             free(nodes_x); free(nodes_y); free(nodes_z); free(triangle_conn);
//             fclose(fp_in);
//             return -1;
//         }
        
//         // Skip tags
//         for (int j = 0; j < num_tags; j++) {
//             fscanf(fp_in, "%d", &tag1);
//         }
        
//         if (elem_type == 2) {  // Triangle
//             if (fscanf(fp_in, "%d %d %d", &n1, &n2, &n3) != 3) {
//                 fprintf(stderr, "Error reading triangle nodes\n");
//                 free(nodes_x); free(nodes_y); free(nodes_z); free(triangle_conn);
//                 fclose(fp_in);
//                 return -1;
//             }
//             triangle_conn[num_triangles * 3 + 0] = n1 - 1;
//             triangle_conn[num_triangles * 3 + 1] = n2 - 1;
//             triangle_conn[num_triangles * 3 + 2] = n3 - 1;
//             num_triangles++;
//         } else {
//             // Skip other element types
//             fgets(line, sizeof(line), fp_in);
//         }
//     }
    
//     fclose(fp_in);
    
//     // Allocate arrays for interpolated field values at mesh nodes
//     double *u_nodes = (double*)malloc(num_nodes * sizeof(double));
//     double *v_nodes = (double*)malloc(num_nodes * sizeof(double));
//     double *w_nodes = (double*)malloc(num_nodes * sizeof(double));
//     double *p_nodes = (double*)malloc(num_nodes * sizeof(double));
    
//     if (!u_nodes || !v_nodes || !w_nodes || !p_nodes) {
//         fprintf(stderr, "Error: Memory allocation failed for field arrays\n");
//         free(nodes_x); free(nodes_y); free(nodes_z); free(triangle_conn);
//         return -1;
//     }
    
//     // Interpolate field values at each mesh node
//     printf("Interpolating field data at mesh nodes...\n");
//     for (i = 0; i < num_nodes; i++) {
//         // Use nearest neighbor interpolation (faster)
//         interpolate_field_data(nodes_x[i], nodes_y[i], nodes_z[i],
//                               myPS, field,
//                               &u_nodes[i], &v_nodes[i], &w_nodes[i], &p_nodes[i]);
        
//         // Alternative: Use IDW with 4 nearest neighbors (more accurate but slower)
//         // interpolate_field_data_idw(nodes_x[i], nodes_y[i], nodes_z[i],
//         //                           myPS, field,
//         //                           &u_nodes[i], &v_nodes[i], &w_nodes[i], &p_nodes[i], 4);
//     }
    
//     // Write VTK file
//     fp_out = fopen(vtk_filename, "w");
//     if (!fp_out) {
//         fprintf(stderr, "Error: Cannot create file %s\n", vtk_filename);
//         free(nodes_x); free(nodes_y); free(nodes_z); free(triangle_conn);
//         free(u_nodes); free(v_nodes); free(w_nodes); free(p_nodes);
//         return -1;
//     }
    
//     // Write VTK header
//     fprintf(fp_out, "# vtk DataFile Version 3.0\n");
//     fprintf(fp_out, "Gmsh mesh with interpolated velocity and pressure fields\n");
//     fprintf(fp_out, "ASCII\n");
//     fprintf(fp_out, "DATASET UNSTRUCTURED_GRID\n");
    
//     // Write points
//     fprintf(fp_out, "POINTS %d double\n", num_nodes);
//     for (i = 0; i < num_nodes; i++) {
//         fprintf(fp_out, "%.16e %.16e %.16e\n", nodes_x[i], nodes_y[i], nodes_z[i]);
//     }
    
//     // Write cells (triangles)
//     fprintf(fp_out, "\nCELLS %d %d\n", num_triangles, num_triangles * 4);
//     for (i = 0; i < num_triangles; i++) {
//         fprintf(fp_out, "3 %d %d %d\n", 
//                 triangle_conn[i * 3 + 0],
//                 triangle_conn[i * 3 + 1],
//                 triangle_conn[i * 3 + 2]);
//     }
    
//     // Write cell types (5 = triangle)
//     fprintf(fp_out, "\nCELL_TYPES %d\n", num_triangles);
//     for (i = 0; i < num_triangles; i++) {
//         fprintf(fp_out, "5\n");
//     }
    
//     // Write point data
//     fprintf(fp_out, "\nPOINT_DATA %d\n", num_nodes);
    
//     // Write velocity vector field (interpolated values)
//     fprintf(fp_out, "VECTORS velocity double\n");
//     for (i = 0; i < num_nodes; i++) {
//         fprintf(fp_out, "%.16e %.16e %.16e\n", u_nodes[i], v_nodes[i], w_nodes[i]);
//     }
    
//     // Write pressure scalar field (interpolated values)
//     fprintf(fp_out, "\nSCALARS pressure double 1\n");
//     fprintf(fp_out, "LOOKUP_TABLE default\n");
//     for (i = 0; i < num_nodes; i++) {
//         fprintf(fp_out, "%.16e\n", p_nodes[i]);
//     }
    
//     fclose(fp_out);
    
//     // Cleanup
//     free(nodes_x);
//     free(nodes_y);
//     free(nodes_z);
//     free(triangle_conn);
//     free(u_nodes);
//     free(v_nodes);
//     free(w_nodes);
//     free(p_nodes);
    
//     printf("Successfully converted %s to %s\n", gmsh_filename, vtk_filename);
//     printf("Nodes: %d, Triangles: %d\n", num_nodes, num_triangles);
    
//     return 0;
// }