// Authors: Dr. Akash Unnikrishnan (1,*) and Prof. Surya Pratap Vanka (2)
//
// Affiliations:
// (1) Faculty of Physics, University of Warsaw, Poland
// (2) University of Illinois at Urbana-Champaign, USA
//
// Contribution note:
// Initial development during PhD research at
// Indian Institute of Technology Gandhinagar, India
//
// Date: February 2026
// Version: 2.5
//
// License: MIT License
// Contact: akash.unnikrishnan@iitgn.ac.in

///////////////////////////////////////////////////////////////////////////////
// NOTES: To run on CPU compile with the command: gcc @sources.txt -o a.out  (Linux/MacOS) or gcc @sources.txt -o a.exe  (Windows)
//        To run on GPU compile with the command: gcc -fopenacc @sources.txt -o a.out  (Linux/MacOS) or gcc -fopenacc @sources.txt -o a.exe  (Windows)
//        For Nvidia GPU, make sure to have the CUDA toolkit installed and properly configured. For AMD GPU, make sure to have the ROCm toolkit installed and properly configured.
//        For Nvidia GPU, it is recommended to compile with nvc for better performance, but it is not mandatory. For AMD GPU, it is recommended to compile with hipcc for better performance, but it is not mandatory.
//        Compile command: nvc -acc @sources.txt -o a.out  (Linux/MacOS) or nvc -acc @sources.txt -o a.exe  (Windows) for Nvidia GPU
//        Run command: ./a.out or ./a.exe (depending on the OS)
//        The code solves the incompressible Navier-Stokes equations using a fractional step method or time implicit method.
//        The spatial discretization is done using polyharmonic spline radial basis functions (PHS-RBF) with appended polynomial basis functions.
//        The linear systems are solved using a geometric multigrid method with Jacobi smoothing.
//        This code is developed for educational and research purposes only.
//        Please acknowledge the use of this code in any publications or presentations.
//        Please send your feedbacks and suggestions to akash.unnikrishnan@iitgn.ac.in
///////////////////////////////////////////////////////////////////////////////

#include "src/lib/functions.h"
#if defined(_WIN32) || defined(_WIN64)
    #include <direct.h>
    #define chdir _chdir
#else
    #include <unistd.h>
#endif

struct parameters parameters;
struct timer timer;

static double elapsed_seconds(struct timespec start, struct timespec end){
    return (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
}

int main(int argc, char *argv[])
{
    // ./a.out <case_dir> runs inside case_dir (inputs and results live there); no argument = current folder
    if (argc > 1){
        if (chdir(argv[1]) != 0){
            printf("ERROR: cannot open case folder '%s'\n", argv[1]);
            exit(1);
        }
        printf("Case folder: %s\n", argv[1]);
    }

    struct timespec clock_start, clock_end, clock_program_begin;
    clock_gettime(CLOCK_MONOTONIC, &clock_program_begin);
    PointStructure* myPointStruct;
    FieldVariables* field;
    
////////////// Read Parameters, grid filenames and mesh data 
    read_parameters_gridfilenames_and_meshdata(&myPointStruct, "flow_parameters.csv", "grid_filenames.csv");
    clock_gettime(CLOCK_MONOTONIC, &clock_end);
    timer.time_read = elapsed_seconds(clock_program_begin, clock_end);
    printf("Time taken to read the grids and flow parameters: %lf\n", timer.time_read);
    AllocateMemoryFieldVariables(&field, myPointStruct, parameters.num_levels);
    parameters.dt = calculate_dt(&myPointStruct[0]);

    if (parameters.write_grid_data) write_processed_grid_data(myPointStruct, 1); // mesh coordinates, normals, boundary tags, corner tags, etc. are written to files for verification

    clock_gettime(CLOCK_MONOTONIC, &clock_start);
    for (int ii = 0; ii<parameters.num_levels ; ii = ii +1)
        create_derivative_matrices_vectorised(&myPointStruct[ii]);
    if(parameters.test>0) test_derivatives(myPointStruct, parameters.num_levels, parameters.dimension);
    clock_gettime(CLOCK_MONOTONIC, &clock_end);
    timer.time_derivatives = elapsed_seconds(clock_start, clock_end);
    printf("Time taken to create derivative matrices: %lf\n", timer.time_derivatives);

////////////// Setting up the boudary condition 
    FILE *bcf = fopen("bc.csv", "r");
    if (bcf != NULL){
        fclose(bcf);  // close immediately, we just tested existence
        printf("Boundary conditions applied from bc.csv file\n");
        for (int i = 0; i < myPointStruct[0].num_boundary_types; i++)
            if (myPointStruct[0].boundary_map[i].bc.type == BC_PRESSURE_OUTLET)
                printf("  %s: pressure outlet, U_c = %g\n",
                       myPointStruct[0].boundary_map[i].name, myPointStruct[0].boundary_map[i].bc.U_c);
    }
    else{
        printf("bc.csv not found: using default init.c file for boundary conditions\n");
        initial_conditions(myPointStruct, field, 1);
        boundary_conditions(myPointStruct, field, 1);
    }
    // check_restart_file(&myPointStruct[0], &field[0]);
    check_restart_file(&myPointStruct[0], &field[0]);
    clock_gettime(CLOCK_MONOTONIC, &clock_end);
    timer.time_initialisation = elapsed_seconds(clock_start, clock_end);
    printf("Time taken for initialisation: %lf\n", timer.time_initialisation);

    // Coloring nodes for parallelization in Gauss-Seidel solver (only for Poisson solver type 2)
    if (parameters.poisson_solver_type == 2)
        for (int ilev = 0; ilev < parameters.num_levels; ilev++)
            setup_point_cloud_multicoloring(&myPointStruct[ilev], &field[ilev]);

    int k = 2;
    for (int i = 0; i<myPointStruct[0].num_nodes; i++){
        field[0].res[i] = sin(2 * 3.14 * k * myPointStruct[0].x[i]) * 
                  sin(2 * 3.14 * k * myPointStruct[0].y[i]);
        field[0].p[i] = 0.0;
    }

    for (int i = 0; i<myPointStruct[1].num_nodes; i++){
        field[1].p[i] = sin(2 * 3.14 * k * myPointStruct[1].x[i]) * 
                  sin(2 * 3.14 * k * myPointStruct[1].y[i]);
    }

////////////// Copy initialized data to GPU memory
    clock_gettime(CLOCK_MONOTONIC, &clock_start);
    copy_all_data_to_gpu(myPointStruct, field);
    clock_gettime(CLOCK_MONOTONIC, &clock_end);
    timer.time_copy_to_gpu = elapsed_seconds(clock_start, clock_end);
    printf("Time taken to copy data to GPU: %lf\n", timer.time_copy_to_gpu);

    FILE *restrf = fopen("restr_coarse.csv", "w");
    for (int i = 0; i<myPointStruct[1].num_nodes; i++){
        fprintf(restrf, "%lf, %lf,%lf\n", myPointStruct[1].x[i], myPointStruct[1].y[i], field[1].p[i]);
    }
    fclose(restrf);


    FS_restrict_residuals_vectorised(&myPointStruct[0], &myPointStruct[1], &field[0], &field[1]);
    FS_prolongate_corrections_vectorised(&myPointStruct[0], &myPointStruct[1], &field[0], &field[1]);

    #pragma acc update self(field[1].source[:myPointStruct[1].num_nodes], field[0].p[:myPointStruct[0].num_nodes])

    double restriction_error_squared = 0.0, prolongation_error_squared = 0.0;
    double restriction_max_error = 0.0, prolongation_max_error = 0.0;
    int restriction_count = 0, prolongation_count = 0;
    for (int i = 0; i < myPointStruct[1].num_nodes; i++) {
        if (myPointStruct[1].boundary_tag[i] || myPointStruct[1].corner_tag[i])
            continue;
        double exact = sin(2 * 3.14 * k * myPointStruct[1].x[i]) *
                       sin(2 * 3.14 * k * myPointStruct[1].y[i]);
        double error = field[1].source[i] - exact;
        restriction_error_squared += error * error;
        if (fabs(error) > restriction_max_error)
            restriction_max_error = fabs(error);
        restriction_count++;
    }
    for (int i = 0; i < myPointStruct[0].num_nodes; i++) {
        if (myPointStruct[0].boundary_tag[i] || myPointStruct[0].corner_tag[i])
            continue;
        double exact = sin(2 * 3.14 * k * myPointStruct[0].x[i]) *
                       sin(2 * 3.14 * k * myPointStruct[0].y[i]);
        double error = field[0].p[i] - exact;
        prolongation_error_squared += error * error;
        if (fabs(error) > prolongation_max_error)
            prolongation_max_error = fabs(error);
        prolongation_count++;
    }
    double restriction_rms_error = sqrt(restriction_error_squared / restriction_count);
    double prolongation_rms_error = sqrt(prolongation_error_squared / prolongation_count);
    printf("Restriction RMS/max error: %.6e / %.6e\n",
           restriction_rms_error, restriction_max_error);
    printf("Prolongation RMS/max error: %.6e / %.6e\n",
           prolongation_rms_error, prolongation_max_error);
    int transfer_test_failed =
        !isfinite(restriction_rms_error) || !isfinite(prolongation_rms_error) ||
        restriction_max_error > 0.1 || prolongation_max_error > 0.1;

    FILE *resf = fopen("res_fine.csv", "w");
    for (int i = 0; i<myPointStruct[0].num_nodes; i++){
        fprintf(resf, "%lf, %lf,%lf\n", myPointStruct[0].x[i], myPointStruct[0].y[i], field[0].res[i]);
    }
    fclose(resf);

    FILE *sourf = fopen("source_coarse.csv", "w");
    for (int i = 0; i<myPointStruct[1].num_nodes; i++){
        fprintf(sourf, "%lf, %lf,%lf\n", myPointStruct[1].x[i], myPointStruct[1].y[i], field[1].source[i]);
    }
    fclose(sourf);

    FILE *prolf = fopen("prol_fine.csv", "w");
    for (int i = 0; i<myPointStruct[0].num_nodes; i++){
        fprintf(prolf, "%lf, %lf,%lf\n", myPointStruct[0].x[i], myPointStruct[0].y[i], field[0].p[i]);
    }
    fclose(prolf);    



//////////////////////////////////////
    free_all_data_from_gpu(myPointStruct, field);
    free_all_memory_from_cpu(myPointStruct, field);
    return transfer_test_failed ? EXIT_FAILURE : EXIT_SUCCESS;
} 