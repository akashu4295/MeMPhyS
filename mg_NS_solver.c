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

#include "src/c_header_files/functions.h"
#include <time.h>

struct parameters parameters;
struct timer timer;

int main()
{
    clock_t clock_start, clock_program_begin = clock();
    struct timespec wall_start, wall_end;
    clock_gettime(CLOCK_MONOTONIC, &wall_start);
    PointStructure* myPointStruct;
    FieldVariables *field;
    double steady_state_error=0.0;
    int it = 0;
    FILE *file2;

///// Read Parameters, grid filenames and mesh data 
    read_flow_parameters("flow_parameters.csv");     // Read the parameters from the file
    read_grid_filenames(&myPointStruct, "grid_filenames.csv", &parameters.num_levels);     // Read the parameters from the file
    read_complete_mesh_data(myPointStruct, parameters.num_levels);
    timer.time_read = (double)(clock() - clock_program_begin)/CLOCKS_PER_SEC;
    printf("Time taken to read the grids and flow parameters: %lf\n", timer.time_read);
    AllocateMemoryFieldVariables(&field, myPointStruct, parameters.num_levels);
    parameters.dt = calculate_dt(&myPointStruct[0]);

    // write_processed_grid_data(myPointStruct, 1); // mesh coordinates, normals, boundary tags, corner tags, etc. are written to files for verification

    clock_start = clock();    // Start the clock
    for (int ii = 0; ii<parameters.num_levels ; ii = ii +1)
        create_derivative_matrices_vectorised(&myPointStruct[ii]);
    if(parameters.test>0) test_derivatives(myPointStruct, parameters.num_levels, parameters.dimension);
    timer.time_derivatives = (double)(clock()-clock_start)/CLOCKS_PER_SEC;
    printf("Time taken to create derivative matrices: %lf\n", timer.time_derivatives);

////////////// Setting up the boudary condition 
    FILE *bcf = fopen("bc.csv", "r");
    if (bcf != NULL){
        fclose(bcf);  // close immediately, we just tested existence
        printf("Boundary conditions applied from bc.csv file\n");
    }
    else{
        printf("bc.csv not found: using default init.c file for boundary conditions\n");
        initial_conditions(myPointStruct, field, 1);
        boundary_conditions(myPointStruct, field, 1);
    }
    check_restart_file(&myPointStruct[0], &field[0]);
    apply_boundary_conditions(myPointStruct, field, 1);
    for (int ii = 0; ii<parameters.num_levels ; ii = ii +1)
        create_laplacian_for_Poisson_equation_vectorised(&myPointStruct[ii]);
    timer.time_initialisation = (double)(clock()-clock_start)/CLOCKS_PER_SEC;
    printf("Time taken for initialisation: %lf\n", timer.time_initialisation);

////////////// Copy data to GPU memory /////////////
    clock_start = clock();    // Start the clock
    copy_all_data_to_gpu(myPointStruct, field);
    timer.time_copy_to_gpu = (double)(clock()-clock_start)/CLOCKS_PER_SEC;
    printf("Time taken to copy data to GPU: %lf\n", timer.time_copy_to_gpu);
    
////////////// Time stepping loop start and writing solution files///////////// 
    clock_start = clock();    // Start the clock
    file2 = fopen("Convergence.csv", parameters.restart ? "a" : "w");  // on restart, append to the old history

    if (parameters.fractional_step)
        if (parameters.dimension == 3){
            for (it = parameters.start_step; it<parameters.num_time_steps; it++ ) 
            {
                steady_state_error = fractional_step_explicit_vectorised(myPointStruct, field);
                printf("Time step: %d, Steady state error: %e\n", it, steady_state_error);
                fflush(stdout);
                fprintf(file2,"%d, %e\n", it, steady_state_error);
                fflush(file2);
                if (steady_state_error < parameters.steady_state_tolerance){
                    printf("Converged at time step: %d\n", it);
                    break;
                }
                if ((it % parameters.write_interval == 0) || (it == parameters.num_time_steps-1)){
                    #pragma acc update host(field[0].u[0:num_nodes], field[0].v[0:num_nodes], field[0].w[0:num_nodes], field[0].p[0:num_nodes])
                    write_vtk(myPointStruct[0].mesh_filename, field, myPointStruct, it);	
                }
            }
        }
        else{
            for (it = parameters.start_step; it<parameters.num_time_steps; it++ ) 
            {
                steady_state_error = fractional_step_explicit_vectorised_2d(myPointStruct, field);
                printf("Time step: %d, Steady state error: %e\n", it, steady_state_error);
                fflush(stdout);
                fprintf(file2,"%d, %e\n", it, steady_state_error);
                fflush(file2);
                if (steady_state_error < parameters.steady_state_tolerance){
                    printf("Converged at time step: %d\n", it);
                    break;
                }
                if ((it % parameters.write_interval == 0) || (it == parameters.num_time_steps-1)){
                    #pragma acc update host(field[0].u[0:num_nodes], field[0].v[0:num_nodes], field[0].p[0:num_nodes])
                    write_vtk(myPointStruct[0].mesh_filename, field, myPointStruct, it);	
                }
            }
        }
    else{
        if (parameters.dimension == 3){
            for (it = parameters.start_step; it<parameters.num_time_steps; it++ ) 
            {
                steady_state_error = time_implicit_solver_vectorised(myPointStruct, field);
                printf("Time step: %d, Steady state error: %e\n", it, steady_state_error);
                fflush(stdout);
                fprintf(file2,"%d, %e\n", it, steady_state_error);
                fflush(file2);
                if (steady_state_error < parameters.steady_state_tolerance){
                    printf("Converged at time step: %d\n", it);
                    break;
                }
                if ((it % parameters.write_interval == 0) || (it == parameters.num_time_steps-1)){
                    #pragma acc update host(field[0].u[0:num_nodes], field[0].v[0:num_nodes], field[0].w[0:num_nodes], field[0].p[0:num_nodes])
                    write_vtk(myPointStruct[0].mesh_filename, field, myPointStruct, it);    
                }
            } 
        }
        else{
            for (it = parameters.start_step; it<parameters.num_time_steps; it++ ) 
                {
                    steady_state_error = time_implicit_solver_vectorised_2d(myPointStruct, field);
                    printf("Time step: %d, Steady state error: %e\n", it, steady_state_error);
                    fflush(stdout);
                    fprintf(file2,"%d, %e\n", it, steady_state_error);
                    fflush(file2);
                    if (steady_state_error < parameters.steady_state_tolerance){
                        printf("Converged at time step: %d\n", it);
                        break;
                    }
                    if ((it % parameters.write_interval == 0) || (it == parameters.num_time_steps-1)){
                        #pragma acc update host(field[0].u[0:num_nodes], field[0].v[0:num_nodes], field[0].p[0:num_nodes])
                        write_vtk(myPointStruct[0].mesh_filename, field, myPointStruct, it);	
                    }
                }
            }
    }
    fclose(file2); 
    timer.time_total = (double)(clock()-clock_start)/CLOCKS_PER_SEC;
    printf("Time taken for the solver: %lf\n", timer.time_total);

    ////////////// Time stepping loop end ///////////// 
    printf("Time_step, dt : %lf\n",parameters.dt);
    printf("Average distance between nodes: %lf\n",myPointStruct[0].d_avg);
    clock_gettime(CLOCK_MONOTONIC, &wall_end);
    timer.time_total = (wall_end.tv_sec - wall_start.tv_sec)
                        + (wall_end.tv_nsec - wall_start.tv_nsec) / 1e9;
    printf("Time for execution (total, wall-clock): %lf\n", timer.time_total);
    write_solver_data(myPointStruct, steady_state_error, it);
    free_PointStructure(myPointStruct, parameters.num_levels);
    free_field(field, parameters.num_levels);
    return 0;
} 