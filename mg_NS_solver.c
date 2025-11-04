// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar
//                University of Illinois at Urbana Champaign
//                 Faculty of Physics, University of Warsaw
// Date : June 2025
// Version : 2.0

///////////////////////////////////////////////////////////////////////////////
// NOTES: This code is work in progress for parallelization using OpenACC directives.
//        However it runs on cpu as well.
//        compile command: gcc @sources.txt -o a.out  (Linux/MacOS)
//        compile command: gcc @sources.txt -o a.exe  (Windows)
//        Run command: ./a.out or ./a.exe (depending on the OS)
//        The code solves the incompressible Navier-Stokes equations using a fractional step method or time implicit method.
//        The spatial discretization is done using polyharmonic spline radial basis functions (PHS-RBF) with appended polynomial basis functions.
//        The linear systems are solved using a geometric multigrid method with Jacobi smoothing.
//        This code is developed for educational and research purposes only.
//        Please acknowledge the use of this code in any publications or presentations.
//        Please send your feedbacks and suggestions to akash.unnikrishnan@iitgn.ac.in
///////////////////////////////////////////////////////////////////////////////

////////////// Header files

#include "header_files/functions.h"
#include "header_files/functions_supplementary.h"
#include "header_files/solvers.h"
#include "header_files/openACC_functions.h"

////////////// Main Program

struct parameters parameters;

int main()
{
    clock_t clock_start = clock(), clock_program_begin = clock();
    PointStructure* myPointStruct;
    FieldVariables *field;
    double steady_state_error;
    FILE *file1;
    FILE *file2;

///// Read Input File Name and Parameters
    read_flow_parameters("flow_parameters.csv");     // Read the parameters from the file
    read_grid_filenames(&myPointStruct, "grid_filenames.csv", &parameters.num_levels);     // Read the parameters from the file

///// Read Mesh data on all levels and compute derivative matrices
    read_complete_mesh_data(myPointStruct, parameters.num_levels);
    printf("Time taken to read the grids and flow parameters: %lf\n", (double)(clock()-clock_start)/CLOCKS_PER_SEC);

    clock_start = clock();    // Start the clock
    for (int ii = 0; ii<parameters.num_levels ; ii = ii +1)
        create_derivative_matrices_vectorised(&myPointStruct[ii]);
    printf("Time taken to create derivative matrices: %lf\n", (double)(clock()-clock_start)/CLOCKS_PER_SEC);

    clock_start = clock();    // Start the clock
    if(parameters.test>0){
        test_derivatives (myPointStruct, parameters.num_levels, parameters.dimension);
        printf("Time taken to test the derivatives: %lf\n", (double)(clock()-clock_start)/CLOCKS_PER_SEC);
    }

////////////// Setting up the boudary condition (Lid driven cavity) /////////////
    
    AllocateMemoryFieldVariables(&field, myPointStruct, parameters.num_levels);
    initial_conditions(myPointStruct, field, 1);
    boundary_conditions(myPointStruct, field, 1);
    parameters.dt = calculate_dt(&myPointStruct[0]);
    printf("Time to setup the problem in cpu: %lf\n", (double)(clock()-clock_program_begin)/CLOCKS_PER_SEC);

////////////// Copy data to GPU memory /////////////

    clock_start = clock();    // Start the clock
    copyin_essentials_to_gpu(myPointStruct);
    copyin_field_to_gpu(field, myPointStruct);
    copyin_parameters_to_gpu();
    printf("Time taken to copy data to GPU: %lf\n", (double)(clock()-clock_start)/CLOCKS_PER_SEC);
    
////////////// Time stepping loop start and writing solution files///////////// 
    clock_start = clock();    // Start the clock
    file2 = fopen("Convergence.csv", "w"); // Write data to a file
    if (parameters.fractional_step)
        if (parameters.dimension == 3){
            for (int it = 0; it<parameters.num_time_steps; it++ ) 
            {
                steady_state_error = fractional_step_explicit_vectorised(myPointStruct, field);
                printf("Time step: %d, Steady state error: %e\n", it, steady_state_error);
                # pragma acc update host(field[0])
                fprintf(file2,"%d, %e\n", it, steady_state_error);
                fflush(file2);
                if (steady_state_error < parameters.steady_state_tolerance){
                    printf("Converged at time step: %d\n", it);
                    break;
                }
                if ((it % parameters.write_interval == 0) || (it == parameters.num_time_steps-1)){
                    file1 = fopen("Solution.csv", "w"); // Write data to a file
                    for (int i = 0; i < myPointStruct[0].num_nodes; i++)
                        fprintf(file1, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", myPointStruct[0].x[i], myPointStruct[0].y[i], myPointStruct[0].z[i], field[0].u[i], field[0].v[i], field[0].w[i], field[0].p[i]);
                    fflush(file1);	
                    fclose(file1);
                }
            }
        }
        else{
            for (int it = 0; it<parameters.num_time_steps; it++ ) 
            {
                steady_state_error = fractional_step_explicit_vectorised_2d(myPointStruct, field);
                printf("Time step: %d, Steady state error: %e\n", it, steady_state_error);
                # pragma acc update host(field[0])
                fprintf(file2,"%d, %e\n", it, steady_state_error);
                fflush(file2);
                if (steady_state_error < parameters.steady_state_tolerance){
                    printf("Converged at time step: %d\n", it);
                    break;
                }
                if ((it % parameters.write_interval == 0) || (it == parameters.num_time_steps-1)){
                    file1 = fopen("Solution.csv", "w"); // Write data to a file
                    for (int i = 0; i < myPointStruct[0].num_nodes; i++)
                        fprintf(file1, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", myPointStruct[0].x[i], myPointStruct[0].y[i], myPointStruct[0].z[i], field[0].u[i], field[0].v[i], field[0].w[i], field[0].p[i]);
                    fflush(file1);	
                    fclose(file1);
                }
            }
        }
    else{
        if (parameters.dimension == 3){
            for (int it = 0; it<parameters.num_time_steps; it++ ) 
            {
                steady_state_error = time_implicit_solver_vectorised(myPointStruct, field);
                printf("Time step: %d, Steady state error: %e\n", it, steady_state_error);
                # pragma acc update host(field[0])
                fprintf(file2,"%d, %e\n", it, steady_state_error);
                fflush(file2);
                if (steady_state_error < parameters.steady_state_tolerance){
                    printf("Converged at time step: %d\n", it);
                    break;
                }
                if ((it % parameters.write_interval == 0) || (it == parameters.num_time_steps-1)){
                    file1 = fopen("Solution.csv", "w"); // Write data to a file
                    for (int i = 0; i < myPointStruct[0].num_nodes; i++)
                        fprintf(file1, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", myPointStruct[0].x[i], myPointStruct[0].y[i], myPointStruct[0].z[i], field[0].u[i], field[0].v[i], field[0].w[i], field[0].p[i]);
                    fflush(file1);	
                    fclose(file1);
                }
            } 
        }
        else{
            for (int it = 0; it<parameters.num_time_steps; it++ ) 
                {
                    steady_state_error = time_implicit_solver_vectorised_2d(myPointStruct, field);
                    printf("Time step: %d, Steady state error: %e\n", it, steady_state_error);
                    # pragma acc update host(field[0])
                    fprintf(file2,"%d, %e\n", it, steady_state_error);
                    fflush(file2);
                    if (steady_state_error < parameters.steady_state_tolerance){
                        printf("Converged at time step: %d\n", it);
                        break;
                    }
                    if ((it % parameters.write_interval == 0) || (it == parameters.num_time_steps-1)){
                        file1 = fopen("Solution.csv", "w"); // Write data to a file
                        for (int i = 0; i < myPointStruct[0].num_nodes; i++)
                            fprintf(file1, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", myPointStruct[0].x[i], myPointStruct[0].y[i], myPointStruct[0].z[i], field[0].u[i], field[0].v[i], field[0].w[i], field[0].p[i]);
                        fflush(file1);	
                        fclose(file1);
                    }
                }
            }
    }
    fclose(file2); 
    printf("Time taken for the solver: %lf\n", (double)(clock()-clock_start)/CLOCKS_PER_SEC);

////////////// Time stepping loop end ///////////// 
    // Write final solution in VTK format
    clock_start = clock();    // Start the clock
    write_vtk(myPointStruct[0].mesh_filename,field,myPointStruct);
    printf("Time taken for VTK writing: %lf\n", (double)(clock()-clock_start)/CLOCKS_PER_SEC);
    printf("Time for execution (total): %lf\n", (double)(clock()-clock_program_begin)/CLOCKS_PER_SEC);

    printf("Time_step, dt : %lf\n",parameters.dt);
    printf("Average distance between nodes: %lf\n",myPointStruct[0].d_avg);
    free_PointStructure(myPointStruct, parameters.num_levels);
    free_field(field, parameters.num_levels);
    return 0;
} 