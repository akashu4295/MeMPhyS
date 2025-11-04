// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign

#ifndef SOLVERS_H
#define SOLVERS_H

#include <time.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "structures.h"

////////////////////////////////////////////////////////////////////////////////////
// Function Declarations
//////////////////////////////////////////////////////////////////////////////////

// Time Implicit Solver Modules

double time_implicit_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void calculate_intermediate_velocity_implicit_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void calculate_mass_residual_implicit_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void multigrid_Poisson_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void update_velocity_implicit_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void relaxation_vectorised(PointStructure* mypointstruct, FieldVariables* field);
void calculate_residuals_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void restrict_residuals_vectorised(PointStructure* myPointStruct_f,PointStructure* myPointStruct, FieldVariables* field_f, FieldVariables* field);
void prolongate_corrections_vectorised(PointStructure* myPointStruct_f,PointStructure* myPointStruct, FieldVariables* field_f, FieldVariables* field);
void update_boundary_pressure_vectorised(PointStructure* mypointstruct, FieldVariables* field);
void update_boundary_pprime_vectorised(PointStructure* mypointstruct, FieldVariables* field);




double time_implicit_solver_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void calculate_intermediate_velocity_implicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void calculate_mass_residual_implicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void multigrid_Poisson_solver_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void update_velocity_implicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void relaxation_vectorised_2d(PointStructure* mypointstruct, FieldVariables* field);
void calculate_residuals_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void restrict_residuals_vectorised_2d(PointStructure* myPointStruct_f,PointStructure* myPointStruct, FieldVariables* field_f, FieldVariables* field);
void prolongate_corrections_vectorised_2d(PointStructure* myPointStruct_f,PointStructure* myPointStruct, FieldVariables* field_f, FieldVariables* field);
void update_boundary_pressure_vectorised_2d(PointStructure* mypointstruct, FieldVariables* field);
void update_boundary_pprime_vectorised_2d(PointStructure* mypointstruct, FieldVariables* field);



// Fractional Step Explicit Solver Modules

double fractional_step_explicit_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_intermediate_velocity_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_mass_residual_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_boundary_dpdn_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_multigrid_Poisson_solver_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_relaxation_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_residuals_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_restrict_residuals_vectorised(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, FieldVariables* field_f, FieldVariables* field_c);
void FS_prolongate_corrections_vectorised(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, FieldVariables* field_f, FieldVariables* field_c);
void FS_update_velocity_vectorised(PointStructure* myPointStruct, FieldVariables* field);
void FS_update_boundary_pressure_vectorised(PointStructure* myPointStruct, FieldVariables* field);

double fractional_step_explicit_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_intermediate_velocity_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_mass_residual_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_boundary_dpdn_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_multigrid_Poisson_solver_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_relaxation_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_calculate_residuals_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_restrict_residuals_vectorised_2d(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, FieldVariables* field_f, FieldVariables* field_c);
void FS_prolongate_corrections_vectorised_2d(PointStructure* myPointStruct_f, PointStructure* myPointStruct_c, FieldVariables* field_f, FieldVariables* field_c);
void FS_update_velocity_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);
void FS_update_boundary_pressure_vectorised_2d(PointStructure* myPointStruct, FieldVariables* field);


#endif