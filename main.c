#include "header_files/get_nodes.c"

// Main function 
int main()
{   
    // Define the parameters
    read_parameters("parameters.csv", &parameters);
    parameters.n_cloud_points = 10; //number of cloud points in the domain
    struct mesh* m1;
    char mesh_filename[50] = "mesh/Square_n_10_unstruc.msh";
    m1 = readmesh(mesh_filename);

    free_mesh(m1);
    return(0);
}