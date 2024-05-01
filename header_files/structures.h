// Author :  Akash Unnikrishnan
// Affiliation : Indian Institute of Technology Gandhinagar

#define DIMENSION_MAX 3 // Maximum dimensionality

double dtemp; int itemp; char temp[50]; // temporary variables used multiple times in the code

// Structure to represent a point in DIMENSION_MAX-dimensional space
typedef struct Point {
    double coords[DIMENSION_MAX];
    int index; // Index of the point in the original dataset
} Point;

// Structure to represent the parameters of the problem
struct parameters
    {
    int dimension; //dimension of the problem, default is 3
    int n_cloud_points; //number of cloud points in the domain
    int poly_degree; //degree of the polynomial basis functions
    int phs_degree; //degree of the PHS basis functions
    int cloud_size_multiplier; //multiplier for the number of cloud points
    int n_poly_terms; //number of polynomial terms
    int* pow_x; //power of x in the polynomial basis functions
    int* pow_y; //power of y in the polynomial basis functions
    int* pow_z; //power of z in the polynomial basis functions
    }parameters; // child created only once and used globally

// Structure to represent the mesh data
typedef struct mesh {   
    char* mesh_filename; // name of the mesh file
    int* n_count; // number of nodes
    int* e_count; // number of elements
    Point* points; // coordinates of the nodes and index
    double* nx; // x component of the normal vector
    double* ny; // y component of the normal vector
    double* nz; // z component of the normal vector
    bool* boundary_tag; // Boolean tag for boundary nodes
    bool* corner_tag; // Boolean tag for corner nodes
    int** cloud_index; // indices of the cloud points
    int x_min, x_max, y_min, y_max, z_min, z_max; // bounding box of the mesh
    double d_avg; // average distance between nodes
}mesh;