#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct parameters
    {
    int dimension; //dimension of the problem, default is 3
    int n_cloud_points; //number of cloud points in the domain
    int poly_degree; //degree of the polynomial basis functions
    int phs_degree; //degree of the PHS basis functions
    int cloud_size_multiplier; //multiplier for the number of cloud points
    int n_poly_terms; //number of polynomial terms
    }parameters;

void read_parameters(char *filename, struct parameters *parameters) {
    FILE *file;
    char ctemp[100];
    int itemp;
    double dtemp;
    file = fopen(filename, "r");
    if (file == NULL)
    {
        printf("Error: Unable to open the file\n");
        exit(1);
    }
    fscanf(file, "%[^,],%d\n", ctemp, &parameters->poly_degree);
    printf("PARAMETERS: %s = %d\n", ctemp, parameters->poly_degree);
    fscanf(file, "%[^,],%d\n", ctemp, &parameters->phs_degree);
    printf("PARAMETERS: %s = %d\n", ctemp, parameters->phs_degree);
    fscanf(file, "%[^,],%d\n", ctemp, &parameters->cloud_size_multiplier);
    printf("PARAMETERS: %s = %d\n", ctemp, parameters->cloud_size_multiplier);
    parameters->n_poly_terms = (parameters->poly_degree + 1) * (parameters->poly_degree + 2) / 2;
    parameters->n_cloud_points = parameters->cloud_size_multiplier * parameters->n_poly_terms;
}
