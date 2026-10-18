#include "functions.h"


// functions to monitor

/* max |u| and max |p| over the solver nodes: out[0] = max |u|, out[1] = max |p| */
void macro_max_velocity(PointStructure* myPointStruct, FieldVariables* field, double* out)
{
    int N = myPointStruct->num_nodes;
    double umax2 = 0.0, pmax = 0.0;

    if (parameters.dimension == 3){
        # pragma acc parallel loop default(present) reduction(max:umax2, pmax)
        for (int i = 0; i < N; i++){
            double u2 = field->u[i]*field->u[i] + field->v[i]*field->v[i] + field->w[i]*field->w[i];
            umax2 = fmax(umax2, u2);
            pmax  = fmax(pmax, fabs(field->p[i]));
        }
    }
    else{
        # pragma acc parallel loop default(present) reduction(max:umax2, pmax)
        for (int i = 0; i < N; i++){
            double u2 = field->u[i]*field->u[i] + field->v[i]*field->v[i];
            umax2 = fmax(umax2, u2);
            pmax  = fmax(pmax, fabs(field->p[i]));
        }
    }

    out[0] = sqrt(umax2);
    out[1] = pmax;
}