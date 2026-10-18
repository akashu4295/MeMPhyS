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

/////////////////////////////////////////////////////////////////////////////
// Forces on monitor_body: Cd, Cl (pressure + viscous)
/////////////////////////////////////////////////////////////////////////////

static int     body_n    = 0;      // number of solver nodes on monitor_body
static int    *body_node = NULL;   // their solver index (RCM order)
static double *body_ds   = NULL;   // their share of the boundary length

/* build the wall node list and ds from the body edges (once, after RCM); corner shares go to the other end of the edge */
void macro_forces_init(PointStructure* myPointStruct)
{
    if (parameters.U_ref <= 0.0 || parameters.L_ref <= 0.0){
        printf("ERROR: monitor_forces needs U_ref > 0 and L_ref > 0 in flow_parameters.csv\n");
        exit(1);
    }

    int N = myPointStruct->num_nodes;
    double *ds = (double*)calloc(N, sizeof(double));   // per solver node, RCM order

    for (int e = 0; e < myPointStruct->num_body_edges; e++){
        int g1 = myPointStruct->body_edges[2*e], g2 = myPointStruct->body_edges[2*e + 1];
        double dx = myPointStruct->x_gmsh[g2] - myPointStruct->x_gmsh[g1];
        double dy = myPointStruct->y_gmsh[g2] - myPointStruct->y_gmsh[g1];
        double L  = sqrt(dx*dx + dy*dy);

        int s1 = myPointStruct->orig_to_solver[g1], s2 = myPointStruct->orig_to_solver[g2];
        int ns = (s1 >= 0) + (s2 >= 0);                // ends that are solver nodes (corners are not)
        if (ns == 0) continue;
        if (s1 >= 0) ds[myPointStruct->rcm_order[s1]] += L / ns;
        if (s2 >= 0) ds[myPointStruct->rcm_order[s2]] += L / ns;
    }

    body_n = 0;
    for (int i = 0; i < N; i++) if (ds[i] > 0.0) body_n++;
    body_node = (int*)malloc(body_n * sizeof(int));
    body_ds   = (double*)malloc(body_n * sizeof(double));

    double sum_ds = 0.0, sum_nx = 0.0, sum_ny = 0.0;
    for (int i = 0, b = 0; i < N; i++){
        if (ds[i] <= 0.0) continue;
        body_node[b] = i;
        body_ds[b]   = ds[i];
        sum_ds += ds[i];
        sum_nx += ds[i] * myPointStruct->x_normal[i];
        sum_ny += ds[i] * myPointStruct->y_normal[i];
        b++;
    }
    free(ds);

    printf("Monitor body '%s': %d nodes, %d edges, length = %.6f, sum(n ds) = (%.2e, %.2e)\n",
           parameters.monitor_body, body_n, myPointStruct->num_body_edges, sum_ds, sum_nx, sum_ny);

    #pragma acc enter data copyin(body_node[0:body_n], body_ds[0:body_n])
}

/* out[0..5] = Cd, Cd_p, Cd_v, Cl, Cl_p, Cl_v on monitor_body (normals point out of the fluid, into the body) */
void macro_forces(PointStructure* myPointStruct, FieldVariables* field, double* out)
{
    int     n    = body_n;
    int    *node = body_node;
    double *dsb  = body_ds;
    int     nc   = myPointStruct->num_cloud_points;
    double  mu   = parameters.mu;
    double  fxp = 0.0, fxv = 0.0, fyp = 0.0, fyv = 0.0;

    # pragma acc parallel loop default(present) reduction(+:fxp, fxv, fyp, fyv)
    for (int b = 0; b < n; b++){
        int i = node[b];

        /* velocity gradients at the wall node from its RBF cloud (current u, v) */
        double ux = 0.0, uy = 0.0, vx = 0.0, vy = 0.0;
        for (int j = 0; j < nc; j++){
            int k = i*nc + j;
            int c = myPointStruct->cloud_index[k];
            ux += myPointStruct->Dx[k] * field->u[c];
            uy += myPointStruct->Dy[k] * field->u[c];
            vx += myPointStruct->Dx[k] * field->v[c];
            vy += myPointStruct->Dy[k] * field->v[c];
        }

        double nx = myPointStruct->x_normal[i], ny = myPointStruct->y_normal[i];
        double p  = field->p[i], ds = dsb[b];

        fxp += p * nx * ds;
        fyp += p * ny * ds;
        fxv -= mu * (2.0*ux*nx + (uy + vx)*ny) * ds;
        fyv -= mu * ((uy + vx)*nx + 2.0*vy*ny) * ds;
    }

    double q = 0.5 * parameters.rho * parameters.U_ref * parameters.U_ref * parameters.L_ref;
    out[0] = (fxp + fxv) / q;  out[1] = fxp / q;  out[2] = fxv / q;
    out[3] = (fyp + fyv) / q;  out[4] = fyp / q;  out[5] = fyv / q;
}

void macro_forces_free(void)
{
    #pragma acc exit data delete(body_node[0:body_n], body_ds[0:body_n])
    free(body_node);
    free(body_ds);
    body_node = NULL; body_ds = NULL; body_n = 0;
}