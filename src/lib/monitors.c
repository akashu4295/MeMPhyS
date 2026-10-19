#include "functions.h"

// Monitors: monitor file handling

static FILE *monitor_file = NULL;   // stays NULL when monitoring is off

void monitors_init(PointStructure* myPointStruct)
{
    if (parameters.monitor_interval <= 0) return;

    monitor_file = fopen("Monitor.csv", parameters.restart ? "a" : "w");   // on restart, append to the old history
    if (!monitor_file){
        printf("WARNING: cannot open Monitor.csv, monitoring disabled\n");
        return;
    }

    if (parameters.monitor_cd || parameters.monitor_cl)
        macro_forces_init(&myPointStruct[0]);

    if (parameters.num_probes > 0)
        macro_probes_init(&myPointStruct[0]);
    
    fseek(monitor_file, 0, SEEK_END);
    if (ftell(monitor_file) > 0) return;   // appending to an existing file: header already there

    fprintf(monitor_file, "time_step, flow_time");
    if (parameters.monitor_max_velocity) fprintf(monitor_file, ", max_velocity, max_abs_pressure");
    if (parameters.monitor_cd) fprintf(monitor_file, ", Cd, Cd_p, Cd_v");
    if (parameters.monitor_cl) fprintf(monitor_file, ", Cl, Cl_p, Cl_v");
    for (int k = 0; k < parameters.num_probes; k++){
        fprintf(monitor_file, ", probe%d_u, probe%d_v", k + 1, k + 1);
        if (parameters.dimension == 3) fprintf(monitor_file, ", probe%d_w", k + 1);
    }
    fprintf(monitor_file, "\n");
    fflush(monitor_file);
}

void monitors_write(PointStructure* myPointStruct, FieldVariables* field, int it)
{
    if (!monitor_file || it % parameters.monitor_interval != 0) return;

    fprintf(monitor_file, "%d, %.10e", it, it * parameters.dt);

    if (parameters.monitor_max_velocity){
        double out[2];
        macro_max_velocity(&myPointStruct[0], &field[0], out);
        fprintf(monitor_file, ", %.10e, %.10e", out[0], out[1]);
    }

    if (parameters.monitor_cd || parameters.monitor_cl){
        double f[6];
        macro_forces(&myPointStruct[0], &field[0], f);
        if (parameters.monitor_cd) fprintf(monitor_file, ", %.10e, %.10e, %.10e", f[0], f[1], f[2]);
        if (parameters.monitor_cl) fprintf(monitor_file, ", %.10e, %.10e, %.10e", f[3], f[4], f[5]);
    }

    if (parameters.num_probes > 0){
        double pv[3*MAX_PROBES];
        int nv = (parameters.dimension == 3) ? 3 : 2;
        macro_probes(&myPointStruct[0], &field[0], pv);
        for (int j = 0; j < nv * parameters.num_probes; j++)
            fprintf(monitor_file, ", %.10e", pv[j]);
    }
    
    fprintf(monitor_file, "\n");
    fflush(monitor_file);
}

void monitors_close(void)
{
    if (monitor_file && (parameters.monitor_cd || parameters.monitor_cl))
        macro_forces_free();
    if (monitor_file) fclose(monitor_file);
    monitor_file = NULL;
}