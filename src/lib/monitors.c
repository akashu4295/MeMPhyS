#include "functions.h"

// Monitors: monitor file handling

static FILE *monitor_file = NULL;   // stays NULL when monitoring is off

void monitors_init(void)
{
    if (parameters.monitor_interval <= 0) return;

    monitor_file = fopen("Monitor.csv", parameters.restart ? "a" : "w");   // on restart, append to the old history
    if (!monitor_file){
        printf("WARNING: cannot open Monitor.csv, monitoring disabled\n");
        return;
    }

    fseek(monitor_file, 0, SEEK_END);
    if (ftell(monitor_file) > 0) return;   // appending to an existing file: header already there

    fprintf(monitor_file, "time_step, flow_time");
    if (parameters.monitor_max_velocity) fprintf(monitor_file, ", max_velocity, max_abs_pressure");
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

    fprintf(monitor_file, "\n");
    fflush(monitor_file);
}

void monitors_close(void)
{
    if (monitor_file) fclose(monitor_file);
    monitor_file = NULL;
}