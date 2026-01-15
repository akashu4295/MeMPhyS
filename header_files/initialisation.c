#include "structures.h"
#include "functions.h"
#include "mat_lib.h"
#include "kdtree.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* trim whitespace in-place */
static char* trim(char* s)
{
    while (isspace((unsigned char)*s)) s++;
    if (*s == 0) return s;

    char* end = s + strlen(s) - 1;
    while (end > s && isspace((unsigned char)*end)) end--;
    end[1] = '\0';
    return s;
}

BCType parse_bc_type(const char* s)
{
    if (!strcmp(s, "wall")) return BC_WALL;
    if (!strcmp(s, "velocity_inlet")) return BC_VELOCITY_INLET;
    if (!strcmp(s, "velocity_outlet")) return BC_VELOCITY_OUTLET;
    if (!strcmp(s, "pressure_outlet")) return BC_PRESSURE_OUTLET;
    return BC_INTERIOR;
}

void read_physical_names(char* meshfile, PointStructure* ps)
{
    FILE* file = fopen(meshfile, "r");
    char temp[128];
    int n;
    while (fscanf(file, "%s", temp) != EOF) {
        if (strcmp(temp, "$PhysicalNames") == 0)
            break;
    }
    fscanf(file, "%d", &n);
    printf("Number of physical names: %d\n", n);
    ps->num_boundary_types = n;
    for (int i = 0; i < n; i++) {
        int dim;
        char raw[64];
        fscanf(file, "%d %hd %63s", &dim, &ps->boundary_map[i].physical_id, raw);

        int len = strlen(raw);
        if (raw[0] == '"' && raw[len-1] == '"') {
            raw[len-1] = '\0';
            strcpy(ps->boundary_map[i].name, raw + 1);
        } else {
            strcpy(ps->boundary_map[i].name, raw);
        }
        ps->boundary_map[i].bc.type = BC_INTERIOR; // default
    }
    fclose(file);
}

void read_boundary_conditions_file(char* bcfile, PointStructure* ps)
{
    FILE* file = fopen(bcfile, "r");
    char line[512];

    if (!file) {
        perror("Cannot open BC file");
        exit(1);
    }

    int lineno = 0;

    while (fgets(line, sizeof(line), file)) {
        lineno++;

        /* strip newline */
        line[strcspn(line, "\r\n")] = 0;

        char* s = trim(line);

        /* skip empty lines and comments */
        if (*s == '\0' || *s == '#')
            continue;

        BCValue bc = {0};
        bc.type = BC_INTERIOR;

        char* tokens[16];
        int ntok = 0;

        /* tokenize by commas */
        char* tok = strtok(s, ",");
        while (tok && ntok < 16) {
            tokens[ntok++] = trim(tok);
            tok = strtok(NULL, ",");
        }

        if (ntok < 2) {
            printf("BC CSV error (line %d): need at least name,type\n", lineno);
            continue;
        }

        const char* name = tokens[0];
        const char* type = tokens[1];

        bc.type = parse_bc_type(type);
        if (bc.type == BC_INTERIOR) {
            printf("BC CSV error (line %d): unknown BC type '%s'\n",
                   lineno, type);
            continue;
        }

        int have_u = 0, have_v = 0, have_w = 0, have_p = 0;

        /* parse key=value pairs */
        for (int i = 2; i < ntok; i++) {
            char* eq = strchr(tokens[i], '=');
            if (!eq) {
                printf("BC CSV warning (line %d): ignoring token '%s'\n",
                       lineno, tokens[i]);
                continue;
            }

            *eq = 0;
            char* key = trim(tokens[i]);
            char* val = trim(eq + 1);

            if (!strcmp(key, "u")) {
                bc.u = atof(val); have_u = 1;
            }
            else if (!strcmp(key, "v")) {
                bc.v = atof(val); have_v = 1;
            }
            else if (!strcmp(key, "w")) {
                bc.w = atof(val); have_w = 1;
            }
            else if (!strcmp(key, "p")) {
                bc.p = atof(val); have_p = 1;
            }
            else {
                printf("BC CSV warning (line %d): unknown key '%s'\n",
                       lineno, key);
            }
        }

        /* semantic validation */
        if ((bc.type == BC_VELOCITY_INLET ||
             bc.type == BC_VELOCITY_OUTLET ||
             bc.type == BC_WALL) &&
            !(have_u && have_v && have_w))
        {
            printf("BC CSV error (line %d): velocity BC requires u,v,w\n", lineno);
            continue;
        }

        if (bc.type == BC_PRESSURE_OUTLET && !have_p) {
            printf("BC CSV error (line %d): pressure BC requires p\n", lineno);
            continue;
        }

        /* map BC to physical name */
        int found = 0;
        for (int i = 0; i < ps->num_boundary_types; i++) {
            if (!strcmp(ps->boundary_map[i].name, name)) {
                ps->boundary_map[i].bc = bc;
                found = 1;
                break;
            }
        }

        if (!found) {
            printf("BC CSV warning (line %d): physical name '%s' not found in mesh\n",
                   lineno, name);
        }
    }

    fclose(file);
}

int bc_priority(BCType t)
{
    if (t == BC_VELOCITY_INLET) return 3;
    if (t == BC_WALL) return 2;
    if (t == BC_PRESSURE_OUTLET) return 1;
    return 0;
}

void assign_node_bc(PointStructure* ps, int node, BCValue new_bc)
{
    BCValue* old = &ps->node_bc[node];

    if (bc_priority(new_bc.type) > bc_priority(old->type)) {
        *old = new_bc;
    }
}


void apply_boundary_conditions_from_file(PointStructure* myPointStruct,
                                   FieldVariables* field,
                                   int numlevels)
{
    for (int ii = 0; ii < numlevels; ii++)
    {
        for (int i = 0; i < myPointStruct[ii].num_nodes; i++)
        {
            if (!myPointStruct[ii].boundary_tag[i]) continue;

            switch (myPointStruct[ii].node_bc[i].type)
            {
                case BC_VELOCITY_INLET:
                case BC_VELOCITY_OUTLET:
                    field[ii].u[i] = myPointStruct[ii].node_bc[i].u;
                    field[ii].v[i] = myPointStruct[ii].node_bc[i].v;
                    field[ii].w[i] = myPointStruct[ii].node_bc[i].w;
                    break;

                case BC_PRESSURE_OUTLET:
                    field[ii].p[i] = myPointStruct[ii].node_bc[i].p;
                    myPointStruct[ii].flag_outlets = true;
                    break;

                case BC_WALL:
                    field[ii].u[i] = 0.0;
                    field[ii].v[i] = 0.0;
                    field[ii].w[i] = 0.0;
                    break;

                default:
                    break;
            }
        }
    }
}
