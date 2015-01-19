#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <new>
#include <iomanip>
#include <ctime>
#include <omp.h>

using namespace std;

#include "macro_mesh.h"
#include "micro_mesh.h"
#include "Keff.h"
#include "soil.h"
#include "user_data.h"
#include "exprem_options.h"
#include "exprem_stats.h"

int Gfunc_serial(double* g, double* u, user_data_struct user_data);
int Gfunc_parallel(double* g, double* u, user_data_struct user_data);
int macro_mesh_properties(const char* mesh_file, \
        macro_mesh_struct &macro_mesh, int* macro_bc);
int micro_mesh_properties(const char* mesh_file, micro_mesh_struct &micro_mesh);
int read_effective_conductivity(const char* Keff_file, \
        Keff_struct &effective_conductivity);
int set_exprem_options(exprem_options_struct &options);
int exprem(int (*Gfunc)(double* g, double* u, user_data_struct user_data), \
        int N, double &t, double tend, double* tspan, double* tspan_solns, \
        double* u0, user_data_struct user_data, exprem_options_struct options, \
        exprem_stats_struct &stats);
void save_solution(string solution_file, double* tspan, double* tspan_solns, 
        int no_solns, user_data_struct user_data);
