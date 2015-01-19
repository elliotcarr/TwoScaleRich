#include <iostream>
#include <string>
#include <new>
#include <stdlib.h>
#include <math.h>

#include "BLAS_LAPACK_config.h"

#ifdef USE_ACCELERATE_FRAMEWORK
#include <Accelerate/Accelerate.h>
#endif

#ifdef USE_INTEL_MKL
#include <mkl.h>
#endif

using namespace std;

#include "macro_mesh.h"
#include "micro_mesh.h"
#include "Keff.h"
#include "soil.h"
#include "user_data.h"
#include "exprem_options.h"
#include "exprem_stats.h"

double weighted_norm(double* PhiError, double* u, double AbsTol, double RelTol, int N);
int phipade(double* phiH, double* H, double alpha, int m, int MaxKrylov);
void check_step(double dt, double MinStep);
double check_flag_Gfunc(int flag, double dt);