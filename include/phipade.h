#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "BLAS_LAPACK_config.h"

#ifdef USE_ACCELERATE_FRAMEWORK
#include <Accelerate/Accelerate.h>
#endif

#ifdef USE_INTEL_MKL
#include <mkl.h>
#endif

using namespace std;

#define max(a,b) (((a) > (b)) ? (a) : (b))