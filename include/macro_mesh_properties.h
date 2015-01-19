#include <iostream>
#include <string>
#include <fstream>
#include <new>
#include <stdlib.h>
#include <vector>
#include <math.h>

using namespace std;

#include "macro_mesh.h"

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

double two_norm(double x1, double y1, double x2, double y2);
