// -----------------------------------------------------------------------------
// Ffunc2.c
//
// Written by Dr Elliot Carr (2013-2014)
// Ecole Centrale Paris and Queensland University of Technology
// 
// This code is part of TwoScalRich.
//
// This Matlab MEX .c file  computes the cell-average over the unit cell
// (see solve_periodic_problem.m).
// -----------------------------------------------------------------------------

#include "mex.h"
#include <math.h>
#include <stdlib.h>

// Matlab gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Loop counter
    int k;
    
    // Input 1
    double* u = mxGetPr(prhs[0]);
    
    // Input 2
    const mxArray* user_data = prhs[1];   
    
    // Mesh properties
    const mxArray* micro_mesh = mxGetField(user_data,0,"micro_mesh");
    int     no_nodes  = mxGetScalar(mxGetField(micro_mesh,0,"no_nodes"));
    double* CV_area   = mxGetPr(mxGetField(micro_mesh,0,"CV_area"));
    double* var_num   = mxGetPr(mxGetField(micro_mesh,0,"variable_number"));
    double  cell_area = mxGetScalar(mxGetField(micro_mesh,0,"cell_area"));
    
    // Loop over nodes and cell-average       
    double f = 0.0;
    for (k=0; k<no_nodes; k++)
	{
        double h = u[((int)var_num[k]-1)];
        f = f + h * CV_area[k];
    }   
    
    // Copy results to output
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    double* out = mxGetPr(plhs[0]);
    out[0] = f / cell_area;
    
    return;
}




