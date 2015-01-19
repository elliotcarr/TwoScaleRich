// -----------------------------------------------------------------------------
// compute_Keff.c
//
// Written by Dr Elliot Carr (2013-2014)
// Ecole Centrale Paris and Queensland University of Technology
// 
// This code is part of TwoScalRich.
//
// This Matlab MEX .c file computes a column of the effective conductivity using
// a cell average of the flux and the solution of the periodic cell problem.
// -----------------------------------------------------------------------------


#include "mex.h"
#include <math.h>
#include <stdlib.h>

// Matlab gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Loop counter
    int i; 
    
    // Input 1
    double* u = mxGetPr(prhs[0]);
    
    // Input 2
    const mxArray* user_data = prhs[1];   
    
    // Mesh properties
    const mxArray* micro_mesh = mxGetField(user_data,0,"micro_mesh");
    int no_nodes         = mxGetScalar(mxGetField(micro_mesh,0,"no_nodes"));
    int no_elements      = mxGetScalar(mxGetField(micro_mesh,0,"no_elements"));
    double* elements     = mxGetPr(mxGetField(micro_mesh,0,"elements"));
    double* soil         = mxGetPr(mxGetField(micro_mesh,0,"soil"));
    double* shape_funcs  = mxGetPr(mxGetField(micro_mesh,0,"shape_funcs"));
    double* var_num      = mxGetPr(mxGetField(micro_mesh,0,"variable_number"));
    double* element_area = mxGetPr(mxGetField(micro_mesh,0,"element_area"));
    double cell_area     = mxGetScalar(mxGetField(micro_mesh,0,"cell_area"));
    
    // KA: Conductivity in Soil A
    // KB: Conductivity in Soil B
    // column: index j in periodic_cell_problem.m
    double KA  = mxGetScalar(mxGetField(user_data,0,"KA"));
    double KB  = mxGetScalar(mxGetField(user_data,0,"KB"));    
    int column = mxGetScalar(mxGetField(user_data,0,"column"));
 
    // Allocate memory for arrays
    double* h = mxCalloc(no_nodes,sizeof(double));
    double* q_av = mxCalloc(2,sizeof(double));
   
    for (i=0; i<no_nodes; i++)
	{
        h[i] = u[((int)var_num[i]-1)];
    }
    
    //--------------------------------------------------------------------------
    // Fluxes: Loop over elements, compute flux and perform cell-average
    for (i=0; i<no_elements; i++)
    {
        unsigned int icv[3];
        double q[2];
        double K_av;
        
        icv[0] = elements[0*no_elements+i]-1;
        icv[1] = elements[1*no_elements+i]-1;
        icv[2] = elements[2*no_elements+i]-1;
        
        // Build Gradient Approximations
        double grad_hx = h[icv[0]] * shape_funcs[i] + \
                h[icv[1]] * shape_funcs[no_elements+i] + \
                h[icv[2]] * shape_funcs[2*no_elements+i];
        double grad_hy = h[icv[0]] * shape_funcs[3*no_elements+i] + \
                h[icv[1]] * shape_funcs[4*no_elements+i] + \
                h[icv[2]] * shape_funcs[5*no_elements+i];
        
        switch ((int)soil[i])
        {
            case 1:
                K_av = KA;
                break;
            case 2:
                K_av = KB;
                break;
            default:
                break;
        }
        
        // Flux vector (depends on column)
        switch (column)
        {
            case 1: 
                q[0] = K_av * (grad_hx + 1.0);
                q[1] = K_av * grad_hy;
                break;
            case 2:
                q[0] = K_av * grad_hx;
                q[1] = K_av * (grad_hy + 1.0);
                break;
            default:
                break;
        }
        
        // Add contribution to average
        q_av[0] = q_av[0] + q[0] * element_area[i];
        q_av[1] = q_av[1] + q[1] * element_area[i];        
    }    
    //--------------------------------------------------------------------------
    
    // Copy results to output
    plhs[0] = mxCreateDoubleMatrix(2,1,mxREAL);
    double* out = mxGetPr(plhs[0]);
    for (i=0; i<2; i++)
    {
        out[i] = q_av[i] / cell_area;
    }
    
    // Free memory allocations
    mxFree(h);
    mxFree(q_av);
    
    return;
}




