// -----------------------------------------------------------------------------
// Ffunc1.c
//
// Written by Dr Elliot Carr (2013-2014)
// Ecole Centrale Paris and Queensland University of Technology
// 
// This code is part of TwoScalRich.
//
// This Matlab MEX .c file computes all the spatial approximation (computation 
// and assembly of fluxes, etc.) for the periodic cell problem.
// -----------------------------------------------------------------------------

#include "mex.h"
#include <math.h>
#include <stdlib.h>

#define min(a,b) (((a) < (b)) ? (a) : (b))

// Matlab gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    // Loop counters
    int i,j,k;
    
    // Input 1
    double* u = mxGetPr(prhs[0]);
    
    // Input 2
    const mxArray* user_data = prhs[1];    
    
    // Mesh properties
    const mxArray* micro_mesh = mxGetField(user_data,0,"micro_mesh");
    int     no_nodes     = mxGetScalar(mxGetField(micro_mesh,0,"no_nodes"));
    int     no_elements  = mxGetScalar(mxGetField(micro_mesh,0,"no_elements"));
    int     no_variables = mxGetScalar(mxGetField(micro_mesh,0,"no_variables"));
    double* elements     = mxGetPr(mxGetField(micro_mesh,0,"elements"));
    double* soil         = mxGetPr(mxGetField(micro_mesh,0,"soil"));
    double* normals      = mxGetPr(mxGetField(micro_mesh,0,"normals"));
    double* shape_funcs  = mxGetPr(mxGetField(micro_mesh,0,"shape_funcs"));
    double* CV_area      = mxGetPr(mxGetField(micro_mesh,0,"CV_area"));
    double* var_num      = mxGetPr(mxGetField(micro_mesh,0,"variable_number"));
    double* variable     = mxGetPr(mxGetField(micro_mesh,0,"variable"));
    
    // KA: Conductivity in Soil A
    // KB: Conductivity in Soil B
    // column: index j in periodic_cell_problem.m
    double KA  = mxGetScalar(mxGetField(user_data,0,"KA"));
    double KB  = mxGetScalar(mxGetField(user_data,0,"KB"));
    int column = mxGetScalar(mxGetField(user_data,0,"column"));
    
    // Allocate memory for arrays
    double* g = mxCalloc(no_variables,sizeof(double));
    double* h = mxCalloc(no_nodes,sizeof(double));
    
    for (k = 0; k<no_nodes; k++)
	{
        h[k] = u[((int)var_num[k]-1)];
    }
    
    //--------------------------------------------------------------------------
    // Fluxes: Loop over elements and compute and assemble fluxes
    for (i=0; i<no_elements; i++)
    {
        unsigned int icv[3];
        unsigned int varj;
        unsigned int varjnb;
        double n[2];
        double q[2];
        double qn;
        unsigned int indx;
        unsigned int jnb;
        double Kvec[3];
        double Cvec[3];
        double Se;
        double alpha_Kvec;
        double beta_Kvec;
        double gamma_Kvec;
        double midx;
        double midy;
        double K_av;
        double grad_hx, grad_hy;
        
        icv[0] = elements[0*no_elements+i]-1;
        icv[1] = elements[1*no_elements+i]-1;
        icv[2] = elements[2*no_elements+i]-1;
        
        // Build Gradient Approximations
        grad_hx = h[icv[0]] * shape_funcs[i] + \
                h[icv[1]] * shape_funcs[no_elements+i] + \
                h[icv[2]] * shape_funcs[2*no_elements+i];
        grad_hy = h[icv[0]] * shape_funcs[3*no_elements+i] + \
                h[icv[1]] * shape_funcs[4*no_elements+i] + \
                h[icv[2]] * shape_funcs[5*no_elements+i];
        
        // Loop over edges within element
        for (j=0; j<3; j++)
        {            
            indx = j*2;
            switch(j)
            {
                case 0:
                    jnb = 1;
                    break;
                case 1:
                    jnb = 2;
                    break;
                case 2:
                    jnb = 0;
                    break;
                default:
                    break;
            }
            
            varj   = (int)var_num[icv[j]]-1;
            varjnb = (int)var_num[icv[jnb]]-1;
            
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
            
            // Vector normal to edge j (directed towards node icv[jnb])
            n[0] = normals[indx*no_elements+i];
            n[1] = normals[(indx+1)*no_elements+i];
            
            // Dot Product
            qn = n[0] * q[0] + n[1] * q[1];
            
            // Distribute flux at edge j to appropriate nodes
            g[varj]   = g[varj]   - qn;
            g[varjnb] = g[varjnb] + qn;
        }
    }
    //--------------------------------------------------------------------------

    // Copy results to output
    plhs[0] = mxCreateDoubleMatrix(no_variables,1,mxREAL);
    double* out = mxGetPr(plhs[0]);
    for (k=0;k<no_variables;k++)
    {
        out[k] = g[k] / (min(KA,KB));
    }
    
    // Free memory allocations
    mxFree(g);
    mxFree(h);
    
    return;
}




