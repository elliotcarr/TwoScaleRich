//------------------------------------------------------------------------------
// read_effective_conductivity.cpp
//
// Elliot Carr, Queensland University of Technology
// 
// This code is part of TwoScaleRich.
//
// This file reads the effective conductivity values computed in 
// Effective Conductivity/effective_conductivity.m
// 
//------------------------------------------------------------------------------

#include "read_effective_conductivity.h"

int read_effective_conductivity(const char* Keff_file, \
        Keff_struct &effective_conductivity)
{
    FILE * myfile;
    
    unsigned int i;
    double epsilonA;
    int no_hvalues;
    double logr;
    double logh_first;
    double* logh_values;
    double* Keff_values;
    long double a;
    
    myfile = fopen(Keff_file,"r");
    
    if (myfile == NULL)
    {
        return 1;
    }
    
    fscanf(myfile, "%Lf", &a);
    epsilonA = a;
    fscanf(myfile, "%i" , &no_hvalues);
    fscanf(myfile, "%Lf", &a);
    logr = a;
    fscanf(myfile, "%Lf", &a);
    logh_first = a;
    
    logh_values = new double [no_hvalues];
    
    for (i=0; i<no_hvalues; i++)
    {
        fscanf(myfile, "%Lf", &a);
        logh_values[i] = a;
    }
    
    Keff_values = new double [4*no_hvalues];
    
    for (i=0; i<no_hvalues; i++)
    {
        fscanf(myfile, "%Lf", &a);
        Keff_values[i] = a;
        fscanf(myfile, "%Lf", &a);
        Keff_values[no_hvalues+i] = a;
        fscanf(myfile, "%Lf", &a);
        Keff_values[2*no_hvalues+i] = a;
        fscanf(myfile, "%Lf", &a);
        Keff_values[3*no_hvalues+i] = a;
    }
    
    effective_conductivity.epsilonA    = epsilonA;
    effective_conductivity.no_hvalues  = no_hvalues;
    effective_conductivity.logr        = logr;
    effective_conductivity.logh_first  = logh_first;
    effective_conductivity.logh_values = logh_values;
    effective_conductivity.Keff_values = Keff_values;
    
    return 0;
}