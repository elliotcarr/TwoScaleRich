// -----------------------------------------------------------------------------
// set_exprem_options.cpp
//
// Elliot Carr, Queensland University of Technology
// 
// This code is part of TwoScaleRich.
//
// This file sets the options for the time stepping solver "exprem.cpp"
//
// -----------------------------------------------------------------------------

#include "set_exprem_options.h"

int set_exprem_options(exprem_options_struct &options)
{
    options.AbsTol        = 1.0e-5; // Absolute error tolerance
    options.RelTol        = 1.0e-5; // Relative error tolerance
    options.MinKrylov     = 2; // Min Krylov subspace dimension
    options.MaxKrylov     = 100; // Max Krylov subspace dimension
    options.MinStep       = 1.0e-10; // Min stepsize
    options.MaxStep       = 1000.0; // Max stepsize
    options.MaxNumSteps   = 100000; // Max Number of time steps permitted
    options.SafetyFac     = 0.8; // Safety factor for local error test
    options.MinStepDecFac = 0.1; // Min factor stepsize can decrease by
    options.MaxStepIncFac = 1.1; // Max factor stepsize can increase by
    options.InitStep      = 1.0e-3; // Initial stepsize
    return 0;
}