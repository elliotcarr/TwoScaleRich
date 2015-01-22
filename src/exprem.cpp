// -----------------------------------------------------------------------------
// exprem.cpp
//
// Elliot Carr, Queensland University of Technology
// 
// This code is part of TwoScaleRich.
//
// This file is a C++ implementation of the EEM time-stepping solver, see, e.g.
// - E. J. Carr, T. J. Moroney and I. W. Turner (2011), Applied Mathematics and 
// Computation, 30(14), pp. 6587-6596.
// - E. J. Carr, I. W. Turner and P. Perre (2013), Journal of Computational 
// Physics, 233, pp. 66-82.
//
// -----------------------------------------------------------------------------

#include "exprem.h"

int exprem(int (*Gfunc)(double* g, double* u, user_data_struct user_data), \
        int N, double &t, double tend, double* tspan, double* tspan_solns, \
        double* u0, user_data_struct user_data, exprem_options_struct options, \
        exprem_stats_struct &stats)
{
    
    // Loop counters
    unsigned int i, j, k;
    
    double dtnew;    
    int inc_one = 1;
    static double* u = new double [N];    
    static double dt;
    static int NumSteps;
    static int FuncEvals;
    static int tspan_count;
    
    // EXPREM solver options (see "set_exprem_options.cpp")
    double AbsTol        = options.AbsTol;
    double RelTol        = options.RelTol;
    int    MinKrylov     = options.MinKrylov;
    int    MaxKrylov     = options.MaxKrylov;
    double MinStep       = options.MinStep;
    double MaxStep       = options.MaxStep;
    int    MaxNumSteps   = options.MaxNumSteps;
    double SafetyFac     = options.SafetyFac;
    double MinStepDecFac = options.MinStepDecFac;
    double MaxStepIncFac = options.MaxStepIncFac;  
    
    int Order = 2.0;
    
    // u = u0
    if (t == 0)
    {
        cblas_dcopy(N, u0, inc_one, u, inc_one);
        dt = options.InitStep;
        
        if (dt < MinStep)
        {
            cout << "ERROR: Intialised stepsize (InitStep = " << dt << ") is " << \
                    "smaller than the minimum allowable stepsize (MinStep = " <<  \
                    MinStep << ").\n";
            exit(EXIT_FAILURE);
        }
        
        NumSteps = 0;
        FuncEvals = 0;
        tspan_count = 1;
        cblas_dcopy(N, u0, inc_one, tspan_solns, inc_one);
    }
    
    if (dt < MinStep)
    {
        cout << "ERROR: Stepsize " << dt << ") is " << \
                "smaller than the minimum allowable stepsize (MinStep = " <<  \
                MinStep << ").\n";
        exit(EXIT_FAILURE);   
    }
    
    // Initialise arrays
    static double* g          = new double [N];
    static double* V          = new double [N*(MaxKrylov+1)];
    static double* H          = new double [(MaxKrylov+1)*MaxKrylov];
    static double* phiH       = new double [MaxKrylov*MaxKrylov];
    static double* upert      = new double [N];
    static double* gpert      = new double [N];
    static double* PhiError   = new double [N];
    static double* LocalError = new double [N];
    static double* ghalf      = new double [N];
    static double* utemp      = new double [N];
    static double* uhalf      = new double [N];
    static double* uinterp    = new double [N];
    
    int flag;    
    flag = Gfunc(g, u, user_data);
    FuncEvals = FuncEvals + 1;
    dt = check_flag_Gfunc(flag, dt);
    
    // Initialise V, H and phiH to zero arrays
    for (j=0; j<MaxKrylov; j++)
    {
        for (i=0; i<MaxKrylov; i++)
        {
            phiH[j*MaxKrylov+i] = 0.0;
        }
        for (i=0; i<MaxKrylov+1; i++)
        {
            H[j*(MaxKrylov+1)+i] = 0.0;
        }
    }
    for (j=0; j<MaxKrylov+1; j++)
    {
        for (i=0; i<N; i++)
        {
            V[j*N + i] = 0.0;
        }
    }
    
    double beta, normu, epsilon;
    bool phi_converged = false;
    int m;
    
    // Begin Arnoldi iteration
    cblas_dcopy(N, g, inc_one, V, inc_one); // V(:,1) = g
    beta = cblas_dnrm2(N, g, inc_one); // beta = norm(g,2)
    cblas_dscal(N, 1.0/beta, V, inc_one); // V(:,1) = V(:,1) / beta
    normu = cblas_dnrm2(N, u, inc_one); // normu = norm(u,2)
    
    if (normu == 0.0)
    {
        epsilon = 1.490116119384766e-08;
    }
    else
    {
        epsilon = 1.490116119384766e-08 * normu;
    }
    
    for (k=0; k<MaxKrylov; k++)
    {
        m = k+1;
        
        // upert = u + epsilon*V(:,k)
        cblas_dcopy(N, u, inc_one, upert, inc_one);
        cblas_daxpy(N, epsilon, V+k*N, inc_one, upert, inc_one);
        
        // gpert = g(upert)
        flag = Gfunc(gpert, upert, user_data);        
        FuncEvals = FuncEvals + 1;
        dt = check_flag_Gfunc(flag, dt);
        
        // V(:,k+1) = (gpert - g)/epsilon
        cblas_dcopy(N, gpert, inc_one, V+(k+1)*N, inc_one);
        cblas_daxpy(N, -1.0, g, inc_one, V+(k+1)*N, inc_one);
        cblas_dscal(N, 1.0/epsilon, V+(k+1)*N, inc_one);
        
        // Orthogonalise
        for (j=0; j<(k+1); j++)
        {
            // H(j,k) = V(:,j)' * V(:,k+1)
            H[k*(MaxKrylov+1)+j] = cblas_ddot(N, V+j*N, inc_one, V+(k+1)*N, \
                    inc_one);
            
            // V(:,k+1) = V(:,k+1) - H(j,k)*V(:,j)
            cblas_daxpy(N, -1.0 * H[k*(MaxKrylov+1)+j] , V+j*N, inc_one, \
                    V+(k+1)*N, inc_one);
        }
        
        // Reorthoganlise
        for (j=0; j<(k+1); j++)
        {
            // theta = V(:,j)' * V(:,k+1);
            double theta = cblas_ddot(N, V+j*N, inc_one, V+(k+1)*N, inc_one);
            
            // V(:,k+1) = V(:,k+1) - theta*V(:,j)
            cblas_daxpy(N, -theta, V+j*N, inc_one, V+(k+1)*N, inc_one);
            H[k*(MaxKrylov+1)+j] = H[k*(MaxKrylov+1)+j] + theta;
        }
        
        // H(k+1,k) = norm(V(:,k+1),2)
        H[k*(MaxKrylov+1)+k+1] = cblas_dnrm2(N, V+(k+1)*N, inc_one);
        
        if (H[k*(MaxKrylov+1)+k+1] <= 1.0e-12)
        {
            // Invariant subspace detected
            phipade(phiH, H, dt, m, MaxKrylov);
            phi_converged = true;
            break;
        }
        else
        {
            // V(:,k+1) = V(:,k+1) / H(k+1,k)
            cblas_dscal(N, 1.0/H[k*(MaxKrylov+1)+k+1], V+(k+1)*N, inc_one);
        }
        
        phipade(phiH, H, dt, m, MaxKrylov);
        
        // PhiError = V(:,k+1)
        cblas_dcopy(N, V+(k+1)*N, inc_one, PhiError, inc_one);
        
        // PhiError = beta*dt*H(k+1,k)*phiHm(k,1)*V(:,k+1)
        cblas_dscal(N, beta*dt*H[k*(MaxKrylov+1)+k+1]*phiH[k], PhiError, \
                inc_one);
        
        double PhiErrorNorm = weighted_norm(PhiError, u, AbsTol, RelTol, N);

        if ((dt*PhiErrorNorm) <= 1.0 && k > (MinKrylov-2))
        {
            phi_converged = true;
            break;
        }
        
    }
    
    // Krylov approximation failed to converge (halve time step)  
    double dtold = dt;    
    while (!phi_converged)
    {
        cout << "Krylov subspace approximation failed to converge. " << \
                "Time step halved from " << dt << " to ";
        dt = 0.5 * dt;
        cout << dt << "\n";
        
        check_step(dt, MinStep);
        
        // phiH(1:m,1:m) = phi(dt*H(1:m,1:m))
        phipade(phiH, H, dt, m, MaxKrylov);
        
        // PhiError = V(:,k+1)
        cblas_dcopy(N, V+m*N, inc_one, PhiError, inc_one);

        // PhiError = beta*dt*H(k+1,k)*phiHm(k,1)*V(:,k+1)
        cblas_dscal(N, beta*dt*H[(m-1)*(MaxKrylov+1)+m]*phiH[m-1], PhiError, \
                inc_one);
        
        double PhiErrorNorm = weighted_norm(PhiError, u, AbsTol, RelTol, N);
        if (dt*PhiErrorNorm <= 1.0)
        {
            cout << "Krylov subspace exceeded maxixmum dimension size " << \
                    "(MaxKrylov = " << m << ").\n" << "Stepsize reduced from " \
                    << dtold << " to " << dt << ".\n";
            break;  
        }
    }
    
    while (true)
    {       
        //----------------------------------------------------------------------
        // Full step solution (called utemp)
        // utemp = u + dt*beta*V(:,1:k)*phi(dt*Hk)*e1
        //        
        // utemp = u
        cblas_dcopy(N, u, inc_one, utemp, inc_one);
        
        // utemp = utemp + dt*beta*V(:,1:k)*phiH(1:k,1)
        cblas_dgemv(CblasColMajor, CblasNoTrans, N, m, dt*beta, V, N, phiH, \
                inc_one, 1.0, utemp,inc_one);
        //----------------------------------------------------------------------
        
        //----------------------------------------------------------------------
        // Compute approximate two half step solution (call uhalf)
        //
        // uhalf = u     + dt/2*beta*V(:,1:k)*phi(dt/2*Hk)*e1
        // uhalf = uhalf + dt/2*V(:,1:k)*phi(dt/2*Hk)*V(:,1:k)'*ghalf
        //
        //
        // phiH(1:m,1:m) = phi(0.5*dt*H(1:m,1:m))
        phipade(phiH, H, 0.5*dt, m, MaxKrylov);
        
        // uhalf = u
        cblas_dcopy(N, u, inc_one, uhalf, inc_one);
        // uhalf = uhalf + dt/2*beta*V(:,1:k)*phiH(1:k,1)
        cblas_dgemv(CblasColMajor, CblasNoTrans, N, m, 0.5*dt*beta, V, N, \
                phiH, inc_one, 1.0, uhalf, inc_one);
        
        // ghalf = Gfunc(uhalf)
        flag = Gfunc(ghalf, uhalf, user_data);
        FuncEvals = FuncEvals + 1;
        dt = check_flag_Gfunc(flag, dt);
        
        double* temp1 = new double [m];
        double* temp2 = new double [m];
        
        for (i=0; i<m; i++)
        {
            temp1[i] = 0.0;
            temp2[i] = 0.0;
        }
        
        // temp = Vk'*ghalf
        cblas_dgemv(CblasColMajor, CblasTrans, N, m, 1.0, V, N, ghalf, \
                inc_one, 0.0, temp1, inc_one);
        
        // temp = phiH(1:k,1:k)*temp
        cblas_dgemv(CblasColMajor, CblasNoTrans, m, m, 1.0, phiH, MaxKrylov, \
                temp1, inc_one, 0.0, temp2, inc_one);
        
        // uhalf = uhalf + dt/2*Vk*temp
        cblas_dgemv(CblasColMajor, CblasNoTrans, N, m, 0.5*dt, V, N, temp2, \
                inc_one, 1.0, uhalf, inc_one);        
        //----------------------------------------------------------------------
        
        delete[] temp1;
        delete[] temp2;
        
        //----------------------------------------------------------------------
        // Estimate local error
        //
        // LocalError = uhalf
        cblas_dcopy(N, uhalf, inc_one, LocalError, inc_one);
        // LocalError = uhalf - utemp
        cblas_daxpy(N, -1.0, utemp, inc_one, LocalError, inc_one);
        
        double LocalErrorNorm = weighted_norm(LocalError, u, AbsTol, RelTol, N);
        //----------------------------------------------------------------------
       
        //----------------------------------------------------------------------        
        // Local Error Test
        if (LocalErrorNorm > 1.0)
        {
            // Local error test failed. Reduce time step and repeat.
            double eta = max(SafetyFac * pow(1.0/LocalErrorNorm,\
                    (1.0/(Order+1.0))),MinStepDecFac);
            //cout << "Local error test failed. " << \
            //    "Time step reduced from " << dt << " to ";
            dt = eta * dt;
            //cout << dt << "\n";
            check_step(dt, MinStep);
            phipade(phiH, H, dt, m, MaxKrylov);
            continue;
        }
        else
        {
            // Local error test successful. 
            // Accept computed solution and adjust time step for next step.
            double eta = max(min(SafetyFac * pow(1.0/LocalErrorNorm,\
                    (1.0/(Order+1.0))),MaxStepIncFac),1.0);
            dtnew = eta * dt;
            dtnew = min(dtnew,tend - (t+dt));
            dtnew = min(dtnew,MaxStep);
            break;
        }
        //----------------------------------------------------------------------
        
    }
    
    // Accept solution
    t = t + dt;
    
    // Build solution at time of interest if required
    if (t >= tspan[tspan_count])
    {
        double dt_interp = tspan[tspan_count] - (t - dt);
        phipade(phiH, H, dt_interp, m, MaxKrylov);
        
        // Compute solution at tspan[tspan_count]
        // uinterp = u
        cblas_dcopy(N, u, inc_one, uinterp, inc_one);
        // uinterp = uinterp + dt_interp*beta*V(:,1:m)*phiH(1:m,1)
        cblas_dgemv(CblasColMajor, CblasNoTrans, N, m, dt_interp*beta, V, N, \
                phiH, inc_one, 1.0, uinterp, inc_one);
        // Store utemp
        cblas_dcopy(N, uinterp, inc_one, tspan_solns+(tspan_count)*N, inc_one);
        
        tspan_count = tspan_count + 1;
    }
    
    // Accept solution
    cblas_dcopy(N, utemp, inc_one, u, inc_one);
    NumSteps = NumSteps + 1;
    
    // Update solver statistics
    stats.CurrentStep = dt;
    stats.CurrentTime = t;
    stats.CurrentSoln = u;
    stats.CurrentKryDim = m;
    stats.NumSteps = NumSteps;
    stats.FuncEvals = FuncEvals;
    
    // Update stepsize for next time step
    dt = dtnew;
    
    return 0;
}

double weighted_norm(double* y, double* u, double AbsTol, double RelTol, int N)
{
    // Computes a weighted norm of y (see page 1566 of Hochbruck, Lubich, 
    // Selhofer (1998), SIAM J. Sci. Comput, 19(5), pg. 1552-1574)
    
    unsigned int i;
    double wnorm = 0.0;
    for (i=0; i<N; i++)
    {
        wnorm = wnorm + pow(y[i] / ((RelTol * fabs(u[i])) + AbsTol),2.0);
    }
    return wnorm = sqrt((1/(double)N)*wnorm);
}

void check_step(double dt, double MinStep)
{
    if (dt < MinStep)
    {
        cout << "Stepsize smaller than minimum allowable stepsize of " << \
                MinStep << "\n";
        exit(EXIT_FAILURE);
    }
}

double check_flag_Gfunc(int flag, double dt)
{
    if (flag == -1)
    {
        cout << "Solution non-physical. Time step halved from " << dt << " to ";
        dt = 0.5*dt;
        cout << dt << "\n";
    }
    
    return dt;
}
