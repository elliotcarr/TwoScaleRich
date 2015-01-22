// -----------------------------------------------------------------------------
// phipade.cpp
//
// Elliot Carr, Queensland University of Technology
// 
// This code is part of TwoScaleRich.
//
// This file computes computes the (6,6) Pade approximation for the 'phi' 
// function required in the EEM time-stepping solver exprem.cpp.
//  
// -----------------------------------------------------------------------------

#include "phipade.h"

int phipade(double* phiH, double* H, double alpha, int m, int MaxKrylov)
{    
    
	if (m == 1)
	{
        // When m = 1, H(1:m,1:m) is a scalar
        
        double A = alpha*H[0];
		double expA;
		double num;
		
		expA = exp(A);
		
		if (expA == 0.0)
		{
			num = -1.0;
		}
		else if (expA == 1.0)
		{		
			num = 0.0;
		}
		else 
		{		
			num = (expA - 1.0) * A / log(expA);
		}
		phiH[0] = num / A;

	}
	else
	{
        
        unsigned int i, j;
        double a[6];
        double b[6];
        
        //----------------------------------------------------------------------
        // Rational Pade approximation takes the form:
        // 
        //  phi(A) =   1 + sum_{j=0}^{5} a[i] A^i
        //            ----------------------------
        //             1 + sum_{j=0}^{5} b[i] A^i
        //  
        // where A = alpha*H(1:m,1:m). 
        //
        // See page 1567 of Hochbruck, Lubich, Selhofer (1998), SIAM J. Sci. 
        // Comput, 19(5), pg. 1552-1574.
        // 
        //----------------------------------------------------------------------
        
        // Pade approximant coefficients
        a[0] = 1.0/26.0;
		a[1] = 5.0/156.0;
		a[2] = 1.0/858.0;
		a[3] = 1.0/5720.0;
		a[4] = 1.0/205920.0;
		a[5] = 1.0/8648640.0;
		
		b[0] = -6.0/13.0;
		b[1] =  5.0/52.0;
		b[2] = -5.0/429.0;
		b[3] =  1.0/1144.0;
		b[4] = -1.0/25740.0;
		b[5] =  1.0/1235520.0;       
        
        int m_sqrd = m*m;
        int inc_one = 1;
        
        double* A    = new double [m_sqrd];
        double* Apow = new double [m_sqrd];
        double* Den  = new double [m_sqrd];
        double* phiA = new double [m_sqrd];
        double* expA = new double [m_sqrd];
        double* temp = new double [m_sqrd];
        
        // Copy H(1:m,1:m) to A
        for (i=0; i<m; i++)
        {
            // Copy H(1:m,m) to A(:,m)
            cblas_dcopy(m, H+i*(MaxKrylov+1), inc_one, A+i*m, inc_one);
        }
        
        // A = alpha*A
        cblas_dscal(m_sqrd, alpha, A, inc_one);
        
        // Compute inf norm of A
        char norm_type = 'i';
        
        
        int dim = m;
        int lda = m;
        double work[m];
        double normA = dlange_(&norm_type, &dim, &dim, A, &lda, work);
        int scale = max(ceil(1+log(normA)/log(2.0)),0);
        
        // A = A / 2^scale
        cblas_dscal(m_sqrd, 1.0/pow(2.0,(double)scale), A, inc_one);         
        
        // Apow = A
        cblas_dcopy(m_sqrd, A, inc_one, Apow, inc_one); 
        
        // Initialise Den and phiA to be the identity matrix
        // Initialise expA to be the idenity matrix to use dgemm later on
        for (j=0; j<m; j++)
        {
            for (i=0; i<m; i++)
            {
                if (i == j)
                {
                    Den[j*m+i]  = 1.0;
                    phiA[j*m+i] = 1.0;
                    expA[j*m+i] = 1.0;
                }
                else
                {
                    Den[j*m+i]  = 0.0;
                    phiA[j*m+i] = 0.0;
                    expA[j*m+i] = 0.0;
                }
                temp[j*m+i] = 0.0;
            }
        }
        
        // phiA = phiA + a[0]*Apow
        cblas_daxpy(m_sqrd, a[0], Apow, inc_one,  phiA, inc_one); 
        // Den = Den + b[0]*Apow
        cblas_daxpy(m_sqrd, b[0], Apow, inc_one,   Den, inc_one); 
        
        for (i=1; i<6; i++)
        {
            // Apow = Apow*A
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, \
                    1.0, Apow, m, A, m, 0.0, temp, m); 
            cblas_dcopy(m_sqrd, temp, inc_one, Apow, inc_one);
            // Num = Num + a[i]*Apow
            cblas_daxpy(m_sqrd, a[i], Apow, inc_one, phiA, inc_one); 
            // phiH = phiH + b[i]*Apow
            cblas_daxpy(m_sqrd, b[i], Apow, inc_one, Den, inc_one); 
        }
        
        // Solve linear system: phiA = inv(Den)*phiA
        int nrhs = m;
        int ldNum = m;
        int ipiv[m];
        int ldb = m;
        int info;
        dgesv_(&dim, &nrhs, Den, &dim, &ipiv[0], phiA, &dim, &info); 
        
        if (info != 0)
        {
            cout << "Error: Solution unsuccessful in Pade approximation at "
                    << "Line 160." << "\n";
            exit(EXIT_FAILURE);
        }               
        
        // Undo scaling
        if (scale > 0)
        {
            // expA = A*phiA + I
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, \
                    1.0,    A, m, phiA, m, 1.0, expA, m); 
            
            // phiA = (expA + I) * phiA/2.0  = 0.5*expA*phiA + 0.5*phiA
            cblas_dcopy(m_sqrd, phiA, inc_one, temp, inc_one);
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m ,m, m, \
                    0.5, expA, m, phiA, m, 0.5, temp, m);
            cblas_dcopy(m_sqrd, temp, inc_one, phiA, inc_one);
 
            for (i=1; i<scale; i++)
            {
                // expA = expA*expA
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m ,m, \
                        m, 1.0, expA, m, expA, m, 0.0, temp, m); 
                cblas_dcopy(m_sqrd, temp, inc_one, expA, inc_one);
                
                // phiA = (expA + I) * phiA/2.0  = 0.5*expA*phiA + 0.5*phiA
                cblas_dcopy(m_sqrd, phiA, inc_one, temp, inc_one);
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m ,m, \
                        m, 0.5, expA, m, phiA, m, 0.5, temp, m);
                cblas_dcopy(m_sqrd, temp, inc_one, phiA, inc_one);
            }
        }
        
        // Copy phiA to phiH(1:m,1:m)
        for (i=0; i<m; i++)
        {
            // Copy phiA(:,m) to phiH(1:m,m)
            cblas_dcopy(m, phiA+i*m, inc_one, phiH+i*MaxKrylov, inc_one);
        }
        
        
        // Free memory allocations
        delete[] expA;
        delete[] phiA;
        delete[] Den;
        delete[] Apow;
        delete[] A;
        delete[] temp;
		
	}
    
    return 0;
    
}