//------------------------------------------------------------------------------
// Gfunc_parallel.cpp
//
// Elliot Carr, Queensland University of Technology
// 
// This code is part of TwoScaleRich.
//
// This file is a parallel (OpenMP) version of Gfunc_serial.cpp
//
//------------------------------------------------------------------------------

#include "Gfunc_parallel.h"

int Gfunc_parallel(double* g, double* u, user_data_struct user_data)
{
    // Loop counters
    int i, j, k, p;
    
    // Macro mesh struct
    macro_mesh_struct macro_mesh;
    macro_mesh = user_data.macro_mesh;
    
    int M_no_verts                  = macro_mesh.no_verts;
    int M_no_nodes                  = macro_mesh.no_nodes;
    int M_no_elements               = macro_mesh.no_elements;
    int *M_elements                 = macro_mesh.elements;
    double *M_shape_funcs           = macro_mesh.shape_funcs;
    double *M_normals               = macro_mesh.normals;
    double *M_boundary_edge_midpts  = macro_mesh.boundary_edge_midpts;
    double *M_boundary_edge_normals = macro_mesh.boundary_edge_normals;
    double *M_boundary_edge_lens    = macro_mesh.boundary_edge_lens;
    int *M_boundary_condition_edges = macro_mesh.boundary_condition_edges;
    double *M_CV_area               = macro_mesh.CV_area;
    double *M_SCV_area              = macro_mesh.SCV_area;
    double *M_midx_mat              = macro_mesh.midx;
    double *M_midy_mat              = macro_mesh.midy;
    double *M_centroid_elements     = macro_mesh.centroid_elements;
    
    // Micro mesh struct
    micro_mesh_struct micro_mesh;
    micro_mesh = user_data.micro_mesh;
   
    int m_no_verts                  = micro_mesh.no_verts;
    int m_no_nodes                  = micro_mesh.no_nodes;
    int m_no_elements               = micro_mesh.no_elements;
    int *m_elements                 = micro_mesh.elements;
    int *m_boundary_nodes           = micro_mesh.boundary_nodes;
    double *m_shape_funcs           = micro_mesh.shape_funcs;
    double *m_normals               = micro_mesh.normals;
    int    *m_boundary_edges        = micro_mesh.boundary_edges;
    double *m_boundary_edge_midpts  = micro_mesh.boundary_edge_midpts;
    double *m_boundary_edge_normals = micro_mesh.boundary_edge_normals;
    double *m_boundary_edge_lens    = micro_mesh.boundary_edge_lens;
    double *m_CV_area               = micro_mesh.CV_area;
    double *m_midx_mat              = micro_mesh.midx;
    double *m_midy_mat              = micro_mesh.midy;
    double m_cell_area              = micro_mesh.cell_area;
    double *m_vectors               = micro_mesh.vectors;
    double *m_centroid_elements     = micro_mesh.centroid_elements;
    double *m_element_area          = micro_mesh.element_area;

    // Effective conductivity struct
    double epsilonA     = user_data.effective_conductivity.epsilonA;
    int no_hvalues      = user_data.effective_conductivity.no_hvalues;
    double logr         = user_data.effective_conductivity.logr;
    double logh_first   = user_data.effective_conductivity.logh_first;
    double* logh_values = user_data.effective_conductivity.logh_values;
    double* Keff_values = user_data.effective_conductivity.Keff_values;
    
    // Model formulation
    int model_type = user_data.model_type;
    
    // Connected sub-domain A hydraulic properties
    double KsatA   = user_data.soilA.Ksat;
    double thetarA = user_data.soilA.thetar;
    double thetasA = user_data.soilA.thetas;
    double alphaA  = user_data.soilA.alpha;
    double nA      = user_data.soilA.n;
    double mA      = user_data.soilA.m;
    
    // Disconnected sub-domain B hydraulic properties (inclusions)
    double KsatB   = user_data.soilB.Ksat;
    double thetarB = user_data.soilB.thetar;
    double thetasB = user_data.soilB.thetas;
    double alphaB  = user_data.soilB.alpha;
    double nB      = user_data.soilB.n;
    double mB      = user_data.soilB.m;
    
    double grad_hx;
    double grad_hy;    
    
    double boundary_flux = user_data.boundary_flux;
    
    double epsilonB = 1.0 - epsilonA;
    
    // Total number of unknowns
    int N = user_data.N;

    // Set number of OpenMP threads        
    unsigned int NUM_THREADS = user_data.NUM_THREADS;
    omp_set_num_threads(NUM_THREADS);
    
    // Allocate memory for arrays
    double *h = new double [N];   
    double* f = new double [N];
    // Promote g to an array indexed by thread number to avoid race conditions
    double* garray = new double [N*NUM_THREADS];    
    
    #pragma omp parallel for private(i) 
    for (j=0; j<NUM_THREADS; j++) // Only outer loop is made private by default
    {
        for (i=0; i<N; i++)
        {
            garray[j*N+i] = 0.0;
            g[i] = 0.0;
            f[i] = 0.0;
        }
    }
    
    // Loop over macro nodes and compute specific water capacity term
    #pragma omp parallel for 
    for (k=0; k<M_no_nodes; k++)
    {
        double h_node;
        double Se;
        double C_node;
        
        h_node = u[k];
        
        // Macro node
        h[k] = h_node;
        
        Se     = pow(1.0 + pow(-alphaA*h_node,nA),-1.0 * mA);
        C_node = (thetarA - thetasA) * Se * mA * pow(-alphaA*h_node,nA) * nA / \
                (h_node * (1.0 + pow(-alphaA*h_node,nA)));
        f[k]   = epsilonA * C_node * M_CV_area[k];        
    }
    
    //--------------------------------------------------------------------------
    // Macroscopic fluxes: Loop over macro elements and compute fluxes
    #pragma omp parallel for private (i,j,p)
    for (k=0; k<M_no_elements; k++)
    {
        unsigned int M_icv[M_no_verts];
        unsigned int varj;
        unsigned int varjnb;
        double n[2];
        double q[2];
        double qn;
        unsigned int indx;
        unsigned int jnb;
        double Kvec11[M_no_verts];
        double Kvec12[M_no_verts];
        double Kvec21[M_no_verts];
        double Kvec22[M_no_verts];
        double alpha_h[M_no_verts];
        double alpha_Kvec11[M_no_verts];
        double alpha_Kvec12[M_no_verts];
        double alpha_Kvec21[M_no_verts];
        double alpha_Kvec22[M_no_verts];
        double midx;
        double midy;
        double K_av11;
        double K_av12;
        double K_av21;
        double K_av22;
        double h_av;
        double centx;
        double centy;
        double x_coord;
        double y_coord;
        double Kb11;
        double Kb12;
        double Kb21;
        double Kb22;
        double bflux;
        double Se;
        double h_node;
        double C_node;
        double grad_hx;
        double grad_hy;
        double hb;
        double logh;
        double n_value;
        double logh_diff;
        double Keff_left;
        double Keff_right;
        
        // OpenMP thread number
        int thread_num = omp_get_thread_num();
        
        // Effective conductivity
        // Keff = [Kvec11, Kvec12; Kvec21, Kvec22];
        for (j=0; j<M_no_verts; j++)
        {
            M_icv[j]  = M_elements[j*M_no_elements+k];
            h_node    = h[M_icv[j]];
            
            logh = log10(-h_node);
            n_value = 1.0 + (logh - logh_first) / logr;
            indx = floor(n_value);
            logh_diff = logh_values[indx] - logh_values[indx+1];
            
            Keff_left  = Keff_values[indx+1];
            Keff_right = Keff_values[indx];
            Kvec11[j]  = Keff_left + (Keff_right - Keff_left) * \
                    (logh - logh_values[indx]) / logh_diff;
            
            Keff_left  = Keff_values[no_hvalues + indx+1];
            Keff_right = Keff_values[no_hvalues + indx];
            Kvec12[j]  = Keff_left + (Keff_right - Keff_left) * \
                    (logh - logh_values[indx]) / logh_diff;
            
            Keff_left  = Keff_values[2*no_hvalues + indx+1];
            Keff_right = Keff_values[2*no_hvalues + indx];
            Kvec21[j]  = Keff_left + (Keff_right - Keff_left) * \
                    (logh - logh_values[indx]) / logh_diff;
            
            Keff_left  = Keff_values[3*no_hvalues + indx+1];
            Keff_right = Keff_values[3*no_hvalues + indx];
            Kvec22[j]  = Keff_left + (Keff_right - Keff_left) * \
                    (logh - logh_values[indx]) / logh_diff;
        }
        
        //----------------------------------------------------------------------
        // Shape function interpolation (Macro)
        //
        // Linear interpolation for triangular elements elements, e.g., for a 
        // variable u, interpolant takes the form: 
        //
        //             alpha_u[0]*x + alpha_u[1]*y + alpha_u[2]
        //
        // where x and y are the coordinates.
        //
        // Bilinear interpolation for quad elements elements, e.g., for a 
        // variable u, interpolant takes the form: 
        //
        //      alpha_u[0]*x + alpha_u[1]*y + alpha_u[2]*x*y + alpha_u[3]
        //
        // where x and y are the coordinates.
        //
        // h - Pressure head
        for (j=0; j<M_no_verts; j++)
        {
            alpha_h[j] = 0.0;
            
            for (i=0; i<M_no_verts; i++)
            {
                alpha_h[j] += h[M_icv[i]] * \
                        M_shape_funcs[(j*M_no_verts+i)*M_no_elements+k];
            }
        }
        // Keff - Effective conductivity
        for (j=0; j<M_no_verts; j++)
        {
            alpha_Kvec11[j] = 0.0;
            alpha_Kvec12[j] = 0.0;
            alpha_Kvec21[j] = 0.0;
            alpha_Kvec22[j] = 0.0;
            
            for (i=0; i<M_no_verts; i++)
            {
                alpha_Kvec11[j] += Kvec11[i] * \
                        M_shape_funcs[(j*M_no_verts+i)*M_no_elements+k];
                alpha_Kvec12[j] += Kvec12[i] * 
                        M_shape_funcs[(j*M_no_verts+i)*M_no_elements+k];
                alpha_Kvec21[j] += Kvec21[i] * \
                        M_shape_funcs[(j*M_no_verts+i)*M_no_elements+k];
                alpha_Kvec22[j] += Kvec22[i] * \
                        M_shape_funcs[(j*M_no_verts+i)*M_no_elements+k];
            }
        }
        //----------------------------------------------------------------------
        
        // Loop over edges within element
        for (j=0; j<M_no_verts; j++)
        {
            
            indx = j*2;
            jnb  = (j+1) % M_no_verts;
            
            //------------------------------------------------------------------
            // Finite Element Linear Interpolation
            //
            // Midpoints of edge j
            midx = M_midx_mat[j*M_no_elements+k];
            midy = M_midy_mat[j*M_no_elements+k];
            // Triangular elements
            if (M_no_verts == 3)
            {
                // Keff approximated at midpoint of edge j
                K_av11 = alpha_Kvec11[0]*midx + alpha_Kvec11[1]*midy \
                        + alpha_Kvec11[2];
                K_av12 = alpha_Kvec12[0]*midx + alpha_Kvec12[1]*midy \
                        + alpha_Kvec12[2];
                K_av21 = alpha_Kvec21[0]*midx + alpha_Kvec21[1]*midy \
                        + alpha_Kvec21[2];
                K_av22 = alpha_Kvec22[0]*midx + alpha_Kvec22[1]*midy \
                        + alpha_Kvec22[2];
                // Gradient of pressure head (h)
                grad_hx = alpha_h[0];
                grad_hy = alpha_h[1];
            }
            // Quad elements
            else if (M_no_verts == 4)
            {
                // Keff approximated at midpoint of edge j
                K_av11 = alpha_Kvec11[0]*midx + alpha_Kvec11[1]*midy \
                        + alpha_Kvec11[2]*midx*midy + alpha_Kvec11[3];
                K_av12 = alpha_Kvec12[0]*midx + alpha_Kvec12[1]*midy \
                        + alpha_Kvec12[2]*midx*midy + alpha_Kvec12[3];
                K_av21 = alpha_Kvec21[0]*midx + alpha_Kvec21[1]*midy \
                        + alpha_Kvec21[2]*midx*midy + alpha_Kvec21[3];
                K_av22 = alpha_Kvec22[0]*midx + alpha_Kvec22[1]*midy \
                        + alpha_Kvec22[2]*midx*midy + alpha_Kvec22[3];                
                // Gradient of pressure head (h)
                grad_hx = alpha_h[0] + alpha_h[2]*midy;
                grad_hy = alpha_h[1] + alpha_h[2]*midx;
            }
            
            // Macro flux vector
            q[0] = -(K_av11 * grad_hx + K_av12 * (grad_hy+1.0));
            q[1] = -(K_av21 * grad_hx + K_av22 * (grad_hy+1.0));
            
            // Vector normal to edge j (directed towards node M_icv[jnb])
            n[0] = M_normals[indx*M_no_elements+k];
            n[1] = M_normals[(indx+1)*M_no_elements+k];
            
            // Dot product
            qn = n[0] * q[0] + n[1] * q[1];
            
            // Distribute flux at edge j to appropriate nodes
            varj   = M_elements[j*M_no_elements+k];
            varjnb = M_elements[jnb*M_no_elements+k];
            garray[thread_num*N+varj] += - qn;
            garray[thread_num*N+varjnb] += qn;
            
            //------------------------------------------------------------------
            // Process Boundary conditions
            //
            // Flux specified
            if (M_boundary_condition_edges[j*M_no_elements+k] == 1.0)
            {
                garray[thread_num*N+varj]   += - boundary_flux * \
                        M_boundary_edge_lens[indx*M_no_elements+k];
                garray[thread_num*N+varjnb] += - boundary_flux * \
                        M_boundary_edge_lens[(indx+1)*M_no_elements+k];
            }
            
            // No flux
            // *Nothing to do*    
            
            // Free drainage (dh/dy = 0 or dh/dx = 0)
            if (M_boundary_condition_edges[j*M_no_elements+k] == 3.0 || 
                    M_boundary_condition_edges[j*M_no_elements+k] == 4.0)
            {
                midx = M_boundary_edge_midpts[4*j*M_no_elements+k];
                midy = M_boundary_edge_midpts[(4*j+1)*M_no_elements+k];
                if (M_no_verts == 3)
                {
                    K_av11 = alpha_Kvec11[0]*midx + alpha_Kvec11[1]*midy \
                            + alpha_Kvec11[2];
                    K_av12 = alpha_Kvec12[0]*midx + alpha_Kvec12[1]*midy \
                            + alpha_Kvec12[2];
                    K_av21 = alpha_Kvec21[0]*midx + alpha_Kvec21[1]*midy \
                            + alpha_Kvec21[2];
                    K_av22 = alpha_Kvec22[0]*midx + alpha_Kvec22[1]*midy \
                            + alpha_Kvec22[2];
                }
                else if (M_no_verts == 4)
                {
                    K_av11 = alpha_Kvec11[0]*midx + alpha_Kvec11[1]*midy \
                            + alpha_Kvec11[2]*midx*midy + alpha_Kvec11[3];
                    K_av12 = alpha_Kvec12[0]*midx + alpha_Kvec12[1]*midy \
                            + alpha_Kvec12[2]*midx*midy + alpha_Kvec12[3];
                    K_av21 = alpha_Kvec21[0]*midx + alpha_Kvec21[1]*midy \
                            + alpha_Kvec21[2]*midx*midy + alpha_Kvec21[3];
                    K_av22 = alpha_Kvec22[0]*midx + alpha_Kvec22[1]*midy \
                            + alpha_Kvec22[2]*midx*midy + alpha_Kvec22[3];
                    grad_hx = alpha_h[0] + alpha_h[2]*midy;
                    grad_hy = alpha_h[1] + alpha_h[2]*midx;
                }
                if (M_boundary_condition_edges[j*M_no_elements+k] == 3.0)
                {
                    // grad_hx = 0
                    q[0] = -(K_av11 * 0.0 + K_av12 * (grad_hy + 1.0));
                    q[1] = -(K_av21 * 0.0 + K_av22 * (grad_hy + 1.0));
                }
                else if (M_boundary_condition_edges[j*M_no_elements+k] == 4.0)
                {
                    // grad_hy = 0
                    q[0] = -(K_av11 * grad_hx + K_av12 * (0.0 + 1.0));
                    q[1] = -(K_av21 * grad_hx + K_av22 * (0.0 + 1.0));
                }
                n[0] = M_boundary_edge_normals[(4*j)*M_no_elements+k];
                n[1] = M_boundary_edge_normals[(4*j+1)*M_no_elements+k];
                bflux = q[0] * n[0] + q[1] * n[1];
                garray[thread_num*N+varj] += - bflux;
                
                midx = M_boundary_edge_midpts[(4*j+2)*M_no_elements+k];
                midy = M_boundary_edge_midpts[(4*j+3)*M_no_elements+k];
                if (M_no_verts == 3)
                {
                    K_av11 = alpha_Kvec11[0]*midx + alpha_Kvec11[1]*midy \
                            + alpha_Kvec11[2];
                    K_av12 = alpha_Kvec12[0]*midx + alpha_Kvec12[1]*midy \
                            + alpha_Kvec12[2];
                    K_av21 = alpha_Kvec21[0]*midx + alpha_Kvec21[1]*midy \
                            + alpha_Kvec21[2];
                    K_av22 = alpha_Kvec22[0]*midx + alpha_Kvec22[1]*midy \
                            + alpha_Kvec22[2];
                }
                else if (M_no_verts == 4)
                {
                    K_av11 = alpha_Kvec11[0]*midx + alpha_Kvec11[1]*midy \
                            + alpha_Kvec11[2]*midx*midy + alpha_Kvec11[3];
                    K_av12 = alpha_Kvec12[0]*midx + alpha_Kvec12[1]*midy \
                            + alpha_Kvec12[2]*midx*midy + alpha_Kvec12[3];
                    K_av21 = alpha_Kvec21[0]*midx + alpha_Kvec21[1]*midy \
                            + alpha_Kvec21[2]*midx*midy + alpha_Kvec21[3];
                    K_av22 = alpha_Kvec22[0]*midx + alpha_Kvec22[1]*midy \
                            + alpha_Kvec22[2]*midx*midy + alpha_Kvec22[3];
                    grad_hx = alpha_h[0] + alpha_h[2]*midy;
                    grad_hy = alpha_h[1] + alpha_h[2]*midx;
                }
                if (M_boundary_condition_edges[j*M_no_elements+k] == 3.0)
                {
                    // grad_hx = 0
                    q[0] = -(K_av11 * 0.0 + K_av12 * (grad_hy + 1.0));
                    q[1] = -(K_av21 * 0.0 + K_av22 * (grad_hy + 1.0));
                }
                else if (M_boundary_condition_edges[j*M_no_elements+k] == 4.0)
                {
                    // grad_hy = 0
                    q[0] = -(K_av11 * grad_hx + K_av12 * (0.0 + 1.0));
                    q[1] = -(K_av21 * grad_hx + K_av22 * (0.0 + 1.0));
                }
                n[0] = M_boundary_edge_normals[(4*j+2)*M_no_elements+k];
                n[1] = M_boundary_edge_normals[(4*j+3)*M_no_elements+k];
                bflux = q[0] * n[0] + q[1] * n[1];
                garray[thread_num*N+varjnb] += - bflux;
                
            }            
        }
        
        //----------------------------------------------------------------------
        // Values of h at micro nodes for macro element k
        //
        // Centroid of macro element k
        centx = M_centroid_elements[k];
        centy = M_centroid_elements[M_no_elements + k];                
        int node_indx;

        for (j=0; j<m_no_nodes; j++)
        {             
            node_indx = M_no_nodes + k*m_no_nodes + j;
            
            if (m_boundary_nodes[j] == 1)
            {
                h[node_indx] = u[node_indx];

                if (model_type == 1)
                {
                    x_coord = centx;
                    y_coord = centy;
                }
                else
                {
                    x_coord = centx + m_vectors[j];
                    y_coord = centy + m_vectors[m_no_nodes + j];                    
                }                
                if (M_no_verts == 3)
                {
                    hb = alpha_h[0]*x_coord + alpha_h[1]*y_coord \
                            + alpha_h[2];
                }
                else if (M_no_verts == 4)
                {
                    hb = alpha_h[0]*x_coord + alpha_h[1]*y_coord \
                            + alpha_h[2]*x_coord*y_coord + alpha_h[3];
                }
                
                f[node_indx] = 1.0;
                garray[thread_num*N+node_indx] = 0.0;
                h[node_indx] = hb;
            }
            else
            {
                h[node_indx] = u[node_indx];
                h_node = h[node_indx];
                Se     = pow(1.0 + pow(-alphaB*h_node,nB),-1.0 * mB);
                C_node = (thetarB - thetasB) * Se * mB * \
                        pow(-alphaB*h_node,nB) * nB / \
                        (h_node * (1.0 + pow(-alphaB*h_node,nB)));
                f[node_indx] = C_node * m_CV_area[j];
            }
        }
        //----------------------------------------------------------------------         
        
        double inclusion_flux[2];
        inclusion_flux[0] = 0.0;
        inclusion_flux[1] = 0.0;
        
        //----------------------------------------------------------------------
        // Microscopic fluxes: Loop over macro elements and compute fluxes
        for (i=0; i<m_no_elements; i++)
        {
            unsigned int m_icv[m_no_verts];
            unsigned int varj;
            unsigned int varjnb;
            double n[2];
            double q[2];
            double qn;
            unsigned int indx;
            unsigned int jnb;
            double Kvec[m_no_verts];
            double Cvec[m_no_verts];
            double Se;
            double alpha_h[m_no_verts];
            double alpha_Kvec[m_no_verts];
            double midx;
            double midy;
            double K_av;
            
            // Loop over micro nodes and compute properties
            for (j=0; j<m_no_verts; j++)
            {
                m_icv[j] = M_no_nodes + \
                        k*m_no_nodes + m_elements[j*m_no_elements+i];
                h_node  = h[m_icv[j]];
                Se      = pow(1.0 + pow(-alphaB*h_node,nB),-1.0 * mB);
                Kvec[j] = KsatB * pow(Se,0.5) * \
                        pow(1.0 - pow(1.0 - pow(Se,1.0/mB),mB),2.0); 
            }
            
            //------------------------------------------------------------------
            // Shape function interpolation (Micro)
            //
            // Linear interpolation for triangular elements elements, e.g., for
            // a variable u, interpolant takes the form:
            //
            //             alpha_u[0]*x + alpha_u[1]*y + alpha_u[2]
            //
            // where x and y are the coordinates.
            //
            // Bilinear interpolation for quad elements elements, e.g., for a
            // variable u, interpolant takes the form:
            //
            //      alpha_u[0]*x + alpha_u[1]*y + alpha_u[2]*x*y + alpha_u[3]
            //
            // where x and y are the coordinates.
            //
            // h - Pressure head
            for (j=0; j<m_no_verts; j++)
            {
                alpha_Kvec[j] = 0.0;
                alpha_h[j] = 0.0;
                
                for (p=0; p<m_no_verts; p++)
                {
                    alpha_Kvec[j] += Kvec[p] * \
                            m_shape_funcs[(j*m_no_verts+p)*m_no_elements+i];
                    alpha_h[j] += h[m_icv[p]] * \
                            m_shape_funcs[(j*m_no_verts+p)*m_no_elements+i];
                }
            }
            
            //------------------------------------------------------------------
            // Calculate micro flux at centroid of micro element
            midx = m_centroid_elements[i];
            midy = m_centroid_elements[m_no_elements+i];
            if (m_no_verts == 3)
            {
                K_av = alpha_Kvec[0]*midx + alpha_Kvec[1]*midy + alpha_Kvec[2];
                grad_hx = alpha_h[0];
                grad_hy = alpha_h[1];
            }
            else if (m_no_verts == 4)
            {
                K_av = alpha_Kvec[0]*midx + alpha_Kvec[1]*midy + \
                        alpha_Kvec[2]*midx*midy + alpha_Kvec[3];
                grad_hx = alpha_h[0] + alpha_h[2]*midy;
                grad_hy = alpha_h[1] + alpha_h[2]*midx;
            }
            
            double qmicro[2];
            qmicro[0] = -K_av * grad_hx;
            qmicro[1] = -K_av * (grad_hy + 1.0);
            inclusion_flux[0] += qmicro[0] * m_element_area[i];
            inclusion_flux[1] += qmicro[1] * m_element_area[i];
            //------------------------------------------------------------------
            
            // Loop over edges within element
            for (j=0; j<m_no_verts; j++)
            {
                
                indx = j*2;
                jnb  = (j+1) % m_no_verts;
                
                // Finite Element Linear Interpolation
                midx = m_midx_mat[j*m_no_elements+i];
                midy = m_midy_mat[j*m_no_elements+i];
                
                if (m_no_verts == 3)
                {
                    K_av = alpha_Kvec[0]*midx + alpha_Kvec[1]*midy + \
                            alpha_Kvec[2];
                    grad_hx = alpha_h[0];
                    grad_hy = alpha_h[1];
                }
                else if (m_no_verts == 4)
                {
                    K_av = alpha_Kvec[0]*midx + alpha_Kvec[1]*midy + \
                            alpha_Kvec[2]*midx*midy + alpha_Kvec[3];
                    grad_hx = alpha_h[0] + alpha_h[2]*midy;
                    grad_hy = alpha_h[1] + alpha_h[2]*midx;
                }
                
                // Micro flux vector
                q[0] = -K_av * grad_hx;
                q[1] = -K_av * (grad_hy+1.0);

                // Vector normal to edge j (directed towards node m_icv[jnb])                
                n[0] = m_normals[indx*m_no_elements+i];
                n[1] = m_normals[(indx+1)*m_no_elements+i];
                
                // Dot product
                qn   = n[0] * q[0] + n[1] * q[1];
                
                varj   = m_elements[j*m_no_elements+i];
                varjnb = m_elements[jnb*m_no_elements+i];
                
                unsigned int varj_indx = m_icv[j]; 
                unsigned int varjnb_indx = m_icv[jnb];
                
                // Distribute flux at edge j to appropriate nodes
                // (Interior nodes only)                
                if (m_boundary_nodes[varj] == 0)
                {
                    garray[thread_num*N+varj_indx] += - qn;
                }
                
                if (m_boundary_nodes[varjnb] == 0)
                {
                    garray[thread_num*N+varjnb_indx] += qn;
                }
                
                // Process interface elements for macro source term
                if (m_boundary_edges[j*m_no_elements+i] == 1.0)
                {
                    
                    midx = m_boundary_edge_midpts[4*j*m_no_elements+i];
                    midy = m_boundary_edge_midpts[(4*j+1)*m_no_elements+i];
                    if (m_no_verts == 3)
                    {
                        K_av = alpha_Kvec[0]*midx + alpha_Kvec[1]*midy \
                                + alpha_Kvec[2];
                        grad_hx = alpha_h[0];
                        grad_hy = alpha_h[1];
                    }
                    else if (m_no_verts == 4)
                    {
                        K_av = alpha_Kvec[0]*midx + alpha_Kvec[1]*midy \
                                + alpha_Kvec[2]*midx*midy + alpha_Kvec[3];
                        grad_hx = alpha_h[0] + alpha_h[2]*midy;
                        grad_hy = alpha_h[1] + alpha_h[2]*midx;
                    }
                    q[0] = -K_av * grad_hx;
                    q[1] = -K_av * (grad_hy+1.0);
                    n[0] = m_boundary_edge_normals[4*j*m_no_elements+i];
                    n[1] = m_boundary_edge_normals[(4*j+1)*m_no_elements+i];
                    qn   = n[0] * q[0] + n[1] * q[1];
                    for (p=0; p<M_no_verts; p++)
                    {
                        garray[thread_num*N+M_icv[p]] += qn * \
                                M_SCV_area[p*M_no_elements+k] / m_cell_area;
                    }
                    
                    midx = m_boundary_edge_midpts[(4*j+2)*m_no_elements+i];
                    midy = m_boundary_edge_midpts[(4*j+3)*m_no_elements+i];
                    if (m_no_verts == 3)
                    {
                        K_av = alpha_Kvec[0]*midx + alpha_Kvec[1]*midy + \
                                alpha_Kvec[2];
                        grad_hx = alpha_h[0];
                        grad_hy = alpha_h[1];
                    }
                    else if (m_no_verts == 4)
                    {
                        K_av = alpha_Kvec[0]*midx + alpha_Kvec[1]*midy + \
                                alpha_Kvec[2]*midx*midy + alpha_Kvec[3];
                        grad_hx = alpha_h[0] + alpha_h[2]*midy;
                        grad_hy = alpha_h[1] + alpha_h[2]*midx;
                    }
                    q[0] = -K_av * grad_hx;
                    q[1] = -K_av * (grad_hy+1.0);
                    n[0] = m_boundary_edge_normals[(4*j+2)*m_no_elements+i];
                    n[1] = m_boundary_edge_normals[(4*j+3)*m_no_elements+i];
                    qn   = n[0] * q[0] + n[1] * q[1];
                    
                    for (p=0; p<M_no_verts; p++)
                    {
                        garray[thread_num*N+M_icv[p]] += qn * \
                                M_SCV_area[p*M_no_elements+k] / (m_cell_area);
                    }
                    
                }                
            }
        }
        
        // Compute addition contribution to macro flux (inclusion flux)        
        if (model_type == 3)
        {
            inclusion_flux[0] = inclusion_flux[0] / m_cell_area;
            inclusion_flux[1] = inclusion_flux[1] / m_cell_area;
            
            for (j=0; j<M_no_verts; j++)
            {
                indx      = j*2;
                jnb       = (j+1) % M_no_verts;
                n[0]      = M_normals[indx*M_no_elements+k];
                n[1]      = M_normals[(indx+1)*M_no_elements+k];
                qn        = n[0] * inclusion_flux[0] + n[1] * inclusion_flux[1];
                varj      = M_elements[j*M_no_elements+k];
                varjnb    = M_elements[jnb*M_no_elements+k];
                garray[thread_num*N+varj] += - qn;
                garray[thread_num*N+varjnb] += qn;
            }
        }
        
    }
    //--------------------------------------------------------------------------
    
    // Add up g up across threads
    for (j=0; j<NUM_THREADS; j++)
    {
        for (i=0; i<N; i++)
        {
            g[i] += garray[j*N+i];
        }
    }
    
    // Scaling (chain rule to get in required form of du/dt = g(u)
    #pragma omp parallel for
    for (i=0; i<N; i++)
    {
        g[i] = g[i] / f[i];
    }    
    
    // Free memory allocations
    delete[] h;
    delete[] f;
    delete[] garray;
    
    return 0;

}
