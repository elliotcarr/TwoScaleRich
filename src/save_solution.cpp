//------------------------------------------------------------------------------
// save_solution.cpp
//
// Written by Dr Elliot Carr (2013-2014)
// Ecole Centrale Paris and Queensland University of Technology
//
// This code is part of TwoScalRich.
//
// This file saves the solution to file.
//
//------------------------------------------------------------------------------

#include "save_solution.h"

void save_solution(string solution_file, double* tspan, double* tspan_solns, 
        int no_solns, user_data_struct user_data)
{
    int i, j, k, m;
    
    macro_mesh_struct macro_mesh;
    macro_mesh = user_data.macro_mesh;
    
    int     M_no_nodes          = macro_mesh.no_nodes;
    int     M_no_elements       = macro_mesh.no_elements;
    int     M_no_verts          = macro_mesh.no_verts;
    int*    M_elements          = macro_mesh.elements;
    double* M_shape_funcs       = macro_mesh.shape_funcs;
    double* M_centroid_elements = macro_mesh.centroid_elements;
    
    micro_mesh_struct micro_mesh;
    micro_mesh = user_data.micro_mesh;
    int     m_no_nodes       = micro_mesh.no_nodes;
    int*    m_boundary_nodes = micro_mesh.boundary_nodes;
    double* m_vectors        = micro_mesh.vectors;
    
    // Total number of unknowns
    int N = user_data.N;
    
    int model_type = user_data.model_type;
    
    // Calibrate solution for Dirichlet condition
    for (j=0; j<(no_solns+1); j++)
    {
        for (k=0; k<M_no_elements; k++)
        {
            int M_icv[M_no_verts];
            double alpha_h[M_no_verts];
            
            double centx = M_centroid_elements[k];
            double centy = M_centroid_elements[M_no_elements + k];
            
            for (m=0; m<M_no_verts; m++)
            {
                M_icv[m]  = M_elements[m*M_no_elements+k];
            }
            
            /* Build interpolant for h for gradient approximations */
            for (m=0; m<M_no_verts; m++)
            {
                alpha_h[m] = 0.0;
                
                for (i=0; i<M_no_verts; i++)
                {
                    alpha_h[m] = alpha_h[m] + tspan_solns[j*N + M_icv[i]] * \
                            M_shape_funcs[(m*M_no_verts+i)*M_no_elements+k];
                    
                }
            }
            
            for (i=0; i<m_no_nodes; i++)
            {
                int node_indx = M_no_nodes + k*m_no_nodes + i;
                
                if (m_boundary_nodes[i] == 1)
                {
                    double x_coord;
                    double y_coord;
                    
                    if (model_type == 1)
                    {
                        x_coord = centx;
                        y_coord = centy;
                    }
                    else
                    {
                        x_coord = centx + m_vectors[i];
                        y_coord = centy + m_vectors[m_no_nodes + i];
                    }
                    if (M_no_verts == 3)
                    {
                        tspan_solns[j*N + node_indx] = alpha_h[0]*x_coord + \
                                alpha_h[1]*y_coord + alpha_h[2];
                    }
                    else if (M_no_verts == 4)
                    {
                        tspan_solns[j*N + node_indx] = alpha_h[0]*x_coord + \
                                alpha_h[1]*y_coord + \
                                alpha_h[2]*x_coord*y_coord + alpha_h[3];
                    }
                }
            }
        }
    }
    
    //--------------------------------------------------------------------------
    // Save solution at tspan values
    ofstream temp_file;
    temp_file.open(solution_file.c_str());
    
    if (!temp_file.is_open())
    {
        cout << "Error opening solution file\n";
    }
    
    temp_file << no_solns+1 << "\n";
    temp_file << scientific;
    temp_file << std::setprecision(20);
    
    for (j=0; j<(no_solns+1); j++)
    {
        if (j == 0)
        {
            temp_file << 0.0;
        }
        else
        {
            temp_file << tspan[j];
        }
        
        for (k=0; k<N; k++)
        {
            temp_file << " ";
            temp_file << tspan_solns[j*N + k];
        }
        
        temp_file << "\n";
        
    }
    temp_file.close();
    //--------------------------------------------------------------------------
    
}