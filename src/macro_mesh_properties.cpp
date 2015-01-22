//------------------------------------------------------------------------------
// macro_mesh_properties.cpp
//
// Written by Elliot Carr
// Queensland University of Technology
//
// This code is part of TwoScalRich.
//
// This file reads in the macro mesh and generates properties such as CV_areas,
// edge lengths, etc.
// 
//------------------------------------------------------------------------------

#include "macro_mesh_properties.h"

int macro_mesh_properties(const char* mesh_file, \
        macro_mesh_struct &macro_mesh, int* macro_bc)
{
    
    FILE *myfile;
    
    unsigned int i, j;
    
    int    no_nodes, no_elements, no_objects, no_boundary_edges;
    double *nodes;
    int    *elements;
    double *normals;
    double *midx;
    double *midy;
    double *boundary_edge_lens;
    double *shape_funcs;
    int    *boundary_nodes_type;
    int    *boundary_edges;
    int    *boundary_condition_edges;
    double *CV_area;
    double *SCV_area;
    double *boundary_edge_midpts;
    double *boundary_edge_normals;
    int    *variable_number;
    double *centroid_elements;
    
    myfile = fopen(mesh_file,"r");
    
    if (myfile == NULL)
    {
        return 1;
    }
    
    char str_temp[80];
    
    int node_number;
    long double node_x;
    long double node_y;
    int no_verts;
    
    // mesh type (3 for triangles, 4 for rectangles)
    fscanf(myfile, "%i", &no_verts); 
    
    // Number of nodes
    fscanf(myfile, "%i", &no_nodes); 
    
        
    //--------------------------------------------------------------------------
    // Build nodes array
    nodes = new double [2*no_nodes];   
    for (i=0; i<no_nodes; i++)
    {        
        fscanf(myfile, "%Lf", &node_x);
        fscanf(myfile, "%Lf", &node_y);
        
        // Node coordinates array
        nodes[i] = node_x;
        nodes[no_nodes+i] = node_y;        
    }
    //--------------------------------------------------------------------------
    
    boundary_nodes_type = new int [5*no_nodes];
    for (i=0; i<no_nodes; i++)
    {        
        for (j=0; j<5; j++)
        {
            boundary_nodes_type[j*no_nodes+i] = 0;
        }
    }
    
    // Number of boundary edges    
    fscanf(myfile, "%i", &no_boundary_edges); 
    
    for (i=0; i<no_boundary_edges; i++)
    {
        int temp1, temp2, temp3;
        fscanf(myfile, "%i", &temp1);
        fscanf(myfile, "%i", &temp2);
        fscanf(myfile, "%i", &temp3);
        
        for (j=0; j<5; j++)
        {
            if (temp3 == -(j+1))
            {
                boundary_nodes_type[j*no_nodes + temp1-1] = 1;
                boundary_nodes_type[j*no_nodes + temp2-1] = 1;
            }
        }
        
    }
    
    //--------------------------------------------------------------------------
    // Build elements array
    fscanf(myfile, "%i", &no_elements); // Number of elements
    elements = new int [no_verts*no_elements];
    for (i=0; i<no_elements; i++)
    {
        if (no_verts == 3)
        {
            int temp1, temp2, temp3;
            fscanf(myfile, "%i", &temp1);
            fscanf(myfile, "%i", &temp2);
            fscanf(myfile, "%i", &temp3);
            
            elements[i] = temp1-1;
            elements[no_elements+i] = temp2-1;
            elements[2*no_elements+i] = temp3-1;
        }
        else if (no_verts == 4)
        {
            int temp1, temp2, temp3, temp4;
            fscanf(myfile, "%i", &temp1);
            fscanf(myfile, "%i", &temp2);
            fscanf(myfile, "%i", &temp3);
            fscanf(myfile, "%i", &temp4);
            
            elements[i] = temp1-1;
            elements[no_elements+i] = temp2-1;
            elements[2*no_elements+i] = temp3-1;
            elements[3*no_elements+i] = temp4-1;
        }      
    }
    //--------------------------------------------------------------------------
    
    fclose(myfile);
    
    //--------------------------------------------------------------------------
    // Boundary edge arrays
    //
    // boundary_edges[p*no_elements+i] = 1 if "i"th edge of element p is a 
    // boundary face (0 otherwise)
    // boundary_condition_edges[p*no_elements+i] = j if "i"th edge of element p
    // is a boundary edge with boundary condition of type j = 1,2,3,4, where:
    // 1: Flux specified (q.n = boundary_flux)
    // 2: No flux (q.n = 0)
    // 3: dh/dx = 0 (Free drainage condition for East and West boundaries)
    // 4: dh/dy = 0 (Free drainage condition for North and South boundaries)
    boundary_edges = new int [no_verts*no_elements];
    boundary_condition_edges = new int [no_verts*no_elements];
    
    for (i=0; i<no_elements; i++)
    {
        
        for (j=0; j<no_verts; j++)
        {              
            boundary_edges[j*no_elements+i] = 0;
        }
        
        for (j=0; j<no_verts; j++)
        {              
            boundary_condition_edges[j*no_elements+i] = 0;
        }
        
    }    

    for (i=0; i<no_elements; i++)
    {
        for (j=0; j<5; j++)
        {
            for (int p=0; p<no_verts; p++)
            {
                int pnb = (p+1) % no_verts;
                int varj = j*no_nodes + elements[p*no_elements+i];
                int varjnb = j*no_nodes + elements[pnb*no_elements+i];
                
                if (boundary_nodes_type[varj] == 1 && \
                        boundary_nodes_type[varjnb] == 1)
                {
                    // Type of BC for boundary edge
                    boundary_condition_edges[p*no_elements+i] = macro_bc[j]; 
                    boundary_edges[p*no_elements+i] = 1;
                }
            }
        }
        
    }
    //--------------------------------------------------------------------------
    
    boundary_edge_lens    = new double [2*no_verts*no_elements];
    boundary_edge_normals = new double [4*no_verts*no_elements];
    boundary_edge_midpts  = new double [4*no_verts*no_elements];
    
    for (i=0; i<no_elements; i++)
    {        
        for (j=0; j<2*no_verts; j++)
        {              
            boundary_edge_lens[j*no_elements+i] = 0.0;
        }
        for (j=0; j<4*no_verts; j++)
        {
            boundary_edge_midpts[j*no_elements+i] = 0.0;
            boundary_edge_normals[j*no_elements+i] = 0.0;
        }        
    }
    
    //--------------------------------------------------------------------------
    // Generate mesh properties such as CV_areas, edge lengths, etc.
    int vert[no_verts];
    double xv[no_verts];
    double yv[no_verts];
    double centroid[2];
    double* midpts = new double [2*no_verts];
    double* p      = new double [2*no_verts];
    double* q      = new double [2*no_verts];
    double* f      = new double [2*no_verts];
    double detA;
    double ordering;
    int jnb;
    double* SCV_area_element = new double [no_verts];
    double a;
    
    normals			     = new double [2*no_verts*no_elements];
    CV_area		         = new double [no_nodes];
    midx			     = new double [no_verts*no_elements];
    midy			     = new double [no_verts*no_elements];
    shape_funcs		     = new double [no_verts*no_verts*no_elements];
    SCV_area             = new double [no_verts*no_elements];
    centroid_elements    = new double [2*no_elements];
    
    int no_variables = 0;
    
    for (i=0; i<no_nodes; i++)
    {
        CV_area[i] = 0.0;        
    }
    
    for (i=0; i<no_elements; i++)
    {
        // Centroid of element i
        centroid[0] = 0.0;
        centroid[1] = 0.0;
        for (j=0; j<no_verts; j++)
        {
            vert[j] = elements[j*no_elements+i];
            xv[j]   = nodes[vert[j]];
            yv[j]   = nodes[no_nodes+vert[j]];
            centroid[0] = centroid[0] + xv[j];
            centroid[1] = centroid[1] + yv[j];            
        }
        
        centroid[0] = centroid[0] / ((double) no_verts);
        centroid[1] = centroid[1] / ((double) no_verts);        
        centroid_elements[i] = centroid[0];
        centroid_elements[no_elements+i] = centroid[1];
        
        ordering = 0.0;
        for (j=0; j<no_verts; j++)
        {
            jnb = (j+1) % no_verts;
            ordering = ordering + (xv[jnb]-xv[j])*(yv[jnb]+yv[j]);
            midpts[2*j+0] = (xv[j] + xv[jnb]) / 2.0;
            midpts[2*j+1] = (yv[j] + yv[jnb]) / 2.0;            
        }
        
        // Control volume and sub-control volume areas
        if (no_verts == 3)
        {
            for (j=0; j<no_verts; j++)
            {
                jnb = (j+1) % no_verts;
                
                p[2*j+0] = centroid[0]   - xv[j];
                p[2*j+1] = centroid[1]   - yv[j];
                q[2*j+0] = midpts[2*j+0] - midpts[2*jnb+0];
                q[2*j+1] = midpts[2*j+1] - midpts[2*jnb+1];
            }
        }
            
        for (j=0; j<no_verts; j++)    
        {
            f[2*j+0] = centroid[0] - midpts[2*j+0];
            f[2*j+1] = centroid[1] - midpts[2*j+1];
        }
        
        for (j=0; j<no_verts; j++)
        {             
            if (no_verts == 3)
            {
                SCV_area_element[j] = 0.5*fabs(p[2*j+0]*q[2*j+1] -\
                        p[2*j+1]*q[2*j+0]);
                if (j == 0)
                {
                    a = p[2*j+0]*q[2*j+1] - p[2*j+1]*q[2*j+0];
                }
            }
            else if (no_verts == 4)
            {
                SCV_area_element[j] = 0.25*fabs((xv[2]-xv[0])*(yv[2]-yv[0]));
            }

            SCV_area[j*no_elements+i] = SCV_area_element[j];
            CV_area[vert[j]] = CV_area[vert[j]] + SCV_area_element[j];            
        }

        // Normals and midpoints of internal edges
        for (j=0; j<no_verts; j++)
        {
            if (no_verts == 3)
            {
                normals[(2*j+0)*no_elements+i] = -sign(a) * f[2*j+1];
                normals[(2*j+1)*no_elements+i] = -sign(a) * -f[2*j+0];
            }
            else if (no_verts == 4)
            {
                normals[(2*j+0)*no_elements+i] = f[2*j+1];
                normals[(2*j+1)*no_elements+i] = -f[2*j+0];
            }
            midx[j*no_elements+i] = (midpts[2*j+0] + centroid[0]) / 2.0;
            midy[j*no_elements+i] = (midpts[2*j+1] + centroid[1]) / 2.0;
        }
        
        // Midpoints, lengths and normals of boundary edges
        for (j=0; j<no_verts; j++)
        {
            if (boundary_edges[j*no_elements+i] == 1)
            {
                jnb = (j+1) % no_verts;
                
                boundary_edge_lens[(2*j+0)*no_elements+i] = \
                        two_norm(xv[j],yv[j],midpts[2*j+0],midpts[2*j+1]);
                boundary_edge_lens[(2*j+1)*no_elements+i] = \
                        two_norm(xv[jnb],yv[jnb],midpts[2*j+0],midpts[2*j+1]);
                
                boundary_edge_midpts[(4*j+0)*no_elements+i] = \
                        (nodes[vert[j]] + midpts[2*j+0]) / 2.0;
                boundary_edge_midpts[(4*j+1)*no_elements+i] = \
                        (nodes[no_nodes+vert[j]] + midpts[2*j+1]) / 2.0;
                boundary_edge_midpts[(4*j+2)*no_elements+i] = \
                        (nodes[vert[jnb]] + midpts[2*j+0]) / 2.0;
                boundary_edge_midpts[(4*j+3)*no_elements+i] = \
                        (nodes[no_nodes+vert[jnb]] + midpts[2*j+1]) / 2.0; 
                
                double bf1[2];
                double bf2[2];
                bf1[0] = midpts[2*j+0] - xv[j];
                bf1[1] = midpts[2*j+1] - yv[j];
                bf2[0] = xv[jnb] - midpts[2*j+0];
                bf2[1] = yv[jnb] - midpts[2*j+1];
                
                boundary_edge_normals[(4*j+0)*no_elements+i] = \
                        -sign(ordering) * bf1[1];
                boundary_edge_normals[(4*j+1)*no_elements+i] = \
                        -sign(ordering) * -bf1[0];
                boundary_edge_normals[(4*j+2)*no_elements+i] = \
                        -sign(ordering) * bf2[1];
                boundary_edge_normals[(4*j+3)*no_elements+i] = \
                        -sign(ordering) * -bf2[0];
            }
        }
        
        // Shape function interpolation coefficients
        if (no_verts == 3)
        {
            detA = xv[1]*yv[2] - xv[2]*yv[1] - xv[0]*yv[2] + xv[2]*yv[0] + \
                    xv[0]*yv[1] - xv[1]*yv[0];
            shape_funcs[i]               = (yv[1]-yv[2]) / detA;
            shape_funcs[no_elements+i]   = (yv[2]-yv[0]) / detA;
            shape_funcs[2*no_elements+i] = (yv[0]-yv[1]) / detA;
            shape_funcs[3*no_elements+i] = (xv[2]-xv[1]) / detA;
            shape_funcs[4*no_elements+i] = (xv[0]-xv[2]) / detA;
            shape_funcs[5*no_elements+i] = (xv[1]-xv[0]) / detA;
            shape_funcs[6*no_elements+i] = (xv[1]*yv[2] - yv[1]*xv[2]) / detA;
            shape_funcs[7*no_elements+i] = (yv[0]*xv[2] - xv[0]*yv[2]) / detA;
            shape_funcs[8*no_elements+i] = (xv[0]*yv[1] - yv[0]*xv[1]) / detA;            
        }
        else if (no_verts == 4)
        {
            detA = -xv[0]*xv[3]*yv[1]*yv[3] + xv[0]*xv[3]*yv[0]*yv[1] - \
                    xv[0]*xv[3]*yv[0]*yv[2] - xv[1]*xv[2]*yv[1]*yv[3] + \
                    xv[1]*xv[2]*yv[2]*yv[3] - xv[1]*xv[3]*yv[0]*yv[1] + \
                    xv[1]*xv[3]*yv[0]*yv[3] + xv[1]*xv[3]*yv[1]*yv[2] - \
                    xv[1]*xv[3]*yv[2]*yv[3] + xv[2]*xv[3]*yv[0]*yv[2] - \
                    xv[2]*xv[3]*yv[0]*yv[3] - xv[2]*xv[3]*yv[1]*yv[2] + \
                    xv[2]*xv[3]*yv[1]*yv[3] + xv[0]*xv[3]*yv[2]*yv[3] + \
                    xv[1]*xv[2]*yv[0]*yv[1] - xv[1]*xv[2]*yv[0]*yv[2] + \
                    xv[0]*xv[1]*yv[0]*yv[2] - xv[0]*xv[1]*yv[0]*yv[3] - \
                    xv[0]*xv[1]*yv[1]*yv[2] + xv[0]*xv[1]*yv[1]*yv[3] - \
                    xv[0]*xv[2]*yv[0]*yv[1] + xv[0]*xv[2]*yv[0]*yv[3] + \
                    xv[0]*xv[2]*yv[1]*yv[2] - xv[0]*xv[2]*yv[2]*yv[3];
            shape_funcs[i] = (xv[2]*yv[1]*yv[2] - xv[3]*yv[1]*yv[3] + \
                    xv[3]*yv[2]*yv[3] - xv[1]*yv[1]*yv[2] + \
                    xv[1]*yv[1]*yv[3] - xv[2]*yv[2]*yv[3]) / detA;
            shape_funcs[no_elements+i] = (-xv[2]*yv[0]*yv[2] + \
                    xv[3]*yv[0]*yv[3] - xv[3]*yv[2]*yv[3] + \
                    xv[0]*yv[0]*yv[2] - xv[0]*yv[0]*yv[3] + \
                    xv[2]*yv[2]*yv[3]) / detA;
            shape_funcs[2*no_elements+i] = (xv[1]*yv[0]*yv[1] - \
                    xv[3]*yv[0]*yv[3] + xv[3]*yv[1]*yv[3] - \
                    xv[0]*yv[0]*yv[1] + xv[0]*yv[0]*yv[3] - \
                    xv[1]*yv[1]*yv[3]) / detA;
            shape_funcs[3*no_elements+i]  = (-xv[1]*yv[0]*yv[1] + \
                    xv[2]*yv[0]*yv[2] - xv[2]*yv[1]*yv[2] + \
                    xv[0]*yv[0]*yv[1] - xv[0]*yv[0]*yv[2] + \
                    xv[1]*yv[1]*yv[2]) / detA;
            shape_funcs[4*no_elements+i]  = (-xv[2]*xv[3]*yv[3] + \
                    xv[2]*xv[3]*yv[2] + xv[1]*xv[3]*yv[3] - \
                    xv[3]*yv[1]*xv[1] - xv[1]*xv[2]*yv[2] + \
                    xv[2]*yv[1]*xv[1]) / detA;
            shape_funcs[5*no_elements+i]  = (xv[2]*xv[3]*yv[3] - \
                    xv[2]*xv[3]*yv[2] - xv[0]*xv[3]*yv[3] + \
                    xv[3]*yv[0]*xv[0] + xv[0]*xv[2]*yv[2] - \
                    xv[2]*yv[0]*xv[0]) / detA;
            shape_funcs[6*no_elements+i]  = (-xv[1]*xv[3]*yv[3] + \
                    xv[3]*yv[1]*xv[1] + xv[0]*xv[3]*yv[3] - \
                    xv[3]*yv[0]*xv[0] - xv[0]*yv[1]*xv[1] + \
                    xv[1]*yv[0]*xv[0]) / detA;
            shape_funcs[7*no_elements+i]  = (-xv[1]*yv[0]*xv[0] + \
                    xv[2]*yv[0]*xv[0] - xv[0]*xv[2]*yv[2] - \
                    xv[2]*yv[1]*xv[1] + xv[0]*yv[1]*xv[1] + \
                    xv[1]*xv[2]*yv[2]) / detA;
            shape_funcs[8*no_elements+i]  = (-xv[2]*yv[1] + xv[1]*yv[2] - \
                    xv[1]*yv[3] + xv[2]*yv[3] + xv[3]*yv[1] - xv[3]*yv[2]) \
                    / detA;
            shape_funcs[9*no_elements+i]  = (xv[0]*yv[3] - xv[2]*yv[3] - \
                    xv[3]*yv[0] + xv[2]*yv[0] - xv[0]*yv[2] + xv[3]*yv[2]) \
                    / detA;
            shape_funcs[10*no_elements+i] = (-xv[0]*yv[3] + xv[1]*yv[3] - \
                    xv[1]*yv[0] + xv[3]*yv[0] - xv[3]*yv[1] + xv[0]*yv[1]) \
                    / detA;
            shape_funcs[11*no_elements+i] = (-xv[0]*yv[1] + xv[0]*yv[2] + \
                    xv[1]*yv[0] - xv[1]*yv[2] - xv[2]*yv[0] + xv[2]*yv[1]) \
                    / detA;
            shape_funcs[12*no_elements+i] = (-xv[1]*xv[2]*yv[1]*yv[3] + \
                    xv[1]*xv[2]*yv[2]*yv[3] - xv[1]*xv[3]*yv[2]*yv[3] + \
                    xv[1]*xv[3]*yv[1]*yv[2] + xv[2]*xv[3]*yv[1]*yv[3] - \
                    xv[2]*xv[3]*yv[1]*yv[2]) / detA;
            shape_funcs[13*no_elements+i] = (xv[0]*xv[2]*yv[0]*yv[3] - \
                    xv[0]*xv[2]*yv[2]*yv[3] + xv[0]*xv[3]*yv[2]*yv[3] - \
                    xv[2]*xv[3]*yv[0]*yv[3] + xv[2]*xv[3]*yv[0]*yv[2] - \
                    xv[0]*xv[3]*yv[0]*yv[2]) / detA;
            shape_funcs[14*no_elements+i] = (xv[0]*xv[1]*yv[1]*yv[3] - \
                    xv[0]*xv[1]*yv[0]*yv[3] + xv[1]*xv[3]*yv[0]*yv[3] - \
                    xv[1]*xv[3]*yv[0]*yv[1] - xv[0]*xv[3]*yv[1]*yv[3] + \
                    xv[0]*xv[3]*yv[0]*yv[1]) / detA;
            shape_funcs[15*no_elements+i] = (xv[0]*xv[1]*yv[0]*yv[2] - \
                    xv[0]*xv[1]*yv[1]*yv[2] + xv[0]*xv[2]*yv[1]*yv[2] - \
                    xv[0]*xv[2]*yv[0]*yv[1] - xv[1]*xv[2]*yv[0]*yv[2] + \
                    xv[1]*xv[2]*yv[0]*yv[1]) / detA;
            
        }
        
    }
    //--------------------------------------------------------------------------
    
    macro_mesh.no_verts                 = no_verts;
    macro_mesh.no_nodes                 = no_nodes;
    macro_mesh.no_elements              = no_elements;
    macro_mesh.elements                 = elements;    
    macro_mesh.nodes                    = nodes;
    macro_mesh.shape_funcs              = shape_funcs;
    macro_mesh.normals                  = normals;
    macro_mesh.boundary_edge_midpts     = boundary_edge_midpts;
    macro_mesh.boundary_edge_normals    = boundary_edge_normals;
    macro_mesh.boundary_edge_lens       = boundary_edge_lens;
    macro_mesh.CV_area                  = CV_area;
    macro_mesh.SCV_area                 = SCV_area;
    macro_mesh.midx                     = midx;
    macro_mesh.midy                     = midy;
    macro_mesh.boundary_condition_edges = boundary_condition_edges;
    macro_mesh.centroid_elements        = centroid_elements;
    
    // Free memory allocations
    delete[] SCV_area_element;
    delete[] boundary_edges;
    
    return 0;

}

double two_norm(double x1, double y1, double x2, double y2)
{
    double norm;
    norm = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
    return norm;
}
