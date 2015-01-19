struct micro_mesh_struct
{
    int    no_nodes;
    int    no_elements;
    int    *elements;
    int    *boundary_nodes;
    double *nodes;
    double *shape_funcs;
    double *normals;
    int    *boundary_edges;
    double *boundary_edge_midpts;
    double *boundary_edge_normals;
    double *boundary_edge_lens;
    double *centroid_elements;
    double *CV_area;
    double *midx;
    double *midy;
    double cell_area;
    double *vectors;
    int    no_verts;
    double *element_area;
};