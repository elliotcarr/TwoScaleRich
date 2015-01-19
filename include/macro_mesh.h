struct macro_mesh_struct
{
    int    no_verts;
    int    no_nodes;
    int    no_elements;
    int    *elements;
    double *nodes;
    double *shape_funcs;
    double *normals;
    double *boundary_edge_midpts;
    double *boundary_edge_normals;
    double *boundary_edge_lens;
    int    *boundary_condition_edges;
    double *CV_area;
    double *SCV_area;
    double *midx;
    double *midy;
    double *centroid_elements;
};