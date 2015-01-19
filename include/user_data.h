struct user_data_struct
{
    macro_mesh_struct macro_mesh;
    micro_mesh_struct micro_mesh;
    Keff_struct effective_conductivity;
    soil_struct soilA;
    soil_struct soilB;
    double boundary_flux;
    int N;
    int model_type;
    unsigned int NUM_THREADS;
};