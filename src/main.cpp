//------------------------------------------------------------------------------
// main.cpp
//
// Written by Dr Elliot Carr (2013-2014)
// École Centrale Paris and Queensland University of Technology
//
// This code is part of TwoScalRich.
//
// This is the main program file.
//
//------------------------------------------------------------------------------

#include "main.h"

int main (int argc, char * const argv[]) {
    
    int i, j, k; 
    int flag;
    time_t tic, toc;
    
    //--------------------------------------------------------------------------
    // START: User changeable parameters
    //--------------------------------------------------------------------------
    
    cout << "There are " << argc << " arguments:" << endl;
    
    // Loop through each argument and print its number and value
    for (int nArg=0; nArg < argc; nArg++)
    {
        cout << nArg << " " << argv[nArg] << endl;
    }
    
    // Problem configuration file
    const char* config_file = argv[1];
    
    // Read configuration file
    FILE * myfile;
    myfile = fopen(config_file,"r");
    
    if (myfile == NULL)
    {
        return 1;
    }
    
    // Mesh files
    char macro_mesh_file[90];
    fscanf(myfile, "%s", macro_mesh_file);
    char micro_mesh_file[90];
    fscanf(myfile, "%s", micro_mesh_file);
    
    // Solution output file
    char solution_file[90];
    fscanf(myfile, "%s", solution_file);
    
    // Simulation end time (hours)
    long double tend;
    fscanf(myfile, "%Lf", &tend);
         
            
       
    
    
    //const char* macro_mesh_file = argv[1]; //"macro_structured_20by20.msh";
    //const char* micro_mesh_file = "microB_unstructured_20by20.msh";      
    
    // Stores solution at time values specified in tspan array
    //string solution_file = "solution.txt";
    int no_solns = (int) tend;
    double* tspan = new double [no_solns+1];
    tspan[0] = no_solns; // Store number as first element
    for (i=1; i<no_solns+1; i++)
    {
        tspan[i] = (double)i*3600.0;
    }    
    
    // Solver statistics
    string statistics_file = "statistics.txt";
    
    // Uniform initial condition [m]
    long double h0;
    fscanf(myfile, "%Lf", &h0);
    
    // Print Solver information (true or false)
    bool PrintStepInfo = true;
    
    // Print Solver information every 'PrintStep' steps
    int PrintStep = 1;
    
    // Soil properties
    // A: connected sub-domain
    // B: disconnected sub-domin (inclusions)
    // ***Important***
    // Changing these parameters requires the effective conductivity to be
    // computed again (see Keff/run.m)
    soil_struct soilA;
    soil_struct soilB;
    
    long double a;
    fscanf(myfile, "%Lf", &a); // hydraulic conductivity  
    soilA.Ksat   = a; //0.044 / 3600.0;  
    fscanf(myfile, "%Lf", &a);
    soilA.thetar = a; // residual moisture content
    fscanf(myfile, "%Lf", &a);
    soilA.thetas = a; // saturated moisture content
    fscanf(myfile, "%Lf", &a);
    soilA.alpha = a; // alpha paramter (van Genuchten relations)
    fscanf(myfile, "%Lf", &a);
    soilA.n = a; // n parameter (van Genuchten relations)

    soilA.m = 1.0 - 1.0 / soilA.n; // m = 1 - 1/n (van Genuchten relations)
    
    fscanf(myfile, "%Lf", &a); // hydraulic conductivity
    soilB.Ksat   = a; //0.044 / 3600.0;  
    fscanf(myfile, "%Lf", &a);
    soilB.thetar = a; // residual moisture content
    fscanf(myfile, "%Lf", &a);
    soilB.thetas = a; // saturated moisture content
    fscanf(myfile, "%Lf", &a);
    soilB.alpha = a; // alpha paramter (van Genuchten relations)
    fscanf(myfile, "%Lf", &a);
    soilB.n = a; // n parameter (van Genuchten relations)
    
    fclose(myfile);

    soilB.m = 1.0 - 1.0 / soilB.n; // m = 1 - 1/n (van Genuchten relations)
    
    // Effective conductivity file
    const char* Keff_file = "Keff.txt";
    
    // Macro boundary conditions
    int *macro_bc = new int [5];
    
    //--------------------------------------------------------------------------
    // Set Macro boundary conditions (boundary faces labelled below)
    // 
    //              |- Influx -| North         
    //          **************************
    //          *                        *
    //          *                        *
    //          *                        * 
    //          *                        *
    //          *                        *
    //          *                        *
    //     West *                        * East
    //          *                        *
    //          *                        *
    //          *                        *
    //          *                        *
    //          *                        *    
    //          *                        *    
    //          **************************
    //                    South
    //
    // 1: Flux specified (q.n = boundary_flux)
    // 2: No flux (q.n = 0)
    // 3: dh/dx = 0 (Free drainage condition for East and West boundaries)
    // 4: dh/dy = 0 (Free drainage condition for North and South boundaries)
    
    // Boundary flux [m/s]
    double boundary_flux  = -0.01 / 3600.0;
    macro_bc[0] = 1; // Influx
    macro_bc[1] = 3; // South
    macro_bc[2] = 2; // East
    macro_bc[3] = 2; // West
    macro_bc[4] = 2; // North (excluding Influx segment)
    //--------------------------------------------------------------------------
    
    // Model type
    // 1: uniform micro boundary condition
    // 2: non-uniform micro boundary condition
    // 3: non-uniform micro boundary condition, inclusion flux included at 
    //    macro-scale    
    int model_type = 3;
    
    //--------------------------------------------------------------------------
    // Get environment variables defined in script files (.sh files)
    
    // Parallel or serial implementation of right-hand side function (Gfunc)
    int parallel = atoi(getenv("parallel"));
    int (*Gfunc_type)(double* g, double* u, user_data_struct user_data);
    if (parallel == 0)
    {
        Gfunc_type = Gfunc_serial;
    }
    else if (parallel == 1)
    {
        Gfunc_type = Gfunc_parallel;
    }
    else
    {
        cout << "Variable `parallel' defined in script .sh file must be " << \
                "either 0 or 1.\n";
        exit(EXIT_FAILURE);
    }    
        
    // Number of OpenMP Threads
    unsigned int NUM_THREADS = max(atoi(getenv("OMP_NUM_THREADS")), 1);    
    //--------------------------------------------------------------------------
    
    string mesh_type      = "structured_";
    string d = "_10by10.txt";
    
    //--------------------------------------------------------------------------
    // Read effective conductivity values 
    Keff_struct effective_conductivity;
    flag = read_effective_conductivity(Keff_file, effective_conductivity);
    //--------------------------------------------------------------------------
    
    //--------------------------------------------------------------------------
    // Process macro and micro meshes and compute mesh properties (normals, 
    // control volume areas, etc.) 
    
    // Macro mesh
    macro_mesh_struct macro_mesh;    
    flag = macro_mesh_properties(macro_mesh_file, macro_mesh, macro_bc);
    
    if (flag == 1)
    {
        cout << "Error opening macro mesh file. Simulation aborted.\n";
        return 0;
    }
    
    // Micro mesh
    micro_mesh_struct micro_mesh;
    flag = micro_mesh_properties(micro_mesh_file, micro_mesh);
    
    if (flag == 1)
    {
        cout << "Error opening micro mesh file. Simulation aborted.\n";
        return 0;
    }
    //--------------------------------------------------------------------------
    
    // EXPREM solver options
    exprem_options_struct exprem_options;    
    flag = set_exprem_options(exprem_options);
    
    // EXPREM solver stats
    exprem_stats_struct exprem_stats;
    
    // Number of unknowns
    int N = macro_mesh.no_nodes + macro_mesh.no_elements * micro_mesh.no_nodes;
    
    // Initial solution
    double* u0 = new double [N];    
    for (i=0; i<N; i++)
    {
        u0[i] = h0;       
    }
    
    // Store initial solution
    double* tspan_solns = new double [N*(no_solns+1)];
    for (i=0; i<N; i++)
    {
        tspan_solns[i] = u0[i];
    }
    
    int num_sims;
    string folder;
    
    num_sims = 1;
    folder = "Results/";

    // User data
    user_data_struct user_data;
    user_data.macro_mesh             = macro_mesh;
    user_data.micro_mesh             = micro_mesh;
    user_data.soilA                  = soilA;
    user_data.soilB                  = soilB;
    user_data.boundary_flux          = boundary_flux;
    user_data.NUM_THREADS            = NUM_THREADS;
    user_data.effective_conductivity = effective_conductivity;
    user_data.model_type             = model_type;
    user_data.N                      = N;
    
    ofstream stats_file;
    stats_file.open(statistics_file.c_str());
    
    if (!stats_file.is_open())
    {
        cout << "Error opening solution file\n";
        return 0;
    }
    
    // Simulation wall time
    double dtime = omp_get_wtime();
    
    cout << scientific;
    cout << std::setprecision(2);  
    
    //--------------------------------------------------------------------------
    // Time stepping         
    double t = 0.0;
    tend = tend * 3600.0; // Hours to seconds
    while (t < tend)
    {
        // exprem takes one internal step and then exits
        flag = exprem(Gfunc_type, N, t, tend, tspan, tspan_solns, u0, \
                user_data, exprem_options, exprem_stats);
        
        if (flag == -1)
        {
            return 0;
        }
        
        stats_file << setw(5) << exprem_stats.NumSteps << " ";
        stats_file << scientific;
        stats_file << std::setprecision(4);
        stats_file << t/3600.0 << " ";
        stats_file << exprem_stats.CurrentStep << " ";
        stats_file << exprem_stats.CurrentKryDim << "\n";
        
        // Print every 'PrintStep' steps
        if (exprem_stats.NumSteps % PrintStep == 0 && PrintStepInfo == true) 
        {
            cout << "NSTEPS: " << setw(5) << exprem_stats.NumSteps << " | ";
            cout << scientific;
            cout << std::setprecision(4) << "TIME [H]: " << t/3600.0 << " | ";
            cout << "STEPSIZE [S]: " << exprem_stats.CurrentStep << " | ";
            cout << "KRYDIM: " << exprem_stats.CurrentKryDim << "\n";
        }
        
    }
    //--------------------------------------------------------------------------
    
    stats_file.close();      
    
    // Simulation wall time
    dtime = omp_get_wtime() - dtime;
    
    cout << scientific;
    cout << fixed << "Runtime: " << dtime << " secs" << "  N: " << N << \
            "   FuncEvals: " << exprem_stats.FuncEvals << "\n";
    cout << "Simulation completed successfully.\n";
    
    // Save solution to file
    save_solution(solution_file, tspan, tspan_solns, no_solns, user_data);
    
    // Free memory allocations
    delete[] tspan;
    delete[] macro_bc;
    delete[] u0;
    delete[] tspan_solns;
    
    return 0;    
}
