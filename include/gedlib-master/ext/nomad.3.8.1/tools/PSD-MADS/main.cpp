/*-------------------------------------------------------------*/
/*                            PSD-MADS                         */
/*-------------------------------------------------------------*/
/*                                                             */
/*  usage:                                                     */
/*                                                             */
/*    "mpirun -np p psdmads param_file bbe ns"                 */
/*    with p > 2, bbe > 0, and 1 <= ns <= number of variables  */
/*    . ns is the number of free variables for each process    */
/*    . bbe is the max number of evaluations for each process  */
/*                                                             */
/*-------------------------------------------------------------*/
/*                                                             */
/*  processes:                                                 */
/*                                                             */
/*            0: master                                        */
/*            1: pollster slave (1 direction)                  */
/*     2 to p-2: regular slaves (2ns directions)               */
/*          p-1: cache server                                  */
/*                                                             */
/*-------------------------------------------------------------*/
/*  See the user guide for other details and the description   */
/*  of the algorithm                                           */
/*-------------------------------------------------------------*/

/*-----------------------------------------------------------*/
#include "Master_Slaves.hpp"
using namespace std;
using namespace NOMAD;

const bool DEBUG = false;

/*-----------------------------------*/
/*           main function           */
/*-----------------------------------*/
int main ( int argc , char ** argv )
{
    
    // MPI initialization:
    MPI_Init ( &argc , &argv );
    int rank , np;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &np   );
    
    // check the arguments and the number of processes:
    if ( np <= 2 || argc < 4 || argc > 5 )
    {
        if ( rank == 0 )
            cerr << "usage: mpirun -np p " << argv[0]
            << " param_file bbe ns [history_file(txt)], with p>2,"
            << " bbe>1, and 1<=ns<=n."
            << endl;
        MPI_Finalize();
        return 1;
    }
    
    // display:
    Display out ( cout );
    out.precision ( 16 );
    
    // parameters:
    NOMAD::Parameters p ( out );
    int bbe = atoi ( argv[2] );
    int ns  = atoi ( argv[3] );
    
    // The history file is handled by the cache server only (rank==np-1)
    string history_file;
    if ( argc == 5 && rank ==np-1 )
    {
        
        history_file = argv[5];
        
        std::ofstream fout;
        fout.open ( history_file.c_str() );
        
        if ( !fout.fail() )
        {
            fout << std::endl;
            fout.close();
        }
        else
        {
            cout << "Warning ( main.cpp, line " << __LINE__
            << "): could not open the history"
            << " in \'"
            << history_file << "\'" << std::endl << std::endl;
            
            history_file.clear();
        }
        
    }
    
    string init_cache_file,save_cache_file;
    bool success_load_init_cache = false;
    Cache init_cache( out );
    
    try
    {
        
        // read the parameters file:
        p.read ( argv[1] );
        
        // Force anisotropy(false) to have a single mesh_index (XMesh)
        p.set_ANISOTROPIC_MESH ( false );
        
        // check the parameters:
        p.check();
        
        // manage the cache file for initialization
        init_cache_file = p.get_x0_cache_file();
        
        if ( ! init_cache_file.empty() )
            throw Exception ( __FILE__ , __LINE__ ,
                             "PSDMADS cannot accept a X0 cache file. Set a regular cache file to be used for initialization instead." );
        
        
        init_cache_file = p.get_cache_file();
        
        // The init cache file is handled by the pollster (rank==1)
        if ( ! init_cache_file.empty() )
        {
            if ( rank == 1 )
            {
                success_load_init_cache = init_cache.load ( init_cache_file , NULL , DEBUG );
                
                cout << "Warning ( main.cpp, line " << __LINE__ << "): " << init_cache_file << " is used for initialization in addition of X0." << std::endl << std::endl;
                
                if ( ! success_load_init_cache || init_cache.size()==0 )
                {
                    init_cache_file.clear();
                }
            }
            
            // The cache file is not loaded by others than pollster
            p.set_CACHE_FILE("");
            p.check();
        }
        
        
        if ( ns < 1 || ns > p.get_dimension() )
            throw Exception ( __FILE__ , __LINE__ ,
                             "Bad value for ns the number of free variables for each process" );
        
        if ( p.get_nb_obj() > 1 )
            throw Exception ( __FILE__ , __LINE__ ,
                             "PSD-MADS is not designed for multi-objective optimization" );
    }
    catch ( exception & e )
    {
        if ( rank == 0 )
            cerr << "error with parameter: " << e.what()<< endl;
        MPI_Finalize();
        return 1;
    }
    
    
    // start the master:
    Master_Slaves master_slaves ( rank , np , bbe , ns , p , DEBUG );
    
    // cache server:
    Cache_Server cache ( out                  ,
                        rank                 ,
                        np                   ,
                        p.get_h_min()        ,
                        p.get_max_bb_eval()  ,
                        history_file         ,
                        false                ,  // ALLOW_MULTIPLE_EVALS
                        DEBUG                  );
    
    
    master_slaves.start( );
    
    
    // Cache server:
    if ( rank == np-1 )
    {
        if ( !DEBUG )
            out << endl << "TIME\tBBE\tOBJ" << endl << endl;
        
        // Start the cache server
        cache.start();
        
        // Save the cache when done (after start is completed)
        if ( ! init_cache_file.empty() )
        {
            bool success_load_cache_file = cache.load ( init_cache_file , NULL , DEBUG );
            
            if ( ! success_load_cache_file )
            {
                cout << "Warning ( main.cpp, line " << __LINE__
                << "): cannot load the cache file " << init_cache_file << " for saving." << std::endl << std::endl;
                
            }
            else
                cache.save( true , DEBUG );
        }
    }
    
    // slaves: algorithm creation and execution:
    if ( rank != 0 && rank != np-1 )
    {
        // MADS run:
        master_slaves.mads_run ( cache , init_cache );
        
        // stop the master:
        master_slaves.stop();
    }
    
    
    
    // stop the cache server:
    cache.stop();
    
    // display the final solution:
    if ( !DEBUG && rank == np-1 )
    {
        cache.display_best_points ( out );
    }
    
    // MPI finalization:
    MPI_Finalize();
    return 0;
}
