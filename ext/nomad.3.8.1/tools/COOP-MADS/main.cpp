/*-----------------------------------------------------*/
/*                      COOP-MADS                      */
/*  see user guide for a description of the algorithm  */
/*-----------------------------------------------------*/
#include "Cache_Server.hpp"
using namespace std;
using namespace NOMAD;

// COOP-MADS parameters:
const bool USE_CACHE_SEARCH     = true;
const bool OPPORT_CACHE_SEARCH  = false;
const bool ALLOW_MULTIPLE_EVALS = false;

/*-----------------------------------*/
/*           main function           */
/*-----------------------------------*/
int main ( int argc , char ** argv ) {

  // MPI initialization:
  MPI_Init ( &argc , &argv );
  int rank , np;
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &np   );

  // check the arguments and the number of processes:
  if ( np <= 1 || argc != 2 ) {
    if ( rank==np-1 )
      cerr << "usage: mpirun -np p " << argv[0]
	   << " param_file, with p>1" << endl;
    MPI_Finalize();
    return 1;
  }

  // display:
  Display out ( cout );
  out.precision ( 16 );

  // parameters:
  Parameters p ( out );

  try {

    // read the parameters file:
    p.read ( argv[1] );
       
    // modify parameters:
    p.set_DISPLAY_DEGREE ( NO_DISPLAY );
    p.set_DISPLAY_STATS ( "process #" + itos(rank) + " BBE OBJ" );
    p.set_SEED ( get_pid() );
      
    p.set_ANISOTROPIC_MESH(false);

    // cache search:
    p.set_CACHE_SEARCH               ( USE_CACHE_SEARCH    );
    p.set_OPPORTUNISTIC_CACHE_SEARCH ( OPPORT_CACHE_SEARCH );  

    // check the parameters:
    p.check();

    if ( p.get_nb_obj() > 1 )
      throw Exception ( __FILE__ , __LINE__ ,
      "COOP-MADS is not designed for multi-objective optimization" );
  }
  catch ( exception & e ) {
    if ( rank==np-1 )
      cerr << "error with parameters" << endl;
    MPI_Finalize();
    return 1;
  }

  // cache server:
  Cache_Server cache ( out                  ,
		       rank                 ,
		       np                   ,
		       p.get_h_min()        ,
		       p.get_max_bb_eval()  ,
		       ALLOW_MULTIPLE_EVALS   );

  if ( rank == np-1 )
    out << endl << "TIME\tBBE\tOBJ" << endl << endl;

  // start the cache server:
  cache.start();

  // clients: algorithm creation and execution:
  if ( rank != np-1 ) {
    Mads mads ( p , NULL , NULL , &cache , NULL );
    mads.run();
  }

  // stop the cache server:
  cache.stop();

  // display the best solutions:
  cache.display_current_solution();
  if ( rank == np-1 )
    out << endl;
  cache.display_best_points();

  // MPI finalization:
  MPI_Finalize();
  return 0;
}
