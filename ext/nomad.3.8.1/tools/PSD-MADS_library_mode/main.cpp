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



/*---------------------------------------------------*/
/*      The evaluator for G2_50 problem              */
/*---------------------------------------------------*/
class My_Evaluator : public NOMAD::Evaluator {
public:
    My_Evaluator  ( const NOMAD::Parameters & p ) :
    NOMAD::Evaluator ( p ) {}
    
    ~My_Evaluator ( void ) {}
    
    bool eval_x ( NOMAD::Eval_Point   & x          ,
                 const NOMAD::Double & h_max      ,
                 bool                & count_eval   ) const
	{
        
        int N=50;
        long double sum1 = 0.0 , sum2 = 0.0 , sum3 = 0.0 , prod1 = 1.0 , prod2 = 1.0, g1=0, g2=0;
        long double Xi;
        
        for (int i = 0 ; i < N ; i++ )
        {
            
            Xi=x[i].value();
            
            sum1  += pow ( cos(Xi) , 4 );
            sum2  += Xi;
            sum3  += (i+1)*Xi*Xi;
            prod1 *= pow ( cos(Xi) , 2 );
            prod2 *= Xi;
        }
        
        
        g1 = -prod2+0.75;
        g2 = sum2 -7.5*N;
        
        long double z  = - fabs ( ( sum1 - 2 * prod1 ) / sqrt(sum3) );
        
        x.set_bb_output  ( 0 , g1 ); // constraint 1
        x.set_bb_output  ( 1 , g2 ); // constraint 2
        x.set_bb_output  ( 2 , z  ); // objective value
        
        
        count_eval = true; // count a black-box evaluation
        
        
        return true;       // the evaluation succeeded
	}
	
	
};



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
    if ( np <= 2 || argc != 4 ) {
        if ( rank == 0 )
            cerr << "usage: mpirun -np p " << argv[0]
            << " param_file bbe ns, with p>2,"
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
    
    try {
        
        // read the parameters file:
        p.read ( argv[1] );
        
        // This option is needed to have the same mesh index for all variables
        p.set_ANISOTROPIC_MESH ( false );
                
        // check the parameters:
        p.check();

        
        if ( ns < 1 || ns > p.get_dimension() )
            throw Exception ( __FILE__ , __LINE__ ,
                             "Bad value for ns the number of free variables for each process" );
        
        if ( p.get_nb_obj() > 1 )
            throw Exception ( __FILE__ , __LINE__ ,
                             "PSD-MADS is not designed for multi-objective optimization" );
    }
    catch ( exception & e ) {
        if ( rank == 0 )
            cerr << "error with parameters" << endl;
        MPI_Finalize();
        return 1;
    }
    
    // custom evaluator creation:
    My_Evaluator ev   ( p );
    
    // start the master:
    Master_Slaves master_slaves ( rank , np , bbe , ns , p , DEBUG ,ev);
    master_slaves.start();
    
    // cache server:
    Cache_Server cache ( out                  ,
                        rank                 ,
                        np                   ,
                        p.get_h_min()        ,
                        p.get_max_bb_eval()  ,
                        false                ,  // ALLOW_MULTIPLE_EVALS
                        DEBUG                  );
    
    // start the cache server:
    if ( rank == np-1 ) {
        if ( !DEBUG )
            out << endl << "TIME\tBBE\tOBJ" << endl << endl;
        cache.start();
    }
    
    // slaves: algorithm creation and execution:
    if ( rank != 0 && rank != np-1 ) {
        
        // MADS run:
        master_slaves.mads_run ( cache );
        
        // stop the master:
        master_slaves.stop();
    }
    
    // stop the cache server:
    cache.stop();
    
    // display the final solution:
    if ( !DEBUG && rank == np-1 )
        cache.display_best_points ( out );
    
    // MPI finalization:
    MPI_Finalize();
    return 0;
}
