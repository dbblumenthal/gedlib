/*-------------------------------------------------------------*/
/*  how to use the NOMAD library with a FORTRAN user function  */
/*-------------------------------------------------------------*/
#include "nomad.hpp"
using namespace std;
using namespace NOMAD;

extern"C" {
  void bb_ ( double x[5] , double fx[3] );
}

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
class My_Evaluator : public Evaluator {
public:
  My_Evaluator  ( const Parameters & p ) :
    Evaluator ( p ) {}

  ~My_Evaluator ( void ) {}

  bool eval_x ( Eval_Point   & x          ,
		const Double & h_max      ,
		bool         & count_eval   ) const {

    double xx[5];
    double fx[3];
    int    i;

    for ( i = 0 ; i < 5 ; ++i )
      xx[i] = x[i].value();

    // call the FORTRAN routine:
    bb_ ( xx , fx );

    for ( i = 0 ; i < 3 ; ++i )
      x.set_bb_output ( i , fx[i] );

    count_eval = true; // count a black-box evaluation

    return true;       // the evaluation succeeded
  }
};

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv ) {

  // display:
  Display out ( std::cout );
  out.precision ( DISPLAY_PRECISION_STD );

  // NOMAD initializations:
  begin ( argc , argv );

  try {

    // parameters creation:
    Parameters p ( out );

    p.set_DIMENSION (5);             // number of variables

    vector<bb_output_type> bbot (3); // definition of
    bbot[0] = OBJ;                 // output types
    bbot[1] = PB;
    bbot[2] = EB;
    p.set_BB_OUTPUT_TYPE ( bbot );

    p.set_X0 ( Point ( 5 , 0.0 ) );  // starting point

    p.set_LOWER_BOUND ( Point ( 5 , -6.0 ) ); // all var. >= -6
    Point ub ( 5 );                           // x_4 and x_5 have no bounds
    ub[0] = 5.0;                              // x_1 <= 5
    ub[1] = 6.0;                              // x_2 <= 6
    ub[2] = 7.0;                              // x_3 <= 7
    p.set_UPPER_BOUND ( ub );

    p.set_MAX_BB_EVAL (100);     // the algorithm terminates after
                                 // 100 black-box evaluations

    // p.set_TMP_DIR ("/tmp");      // repertory for
                                   // temporary files

    // parameters validation:
    p.check();

    // custom evaluator creation:
    My_Evaluator ev   ( p );

    // algorithm creation and execution:
    Mads mads ( p , &ev );
    mads.run();
  }
  catch ( exception & e ) {
    cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
  }

  Slave::stop_slaves ( out );
  end();

  return EXIT_SUCCESS;
}
