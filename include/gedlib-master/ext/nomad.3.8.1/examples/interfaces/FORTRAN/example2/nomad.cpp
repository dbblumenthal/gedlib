/*---------------------------------*/
/*  FORTRAN and the NOMAD library  */
/*        (see readme.txt)         */
/*---------------------------------*/
#include "nomad.hpp"
using namespace std;
using namespace NOMAD;

extern "C" {

  void nomad_ ( int    * n              ,
		int    * m              ,
		double * x              ,
		double * lb             ,
		double * ub             ,
		int    * max_bbe        ,
		int    * display_degree   );

  void bb_ ( double x[5] , double fx[3] );
}

/*----------------------------------------*/
/*  The problem: the evaluation is made   */
/*  by calling a FORTRAN routine          */
/*----------------------------------------*/
class My_Evaluator : public Evaluator {
public:
  My_Evaluator  ( const Parameters & p ) :
    Evaluator ( p ) {}

  ~My_Evaluator ( void ) {}

  bool eval_x ( Eval_Point   & x          ,
                const Double & h_max      ,
                bool         & count_eval   ) const {

    int n = x.size();
    int m = x.get_bb_outputs().size();

    int i;

    double * xx = new double [n];
    double * fx = new double [m];

    for ( i = 0 ; i < x.size() ; ++i )
      xx[i] = x[i].value();

    // call the FORTRAN routine:
    bb_ ( xx , fx );

    for ( i = 0 ; i < m ; ++i )
      x.set_bb_output ( i , fx[i] );

    count_eval = true; // count a black-box evaluation

    return true;       // the evaluation succeeded
  }
};

/*---------------------------------------*/
/*  this routine launches NOMAD and can  */
/*  be called from a FORTRAN program     */
/*---------------------------------------*/
void nomad_ ( int    * n       ,          // # of variables
	      int    * m       ,          // # of outputs (obj + m-1 constraints)
	      double * x       ,          // starting point (IN) / solution (OUT)
	      double * lb      ,          // lower bounds for each variable
	      double * ub      ,          // upper bounds
	      int    * max_bbe ,          // max # of evaluations (-1: not considered)
	      int    * display_degree ) { // display_degree (0-4; 0: no display)

  // display:
  Display out ( std::cout );
  out.precision ( DISPLAY_PRECISION_STD );

  try {

    int i;

    // parameters creation:
    Parameters p ( out );

    p.set_DIMENSION (*n);             // number of variables

    vector<bb_output_type> bbot (*m); // definition of output types:
    bbot[0] = OBJ;                    // first output : objective value to minimize
    for ( i = 1 ; i < *m ; ++i )      // other outputs: constraints cj <= 0
      bbot[i] = PB;
    p.set_BB_OUTPUT_TYPE ( bbot );

    // starting point and bounds:
    Point px0 ( *n );
    Point plb ( *n );
    Point pub ( *n );
    for ( i = 0 ; i < *n ; ++i ) {
      px0[i] = x [i];
      if ( lb[i] > -1e20 )
	plb[i] = lb[i];
      if ( ub[i] < 1e20 )
	pub[i] = ub[i];
    }
    p.set_X0          ( px0 );
    p.set_LOWER_BOUND ( plb );
    p.set_UPPER_BOUND ( pub );

    // maximum number of black-box evaluations:
    if ( *max_bbe > 0 )
      p.set_MAX_BB_EVAL ( *max_bbe );

    // display degree:
    p.set_DISPLAY_DEGREE ( *display_degree );

    // parameters validation:
    p.check();

    // custom evaluator creation:
    My_Evaluator ev ( p );

    // algorithm creation and execution:
    Mads mads ( p , &ev );
    mads.run();

    // get the solution:
    const Eval_Point * bf = mads.get_best_feasible();
    if ( bf )
      for ( i = 0 ; i < *n ; ++i )
	x[i] = (*bf)[i].value();
  }
  catch ( exception & e ) {
    cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
  }
}
