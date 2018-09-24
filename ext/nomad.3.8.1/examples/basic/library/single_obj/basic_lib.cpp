/*-----------------------------------------------------*/
/*  how to use the NOMAD library with a user function  */
/*-----------------------------------------------------*/
#include "nomad.hpp"
using namespace std;
// using namespace NOMAD; avoids putting NOMAD:: everywhere

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
class My_Evaluator : public NOMAD::Evaluator {
public:
  My_Evaluator  ( const NOMAD::Parameters & p ) :
    NOMAD::Evaluator ( p ) {}

  ~My_Evaluator ( void ) {}

  bool eval_x ( NOMAD::Eval_Point   & x          ,
		const NOMAD::Double & h_max      ,
		bool                & count_eval   ) const
	{

    NOMAD::Double c1 = 0.0 , c2 = 0.0;
    for ( int i = 0 ; i < 5 ; i++ ) 
	{
      c1 += (x[i]-1).pow2();
      c2 += (x[i]+1).pow2();
    }
    x.set_bb_output  ( 0 , x[4]  ); // objective value
    x.set_bb_output  ( 1 , c1-25 ); // constraint 1
    x.set_bb_output  ( 2 , 25-c2 ); // constraint 2

    count_eval = true; // count a black-box evaluation

    return true;       // the evaluation succeeded
}
	
	
};

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv ) {

  // display:
  NOMAD::Display out ( std::cout );
  out.precision ( NOMAD::DISPLAY_PRECISION_STD );

  try {

    // NOMAD initializations:
    NOMAD::begin ( argc , argv );

    // parameters creation:
    NOMAD::Parameters p ( out );

    p.set_DIMENSION (5);             // number of variables

    vector<NOMAD::bb_output_type> bbot (3); // definition of
    bbot[0] = NOMAD::OBJ;                   // output types
    bbot[1] = NOMAD::PB;
    bbot[2] = NOMAD::EB;
    p.set_BB_OUTPUT_TYPE ( bbot );

//    p.set_DISPLAY_ALL_EVAL(true);   // displays all evaluations.
    p.set_DISPLAY_STATS ( "bbe ( sol ) obj" );

    p.set_X0 ( NOMAD::Point(5,0.0) );  // starting point

    p.set_LOWER_BOUND ( NOMAD::Point ( 5 , -6.0 ) ); // all var. >= -6
    NOMAD::Point ub ( 5 );                    // x_4 and x_5 have no bounds
    ub[0] = 5.0;                              // x_1 <= 5
    ub[1] = 6.0;                              // x_2 <= 6
    ub[2] = 7.0;                              // x_3 <= 7
    p.set_UPPER_BOUND ( ub );

    p.set_MAX_BB_EVAL (100);     // the algorithm terminates after
                                 // 100 black-box evaluations
    p.set_DISPLAY_DEGREE(2);
    p.set_SOLUTION_FILE("sol.txt");

    // parameters validation:
    p.check();

    // custom evaluator creation:
    My_Evaluator ev   ( p );

    // algorithm creation and execution:
    NOMAD::Mads mads ( p , &ev );
    mads.run();
  }
  catch ( exception & e ) {
    cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
  }

  NOMAD::Slave::stop_slaves ( out );
  NOMAD::end();

  return EXIT_SUCCESS;
}
