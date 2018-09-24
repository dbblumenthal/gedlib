/*-------------------------------------------------------------------------------------*/
/*           example of a user search for the periodic function f(x)=sin(2x)           */
/*-------------------------------------------------------------------------------------*/
#include "nomad.hpp"
using namespace std;
using namespace NOMAD;

const double PI = 3.141592654;

/*------------------------------------------------*/
/*               The problem                      */
/*------------------------------------------------*/
/*       n=1, m=1                                 */
/*       the periodic function f(x)=sin(2x)       */
/*------------------------------------------------*/
class My_Evaluator : public Evaluator {

public:

  // ctor:
  My_Evaluator  ( const Parameters & p ) :
    Evaluator ( p ) {}

  // dtor:
  ~My_Evaluator ( void ) {}

  // evaluation of a point:
  bool eval_x ( Eval_Point          & x          ,
		const NOMAD::Double & h_max      ,
		bool                & count_eval   ) const {
    x.set_bb_output ( 0 , sin ( 2*x[0].value() ) );
    count_eval = true;
    return true;
  }
};

/*------------------------------------------------*/
/*                   user search                  */
/*------------------------------------------------*/
class My_Search : public Search {

public:

  // ctor:
  My_Search ( Parameters & p )
    : Search ( p , USER_SEARCH ) {}

  // dtor:
  ~My_Search ( void ) {}

  // the search:
  void search (  Mads              & mads           ,
		 int               & nb_search_pts  ,
		 bool              & stop           ,
		 stop_type         & stop_reason    ,
		 success_type      & success        ,
		 bool              & count_search   ,
		 const Eval_Point *& new_feas_inc   ,
		 const Eval_Point *& new_infeas_inc   );
};

/*-------------------*/
/*  the user search  */
/*-------------------*/
void My_Search::search ( Mads              & mads           ,
			 int               & nb_search_pts  ,
			 bool              & stop           ,
			 stop_type         & stop_reason    ,
			 success_type      & success        ,
			bool              & count_search   ,
			 const Eval_Point *& new_feas_inc   ,
			 const Eval_Point *& new_infeas_inc   ) {

  nb_search_pts = 0;
  success       = UNSUCCESSFUL;
  count_search  = false;

  // current feasible incumbent:
  const Eval_Point * feas_inc = mads.get_best_feasible();

  if ( !feas_inc )
    return;

  // xk:
  double xk = (*feas_inc)[0].value();
  if ( xk < 0 )
    return;

  // get a signature:
  Signature * signature = feas_inc->get_signature();
  if ( !signature )
    return;

  // count the search:
  count_search  = true;

  // construct the search point (tk = qk.xk - 1):
  Eval_Point * tk = new Eval_Point;
  tk->set ( 1 , 1 );
  tk->set_signature  ( signature );

  (*tk)[0] = static_cast<int>(ceil(1.0/xk)) * xk - 1.0;
    
    
    // Projection maybe needed
    //    const NOMAD::Display & out= _p.out();
    //    NOMAD::dd_type display_degree = out.get_search_dd();
    //  if ( display_degree == NOMAD::FULL_DISPLAY )
    //  {
    //		out << "candidate";
    //        out << " (before projection)";
    //		out << ": ( " << *tk << " )" << std::endl;
    //  }
    //
    //  // Project to the mesh
    //  tk->project_to_mesh(*feas_inc,signature->get_mesh()->get_delta(),signature->get_lb(),signature->get_ub() );
    //
    //    if ( display_degree == NOMAD::FULL_DISPLAY )
    //    {
    //		out << "candidate";
    //        out << " (after projection)";
    //		out << ": ( " << *tk << " )" << std::endl;
    //    }
    
    

  // Evaluator_Control:
  Evaluator_Control & ev_control = mads.get_evaluator_control();

  // add the new point to the ordered list of search trial points:
  ev_control.add_eval_point ( tk                       ,
			      _p.out().get_search_dd() ,
			      false                    ,
			      Double()                 ,
			      Double()                  ,
				  Double()                 ,
				  Double()                 );

  nb_search_pts = 1;

  // evaluation:
  new_feas_inc = new_infeas_inc = NULL;
  ev_control.eval_list_of_points ( _type                   ,
				   mads.get_true_barrier() ,
				   mads.get_sgte_barrier() ,
				   mads.get_pareto_front() ,
				   stop                    ,
				   stop_reason             ,
				   new_feas_inc            ,
				   new_infeas_inc          ,
				   success                   );
}

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv ) {

  // NOMAD initializations:
  begin ( argc , argv );

  // display:
  Display out ( std::cout );
  out.precision ( DISPLAY_PRECISION_STD );

  //  parameters creation:
  Parameters p ( out );

  p.set_DIMENSION (1);             // number of variables

  vector<bb_output_type> bbot (1); // definition of
  bbot[0] = OBJ;                   // output types
  p.set_BB_OUTPUT_TYPE ( bbot );

  p.set_X0 ( Point ( 1 , PI-3.0 ) );  // starting point

  p.set_LOWER_BOUND ( Point ( 1 , -PI/2.0 ) );
  p.set_UPPER_BOUND ( Point ( 1 ,  PI/2.0 ) );

  p.set_SPECULATIVE_SEARCH ( false );

  p.set_INITIAL_MESH_SIZE        ( 1.0 );

  p.set_DISPLAY_DEGREE ( "0300" ); // display only the search step

  p.set_DISPLAY_STATS ( "bbe sol obj" );

  p.set_MAX_BB_EVAL ( 10 );

  // parameters validation:
  p.check();

  // custom evaluator creation:
  My_Evaluator ev ( p );

  // algorithm creation:
  Mads mads ( p , &ev );

  // user search:
  My_Search my_search  ( p );
  mads.set_user_search ( &my_search );

  // algorithm execution:
  mads.run();

  Slave::stop_slaves ( out );
  end();

  return EXIT_SUCCESS;
}
