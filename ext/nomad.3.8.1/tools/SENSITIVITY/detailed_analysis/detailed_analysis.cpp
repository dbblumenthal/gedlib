/*----------------------------------------------------------------------------------------*/
/*                    Sensitivity analysis with BiMADS (detailed method)                  */
/*----------------------------------------------------------------------------------------*/
/*                                                                                        */
/*  This program performs a sensitivity analysis by running a biobjective optimization    */
/*                                                                                        */
/*  The two objective functions are:                                                      */
/*             1. the original objective                                                  */
/*             2. the constraint of interest that may be cj<=0 or a bound lb<=x or x<=ub  */
/*                                                                                        */
/*  The output is a cache file which may be opened by the 'cache_inspect' program         */
/*                                                                                        */
/*  usage:                                                                                */
/*                                                                                        */
/*  ./detailed_analysis.exe param_file cstr_index cache_file lj uj bb_eval init_bb_eval   */
/*                                                                                        */
/*         param_file  : parameters file                                                  */
/*         cache_file  : cache file to be treated with cache_inspect                      */
/*         cstr_index  : constraint index j in [0;m+2n-1]                                 */
/*                       m corresponds to the number of outputs                           */
/*                       n corresponds to the number of variables                         */
/*                       j<m corresponds to constraint cj<=0                              */
/*                       j=m+2i corresponds to the lower bound of variable i in [0;n-1]   */
/*                       j=m+2i+1 corresponds to the upper bound of variable i in [0;n-1] */
/*         lj, uj      : optional, may be - or INF                                        */
/*                       constraint cj<=0 is replaced by cj<=e with lj<=e<=uj             */
/*                       lower bound lb-xi<=0 is replaced by lb-xi<=e with lj<=e<=uj      */
/*                       upper bound xi-ub<=0 is replaced by xi-ub<=e with lj<=e<=uj      */
/*         bb_eval     : optional, max number of blackbox evaluations                     */
/*                       if none specified, param_file is considered                      */
/*         init_bb_eval: optional, max number of blackbox evaluations for the two         */
/*                       first MADS runs (min f1 and min f2)                              */
/*                       if none specified, a default of 10% of the max number of         */
/*                       evaluations is considered                                        */
/*                                                                                        */
/*----------------------------------------------------------------------------------------*/
#include "nomad.hpp"
#include <cstring>

using std::endl;

/*-------------------------------------*/
/*               prototypes            */
/*-------------------------------------*/
bool check_args_1 ( int                    ,
		    char                ** ,
		    const NOMAD::Display & ,
		    std::string          & ,
		    int                  & ,
		    std::string          & ,
		    NOMAD::Double        & ,
		    NOMAD::Double        & ,
		    int                  & ,
		    int                  &   );

void check_args_2 ( const NOMAD::Parameters & ,
		    int                       ,
		    const NOMAD::Double     & ,
		    const NOMAD::Double     & ,
		    std::string             &   );

bool check_bounds ( const NOMAD::Point & ,
		    const NOMAD::Point & ,
		    const NOMAD::Point &   );

void display_usage ( const NOMAD::Display & , char ** );

void display_objectives ( const NOMAD::Display    & ,
			  const NOMAD::Parameters & ,
			  const NOMAD::Parameters & ,
			  int                       ,
			  const NOMAD::Double     & ,
			  const NOMAD::Double     &   );

void create_new_parameters ( const NOMAD::Parameters & ,
			     NOMAD::Parameters       & ,
			     const NOMAD::Cache      & ,
			     int                       ,
			     const NOMAD::Double     & ,
			     const NOMAD::Double     & ,
			     int                       ,
			     int                       ,
			     int                     &   );

const NOMAD::Eval_Point * select_x0 ( const NOMAD::Parameters & ,
				      const NOMAD::Cache      & ,
				      const NOMAD::Point      & ,
				      const NOMAD::Point      & ,
				      const NOMAD::Double     & ,
				      const NOMAD::Double     & ,
				      int                       ,
				      bool                        );

/*----------------------------------------*/
/*            multi-obj evaluator         */
/*----------------------------------------*/
class Sensitivity_Evaluator : public NOMAD::Multi_Obj_Evaluator {
private:
  NOMAD::Evaluator    _ev0;
  NOMAD::Cache      * _cache;
  int                 _cstr_index;
  int                 _m;
  const NOMAD::Double _lj;
  const NOMAD::Double _uj;
  int                 _auto_max_bbe;
public:
  Sensitivity_Evaluator ( const NOMAD::Parameters & p            ,
			  const NOMAD::Parameters & p0           ,
			  NOMAD::Cache            & cache        ,
			  int                       cstr_index   ,
			  const NOMAD::Double     & lj           ,
			  const NOMAD::Double     & uj           ,
			  int                       auto_max_bbe   ) :
    Multi_Obj_Evaluator ( p                      ) ,
    _ev0                ( p0                     ) ,
    _cache              ( &cache                 ) ,
    _cstr_index         ( cstr_index             ) ,
    _m                  ( p0.get_bb_nb_outputs() ) ,
    _lj                 ( lj                     ) ,
    _uj                 ( uj                     ) ,
    _auto_max_bbe       ( auto_max_bbe           )   {}

  ~Sensitivity_Evaluator ( void ) {}

  virtual bool eval_x ( NOMAD::Eval_Point   & x          ,
			const NOMAD::Double & h_max      ,
			bool                & count_eval   ) const;
  
  virtual void update_mads_run ( const NOMAD::Stats             & stats        ,
				 const NOMAD::Evaluator_Control & ev_control   ,
				 const NOMAD::Barrier           & true_barrier ,
				 const NOMAD::Barrier           & sgte_barrier ,
				 const NOMAD::Pareto_Front      & pareto_front   );
};

/*-------------------------------------*/
/*                eval_x()             */
/*-------------------------------------*/
bool Sensitivity_Evaluator::eval_x
( NOMAD::Eval_Point   & x          ,
  const NOMAD::Double & h_max      ,
  bool                & count_eval   ) const {

  bool eval_ok;
  int  i , n = x.size();

  NOMAD::Eval_Point * y = new NOMAD::Eval_Point ( n , _m );
  for ( i = 0 ; i < n ; ++i )
    (*y)[i] = x[i];

  // search in cache:
  const NOMAD::Eval_Point * cache_x = _cache->find ( *y );

  // point in cache:
  if ( cache_x ) {
    count_eval = false;
    delete y;
  }

  // point not in cache:
  else {
    eval_ok = _ev0.eval_x ( *y , h_max , count_eval );
    _cache->insert ( *y );
    cache_x = y;
  }
 
  // set bb outputs:
  for ( i = 0 ; i < _m ; ++i )
    x.set_bb_output ( i , (cache_x->get_bb_outputs())[i] );

  // constraint index j indicates a bound constraint:
  if ( _cstr_index >= _m ) {

    // index i of variable xi:
    i = ( _cstr_index - _m ) / 2;

    // lower bound:
    if ( ( _cstr_index - _m ) % 2 == 0 )
      x.set_bb_output ( _m , -x[i] );

    // upper bound:
    else
      x.set_bb_output ( _m ,  x[i] );
  }

  // regular constraint cj <= 0 replaced by lj <= cj <= uj:
  else {
    const NOMAD::Double & cj = (cache_x->get_bb_outputs())[_cstr_index];
    if ( ( _uj.is_defined() &&  cj > _uj ) ||
	 ( _lj.is_defined() &&  cj < _lj )    ) {
      x.set_bb_output   ( _cstr_index , NOMAD::INF );
      x.set_eval_status ( NOMAD::EVAL_FAIL         );
      return false;
    }
  }

  x.set_eval_status ( cache_x->get_eval_status() );
  eval_ok = cache_x->is_eval_ok();

  return eval_ok;
}

/*-------------------------------------*/
/*     code called between MADS runs   */
/*-------------------------------------*/
void Sensitivity_Evaluator::update_mads_run
( const NOMAD::Stats             & stats         ,
  const NOMAD::Evaluator_Control & ev_control    ,
  const NOMAD::Barrier           & true_barrier  ,
  const NOMAD::Barrier           & sgte_barrier  ,
  const NOMAD::Pareto_Front      & pareto_front    ) {

  if ( _auto_max_bbe > 0 && stats.get_mads_runs() == 2 ) {
    NOMAD::Parameters * p = const_cast<NOMAD::Parameters *> ( &_p );
    p->set_MAX_BB_EVAL ( _auto_max_bbe );
    p->check();
  }
}

/*-------------------------------------*/
/*-------------------------------------*/
/*            main function            */
/*-------------------------------------*/
/*-------------------------------------*/
int main ( int argc , char ** argv ) {

  NOMAD::begin ( argc , argv );

  // display:
  NOMAD::Display out ( std::cout );
  out.precision ( NOMAD::DISPLAY_PRECISION_STD );

  // check and interpret the arguments:
  std::string   err;
  std::string   param_file;
  std::string   cache_file;
  int           cstr_index;
  int           command_line_bbe;
  int           command_line_bbe0;
  NOMAD::Double lj , uj;

  if ( !check_args_1 ( argc              ,
		       argv              ,
		       out               ,
		       param_file        ,
		       cstr_index        ,
		       cache_file        ,
		       lj                ,
		       uj                ,
		       command_line_bbe  ,
		       command_line_bbe0   ) ) {
    NOMAD::end();
    return EXIT_FAILURE;
  }

  // read the cache file:
  NOMAD::Cache cache ( out );
  if ( !cache.load ( cache_file , NULL , false ) )
    return EXIT_FAILURE;
  
  try {

    // read the original parameters:
    NOMAD::Parameters p0 ( out );
    p0.read ( param_file );
    p0.check();

    // check the arguments:
    check_args_2 ( p0 , cstr_index , lj , uj , err );
    if ( !err.empty() )
      throw NOMAD::Exception ( "detailed_analysis.cpp" , __LINE__ , err );

    // create new parameters:
    int auto_max_bbe = -1;
    NOMAD::Parameters p  ( out );
    create_new_parameters ( p0                ,
			    p                 ,
			    cache             ,
			    cstr_index        ,
			    lj                ,
			    uj                ,
			    command_line_bbe  ,
			    command_line_bbe0 ,
			    auto_max_bbe        );
    p.check();
    // out << p << endl;

    // display the objectives:
    display_objectives ( out , p0 , p , cstr_index , lj , uj );

    // sensitivity evaluator:
    Sensitivity_Evaluator ev ( p , p0 , cache , cstr_index , lj , uj , auto_max_bbe );

    // biobjective optimization:
    NOMAD::Mads mads ( p , &ev );
    mads.multi_run();

    // save the cache file:
    cache.save ( true , true );
    out << endl;
  }
  catch ( std::exception & e ) {
    if ( NOMAD::Slave::is_master() ) {
      err = std::string ( "NOMAD has been interrupted: " ) + e.what();
      std::cerr << endl << err << endl << endl;
    }
  }

  NOMAD::Slave::stop_slaves ( out );
  NOMAD::end();

  return EXIT_SUCCESS;
}

/*-------------------------------------*/
/*        create new parameters        */
/*-------------------------------------*/
void create_new_parameters ( const NOMAD::Parameters & p0                ,
			     NOMAD::Parameters       & p                 ,
			     const NOMAD::Cache      & cache             ,
			     int                       cstr_index        ,
			     const NOMAD::Double     & lj                ,
			     const NOMAD::Double     & uj                ,
			     int                       command_line_bbe  ,
			     int                       command_line_bbe0 ,
			     int                     & auto_max_bbe        ) {
  // copy parameters:
  p.set_DIMENSION                    ( p0.get_dimension()                    );
  p.set_BB_INPUT_TYPE                ( p0.get_bb_input_type()                );
  p.set_H_MIN                        ( p0.get_h_min()                        );
  p.set_H_NORM                       ( p0.get_h_norm()                       );
  p.set_ADD_SEED_TO_FILE_NAMES       ( p0.get_add_seed_to_file_names()       );
  p.set_DIRECTION_TYPE               ( p0.get_direction_types()              );
  p.set_SEC_POLL_DIR_TYPE            ( p0.get_sec_poll_dir_types()           );
  p.set_OPPORTUNISTIC_EVAL           ( p0.get_opportunistic_eval()           );
  p.set_OPPORTUNISTIC_MIN_NB_SUCCESS ( p0.get_opportunistic_min_nb_success() );
  p.set_OPPORTUNISTIC_MIN_EVAL       ( p0.get_opportunistic_min_eval()       );
  p.set_OPPORTUNISTIC_MIN_F_IMPRVMT  ( p0.get_opportunistic_min_f_imprvmt()  );
  p.set_OPPORTUNISTIC_LUCKY_EVAL     ( p0.get_opportunistic_lucky_eval()     );
  p.set_SPECULATIVE_SEARCH           ( p0.get_speculative_search()           );
  p.set_VNS_SEARCH                   ( p0.get_VNS_trigger()                  );
  p.set_CACHE_SEARCH                 ( p0.get_cache_search()                 );
  p.set_OPPORTUNISTIC_CACHE_SEARCH   ( p0.get_opportunistic_cache_search()   );
  p.set_SNAP_TO_BOUNDS               ( p0.get_snap_to_bounds()               );
  p.set_SEED                         ( p0.get_seed()                         );
  p.set_HISTORY_FILE                 ( p0.get_history_file()                 );
  p.set_MULTI_USE_DELTA_CRIT         ( p0.get_multi_use_delta_crit()         );
  p.set_MULTI_FORMULATION            ( p0.get_multi_formulation()            );

  if ( p0.get_max_time() > 0 )
    p.set_MAX_TIME ( p0.get_max_time() );

  if ( p0.get_max_cache_memory() > 0 )
    p.set_MAX_CACHE_MEMORY ( p0.get_max_cache_memory() );

  if ( p0.get_multi_nb_mads_runs() > 0 )
    p.set_MULTI_NB_MADS_RUNS ( p0.get_multi_nb_mads_runs() );

  if ( p0.get_multi_f_bounds().is_defined() )
    p.set_MULTI_F_BOUNDS ( p0.get_multi_f_bounds() );

  // LH search is disabled:
  p.set_LH_SEARCH ( 0 , 0 );
  //   int n0 = p0.get_LH_search_p0();
  //   int ni = p0.get_LH_search_pi();
  //   if ( n0 >= 0 || ni >= 0 ) {
  //     p.set_LH_SEARCH ( n0 , ni );
  //     p.set_OPPORTUNISTIC_LH ( p0.get_opportunistic_LH() );
  //   }

  // DISPLAY_DEGREE:
  p.set_DISPLAY_DEGREE ( NOMAD::NORMAL_DISPLAY );

  // BB_OUTPUT_TYPE:
  std::vector<NOMAD::bb_output_type> bbo = p0.get_bb_output_type();
  int                                m   = static_cast<int>(bbo.size());
  if ( cstr_index < m )
    bbo[cstr_index] = NOMAD::OBJ;
  else
    bbo.push_back ( NOMAD::OBJ );
  p.set_BB_OUTPUT_TYPE ( bbo );
  
  // bounds and rescaling:
  NOMAD::Point lb = p0.get_lb();
  NOMAD::Point ub = p0.get_ub();
  {
    
    int m  = p0.get_bb_nb_outputs();

    if ( cstr_index >= m ) {

      NOMAD::Point r1 = ub - lb;

      int i = ( cstr_index - m ) / 2;

      // lower bound:
      if ( ( cstr_index - m ) % 2 == 0 ) {
	if ( lj.is_defined() )
	  ub[i]  = lb[i]-lj;
	if ( uj.is_defined() )
	  lb[i] -= uj;
      }
    
      // upper bound:
      else {
	if ( lj.is_defined() )
	  lb[i]  = lj+ub[i];
	if ( uj.is_defined() )
	  ub[i] += uj;
      } 

      // rescaling via initial mesh size:
      if ( r1[i].is_defined() &&
	   r1[i] != 0.0       &&
	   ub[i].is_defined() &&
	   lb[i].is_defined()    ) {
	NOMAD::Point delta_0 = p0.get_initial_mesh_size();
	delta_0 *= ((ub[i]-lb[i]) / r1[i]);
	p.set_INITIAL_MESH_SIZE ( delta_0 );
      }
    }
  
    p.set_LOWER_BOUND ( lb );
    p.set_UPPER_BOUND ( ub );
  }

  // starting point: we select a feasible starting point from the cache:
  const NOMAD::Eval_Point * x0 = select_x0 ( p0         ,
					     cache      ,
					     lb         ,
					     ub         ,
					     lj         ,
					     uj         ,
					     cstr_index ,
					     true         );
  // set the starting point:
  if ( x0 )
    p.set_X0 ( *x0 );

  // if no starting point has been chosen,
  // consider the user-defined starting points:
  else {
    bool no_x0 = true;
    std::vector<NOMAD::Point *> x0s = p0.get_x0s();
    for ( size_t k = 0 ; k < x0s.size() ; ++k ) {
      if ( check_bounds ( *x0s[k] , lb , ub ) ) {
	p.set_X0 ( *x0s[k] );
	no_x0 = false;
      }
    }
    if ( no_x0 ) {

      // if still no starting point could be selected,
      // choose an infeasible starting points from the cache:
      x0 = select_x0 ( p0         ,
		       cache      ,
		       lb         ,
		       ub         ,
		       lj         ,
		       uj         ,
		       cstr_index ,
		       false        );

      if ( x0 )
	p.set_X0 ( *x0 );
      else
	std::cerr << "Warning: could not find a valid starting point" << endl;
    }
  }

  // max number of evaluations:
  auto_max_bbe    = -1;
  int max_bbe     = p0.get_max_bb_eval();
  int overall_bbe = p0.get_multi_overall_bb_eval();

  if ( command_line_bbe < 0 ) {

    // max_bbe and overall_bbe are defined:
    if ( max_bbe > 0 && overall_bbe > 0 ) {
      p.set_MULTI_OVERALL_BB_EVAL ( overall_bbe );
      p.set_MAX_BB_EVAL           ( max_bbe     );
    }
    else {

      // max_bbe is defined and overall_bbe is undefined:
      if ( overall_bbe <= 0 && max_bbe > 0 ) {

	p.set_MULTI_OVERALL_BB_EVAL ( max_bbe );

	p.check();
	auto_max_bbe = p.get_max_bb_eval();

	max_bbe = max_bbe * 10 / 100;

	if ( max_bbe <= auto_max_bbe )
	  auto_max_bbe = -1;
	else
	  p.set_MAX_BB_EVAL ( max_bbe );
      }
    
      // max_bbe is undefined and overall_bbe is defined:
      else if ( overall_bbe > 0 && max_bbe <= 0 ) {
	
	p.set_MULTI_OVERALL_BB_EVAL ( overall_bbe );

	p.check();
	auto_max_bbe = p.get_max_bb_eval();

	max_bbe = overall_bbe * 10 / 100;

	if ( max_bbe <= auto_max_bbe )
	  auto_max_bbe = -1;
	else
	  p.set_MAX_BB_EVAL ( max_bbe );

      }

      // nothing is defined:
      else {
	p.set_MAX_BB_EVAL           ( -1 );
	p.set_MULTI_OVERALL_BB_EVAL ( -1 );

	p.check();
	auto_max_bbe = p.get_max_bb_eval();

	max_bbe = 2 * auto_max_bbe;
	p.set_MAX_BB_EVAL ( max_bbe );
      }
    }
  }
  else {
    p.set_MULTI_OVERALL_BB_EVAL ( command_line_bbe );

    p.check();
    auto_max_bbe = p.get_max_bb_eval();

    if ( command_line_bbe0 > 0 )
      p.set_MAX_BB_EVAL ( command_line_bbe0 );
    else {
      max_bbe = command_line_bbe * 10 / 100;
      if ( max_bbe <= auto_max_bbe )
	auto_max_bbe = -1;
      else
	p.set_MAX_BB_EVAL ( max_bbe );
    }
  }
}

/*------------------------------------------------------*/
/*  select a starting point from the cache              */
/*    . the point with the best f2 value is considered  */
/*    . the point may be feasible or not                */
/*------------------------------------------------------*/
const NOMAD::Eval_Point * select_x0 ( const NOMAD::Parameters & p0         ,
				      const NOMAD::Cache      & cache      ,
				      const NOMAD::Point      & lb         ,
				      const NOMAD::Point      & ub         ,
				      const NOMAD::Double     & lj         ,
				      const NOMAD::Double     & uj         ,
				      int                       cstr_index ,
				      bool                      check_feas   ) {
  
  NOMAD::Double f1 , f2 , best_f2 = 1e20 , best_f1 = 1e20;
  bool          is_feasible;
  int           i , i1 = *(p0.get_index_obj().begin()) , m = p0.get_bb_nb_outputs();
  const std::vector<NOMAD::bb_output_type> & bbot = p0.get_bb_output_type();
  const NOMAD::Eval_Point * x0 = NULL , * cur = cache.begin();


  while ( cur ) {
      
    if ( cur->get_eval_status() == NOMAD::EVAL_OK &&
	 check_bounds ( *cur , lb , ub ) ) {

      const NOMAD::Point & bbo = cur->get_bb_outputs();

      is_feasible = true;

      if ( check_feas )
	for ( i = 0 ; i < m ; ++i )
	  if ( i != cstr_index && NOMAD::bbot_is_constraint ( bbot[i] ) )
	    if ( bbo[i] > 0.0 ) {
	      is_feasible = false;
	      break;
	    }
	  
      if ( is_feasible ) {

	if ( cstr_index < m ) {
	  f2 = bbo[cstr_index];
	  if ( ( lj.is_defined() && f2 < lj ) ||
	       ( uj.is_defined() && f2 > uj      ) )
	    f2 = 1e20;
	}
	else
	  f2 = ( ( ( cstr_index - m ) % 2 ) ? 1.0 : -1.0 ) *
	    (*cur)[( cstr_index - m ) / 2];
      
	if ( f2 < best_f2 ) {
	  x0      = cur;
	  best_f2 = f2;
	}
	else if ( f2 == best_f2 ) {
	  f1 = bbo[i1];
	  if ( f1 < best_f1 ) {
	    x0      = cur;
	    best_f1 = f1;
	  }
	}
      }
    }
    cur = cache.next();
  }
  
  return x0;
}

/*-------------------------------------*/
/*           check the bounds          */
/*-------------------------------------*/
bool check_bounds ( const NOMAD::Point & x  ,
		    const NOMAD::Point & lb ,
		    const NOMAD::Point & ub   ) {

  int n = x.size();
  if ( n != lb.size() || n != ub.size() )
    return false;
  
  for ( int i = 0 ; i < n ; ++i )
    if ( !x[i].is_defined() ||
	 ( lb[i].is_defined() && x[i] < lb[i] ) ||
	 ( ub[i].is_defined() && ub[i] < x[i] )    )
      return false;

  return true;
}

/*-------------------------------------*/
/*          check the arguments        */
/*-------------------------------------*/
void check_args_2 ( const NOMAD::Parameters  & p          ,
		    int                        cstr_index ,
		    const NOMAD::Double      & lj         ,
		    const NOMAD::Double      & uj         ,
		    std::string              & error        ) {
  error.clear();

  std::vector<NOMAD::bb_output_type>
    bb_output_type = p.get_bb_output_type();
  
  const std::list<int> & l_index_obj = p.get_index_obj();

  if ( l_index_obj.size() != 1 ) {
    error = ": original problem is multiobjective";
    return;
  }

  int index_obj = *l_index_obj.begin();
  if ( cstr_index == index_obj ) {
    error = ": constraint index == objective index";
    return;
  }

  int m = bb_output_type.size();

  // constraint index is a true constraint:
  if ( cstr_index < m ) {
    
    if ( bb_output_type[cstr_index] != NOMAD::EB    &&
	 bb_output_type[cstr_index] != NOMAD::PB    &&
	 bb_output_type[cstr_index] != NOMAD::PEB_P &&
	 bb_output_type[cstr_index] != NOMAD::PEB_E &&
	 bb_output_type[cstr_index] != NOMAD::FILTER   ) {
      error = ": output[" + NOMAD::itos(cstr_index)
	+ "] is not a constraint";
      return;
    }
  }
    
  // constraint index is a bound:
  else {
    if ( cstr_index >= m + 2*p.get_dimension() ) {
      error = ": constraint index >= m+2n";
      return;
    }
    
    // index i of variable xi:
    int i = ( cstr_index - m ) / 2;

    // lower bound:
    if ( ( cstr_index - m ) % 2 == 0 ) {
      if ( !p.get_lb()[i].is_defined() ) {
	error = ": variable x" + NOMAD::itos(i) + " has no lower bound";
	return;
      }
    }
    
    // upper bound:
    else if ( !p.get_ub()[i].is_defined() ) {
      error = ": variable x" + NOMAD::itos(i) + " has no upper bound";
      return;
    }
  } 
}

/*-------------------------------------*/
/*  check and interpret the arguments  */
/*-------------------------------------*/
bool check_args_1 ( int                    argc              ,
		    char                ** argv              ,
		    const NOMAD::Display & out               ,
		    std::string          & param_file        ,
		    int                  & cstr_index        ,
		    std::string          & cache_file        ,
		    NOMAD::Double        & lj                ,
		    NOMAD::Double        & uj                ,
		    int                  & command_line_bbe  ,
		    int                  & command_line_bbe0   ) {

  // number of arguments:
  if ( argc < 4 || argc > 8 ) {
    display_usage ( out , argv );
    return false;
  }

  // parameters file:
  param_file = argv[1];

  // cache file:
  cache_file = argv[2];

  // constraint index:
  int i , n;
  {
    n = strlen ( argv[3] );
    for ( i = 0 ; i < n ; ++i )
      if ( !isdigit ( argv[3][i] ) ) {
	display_usage ( out , argv );
	out << "(error with cstr_index)" << endl << endl;
	return false;
      }
    cstr_index = atoi ( argv[3] );
    if ( cstr_index < 0 ) {
      display_usage ( out , argv );
      out << "(error with cstr_index: < 0)" << endl << endl;
      return false;
    }
  }

  // lj:
  lj.reset();
  if ( argc > 4 )
    if ( !lj.atof ( argv[4] ) || ( lj.is_defined() && lj > 0.0 ) ) {
      display_usage ( out , argv );
      out << "(error with lj)" << endl << endl;
      return false;
    }

  // uj:
  uj.reset();
  if ( argc > 5 )
    if ( !uj.atof ( argv[5] ) || ( uj.is_defined() && uj < 0.0 ) ) {
      display_usage ( out , argv );
      out << "(error with uj)" << endl << endl;
      return false;
    }

  if ( lj.is_defined() && uj.is_defined() && uj <= lj ) {
    display_usage ( out , argv );
    out << "(error: uj <= lj)" << endl << endl;
    return false;
  }

  // command_line_bbe:
  command_line_bbe = -1;
  if ( argc > 6 ) {
    n = strlen ( argv[6] );
    for ( i = 0 ; i < n ; ++i )
      if ( !isdigit ( argv[6][i] ) ) {
	display_usage ( out , argv );
	out << "(error with bb_eval)" << endl << endl;
	return false;
      }
    command_line_bbe = atoi ( argv[6] );
    if ( command_line_bbe <= 0 ) {
      display_usage ( out , argv );
      out << "(error with bb_eval: <= 0)" << endl << endl;
      return false;
    }
  }

  // command_line_bbe0:
  command_line_bbe0 = -1;
  if ( argc > 7 ) {
    n = strlen ( argv[7] );
    for ( i = 0 ; i < n ; ++i )
      if ( !isdigit ( argv[7][i] ) ) {
	display_usage ( out , argv );
	out << "(error with init_bb_eval)" << endl << endl;
	return false;
      }
    command_line_bbe0 = atoi ( argv[7] );
    if ( command_line_bbe0 <= 0 ) {
      display_usage ( out , argv );
      out << "(error with init_bb_eval: <= 0)" << endl << endl;
      return false;
    }
  }

  return true;
}

/*----------------------------------------*/
/*         display the objectives         */
/*----------------------------------------*/
void display_objectives ( const NOMAD::Display    & out        ,
			  const NOMAD::Parameters & p0         ,
			  const NOMAD::Parameters & p          ,
			  int                       cstr_index ,
			  const NOMAD::Double     & lj         ,
			  const NOMAD::Double     & uj           ) {
  out << endl
      << "Sensitivity analysis with biobjective optimization:"
      << endl << endl;

  const std::list<int> & index_obj = p.get_index_obj();
  if ( index_obj.size() != 2 )
    throw NOMAD::Exception ( "detailed_analysis.cpp" , __LINE__ ,
			     ": new Parameters object is not biobjective" );
  int k = 0 , i0 = *p0.get_index_obj().begin() , m = p0.get_bb_nb_outputs();
  
  std::list<int>::const_iterator it , end = index_obj.end();
  
  for ( it = index_obj.begin() ; it != end ; ++it ) {
	
    out << "\tobjective f" << ++k << ": ";
    if ( *it == i0 )
      out << "original objective (output #" << *it << ")";
    else if ( *it < m ) {
      out << "constraint c" << *it << "<=0 ";
      if ( lj.is_defined() && uj.is_defined() )
	out << "replaced by " << lj << "<=c" << *it << "<=" << uj;
      else if ( lj.is_defined() )
	out << "replaced by " << lj << "<=c" << *it;
      else if ( uj.is_defined() )
	out << "replaced by c" << *it << "<=" << uj;
      else
	out << "removed";
    }
    else {
      
      int i = ( cstr_index - m ) / 2;
      
      if ( ( cstr_index - m ) % 2 == 0 ) {
	const NOMAD::Double & v = p0.get_lb()[i];
	out << "-x" << i << ": lower bound " << v << "-x" << i << "<=0 ";
	if ( lj.is_defined() && uj.is_defined() )
	  out << "replaced by " << lj << "<=" << v << "-x" << i << "<=" << uj;
	else if ( lj.is_defined() )
	  out << "replaced by " << lj << "<=" << v << "-x" << i;
	else if ( uj.is_defined() )
	  out << "replaced by " << v << "-x" << i << "<=" << uj;
	else
	  out << "removed";
      }
      else {
	const NOMAD::Double & v = p0.get_ub()[i];
	out << "x" << i << ": upper bound x" << i << "-" << v << "<=0 ";
	if ( lj.is_defined() && uj.is_defined() )
	  out << "replaced by " << lj << "<=x" << i << "-" << v << "<=" << uj;
	else if ( lj.is_defined() )
	  out << "replaced by " << lj << "<=x" << i << "-" << v;
	else if ( uj.is_defined() )
	  out << "replaced by x" << i << "-" << v << "<=" << uj;
	else
	  out << "removed";
      }
    }
    out << endl;
  }
}

/*-----------------*/
/*  display_usage  */
/*-----------------*/
void display_usage ( const NOMAD::Display & out , char ** argv ) {
  out << endl << "usage: " << argv[0]
      << " param_file cache_file cstr_index lj uj bb_eval init_bb_eval"
      << endl << endl
      << "       param_file: parameters file" << endl
      << "       cache_file: cache file to be treated with cache_inspect" << endl
      << "       cstr_index: constraint index j in [0;m+2n-1]" << endl
      << "                   m corresponds to the number of outputs" << endl
      << "                   n corresponds to the number of variables" << endl
      << "                   j<m corresponds to constraint cj<=0" << endl
      << "                   j=m+2i corresponds to the lower bound of variable i"
      << " in [0;n-1]" << endl
      << "                   j=m+2i+1 corresponds to the upper bound of variable i"
      << " in [0;n-1]" << endl
      << "           lj, uj: optional, may be - or INF"                             << endl
      << "                   constraint cj<=0 is replaced by cj<=e with lj<=e<=uj"        << endl
      << "                   lower bound lb-xi<=0 is replaced by lb-xi<=e with lj<=e<=uj" << endl
      << "                   upper bound xi-ub<=0 is replaced by xi-ub<=e with lj<=e<=uj" << endl
      << "          bb_eval: optional, max number of blackbox evaluations"      << endl
      << "                   if none specified, param_file is considered"       << endl
      << "     init_bb_eval: optional, max number of blackbox evaluations for "
      << "the two first MADS runs" << endl
      << "                   if none specified, 10% of the max number of evaluations "
      << "is considered" << endl
      << endl << endl;
}
