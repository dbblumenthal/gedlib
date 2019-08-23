/*--------------------------------------------------*/
/*  NOMAD_PLOT -- single and bi-objective versions  */
/*--------------------------------------------------*/
#include "nomad.hpp"
using namespace std;
using namespace NOMAD;

class Interface_Graphique
{
  public:

  Interface_Graphique(){};
  ~Interface_Graphique(){};
    
    void create_multi()
    {
      ofstream deb("ini.txt", ios::out | ios::trunc);
      {
	deb << "0" <<"\n" ;
	deb.close();
      }
      rename("ini.txt", "init.txt");
    }

    void create_single()
    {
      ofstream deb("ini.txt", ios::out | ios::trunc);
      {
	deb << "1" <<"\n";
	deb.close();
      }
      rename("ini.txt", "init.txt");
    }

    void new_single(int bbe, double f)
    {
      FILE* fp = NULL;
      fp = fopen( "out.txt", "rb" );
      while(!fp)
      {
        fp = fopen( "out.txt", "rb" );
      }
      fclose(fp);

      ofstream deb("out.txt", ios::out | ios::trunc);
      {
	deb << bbe << "\n"<< f;
	deb.close();
      }

      rename("out.txt", "in.txt");
    }

    void new_multi(int nb_cache_pts, double* cache_f1, double* cache_f2,
                    int nb_pareto_pts, double* pareto_f1, double* pareto_f2) const
    {
      FILE* fp = NULL;
      fp = fopen( "out.txt", "rb" );
      while(!fp)
      {
        fp = fopen( "out.txt", "rb" );
      }
      fclose(fp);

      ofstream deb("out.txt", ios::out | ios::trunc);
      {
	deb << nb_cache_pts << "\n";
        for (int i = 0; i < nb_cache_pts; ++i)
	{
          deb << cache_f1[i] << "\n"<< cache_f2[i] << "\n";
        }
	deb << nb_pareto_pts << "\n";
        for (int i = 0; i < nb_pareto_pts; ++i)
	{
          deb << pareto_f1[i] << "\n"<< pareto_f2[i] << "\n";
        }
        deb.close();
      }

      rename("out.txt", "in.txt");
    }
};


/*----------------------------------------------------*/
/*  custom evaluator used to define update_success()  */
/*                   (bi-objective)                   */
/*----------------------------------------------------*/
class My_Multi_Obj_Evaluator : public Multi_Obj_Evaluator {

private:

  Double                  _h_min;
  set<const Eval_Point *> _cache_pts;
  

public:

  Interface_Graphique inter;
  // constructor:
  My_Multi_Obj_Evaluator ( const Parameters & p );

  // destructor:
  virtual ~My_Multi_Obj_Evaluator ( void ) {}

  // plot the points:
  void plot_points ( int      nb_cache_pts  ,
		     double * cache_f1      ,
		     double * cache_f2      ,
		     int      nb_pareto_pts ,
		     double * pareto_f1     ,
		     double * pareto_f2       ) const;

  // function that is invoked after each MADS run:
  virtual void update_mads_run ( const Stats             & stats         ,
 				 const Evaluator_Control & ev_control    ,
 				 const Barrier           & true_barrier  ,
 				 const Barrier           & sgte_barrier  ,
 				 const Pareto_Front      & pareto_front    ); 
};

/*-----------------------*/
/*  constructor (bi-obj) */
/*-----------------------*/
My_Multi_Obj_Evaluator::My_Multi_Obj_Evaluator ( const Parameters & p ) :
  Multi_Obj_Evaluator ( p             ) ,
  _h_min              ( p.get_h_min() ) {
  inter.create_multi();
}

/*------------------------------------------------*/
/*  function that is invoked after each MADS run  */
/*  (bi-obj)                                      */
/*------------------------------------------------*/
void My_Multi_Obj_Evaluator::update_mads_run ( const Stats             & stats         ,
					       const Evaluator_Control & ev_control    ,
					       const Barrier           & true_barrier  ,
					       const Barrier           & sgte_barrier  ,
					       const Pareto_Front      & pareto_front    ) {

  // objective indexes:
  int i1 = Multi_Obj_Evaluator::get_i1();
  int i2 = Multi_Obj_Evaluator::get_i2();

  // browse all the cache points:
  const Cache              & cache = ev_control.get_cache();
  const Eval_Point         * cur   = cache.begin();
  size_t                     size_before;
  vector<const Eval_Point *> cache_pts;
  vector<const Eval_Point *> pareto_pts;

  while (cur) {
    if ( cur->is_feasible(_h_min) && cur->is_eval_ok() ) {
      size_before = _cache_pts.size();     
      _cache_pts.insert ( cur );
      if ( _cache_pts.size() > size_before )
	cache_pts.push_back ( cur );
    }
    cur = cache.next();
  }

  int      nb_cache_pts = cache_pts.size();
  double * cache_f1     = new double [nb_cache_pts];
  double * cache_f2     = new double [nb_cache_pts];

  int k;

  for ( k = 0 ; k < nb_cache_pts ; ++k ) {
    cache_f1[k] = cache_pts[k]->get_bb_outputs()[i1].value();
    cache_f2[k] = cache_pts[k]->get_bb_outputs()[i2].value();
  }

  // browse the pareto points:
  int      nb_pareto_pts = pareto_front.size();
  double * pareto_f1     = new double [nb_pareto_pts];
  double * pareto_f2     = new double [nb_pareto_pts];

  k   = 0;
  cur = pareto_front.begin();
  while ( cur ) {

    pareto_f1[k] = cur->get_bb_outputs()[i1].value();
    pareto_f2[k] = cur->get_bb_outputs()[i2].value();

    cur = pareto_front.next();
    ++k;
  }

  // call the plot function:
  plot_points ( nb_cache_pts  , cache_f1  , cache_f2 ,
		nb_pareto_pts , pareto_f1 , pareto_f2 );

  delete [] cache_f1;
  delete [] cache_f2;
  delete [] pareto_f1;
  delete [] pareto_f2;
}

/*----------------------------*/
/*  plot points (multi-obj)  */
/*----------------------------*/
void My_Multi_Obj_Evaluator::plot_points ( int      nb_cache_pts  ,
					   double * cache_f1      ,
					   double * cache_f2      ,
					   int      nb_pareto_pts ,
					   double * pareto_f1     ,
					   double * pareto_f2       )  const{
  inter.new_multi(nb_cache_pts, cache_f1, cache_f2, nb_pareto_pts, pareto_f1, pareto_f2); 
}

/*----------------------------------------------------*/
/*  custom evaluator used to define update_success()  */
/*                 (single-objective)                 */
/*----------------------------------------------------*/
class My_Evaluator : public Evaluator {

private:
Interface_Graphique inter;

public:

  // constructor:
  My_Evaluator  ( const Parameters & p );

  // destructor:
  virtual ~My_Evaluator ( void ) {}

  // plot a new success:
  void plot_success ( int bbe , const Double & f );

  // function that is invoked at each success:
  virtual void update_success ( const Stats & stats , const Eval_Point & x );

};

/*---------------------------*/
/*  constructor (single-obj) */
/*---------------------------*/
My_Evaluator::My_Evaluator  ( const Parameters & p ) : Evaluator ( p ) {
  inter.create_single();
}

/*--------------------------------------------*/
/*  function that is invoked at each success  */
/*  (single-obj)                              */
/*--------------------------------------------*/
void My_Evaluator::update_success ( const Stats & stats , const Eval_Point & x ) {
  if ( x.is_feasible ( _p.get_h_min() ) )
    plot_success ( stats.get_bb_eval()+1 , x.get_f() );
}

/*-----------------------------------*/
/*  plot a new success (single-obj)  */
/*-----------------------------------*/
void My_Evaluator::plot_success ( int bbe , const Double & f ) {
  inter.new_single ( bbe, f.value() );
}

/*-------------------------*/
/*  tells the GUI to wait  */
/*-------------------------*/
void GUI_wait ( void ) {
  ofstream fin("stop.txt", ios::out | ios::trunc);
  {
    fin.close();
  }
}

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

    // usage:
    if ( argc < 2 ) {
      display_usage ( cerr );
      return EXIT_SUCCESS;
    }

    // parameters file:
    string param_file_name = argv[1];
    string opt             = param_file_name;
    NOMAD::toupper(opt);

    // display version if option '-v' has been specified:
    if ( opt == "-V" ) {
      display_version ( out );
      return EXIT_SUCCESS;
    }

    // display info if option '-i' has been specified:
    if ( opt == "-I" || opt == "-INFO" ) {
      display_info  ( out );
      display_usage ( out );
      return EXIT_SUCCESS;
    }
    
    // parameters creation:
    Parameters p ( out );

    // display help on parameters if option '-h' has been specified:
    if ( opt == "-H" || opt == "-HELP" ) {
      p.help( (argc>2) ? argv[2] : "all" );
      return EXIT_SUCCESS;
    }

    // read parameters file:
    p.read ( param_file_name );

    // parameters check:
    p.check();

    // display NOMAD info:
    if ( p.get_display_degree() > 1 )
      display_info ( out );

    // parameters display:
    if ( p.get_display_degree() > 2 )
      out << endl << "parameters:" << endl << p << endl;

    // single-objective:
    if ( p.get_nb_obj() == 1 ) {

      // custom evaluator:
      My_Evaluator ev ( p );

      // algorithm creation and execution:
      Mads mads ( p , &ev );
      mads.run();

      // plot the last point:
      const Eval_Point * bf = mads.get_best_feasible();
      if ( bf )
	ev.plot_success ( mads.get_stats().get_bb_eval() , bf->get_f() );
    }

    // bi-objective:
    else {

      // custom evaluator:
      My_Multi_Obj_Evaluator ev ( p );

      // algorithm creation and execution:
      Mads mads ( p , &ev );
      mads.multi_run();
    }


  }
  catch ( exception & e ) {
    cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
  }

  GUI_wait();


  Slave::stop_slaves ( out );
  end();

  return EXIT_SUCCESS;
}

/*------------------------------------------*/
/*            display NOMAD version         */
/*------------------------------------------*/
void NOMAD::display_version ( const NOMAD::Display & out )
{
#ifdef USE_MPI
  if ( !NOMAD::Slave::is_master() )
    return;
#endif
  out << std::endl << "NOMAD - version " << NOMAD::VERSION << " - www.gerad.ca/nomad"
      << std::endl << std::endl;
}

/*------------------------------------------*/
/*          display NOMAD information       */
/*------------------------------------------*/
void NOMAD::display_info ( const NOMAD::Display & out )
{
#ifdef USE_MPI
  if ( !NOMAD::Slave::is_master() )
    return;
#endif
  NOMAD::display_version ( out );
  out << NOMAD::open_block ( "Copyright (C) 2001-2010" )
      << "Mark A. Abramson     - The Boeing Company"                 << std::endl
      << "Charles Audet        - Ecole Polytechnique de Montreal"    << std::endl
      << "Gilles Couture       - Ecole Polytechnique de Montreal"    << std::endl
      << "John E. Dennis, Jr.  - Rice University"                    << std::endl
      << "Sebastien Le Digabel - Ecole Polytechnique de Montreal"    << std::endl
      << NOMAD::close_block()
      << std::endl
      << "Funded in part by AFOSR and Exxon Mobil."                  << std::endl
      << std::endl
      << "License   : \'" << NOMAD::LGPL_FILE       << "\'" << std::endl
      << "User guide: \'" << NOMAD::USER_GUIDE_FILE << "\'" << std::endl
      << "Examples  : \'" << NOMAD::EXAMPLES_DIR    << "\'" << std::endl
      << "Tools     : \'" << NOMAD::TOOLS_DIR       << "\'" << std::endl
      << std::endl
      << "Please report bugs to nomad@gerad.ca"
      << std::endl;
}

/*------------------------------------------*/
/*             display NOMAD usage          */
/*------------------------------------------*/
void NOMAD::display_usage ( const NOMAD::Display & out )
{
#ifdef USE_MPI
  if ( !NOMAD::Slave::is_master() )
    return;
  out << std::endl
      << "Run NOMAD.MPI: mpirun -np p nomad.MPI parameters_file" << std::endl
      << "Info         : nomad -i"                    << std::endl
      << "Help         : nomad -h keyword (or 'all')" << std::endl
      << "Version      : nomad -v"                    << std::endl
      << std::endl;  
#else
  out << std::endl
      << "Run NOMAD: nomad parameters_file"          << std::endl
      << "Info     : nomad -i"                       << std::endl
      << "Help     : nomad -h keyword(s) (or 'all')" << std::endl
      << "Version  : nomad -v"                       << std::endl
      << std::endl; 
#endif
}
