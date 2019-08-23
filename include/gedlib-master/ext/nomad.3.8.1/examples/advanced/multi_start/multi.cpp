// This program is not intended to be linked with the MPI NOMAD library.

#include "nomad.hpp"
using namespace std;
using namespace NOMAD;

const double LB = 0.0;  // lower bounds
const double UB = 4.5;  // upper bounds

// LH search used to generate a list of starting points:
void LH_x0 ( int n , int p , vector<Point *> & x0_pts );
void LH_values_for_var_i ( int     ind ,
			   int     p   ,
			   Point & x     );

/*--------------------------------------*/
/*            custom evaluator          */
/*--------------------------------------*/
class My_Evaluator : public Evaluator {

private:

  int      _NBC;
  int      _M;
  double * _x;
  double * _y;

public:

  // constructor:
  My_Evaluator ( const Parameters & param , int NBC , int M )
    : Evaluator ( param            ) ,
      _NBC      ( NBC              ) ,
      _M        ( M                ) ,
      _x        ( new double [NBC] ) ,
      _y        ( new double [NBC] )   {
    _y[0] = _y[1] = _x[2] = 0.0;
  }

  // destructor:
  ~My_Evaluator ( void ) { delete [] _x; delete [] _y; }

  // eval_x:
  bool eval_x ( Eval_Point   & x          ,
		const Double & h_max      ,
		bool         & count_eval   ) const;
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

    // usage:
    if ( argc != 4 ) {
      cerr << "\nusage: multi param.txt nb_circles nb_mads_runs\n\n";
      return 1;
    }

    int NBC = atoi(argv[2] );

    if ( NBC < 4 ) {
      cerr << "\nthe number of circles must be > 3\n\n";
      return 1;
    }

    int   N = 2*NBC-3;
    //int   M = ( NBC * ( NBC + 1 ) - 8 ) / 2;
    int   M = 2;

    // list of x0 points (LH strategy is used):
    // -----------------

    vector<Point *> x0_pts;

    // srand ( static_cast<int> ( time(0) ) );
    srand(0);
    rand();

    int i;
    int nb_mads_runs = atoi ( argv[3] );

    LH_x0 ( N , nb_mads_runs , x0_pts );

    // read best_x.txt:
    ifstream fin ( "best_x.txt");

    if ( !fin.fail() )
      for ( int i = 0 ; i < N ; ++i )
        fin >> (*x0_pts[0])[i];

    fin.close();

    // display all starting points:
    out << endl;
    for ( int j = 0 ; j < nb_mads_runs ; ++j )
      out << "starting point # " << j << ": ( " << *x0_pts[j] << " )" << endl;
    out << endl;

    // parameters creation:
    // --------------------
    Parameters param ( out );

    param.set_DIMENSION (N);
    vector<bb_output_type> bbot (M+1);
    for ( i = 0 ; i < M ; ++i )
      bbot[i] = PB;
    bbot[M] = OBJ;
    param.set_BB_OUTPUT_TYPE ( bbot );

    Point lb ( N , LB );
    Point ub ( N , UB );
    param.set_LOWER_BOUND ( lb );
    param.set_UPPER_BOUND ( ub );

    param.read ( argv[1] );

    param.set_DISPLAY_DEGREE (0);

    param.set_X0 ( *x0_pts[0] );

    // parameters check:
    param.check();

    // out << param << endl;

    // custom evaluator creation:
    My_Evaluator ev ( param , NBC , M );

    const Eval_Point * cur_x;
    Point              best_x (N);
    Double             best_f = INF , worst_f = 0.0 , avg_f = 0.0;

    // MADS runs:
    // ----------
    int bbe = 0;
    i = 0;
    while ( true ) {

      // algorithm creation:
      Mads mads ( param , &ev );
      mads.run();

      bbe += mads.get_cache().size();

      // displays and remember the best point:
      out << "run #" << setw(2) << i << ": ";
      cur_x = mads.get_best_feasible();
      if ( cur_x ) {

	out << "f=" << cur_x->get_f() << endl;

        if ( cur_x->get_f() < best_f ) {

	  best_f = cur_x->get_f();
	  best_x = *cur_x;
	}

	if ( cur_x->get_f() > worst_f )
	  worst_f = cur_x->get_f();

	avg_f += cur_x->get_f();

      }
      else
	out << "NULL" << endl;

      if ( ++i == nb_mads_runs )
	break;

      // preparation of next run:
      mads.reset();
      param.reset_X0();
      param.set_X0 ( *x0_pts[i] );
      param.check();
    }

    // display the solution:
    out << endl << "bb eval : " << bbe << endl
	<< "best    : " << best_f;
    out << endl
	<< "worst   : " << worst_f << endl
	<< "solution: ";
    out << "x = ( ";
    best_x.display ( out , " " , -1 , -1 );
    out << " ) f(x) = " << best_f.value();
    out << endl << endl;

    ofstream fout ( "best_x.txt" );
    fout << setprecision(32);
    best_x.display ( fout , " " , -1 , -1 );
    fout.close();

    // delete x0 points:
    for ( i = 0 ; i < nb_mads_runs ; ++i )
      delete x0_pts[i];
  }
  catch ( exception & e ) {
    cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
  }

  end();

  return EXIT_SUCCESS;
}


/*---------------------------------*/
/*  methods of class My_Evaluator  */
/*---------------------------------*/

// eval_x:
bool My_Evaluator::eval_x ( Eval_Point   & ep         ,
			    const Double & h_max      ,
			    bool         & count_eval   ) const {
	
  count_eval = true;

  int epi = 0;

  _x[0] = ep[epi++].value();  _y[0] = 0.0;
  _x[1] = ep[epi++].value();  _y[1] = 0.0;
  _y[2] = ep[epi++].value();  _x[2] = 0.0;
  _x[3] = ep[epi++].value();
  _y[3] = ep[epi++].value();

  int ic = 0 , i , j;
  double constraint = 0.0;

  for ( i = 4 ; i < _NBC ; ++i ) {
    _x[i] = ep[epi++].value();
    _y[i] = ep[epi++].value();
  }

  ep.set_bb_output ( ic++ , constraint );

  double dist , distmax = 0.0, avg_dist = 0.0;
  constraint = 1.0;

  /*for ( i = 0 ; i < _NBC ; ++i ) {
    for( j = i+1 ; j < _NBC ; ++j ) {
    dist = pow(_x[i]-_x[j],2) + pow(_y[i]-_y[j],2);
    avg_dist += sqrt(dist);
    if ( dist > distmax )
    distmax = dist;
    if( dist < 1 ) constraint *= dist;
    }
    }
    ep.set_bb_output ( ic++ , 1.0-constraint );
    ep.set_bb_output ( _M , sqrt(distmax) );
  */

  for ( i = 0 ; i < _NBC-1 ; ++i ) {
    // out << "  ( "<< _x[i] << ",  " << _y[i] << " ) \n";
    for( j = i+1 ; j < _NBC ; ++j ) {
      dist = pow(_x[i]-_x[j],2) + pow(_y[i]-_y[j],2);
      if ( dist > distmax )
	distmax = dist;
      if( dist < 1 ) {
	constraint *= sqrt(dist);
      }
      else if( dist > 14)
	avg_dist += sqrt(dist);
      // out << i << " " << j << " . "  << sqrt(dist) << endl;
    }
  }
  // out << distmax << " " <<  avg_dist << " " << constraint << "                 ";

  if (constraint < 0.9999)
    constraint = 1001.0-constraint;
  else
    constraint = sqrt(distmax) + avg_dist/(10.0*_NBC);

  // out << constraint << "\n";

  ep.set_bb_output ( ic++ , 0.0 );
  ep.set_bb_output ( _M , constraint);

  return true;
}

/*----------------------------------------*/
/*  LH search used to generate x0 points  */
/*----------------------------------------*/
void LH_x0 ( int n , int p , vector<Point *> & x0_pts ) {

  // pts contains n points of dimension p: each of these points contains
  // p different values for each variable:
  Point ** pts = new Point * [n] , * x;

  int pm1 = p - 1 , i;

  // creation of p search points:
  for ( int k = 0 ; k < p ; ++k ) {

    x = new Point (n);

    for ( i = 0 ; i < n ; ++i ) {

      if ( k==0 ) {
	pts[i] = new Point(p);

	LH_values_for_var_i ( i , p , *pts[i] );

      }

      (*x)[i] = (*pts[i])[k];

      if ( k == pm1 )
	delete pts[i];
    }

    x0_pts.push_back ( x );

  }

  delete [] pts;
}


void LH_values_for_var_i ( int     ind ,
			   int     p   ,
			   Point & x     ) {

  Random_Pickup rp(p);
  int    i;
  double w = (UB - LB)/p;
  Double v;

  for ( i = 0 ; i < p ; ++i ) {
    v = LB + ( i + rand()/D_INT_MAX ) * w;
    x[rp.pickup()] = v;
  }
}
