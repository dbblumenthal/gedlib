/*---------------------------------------------------------------*/
/*                    NOMAD for a DLL-black-box                  */
/*---------------------------------------------------------------*/
#include "nomad.hpp"
#include <windows.h>
using namespace std;
using namespace NOMAD;

// black-box outputs (m):
// ----------------------
//    . the first one is the objective function
//    . all the others (1 to m-1) are PB constraints with format g(x) <= 0


/*------------------------------------------*/
/*               DLL functions              */
/*------------------------------------------*/
typedef int  (_stdcall *GET_NM)( void );
typedef void (_stdcall *EVAL_X)( double * , double * , bool * );
typedef void (_stdcall *INIT  )( void );

GET_NM get_n_dll;
GET_NM get_m_dll;
EVAL_X eval_x_dll;
//INIT init_dll;
//INIT clear_dll;

// simulate the DLL functions (for tests):
/*
void init  ( void ) {}
void clear ( void ) {}
int  get_n_dll  ( void ) { return 5; }
int  get_m_dll  ( void ) { return 3; }
void eval_x_dll ( double * x , double * fx , bool & cnt_eval ) {
  double c1 = 0.0 , c2 = 0.0;
  for ( int i = 0 ; i < 5 ; ++i ) {
    c1 += pow ( x[i]-1 , 2 );
    c2 += pow ( x[i]+1 , 2 );
  }
  fx[0] = x[4];
  fx[1] = c1 - 25;
  fx[2] = 25 - c2;
  cnt_eval = true;
}
*/

/*--------------------------------------*/
/*            custom evaluator          */
/*--------------------------------------*/
class My_Evaluator : public Multi_Obj_Evaluator {

private:

  int      _n;
  int      _m;
  double * _px;
  double * _fx;

public:

  // ctor:
  My_Evaluator ( const Parameters & p , int n , int m )
    : Multi_Obj_Evaluator ( p ) ,
      _n        ( n               ) ,
      _m        ( m               ) ,
      _px       ( new double [_n] ) ,
      _fx       ( new double [_m] )   {}

  // dtor:
  ~My_Evaluator ( void ) { delete [] _px; delete [] _fx; }

  // eval_x:
  bool eval_x ( Eval_Point          & x        ,
		const NOMAD::Double & h_max    ,
		bool                & cnt_eval   ) const;
};

// eval_x:
bool My_Evaluator::eval_x ( Eval_Point          & x        ,
			    const NOMAD::Double & h_max    ,
			    bool                & cnt_eval   ) const {
  int i;
  for ( i = 0 ; i < _n ; ++i )
    _px[i] = x[i].value();   
  eval_x_dll ( _px , _fx , &cnt_eval );
  for ( i = 0 ; i < _m ; ++i )
    x.set_bb_output ( i , _fx[i] );
  return true;
}

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv ) {

  // use:
  if ( argc != 2 ) {
    cerr << "\nuse: nomad_for_DLL.exe param.txt\n\n";
    return 1;
  } 

  // NOMAD initializations:
  begin ( argc , argv );

  // display:
  NOMAD::Display out ( std::cout );
  out.precision ( NOMAD::DISPLAY_PRECISION_STD );

  // random seed:
  unsigned seed = 0;
  srand(seed);

  // define DLL functions:
  // get handle to dll:
  HINSTANCE hdll = LoadLibrary ( TEXT("bb.dll") );
  if ( !hdll ) {
    cerr << "ERROR: unable to open bb.dll" << endl;
	FreeLibrary ( hdll );
	return 1;
  }

  // get function pointers:
  get_n_dll = (GET_NM)GetProcAddress(hdll, "GET_N");
  if ( !get_n_dll ) {
    cerr << "ERROR: unable to find GET_N function in dll" << endl;
	FreeLibrary ( hdll );
	return 1;
  }

  get_m_dll = (GET_NM)GetProcAddress(hdll, "GET_M");
  if ( !get_m_dll ) {
    cerr << "ERROR: unable to find GET_M function in dll" << endl;
	FreeLibrary ( hdll );
	return 1;
  }

  eval_x_dll = (EVAL_X)GetProcAddress(hdll, "EVAL_X");
  if ( !eval_x_dll ) {
    cerr << "ERROR: unable to find EVAL_X function in dll" << endl;
	FreeLibrary ( hdll );
	return 1;
  }

 /* init_dll = (INIT)GetProcAddress(hdll, "INIT");
  if ( !init_dll ) {
    cerr << "ERROR: unable to find INIT function in dll" << endl;
	FreeLibrary ( hdll );
	return 1;
  }

  clear_dll = (INIT)GetProcAddress(hdll, "CLEAR");
  if ( !clear_dll ) {
    cerr << "ERROR: unable to find CLEAR function in dll" << endl;
	FreeLibrary ( hdll );
	return 1;
  }*/

  // MADS:
  try {

    /* init_dll(); */

    int n = get_n_dll();
    int m = get_m_dll();

    // parameters creation:
    // --------------------
    Parameters p ( out );

    // dimension:
    p.set_DIMENSION ( n );

    // definition of output types:
    vector<bb_output_type> bbot (m);
    bbot[0] = bbot[1] = OBJ;
    for ( int i = 2 ; i < m ; ++i )
      bbot[i] = PB;
    p.set_BB_OUTPUT_TYPE ( bbot );

    // read parameters file:
    p.read ( argv[1] );

    // parameters check:
    p.check();

    // display parameters:
    // out << p << endl;

    // custom evaluator creation:
    My_Evaluator ev ( p , n , m );

    // algorithm creation and execution:
    Mads mads ( p , &ev );
    mads.multi_run();

    // algorithm display:
    // out << mads << endl;

    /* clear_dll(); */
  }
  catch ( exception & e ) {
    cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
  }

  // Release the Dll:
  FreeLibrary ( hdll );

  NOMAD::Slave::stop_slaves ( out );
  NOMAD::end();

  return EXIT_SUCCESS;
}
