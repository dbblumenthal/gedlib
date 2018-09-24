#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
using namespace std;

const double GAMMA[13] =
  { 0.0137 ,
    0.0274 ,
    0.0434 ,
    0.0866 ,
    0.137  , 
    0.274  ,
    0.434  , 
    0.866  ,
    1.37   ,
    2.74   ,
    4.34   ,
    5.46   ,
    6.88     };

const double ETA[13] =
  { 3220 ,
    2190 ,
    1640 ,
    1050 ,
    766 ,
    490 ,
    348 ,
    223 ,
    163 ,
    104 ,
    76.7 ,
    68.1 ,
    58.2   };

int main ( int argc , char ** argv ) {

  double f = 1e20;
  double eta_0 , lambda , beta;

  if ( argc >= 2 ) {
    
    ifstream in ( argv[1] );

    in >> eta_0 >> lambda >> beta;
    
    eta_0  *= 955.4;
    lambda *= 38.5;
    beta   *= 0.035;

    if ( in.fail() )
      f = 1e20;
    else {
      f = 0.0;
      for ( int i = 0 ; i < 13 ; ++i )
        f += fabs ( ETA[i] -
		    eta_0 *
		    pow ( 1 + lambda*lambda*GAMMA[i]*GAMMA[i] ,
			  (beta-1)/2.0 ) );
    }
    
    in.close();
  }

/* TOTO  A VIRER 
  srand(getpid()+time(0));
 
  const int N = rand()%1000;

 
  for ( int i = 0 ; i < N ; ++i )
	for ( int j = 0 ; j < N ; ++j )
	double r = pow(rand()%1000 , 3 );
*/
 
  cout << setprecision(15);
  cout << f << endl; // " " << -f /* -275*/ << endl;
  
  return 0;
}
