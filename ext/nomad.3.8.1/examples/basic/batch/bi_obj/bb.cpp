// Deb Discontinuous, alpha=2, q=4
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

#define ALPHA  2.0
#define Q      4.0
#define PIT2   6.2831853071796

int main ( int argc , char ** argv ) {

  double x1 , x2;

  if ( argc == 2 ) {

    ifstream in ( argv[1] );

    in >> x1 >> x2;

    if ( in.fail() ) {
      cout << 1e20 << " " << 1e20 << endl;
      in.close();
      return 1;
    }

    in.close();
  }

  else if ( argc == 3 ) {
    x1 = atof(argv[1]);
    x2 = atof(argv[2]);
  }
  else  {
    cout << 1e20 << " " << 1e20 << endl;
    return 1;
  }

  double f1  = x1;
  double g   = 1 + 10.0 * x2;
  double f1g = f1 / g;
  double f2  = g - g * pow(f1g,ALPHA) - f1 * sin ( PIT2 * Q * f1 );


  cout.precision(12);

  cout << f1 << " " << f2 << endl;
  
  return 0;
}

