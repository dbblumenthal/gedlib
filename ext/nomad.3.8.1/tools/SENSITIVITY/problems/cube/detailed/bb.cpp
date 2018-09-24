#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
using namespace std;

const double X[] = { 5 , 6 , 2 };
//const double X[] = { 10 , 15 , 4 };

int main ( int argc , char ** argv ) {

  double z = 1e+20 , c0 = 1e+20; //, c2 = 1e+20 , c3 = 1e+20 , c0 = 1e+20;

  if ( argc != 2 ) {
	  cout << c0 << " " << z << endl;
    return 1;
  }

  double x[3];

  ifstream in ( argv[1] );
  if ( in.fail() ) {
    cout << c0 << " " << z << endl;
    return 1;
  }

  int i;

  for ( i = 0 ; i < 3 ; i++ )
    in >> x[i];

  if ( in.fail() ) {
	  cout << c0 << " " << z << endl;
    in.close();
    return 1;
  }

  in.close();

  /*-------------------------------------------------------------------*/

  /*c1 = x[0] - 2.0;
  c2 = x[1] - 2.0;
  c3 = x[2] - 2.0;*/

  c0 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 9.0; 

  z = pow (x[0]-X[0] , 2.0 ) + pow (x[1]-X[1] , 2.0 ) + pow (x[2]-X[2] , 2.0 );

  cout << setprecision(16);
//    cout << c0 << " " << c1 << " " << c2 << " " << c3 << " " << z << endl;
    // cout << c0 << " " << c2 << " " << z << endl;
  cout << c0 << " " << z << endl;

  return 0;
}
