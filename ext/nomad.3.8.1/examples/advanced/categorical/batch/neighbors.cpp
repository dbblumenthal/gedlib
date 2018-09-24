/*--------------------------------------*/
/*  construct the extended poll points  */
/*      (categorical neighborhoods)     */
/*--------------------------------------*/
#include <iostream>
#include <fstream>
using namespace std;

int main ( int argc , char ** argv ) {

  if ( argc != 2 )
    return 1;

  ifstream in ( argv[1] );
  
  double t0 , v0 , t1 , v1;

  in >> t0 >> v0 >> t1 >> v1;

  in.close();

  if ( in.fail() )
    return 1;

  int t2 = static_cast<int> ( 3 - t0 - t1 );

  // neighbor #1:
  cout << t2 << " " << v0 << " " << t1 << " " << v1 << endl;

  // neighbor #2:
  cout << t0 << " " << v0 << " " << t2 << " " << v1 << endl;


  return 0;
}
