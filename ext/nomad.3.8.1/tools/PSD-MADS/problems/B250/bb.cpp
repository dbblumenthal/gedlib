/*--------------------------------------------------------*/
/* description : black-box described in the MoVars paper  */
/*               (version G >= 250)                       */
/* author      : Sebastien Le Digabel                     */
/* date        : 2007-10-31                               */
/*--------------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const double SQRT_A = pow(10.0,-2.5);
const double E1     = exp(-0.1);

const double H  = 0.0625;
const double H2 = 0.00390625;

const int J[15][6] = {
  { 2 } ,
  { 1  , 3 } ,
  { 1  , 2  , 4 } ,
  { 1  , 2  , 3  , 5 } ,
  { 1  , 2  , 3  , 4  , 6 } ,
  { 1  , 2  , 3  , 4  , 5  , 7 } ,
  { 2  , 3  , 4  , 5  , 6  , 8 } ,
  { 3  , 4  , 5  , 6  , 7  , 9 } ,
  { 4  , 5  , 6  , 7  , 8  , 10 } ,
  { 5  , 6  , 7  , 8  , 9  , 11 } ,
  { 6  , 7  , 8  , 9  , 10 , 12 } ,
  { 7  , 8  , 9  , 10 , 11 , 13 } ,
  { 8  , 9  , 10 , 11 , 12 , 14 } ,
  { 9  , 10 , 11 , 12 , 13 , 15 } ,
  { 10 , 11 , 12 , 13 , 14 }
};

const int L [15] = { 1 , 2 , 3 , 4 , 5 , 6 , 6 , 6 , 6 , 6 , 6 , 6 , 6 , 6 , 5 };

/*---------------------------------------*/
/*           penalty function II         */
/*---------------------------------------*/
/*  S1 variables : 1 --> 15              */
/*    (n=15 variables)                   */
/*---------------------------------------*/
double f_penalty ( double x[60] ) {

  double f = pow ( x[0] - 0.2 , 2 );
  int    i;

  for ( i = 2 ; i <= 15 ; i++ )
    f += pow ( SQRT_A * ( exp ( x[i-1] / 10 ) + exp ( x[i-2] / 10 )
			  - exp ( i/10.0 ) - exp ( (i-1)/10.0 ) ) , 2 );

  for ( i = 16 ; i < 30 ; i++ )
    f += pow ( SQRT_A * ( exp ( x[i-15] / 10 ) - E1 ) , 2 );

  double s = -1.0;
  for ( i = 1 ; i <= 15 ; i++ )
    s += (16-i) * pow ( x[i-1] , 2 );

  return f + pow(s,2);
}

/*---------------------------------------*/
/*             trigo. function           */
/*---------------------------------------*/
/*  S2 variables : 16 --> 45             */
/*    (n=30 variables)                   */
/*---------------------------------------*/
double f_trigo ( double x[60] ) {

  int    i;
  double s = cos(x[15]);
  for ( i = 16 ; i < 45 ; i++ )
    s += cos(x[i]);
  
  double f = 0.0;
  for ( i = 1 ; i <= 30 ; i++ )
    f += pow ( 30 - s + i * ( 1 - cos(x[i+14]) ) - sin (x[i+14]) , 2 );

  return f;
}

/*---------------------------------------*/
/*       Brown almost-linear function    */
/*---------------------------------------*/
/*  S3 variables : 46 --> 60             */
/*    (n=15 variables)                   */
/*---------------------------------------*/
double f_brown ( double x[60] ) {

  int i;

  double s = x[45];
  for ( i = 46 ; i < 60 ; i++ )
    s += x[i];

  double f = 0.0;
  for ( i = 1 ; i < 15 ; i++ )
    f += pow ( x[i+44] + s - 16 , 2 );

  s = x[45];
  for ( i = 46 ; i < 60 ; i++ )
    s *= x[i];
  
  return f + pow ( s-1 , 2 );
}

/*---------------------------------------*/
/*        Broyden banded function        */
/*---------------------------------------*/
/*  S1 variables : 1 --> 15              */
/*    (n=15 variables)                   */
/*---------------------------------------*/
double f_broyden_banden ( double x[60] ) {
  int i , j;
  double s , f = 0.0;
  for ( i = 1 ; i <= 15 ; i++ ) {
    s = 0.0;
    for ( j = 0 ; j < L[i-1] ; j++ )
      s += x[ J[i-1][j] - 1 ] * ( 1 + x[ J[i-1][j] - 1 ] );
    f += pow ( x[i-1] * ( 2 + 5 * pow ( x[i-1] , 2 ) ) + 1 - s , 2 );
  }
  return f;
}

/*---------------------------------------*/
/*        Broyden tridiagonal function   */
/*---------------------------------------*/
/*  S2 variables : 16 --> 45             */
/*    (n=30 variables)                   */
/*---------------------------------------*/
double f_broyden_tridiag ( double x[60] ) {
  double f = pow ( ( 3 - 2 * x[15] ) * x[15] - 2 * x[16] + 1 , 2 );
  for ( int i = 16 ; i <= 43 ; i++ )
    f += pow ( ( 3 - 2 * x[i] ) * x[i] - x[i-1] - 2 * x[i+1] + 1 , 2 );
  return f + pow ( (3-2*x[44]) * x[44] - x[43] + 1 , 2 );
}

/*---------------------------------------*/
/*    discrete boundary value function   */
/*---------------------------------------*/
/*  S3 variables : 46 --> 60             */
/*    (n=15 variables)                   */
/*---------------------------------------*/
double f_discrete_boundary ( double x[60] ) {
  double f = pow ( 2*x[45] - x[46] + H2 * pow( x[45] + H + 1 , 3 ) / 2.0 , 2 );
  for ( int i = 2 ; i <= 14 ; i++ )
    f += pow ( 2*x[i+44] - x[i+43] - x[i+45] + H2 * pow( x[i+44] + i * H + 1 , 3 ) / 2.0 , 2 );
  return f + pow ( 2*x[59] - x[58] + H2 * pow( x[59] + 15 * H + 1 , 3 ) / 2.0 , 2 );
}

/*---------------------------------------*/
/*               main function           */
/*---------------------------------------*/
int main ( int argc , char ** argv ) {

  if ( argc != 2 ) {
    cout << 1e20 << " " << 1e20 << endl;
    return 0;
  }
  
  // input file :
  ifstream in (argv[1]);
  double   x[60];
  int      i;

  for ( i = 0 ; i < 60 ; i++ )
    in >> x[i];
  
  if ( in.fail() ) {
    cout << 1e20 << " " << 1e20 << endl;
    return 0;
  }
   
  in.close();

//   cout << "f_penalty=" << f_penalty(x) << endl;
//   cout << "f_trigo=" << f_trigo(x) << endl;
//   cout << "f_brown=" << f_brown(x) << endl;

//   cout << "f_broyden_banden=" << f_broyden_banden(x) << endl;
//   cout << "f_broyden_tridiag=" << f_broyden_tridiag(x) << endl;
//   cout << "f_discrete_boundary=" << f_discrete_boundary(x) << endl;


  cout << /* "f(x)=" << */ f_penalty(x)+f_trigo(x)+f_brown(x) << " "
       << /* "g(x)=" << */ -f_broyden_banden(x)-f_broyden_tridiag(x)-f_discrete_boundary(x) + 250
       << endl;

  return 0;
}
