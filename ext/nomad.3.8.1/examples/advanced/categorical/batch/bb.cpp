/*-------------------------------------------------------------------*/
/*            Example of a problem with categorical variables        */
/*-------------------------------------------------------------------*/
/*                                                                   */
/*  . portfolio problem with 2 assets                                */
/*                                                                   */
/*  . NOMAD is used in batch mode                                    */
/*                                                                   */
/*  . the number of variables is fixed to 4 (2 assets)               */
/*                                                                   */
/*  . variables are of the form (t0 v0 t1 v1) where ti is the type   */
/*    of an asset and vi is the money invested in this asset         */
/*                                                                   */
/*  . categorical variables are t0 and t1                            */
/*                                                                   */
/*  . with a $10,000 budget, the problem consists in minimizing      */
/*    some measure of the risk and of the revenue                    */
/*                                                                   */
/*-------------------------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

int main ( int argc , char ** argv ) {


  double g1 = 1e20 , g2 = 1e20 , f = 1e20;

  if ( argc == 2 ) {
    ifstream in ( argv[1] );
    if ( !in.fail() ) {

      double t0 , v0 , t1 , v1;

      in >> t0 >> v0 >> t1 >> v1;

      if ( !in.fail() ) {

	double v[3] , vmin = 10000;
	v[0] = v[1] = v[2] = 0.0;

	v[static_cast<int>(t0)] = v0;
	v[static_cast<int>(t1)] = v1;

	if ( v0 < vmin )
	  vmin = v0;
	if ( v1 < vmin )
	  vmin = v1;
      
	// constraints (budget and each asset is considered with at least 1$):
	double vt = v[0] + v[1] + v[2]; 
	double h  = vt - 10000;
	g1 = h;
	g2 = 1-vmin;

	if ( h <= 0 && vmin >= 1 ) {

	  // compute the risk and revenue:
	  double vt2  = vt*vt;
	  double rev  = v[0] * 0.0891 + v[1] * 0.2137 + v[2] * 0.2346;
	  double risk = 0.01 * (v[0]/vt)*(v[0]/vt) +
	                0.05 * (v[1]/vt)*(v[1]/vt) +
	                0.09 * (v[2]/vt)*(v[2]/vt) +
	                0.02 * (v[0]*v[1]/vt2)     +
	                0.02 * (v[0]*v[2]/vt2)     +
	                0.10 * (v[1]*v[2]/vt2);

	  // the objective is taken as a scaled distance
	  // between (risk,revenue) and (risk_best,rev_best):
	  double a = ( risk - 0.01 ) * 100 / 0.08;
	  double b = ( rev  - 891  ) * 100 / 1455;

	  f = sqrt ( a*a + (100-b)*(100-b) );
	}
	else
	  f = 145;
      }
      in.close();
    }
  }
  cout.precision(18);
  cout <<  g1 << " " << g2 << " " << f << endl;

  return 0;
}
