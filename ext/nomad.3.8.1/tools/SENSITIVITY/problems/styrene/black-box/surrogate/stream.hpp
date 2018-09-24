#ifndef STREAM_H
#define STREAM_H
#include "chemical.hpp"
#include <iomanip>
using namespace std;

class stream {

private:

  stream ( void ) {}
  stream ( const stream & s ) {}

public:
  void check_error();

  double P, T, m, v;  // Pressure in atm, temperature in K,
                      // total mass flow in kg/s, total volume flow in m3/s

  int i, j, k, error, warning;

  string name;

  int nb;		// number of chemicals to store;

  chemical ** chem;     // list of stored chemicals

  thermolib *thermo;

  double tmp, tmp1,tmp2, *tab1, *tab2, *tab3, *tab4;
  double step;
  void bubble(), dew();
  
public:





  double dp, bp;
  double rho() {if(v!=0) return m/v; else exit(0);} // Apparent density in kg/m3
  double n();		// total mole flow (mol/s)
  double quality();	// returns the vapor fraction of stream (1=all vapor, 0 = all liquid)
  double K(int);        // returns the vapor-liquid equilibrium constant of compound i
   
  // constructor :
  stream ( const string _name , int n , chemical ** list ) : P      (0                 ) ,
							     T      (0                 ) ,
							     m      (0                 ) ,
							     v      (0                 ) ,
							     error  (0                 ) ,
							     warning(0                 ) ,
							     name   (_name             ) ,
							     nb     (n                 ) ,
							     chem   (new chemical *[nb]) ,
							     thermo (new thermolib (nb)) ,
							     tab1   (new double    [nb]) ,
							     tab2   (new double    [nb]) ,
							     tab3   (new double    [nb]) ,
							     tab4   (new double    [nb])   {
    for  ( i = 0 ; i < nb ; i++ )
      chem[i] = new chemical(*list[i]);
  }


  // affectation :
  stream & operator = ( const stream & s );

  void set_thermo ( thermolib * t ) { (*thermo) = *t; }


  void set(double, double); // Sets P and T
  void set(double*);        // Sets mass flows
  void set( const string & n ) { name = n;}
  void set ( int _nb , chemical ** _chem);


  // destructor :
  ~stream ( void );

  void write();
  void purge();
};


#endif
