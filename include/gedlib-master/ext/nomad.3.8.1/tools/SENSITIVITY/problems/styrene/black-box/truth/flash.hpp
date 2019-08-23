/*
This unit takes operating P and T, and indexes of streams. It applies
the Rachford-Rice procedure to solve the isothermal flash problem
(ref : Seader & Henley).

Structure in the .process file:
flash {name} {pressure} {temperature} {index of input stream} {index of output liquid and output vapor}

How to use:
   1- Call the constructor : flash1 = new flash(in, out_L, out_V);			//in is the feed, out_L is the liquid output and out_V is the vapor output 
   2- Set P and T: flash1->set(P,T);
   3- Set the name: flash1->set(name);
   4a- Perform an isothermal flash : flash1->solve();
   4b- Perform an adiabatic flash: flash1->adiabatic();
*/
#ifndef FLASH_H
#define FLASH_H

#include "stream.hpp"
#include "bissection.hpp"
using namespace std;

class flash {
private:
  bool success;
  bissection<flash> *solver;
  // ofstream log, results;
//   char name[31], filename[41];  //name of the unit

  string name;
  
  int i, task;          //task=0: isothermal flash;   task=1:adiabatic flash
  stream *F, *Fcopy;		   //pointer to the input stream
  stream *L, *V;	//pointers to liquid and vapor output streams
  double f_x, x, *K, *z;			//pressure (given) and temperature (given)
  double Q, Tin, step, vol;			//required power, in kW
  
public:
  flash(){P=0.0; T=0.0;}
  flash(stream*, stream*, stream*);		//defines the connectivities of this unit
  ~flash(){delete Fcopy; delete [] K; delete []  z; delete solver;}
  double P ,T;
  void set(double, double);
  void set( const string & n ) { name = n; }
  bool solve();								//applies the Rachford-Rice procedure
  bool adiabatic();						//adiabatic isobaric flash
  double f(double);						//returns the function to the solver
  void write();
  double get_water ( void );
  double get_cost ( void );

  double get_power ( void ) const { return Q; }
  void cost(), water(), power();
};
#endif
