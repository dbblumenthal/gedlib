/*
This unit takes more than two input streams and merge them in one
output stream. The pressure must be specified by the user, and the
temperature of the output stream is computed.

Structure in the .process file:
mix {name} {pressure} {nb_in} {indexes of input streams} {index of output stream}

How to use:
   1- Call the constructor: mix1 = new mix(nb_in, list1_in, out);
   2- Set the operating pressure : mix1->set(P);
   3- Set the name of the unit: mix1->set(name);
   4- Solve: bool=mix1->solve();
*/
#ifndef MIX_H
#define MIX_H

#include "stream.hpp"
#include "bissection.hpp"
using namespace std;

class mix {
 private:
  int i, j;
  bool success;
  bissection<mix> *solver;
  string name;
  int nb_in;			//number of input streams
  stream **in;		//list pointers to input streams
  stream *out;		//pointer to output stream
  // double min, max;
  
 public:
  double P, T;			//pressure (given) and temperature (unknown)   
  mix(){P=0.0;}
  mix(int, stream**, stream*);		//defines the connectivities of this unit
  ~mix(){delete solver;}
  void set(double p) {P=p;}
  void set ( const string & n ) { name = n; }
  double f(double); //returns the function to the solver
  bool solve();     //finds the temperature and computes mass balance
  void write();
};
#endif
