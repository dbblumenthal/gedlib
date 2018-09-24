/*
This unit simulates a chemical reactor. Actually, only the pfr can be used.
(ref : Fogler).

Structure in the .process file:
reactor {name} {pfr or cstr} {index of input stream} {index of output stream} {length in m} {diameter in m} {nb_react, list of reactions occuring} {U in kW/m2.K}{Ta in K}

How to use:
   1- Call the constructor : react = new reactor<pfr or cstr>(in, out);		
   2- Set dimensions and reactions : react->set(length, diameter, nb_react, list_react);	//list_react is the list of reactions names
   3- Set cooling parameters : react->set(U, Ta);
   4- Set the name : react->set(name);
   5- Run the model: react->solve();
*/
#ifndef REACTOR_H
#define REACTOR_H

#include "pfr.hpp"
using namespace std;

template<class TYPE>
class reactor {
private:
  // ofstream log;
  bool success;
  string name;
  int i ,j, m, n;
  double V, L, D, U, Ta;
  stream *in, *out;
  TYPE *model;
  reaction ** rx;
  double ** table;

public:
  // reactor(){};
  reactor(stream*, stream*);
  void set( const string & n) { name = n; }
  void set(double, double, int, const string * );
  void set(double u, double ta) {U=u;Ta=ta;}
  bool solve();
  void write();

  double get_cost  ( void ) const { return model->get_cost() ; }
  double get_water ( void ) const { return model->get_water(); }

  ~reactor();
};
#endif
