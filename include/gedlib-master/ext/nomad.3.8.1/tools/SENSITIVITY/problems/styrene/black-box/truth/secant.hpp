/*
To use the secant solver to find the root of a scalar function:
   (the parametric object E must have public function such as E->f(x),
      where x is the point at which evaluate f.)
   1- construct the solver :    solver = new secant<E>();   
   2- set the solver :    solver->set(unit, x0, x1);   //unit is usually the pointer *this, and x0 and x1 are two required initial points
   3- launch the solver : bool = solver->run();   //will return true is success, false if the solver failed
*/
#ifndef SECANT_H
#define SECANT_H

#include "defines.hpp"
using namespace std;


template <class E>
class secant {
private:
  double x_last, x_now, x_next;
  double f_last, f_now, error;
  int i;
  bool OK;
  E *unit;
	  
public:
  secant();
  void set(E*, double, double);
  bool run();
  ~secant(){}
};
#endif
