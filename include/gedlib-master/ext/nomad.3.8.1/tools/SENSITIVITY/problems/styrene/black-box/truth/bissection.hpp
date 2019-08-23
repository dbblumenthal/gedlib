/*
To use the bissection solver to find the root of a scalar function
when the secant solver fails.
   (the parametric object E must have public function such as E->f(x),
      where x is the point at which evaluate f.)
   1- construct the solver :    solver = new bissect<E>();   
   2- set the solver :    solver->set(unit, x0, x1);   //unit is usually the pointer *this, and x0 and x1 are the interval's bounds
   3- launch the solver : bool = solver->run();   //will return true is success, false if the solver failed
*/
#ifndef BISSECTION_H
#define BISSECTION_H

#include "defines.hpp"
using namespace std;

template <class E>
class bissection {
private:
  double x1, xm, x2;
  double f1, fm, f2;
  int i;
  bool OK;
  E *unit;
	  
public:
  bissection(){x1=xm=x2=f1=fm=f2=0; OK=false;}
  void set(E* tmp, double xx1, double xx2) {unit=tmp; x1=xx1; x2=xx2;}
  bool run();
  ~bissection(){}
};
#endif
