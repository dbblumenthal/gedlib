/*
To use the Runge-Kutta solver for systems of first order differential equations,
   (the parametric object E must have public function such as E->f(i,x,y),
      where i is the index of the function to evaluate, x is the time and
	  y is a point such as y(x), returns values of differential equation i)
	  
   1- construct the solver :    solver = new RungeKutta<objet>(int);   //the integer is the dimension of x
   2- set the solver :    solver->set(unit, y0, x0, xn);   //unit is usually the pointer *this, y0 are the initial conditions, and x0 and xn is the time interval
   3- launch the solver : bool = solver->run();   //will return true is success, false if the solver failed????????????????????????????
   4- delete the solver : delete solver;
(ref  :Fortin)
*/
#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "defines.hpp"
using namespace std;

template <class E>
class RungeKutta
{
private:
  double *k1, *k2, *k3, *k4, *y_tmp, *y;
  // double k1[MAX_DIM], k2[MAX_DIM], k3[MAX_DIM], k4[MAX_DIM], y_tmp[MAX_DIM], y[MAX_DIM];
  double h, x0, xn, x;
  int i, j, m;
  bool success;
  E *unit;
	  
public:
  RungeKutta(int);
  ~RungeKutta();
  void set(E*, double*, double, double);
  double dx(){return h;}
  bool run();
};
#endif
