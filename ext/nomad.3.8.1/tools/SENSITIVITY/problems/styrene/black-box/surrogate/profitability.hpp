#ifndef PROFITABILITY_H
#define PROFITABILITY_H
#include "cashflow.hpp"
#include "secant.hpp"

#include <iomanip>

using namespace std;

class profitability
{
private:
  cashflow *C;
//   ofstream out;
//   char name[41]; 
  bool OK;
  double ROI(), RR(), DFR();
  double PT(), AEC(), NPV();
  int i;
  double den, num, sum;
  secant<profitability> *solver;
  
public:
  profitability(cashflow* c){C=c;}
  ~profitability(){delete solver;};
  // void set(char n[31]) {strcpy(name, n); strcat(name, ".econo");}
  bool run ( double * y );
  double f(double);
};
#endif
