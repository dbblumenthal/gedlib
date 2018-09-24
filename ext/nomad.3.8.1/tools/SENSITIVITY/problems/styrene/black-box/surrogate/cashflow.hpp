#ifndef CASHFLOW_H
#define CASHFLOW_H
#include "defines.hpp"
// #include <iomanip>
using namespace std;

class cashflow
{
   private:
  // char name[31], filename[41];
  // ifstream in;
  // ofstream out;
  double Itot, Ctot, Rtot;
  void set_Inv(), set_Amort(), set_C_R();
  double temp;
  int i, j, counter;
  bool OK;
  double yield_tab[15];

  double yield ( int k ) const { return yield_tab[(k==15) ? 14 : k%15]; }
  
public:
  double *Inv, *Coper, *Amort, *Rev, *Flow, *Flowact;
  double  i_rate, a_rate;
  int N;
  cashflow(int);
  ~cashflow();
  // void set(char n[31]) {strcpy(name, n);}
  void set_rates(double d1, double d2){i_rate=d1; a_rate=d2;}
  void set_basics(double d1, double d2, double d3){Itot=d1; Ctot=d2; Rtot=d3;}
  bool run();
};
#endif
