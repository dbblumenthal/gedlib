#ifndef PFR_H
#define PFR_H

#include "reaction.hpp"
#include "RungeKutta.hpp"
#include "stream.hpp"
using namespace std;

class pfr {
private:
  // terminator *term;
  // ofstream results;
  string name;
  bool OK, explode;
  int i ,j, n, m;
  double L,D,dL, U, Ta, m_in, sum, P;
  stream *F;
  double **a, *C, T, *y, *r, tmp, tmp1;
  reaction **rx;
  RungeKutta<pfr> *solver;
  
public:
  // pfr(){};
  pfr ( stream * , stream * , double ** , int , reaction ** , double , double );
  void set ( const string & n ) { name = n; }
  void set(double l, double d) {L=l; D=d;}
  bool run();
  void water();
  void  cost();
  double get_cost ( void );
  double get_water ( void );

  double f(int, double, double*);
  ~pfr();
};
#endif
