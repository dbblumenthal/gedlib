#include "RungeKutta.hpp"

#include <iomanip>

using namespace std;

template <class E>
RungeKutta<E>::RungeKutta(int dim)
{
   m = dim;
   k1 = new double[m];
   k2 = new double[m];
   k3 = new double[m];
   k4 = new double[m];
   y = new double[m];
   y_tmp = new double[m];
}

template <class E>
RungeKutta<E>::~RungeKutta()
{
   delete [] k1;
   delete [] k2;
   delete [] k3;
   delete [] k4;
   delete [] y;
   delete [] y_tmp;
}

template <class E>
void RungeKutta<E>::set ( E * tmp , double * y0 , double beg , double end )
{
  unit=tmp;
  x0=beg; xn=end;
  x=x0;
  h=double(xn-x0)/double(N_INTER);
  for (i=0;i<m;i++) {y[i]=y0[i];}
  success=true;
}

template <class E>
bool RungeKutta<E>::run() {
  for(j=0;j<MAX_ITER_RK;j++) {
    //Avoid going out of x interval
    if (x+h >xn) {
      h = xn-x;
      j = MAX_ITER_RK;
    }

    //Compute k1, k2, k3, k4
    for(i=0;i<m;i++)
      k1[i] = h*unit->f(i, x, y);
    for(i=0;i<m;i++)
      y_tmp[i] = y[i]+k1[i]/2.0;
    for(i=0;i<m;i++)
      k2[i] = h*unit->f(i, x+h/2.0, y_tmp);
    for(i=0;i<m;i++)
      y_tmp[i] = y[i]+k2[i]/2.0;
    for(i=0;i<m;i++)
      k3[i] = h*unit->f(i, x+h/2.0, y_tmp);
    for(i=0;i<m;i++)
      y_tmp[i] = y[i]+k3[i];
    for ( i = 0 ; i < m ; i++ )
      k4[i] = h*unit->f ( i , x+h , y_tmp );
    //Compute the new y
    for(i=0;i<m;i++)
      y[i]+=(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
    x += h;
  }

  if ( x < xn-EPS ) {// MODIF SEB (le EPS)
    success=false;
    
    // cout.setf(ios::fixed);
    // cout << setprecision(12);
    // cout << "x=" << x << " < xn=" << xn << " diff=" << xn-x << endl;
  }
  
  return success;
}
