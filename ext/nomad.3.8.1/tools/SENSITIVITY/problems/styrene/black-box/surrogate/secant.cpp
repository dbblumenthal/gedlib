#include "secant.hpp"
using namespace std;

template <class E>
secant<E>::secant()
{
   x_last=0;
   x_now=0;
   x_next=0;
   f_last=0;
   f_now=0;
   OK=false;
}

template <class E>
void secant<E>::set(E* tmp, double x1, double x2)
{
   unit=tmp;
   x_last=x1;
   x_now=x2;
   OK=false;
}

template <class E>
bool secant<E>::run()
{
  // if(DEBUG) cout<<endl<<"begin solve secant";
   f_last = unit->f(x_last);
   for (i=1; i<MAX_ITER_SECANT; i++)
   {
      f_now = unit->f(x_now);
     // if(DEBUG) cout<<endl<<" x = "<<x_now<<"    f(x) = "<<f_now;
      x_next = x_now - (f_now*(x_now-x_last)/(f_now-f_last));
      if (fabs((x_next-x_now)/x_now)<=TOL_SECANT)
      {
         i=MAX_ITER_SECANT;
         OK=true;
      }
      else
      {
         x_last=x_now;
         f_last=f_now;
         x_now=x_next;
      }
   }
   return OK;
}
