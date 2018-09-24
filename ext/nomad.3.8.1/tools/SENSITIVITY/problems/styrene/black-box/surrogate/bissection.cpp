#include "bissection.hpp"
using namespace std;

template <class E>
bool bissection<E>::run()
{
   for (i=1; i<MAX_ITER_BISSECTION; i++)
   {
      xm=(x1+x2)/2;
      // if(DEBUG) cout<<endl<<x1<<"  "<<xm<<"  "<<x2;
      if (fabs(x1-x2)/fabs(xm) < TOL_BISSECTION)
      {
         i=MAX_ITER_BISSECTION;
         OK=true;
      }
      else
      {
         f1 = unit->f(x1);
         fm = unit->f(xm);
         f2 = unit->f(x2);
         if (f1*fm < 0.0) x2 = xm;
         if (fm*f2 < 0.0) x1 = xm;
      }
   }
   // if (DEBUG) system("pause");
   return OK;
}
