#include "mix.hpp"
#include "bissection.cpp"
using namespace std;

mix::mix ( int n , stream ** s1 , stream * s2 ) {
  nb_in=n;
  in = s1;
  out= s2;
  success = true;
  solver  = new bissection<mix>();
}

bool mix::solve()
{
   //Mass balance
   out->m=0;
   out->purge();
   for (j=0; j<out->nb;j++)
      for (i=0;i<nb_in;i++) if(in[i]->chem[j]->m>EPS)
      {
          out->chem[j]->m+=in[i]->chem[j]->m;
          out->m+=in[i]->chem[j]->m;
      }
   //Find the temperature
   double max=0.0; double min=1e6;
   for(i=0;i<nb_in;i++)
      for(j=0;j<out->nb;j++)
      {
         if(in[i]->chem[j]->T>max && in[i]->chem[j]->m>EPS) max=in[i]->chem[j]->T;
         if(in[i]->chem[j]->T<min && in[i]->chem[j]->m>EPS) min=in[i]->chem[j]->T;
      }
    if(fabs(max-min)<EPS) T=max;
    else
    {
      solver->set(this, min, max);
      success = solver->run();
    }
   out->set(P,T);
   //    if (success==false)
   //    {
   //       log.open(MESSAGES, ios::app);
   //       log<<"   --> Warning <--  Solver of "<<name<<" did not converge.\n";
   //       log.close();
   //    }
   //    min = 0;
   //    for(i=0;i<nb_in;i++)
   //      min+=in[i]->m;
   //    if(fabs(min-out->m)>sqrt(EPS))
   //    {
   //       log.open(MESSAGES, ios::app);
   //       log<<"   --> Warning <--  Block "<<name<<" is not in mass balance ("<<fabs(min-out->m)/min<<").\n";
   //       log.close();
   //    }

   // out->write();  // WRITE TOTO
   return success;
}

double mix::f(double x)
{
   T=x;
   double energy=0.0;   //in kW
   for (j=0; j<out->nb;j++)
      for (i=0;i<nb_in;i++)
          energy += in[i]->chem[j]->dH(in[i]->T, T, P)*in[i]->chem[j]->n()/1000;
   return energy;
}

void mix::write() {
  cout << "WRITE FILE " << RUNTIME << name << ".unit" << " :\n\tBEGIN\n";
  cout << "\t>>         " << name
       << endl << "\t>>           streams in: ";
  for ( int i = 0 ; i < nb_in ; i++ )
    cout << in[i]->name << " ";
  cout << endl << "\t>>           stream out: " << out->name;
  cout <<endl << "\t>>           P = " << P << " atm,  T = " << T;
  if (success)
    cout << " K (converged normally)";
  cout << "\n\tEND\n\n";
}
