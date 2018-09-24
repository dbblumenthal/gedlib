#include "heatx.hpp"
#include "bissection.cpp"

heatx::heatx(bool m, stream* s1, stream* s2)
{
   in = s1;
   out = s2;
   out->m=0;
   for(i=0;i<in->nb;i++)
   {
      out->chem[i]->m=in->chem[i]->m;
      out->m+=out->chem[i]->m;
   }
   solver = new bissection<heatx>();
   mode = m;
   success = true;
}

void heatx::set(double d1, double d2)
{
   if(mode==0) T=d1;
   if(mode==1) Q=d1;
   eta = d2;
}

bool heatx::solve()
{
   if(mode==0)
   {
      Q = 0.0;
      out->set(in->P, T);
      // out->write();  // WRITE TOTO
      for(i=0;i<in->nb;i++)
         Q+=out->chem[i]->dH(in->T, out->T, in->P)*out->chem[i]->n();
      if(eta>EPS)
	Qreal = Q/eta;
      else {
	Qreal=Q;
	success=false;
      }
   }
   if(mode==1)
   {
      Qreal = eta*Q;
      min = in->T;
      max = 2000;
      solver->set(this, min, max);
      success = solver->run();
      out->set(in->P, T);
      //out->write();  // WRITE TOTO
   }
   return success;
}


double heatx::f(double x)
{
   T=x;
   max = Qreal;
   for(i=0;i<in->nb;i++)
      max -= out->chem[i]->dH(in->T, T, in->P)*out->chem[i]->n();
   return max;
}

void heatx::write()
{

  cout << "WRITE FILE " << RUNTIME << name << ".unit" << " :\n\tBEGIN\n";
  cout << "\t>>         " << name;
  cout << endl << "\t>>           stream in: "  << in->name  << "  T = " << in->T <<" K";
  cout << endl << "\t>>           stream out: " << out->name << "  T = " << out->T <<" K";
  if (mode==0)
    cout<<endl<<"\t>>           Heat duty : "<<Qreal<<" kW";
  if (mode==1)
    cout<<endl<<"\t>>           Heat duty : "<<Q<<" kW";
  if (success)
    cout<<endl<<"\t>>           Heat losses "<<fabs(Qreal-Q)<<" kW  (converged normally)";
  cout << "\n\tEND\n\n";

  cost();
  power();
  water();
}



double heatx::get_cost ( void ) {
  if(mode==1) min=fabs(Q)/0.225/(eta)/fabs(out->T-in->T);
  if(mode==0) min=fabs(Qreal)/0.25/(eta)/fabs(out->T-in->T);
  if(min<10) min=10; if(min>1000) min=1000;
  max = 4.3247-0.303*log10(min)+0.1634*pow(log10(min),2);
  T=in->P;
  T = (T-1)*1.01325;
  if (fabs(T)<EPS) T=0.1; if(T>100) T=100;
  min=0.03881-0.11272*log10(T)+0.08183*pow(log10(T),2);
  min=pow(10, min);
  max = (1.63+1.66*2.5*min)*pow(10, max);
  max = max*MS_YEAR/MS_2001;
  return max;
}


void heatx::cost()
{
  cout << setprecision(5);
  cout << "WRITE FILE " << RUNTIME << name << ".cost" << " :\n\tBEGIN\n";
  cout << "\t>>" << get_cost();
  cout << "\n\tEND\n\n";
}



double heatx::get_water ( void )
{
  max = (Q<0.0) ? fabs(Q)/(4.185*0.10*(out->T-298)) : 0.0;
  return max;
}


void heatx::water()
{
  if(Q<0.0) max = fabs(Q)/(4.185*0.10*(out->T-298));
  else max = 0.0;
  cout << "WRITE FILE " << RUNTIME << name << ".water" << " :\n\tBEGIN\n";
  cout << "\t>>" << max;
  cout << "\n\tEND\n\n";

}

double heatx::get_power ( void ) {
  max = (mode) ? Q : Qreal;
  if (max>EPS)
    return max;
  return 0.0;
}

void heatx::power()
{
  if(mode==0) max = Qreal;
  if(mode==1) max = Q;

  cout << "WRITE FILE " << RUNTIME << name << ".power" << " :\n\tBEGIN\n";
  if(max>EPS) cout<< "\t>>" << max;
  else cout<< "\t>>" << 0;
  cout << "\n\tEND\n\n";
}
