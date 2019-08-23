#include "flash.hpp"
#include "bissection.cpp"
using namespace std;

flash::flash ( stream * in , stream * out_L , stream * out_V ) {
  F = in;
  Fcopy = new stream("Fcopy", F->nb, F->chem);
  Tin = F->T;
  z = new double[F->nb];
  for ( i = 0 ; i < F->nb ; i++ )
    z[i] = F->chem[i]->n()/F->n();
  
  L = out_L;
  V = out_V;
  success = true;
  K = new double[F->nb];
  task=0;
  solver = new bissection<flash>();
}

void flash::set(double p, double t)
{
   P=p;
   T=t;
   for (i=0;i<F->nb;i++)
   {
     if(F->chem[i]->Tc<T) K[i] = F->chem[i]->Psat(T)/P;
     else K[i]=1;
   }
   F->set(P,T);
}

bool flash::solve()
{
   L->purge();
   V->purge();
   f_x=F->quality();

   if( 0.0 < f_x && f_x < 1.0)
   {

     // TOTO
     for ( i = 0 ; i < F->nb ; i++ ) {
       if ( F->chem[i]->Tc < T ) {
	 F->m -= F->chem[i]->m;
	 F->chem[i]->m = 0;
       }
     }
     
     for ( i = 0 ; i < F->nb ; i++ )
       z[i] = F->chem[i]->n()/F->n();

     solver->set(this, 0.0, 1.0);

     success=solver->run();

      if (!success)
      {
//          if(task==0){
//          log.open(MESSAGES, ios::app);
//          log<<"   --> Warning <--  Solver of FLASH "<<name<<" did not converge.\n";
//          log.close();}
         for (i=0;i<F->nb;i++)
         {
             if (T<F->chem[i]->Tc && T>F->chem[i]->Tboil(P)) {V->chem[i]->m=F->chem[i]->m; V->m+=V->chem[i]->m;}
             if (T<F->chem[i]->Tc && T<=F->chem[i]->Tboil(P)) {L->chem[i]->m=F->chem[i]->m; L->m+=L->chem[i]->m;}
         }
      }
      else
      {

         V->m = x*F->n();
         L->m = F->n() - V->m;
         //Distribute liquid components
         for (i=0;i<L->nb;i++)
         {
            L->chem[i]->m = (L->m*z[i])/(1+x*(K[i]-1))*L->chem[i]->M/1000.0;
            L->chem[i]->state=0;
         }


         L->m=0.0; for(i=0;i<L->nb;i++) L->m+=L->chem[i]->m;
         //Distribute vapor components
         for (i=0;i<V->nb;i++)
         {
            V->chem[i]->m = V->m*L->chem[i]->n()*K[i]/L->n()*V->chem[i]->M/1000.0;
            V->chem[i]->state=1;
         }
         V->m=0.0; for(i=0;i<V->nb;i++) V->m+=V->chem[i]->m;
      }
      for(i=0;i<F->nb;i++)
      if(F->chem[i]->Tc<T){V->chem[i]->m=Fcopy->chem[i]->m; V->m+=Fcopy->chem[i]->m;}

   }
   else
   {
 /*     if(task==0)
      {
         log.open("runtime\\messages.r", ios::app);
         if (T<F->dp) log<<"   --> Warning <--  Mixture in "<<name<<" can't be flashed (bp="<<F->bp<<" dp="<<F->dp<<").\n";
         log.close();
      }     */
      for (i=0;i<F->nb;i++)
      {
         if (F->chem[i]->Tc<T || f_x>=1) {V->chem[i]->m=Fcopy->chem[i]->m; V->m+=V->chem[i]->m; }
         else {L->chem[i]->m=Fcopy->chem[i]->m; L->m+=L->chem[i]->m;}

      }
      success = true;
   }
   L->set(P,T);
   V->set(P,T);
   Q = 0.0;
   for (i=0;i<F->nb;i++)
   {
      Q += L->chem[i]->dH(Tin, T, P)*L->chem[i]->n();
      Q += V->chem[i]->dH(Tin, T, P)*V->chem[i]->n();
   }
   F->m=0;
   for(i=0;i<Fcopy->nb;i++) {F->chem[i]->m = Fcopy->chem[i]->m; F->m+=F->chem[i]->m;}
   F->set(F->P,Tin);
//    if(fabs(V->m+L->m-F->m)>sqrt(EPS))
//    {
//       log.open(MESSAGES, ios::app);
//       log<<"   --> Warning <--  Block "<<name<<" is not in mass balance ("<<fabs(V->m+L->m-F->m)/F->m<<").\n";
//       log.close();
//    }
//   V->write();//  TOTO
 //  L->write(); // TOTO

   return success;
}

double flash::f(double psy)
{
   x=psy;
   f_x=0.0;
   for(i=0;i<F->nb;i++) f_x += (z[i]*(1-K[i]))/(1+psy*(K[i]-1));
   return f_x;
}

bool flash::adiabatic()
{
   task=1;
   F->set(P,T); T=F->dp;
   step=-5;
   Q=1;
   

   while (fabs(step)>0.01 && fabs(Q)>0.1)
   {
      T+=step;
      F->set(P,T);

      for (i=0;i<F->nb;i++)
	K[i] = F->chem[i]->Psat(T)/P;

      success=solve();


      if (Q<0 && step<0) step*=-0.5;
      if (Q>0 && step>0) step*=-0.5;
   }
   if (fabs(Q)<0.1) return true;
   else return false;
}

void flash::write() {
  cout << "WRITE FILE " << RUNTIME << name << ".unit" << " :\n\tBEGIN\n";
  cout << "\t>>         "<<name;
  cout << endl << "\t>>           stream in : "<<F->name;
  cout <<endl<<"\t>>           streams out : "<<L->name<<" (liq.)  "<<V->name<<" (vap.)";
  cout <<endl<<"\t>>           P = "<<P<<" atm,  T = "<<T<<" K";
  cout <<endl<<"\t>>           Heat duty = "<<Q;
  if (success==true) cout <<" kW (converge normally)";
  cout << "\n\tEND\n\n";
  cost();
  power();
  water();
}


double flash::get_cost ( void ) {
  vol=15.0*(L->v+V->v);
  if(vol<0.3) vol=0.3; if(vol>520)vol=520;
  step = 3.4974+0.4485*log10(vol)+0.1074*pow(log10(vol),2);
  step = pow(10, step);
  P= (P-1)*101.325/100;
  f_x=pow(2.0*vol/pi, 1.0/3.0);
  vol=(P+1)*f_x/(317.46*(850-0.6*(P+1)))+0.0315;
  step *=(2.25+ 1.82*vol*2.2);
  step = step*MS_YEAR/MS_2001;
  return step;
}

void flash::cost() {
  cout << "WRITE FILE " << RUNTIME << name << ".cost" << " :\n\tBEGIN\n";
  cout << "\t>>" << get_cost();
  cout << "\n\tEND\n\n";
}
void flash::power() {
  cout << "WRITE FILE " << RUNTIME << name << ".power" << " :\n\tBEGIN\n";
  cout << "\t>>" << Q;
  cout << "\n\tEND\n\n";
}

double flash::get_water ( void ) {
  step = (Q<0.0) ? fabs(Q)/(4.185*0.10*(Tin-298)) : 0.0;
  return step;
}

void flash::water() {
  cout << "WRITE FILE " << RUNTIME << name << ".water" << " :\n\tBEGIN\n";
  if(Q<0.0)
    step= (fabs(Q)/(4.185*0.10*(Tin-298)));
  else
    step= 0.0;
  cout << "\t>>" << step;
  cout << "\n\tEND\n\n";
}
