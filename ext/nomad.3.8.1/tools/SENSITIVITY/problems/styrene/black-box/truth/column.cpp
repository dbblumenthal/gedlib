#include "column.hpp"
using namespace std;

column::column(stream* in, stream* out_B, stream* out_D)
{
  F = in;
  B = out_B;
  D = out_D;
  L=new stream("columnL", F->nb, F->chem);
  V=new stream("columnV", F->nb, F->chem);
  flasher = new flash(F, L, V);
  alpha_1 = new double[F->nb];
  alpha_f = new double[F->nb];
  alpha_N = new double[F->nb];
  alpha_m = new double[F->nb];
  T_d=0; T_b=0; T_f=F->T;
//   for(i=0;i<F->nb;i++)
//     if(F->chem[i]->Tc<F->T && F->chem[i]->m>EPS)
//       {
// 	logf.open(MESSAGES,ios::app);
// 	logf<<"   --> Warning <--  Presence of gas in column.\n";
// 	logf.close();
// 	i=F->nb;
//       }
}

column::~column()
{
   delete L; delete V; delete flasher;
   delete []  alpha_1;
   delete [] alpha_f;
   delete [] alpha_N;
   delete [] alpha_m;
}

void column::set(double p, int lk, double xd, int hk, double xb)
{
   P=p;
   HK=hk-1; x_B=xb;
   LK=lk-1; x_D=xd;
   //logf.open(MESSAGES, ios::app);
 //  F->write(); system("pause");
 //  if (F->chem[LK]->m<EPS) logf<<"   ==> Error <==  Flow of light key in column "<<name<<" is zero.\n";
//   if (F->chem[HK]->m<EPS) logf<<"   ==> Error <==  Flow of heavy key in column "<<name<<" is zero.\n";
   //logf.close();
}

bool column::solve()
{

   OK=true;
 //  B->thermo=F->thermo; D->thermo=F->thermo;
   //flash once the feed stream


   
   flasher->set(P,F->T);

   flasher->adiabatic();


   T_f=flasher->T;
   L->set(P, T_f); V->set(P, T_f);
   // L->write(); V->write();  TOTO
   //check if a column is needed; if not, bypass block
   if(F->chem[LK]->n()/F->n()<0.001) OK=false;
   if(F->chem[HK]->n()/F->n()<0.001) OK=false;
   if(!OK)
     return false;
//    {
//       strcpy(filename, B->name); *B=*L; B->set(filename); B->write();
//       strcpy(filename, D->name); *D=*V; D->set(filename); D->write();
//    }
   else
   {





      //apply the FUG method
      first_split();
      Nmin = Fenske();
      N=Nmin+1;
      while (fabs(N-Nmin)>0.1)
      {
         N=Nmin;
         D->set(P, T_f); T_d=D->bp;
         B->set(P, T_f); T_b=B->bp;
         set_alpha();
         distribute();
         Nmin = Fenske();
         if (Nmin<1) Nmin=1;
      }
      D->set(P, T_d);
      B->set(P, T_b);
      if(fabs(Nmin)<=MIN_PLATES || fabs(Nmin)>MAX_PLATES)  OK=false;
      else
      {
         Rmin = Underwood();
         if(Rmin>100) Rmin=100;
         if(L->chem[HK]->m+L->chem[LK]->m<EPS) Rmin=10.0;
         if (Nmin<5) Ract = 1.5*Rmin;
         if (5<Nmin && Nmin<15) Ract = 1.3*Rmin;
         if (15<=Nmin) Ract = 1.1*Rmin;
         N = Gilliland();
         feed = Kirkbride();
         condense();
         reboil();
      }
      // B->write();  TOTO
      // D->write(); TOTO
   }
   return OK;
}

void column::first_split()
{
   B->purge(); D->purge();
   set_alpha();
   //Check if LK is really lighter than HK
   if (alpha_m[LK]<1)
   {
//       logf.open(MESSAGES,ios::app);
//       logf<<"   --> Warning <--  Swapping keys in column "<<name<<".\n";
//       logf.close();
      feed=LK; LK=HK; HK=feed; set_alpha();
   }
   for(i=0;i<F->nb;i++)
   {
      if (i!=LK && i!=HK && F->chem[i]->m>EPS)
      {
         if(alpha_f[i] > alpha_f[LK]) //volatile
         {
            D->chem[i]->m = (alpha_f[i]-alpha_f[LK])/alpha_f[i]*F->chem[i]->m;
            D->m += D->chem[i]->m;
            B->chem[i]->m = F->chem[i]->m-D->chem[i]->m;
            B->m+=B->chem[i]->m;
         }
         if(alpha_f[i] < 1)           //not volatile
         {
            B->chem[i]->m = (alpha_f[HK]-alpha_f[i])/alpha_f[i]*F->chem[i]->m;
            B->m += B->chem[i]->m;
            D->chem[i]->m = F->chem[i]->m-B->chem[i]->m;
            D->m+=D->chem[i]->m;
         }
         if(1 <= alpha_f[i] && alpha_f[i]<=alpha_f[LK]) //ambiguous volatility
         {
            D->chem[i]->m = (alpha_f[i]-1)/(alpha_f[LK]-1)*F->chem[i]->m;
            B->chem[i]->m = F->chem[i]->m-D->chem[i]->m;
            D->m+=D->chem[i]->m;
            B->m+=B->chem[i]->m;
         }
      }
   }
   D->chem[HK]->m = D->n()*x_D/(1-x_D)*D->chem[HK]->M/1000.0;
   if(D->chem[HK]->m<EPS) D->chem[HK]->m=0.01*F->chem[HK]->m;
   B->chem[LK]->m = B->n()*x_B/(1-x_B)*B->chem[LK]->M/1000.0;
   if(B->chem[LK]->m<EPS) B->chem[LK]->m=0.01*F->chem[LK]->m;
   B->chem[HK]->m = F->chem[HK]->m - D->chem[HK]->m;
   D->chem[LK]->m = F->chem[LK]->m - B->chem[LK]->m;
   D->m += (D->chem[LK]->m + D->chem[HK]->m);
   B->m += (B->chem[LK]->m + B->chem[HK]->m);
}
void column::distribute()
{
   D->m=0; B->m=0;
   for(i=0;i<F->nb;i++)
   {
      if (i!=LK && i!=HK  && F->chem[i]->m>EPS)
      {
         if(alpha_m[i] > 1) //volatile and ambiguous
         {
            B->chem[i]->m = F->chem[i]->m/(1+D->chem[HK]->n()/B->chem[HK]->n()*pow(alpha_m[i], Nmin));
            D->chem[i]->m = F->chem[i]->m - B->chem[i]->m;
         }
         if(alpha_m[i] <= 1)           //not volatile
         {
            D->chem[i]->m = F->chem[i]->m*(D->chem[HK]->n()/B->chem[HK]->n()*pow(alpha_m[i], Nmin))/(1+D->chem[HK]->n()/B->chem[HK]->n()*pow(alpha_m[i], Nmin));
            B->chem[i]->m = F->chem[i]->m - D->chem[i]->m;
         }
         D->m+=D->chem[i]->m;
         B->m+=B->chem[i]->m;
      }
   }
   D->m += (D->chem[LK]->m + D->chem[HK]->m);
   B->m += (B->chem[LK]->m + B->chem[HK]->m);
}

void column::set_alpha()
{
   for(i=0;i<F->nb; i++)
   {
      if (T_b>EPS && F->chem[i]->m>EPS) alpha_1[i] = F->chem[i]->Psat(T_b)/F->chem[HK]->Psat(T_b);
      else alpha_1[i]=0;
      if (T_d>EPS&& F->chem[i]->m>EPS) alpha_N[i] = F->chem[i]->Psat(T_d)/F->chem[HK]->Psat(T_d);
      else alpha_N[i]=0;
      if (T_f>EPS&& F->chem[i]->m>EPS) alpha_f[i] = F->chem[i]->Psat(T_f)/F->chem[HK]->Psat(T_f);
      else alpha_f[i]=0;
      alpha_m[i] = pow(alpha_1[i]*alpha_f[i]*alpha_N[i], 1.0/3.0);
   }
   for(i=0;i<F->nb;i++) if(alpha_m[i]<EPS&& F->chem[i]->m>EPS) alpha_m[i] = alpha_f[i];
}

void column::reboil()
{
   Q_reboil = 0.0;
   for (i=0;i<F->nb;i++) if(F->chem[i]->m>EPS)
   {   
      Q_reboil += B->chem[i]->Cp(false)*(T_b-T_f)*B->chem[i]->n()/1000.0; //energy to go from input to bottoms T
      Q_reboil += D->chem[i]->Cp(false)*(T_f-T_d)*D->chem[i]->n()/1000.0; //energy to go from input to tops T
   }
   Q_reboil += Q_condens;
}

void column::condense()
{
   Q_condens = 0.0;
   for (i=0;i<F->nb;i++) if(F->chem[i]->m>EPS)
   {   
      Q_condens += D->chem[i]->Hvap(T_d)*(1+Ract)*D->chem[i]->n();
   }
}

void column::write()
{
  cout << setprecision(11);

  cout << "WRITE FILE " << RUNTIME << name << ".unit" << " :\n\tBEGIN\n";
  cout <<"\t>>         "<<name;
  cout <<endl<<"\t>>           stream in: "<<F->name;
  cout <<endl<<"\t>>           streams out: "<<B->name<<" (bot.)  "<<D->name<<" (top.)";
  cout <<endl<<"\t>>           P = "<<P<<" atm,  T(0) = "<<T_b<<",  T("<<feed<<") = "<<T_f<<",  T("<<int(N)<<") = "<<T_d<<"  K";
  cout <<endl<<"\t>>           Number of stages: "<<int(N)<<" (feeding at stage "<<feed<<")";
  cout <<endl<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(5)<<"\t>>           LK purity = "<<D->chem[LK]->n()/D->n()<<"  HK purity = "<<B->chem[HK]->n()/B->n();
  cout <<endl<<"\t>>           Reboiler duty: "<<Q_reboil<<" kW   Condenser duty: "<<(-1)*Q_condens<<" kW";
  cout << "\n\tEND\n\n";
  cost(); water(); power();
}

double column::get_cost()
{
  //cost of vessel
  vol=(0.45*N)*(pow(300*D->v, 1.5)/2.4/sqrt(B->v))*sqrt(D->m/B->m);
  if(vol<0.3) vol=0.3; if(vol>520)vol=520;
  money = 3.4974+0.4485*log10(vol)+0.1074*pow(log10(vol),2);
  money = pow(10, money);
  P= (P-1)*101.325/100;
  diam=sqrt(4.0*vol/pi/N/0.45);
  vol=(P+1)*diam/(317.46*(850-0.6*(P+1)))+0.0315;
  money *=(2.25+ 1.82*vol*2.2);
  //cost of trays
  vol = (pow(300*D->v, 1.5)/2.4/sqrt(B->v))*sqrt(D->m/B->m);
  diam = 2.9949+0.4465*log10(vol)+0.3961*pow(log10(vol),2);
  money+=1.5*pow(10, diam);
  //cost of reboiler     U=5250W/m2.K
  vol=fabs(Q_reboil)/0.85/5.25/15.0;
  if(vol<10) vol=10; if(vol>100) vol=100;
  vol = 4.4646-0.5277*log10(vol)+0.3955*pow(log10(vol),2);
  money += (1.63+1.66*2.5)*pow(10, vol);
  //cost of condenser    U=1850W/m2.K
  vol=fabs(Q_condens)/0.85/1.85/(0.5*(T_d-298));
  if(vol<1) vol=1; if(vol>100) vol=100;
  vol = 3.9912+0.0668*log10(vol)+0.243*pow(log10(vol),2);
  money += (1.74+1.55*2.5)*pow(10, vol);
  money = money*MS_YEAR/MS_2001;
  return money;
}


void column::cost()
{
  cout << "WRITE FILE " << RUNTIME << name << ".cost" << " :\n\tBEGIN\n";
  cout << "\t>>" << get_cost();
  cout << "\n\tEND\n\n";
}
void column::power()
{
  cout << "WRITE FILE " << RUNTIME << name << ".power" << " :\n\tBEGIN\n";
  money =(Q_reboil/0.85-Q_condens);
  cout << "\t>>" << money;
  cout << "\n\tEND\n\n";
}
void column::water()
{
  cout << "WRITE FILE " << RUNTIME << name << ".water" << " :\n\tBEGIN\n";
  money = (fabs(Q_condens)/(4.185*0.85*0.25*fabs(T_d-298)));
  cout << "\t>>" << money;
  cout << "\n\tEND\n\n";
}
