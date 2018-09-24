#include "stream.hpp"
using namespace std;

// destructor :
stream::~stream ( void ) {
  delete thermo;
  delete [] tab1;
  delete [] tab2;
  delete [] tab3;
  delete [] tab4;
  for ( i = 0 ; i < nb ; i++ )
    delete chem[i];
  delete [] chem;
}

void stream::check_error()
{
  
   if (error>MAX_ERROR) {
     cout << "ERROR 4\n\n";
     exit(0);
   }
   if (warning>MAX_WARNING) {
     cout << "ERROR 5\n\n";
     exit(0);
   }
}

void stream::set ( int _nb , chemical ** _chem) {
  int i;
  for ( i = 0 ; i < nb ; i++ )
    delete chem[i];
  delete [] chem;
  nb   = _nb;
  chem = new chemical * [nb];
  for  ( i = 0 ; i < nb ; i++ )
    chem[i] = new chemical(*_chem[i]);
}

void stream::set(double pres, double temp) //set P, T, find state and resulting v
{
   if(thermo->get_dim()!=nb)
   {
     delete thermo;
     thermo = new thermolib(nb);
   }
   P=pres;
   T=temp;
   for (i=0;i<nb;i++) {chem[i]->P=pres; chem[i]->T=temp;}
   if(n()>EPS)
   {
      v=quality();
      tmp1=0;
      for (i=0;i<nb;i++)
      {
         if (T>chem[i]->Tc || v==1) chem[i]->state=1;
         else chem[i]->state=0;
      }
      v=0;
      for (i=0;i<nb;i++)
      {
         if (chem[i]->state==0) {if (chem[i]->m>EPS) v+= chem[i]->m/chem[i]->rho(); tab4[i]=0;}
         if (chem[i]->state==1) {tab4[i] = chem[i]->n(); tmp1+=tab4[i];}
      }
      if (tmp1>EPS)
      {
         for (i=0;i<nb;i++)
         {
            tab1[i] = chem[i]->Pc;  //cout<<endl<<tab1[i];
            tab2[i] = chem[i]->Tc; // cout<<"  "<<tab2[i];
            tab3[i] = chem[i]->omega();  //cout<<"  "<<tab3[i];
            tab4[i]/=tmp1;    //          cout<<"  "<<tab4[i];
         }
         thermo->send(tab1,tab2,tab3, tab4);
         thermo->set(P,T,0.0,tmp1);
         v+=thermo->v();
      }
   }
   else v= 0.0;
}

void stream::set(double *list)
{
   m=0;
   for (i=0; i<nb; i++)
   {
      chem[i]->m=list[i];
      m+=list[i];
   }
}

void stream::bubble()
{
   bp=1.1e6;
   for(i=0;i<nb;i++) if(T<chem[i]->Tc && chem[i]->Tboil(P)<bp && chem[i]->m>EPS) bp=chem[i]->Tboil(P);
   if (bp==1.1e6) bp=0.0;
   else
   {
      step=2;
      while (fabs(step)>TOL_BP && fabs(tmp1-1)>TOL_BP)
      {
         //if(DEBUG) cout<<endl<<bp;
         bp+=step;
         tmp1=tmp2=0;
         for (i=0;i<nb;i++) if(T<chem[i]->Tc) {tmp2+=chem[i]->n(); tmp1+=chem[i]->n()*chem[i]->Psat(bp)/P;}
         tmp1/=tmp2;
         step=10*(1.0-tmp1);
      }
   }
}

void stream::dew()
{
   dp=0.0;
   tmp1=10;
   for(i=0;i<nb;i++) if(T<chem[i]->Tc && chem[i]->Tboil(P)>dp && chem[i]->m>EPS) dp=chem[i]->Tboil(P);
   if (dp>EPS)
   {
   dp=bp;
      step=1;
      while (fabs(step)>TOL_DP && fabs(tmp1/tmp2-1)>TOL_DP)
      {
         dp+=step;
         if(dp<bp) dp=bp;
         //if(DEBUG) cout<<endl<<dp;
         tmp1=tmp2=0;
          for (i=0;i<nb;i++) if(T<chem[i]->Tc) {tmp2+=chem[i]->n(); tmp1+=chem[i]->n()/chem[i]->Psat(dp)*P; }
         if (step/fabs(step)*tmp2/tmp1 >1 || step/fabs(step)*tmp1/tmp2 <-1) step*= -0.1;
      }
   }
}

double stream::quality()
{
   if(T>EPS)
   {
      bubble();
      dew();
      //if(DEBUG) cout<<endl<<name<<"  bp="<<bp<<"  dp="<<dp<<"  T="<<T;  system("pause");
      if (bp < dp)
      {
         if (bp < T && T< dp) tmp= (T-bp)/(dp-bp);
         if (T<= bp) tmp= 0.0;
         if (T >= dp) tmp= 1.0;
      }
   }
   else tmp= 0.0;
   return tmp;
}

double stream::K(int i)
{
   for (j=0;j<nb;j++)
     tab4[j] = chem[j]->n()/n();
   thermo->send(tab1,tab2,tab3, tab4);
   if (T>EPS && P>EPS)
   {
      thermo->set(P,T,v,n());
      return thermo->K(i);
   }
   else
   {
	   ofstream logf;
      logf.open(MESSAGES, ios::app);
      logf<<"   --> Warning <--  Cannot compute K of "<<chem[i]->name<<" in stream "<<name<<".\n";
      logf.close();
      warning++;
      check_error();
      return 1.0;
   }
}

double stream::n()
{
   tmp=0.0;
   for ( k = 0 ; k < nb ; k++ )
     tmp+=chem[k]->n();
   return tmp;
}


// affectation :
stream & stream::operator = ( const stream & s ) {
  
  (*thermo) = *(s.thermo);

  delete [] tab1;
  delete [] tab2;
  delete [] tab3;
  delete [] tab4;
  for ( i = 0 ; i < nb ; i++ )
    delete chem[i];
  delete [] chem;
    
  P = s.P;
  T = s.T;
  m = s.m;
  v = s.v;
  i = s.i;
  j = s.j;
  k = s.k;
  error = s.error;
  warning = s.warning;
  name = s.name;
  nb = s.nb;
  chem = new chemical *[nb];
  
  tab1 = new double[nb];
  tab2 = new double[nb];
  tab3 = new double[nb];
  tab4 = new double[nb];

  step = s.step;

  for  ( i = 0 ; i < nb ; i++ ) {
    chem[i] = new chemical(*s.chem[i]);
    tab1[i] = s.tab1[i];
    tab2[i] = s.tab2[i];
    tab3[i] = s.tab3[i];
    tab4[i] = s.tab4[i];
  }

  tmp  = s.tmp;
  tmp1 = s.tmp1;
  tmp2 = s.tmp2;
  
  return *this;
}


void stream::write() {
  
  string file_name = RUNTIME + name + ".stream";

  cout << "WRITE FILE " << file_name << " :\n\tBEGIN\n";

  cout.unsetf(ios::scientific);
  cout << setiosflags(ios::showpoint | ios::fixed)<<setprecision(1);

  cout << "\t>>" << setw(8) << P << "   " << setw(9) << T <<"  ";

  cout << resetiosflags(ios::fixed) << setiosflags(ios::scientific) << setprecision(3);

  if (m>EPS)
    cout << setw(11) << m << setw(11) << v;
  else
    cout << "       x          x   ";
  for  ( i = 0 ; i < nb ; i++ ) {
    if ( chem[i]->m <= EPS )
      cout << "       x   ";
    else
      cout << setw(11) << chem[i]->m;
  }
  cout << "\n\tEND\n\n";

  cout.setf(ios::scientific);

}

void stream::purge()
{
  P=T=m=v=0.0;
  for(i=0;i<nb;i++)
    chem[i]->purge();
}
