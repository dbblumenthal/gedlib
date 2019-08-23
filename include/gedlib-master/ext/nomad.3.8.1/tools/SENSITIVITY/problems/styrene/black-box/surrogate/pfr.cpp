#include "pfr.hpp"
#include "RungeKutta.cpp"
using namespace std;

pfr::pfr ( stream * s1 , stream * s2 , double ** t , int nb_r , reaction ** rr , double u , double ta ) {

  F = s2;
  F->m=0;
  P=s1->P;
  for ( i = 0 ; i < s1->nb ; i++ ) {
    F->chem[i]->m = s1->chem[i]->m;
    F->m+=F->chem[i]->m;
  }
  F->set(s1->P, s1->T);
  m_in=F->m;
  a = t;
  rx = rr;
  n=nb_r;
  m=  F->nb;
  U=u;
  Ta=ta;
  T = F->T;
  C = new double[m];
  y = new double[m+1];
  r=new double[n];
  OK=true;
  explode=true;
  solver = new RungeKutta<pfr>(m+1);
}

pfr::~pfr() {
  delete [] r;
  delete [] C;
  delete [] y;
  delete solver;
}

bool pfr::run() {

  for ( i = 0 ; i < m ; i++ )
    y[i]=F->chem[i]->n();

  y[m]=T;

  solver->set ( this , y , 0.0 , L );

  dL=solver->dx();

  OK=solver->run();

  sum = F->m;
  F->m = 0;

  for ( i = 0 ; i < m ; i++ )
    F->m+=F->chem[i]->m;
  for ( i = 0 ; i < m ; i++ )
    F->chem[i]->m *= sum/F->m;
  // if (OK)
  //mass balance!
  if ( fabs(m_in-F->m) > EPS || !explode )
    OK=false;

  return OK;
}

double pfr::f ( int eq , double l , double * y ) {


  sum=F->m;
  F->m=0;
  for(i=0; i<m;i++)
    {
      if(y[i]<0) y[i]=0;
      F->chem[i]->m = y[i]*F->chem[i]->M/1000.0;
      F->m+=F->chem[i]->m;
    }


  for(i=0; i<m;i++)
    F->chem[i]->m *= sum/F->m;


  F->m=sum;
  T=y[m];
  if(T>MAX_TEMP)
    {
      cout << "ERROR 11\n\n";
      exit(0);
    }


  for(i=0; i<m;i++)
    C[i]=F->chem[i]->n()/F->v;

  for(j=0;j<n;j++)
    r[j] = rx[j]->rate(T,C);

  if(0<=eq && eq<m)  //return dFi/dL
    {
      tmp=0.0;
      for(j=0;j<n;j++) tmp+=a[eq][j]*r[j];
      tmp *= (pi*D*D/4.0);
    }



  if(eq==m)   //return dT/dL
    {


      F->set(F->P,T);

      tmp=0.0;
      for(j=0;j<n;j++)
	tmp -= r[j]*rx[j]->dHr(T);


      tmp *= (pi*D*D/4.0);


      tmp += (pi*D)*U*(Ta-T);


      tmp1=0.0;
      for(i=0;i<m;i++)
	tmp1+= y[i]*F->chem[i]->Cp()*0.001;
      tmp /= tmp1;


      if(fabs(tmp*dL)>500.0)
	{
	  cout << "ERROR 13\n\n";
	  exit(0);                           
	}
    }




  return tmp;
}



double pfr::get_cost ( void ) {
  dL=L*pi*pow(D,2)/4.0;
  if(dL<0.3) dL=0.3; if(dL>520) dL=520;
  sum = 3.4974+0.4485*log10(dL)+0.1074*pow(log10(dL),2);
  sum = pow(10, sum);
  P= (P-1)*101.325/100;
  dL=(P+1)*D/(317.46*(850-0.6*(P+1)))+0.0315;
  sum *=(2.25+ 1.82*dL*4.2);
  sum = sum*MS_YEAR/MS_2001;
  return sum;
}

double pfr::get_water() {
  sum = (U>EPS && T>Ta) ? U*L*pi*pow(D,2)/4*(T-Ta)/4.185/25.0 : 0.0;
  return sum;
}

void pfr::cost() {
  cout << "WRITE FILE " << RUNTIME << name << ".cost" << " :\n\tBEGIN\n";
  cout << "\t>>" << get_cost();
  cout << "\n\tEND\n\n";
}

void pfr::water() {
  cout << "WRITE FILE " << RUNTIME << name << ".water" << " :\n\tBEGIN\n";
  if (U>EPS && T>Ta) sum = (U*L*pi*pow(D,2)/4*(T-Ta)/4.185/25.0);
  else sum = 0.0;
  cout << "\t>>" << sum;
  cout << "\n\tEND\n\n";
}
