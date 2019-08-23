#include "reaction.hpp"
using namespace std;


int reaction::find_chemical ( const string & chem_name ) const {
  for ( int i = 0 ; i < m ; i++ )
    if ( list[i]->CAS == chem_name )
      return i;
  return -1;
}

// donnees hardcodees
reaction::reaction ( const string & in1 , int dim , chemical ** in2 ) {

  m    = dim;  // nbre de chemicals
  list = in2;  // liste des chemicals

  n      = new double[m];
  safe_n = new double[m];
  safe_a = new double[m];
  a      = new double[m];

  int i , j;
  for ( i = 0 ; i < m ; i++ ) {
    a[i]=0.0;
    n[i]=0.0;
  }

  // 1/5 :
  if ( in1 == "eb2sty" ) {
    k0 = 3.525e5;
    E  = 90.85;
    if ( (j = find_chemical ("100-41-4")) < 0 ) {
      cout << "ERROR 10a\n\n";
      exit(0);
    }
    a[j] = -1;
    n[j] =  1;
    if ( (j = find_chemical ("1333-74-0")) < 0 ) {
      cout << "ERROR 10b\n\n";
      exit(0);
    }
    a[j] = 1;
    n[j] = 0;
    if ( (j = find_chemical ("100-42-5")) < 0 ) {
      cout << "ERROR 10c\n\n";
      exit(0);
    }
    a[j] = 1;
    n[j] = 0;
  }

  // 2/5 :
  else if ( in1 == "sty2eb" ) {
    k0 = 2.754e-4;
    E  = -18.653;
    if ( (j = find_chemical ("100-41-4")) < 0 ) {
      cout << "ERROR 10d\n\n";
      exit(0);
    }
    a[j] = 1;
    n[j] = 0;
    if ( (j = find_chemical ("1333-74-0")) < 0 ) {
      cout << "ERROR 10e\n\n";
      exit(0);
    }
    a[j] = -1;
    n[j] =  1;
    if ( (j = find_chemical ("100-42-5")) < 0 ) {
      cout << "ERROR 10f\n\n";
      exit(0);
    }
    a[j] = -1;
    n[j] =  1;
  }

  // 3/5 :
  else if ( in1 == "eb2bz" ) {
    k0 = 9.577e4;
    E  = 111.375;
    if ( (j = find_chemical ("100-41-4")) < 0 ) {
      cout << "ERROR 10g\n\n";
      exit(0);
    }
    a[j] = -1;
    n[j] =  1;
    if ( (j = find_chemical ("71-43-2")) < 0 ) {
      cout << "ERROR 10h\n\n";
      exit(0);
    }
    a[j] = 1;
    n[j] = 0;
    if ( (j = find_chemical ("74-85-1")) < 0 ) {
      cout << "ERROR 10i\n\n";
      exit(0);
    }
    a[j] = 1;
    n[j] = 0;
  }

  // 4/5 :
  else if ( in1 == "eb2tol" ) {
    k0 = 6.077e8;
    E  = 207.850;
    if ( (j = find_chemical ("100-41-4")) < 0 ) {
      cout << "ERROR 10j\n\n";
      exit(0);
    }
    a[j] = -1;
    n[j] =  1;
    if ( (j = find_chemical ("1333-74-0")) < 0 ) {
      cout << "ERROR 10k\n\n";
      exit(0);
    }
    a[j] = -1;
    n[j] =  1;
    if ( (j = find_chemical ("108-88-3")) < 0 ) {
      cout << "ERROR 10l\n\n";
      exit(0);
    }
    a[j] = 1;
    n[j] = 0;
    if ( (j = find_chemical ("74-82-8")) < 0 ) {
      cout << "ERROR 10m\n\n";
      exit(0);
    }
    a[j] = 1;
    n[j] = 0;
  }

  // 5/5 :
  else if ( in1 == "tol2bz" ) {
    k0 = 1;
    E  = 19.038;
    if ( (j = find_chemical ("1333-74-0")) < 0 ) {
      cout << "ERROR 10n\n\n";
      exit(0);
    }
    a[j] = -1;
    n[j] = 0.5;
    if ( (j = find_chemical ("108-88-3")) < 0 ) {
      cout << "ERROR 10o\n\n";
      exit(0);
    }
    a[j] = -1;
    n[j] =  1;
    if ( (j = find_chemical ("71-43-2")) < 0 ) {
      cout << "ERROR 10p\n\n";
      exit(0);
    }
    a[j] = 1;
    n[j] = 0;
    if ( (j = find_chemical ("74-82-8")) < 0 ) {
      cout << "ERROR 10q\n\n";
      exit(0);
    }
    a[j] = 1;
    n[j] = 0;
  }
  else {
    cout << "ERROR 12\n\n";
    exit(0);
  }

  for ( i = 0 ; i < m ; i++ ) {
    safe_n[i]=n[i];
    safe_a[i]=a[i];
  }
}

reaction::~reaction() {
  delete [] a;
  delete [] n;
  delete [] safe_n;
  delete [] safe_a;
}

double reaction::dHr(double T)
{
  int i , j;
  for (i=0;i<m;i++)
    if(safe_a[i]!=a[i])
      {
	if(a[i]>safe_a[i]) a[i]=safe_a[i];
	else safe_a[i]=a[i];
      }
  double tmp=0.0;
  for (i=0;i<m;i++) tmp += a[i]*list[i]->Ho;
  if(fabs(T-298)>EPS)
    for (i=0;i<m;i++)
      for (j=1;j<=4;j++) tmp += a[i]*list[i]->Cp_param[j-1]*(pow(T,j)-pow(298.0,j))/j/1000.0;
  return tmp;
}

double reaction::rate(double T, double* C)
{
  double tmp = k0*exp(-1000*E/8.3144/T);
  for ( int i=0;i<m;i++)
    {
      if(safe_n[i]!=n[i]) n[i]=safe_n[i];
      if(C[i]>EPS && fabs(n[i])>EPS) tmp *= pow(C[i], n[i]);
    }
  return tmp;
}
