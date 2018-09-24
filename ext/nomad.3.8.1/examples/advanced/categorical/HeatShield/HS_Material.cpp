#include "HS_Material.hpp"
using namespace std;

/*------------------------*/
/*  conversion constants  */
/*------------------------*/
const double HS_Material::DENSITY_CC = 2.767974E-2;
const double HS_Material::T_CC       = 1/1.8;
const double HS_Material::K_CC       = 1/57.7882;
const double HS_Material::S_CC       = 689.476;
const double HS_Material::YM_CC      = 6.89476E5;

/*---------------*/
/*  constructor  */
/*---------------*/
HS_Material::HS_Material  ( material_type t ,
			    double        d   ) : _type    ( t              ) ,
						  _density ( d * DENSITY_CC ) ,
						  _nT      ( 0              ) , 
						  _T       ( NULL           ) ,
						  _nk      ( 0              ) ,
						  _k       ( NULL           ) ,
						  _kA      ( NULL           ) ,
						  _kB      ( NULL           ) ,
						  _kC      ( NULL           ) ,
						  _ns      ( 0              ) ,
						  _s       ( NULL           ) ,
						  _sA      ( NULL           ) ,
						  _sB      ( NULL           ) ,
						  _sC      ( NULL           ) ,
						  _nTCI    ( 0              ) ,
						  _TCI     ( NULL           ) ,
						  _ne      ( 0              ) ,
						  _e       ( NULL           ) ,
						  _ekA     ( NULL           ) ,
						  _ekB     ( NULL           ) ,
						  _ekC     ( NULL           ) ,
						  _ksA     ( NULL           ) ,
						  _ksB     ( NULL           ) ,
						  _ksC     ( NULL           ) ,
						  _nTEI    ( 0              ) ,
						  _TEI     ( NULL           ) ,
						  _ny      ( 0              ) ,
						  _y       ( NULL           ) ,
						  _yA      ( NULL           ) ,
						  _yB      ( NULL           ) ,
						  _yC      ( NULL           )   {}
/*--------------*/
/*  destructor  */
/*--------------*/
HS_Material::~HS_Material ( void ) 
 {
  delete [] _T;
  delete [] _k;
  delete [] _kA;
  delete [] _kB;
  delete [] _kC;
  delete [] _s;
  delete [] _sA;
  delete [] _sB;
  delete [] _sC;
  delete [] _TCI;
  delete [] _e;
  delete [] _ekA;
  delete [] _ekB;
  delete [] _ekC;
  delete [] _ksA;
  delete [] _ksB;
  delete [] _ksC;
  delete [] _TEI;
  delete [] _y;
  delete [] _yA;
  delete [] _yB;
  delete [] _yC;
}

/*---------------*/
/*  SET methods  */
/*---------------*/

// set_T:
void HS_Material::set_T ( int n , double * T , bool convert ) 
 {
  delete [] _T;
  double CC = (convert) ? T_CC : 1.0;
  _nT = n;
  _T  = new double [_nT];
  for ( int i = 0 ; i < _nT ; ++i )
    _T[i] = T[i] * CC;
}

// set_y:
void HS_Material::set_y ( int n , double * y , bool convert ) 
 {
  delete [] _y;
  double CC = (convert) ? YM_CC : 1.0;
  _ny = n;
  _y  = new double [_ny];
  for ( int i = 0 ; i < _ny ; ++i )
    _y[i] = y[i] * CC;
}

// set_k:
void HS_Material::set_k ( int n , double * k , bool convert ) 
 {
  delete [] _k;
  double CC = (convert) ? K_CC : 1.0;
  _nk = n;
  _k  = new double [_nk];
  for ( int i = 0 ; i < _nk ; ++i )
    _k[i] = k[i] * CC;
}

// set_yABC:
void HS_Material::set_yABC ( void ) 
 {
  delete [] _yA;
  delete [] _yB;
  delete [] _yC;
  _yA = new double [ _ny-1 ];
  _yB = new double [ _ny-1 ];
  _yC = new double [ _ny-1 ];

  double * T = new double [ _ny ];
  int i1 = _ny-1;
  int i2 = _nT-1;
  while ( i1 >= 0 )
    T[i1--] = _T[i2--];

  CSplineABC ( _ny-1 , T , _y , _yA , _yB , _yC );

  delete [] T;
}

// set_kABC:
void HS_Material::set_kABC ( void ) 
 {
  delete [] _kA;
  delete [] _kB;
  delete [] _kC;
  _kA = new double [ _nk-1 ];
  _kB = new double [ _nk-1 ];
  _kC = new double [ _nk-1 ];   
  CSplineABC ( _nk-1 , _T , _k , _kA , _kB , _kC );  
}

// set_s:
void HS_Material::set_s ( int n , double * s , bool convert ) 
 {
  delete [] _s;
  double CC = (convert) ? S_CC : 1.0;
  _ns = n;
  _s  = new double [_ns];
  for ( int i = 0 ; i < _ns ; ++i )
    _s[i] = s[i] * CC;
}

// set_e:
void HS_Material::set_e ( int n , double * e ) 
 {
  delete [] _e;
  _ne = n;
  _e  = new double [_ne];
  for ( int i = 0 ; i < _ne ; ++i )
    _e[i] = e[i];
}

// set_sABC:
void HS_Material::set_sABC ( void ) 
 {
  delete [] _sA;
  delete [] _sB;
  delete [] _sC;
  _sA = new double [ _ns-1 ];
  _sB = new double [ _ns-1 ];
  _sC = new double [ _ns-1 ];   

  // hypothesis: _ns may be smaller than _nT:
  double * T = new double [ _ns ];
  int i1 = _ns-1;
  int i2 = _nT-1;
  while ( i1 >= 0 )
    T[i1--] = _T[i2--];

  CSplineABC ( _ns-1 , T , _s , _sA , _sB , _sC );  

  delete [] T;
}

// set_ksABC:
void HS_Material::set_ksABC ( void ) 
 {
  delete [] _ksA;
  delete [] _ksB;
  delete [] _ksC;

  _ksA = new double [_nT-1];
  _ksB = new double [_nT-1];
  _ksC = new double [_nT-1];

  double * ks = new double [_nT];

  for ( int i = 0 ; i < _nT ; ++i )
    ks[i] = _k[i] / _s[i];

  CSplineABC ( _nT-1 , _T , ks , _ksA , _ksB , _ksC );

  delete [] ks;
}

// set_ekABC1:
void HS_Material::set_ekABC1 ( void ) 
 {

  delete [] _ekA;
  delete [] _ekB;
  delete [] _ekC;

  _ekA = new double [_ne-1];
  _ekB = new double [_ne-1];
  _ekC = new double [_ne-1];

  double * T = new double [ _ne ];
  int i1 = _ne-1;
  int i2 = _nT-1;
  while ( i1 >= 0 )
    T[i1--] = _T[i2--];
  T[0] = 0.0;

  CSplineABC ( _ne-1 , T , _e , _ekA , _ekB , _ekC );

  delete [] T;
}

// set_ekABC2:
void HS_Material::set_ekABC2 ( void ) 
 {

  delete [] _ekA;
  delete [] _ekB;
  delete [] _ekC;

  _ekA = new double [_nT-1];
  _ekB = new double [_nT-1];
  _ekC = new double [_nT-1];

  double * ek = new double [ _nT ];
  for ( int i = 0 ; i < _nT ; ++i )
    ek[i] = _e[i] * _k[i];

  CSplineABC ( _nT-1 , _T , ek , _ekA , _ekB , _ekC );

  delete [] ek;
}

// set_TCI:
void HS_Material::set_TCI ( void ) 
 {
  
  delete [] _TCI;

  int       i;
  int      nx = _nT - 1;
  double *  x = new double [nx];
  double * fx = new double [nx];

  for ( i = 0 ; i < nx ; ++i )
    x[i] = ( _T[i] + _T[i+1] ) / 2.0;
  
  thermal_conductivity ( nx , x , fx );
  
  _nTCI = nx;
  _TCI  = new double[nx];
  for ( i = 0 ; i < nx ; ++i )
    _TCI[i] = (_k[i]+_k[i+1]+4*fx[i]) * (_T[i+1]-_T[i]) / 6.0;
   
  delete [] x;
  delete [] fx;
}

// set_TEI:
void HS_Material::set_TEI ( void ) 
 {

  delete [] _TEI;

  int       i;
  int      nx  = _nT - 1;
  double *  x  = new double [nx];
  double * fx1 = new double [nx];
  double * fx2 = new double [nx];

  for ( i = 0 ; i < nx ; ++i )
    x[i] = ( _T[i] + _T[i+1] ) / 2.0;

  thermal_conductivity_stress ( nx , x , fx1 );
  thermal_expansion ( nx , x , fx2 );
 
  _nTEI = nx;
  _TEI  = new double[nx];
  for ( i = 0 ; i < nx ; ++i )
    _TEI[i] = ( _e[i]*_k[i] + _e[i+1]*_k[i+1] + 4*fx1[i]*fx2[i] ) *
      (_T[i+1]-_T[i]) / 6;

  delete [] x;
  delete [] fx1;
  delete [] fx2;
}

/*-------------------------------------------------*/
/*  Material: get the T indexes between T1 and T2  */
/*-------------------------------------------------*/
bool HS_Material::get_T_indexes ( double T1       ,
				  double T2       ,
				  int  & indfirst ,
				  int  & indlast    ) const {
  indfirst = -1;
  indlast  = -1;

  if ( T2 <= _T[0] || T1 >= _T[_nT-1] || T2 <= T1 )
    return false;

  indfirst = 0;
  while ( _T[indfirst] <= T1 )
    ++indfirst;
  
  indlast = _nT-1;
  while ( _T[indlast] >= T2 )
    --indlast;

  if ( indfirst > indlast ) 
 {
    indfirst = indlast = -1;
    return false;
  }

  return true;
}

/*---------------------*/
/*  thermal expansion  */
/*---------------------*/
void HS_Material::thermal_expansion ( int            nx ,
				      const double * x  ,
				      double       * y    ) const {
  double w;
  int    ind;
  for ( int i = 0 ; i < nx ; ++i ) 
 {
    if ( x[i] <= _T[0] )
      ind = 0;
    else if ( x[i] >= _T[_nT-1] )
      ind = _nT-2;
    else 
 {
      ind = 0;
      while ( x[i] >= _T[ind+1] )
	++ind;
    }
    w = x[i] - _T[ind];
    y[i] = _e[ind]*_k[ind] + w * ( _ekC[ind]  + w * ( _ekB[ind]  + w * _ekA[ind] ) );
  }
}

/*-------------------------------*/
/*  thermal conductivity/stress  */
/*-------------------------------*/
void HS_Material::thermal_conductivity_stress ( int            nx ,
						const double * x  ,
						double       * y    ) const {
  double w;
  int    ind;
  for ( int i = 0 ; i < nx ; ++i ) 
 {
    if ( x[i] <= _T[0] )
      ind = 0;
    else if ( x[i] >= _T[_nT-1] )
      ind = _nT-2;
    else 
 {
      ind = 0;
      while ( x[i] >= _T[ind+1] )
	++ind;
    }
    w = x[i] - _T[ind];
    y[i] = _k[ind] * _s[ind] + w  * ( _ksC[ind] + w * ( _ksB[ind] + w * _ksA[ind] ) );
  }
}

/*------------------------*/
/*  thermal conductivity  */
/*------------------------*/
void HS_Material::thermal_conductivity ( int            nx ,
					 const double * x  ,
					 double       * y    ) const {
  double w;
  int    ind;
  for ( int i = 0 ; i < nx ; ++i ) 
 {
    if ( x[i] <= _T[0] )
      ind = 0;
    else if ( x[i] >= _T[_nT-1] )
      ind = _nT-2;
    else 
 {
      ind = 0;
      while ( x[i] >= _T[ind+1] )
	++ind;
    }
    w = x[i] - _T[ind];
    y[i] = _k[ind] + w * ( _kC[ind]  + w * ( _kB[ind]  + w * _kA[ind] ) );
  }
}

/*------------------------*/
/*  thermal YoungModulus  */
/*------------------------*/
void HS_Material::thermal_YM ( int            nx ,
			       const double * x  ,
			       double       * y    ) const {
  double w;
  int    ind;
  for ( int i = 0 ; i < nx ; ++i ) 
 {
    if ( x[i] <= _T[0] )
      ind = 0;
    else if ( x[i] >= _T[_nT-1] )
      ind = _nT-2;
    else 
 {
      ind = 0;
      while ( x[i] >= _T[ind+1] )
	++ind;
    }
    w = x[i] - _T[ind];

    y[i] = _y[ind] + w * ( _yC[ind]  + w * ( _yB[ind]  + w * _yA[ind] ) );
  }
}

/*------------------*/
/*  thermal stress  */
/*------------------*/
double HS_Material::thermal_stress ( double T ) const {

  int ind;
  if ( T <= _T[0] )
    ind = 0;
  else if ( T >= _T[_nT-1] )
    ind = _nT-2;
  else 
 {
    ind = 0;
    while ( T >= _T[ind+1] )
      ++ind;
  }

  double w = T - _T[ind];

  return _s[ind] + w * ( _sC[ind]  + w * ( _sB[ind]  + w * _sA[ind] ) );
}

/*---------------------------------------------------------*/
/*  tridiagonal solver to determine cubic spline z-values  */
/*---------------------------------------------------------*/
void HS_Material::CSplineABC ( int            n ,
			       const double * x ,            // size n+1
			       const double * y ,            // size n+1
			       double       * A ,            // size n
			       double       * B ,            // size n
			       double       * C   ) const {  // size n

  double * h = new double [n];
  double * b = new double [n];
  double * u = new double [n-1];
  double * v = new double [n-1];

  int i;

  for ( i = 0 ; i < n ; ++i ) 
 {
    h[i] = x[i+1] - x[i];
    b[i] = 6*(y[i+1]-y[i])/h[i];
    if ( i < n-1 ) 
 {
      u[i] = 2 * ( h[i] + x[i+2] - x[i+1] );
      v[i] = -b[i] +  6 * ( y[i+2]-y[i+1] ) / (x[i+2]-x[i+1]);
    }
  }

  for ( i = 1 ; i < n-1 ; ++i ) 
 {
    u[i] -= h[i]*h[i]  / u[i-1];
    v[i] -= h[i]*v[i-1]/ u[i-1];
  }

  double * z = new double[n+1];
  z[0] = z[n] = 0;

  for ( i = n-1 ; i > 0 ; --i )
    z[i] = (v[i-1] - h[i]*z[i+1])/u[i-1];

  // compute cubic spline coefficients:
  for ( i = 0 ; i < n ; ++i ) 
 {
    A[i] = (z[i+1] - z[i]) / (6*h[i]);
    B[i] = z[i] / 2;
    C[i] = -h[i]*z[i+1]/6 - h[i]*z[i]/3 + (y[i+1]-y[i])/h[i];
  }

  delete [] h;
  delete [] b;
  delete [] u;
  delete [] v;
  delete [] z;
}

/*-------------------------*/
/*  display material type  */
/*-------------------------*/
ostream & operator << ( ostream & out , material_type mt ) 
 {

  switch ( mt ) 
 {
  case _NYLON_:
    out << "nylon";
    break;
  case _TEFLON_:
    out << "teflon";
    break;
  case _EPOXY_NORMAL_:
    out << "epoxy_normal";
    break;
  case _EPOXY_PLANE_:
    out << "epoxy_plane";
    break;
  case _ALUMINIUM_:
    out << "aluminum";
    break;
  case _STEEL_:
    out << "steel";
    break;
  case _CARBON_STEEL_:
    out << "carbon_steel";
    break;
  default:
    out << "undefined";
  }
  return out;
}
