#ifndef __HS_MATERIAL__
#define __HS_MATERIAL__

#include <iostream>
#include <iomanip>
#include <fstream>

/*----------------------*/
/*  types of materials  */
/*----------------------*/
enum material_type { _NYLON_        ,     // 0
		     _TEFLON_       ,     // 1
		     _EPOXY_NORMAL_ ,     // 2
		     _EPOXY_PLANE_  ,     // 3
		     _ALUMINIUM_    ,     // 4
		     _STEEL_        ,     // 5
		     _CARBON_STEEL_   };  // 6

// display material type:	
std::ostream & operator << ( std::ostream & , material_type );

/*-----------------*/
/*  material type  */
/*-----------------*/
class HS_Material {

private:

  material_type _type;
  double        _density;
  int           _nT;
  double      * _T;
  int           _nk;
  double      * _k;
  double      * _kA;
  double      * _kB;
  double      * _kC;
  int           _ns;
  double      * _s;
  double      * _sA;
  double      * _sB;
  double      * _sC;
  int           _nTCI;
  double      * _TCI;
  int           _ne;
  double      * _e;
  double      * _ekA;
  double      * _ekB;
  double      * _ekC;
  double      * _ksA;
  double      * _ksB;
  double      * _ksC;
  int           _nTEI;
  double      * _TEI;
  int           _ny;
  double      * _y;
  double      * _yA;
  double      * _yB;
  double      * _yC;

  // tridiagonal solver to determine cubic spline z-values:
  void CSplineABC ( int            n ,
		    const double * x ,
		    const double * y ,
		    double       * A ,
		    double       * B ,
		    double       * C   ) const;
public:

  // conversion constants:
  static const double DENSITY_CC;
  static const double T_CC;
  static const double K_CC;
  static const double S_CC;
  static const double YM_CC;

  // constructor:
  HS_Material  ( material_type t , double d );

  // destructor:
  ~HS_Material ( void );

  // SET methods:
  void set_T      ( int n , double * T , bool convert = true );
  void set_y      ( int n , double * y , bool convert = true );
  void set_k      ( int n , double * k , bool convert = true );
  void set_s      ( int n , double * s , bool convert = true );
  void set_e      ( int n , double * e );
  void set_yABC   ( void );
  void set_kABC   ( void );
  void set_sABC   ( void );
  void set_ekABC1 ( void );
  void set_ekABC2 ( void );
  void set_TCI    ( void );
  void set_ksABC  ( void );
  void set_TEI    ( void );

  // get the T indexes between T1 and T2:
  bool get_T_indexes ( double T1       ,
		       double T2       ,
		       int  & indfirst ,
		       int  & indlast    ) const;

  // other GET methods:
  double get_density ( void    ) const { return _density;  }
  int    get_nT      ( void    ) const { return _nT;       }
  double get_T       ( int ind ) const { return _T  [ind]; }
  double get_k       ( int ind ) const { return _k  [ind]; }
  double get_e       ( int ind ) const { return _e  [ind]; }
  double get_TCI     ( int ind ) const { return _TCI[ind]; }
  double get_TEI     ( int ind ) const { return _TEI[ind]; }

  // thermal expansion:
  void thermal_expansion ( int            nx ,
			   const double * x  ,
			   double       * y    ) const;

  // thermal conductivity:
  void thermal_conductivity ( int            nx ,
			      const double * x  ,
			      double       * y    ) const;

  // thermal stress:
  double thermal_stress ( double T ) const;

  // thermal conductivity/stress:
  void thermal_conductivity_stress ( int            nx ,
				     const double * x  ,
				     double       * y    ) const;

  // thermal YoungModulus:
  void thermal_YM ( int            nx ,
		    const double * x  ,
		    double       * y    ) const;

};

#endif
