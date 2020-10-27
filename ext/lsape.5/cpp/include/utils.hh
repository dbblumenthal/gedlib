// -----------------------------------------------------------   
/** \file utils.h
 *  \brief Manipulation of assignments/matchings
 * \author Sebastien Bougleux (Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, Caen, France)
*/
/* -----------------------------------------------------------
   This file is part of LIBLSAP.
   
   LIBLSAP is free software: you can redistribute it and/or modify
   it under the terms of the CeCILL-C License. See README for more
   details.

   -----------------------------------------------------------
   Creation: December 5 2015
   Last modif: March 2018
*/

#ifndef __UTILS_HH__
#define __UTILS_HH__

#define ZERO_EPS 1e-12

#include <random>
#include<iostream>
namespace liblsap {

  // -----------------------------------------------------------
  /*template <typename DT>
  inline bool isZero(const DT &e)
  { return (std::abs(e) < (DT)1e-12); }
  */
  // -----------------------------------------------------------
  int randInt(unsigned int mn, unsigned int mx)
  { return rand()%(mx-mn+1) + mn; }

  // -----------------------------------------------------------
  template <typename DT, typename IndexType = int>
  void randVecStoch(const IndexType &n, DT *v)
  {
    std::random_device rd;
    std::mt19937 e2(rd());
    std::normal_distribution<DT> distribution(0,1.0);
    for (int i = 0; i < n; i++) v[i] = distribution(e2);
  }

  // -----------------------------------------------------------
  template <typename IT>
  IT factorial(const IT &n)
  {
    IT f = 1;
    for (IT i = 1; i <= n; i++) f *= i;
    return f;
  }

  // -----------------------------------------------------------
  inline int compute_s(double f, double C)
  {
    int exp;
    std::frexp(f*C,&exp);
    return 53-exp;
  }

  // -----------------------------------------------------------
  inline double compute_S(double f, double C, double &one_over_S)
  {
    int exp;
    std::frexp(f*C,&exp);
    one_over_S = std::ldexp(1,exp-53);
    return std::ldexp(1,53-exp);
  }

  // -----------------------------------------------------------
  inline double scale_data(double d, double S, double one_over_S)
  {
    if (d == 0) return 0;
    int sign_d = 1;
    if (d < 0) { sign_d = -1; d = -d; }
    return sign_d * std::floor(d*S) * one_over_S;
  }

  // -----------------------------------------------------------
  inline double scale_data(double d, double S)
  {
    if (d == 0) return 0;
    int sign_d = 1;
    if (d < 0) { sign_d = -1; d = -d; }
    return sign_d * std::floor(d*S);
  }

  // -----------------------------------------------------------
  bool arith_scale(int nr, int nc, double *d, double *dres, double f, double dmax)
  {
    int n = nr*nc;
    // compute scaling parameters
    double one_over_S;
    double S = compute_S(f,dmax,one_over_S);
    // scale values in d
    bool scaling = false;
    double dtmp = 0;
    for (double *dit = d, *ditres = dres, *dend = d+n; dit != dend; dit++, ditres++) {
      dtmp = *dit;
      *ditres = scale_data(dtmp,S,one_over_S);
      if (*ditres != dtmp) scaling = true;
    }
    return scaling;
  }

  // -----------------------------------------------------------
  bool arith_scale(int nr, int nc, double *d, long int *dres, double f, double dmax)
  {
    int n = nr*nc;
    // compute scaling parameters
    double one_over_S;
    double S = compute_S(f,dmax,one_over_S);
    // scale values in d
    bool scaling = false;
    double dtmp = 0, dtmp2 = 0;
    long int *ditres = dres;
    for (double *dit = d, *dend = d+n; dit != dend; dit++, ditres++) {
      dtmp = *dit;
      dtmp2 = scale_data(dtmp,S);
      *ditres = static_cast<long int>(dtmp2);
      if (dtmp != dtmp2) scaling = true;
    }
    return scaling;
  }

  // -----------------------------------------------------------
  bool arith_scale(int nr, int nc, double *d, double *dres, double f)
  {
    int n = nr*nc;
    // compute max value in d
    double dmax = std::numeric_limits<double>::max();
    for (double *dit = d, *dend = d+n; dit != dend; dit++)
      if (*dit > dmax) dmax = *dit;
    return arith_scale(nr,nc,d,dres,dmax);
  }
  
} // end namespace

#endif
