// -----------------------------------------------------------   
/** \file lsape.hh
 *  \brief Main file for solving the LSAPE / minimal-cost error-correcting bipartite graph matching problem
 * \author Sebastien Bougleux (Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, Caen, France)
*/
   
/* -----------------------------------------------------------
   
   This file is part of LSAPE.
   
   LSAPE is free software: you can redistribute it and/or modify
   it under the terms of the CeCILL-C License. See README for more
   details.

   -----------------------------------------------------------
   
   Creation: Oct. 2017
   Last modif: Oct. 2017
   
*/

#ifndef __LSAPE__
#define __LSAPE__

#include <cstring>
#include <limits>
#include "hungarian-lsap.hh"
#include "bistoch.hh"

// -----------------------------------------------------------
template <typename DT, typename IT>
DT ecPermCost(const DT *C, const IT &nr, const IT &nc, const IT *rho, const IT *varrho)
{
  DT sumcost = 0;
  const IT n = nr-1, m = nc-1;
  for (int i = 0; i < n; i++) sumcost += C[rho[i]*nr+i];
  for (int j = 0; j < m; j++)
    if (varrho[j] == n) sumcost += C[j*nr+varrho[j]];
  return sumcost;
}

// -----------------------------------------------------------
template <typename DT, typename IT>
bool triangularInequality(const DT *C, const IT &nr, const IT &nc)
{
  IT m = nc-1, n = nr-1, i, j;
  IT n1m = (n+1)*m;
  const DT *ptC = C, *ptIns;

  for (j = 0; j < m; j++, ++ptC)
  {
    ptIns = ptC + n;
    for (i = 0; i < n; i++, ++ptC)
    {
      if (*ptC > *(ptC+n1m) + *ptIns) return false;
    }
  }
  return true;
}

// -----------------------------------------------------------
template <typename DT, typename IT>
bool triangularInequality(const DT *C, const IT &nr, const IT &nc, IT i, IT j)
{ return (C[j*nr+i] <= C[(nc-1)*nr+i] + C[j*nr+(nr-1)]); }

// -----------------------------------------------------------
template <typename DT, typename IT>
void ecMtxStochasticBarycenter(IT nr, IT nc, DT *J)
{
    const DT val = ((DT)2) / (nr + nc);
    std::memset(J,val,sizeof(DT)*(nr*nc-1));
    J[nr*nc-1] = 0;
}

// -----------------------------------------------------------
template <typename DT, typename IT>
DT* randBiStochExt(const IT &nr, const IT &nc, DT tol = 1e-8, int nbit = 50)
{
  IT nm = nr + nc, i, j, n1 = nr+1, m1 = nc+1;
  DT *X = new DT[nm*nm];
  DT *A = new DT[std::max(n1,m1)];
  DT *Z = new DT[std::max(n1,m1)], sum;
  DT *ptA = A, *ptZ = Z, *ptX = X;
  
  for (j = 0; j < nc; j++)
  {
    randVecStoch(n1,A);
    randVecStoch(n1,Z);
    ptA = A; ptZ = Z;
    sum = 0;
    for (i = 0; i < n1; i++, ++ptA, ++ptZ)
    {
      *ptA = (*ptA) * (*ptA) + (*ptZ) * (*ptZ);
      sum += *ptA;
    }
    ptA = A;
    for (i = 0; i < nr; i++, ++ptA, ++ptX) *ptX = *ptA / sum;
    for (i = nr; i < nm; i++, ++ptX) *ptX = 0;
    X[j*nm+j+nr] = *ptA / sum;
  }

  for (j = nc; j < nm; j++)
  {
    randVecStoch(m1,A);
    randVecStoch(m1,Z);
    ptA = A; ptZ = Z;
    sum = 0;
    for (i = 0; i < m1; i++, ++ptA, ++ptZ)
    {
      *ptA = (*ptA) * (*ptA) + (*ptZ) * (*ptZ);
      sum += *ptA;
    }
    ptA = A;
    for (i = 0; i < nr; i++, ++ptX) *ptX = 0;
    for (i = nr; i < nm; i++, ++ptA, ++ptX) *ptX = *ptA / sum;
    X[j*nm+j-nc] = *ptA / sum;
  }

  sinkhornKnopp(X,nm,nm,tol,nbit);

  delete[] A, delete[] Z;
  return X;
}
// -----------------------------------------------------------
template <typename DT,typename IT>
void reduceExt(const DT *E, const IT &nr, const IT &nc, DT *R)
{
  // copy sub
  IT i, j, nm = nr+nc;
  for (j = 0; j < nc; j++)
  {
    for (i = 0; i < nr; i++)
      R[j*(nr+1)+i] = E[j*nm+i];
    R[j*(nr+1)+nr] = E[j*nm+j+nr];
  }
  for (i = 0; i < nr; i++)
    R[nc*(nr+1)+i] = E[(i+nc)*nm+i];
  R[(nr+1)*(nc+1)-1] = 0;
}

// -----------------------------------------------------------
/** \brief Extend a (n+1)x(m+1) LSAPE instance to an equivalent (n+m)x(m+n) LSAP instance
 *  \param[in] C nxm edit cost matrix
 *  \param[in] n number of rows (last row correspond to the null element)
 *  \param[in] m number of colums (last column correspond to the null element)
 *  \param[out] Cext (n-1+m-1)x(m-1+n-1) extended cost matrix (must be previously allocated)
 *  \param[in] valmx value given to outer diagonal elements of Cext for removals and insertions
 */
template <typename DT,typename IT>
void EBPinst(const DT *C, const IT &n, const IT &m, DT *Cext, DT valmx = std::numeric_limits<DT>::max())
{
  IT msub = m-1, nsub = n-1, i, j;
  IT mpn = msub+nsub;
  // copy subsitutions
  for (j = 0; j < msub; j++)
    std::memcpy(Cext+j*mpn,C+j*n,sizeof(DT)*nsub);
  // copy insertions
  for (j = 0; j < msub; j++)
    for (i = nsub; i < mpn; i++)
      if (i != j+nsub) Cext[j*mpn+i] = valmx;
      else Cext[j*mpn+i] = C[j*n+nsub];
  // copy removals
  for (j = msub; j < mpn; j++)
    for (i = 0; i < nsub; i++)
      if (i+msub != j) Cext[j*mpn+i] = valmx;
      else Cext[j*mpn+i] = C[msub*n+i];
  // set completness
  for (i = nsub; i < mpn; i++)
    for (j = msub; j < mpn; j++)
      Cext[j*mpn+i] = 0;
}

// -----------------------------------------------------------
template <typename DT, typename IT>
void EBP(const DT *C, const IT &nr, const IT &nc, IT *rho, DT *u = NULL, DT *v = NULL, IT *varrho = NULL, unsigned short init_type = 1)
{
  IT npm = nr+nc-2, i = 0, j;
  IT npm2 = npm*npm;
  DT *Cext = new DT[npm2];

  // find max of C
  DT mx = std::numeric_limits<DT>::min();
  for (const DT *ptC = C; i < npm2; i++, ++ptC) if (*ptC > mx) mx = *ptC;

  // construct the instance of the symmetric LSAP
  EBPinst(C,nr,nc,Cext,mx+10);
  
  DT *ux = new DT[npm], *vx = new DT[npm];
  IT *rhox = new IT[npm], *varrhox = NULL;
  if (varrho) varrhox = new IT[npm];

  hungarianLSAP(Cext,npm,npm,rhox,ux,vx,varrhox,init_type);
  
  // reconstruction
  IT n = nr-1, m = nc-1;
  for (i = 0; i < n; i++)
    if (rhox[i] >= m) rho[i] = m;
    else rho[i] = rhox[i];
  if (varrho)
    for (j = 0; j < m; j++)
      if (varrhox[j] >= n) varrho[j] = n;
      else varrho[j] = varrhox[j];
  
  delete[] Cext; delete[] rhox; if (varrho) delete[] varrhox;
  delete[] ux; delete[] vx;
}

// -----------------------------------------------------------
/** \brief Reduce a (n+1)x(m+1) LSAPE instance to an equivalent nxm LSAP instance for RBP
 *  \param[in] C nxm edit cost matrix
 *  \param[in] nr number of rows (last row correspond to the null element)
 *  \param[in] nc number of colums (last column correspond to the null element)
 *  \param[out] Cred (n-1)x(m-1) reduced cost matrix (must be previously allocated)
 */

#include<iostream>
template <typename DT, typename IT>
void RBPinst(const DT *C, const IT &nr, const IT &nc, DT *Cred)
{
  IT m = nc-1, n = nr-1, i, j;
  IT nm = (n+1)*m;
  const DT *ptC = C, *ptIns = C;
  DT *ptR = Cred, cd;

  if (n < m)
  {
    for (j = 0; j < m; j++, ++ptC)
    {
      ptIns = ptC + n;
      for (i = 0; i < n; i++, ++ptC, ++ptR)
      {
        cd = *ptC - *ptIns;
        *ptR = std::min(cd,C[m*nr+i]);
      }
    }
  }
  else if (n == m)
  {
    for (j = 0; j < m; j++, ++ptC)
    {
      ptIns = ptC + n;
      for (i = 0; i < n; i++, ++ptC, ++ptR)
      {
        cd = C[m*nr+i] + *ptIns;
      	*ptR = std::min(*ptC, cd);
      }
    }
  }
  else
  {
    for (j = 0; j < m; j++, ++ptC)
    {
      ptIns = ptC + n;
      for (i = 0; i < n; i++, ++ptC, ++ptR)
      {
        cd = *ptC - C[m*nr+i];  // C[j*nr+i]-C[m*nr+i];
        *ptR = std::min(cd,*ptIns); // C[j*nr+n]);
      }
    }
  }
}
// -----------------------------------------------------------
template <typename DT, typename IT>
void RBP(const DT *C, const IT &nr, const IT &nc, IT *rho, DT *u = NULL, DT *v = NULL, IT *varrho = NULL, unsigned short init_type = 1)
{
  IT n = nr-1, m = nc-1, i, j;
  DT *Cred = new DT[n*m];
  RBPinst(C,nr,nc,Cred);
  
  // Cred contains negative values -> positive
  if (n != m)
  {
    DT mn = std::numeric_limits<DT>::max();
    for (j = 0; j < m; j++)
      for (i = 0; i < n; i++)
      {
        const DT &cr = Cred[j*n+i];
        if (cr < mn) mn = cr;
      }
    if (mn < 0)
    {
        for (j = 0; j < m; j++)
            for (i = 0; i < n; i++)
            Cred[j*n+i] -= mn;
    }

  }

  hungarianLSAP(Cred,n,m,rho,u,v,varrho,init_type);
  
  // reconstruction

  if (n > m)
  {
    for (j = 0; j < m; j++)
    {
        i = varrho[j];
        if (varrho[j] == -1) std::cerr << "bizzzz!!!!!!!\n";
        if (!triangularInequality(C,nr,nc,i,j))
        { rho[i] = m; varrho[j] = n; }
    }
    for (i = 0; i < n; i++) if (rho[i] == -1) rho[i] = m;
  }
  /*
  if (n < m)
  {
    if (varrho)
      for (j = 0; j < m; j++)
	if (varrho[j] == -1) varrho[j] = n;
  }
  else if (n > m)
    for (i = 0; i < n; i++) if (rho[i] == -1) rho[i] = m;
  */
  delete[] Cred;
}
// -----------------------------------------------------------
/** \brief Reduce a (n+1)x(m+1) LSAPE instance to an equivalent nxm LSAP instance for RBP
 *  \param[in] C nxm edit cost matrix
 *  \param[in] nr number of rows (last row correspond to the null element)
 *  \param[in] nc number of colums (last column correspond to the null element)
 *  \param[out] Cred (n-1)x(m-1) reduced cost matrix (must be previously allocated)
 */
template <typename DT, typename IT>
void ARBPinst(const DT *C, const IT &nr, const IT &nc, DT *Cred)
{
  IT m = nc-1, n = nr-1, i, j;
  IT nm = (n+1)*m;
  const DT *ptC = C, *ptIns = C;
  DT *ptR = Cred, cd;

  if (n < m)
  {
    for (j = 0; j < m; j++, ++ptC)
    {
      ptIns = ptC + n;
      for (i = 0; i < n; i++, ++ptC, ++ptR)
      {
        cd = *ptC - *ptIns;
        *ptR = std::min(cd,C[m*nr+i]);
      }
    }
  }
  else if (n == m)
  {
    for (j = 0; j < m; j++, ++ptC)
    {
      ptIns = ptC + n;
      for (i = 0; i < n; i++, ++ptC, ++ptR)
      {
        cd = C[m*nr+i] + *ptIns;
      	*ptR = std::min(*ptC, cd);
      }
    }
  }
  else
  {
    for (j = 0; j < m; j++, ++ptC)
    {
      ptIns = ptC + n;
      for (i = 0; i < n; i++, ++ptC, ++ptR)
      {
        //cd = *ptC - C[m*nr+i];  // C[j*nr+i]-C[m*nr+i];
        *ptR = *ptC - C[m*nr+i];
      }
    }
  }
}
// -----------------------------------------------------------
template <typename DT, typename IT>
void ARBP(const DT *C, const IT &nr, const IT &nc, IT *rho, DT *u = NULL, DT *v = NULL, IT *varrho = NULL, unsigned short init_type = 1)
{
  IT n = nr-1, m = nc-1, i, j;
  DT *Cred = new DT[n*m];
  ARBPinst(C,nr,nc,Cred);
  
  // Cred contains negative values -> positive
  if (n != m)
  {
    DT mn = std::numeric_limits<DT>::max();
    for (j = 0; j < m; j++)
      for (i = 0; i < n; i++)
      {
        const DT &cr = Cred[j*n+i];
        if (cr < mn) mn = cr;
      }
    if (mn < 0)
    {
        for (j = 0; j < m; j++)
            for (i = 0; i < n; i++)
            Cred[j*n+i] -= mn;
    }

  }

  hungarianLSAP(Cred,n,m,rho,u,v,varrho,init_type);
  
  // reconstruction

  if (n > m)
  {
    for (i = 0; i < n; i++) if (rho[i] == -1) rho[i] = m;
  }
  else if (n < m)
  {
    if (varrho)
      for (j = 0; j < m; j++)
        if (varrho[j] == -1) varrho[j] = n;
  }
  delete[] Cred;
}

// -----------------------------------------------------------
/** \brief Reduce a (n+1)x(m+1) LSAPE instance to an equivalent nxm LSAP instance for FBP
 *  \param[in] C nxm edit cost matrix
 *  \param[in] nr number of rows (last row correspond to the null element)
 *  \param[in] nc number of colums (last column correspond to the null element)
 *  \param[out] Cred (n-1)x(m-1) reduced cost matrix (must be previously allocated)
 */
template <typename DT, typename IT>
void FBPinst(const DT *C, const IT &nr, const IT &nc, DT *Cred)
{
  IT m = nc-1, n = nr-1, i, j;
  IT nm = (n+1)*m;
  const DT *ptC = C, *ptIns;
  DT *ptR = Cred;

  for (j = 0; j < m; j++, ++ptC)
  {
    ptIns = ptC + n;
    for (i = 0; i < n; i++, ++ptC, ++ptR)
    {
      *ptR = *ptC - C[m*nr+i] - *ptIns;
    }
  }
}

// -----------------------------------------------------------
template <typename DT, typename IT>
void FBP(const DT *C, const IT &nr, const IT &nc, IT *rho, DT *u = NULL, DT *v = NULL, IT *varrho = NULL, unsigned short init_type = 1)
{
  IT n = nr-1, m = nc-1, i, j;
  DT *Cred = new DT[n*m];
  FBPinst(C,nr,nc,Cred);

  // Cred contains negative values -> positive
  DT mn = std::numeric_limits<DT>::max();
  for (j = 0; j < m; j++)
    for (i = 0; i < n; i++)
    {
      const DT &cr = Cred[j*n+i];
      if (cr < mn) mn = cr;
    }
  if (mn < 0)
  {
    for (j = 0; j < m; j++)
      for (i = 0; i < n; i++)
	Cred[j*n+i] -= mn;
  }

  hungarianLSAP(Cred,n,m,rho,u,v,varrho,init_type);
  
  // reconstruction
  if (n < m) // no removal
  {
    if (varrho)
      for (i = 0; i < m; i++) if (varrho[i] == -1) varrho[i] = n;
  }
  else if (n > m) // no insertion
  {
    for (i = 0; i < n; i++) if (rho[i] == -1) rho[i] = m;
  }
  
  delete[] Cred;
}

// -----------------------------------------------------------
template <typename DT, typename IT>
void FBP0inst(const DT *C, const IT &nr, const IT &nc, DT *Cext, const IT &nrc)
{
  IT m = nc-1, n = nr-1, i, j, mn;
  const DT *ptC = C, *ptIns = NULL;
  DT *ptR = Cext, *ptRem = NULL;

  if (n < m)
  {
    mn = m-n;
    // substitutions
    for (j = 0; j < m; j++, ++ptC, ptR+=mn)
    {
      const DT &ins = C[j*nr+n];
      for (i = 0; i < n; i++, ++ptC, ++ptR) 
	*ptR = *ptC - ins - C[m*nr+i];
    }
    // insertions
    ptR = Cext+n;
    for (j = 0; j < m; j++, ptR+=n)
      for (i = 0; i < mn; i++, ptR++) *ptR = 0;
  }
  else 
    if (n > m)
    {
      mn = n-m;
      // substitutions
      for (j = 0; j < m; j++, ++ptC)
      {
	const DT &ins = C[j*nr+n];
	for (i = 0; i < n; i++, ++ptC, ++ptR)
	  *ptR = *ptC - ins - C[m*nr+i];
      }
      // removals
      ptRem = Cext+n*m;
      for (j = 0; j < mn; j++, ptRem+=n) std::memcpy(ptRem,0,sizeof(DT)*n);
    }
    else // only substitutions
    {
       for (j = 0; j < m; j++, ++ptC)
       {
	 const DT &ins = C[j*nr+n];
	 for (i = 0; i < n; i++, ++ptC, ++ptR)
	   *ptR = *ptC - ins - C[m*nr+i];
       }
    }
}

// -----------------------------------------------------------
template <typename DT, typename IT>
void FBP0(const DT *C, const IT &nr, const IT &nc, IT *rho, DT *u, DT *v, IT *varrho = NULL, unsigned short init_type = 1)
{
  IT n = nr-1, m = nc-1, i, j;
  IT mxnm = std::max(n,m);
  DT *Cext = new DT[mxnm*mxnm];

  FBP0inst(C,nr,nc,Cext,mxnm);
  
  // the instance is negative -> find min and translate
  DT mn = std::numeric_limits<DT>::max();
  for (j = 0; j < m; j++)
    for (i = 0; i < n; i++)
    {
      const DT &c = Cext[j*mxnm+i];
      if (c < mn) mn = c;
    }
  if (mn < 0)
  {
    DT *pt = Cext;
    for (i = 0; i < mxnm*mxnm; i++, ++pt) *pt -= mn;
  }
  
  if (n <=m) // no removal
  {
    IT *rhox = new IT[mxnm];
    DT *ux = new DT[mxnm];

    hungarianLSAP(Cext,mxnm,mxnm,rhox,ux,v,varrho,init_type);
  
    for (i = 0; i < n; i++) { rho[i] = rhox[i]; u[i] = ux[i]; }
    if (varrho) for (j = 0; j < m; j++) if (varrho[j]>=n) varrho[j]=n;
    
    delete[] rhox; delete[] ux;
  }
  else  // no insertion
  {
    IT *varrhox = NULL;
    if (varrho) varrhox = new IT[mxnm];
    DT *vx = new DT[mxnm];
    
    hungarianLSAP(Cext,mxnm,mxnm,rho,u,vx,varrhox,init_type);

    for (i = 0; i < n; i++) if (rho[i]>=m) rho[i] = m;
    if (varrho) for (j = 0; j < m; j++) varrho[j] = varrhox[j];
    for (j = 0; j < m; j++) v[j] = vx[j];

    if (varrho) delete[] varrhox; delete[] vx;
  }
  
  delete[] Cext;
}

// -----------------------------------------------------------
/** \brief Reduce a (n+1)x(m+1) LSAPE instance to an equivalent nxm LSAP instance for FBP
 *  \param[in] C nxm edit cost matrix
 *  \param[in] nr number of rows (last row correspond to the null element)
 *  \param[in] nc number of colums (last column correspond to the null element)
 *  \param[out] Cred (n-1)x(m-1) reduced cost matrix (must be previously allocated)
 */
template <typename DT, typename IT>
void SFBPinst(const DT *C, const IT &nr, const IT &nc, DT *Cext, const IT &nrc)
{
  IT m = nc-1, n = nr-1, i, j, mn;
  const DT *ptC = C, *ptIns = NULL;
  DT *ptR = Cext, *ptRem = NULL;

  if (n < m)
  {
    mn = m-n;
    // substitutions
    for (j = 0; j < m; j++, ++ptC, ptR+=mn)
      for (i = 0; i < n; i++, ++ptC, ++ptR) *ptR = *ptC;
    // insertions
    ptIns = C+n;
    ptR = Cext+n;
    for (j = 0; j < m; j++, ptIns+=nr, ptR+=n)
      for (i = 0; i < mn; i++, ptR++) *ptR = *ptIns;
  }
  else 
    if (n > m)
    {
      mn = n-m;
      // substitutions
      for (j = 0; j < m; j++, ++ptC)
        for (i = 0; i < n; i++, ++ptC, ++ptR) *ptR = *ptC;
      // removals
      ptRem = Cext+n*m;
      ptC = C+(n+1)*m;
      for (j = 0; j < mn; j++, ptRem+=n) std::memcpy(ptRem,ptC,sizeof(DT)*n);
    }
    else // only substitutions
    {
       for (j = 0; j < m; j++, ++ptC)
	 for (i = 0; i < n; i++, ++ptC, ++ptR) *ptR = *ptC;
    }  
}

// -----------------------------------------------------------
template <typename DT, typename IT>
void SFBP(const DT *C, const IT &nr, const IT &nc, IT *rho, DT *u, DT *v, IT *varrho = NULL, unsigned short init_type = 1)
{
  IT n = nr-1, m = nc-1;
  IT mxnm = std::max(n,m);
  DT *Cext = new DT[mxnm*mxnm];
  SFBPinst(C,nr,nc,Cext,mxnm);
  
  if (n <=m) // no removal
  {
    IT *rhox = new IT[mxnm];
    DT *ux = new DT[mxnm];

    hungarianLSAP(Cext,mxnm,mxnm,rhox,ux,v,varrho,init_type);
  
    for (int i = 0; i < n; i++) { rho[i] = rhox[i]; u[i] = ux[i]; }
    if (varrho) for (int j = 0; j < m; j++) if (varrho[j]>=n) varrho[j]=n;
    
    delete[] rhox; delete[] ux;
  }
  else  // no insertion
  {
    IT *varrhox = NULL;
    if (varrho) varrhox = new IT[mxnm];
    DT *vx = new DT[mxnm];
    
    hungarianLSAP(Cext,mxnm,mxnm,rho,u,vx,varrhox,init_type);

    for (int i = 0; i < n; i++) if (rho[i]>=m) rho[i] = m;
    if (varrho) for (int j = 0; j < m; j++) varrho[j] = varrhox[j];
    for (int j = 0; j < m; j++) v[j] = vx[j];

    if (varrho) delete[] varrhox; delete[] vx;
  }
  
  delete[] Cext;
}

// -----------------------------------------------------------
#endif
