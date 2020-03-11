// =========================================================================
/** \file greedy-lsape.hh
 *  \brief Greedy algorithms for approximating symmetric and asymmetric Linear Sum Assignment Problems with Error-Correction (LSAPE)
 * \author Sebastien Bougleux (Normandie Univ, CNRS - ENSICAEN - UNICAEN, GREYC, Caen, France)
 */
// =========================================================================
/* 
   This file is part of LSAPE.
   
   LSAPE is free software: you can redistribute it and/or modify
   it under the terms of the CeCILL-C License. See README for more
   details.
   
   You should have received a copy of the CeCILL-C License along with 
   LSAPE. If not, see http://www.cecill.info/index.en.html

   -----------------------------------------------------------
   
   Creation: December 5 2015
   Modifications: July 6 2017
   
   -----------------------------------------------------------
   Main functions:

   cost = greedyLSAPE(C,n,m,rho,varrho,greedy_type=2)
   greedykLSAPE(C,n,m,rho,varrho,k,costs,greedy_type=2)

   -> see at the end of this file or generate doxygen documentation
*/
// =========================================================================

#ifndef _GREEDY_LSAPE_
#define _GREEDY_LSAPE_

#include <iostream>
#include <limits>
#include <algorithm>
#include "utils.hh"

// -----------------------------------------------------------
// Basic greedy LSAP for n >= m only (not checked)
// return the cost of the approximate solution
// -----------------------------------------------------------
template <class DT, typename IT>
DT greedyBasicLSAPE(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho = NULL)
{
  IT i = -1, imin, itmp;
  DT cmin, mxdt = std::numeric_limits<DT>::max(), approx = 0;
  bool deletevarrho = false;
  if (varrho == NULL) { deletevarrho = true; varrho = new IT[m]; }

  IT *unass = new IT[n+1], *pti_unass = NULL;
  IT *pti_unass_beg = unass, *pti_unass_end = unass+n, *pti_min = NULL;
  
  for (i = 0; i < n; i++) { unass[i] = i; rho[i] = -1; }
  
  // augmentation of columns
  for (IT j = 0; j < m; j++)
  {
    // find the min among unassigned rows
    cmin = mxdt;
    for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++)
    {
      const DT &cij = C[j*n+*pti_unass];
      if (cij  < cmin) { cmin = cij; pti_min = pti_unass; }
    }
    // assign the row which provides the minimum and update the approximate solution
    imin = *pti_min; rho[imin] = j; varrho[j] = imin;
    *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
    approx += cmin;
  }
  
  delete[] unass;
  if (deletevarrho) delete[] varrho;
  return approx;
}

// -----------------------------------------------------------
// Refined greedy LSAP for n >= m only (not checked)
// return the cost of the approximate solution
// -----------------------------------------------------------
template <class DT, typename IT>
DT greedyRefinedLSAPE(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho = NULL)
{
  IT nass = 0, i = -1, j = -1, imin, jmin;
  DT cmin, ckmin, mxdt = std::numeric_limits<DT>::max(), approx = 0;
  bool deletevarrho = false;
  if (varrho == NULL) { deletevarrho = true; varrho = new IT[m]; }

  IT *unassi = new IT[n+1], *pti_unass = NULL;
  IT *pti_unass_beg = unassi, *pti_unass_end = unassi+n, *pti_min = NULL;

  IT *ptj_unass = NULL, *unassj = new IT[m+1];
  IT *ptj_unass_beg = unassj, *ptj_unass_end = unassj+m, *ptj_min = NULL;

  IT *ptj_unass1 = NULL;
  
  for (i = 0; i < n; i++) { unassi[i] = i; rho[i] = -1; }
  for (j = 0; j < m; j++) { unassj[j] = j; }
  
  // augmentation of columns
  for (ptj_unass1 = ptj_unass_beg; ptj_unass1 != ptj_unass_end;)
  {
    j = *ptj_unass1;
    // find the min among unassigned rows
    cmin = mxdt;
    for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++)
    {
      const DT &cij = C[j*n+*pti_unass];
      if (cij  < cmin) { cmin = cij; pti_min = pti_unass; }
    }
    // find the min among unassigned columns for imin
    imin = *pti_min;
    ckmin = mxdt;
    for (ptj_unass = ptj_unass_beg; ptj_unass != ptj_unass_end; ptj_unass++)
    {
      const DT &cik = C[*ptj_unass*n+imin];
      if (cik  < ckmin) { ckmin = cik; ptj_min = ptj_unass; }
    }
    // assign the row and column which provides the minimum
    if (cmin <= ckmin)
    {
      rho[imin] = j; varrho[j] = imin;
      *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
      ptj_unass_beg++; 
      approx += cmin;
    }
    else
    {
      jmin = *ptj_min; rho[imin] = jmin; varrho[jmin] = imin;
      *ptj_min = *ptj_unass_beg; *ptj_unass_beg = jmin; ptj_unass_beg++;
      *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
      approx += ckmin;
    }
    ptj_unass1 = ptj_unass_beg;
  }
  
  delete[] unassi; delete[] unassj;
  if (deletevarrho) delete[] varrho;
  return approx;
}

// -----------------------------------------------------------
// Loss greedy LSAP for n >= m only (not checked)
// return the cost of the approximate solution
// -----------------------------------------------------------
template <class DT, typename IT>
DT greedyLossLSAPE(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho = NULL)
{
  IT nass = 0, i = -1, j = -1, imin, jmin, imin2, imin3;
  DT cmin, cij, ckmin, mxdt = std::numeric_limits<DT>::max(), approx = 0, cmin2, cmin3;
  bool deletevarrho = false;
  if (varrho == NULL) { deletevarrho = true; varrho = new IT[m]; }

  IT *unassi = new IT[n+1], *pti_unass = NULL;
  IT *pti_unass_beg = unassi, *pti_unass_end = unassi+n, *pti_min = NULL, *pti_min2 = NULL, *pti_min3 = NULL;

  IT *ptj_unass = NULL, *unassj = new IT[m+1];
  IT *ptj_unass_beg = unassj, *ptj_unass_end = unassj+m, *ptj_min = NULL;

  IT *ptj_unass1 = NULL;
  
  for (i = 0; i < n; i++) { unassi[i] = i; rho[i] = -1; }
  for (j = 0; j < m; j++) { unassj[j] = j; }
  
  // augmentation of columns
  for (ptj_unass1 = ptj_unass_beg; ptj_unass1 != ptj_unass_end;)
  {
    j = *ptj_unass1;
    // find the min among unassigned rows
    cmin = mxdt;
    for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++)
    {
      cij = C[j*n+*pti_unass];
      if (cij  < cmin) { cmin = cij; pti_min = pti_unass; }
    }
    // find the min among unassigned columns for imin
    imin = *pti_min;
    ckmin = mxdt;
    for (ptj_unass = ptj_unass_beg; ptj_unass != ptj_unass_end; ptj_unass++)
    {
      const DT &cik = C[*ptj_unass*n+imin];
      if (cik  < ckmin) { ckmin = cik; ptj_min = ptj_unass; }
    }
    // assign the row and column which provides the minimum
    jmin = *ptj_min;
    if (j == jmin)
    {
      rho[imin] = j; varrho[j] = imin;
      *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
      ptj_unass_beg++;
      approx += cmin;
    }
    else
    {
      // find the min among unassigned rows different to imin, for j => find the 2nd min
      cmin3 = cmin2 = mxdt;
      for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++)
      {
	if (*pti_unass != imin)
	{
	  cij = C[j*n+*pti_unass];
	  if (cij  < cmin2) { cmin2 = cij; pti_min2 = pti_unass; }
	  cij = C[jmin*n+*pti_unass];
	  if (cij  < cmin3) { cmin3 = cij; pti_min3 = pti_unass; }
	}
      }
      imin2 = *pti_min2;
      imin3 = *pti_min3;
      if (cmin + cmin3 < cmin2 + ckmin) // remove j, jmin, imin, imin3
      {
	rho[imin] = j; varrho[j] = imin;
	rho[imin3] = jmin; varrho[jmin] = imin3;
	*pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
	ptj_unass_beg++;
	*ptj_min = *ptj_unass_beg; *ptj_unass_beg = jmin; ptj_unass_beg++;
	*pti_min3 = *pti_unass_beg; *pti_unass_beg = imin3; pti_unass_beg++;
	approx += cmin + cmin3;
      }
      else // remove j, jmin, imin, imin2
      {
	rho[imin2] = j; varrho[j] = imin2;
	rho[imin] = jmin; varrho[jmin] = imin;
	ptj_unass_beg++;
	*ptj_min = *ptj_unass_beg; *ptj_unass_beg = jmin; ptj_unass_beg++;
	*pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
	*pti_min2 = *pti_unass_beg; *pti_unass_beg = imin2; pti_unass_beg++;
	approx += ckmin + cmin2;
      }
    }
    ptj_unass1 = ptj_unass_beg;
  }
  
  delete[] unassi; delete[] unassj;
  if (deletevarrho) delete[] varrho;
  return approx;
}

// -----------------------------------------------------------
// BasicSort greedy LSAP for n >= m only (not checked)
// return the cost of the approximate solution
// -----------------------------------------------------------
template <class DT, typename IT>
class CostSortComp
{
private:
  const DT *_C;
public:
  CostSortComp(const DT *C) : _C(C) { }
  ~CostSortComp() { }
  inline bool operator()(const IT &ij, const IT &kl) { return (_C[ij] < _C[kl]); }
};

template <class DT, typename IT>
DT greedyBasicSortLSAPE(const DT *C, const IT &n, const IT &m, IT *rho, IT *varrho = NULL)
{ 
  // sort the costs
  IT *idxs = new IT[n*m+1];
  IT *pt_idxs_end = idxs+n*m, i = 0, j;
  DT approx = 0;
  for (IT *pt_idx = idxs; pt_idx != pt_idxs_end; pt_idx++, i++) *pt_idx = i;
  
  CostSortComp<DT,IT> comp(C);
  std::sort(idxs,idxs+n*m,comp);

  // assign element in ascending order of the costs
  bool deletevarrho = false;
  if (varrho == NULL) { deletevarrho = true; varrho = new IT[m]; }
  for (IT i = 0; i < n; i++) rho[i] = -1;
  for (IT j = 0; j < m; j++) varrho[j] = -1;
  for (IT *pt_idx = idxs, nbe = 0; pt_idx != pt_idxs_end && nbe < m; pt_idx++)
  {
    ind2sub(*pt_idx,n,i,j);
    if (rho[i] == -1 && varrho[j] == -1)
    {
      rho[i] = j;
      varrho[j] = i;
      approx += C[*pt_idx];
      nbe++;
    }
  }
  if (deletevarrho) delete[] varrho;
  delete[] idxs;
  return approx;
}

// -----------------------------------------------------------
// Counting sort greedy LSAP for n >= m only (not checked)
// and for integer values
// return the cost of the approximate solution
// -----------------------------------------------------------
template <class DT, typename IT>
DT greedyCountingSortLSAPE(const DT *C, const IT &nrows, const IT &ncols, IT *rho, IT *varrho)
{
  //  IT n = nrows-1, m = ncols-1;
  // find min and max values
  DT minc = C[0], maxc = C[0];
  const DT *ite = C+nrows*ncols-1;
  for (const DT *itc = C+1; itc < ite; itc++)
  {
    const DT &vc = *itc;
    if (vc < minc) minc = vc;
    else if (vc > maxc) maxc = vc;
  }

  // construct histogram
  IT nbins = maxc - minc + 1;
  IT *bins = new IT[nbins+1];
  const IT *itbe = bins+nbins;
  for (IT *itb = bins; itb < itbe; itb++) *itb = 0;
  for (const DT *itc = C; itc < ite; itc++) bins[*itc-minc]++;

  // starting index for each cost value
  IT tot = 0, oldcount;
  for (IT i = 0; i < nbins; i++) { oldcount = bins[i]; bins[i] = tot; tot += oldcount; }

  // reoder the costs, preserving order of C with equal keys and favor substitutions
  IT *idxs = new IT[nrows*ncols], k = 0;
  for (const DT *itc = C; itc < ite; itc++, k++) { idxs[bins[*itc-minc]] = k; bins[*itc-minc]++; }

  // assign element in ascending order of the costs
  IT *pt_idxs_end = idxs+nrows*ncols-1, i = 0, j = 0, approx = 0;
  for (; i < nrows-1; i++) rho[i] = -1;
  for (; j < ncols-1; j++) varrho[j] = -1;
  for (IT *pt_idx = idxs, nbc = 0; pt_idx != pt_idxs_end && nbc < nrows+ncols-2; pt_idx++)
  {
    ind2sub(*pt_idx,nrows,i,j);
    if (i == nrows-1) { if (varrho[j] == -1) { varrho[j] = i; approx += C[*pt_idx]; nbc++; } }
    else
      if (rho[i] == -1)
	if (j == ncols-1) { rho[i] = j; approx += C[*pt_idx]; nbc++; }
	else if (varrho[j] == -1) { rho[i] = j; varrho[j] = i; approx += C[*pt_idx]; nbc += 2; }
  }
  delete[] idxs; delete[] bins;
  return approx;
}

template<> float greedyCountingSortLSAPE<>(const float *C, const int &n, const int &m, int *rho, int *varrho) { std::cerr << "case not supported => return MAX cost\n"; return std::numeric_limits<float>::min(); }
template<> double greedyCountingSortLSAPE<>(const double *C, const int &n, const int &m, int *rho, int *varrho) { std::cerr << "case not supported => return MAX cost\n"; return std::numeric_limits<double>::min(); }

// -----------------------------------------------------------
// Counting sort greedy LSAP for n >= m only (not checked)
// and for integer values
// return the cost of the approximate solution
// -----------------------------------------------------------
template <class DT, typename IT>
void greedyCountingSortLSAPE(const DT *C, const IT &nrows, const IT &ncols, IT *rho, IT *varrho, const IT &kcheap, DT *approx)
{
  // find min and max values
  DT minc = C[0], maxc = C[0];
  const DT *ite = C+nrows*ncols-1;
  for (const DT *itc = C+1; itc < ite; itc++)
  {
    const DT &vc = *itc;
    if (vc < minc) minc = vc;
    else if (vc > maxc) maxc = vc;
  }
  
  // construct histogram
  IT nbins = maxc - minc + 1;
  IT *bins = new IT[nbins+1];
  const IT *itbe = bins+nbins;
  for (IT *itb = bins; itb < itbe; itb++) *itb = 0;
  for (const DT *itc = C; itc < ite; itc++) bins[*itc-minc]++;

  // starting index for each cost value
  IT tot = 0, oldcount;
  for (IT i = 0; i < nbins; i++) { oldcount = bins[i]; bins[i] = tot; tot += oldcount; }

  // reoder the costs, preserving order of C with equal keys
  IT *idxs = new IT[nrows*ncols], k = 0;
  for (const DT *itc = C; itc < ite; itc++, k++) { idxs[bins[*itc-minc]] = k; bins[*itc-minc]++; }
  
  // assign element in ascending order of the costs
  IT *pt_idxs_end = idxs+nrows*ncols-1, i, j;
  IT *pt_idx0 = idxs;
  for (IT ncheap = 0; ncheap < kcheap && pt_idx0 != pt_idxs_end; ncheap++, pt_idx0++)
  {
    //approx[ncheap] = 0;
    for (i = 0; i < nrows-1; i++) rho[ncheap*(nrows-1)+i] = -1;
    for (j = 0; j < ncols-1; j++) varrho[j*kcheap+ncheap] = -1;
    for (IT *pt_idx = pt_idx0, nbc = 0; pt_idx != pt_idxs_end && nbc < nrows+ncols-2; pt_idx++)
    {
      ind2sub(*pt_idx,nrows,i,j);
      if (i == nrows-1) { if (varrho[j*kcheap+ncheap] == -1) { varrho[j*kcheap+ncheap] = i; approx[ncheap] += C[*pt_idx]; nbc++; } }
      else
	if (rho[ncheap*(nrows-1)+i] == -1)
	  if (j == ncols-1) { rho[ncheap*(nrows-1)+i] = j; approx[ncheap] += C[*pt_idx]; nbc++; }
	  else if (varrho[j*kcheap+ncheap] == -1) { rho[ncheap*(nrows-1)+i] = j; varrho[j*kcheap+ncheap] = i; approx[ncheap] += C[*pt_idx]; nbc += 2; }
    }
  }
  delete[] idxs; delete[] bins;
}

template<> void greedyCountingSortLSAPE<>(const float *C, const int &n, const int &m, int *rho, int *varrho, const int &k, float *approx) { std::cerr << "case not supported\n"; }
template<> void greedyCountingSortLSAPE<>(const double *C, const int &n, const int &m, int *rho, int *varrho, const int &k, double *approx) { std::cerr << "case not supported\n"; }

// --------------------------------------------------------------------------------
// Main function: greedy algorithms for both square and rectangular cost matrices
/**
 * \brief Compute an approximate solution to the LSAPE as an assignment with low costs, with greedy algorithms given a square or a rectangular cost matrix
 * \param[in] C nxm cost matrix, represented as an array of size \c nm obtained by concatenating its columns (last row and column correspond to null elements)
 * \param[in] n Number of rows of \p C (size of the 1st set)
 * \param[in] m Number of columns of \p C (size of the 2nd set)
 * \param[out] rho Array of size n-1 (must be previously allocated) for the assignment of the rows to the columns, rho[i]=n-1 indicates that i is assigned to the null element, else rho[i]=j indicates that i is assigned to j
 * \param[out] varrho Array of size m-1 (must be previously allocated) for the assignement of the columns to the rows
 * \param[in] greedy_type 0:Basic, 1:Refined, 2: Loss (default), 3: Basic sort, 4: Counting sort (integers only)
 * \return Cost of the assignment
 *
 * \note Adapted from the reference\n
 * Approximate Graph Edit Distance in Quadratic Time. K. Riesen, M. Ferrer, H. Bunke. ACM Transactions on Computational Biology and Bioinformatics, 2015.
 */
// --------------------------------------------------------------------------------
template <class DT, typename IT>
DT greedyLSAPE(const DT *C, const IT &nrows, const IT &ncols, IT *rho, IT *varrho, unsigned short greedy_type = 2)
{
  switch(greedy_type)
  {
  case 0: return greedyBasicLSAPE<DT,IT>(C,nrows,ncols,rho,varrho);
  case 1: return greedyRefinedLSAPE<DT,IT>(C,nrows,ncols,rho,varrho);
  case 3: return greedyBasicSortLSAPE<DT,IT>(C,nrows,ncols,rho,varrho);
  case 4: return greedyCountingSortLSAPE<DT,IT>(C,nrows,ncols,rho,varrho);
  default: return greedyLossLSAPE<DT,IT>(C,nrows,ncols,rho,varrho);
  }
}

// --------------------------------------------------------------------------------
/**
 * \brief Compute an approximate solution to the LSAPE as an assignment with low costs, with greedy algorithms given a square or a rectangular cost matrix
 * \param[in] C nxm cost matrix, represented as an array of size \c nm obtained by concatenating its columns (last row and column correspond to null elements)
 * \param[in] n Number of rows of \p C (size of the 1st set)
 * \param[in] m Number of columns of \p C (size of the 2nd set)
 * \param[out] rho Matrix of size (n-1)xk (must be previously allocated) for the assignment of the rows to the columns, rho[k*(n-1)+i]=m-1 indicates that i is assigned to the null element in the k-st assignment, else rho[k*(n-1)+i]=j indicates that i is assigned to j in the k-st assignment
 * \param[out] varrho Matrix of size kx(m-1) (must be previously allocated) for the assignement of the columns to the rows, varrho[j*kass+k]=n-1 indicates that j is assigned to the null element in the k-st assignment, else varrho[j*kass+k]=i indicates that j is assigned to i in the k-st assignment
 * \param[in] kass Number of requested assignments
 * \param[out] costs Array of size kass so that costs[k] returns the cost of the k-st assignment (must be previously allocated)
 * \param[in] greedy_type 0:Basic, 1:Refined, 2: Loss (default), 3: Basic sort (std::sort), 4: Counting sort (integers only)
 */
template <class DT, typename IT>
void greedykLSAPE(const DT *C, const IT &nrows, const IT &ncols, IT *rho, IT *varrho, const IT &kass, DT* costs, unsigned short greedy_type = 2)
{
  //  switch(greedy_type)
  //{
  //case 0: return greedyBasicLSAPE<DT,IT>(C,nrows,ncols,rho,varrho);
  //case 1: return greedyRefinedLSAPE<DT,IT>(C,nrows,ncols,rho,varrho);
  //case 3: return greedyBasicSortLSAPE<DT,IT>(C,nrows,ncols,rho,varrho);
  //case 4: 
  greedyCountingSortLSAPE<DT,IT>(C,nrows,ncols,rho,varrho,kass,costs);
  //default: return greedyLossLSAPE<DT,IT>(C,nrows,ncols,rho,varrho);
  //}
}

// -----------------------------------------------------------
#endif
