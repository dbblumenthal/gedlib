// -----------------------------------------------------------   
/** \file lsap.hh
 *  \brief Solvers for the Linear Sum Assignment Problem (LSAP) and its dual problem, for balanced and unbalanced instances.
 *  \author Sebastien Bougleux (Normandie Univ, UNICAEN, ENSICAEN, CNRS, GREYC, Image team, Caen, France)
*/
/* -----------------------------------------------------------
   This file is part of LSAPE.
   
   LSAPE is free software: you can redistribute it and/or modify
   it under the terms of the CeCILL-C License. See README for more
   details.

   -----------------------------------------------------------   
   Creation:
   Last modif:
*/

#ifndef __LSAP_HH__
#define __LSAP_HH__

#ifndef sub2idx
#define sub2idx(i,j,n1) (j*n1+i)
#endif

#include <typeinfo>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <utils.hh>
#include <matchings.hh>
#include <graphs.hh>

namespace liblsap {

  template <typename DataType = int, typename IndexType = int,
	    class PrimContainer = NodeMaps<IndexType> >
  class LSAP {

  public:
        
  protected:
    IndexType _n1;
    IndexType _n2;
    DataType *_C;
    bool _C_internal;
    PrimContainer *_primals;
    bool _prims_internal;
    DataType *_dual[2];
    bool _dual_internal;
    static constexpr DataType _zero = 0;
    IndexType _lb_nb_prim;
    IndexType _ub_nb_prim;
    
    // ==============================================================
    // PRIMAL-DUAL HUNGARIAN (AUGMENTING PATHS)
    // ==============================================================
    // --------------------------------------------------------------
    // rows
    IndexType basicInitS1(Matching<IndexType> &prim)
    {
      IndexType i = 0, j, nb_match = 0;
      DataType mn;
      
      // find the min of each row
      for (; i < _n1; i++)
      {
	mn = std::numeric_limits<DataType>::max();
	for (j = 0; j < _n2; j++)
	{
	  const DataType &c = _C[sub2idx(i,j,_n1)];
	  if (c < mn) mn = c;
	}
	dual(0,i) = mn;
	prim(0,i) = _n2;
      }
      
      for (j = 0; j < _n2; j++) { dual(1,j) = 0; prim(1,j) = _n1; }
      
      // assign
      for (i = 0; i < _n1; i++)
	if (prim(0,i) == _n2) {
	  for (j = 0; j < _n2; j++)
	    if (prim(1,j) == _n1 && _C[sub2idx(i,j,_n1)] == dual(0,i))
	      { prim(0,i) = j; prim(1,j) = i; nb_match++; break; }
	}
      return nb_match;
    }
    // --------------------------------------------------------------
    // columns
    IndexType basicInitS2(Matching<IndexType> &prim)
    {
      IndexType i = 0, j = 0, nb_match = 0;
      DataType mn;
      
      // find the min of each row
      for (; j < _n2; j++)
      {
	mn = std::numeric_limits<DataType>::max();
	for (i = 0; i < _n1; i++)
	{
	  const DataType &c = _C[sub2idx(i,j,_n1)];
	  if (c < mn) mn = c;
	}
	dual(1,j) = mn;
	prim(1,j) = _n1;
      }
      
      for (i = 0; i < _n1; i++) { dual(0,i) = 0; prim(0,i) = _n2; }
      
      // assign
      for (j = 0; j < _n2; j++)
      {
	if (prim(1,j) == _n1) {
	  for (i = 0; i < _n1; i++) {
	    if (prim(0,i) == _n2 && _C[sub2idx(i,j,_n1)] == dual(1,j))
	      { prim(0,i) = j; prim(1,j) = i; nb_match++; break; }
	  }
	}
      }
      return nb_match;
    }
    // --------------------------------------------------------------
    IndexType basicInit(Matching<IndexType> &prim)
    {
      IndexType i = 0, j, nb_match = 0;
      DataType mn, val, c;
      
      // find the min of each row
      for (; i < _n1; i++)
      {
	mn = std::numeric_limits<DataType>::max();
	for (j = 0; j < _n2; j++)
	{
	  c = _C[sub2idx(i,j,_n1)];
	  if (c < mn) mn = c;
	}
	dual(0,i) = mn;
	prim(0,i) = _n1;
      }
  
      // find the min of each column
      for (j = 0; j < _n2; j++)
      {
	mn = std::numeric_limits<DataType>::max();
	for (i = 0; i < _n1; i++)
	{
	  val = _C[sub2idx(i,j,_n1)] - dual(0,i);
	  if (val < mn) mn = val;
	}
	dual(1,j) = mn;
	prim(1,j) = _n1;
      }
  
      // assign
      for (i = 0; i < _n1; i++)
      {
	if (prim(0,i) == _n2)
	  for (j = 0; j < _n2; j++)
	    if (prim(1,j) == _n1 && _C[sub2idx(i,j,_n1)] - dual(0,i) == dual(1,j))
	      { prim(0,i) = j; prim(1,j) = i; nb_match++; break; }
      }
      return nb_match;
    }
    // --------------------------------------------------------------
    IndexType initSquare(Matching<IndexType> &prim, unsigned short init_type = 1)
    {
      IndexType nb_assign = 0;
      switch (init_type) {
      case 0:
	{
	  for (IndexType i = 0; i < _n1; i++)
	    { dual(0,i) = dual(1,i) = 0; prim(0,i) = prim(1,i) = _n1; }
	}
	break;
      default:
	nb_assign = basicInit(prim);
      }
      return nb_assign;
    }
    // --------------------------------------------------------------
    IndexType initRect(Matching<IndexType> &prim, unsigned short init_type = 1)
    {
      IndexType nb_assign = 0;
      switch (init_type) {
      case 0:
	{
	  for (IndexType i = 0; i < _n1; i++) { dual(0,i) = 0; prim(0,i) = _n2; }
	  for (IndexType j = 0; j < _n2; j++) { dual(1,j) = 0; prim(1,j) = _n1; }
	}
	break;
      default:
	nb_assign = (_n1 < _n2 ? basicInitS1(prim) : basicInitS2(prim));
      }
      return nb_assign;
    }
    // --------------------------------------------------------------
    void augmentingPathS2(const IndexType &k, Matching<IndexType> &prim, 
			  IndexType *U, IndexType *SV, IndexType *pred, DataType *pi,
			  IndexType &zi)
    {
      IndexType i = 0, j = k, r = 0, *SVptr = SV, *ulutop = U, *uluptr = NULL;
      DataType delta = 0, mx = std::numeric_limits<DataType>::max(), cred = 0;
      const IndexType *svptr = NULL, *luptr = NULL, *uend = U+_n1, *lusutop = U;

      for (i = 0; i < _n1; i++) { pi[i] = mx; U[i] = i; }
      
      while (true)
      {
	*SVptr = j; *(++SVptr) = _n2;
	for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
	{
	  r = *uluptr;
	  //cred = std::abs(_C[sub2idx(r,j,_n1)] - (dual(0,r) + dual(1,j)));
	  //if (isZero(cred) && cred != zero) std::cout << "cred=" << cred << std::endl;
	  cred = _C[sub2idx(r,j,_n1)] - (dual(0,r) + dual(1,j));
	  if (cred < pi[r])
	  {
	    //	    if (isZero(cred)) cred = zero;
	    pred[r] = j;
	    pi[r] = cred;
	    if (cred == _zero)
	    {
	      if (prim(0,r) == _n2) { zi = r; return; }
	      i = *ulutop; *ulutop = r; *uluptr = i; ++ulutop;
	    }
	  }
	}
	
	if (lusutop == ulutop) // dual update
	{
	  delta = mx;
	  for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
	    if (pi[*uluptr] < delta) delta = pi[*uluptr];
	  for (svptr = SV; *svptr != _n2; ++svptr) _dual[1][*svptr] += delta;
	  for (luptr = U; luptr != ulutop; ++luptr) _dual[0][*luptr] -= delta;
	  for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
	  {
	    pi[*uluptr] -= delta;
	    if (pi[*uluptr] == _zero) 
	    {
	      if (prim(0,*uluptr) == _n2) { zi = *uluptr; return; }
	      r = *ulutop; *ulutop = *uluptr; *uluptr = r; ++ulutop;
	    }
	  }
	} // end dual update
	i = *lusutop; ++lusutop; // i is now in SU
	j = prim(0,i);
      }
    }
    // -----------------------------------------------------------
    void augmentingPathS1(const IndexType &k, Matching<IndexType> &prim,
			  IndexType *V, IndexType *SU, IndexType *pred,
			  DataType *pi, IndexType &zj)
    {
      const IndexType *suptr = NULL, *lvptr = NULL, *vend = V+_n2, *lvsvtop = V;
      IndexType i = k, j = 0, c = 0, *SUptr = SU, *vlvtop = V, *vlvptr = NULL;
      DataType delta = 0, mx = std::numeric_limits<DataType>::max(), cred = 0;
      //*SU = -1;
      //zj = -1;
      
      for (j = 0; j < _n2; j++) { pi[j] = mx; V[j] = j; }
      
      while (true)
      {
	*SUptr = i; *(++SUptr) = _n1;
	for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // U\LU
	{
	  c = *vlvptr;
	  cred =  _C[sub2idx(i,c,_n1)] - (dual(0,i) + dual(1,c));
	  if (cred < pi[c])
	  {
	    pred[c] = i;
	    pi[c] = cred;
	    if (cred == 0)
	    {
	      if (prim(1,c) == _n1) { zj = c; return; }
	      j = *vlvtop; *vlvtop = c; *vlvptr = j; ++vlvtop;
	    }
	  }
	}
        
	if (lvsvtop == vlvtop) // dual update
	{
	  delta = mx;
	  for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // V\LV
	    if (pi[*vlvptr] < delta) delta = pi[*vlvptr];
	  for (suptr = SU; *suptr != _n1; ++suptr) dual(0,*suptr) += delta;
	  for (lvptr = V; lvptr != vlvtop; ++lvptr) dual(1,*lvptr) -= delta;
	  for (vlvptr = vlvtop;  vlvptr != vend; ++vlvptr) // V\LV
	  {
	    pi[*vlvptr] -= delta;
	    if (pi[*vlvptr] == 0)
	    {
	      if (prim(1,*vlvptr) == _n1) { zj = *vlvptr; return; }
	      c = *vlvtop; *vlvtop = *vlvptr; *vlvptr = c; ++vlvtop;
	    }
	  }
	} // end dual update
	j = *lvsvtop; ++lvsvtop; // j is now in SV
	i = prim(1,j);
      }
    }
    // --------------------------------------------------------------
    /*void augmentingPathS2_C(const IndexType &k, const Matching<IndexType> &prim, 
			    IndexType *U, IndexType *SV, IndexType *pred, DataType *pi,
			    IndexType &zi)
    {
      IndexType i = 0, j = k, r = 0, *SVptr = SV, *ulutop = U, *uluptr = NULL;
      DataType delta = 0, mx = std::numeric_limits<DataType>::max(), cred = 0;
      const IndexType *svptr = NULL, *luptr = NULL, *uend = U+_n1, *lusutop = U;

      for (i = 0; i < _n1; i++) { pi[i] = mx; U[i] = i; }
      
      while (true)
      {
	*SVptr = j; *(++SVptr) = _n2;
	for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
	{
	  r = *uluptr;
	  cred = _Cred[sub2idx(r,j,_n1)];
	  if (cred < pi[r])
	  {
	    pred[r] = j;
	    pi[r] = cred;
	    if (cred == _zero)
	    {
	      if (prim(0,r) == _n2) { zi = r; return; }
	      i = *ulutop; *ulutop = r; *uluptr = i; ++ulutop;
	    }
	  }
	}	
	if (lusutop == ulutop) // dual update
	{
	  delta = mx;
	  for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
	    if (pi[*uluptr] < delta) delta = pi[*uluptr];
	  for (svptr = SV; *svptr != _n2; ++svptr) {
	    for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU 
	      Cred[sub2idx(*uluptr,*svptr,_n1)] -= delta;
	    _dual[1][*svptr] += delta;
	  }
	  for (luptr = U; luptr != ulutop; ++luptr) {
	  // TODO loop on V\LV ie vlv (must do the same as for ulu and u)
	    _dual[0][*luptr] -= delta;
	  }
	  for (uluptr = ulutop;  uluptr != uend; ++uluptr) // U\LU
	  {
	    pi[*uluptr] -= delta;
	    if (pi[*uluptr] == _zero) 
	    {
	      if (prim(0,*uluptr) == _n2) { zi = *uluptr; return; }
	      r = *ulutop; *ulutop = *uluptr; *uluptr = r; ++ulutop;
	    }
	  }
	} // end dual update
	i = *lusutop; ++lusutop; // i is now in SU
	j = prim(0,i);
      }
      }*/
    // --------------------------------------------------------------
    void hungarianSquare(Matching<IndexType> &prim, unsigned short init_type = 1)
    {
      IndexType i, j, r, k, *SV = NULL, *pred = NULL, *U = NULL;
      DataType *pi = NULL;
      IndexType nass = initSquare(prim,init_type);
      
      // augment columns
      if (nass < _n1)
      {
	U = new IndexType[_n1+1]; SV = new IndexType[_n1+1];
	pi = new DataType[_n1+1]; pred = new IndexType[_n1];
	for (k = 0; k < _n1; k++)
	  if (prim(1,k) == _n1)
	  {
	    augmentingPathS2(k,prim,U,SV,pred,pi,i);
	    for (j = _n1; j != k;)  // update primal solution = new partial assignment
	    {
	      j = pred[i]; prim(0,i) = j;
	      r = prim(1,j); prim(1,j) = i; i = r;
	    }
	    if (!isPrimalDual(prim)) std::cout << "not prim-dual !!!\n";
	  }
	delete[] U; delete[] pred; delete[] pi; delete[] SV;
      }
    }
    // --------------------------------------------------------------
    void hungarianRect(Matching<IndexType> &prim, unsigned short init_type = 1)
    {
      IndexType nass = 0, mass = 0, i, j, r, c, k, *S = NULL, *U = NULL, *V = NULL, *pred = NULL;
      DataType *pi = NULL;
      
      if (_n1 < _n2) // augmentation of rows
      {
	nass = initRect(prim,init_type);
	if (nass < _n1)
	{
	  V = new IndexType[_n2+1]; S = new IndexType[_n1+1];
	  pi = new DataType[_n2]; pred = new IndexType[_n2];
	  for (k = 0; k < _n1; k++)
	    if (prim(0,k) == _n2)
	    {
	      augmentingPathS1(k,prim,V,S,pred,pi,j);
	      for (i = _n1; i != k;)  // update primal solution = new partial assignment
	      {
		i = pred[j]; prim(1,j) = i;
		c = prim(0,i); prim(0,i) = j; j = c;
	      }
	      nass++;
	    }
	  delete[] V; delete[] pred; delete[] pi; delete[] S;
	}
      }
      else // nrows > ncols, augmentation of columns
      {
	mass = initRect(prim,init_type);
	if (mass < _n2)
	{
	  U = new IndexType[_n1+1]; S = new IndexType[_n2+1];
	  pi = new DataType[_n1]; pred = new IndexType[_n1];
	  for (k = 0; k < _n2; k++)
	    if (prim(1,k) == _n1)
	    {
	      augmentingPathS2(k,prim,U,S,pred,pi,i);
	      for (j = _n2; j != k;)  // update primal solution = new partial assignment
	      {
		j = pred[i]; prim(0,i) = j;
		r = prim(1,j); prim(1,j) = i; i = r;
	      }
	      mass++;
	    }
	  delete[] U; delete[] pred; delete[] pi; delete[] S;
	}
      }
    }
    // --------------------------------------------------------------
    void hungarian(Matching<IndexType> &prim, unsigned short init_type = 1)
    { (_n1 == _n2 ? hungarianSquare(prim,init_type) : hungarianRect(prim,init_type)); }
    
    // ==============================================================
    // GREEDY METHODS FOR APPROXIMATION
    // ==============================================================
    // --------------------------------------------------------------
    DataType greedyBasic(IndexType *prim12, IndexType *prim21 = NULL,
			 IndexType *perm1 = NULL, IndexType *perm2 = NULL)
    {
      IndexType i, imin;
      DataType cmin, approx = _zero;
      const DataType mxdt = std::numeric_limits<DataType>::max();
      bool delete21 = false;
      if (prim21 == NULL) { delete21 = true; prim21 = new IndexType[_n2]; }
      IndexType nmx = std::max(_n1,_n2);
      IndexType *unass = new IndexType[nmx+1], *pti_unass = NULL;
      IndexType *pti_unass_beg = unass, *pti_min = NULL;
      const IndexType *pti_unass_end = unass+nmx;
  
      if (_n1 > _n2) // assign columns
      {
	if (perm1 == NULL) for (i = 0; i < _n1; i++) { unass[i] = i; prim12[i] = _n2; }
	else for (i = 0; i < _n1; i++) { unass[i] = perm1[i]; prim12[i] = _n2; }
    
	if (perm2 == NULL)
	{
	  for (IndexType j = 0; j < _n2; j++)
	  {
	    // find the min among unassigned rows
	    cmin = mxdt;
	    for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++)
	    {
	      const DataType &cij = _C[sub2idx(*pti_unass,j,_n1)];
	      if (cij  < cmin) { cmin = cij; pti_min = pti_unass; }
	    }
	    // assign the row which provides the minimum and update the approximate solution
	    imin = *pti_min; prim12[imin] = j; prim21[j] = imin;
	    *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
	    approx += cmin;
	  }
	}
	else // a permutation of the rows is given
	{
	  for (IndexType jj = 0, *permit = perm2, j = 0; jj < _n2; jj++, ++permit)
	  {
	    j = *permit;
	    // find the min among unassigned rows
	    cmin = mxdt;
	    for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++)
	    {
	      const DataType &cij = _C[sub2idx(*pti_unass,j,_n1)];
	      if (cij  < cmin) { cmin = cij; pti_min = pti_unass; }
	    }
	    // assign the row which provides the minimum and update the approximate solution
	    imin = *pti_min; prim12[imin] = j; prim21[j] = imin;
	    *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
	    approx += cmin;
	  }
	}
      }
      else // assign rows for _n1 < _n2
      {
	if (perm2 == NULL) for (i = 0; i < _n1; i++) { unass[i] = i; prim21[i] = _n1; }
	else for (i = 0; i < _n2; i++) { unass[i] = perm2[i]; prim21[i] = _n1; }
	
	if (perm1 == NULL)
	{  
	  for (i = 0; i < _n1; i++)
	  {
	    // find the min among unassigned columns
	    cmin = mxdt;
	    for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++)
	    {
	      const DataType &cij = _C[sub2idx(i,*pti_unass,_n1)];
	      if (cij  < cmin) { cmin = cij; pti_min = pti_unass; }
	    }
	    // assign the column which provides the minimum and update the approximate solution
	    imin = *pti_min; prim21[imin] = i; prim12[i] = imin;
	    *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
	    approx += cmin;
	  }
	}
	else // a permutation of the rows is given
	{
	  for (IndexType ii = 0, *permit = perm1; ii < _n1; ii++, ++permit)
	  {
	    i = *permit;
	    // find the min among unassigned columns
	    cmin = mxdt;
	    for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++)
	    {
	      const DataType &cij = _C[sub2idx(i,*pti_unass,_n1)];
	      if (cij  < cmin) { cmin = cij; pti_min = pti_unass; }
	    }
	    // assign the column which provides the minimum and update the approximate solution
	    imin = *pti_min; prim21[imin] = i; prim12[i] = imin;
	    *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
	    approx += cmin;
	  }
	}
      }
      delete[] unass;
      if (delete21) { delete[] prim21; prim21 = NULL; }
      return approx;
    }
    // -----------------------------------------------------------
    // Refined greedy estimation for _n1 >= _n2
    // return the cost of the estimation
    // -----------------------------------------------------------
    DataType greedyRefined2(IndexType *prim12, IndexType *prim21 = NULL,
			    IndexType *perm1 = NULL, IndexType *perm2 = NULL)
    {
      IndexType nass = 0, i, j, imin, jmin, *unassi = new IndexType[_n1+1], *pti_unass = NULL;
      IndexType *pti_unass_beg = unassi, *pti_min = NULL;
      const IndexType *pti_unass_end = unassi+_n1;
      IndexType *ptj_unass = NULL, *unassj = new IndexType[_n2+1];
      IndexType *ptj_unass_beg = unassj, *ptj_unass_end = unassj+_n2, *ptj_min = NULL;
      IndexType *ptj_unass1 = NULL;
      DataType cmin, ckmin, approx = _zero, c_cur = _zero;
      const DataType mxdt = std::numeric_limits<DataType>::max();
      bool delete21 = false;
      if (prim21 == NULL) { delete21 = true; prim21 = new IndexType[_n2]; }
      
      if (perm1 == NULL) for (i = 0; i < _n1; i++) { unassi[i] = i; prim12[i] = _n2; }
      else for (i = 0; i < _n1; i++) { unassi[i] = perm1[i]; prim12[i] = _n2; }
      if (perm2 == NULL) for (j = 0; j < _n2; j++) { unassj[j] = j; prim21[j] = _n1; }
      else for (j = 0; j < _n2; j++) { unassj[j] = perm2[j]; prim21[j] = _n1; }
      
      // augmentation of columns - 2nd set
      for (ptj_unass1 = ptj_unass_beg; ptj_unass1 != ptj_unass_end;)
      {
	j = *ptj_unass1;
	// find the min among unassigned rows - 1st set
	cmin = mxdt;
	for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++) {
	  c_cur = _C[sub2idx(*pti_unass,j,_n1)];
	  if (c_cur  < cmin) { cmin = c_cur; pti_min = pti_unass; }
	}
	// find the min among unassigned columns for imin
	imin = *pti_min;
	ckmin = mxdt;
	for (ptj_unass = ptj_unass_beg; ptj_unass != ptj_unass_end; ptj_unass++) {
	  c_cur = _C[sub2idx(imin,*ptj_unass,_n1)];
	  if (c_cur  < ckmin) { ckmin = c_cur; ptj_min = ptj_unass; }
	}
	// assign the row and column which provides the minimum
	if (cmin <= ckmin) {
	  prim12[imin] = j; prim21[j] = imin;
	  *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
	  ptj_unass_beg++;
	  approx += cmin;
	}
	else {
	  jmin = *ptj_min; prim12[imin] = jmin; prim21[jmin] = imin;
	  *ptj_min = *ptj_unass_beg; *ptj_unass_beg = jmin; ptj_unass_beg++;
	  *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
	  approx += ckmin;
	}
	ptj_unass1 = ptj_unass_beg;
      }
      
      delete[] unassi; delete[] unassj;
      if (delete21) { delete[] prim21; prim21 = NULL; }
      return approx;
    }
    // -----------------------------------------------------------
    // Refined greedy estimation for _n1 <= _n2
    // return the cost of the estimation
    // -----------------------------------------------------------
    DataType greedyRefined1(IndexType *prim12, IndexType *prim21 = NULL,
			    IndexType *perm1 = NULL, IndexType *perm2 = NULL)
    {
      IndexType nass = 0, i, j, imin, jmin;
      IndexType *unassi = new IndexType[_n1+1], *pti_unass = NULL;
      IndexType *pti_unass_beg = unassi, *pti_unass_end = unassi+_n1, *pti_min = NULL;
      IndexType *ptj_unass = NULL, *unassj = new IndexType[_n2+1];
      IndexType *ptj_unass_beg = unassj, *ptj_min = NULL;
      const IndexType *ptj_unass_end = unassj+_n2;
      IndexType *pti_unass1 = NULL;
      DataType cmin, ckmin, approx = _zero, c_cur = _zero;
      const DataType mxdt = std::numeric_limits<DataType>::max();
      bool deleteprim21 = false;
      if (prim21 == NULL) { deleteprim21 = true; prim21 = new IndexType[_n2]; }
      
      if (perm1 == NULL) for (i = 0; i < _n1; i++) { unassi[i] = i; prim12[i] = _n2; }
      else for (i = 0; i < _n1; i++) { unassi[i] = perm1[i]; prim12[i] = _n2; }
      if (perm2 == NULL) for (j = 0; j < _n2; j++) { unassj[j] = j; prim21[j] = _n1; }
      else for (j = 0; j < _n2; j++) { unassj[j] = perm2[j]; prim21[j] = _n1; }
      
      // augmentation of rows
      for (pti_unass1 = pti_unass_beg; pti_unass1 != pti_unass_end;)
      {
	i = *pti_unass1;
	// find the min among unassigned columns
	cmin = mxdt;
	for (ptj_unass = ptj_unass_beg; ptj_unass != ptj_unass_end; ptj_unass++) {
	  c_cur = _C[sub2idx(i,*ptj_unass,_n1)];
	  if (c_cur  < cmin) { cmin = c_cur; ptj_min = ptj_unass; }
	}
	// find the min among unassigned rows for jmin
	jmin = *ptj_min;
	ckmin = mxdt;
	for (pti_unass = pti_unass_beg; pti_unass != pti_unass_end; pti_unass++){
	  c_cur = _C[sub2idx(*pti_unass,jmin,_n1)];
	  if (c_cur  < ckmin) { ckmin = c_cur; pti_min = pti_unass; }
	}
	// assign the row and column which provides the minimum
	if (cmin <= ckmin) {
	  prim21[jmin] = i; prim12[i] = jmin;
	  *ptj_min = *ptj_unass_beg; *ptj_unass_beg = jmin; ptj_unass_beg++;
	  pti_unass_beg++; 
	  approx += cmin;
	}
	else {
	  imin = *pti_min; prim21[jmin] = imin; prim12[imin] = jmin;
	  *pti_min = *pti_unass_beg; *pti_unass_beg = imin; pti_unass_beg++;
	  *ptj_min = *ptj_unass_beg; *ptj_unass_beg = jmin; ptj_unass_beg++;
	  approx += ckmin;
	}
	pti_unass1 = pti_unass_beg;
      }
  
      delete[] unassi; delete[] unassj;
      if (deleteprim21) delete[] prim21;
      return approx;
    }
    // --------------------------------------------------------------
    template <class BPGraph>
    void equivalenceGraph(Matching<IndexType> &prim, BPGraph &G)
    {
      IndexType i, j = 0, k;
      DataType val = 0;
      for (; j < _n2; j++) {
	k = G.setidx2idx(1,j);
	for (i = 0; i < _n1; i++) {
	  val = std::abs(_C[sub2idx(i,j,_n1)] - (dual(0,i) + dual(1,j)));
	  if (val == 0)
	    if (prim(0,i) == j) G.insertArc(G.setidx2idx(0,i),k);
	    else G.insertArc(k,G.setidx2idx(0,i));
	  else if (prim(0,i) == j) std::cout << "bug equivalenceGraph val=" << val << std::endl;
	}
      }
    }
    // --------------------------------------------------------------
    void _clearDual()
    {
      if (_dual[0]) for (IndexType i = 0; i < _n1; i++) dual(0,i) = 0;
      if (_dual[1]) for (IndexType j = 0; j < _n2; j++) dual(1,j) = 0;
    }
    // --------------------------------------------------------------
    bool isPrimalDual(Matching<IndexType> &prim)
    {
      IndexType i = 0, j = 0;
      DataType a, b, c, d, e;
      for (; j < _n2; j++)
	for (i = 0; i < _n1; i++) {
	  if (prim(0,i) == j && _C[sub2idx(i,j,_n1)]!=(dual(0,i)+dual(1,j))) {
	    a = dual(0,i); b = dual(1,j); c = _C[sub2idx(i,j,_n1)];
	    e = a+b;
	    d = std::abs(c-e);
	    std::cout  << a << " " << b << " " << e << " " << c << " "
		       << d << " " << (e==c) << std::endl;
	    return false;
	  }
	}
      return true;
    }
    // --------------------------------------------------------------
    void findSolution(Matching<IndexType> &prim, unsigned short init_type = 1)
    { hungarian(prim,init_type); }
    
    
  public:
    // --------------------------------------------------------------
    LSAP() : _n1(0), _n2(0), _C(NULL), _C_internal(true), _primals(NULL), 
	     _prims_internal(true), _dual_internal(true)
    { _dual[0] = NULL; _dual[1] = NULL; }
    // --------------------------------------------------------------
    LSAP(const IndexType &n1, const IndexType &n2, DataType *C = NULL)
      : _n1(n1), _n2(n2), _C(C), _C_internal(C == NULL),
	_primals(new PrimContainer), _prims_internal(true),
	_dual_internal(true), _lb_nb_prim(0), _ub_nb_prim(-1) {
      _dual[0] = new DataType[_n1];
      _dual[1] = new DataType[_n2];
    }
    // --------------------------------------------------------------
    LSAP(const IndexType &n1, const IndexType &n2, PrimContainer &prims, DataType *C = NULL)
      : _n1(n1), _n2(n2), _C(C), _C_internal(C == NULL), _primals(&prims), _prims_internal(false),
	_dual_internal(true), _lb_nb_prim(0), _ub_nb_prim(-1) {
      _dual[0] = new DataType[_n1];
      _dual[1] = new DataType[_n2];
    }
    // --------------------------------------------------------------
    /*LSAP(const IndexType &n1, const IndexType &n2, DataType *C, IndexType *prim12,
	 IndexType *prim21 = NULL, DataType *dual1 = NULL, DataType *dual2 = NULL)
      : _n1(n1), _n2(n2), _C(C), _C_internal(false), _primals(new std::list<Matchings<IndexType> >),
	_prim(NULL), _prim_internal(true), _prims_internal(true), _dual_internal(dual1 == NULL) {
      if (dual1 == NULL) _dual[0] = new DataType[_n1];
      else _dual[0] = dual1;
      if (dual2 == NULL) _dual[1] = new DataType[_n2];
      else _dual[1] = dual2;
      _primals.emplace_back(_n1,_n2,prim12,prim21);
      }*/
    // --------------------------------------------------------------
    /*LSAP(const IndexType &n1, const IndexType &n2, DataType *C, Matching<IndexType> &prim,
	 DataType *dual1 = NULL, DataType *dual2 = NULL)
      : _n1(n1), _n2(n2), _C(C), _C_internal(false), _primals(new std::list<Matchings<IndexType> >),
	_prim(NULL), _prim_internal(true), _prims_internal(true), _dual_internal(dual1 == NULL) {
      if (dual1 == NULL) _dual[0] = new DataType[_n1];
      else _dual[0] = dual1;
      if (dual2 == NULL) _dual[1] = new DataType[_n2];
      else _dual[1] = dual2;
      _primals.emplace_back(_n1,_n2,prim12,prim21);
      }*/
    // --------------------------------------------------------------
    ~LSAP()
    {
      if (_C_internal && _C) delete[] _C;
      if (_prims_internal && _primals) delete _primals;
      if (_dual_internal) {
	if (_dual[0]) delete[] _dual[0];
	if (_dual[1]) delete[] _dual[1];
      }
    }
    // --------------------------------------------------------------
    const IndexType& nbRows() { return _n1; }
    const IndexType& nbCols() { return _n2; }
    // --------------------------------------------------------------
    Matching<IndexType>& primal() //{ return _primals->front(); }
    { return _primals->nodeMap(0); }
    // --------------------------------------------------------------
    Matching<IndexType>& primal() const //{ return _primals->front(); }
    { return _primals->nodeMap(0); }
    // --------------------------------------------------------------
    PrimContainer& primals() { return *_primals; }
    // --------------------------------------------------------------
    const PrimContainer& primals() const { return *_primals; }
    // --------------------------------------------------------------
    const IndexType& primal(unsigned short idx_set, const IndexType &idx_elt)
    { return _primals->front()(idx_set,idx_elt); }
    // --------------------------------------------------------------
    bool isPrimal(Matching<IndexType> &nm)
    {
      for (IndexType i = 0, iend = nm.nbS1(); i < iend; i++)
	if (_C[sub2idx(i,nm(i),_n1)] != 0) return false;
      for (IndexType i = 0, iend = nm.nbS2(); i < iend; i++)
	if (_C[sub2idx(nm(1,i),i,_n1)] != 0) return false;
      return true;
    }
    // --------------------------------------------------------------
    bool isPrimal()
    {
      for (typename PrimContainer::iterator it = _primals->begin(), end = _primals->end(); it != end; ++it)
	if (!isPrimal(*it)) return false;
      return true;
    }
    // --------------------------------------------------------------
    DataType& dual(unsigned short idx_set, const IndexType &idx_elt)
    { return _dual[idx_set][idx_elt]; }
    // --------------------------------------------------------------
    const DataType& dual(unsigned short idx_set, const IndexType &idx_elt) const
    { return _dual[idx_set][idx_elt]; }
    // --------------------------------------------------------------
    LSAP<DataType,IndexType,PrimContainer>& clearDual()
    { _clearDual(); return *this; }
    // --------------------------------------------------------------
    /*LSAP<DataType,IndexType>& setPrimal(Matching<IndexType> &primal)
    {
      if (_prim_internal && _prim) delete _prim;
      _prim_internal = false;
      _prim = &primal;
      return *this;
      }*/
    // --------------------------------------------------------------
    /*LSAP<DataType,IndexType>& setPrimals(std::list<Matching<IndexType> > &prims)
    {
      if (_prims_internal && _prim) delete _prim;
      _prims_internal = false;
      _primals = &prims;
      return *this;
      }*/
    // --------------------------------------------------------------
    /*LSAP<DataType,IndexType>& setPrimal(IndexType *prim12, IndexType *prim21 = NULL)
    {
      Matching<IndexType> *prim = new Matching<IndexType>(_n1,_n2,prim12,prim21);
      setPrimal(*prim);
      _prim_internal = true;
      return *this;
      }*/
    // --------------------------------------------------------------
    /*DataType cost() const
    {
      DataType c = _zero;
      if (_prim) {
	for (IndexType i = 0; i < _n1; i++)
	  if (_prim->operator()(0,i) < _n2) c += _C[sub2idx(i,_prim->operator()(0,i),_n1)];
      }
      else return std::numeric_limits<DataType>::min();
      return c;
      }*/
    // --------------------------------------------------------------
    /*DataType cost(const Matching<IndexType> &m) const
    {
      DataType c = _zero;
      if (m.nbS1() != _n1 || m.nbS2() != _n2) return std::numeric_limits<DataType>::min();
      for (IndexType i = 0; i < _n1; i++)
	if (m(0,i) < _n2) c += _C[sub2idx(i,m(0,i),_n1)];
      return c;
      }*/
    // --------------------------------------------------------------
    /*    LSAP<DataType>& translate(const DataType &val)
    {
      for (DataType *pt = _C, *pte = _C+_n1*_n2; pt != pte; ++pt) *pt += val;
      return *this;
      }*/
    // --------------------------------------------------------------
    /*LSAP<DataType,IndexType>& genRandomInstance(DataType alpha = 1)
    {
      DataType mn, mx,
	mnn = std::numeric_limits<DataType>::max(), mnx = std::numeric_limits<DataType>::min();
      for (IndexType j = 0; j < _n2; j++) {
	randVecStoch(_n1,&_C[sub2idx(0,j,_n1)]);
	mn = std::numeric_limits<DataType>::max();
	mx = std::numeric_limits<DataType>::min();
	for (IndexType i = 0; i < _n1; i++) {
	  const DataType &d = _C[sub2idx(i,j,_n1)];
	  if (d < mn) mn = d;
	  if (d > mx) mx = d;
	}
	if (mn < mnn) mnn = mn;
	if (mx > mnx) mnx = mx;
      }
      alpha /= (mnx - mnn);
      for (DataType *pt = _C, *pte = _C+_n1*_n2; pt != pte; ++pt) *pt = alpha * (*pt-mnn);
      return *this;
    }
    // --------------------------------------------------------------
    LSAP<DataType,IndexType>& genMacholWienInstance(DataType c = 1)
    {
      for (IndexType j = 0; j < _n2; j++)
	for (IndexType i = 0; i < _n1; i++)
	  _C[sub2idx(i,j,_n1)] = c*(i+1)*(j+1);
      return *this;
    }
    // --------------------------------------------------------------
    LSAP<DataType,IndexType>& genRandPermInstance(unsigned int nb_perm = 1)
    {
      IndexType i = 0;
      for (IndexType j = 0; j < _n2; j++)
	for (i = 0; i < _n1; i++)
	  _C[sub2idx(i,j,_n1)] = 1;
      std::srand(unsigned(std::time(0)));
      for (; nb_perm > 0; nb_perm--) {
	std::vector<IndexType> v;
	for (i = 0; i < _n1; i++) v.push_back(i);
	std::random_shuffle(v.begin(),v.end());
	for (i = 0; i < _n1; i++) _C[sub2idx(i,v[i],_n1)] = 0;
      }
      return *this;
      }*/
    // --------------------------------------------------------------
    LSAP<DataType,IndexType,PrimContainer>& genRandInstance(IndexType nb_perm = 1,
							    IndexType nb_scc = 3,
							    bool cumul = false)
    {
      IndexType i = 0, i_beg = 0, i_end = 0, step = _n1/nb_scc, diff = _n1%nb_scc, s = 0;
      if (_C == NULL) _C = new DataType[_n1*_n2];
      for (IndexType j = 0; j < _n2; j++)
	for (i = 0; i < _n1; i++)
	  _C[sub2idx(i,j,_n1)] = nb_perm;

      std::vector<IndexType> v;
      if (!cumul) {
	for (; nb_perm > 0; nb_perm--) {
	  for (s = 0; s < nb_scc; s++) {
	    i_beg = s*step;
	    i_end = (s+1)*step + (s == nb_scc-1 ? diff : 0);
	    v.clear();
	    for (i = i_beg; i < i_end; i++) v.push_back(i);
	    std::random_shuffle(v.begin(),v.end());
	    for (i = i_beg; i < i_end; i++) _C[sub2idx(i,v[i-i_beg],_n1)] = 0;
	  }
	}
      }
      else {
	for (; nb_perm > 0; nb_perm--) {
	  for (s = 0; s < nb_scc; s++) {
	    i_beg = s*step;
	    i_end = (s+1)*step + (s == nb_scc-1 ? diff : 0);
	    v.clear();
	    for (i = i_beg; i < i_end; i++) v.push_back(i);
	    std::random_shuffle(v.begin(),v.end());
	    for (i = i_beg; i < i_end; i++) _C[sub2idx(i,v[i-i_beg],_n1)]--;
	  }
	}
      }
      
      return *this;
    }
    // --------------------------------------------------------------
    LSAP<DataType,IndexType,PrimContainer>& solve(unsigned short init_type = 1)
    {
      if (_C == NULL) return *this;
      if (_primals == NULL) { _primals = new PrimContainer; _prims_internal = true; }
      else _primals->clear();
      if (_dual[0] == NULL) { _dual_internal = true; _dual[0] = new DataType[_n1]; _dual[1] = new DataType[_n2]; }
      _primals->push(_n1,_n2);
      findSolution(_primals->nodeMap(0),init_type);
      return *this;
    }
    // --------------------------------------------------------------
    LSAP<DataType,IndexType,PrimContainer>& enumerate(IndexType nbPrim,
						      ENUM_EDG_SELECT option = ENUM_EDG_SELECT_RAND,
						      bool compute_1st = false,
						      unsigned short init_type = 1)
    {
      if (_C == NULL) return *this; // TODO error
      if (_primals == NULL) { _primals = new PrimContainer; _prims_internal = true; }
      else if (_primals->size() > 1) _primals->erase(_primals->begin()+1,_primals->end());
      
      // find a 1st solution
      if (compute_1st || _primals->size() == 0) {
	_primals->push(_n1,_n2);
	findSolution(_primals->nodeMap(0),init_type);
      }
      if (nbPrim == 1) return *this;

      // try find other solutions
      nbPrim--;

      /*if (option == ENUM_EDG_SELECT_BALANCED) {
	BipartiteGraph<NeibArc<DegreeNode<IndexType> > > G(_n1,_n2);
	equivalenceGraph<BipartiteGraph<NeibArc<DegreeNode<IndexType> > > >(_primals->nodeMap(0),G);
	G.enumMaximumMatchings(_primals->nodeMap(0),nbPrim,*_primals,option);
	return *this;
      }

      if (option <= ENUM_EDG_SELECT_RAND) {
	BipartiteGraph<NeibArc<IndexNode<IndexType> > > G(_n1,_n2);
	equivalenceGraph<BipartiteGraph<NeibArc<IndexNode<IndexType> > > >(_primals->nodeMap(0),G);
	G.enumMaximumMatchings(_primals->nodeMap(0),nbPrim,*_primals,option);
	return *this;
	}*/

      BipartiteGraph<WeightedArc<IndexType,DegreeNode<IndexType> > > G(_n1,_n2);
      equivalenceGraph<BipartiteGraph<WeightedArc<IndexType,DegreeNode<IndexType> > > >(_primals->nodeMap(0),G);
      G.enumMaximumMatchings(_primals->nodeMap(0),nbPrim,*_primals,option);
      return *this;
    }
    // --------------------------------------------------------------
    /*LSAP<DataType,IndexType,PrimContainer>& enumerat(IndexType nb_prim,
						     ENUM_EDG_SELECT select_edg = ENUM_EDG_SELECT_RAND)
    {
      if (_C == NULL) return *this; // TODO error
      if (_primals == NULL) { _primals = new PrimContainer; _prims_internal = true; }
      else if (_primals->size() > 1) _primals->erase(_primals->begin()+1,_primals->end());
      // find a 1st solution
      if (_primals->size() == 0) {
	_primals->push(_n1,_n2);
	findSolution(_primals->nodeMap(0));
      }
      if (nb_prim == 1) return *this;
      nb_prim--;
      BipartiteGraph<NeibArc<IndexNode<IndexType> > > G(_n1,_n2);
      equivalenceGraph<BipartiteGraph<NeibArc<IndexNode<IndexType> > > >(_primals->nodeMap(0),G);
      G.enumMatchings(_primals->nodeMap(0),nb_prim,*_primals,select_edg);
      return *this;
      }*/
    // --------------------------------------------------------------
    LSAP<DataType,IndexType,PrimContainer>&
    enumerateDissimilar(IndexType nbPrim, ENUM_DISS_ALGO algo = ENUM_DISS_MXW_DFS,
			ENUM_EDG_SELECT edge_select = ENUM_EDG_SELECT_RAND,
			ENUM_CYCLE_SELECT cycle_select_opt = ENUM_CYCLE_SELECT_MXW, 
			bool several_scc = true, bool compute_1st = false, unsigned short init_type = 1)
    {
      if (_C == NULL) return *this; // TODO error
      if (_primals == NULL) { _primals = new PrimContainer; _prims_internal = true; }
      else _primals->erase(_primals->begin()+1,_primals->end());
      
      // find a 1st solution
      if (compute_1st || _primals->size() == 0) {
	_primals->push(_n1,_n2);
	findSolution(_primals->nodeMap(0),init_type);
      }
      if (nbPrim == 1) return *this;

      // try find other solutions
      nbPrim--;

      BipartiteGraph<WeightedArc<IndexType,IndexNode<IndexType> > > G(_n1,_n2);
      equivalenceGraph<BipartiteGraph<WeightedArc<IndexType,IndexNode<IndexType> > > >
	(_primals->nodeMap(0),G);
      G.enumDissimilarMaximumMatchings(_primals->nodeMap(0),nbPrim,*_primals,algo,edge_select,
				       cycle_select_opt,several_scc);
      return *this;
    }
    // --------------------------------------------------------------
    /*LSAP<DataType,IndexType>& approx(unsigned short approx_type = 0)
    {
      IndexType *perm1 = NULL, *perm2 = NULL;
      if (!_prim) { _prim_internal = true; _prim = new Matching<IndexType>(_n1,_n2); }
      DataType c = std::numeric_limits<DataType>::max();
      switch (approx_type) {
      case 0: c = greedyBasic(_prim->array(0),_prim->array(1),perm1,perm2); break;
      case 1: c = (_n1 < _n2 ? greedyRefined1(_prim->array(0),_prim->array(1),perm1,perm2) :
		   greedyRefined2(_prim->array(0),_prim->array(1),perm1,perm2)); break;
	case 3: break;
      case 4: c = greedySort(_prim->array(0),_prim->array(1),perm1,perm2); break;
      case 5: c = greedyCountingSort(_prim->array(0),_prim->array(1),perm1,perm2); break;
      }
      return *this;
    }*/
    // --------------------------------------------------------------
    /*const IndexType& lbNbPrimals() 
    {
      if (_lb_nb_prim > 0) return _lb_nb_prim;
      IndexType res = 0, dout = 0;
      solve(1);
      BipartiteGraph<NeibArc<DegreeNode<IndexType> > > G(_n1,_n2);
      equivalenceGraph(_primals->front(),G);
      StronglyConnectedComponents<BipartiteGraph<NeibArc<DegreeNode<IndexType> > > > cc(G);
      cc.decompose();
      // trim arcs not in an SCC (not in a matching)
      for (IndexType n = G.minNodeIdx(); n != 0; n++) // nodes of V1
      {
	if (cc.scc(n) == NULL) // not in an SCC => arc in all matchings
	  G.deleteArcs(n);
	else // in an SCC, unconnect arcs to a node in a different SCC
	  for (NeibArc<DegreeNode<IndexType> > *arc = G.outgoingArc(n),
		 *a_next = NULL; arc != NULL; arc = a_next)
	  {
	    a_next = arc->next;
	    if (cc.scc(arc->to->index) != cc.scc(n)) { G.unconnect(arc); delete arc; arc = NULL; }
	  }
      }
      for (IndexType n = 0, mx = G.maxNodeIdx()+1; n != mx; n++) // nodes of V2
      {
	if (cc.scc(n) == NULL) // not in an SCC => arc in all matchings
	  G.deleteArcs(n);
	else // in an SCC, unconnect arcs to a node in a different SCC
	  for (NeibArc<DegreeNode<IndexType> > *arc = G.outgoingArc(n),
		 *a_next = NULL; arc != NULL; arc = a_next)
	  {
	    a_next = arc->next;
	    if (cc.scc(arc->to->index) != cc.scc(n)) { G.unconnect(arc); delete arc; arc = NULL; }
	  }
	dout = G.node(n)->degOut()-1;
	if (dout > 0) res += dout;
      }
      std::cout << res << " " << cc.nbSCC() << std::endl;
      cc.print();
      _lb_nb_prim = res+cc.nbSCC();
      return _lb_nb_prim;
      }*/
    // --------------------------------------------------------------
    // Number of injections from the smallest to the largest set
    // does not consider the cost, ie, for any instance
    // max(_n1,_n2)!/(max(_n1,_n2)-min(_n1,_n2))!
    /*IndexType ubNbPrimalsAnyInstance()
    {
      IndexType res = 1;
      if (_n1 < _n2) for (IndexType i = 0, iend = _n2-_n1; i < iend; i++) res *= _n2-i;
      else
	if (_n1 > _n2) for (IndexType i = 0, iend = _n1-_n2; i < iend; i++) res *= _n1-i;
	else for (IndexType i = 0; i < _n1; i++) res *= _n1-i;
      if (_ub_nb_prim < 0) _ub_nb_prim = res;
      return res;
      }*/
    // --------------------------------------------------------------
    /*const IndexType& ubNbPrimals()
    {
      if (_ub_nb_prim > 0) return _ub_nb_prim;
      return _ub_nb_prim;
      }*/
    // --------------------------------------------------------------
    LSAP<DataType,IndexType,PrimContainer>&  printInstance()
    {
      for (IndexType i = 0; i < _n1; i++)
      {
	for (IndexType j = 0; j < _n2; j++)
	  std::cout << _C[sub2idx(i,j,_n1)] << "  ";
	std::cout << std::endl;
      }
      return *this;
    }
    // --------------------------------------------------------------
    LSAP<DataType,IndexType,PrimContainer>&  loadInstance(const char* filename)
    {
      std::ifstream fin;
      fin.exceptions(std::ifstream::failbit|std::ifstream::badbit);
      try {
	fin.open(filename);

	IndexType nbr = 0, nbc = 0;
	fin >> nbr >> nbc;
	
	if (nbr != _n1 || nbc != _n2 || !_C_internal) {
	  if (_C_internal) delete[] _C;
	  else _C_internal = true;
	  _n1 = nbr; _n2 = nbc;
	  _C = new DataType[_n1*_n2];
	}
	
	for (IndexType j = 0, i = 0; i < _n1; i++)
	  for (j = 0; j < _n2; j++)
	    fin >> _C[sub2idx(i,j,_n1)];
	
	fin.close();
      }
      catch (std::ifstream::failure e) {
	std::cerr << "Exception opening/reading/closing file " << filename << std::endl;
	exit(1);
      }
      
      return *this;
    }

    // --------------------------------------------------------------
    LSAP<DataType,IndexType,PrimContainer>&  saveInstance(const char* filename)
    {
      std::ofstream fout(filename);
      fout << _n1 << " " << _n2 << std::endl;
      for (IndexType i = 0; i < _n1; i++)
      {
	for (IndexType j = 0; j < _n2; j++) fout << _C[sub2idx(i,j,_n1)] << "  ";
	fout << std::endl;
      }
      return *this;
    }
  };

}

#endif
