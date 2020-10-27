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

#ifndef __LSAPE_HH__
#define __LSAPE_HH__

#ifndef sub2idx
#define sub2idx(i,j,n1) (j*n1+i)
#endif

#include <typeinfo>
#include <iomanip>
#include <limits>
#include <algorithm>
#include "utils.hh"
#include "matchings.hh"
#include "graphs.hh"
#include "hungarian-lsape.h"

namespace liblsap {

  template <typename DataType, typename IndexType,
	    class PrimContainer = NodeMaps<IndexType> >
  class LSAPE {

  public:
    typedef typename PrimContainer::NMap NMap;
    
  protected:
    IndexType _n1;
    IndexType _n2;
    IndexType _nr;
    IndexType _nc;
    DataType *_C;
    bool _C_internal;
    DataType *_dual[2];
    bool _dual_internal;
    PrimContainer *_primals;
    bool _prims_internal;
    
    // --------------------------------------------------------------
    void findSolution(NMap &prim, unsigned short init_type = 1)
    { hungarianLSAPE<DataType,IndexType>(_C,_nr,_nc,
					 prim.array(0),prim.array(1),
					 _dual[0],_dual[1],init_type); }
    // --------------------------------------------------------------
    template <class BPGraph>
    void equivalenceGraph(const NMap &prim, BPGraph &G)
    {
      IndexType i, j = 0, k;
      DataType zero = 0, val = 0;
      for (; j < _nc; j++) {
	k = G.setidx2idx(1,j);
	for (i = 0; i < _n1; i++) {
	  val = _C[sub2idx(i,j,_nr)] - (dual(0,i) + dual(1,j));
	  if (val == zero)
	    if (prim(0,i) == j) G.insertArc(G.setidx2idx(0,i),k);
	    else G.insertArc(k,G.setidx2idx(0,i));
	  else if (prim(0,i) == j) std::cout << val <<  " at " << i << " " << j << std::endl;
	}
      }
      // last row: dummy node in S1
      for (j = 0; j < _n2; j++) {
	val = _C[sub2idx(_n1,j,_nr)] - dual(1,j);
	if (val == zero)
	  if (prim(1,j) == _n1) G.insertArc(-1,G.setidx2idx(1,j));
	  else G.insertArc(G.setidx2idx(1,j),-1);
	else if (prim(1,j) == _n1) std::cout << val <<  " at " << _n1 << " " << j << std::endl;
      }
    }
    
  public:
    // --------------------------------------------------------------
    LSAPE() : _n1(0), _n2(0), _nr(0), _nc(0), _C(NULL), _C_internal(true), _dual_internal(true),
	      _primals(NULL), _prims_internal(true) {
      _dual[0] = NULL; _dual[1] = NULL;
    }
    // --------------------------------------------------------------
    LSAPE(const IndexType &n1, const IndexType &n2)
      : _n1(n1), _n2(n2), _nr(n1+1), _nc(n2+1), _C(new DataType[(n1+1)*(n2+1)]),
	_C_internal(true), _dual_internal(true), _primals(NULL), _prims_internal(true)  {
      _dual[0] = new DataType[_nr]; _dual[0][_n1] = 0;
      _dual[1] = new DataType[_nc]; _dual[1][_n2] = 0;
    }
    // --------------------------------------------------------------
    LSAPE(const size_t &n1, const size_t &n2, DataType *C)
      : _n1(n1), _n2(n2), _nr(n1+1), _nc(n2+1), _C(C), _C_internal(false), _dual_internal(true),
	_primals(NULL), _prims_internal(true) {
      _dual[0] = new DataType[_nr]; _dual[0][_n1] = 0;
      _dual[1] = new DataType[_nc]; _dual[1][_n2] = 0;
    }
    // --------------------------------------------------------------
    /*LSAPE(const size_t &n1, const size_t &n2, double *C)
      : _n1(n1), _n2(n2), _nr(n1+1), _nc(n2+1), _C(new DataType[(n1+1)*(n2+1)]), _C_internal(true), _dual_internal(true),
	_primals(NULL), _prims_internal(true) {
      _dual[0] = new DataType[_nr]; _dual[0][_n1] = 0;
      _dual[1] = new DataType[_nc]; _dual[1][_n2] = 0;
      
      }*/
    // --------------------------------------------------------------
    ~LSAPE()
    {
      if (_C_internal && _C) delete[] _C;
      if (_prims_internal && _primals) delete _primals;
      if (_dual_internal) {
	if (_dual[0]) delete[] _dual[0];
	if (_dual[1]) delete[] _dual[1];
      }
    }
    // --------------------------------------------------------------
    NMap& primal() { return _primals->front(); }
    // --------------------------------------------------------------
    const NMap& primal() const { return _primals->front(); }
    // --------------------------------------------------------------
    PrimContainer& primals() { return *_primals; }
    // --------------------------------------------------------------
    const PrimContainer& primals() const { return *_primals; }
    // --------------------------------------------------------------
    bool isPrimal(NMap &nm)
    {
      for (IndexType i = 0, iend = nm.nbS1(); i < iend; i++)
	if (_C[sub2idx(i,nm(i),_nr)] != 0) return false;
      for (IndexType i = 0, iend = nm.nbS2(); i < iend; i++)
	if (_C[sub2idx(nm(1,i),i,_nr)] != 0) return false;
      return true;
    }
    // --------------------------------------------------------------
    bool isPrimal()
    {
      for (typename PrimContainer::iterator it = _primals->begin(), end = _primals->end(); it != end; ++it)
	if (!isPrimal(**it)) return false;
      return true;
    }
    // --------------------------------------------------------------
    DataType& dual(unsigned short idx1or2, const IndexType &idx_elt)
    { return _dual[idx1or2][idx_elt]; }
    const DataType& dual(unsigned short idx1or2, const IndexType &idx_elt) const
    { return _dual[idx1or2][idx_elt]; }
    /*LSAPE<DataType,IndexType>& clearDual()
      { _clearDual(); return *this; }*/
    // --------------------------------------------------------------
    DataType maxInstanceValue()
    {
      IndexType n = _nr*_nc;
      DataType dres = std::numeric_limits<DataType>::min();
      for (DataType *dit = _C, *dend = _C+n; dit != dend; dit++)
	if (*dit > dres) dres = *dit;
      return dres;
    }
    // --------------------------------------------------------------
    LSAPE<DataType,IndexType,PrimContainer>& arithScaleInstance()
    {
      DataType dmax = maxInstanceValue();
      arith_scale(_nr,_nc,_C,_C,_nr+_nc,dmax);
      return *this;
    }
    // --------------------------------------------------------------
    LSAPE<DataType,IndexType,PrimContainer>& genRandInstance(IndexType nb_prim = 1,
							     IndexType ec_percent = 20,
							     bool cumul = false)
    {
      if (_C == NULL) { _C = new DataType[_nr*_nc]; _C_internal = true; }
      IndexType i = 0, j = 0;
      for (j = 0; j < _nc; j++)
	for (i = 0; i < _nr; i++)
	  _C[sub2idx(i,j,_nr)] = nb_prim;
      
      if (!cumul) {
	for (; nb_prim > 0; nb_prim--) {
	  NMap m(_n1,_n2,true);
	  m.rand(ec_percent);
	  for (i = 0; i < _n1; i++) _C[sub2idx(i,m(i),_nr)] = 0;
	  for (j = 0; j < _n2; j++)
	    if (m(1,j) == _n1) _C[sub2idx(m(1,j),j,_nr)] = 0;
	}
      }
      else {
	for (; nb_prim > 0; nb_prim--) {
	  NMap m(_n1,_n2,true);
	  m.rand(ec_percent);
	  for (i = 0; i < _n1; i++) _C[sub2idx(i,m(i),_nr)]--;
	  for (j = 0; j < _n2; j++)
	    if (m(1,j) == _n1) _C[sub2idx(m(1,j),j,_nr)]--;
	}
      }

      return *this;
    }
    // --------------------------------------------------------------
    LSAPE<DataType,IndexType,PrimContainer>& solve(unsigned short init_type = 1)
    {
      if (_C == NULL) return *this;
      if (_primals == NULL) { _primals = new PrimContainer; _prims_internal = true; }
      else _primals->clear();
      if (_dual[0] == NULL)
	{ _dual_internal = true; _dual[0] = new DataType[_nr]; _dual[1] = new DataType[_nc]; }
      _primals->push(_n1,_n2,true);
      findSolution(_primals->nodeMap(0),init_type);
      return *this;
    }
    // --------------------------------------------------------------
    /*LSAPE<DataType,IndexType,PrimContainer>& enumerate(IndexType nbPrim,
						       ENUM_OPTION option = ENUM_RAND_SCC_RAND_EDGE,
						       bool compute_1st = false,
						       unsigned short init_type = 1)
    {
      if (_C == NULL) return *this; // TODO error
      if (_primals == NULL) { _primals = new PrimContainer; _prims_internal = true; }
      else clearPrimals(compute_1st);
      
      // find a 1st solution
      if (compute_1st || _primals->size() == 0) {
	_primals->emplace_back(_n1,_n2,true);
	findSolution(_primals->front(),init_type);
      }
      if (nbPrim == 1) return *this;

      // try find other solutions
      nbPrim--;
      if (option <= ENUM_RAND_SCC_RAND_EDGE) {
	BipartiteGraphEC<NeibArc<IndexNode<IndexType> > > G(_n1,_n2);
	equivalenceGraph<BipartiteGraphEC<NeibArc<IndexNode<IndexType> > > >(_primals->front(),G);
	G.enumMaximumMatchings(_primals->front(),nbPrim,*_primals,option);
	return *this;
      }
      
      BipartiteGraphEC<NeibArc<DegreeNode<IndexType> > > G(_n1,_n2);
      equivalenceGraph<BipartiteGraphEC<NeibArc<DegreeNode<IndexType> > > >(_primals->front(),G);
      G.enumMaximumMatchings(_primals->front(),nbPrim,*_primals,option);
      return *this;
      }*/
    // --------------------------------------------------------------
    LSAPE<DataType,IndexType,PrimContainer>& enumerate(IndexType nbPrim,
						       ENUM_EDG_SELECT option = ENUM_EDG_SELECT_RAND,
						       bool compute_1st = false,
						       unsigned short init_type = 1)
    {
      if (_C == NULL) return *this; // TODO error
      if (_primals == NULL) { _primals = new PrimContainer; _prims_internal = true; }
      else if (_primals->size() > 1) _primals->erase(_primals->begin()+1,_primals->end());
      
      // find a 1st solution
      if (compute_1st || _primals->size() == 0) {
	_primals->push(_n1,_n2,true);
	findSolution(_primals->nodeMap(0),init_type);
      }
      if (nbPrim == 1) return *this;

      // try find other solutions
      nbPrim--;
      
      BipartiteGraphEC<WeightedArc<IndexType,DegreeNode<IndexType> > > G(_n1,_n2);
      equivalenceGraph<BipartiteGraphEC<WeightedArc<IndexType,DegreeNode<IndexType> > > >(_primals->nodeMap(0),G);
      G.enumMaximumMatchings(_primals->nodeMap(0),nbPrim,*_primals,option);
      return *this;
    }
    // --------------------------------------------------------------
    LSAPE<DataType,IndexType,PrimContainer>&
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
	_primals->push(_n1,_n2,true);
	findSolution(_primals->nodeMap(0),init_type);
      }
      if (nbPrim == 1) return *this;

      // try find other solutions
      nbPrim--;

      BipartiteGraphEC<WeightedArc<IndexType,IndexNode<IndexType> > > G(_n1,_n2);
      equivalenceGraph<BipartiteGraphEC<WeightedArc<IndexType,IndexNode<IndexType> > > >
	(_primals->nodeMap(0),G);
      G.enumDissimilarMaximumMatchings(_primals->nodeMap(0),nbPrim,*_primals,algo,edge_select,
				       cycle_select_opt,several_scc);
      return *this;
    }
    // --------------------------------------------------------------
    LSAPE<DataType,IndexType,PrimContainer>&  loadInstance(const char* filename)
    {
      std::ifstream fin;
      fin.exceptions(std::ifstream::failbit|std::ifstream::badbit);
      try {
	fin.open(filename);

	IndexType nbr = 0, nbc = 0;
	fin >> nbr >> nbc;
	
	if (nbr != _nr || nbc != _nc || !_C_internal) {
	  if (_C_internal && _C) delete[] _C;
	  _C_internal = true;
	  _nr = nbr; _nc = nbc; _n1 = nbr-1; _n2 = nbc-1;
	  _C = new DataType[_nr*_nc];
	}
	
	DataType val;
	for (IndexType j = 0, i = 0; i < _nr; i++)
	  for (j = 0; j < _nc; j++) {
	    fin >> val;
	    _C[sub2idx(i,j,_nr)] = val;
	  }

	_C[sub2idx(_n1,_n2,_nr)] = 0;
	fin.close();
      }
      catch (std::ifstream::failure e) {
	std::cerr << "Exception opening/reading/closing file " << filename << std::endl;
	exit(1);
      }
      
      return *this;
    }
    // --------------------------------------------------------------
    LSAPE<DataType,IndexType,PrimContainer>& saveInstance(const char* filename)
    {
      std::ofstream fout(filename);
      fout << _nr << " " << _nc << std::endl;
      for (IndexType i = 0; i < _nr; i++)
      {
	for (IndexType j = 0; j < _nc; j++) fout << _C[sub2idx(i,j,_nr)] << "  ";
	fout << std::endl;
      }
      fout.close();
      return *this;
    }
    // --------------------------------------------------------------
    LSAPE<DataType,IndexType,PrimContainer>&  printC()
    {
      for (IndexType i = 0; i < _nr; i++)
      {
	for (IndexType j = 0; j < _nc; j++)
	  std::cout << _C[sub2idx(i,j,_nr)] << "  ";
	std::cout << std::endl;
      }
      return *this;
    }
  };

}

#endif
