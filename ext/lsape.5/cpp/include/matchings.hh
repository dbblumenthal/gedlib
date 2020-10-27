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
   Creation: October 5 2015
   Last modif: Feb 2020
*/

#ifndef __MATCHINGS_HH__
#define __MATCHINGS_HH__

#include <iostream>
#include <list>
#include <map>
#include <vector>
#include "utils.hh"
#include <algorithm>

#ifdef LSAPE_PARA
#include <omp.h>
#endif

namespace liblsap {
  
  // ==============================================================
  template <typename IndexType = int>
  class Matching {

  protected:
    IndexType _n1;    // number of elements in S1
    IndexType _n2;    // number of elements in S2
    IndexType *_m[2]; // Without dummy nodes:
                      // if (_n1 < _n2) then _m[1][j]=_n1 <=> j unmatched
                      // if (_n1 > _n2) then _m[0][i]=_n2 <=> i unmatched
                      // in other cases _m[0][i]=j <=> _m[1][j]=i
    bool _internal[2];
    bool _ec;
  public:
    // --------------------------------------------------------------
    // null matching
    Matching(bool ec = false) : _n1(0), _n2(0), _ec(ec)
    { //std::cout << "create matching Matching(bool ec = false) \n";
      _m[0] = NULL; _m[1] = NULL; _internal[0] = _internal[1] = true; }
    // --------------------------------------------------------------
    // balanced matching
    Matching(const IndexType &n1, bool ec = false) : _n1(n1), _n2(n1), _ec(ec)
    {
      //std::cout << "create matching Matching(const IndexType &n1, bool ec = false)" << std::endl;
      _m[0] = new IndexType[n1]; _m[1] = NULL;
      _internal[0] = _internal[1] = true;
      for (IndexType* m = _m[0], *end = _m[0]+n1; m != end; m++) *m = n1;
    }
    // --------------------------------------------------------------
    // matching with (init_value=true by default) or withoot (=false)
    // initialization of assignment to _n1 or _n2
    Matching(IndexType n1, IndexType n2, bool ec = false, bool init_values = true)
      : _n1(n1), _n2(n2), _ec(ec)
    {
      //std::cout << "create matching Matching(IndexType n1, IndexType n2, bool ec = false, bool init_values = true)" << std::endl;
      _m[0] = new IndexType[n1]; _m[1] = new IndexType[n2];
      _internal[0] = true; _internal[1] = true;
      if (init_values) clear();
    }
    // --------------------------------------------------------------
    // matching with (init_value=true by default) or withoot (=false)
    // initialization of assignment to _n1 or _n2
    /*Matching(const IndexType &n1, const IndexType &n2,
	     IndexType *m12, IndexType *m21 = NULL, bool ec = false, bool init_values = true)
      : _n1(n1), _n2(n2), _ec(ec)
      {
      _m[0] = m12;
      _internal[0] = false;
      if (m21 == NULL) { _m[1] = new IndexType[n2]; _internal[1] = true; }
      else { _internal[1] = false; _m[1] = m21; }
      if (init_values) clear();
      }*/
    // --------------------------------------------------------------
    // copy constructor
    Matching(const Matching<IndexType> &m)
      : _n1(m._n1), _n2(m._n2), _ec(m._ec)
    {
      //std::cout << "create matching Matching(Matching<IndexType> &m)" << std::endl;
      _internal[0] = true; _internal[1] = true;
      if (m._m[0]) {
      	_m[0] = new IndexType[_n1];
	for (IndexType i = 0; i < _n1; i++) _m[0][i] = m._m[0][i];
	//for (IndexType *i1 = _m[0], *i2 = m._m[0], *end = _m[0]+_n1; i1 != end; i1++, i2++)
	  //*i1 = *i2;
      }
      else _m[0] = NULL;
      if (m._m[1]) {
	_m[1] = new IndexType[_n2];
	for (IndexType i = 0; i < _n2; i++) _m[1][i] = m._m[1][i];
	//for (IndexType *i1 = _m[1], *i2 = m._m[1], *end = _m[1]+_n2; i1 != end; i1++, i2++)
	  //*i1 = *i2;
      }
      else _m[1] = NULL;
    }
    // --------------------------------------------------------------
    ~Matching() {
      //std::cout << "delete matching" << std::endl;
      if (_internal[0] && _m[0]) { delete[] _m[0]; _m[0] = NULL; }
      if (_internal[1] && _m[1]) { delete[] _m[1]; _m[1] = NULL; }
    }
    // --------------------------------------------------------------
    const IndexType& nbS1() const { return _n1; }
    const IndexType& nbS2() const { return _n2; }
    const IndexType& nbS(const unsigned short &idx_set) const { return (idx_set == 0 ? _n1 : _n2); }
    IndexType* array(unsigned short idx_set) { return _m[idx_set]; }
    const IndexType* array(unsigned short idx_set) const { return _m[idx_set]; }
    const bool& ec() const { return _ec; }
    // --------------------------------------------------------------
    // initialize the assignment to _n2 or _n1
    void clear()
    {
      if (_m[0])
	for (IndexType* m = _m[0], *end = _m[0]+_n1-1; m <= end; m++) *m = _n2;
      if (_m[1])
	for (IndexType* m = _m[1], *end = _m[1]+_n2-1; m <= end; m++) *m = _n1;
    }
    // --------------------------------------------------------------
    IndexType& operator()(const IndexType &n) { return _m[0][n]; }
    // --------------------------------------------------------------
    const IndexType& operator()(const IndexType &n) const { return _m[0][n]; }
    // --------------------------------------------------------------
    IndexType& operator()(const unsigned short &from1or2, const IndexType &n)
    { return _m[from1or2][n]; }
    // --------------------------------------------------------------
    const IndexType& operator()(const unsigned short &from1or2, const IndexType &n) const
    { return _m[from1or2][n]; }
    // --------------------------------------------------------------
    bool covered(const IndexType &n) { return (_m[0][n] < _n2); }
    // --------------------------------------------------------------
    bool covered(const unsigned short &from1or2, const IndexType &n)
    { return (_m[from1or2][n] < (from1or2 == 0 ? _n2 : _n1)); }
    // --------------------------------------------------------------
    void transpose()
    {
      IndexType ntmp = _n1;
      _n1 = _n2; _n2 = ntmp;
      IndexType *m = _m[0];
      _m[0] = _m[1];
      _m[1] = m;
      bool inter = _internal[0];
      _internal[0] = _internal[1];
      _internal[1] = inter;
    }
    // --------------------------------------------------------------
    void rand(IndexType percent=0)
    {
      bool n1_inf_n2 = _n1 <= _n2;
      IndexType i, j;
      if (!n1_inf_n2) transpose();
      for (j = 0; j < _n2; j++) { _m[1][j] = _n1; }
      std::vector<IndexType> v;
      for (j = 0; j < _n2; j++) v.push_back(j);
      std::random_shuffle(v.begin(),v.end());
      for (i = 0; i < _n1; i++) {
	_m[0][i] = v[i];
	_m[1][v[i]] = i;
      }
      if (_ec && percent > 0) {
	for (IndexType nbEC = _n1 * percent / 100, rd; nbEC > 0; nbEC--) {
	  rd = randInt(0,_n1-1);
	  j = _m[0][rd];
	  if (j == _n1) continue;
	  _m[1][j] = _n1;
	  _m[0][rd] = _n2;
	}
      }
      if (!n1_inf_n2) transpose();
    }
    // --------------------------------------------------------------
    void print(bool also2to1 = false) const
    {
      std::cout << "matching " << _n1 << " to " << _n2 << ": ";
      for (IndexType i = 0; i < _n1; i++) {
	std::cout << i << "->" << _m[0][i] << ",";
      }
      std::cout << std::endl;
      if (also2to1) {
	std::cout << "matching " << _n2 << " to " << _n1 << ": ";
	for (IndexType i = 0; i < _n2; i++) {
	  std::cout << i << "->" << _m[1][i] << ",";
	}
	std::cout << std::endl;
      }
    }
    // --------------------------------------------------------------
    double nsim(Matching<IndexType> &other)
    {
      if (other.nbS1() != nbS1() || other.nbS2() != nbS2()) return 0;
      unsigned long int res[2];
      for (unsigned short c = 0; c < 2; c++) {
	res[c] = 0; 
	for (const IndexType *im1 = array(c), *im2 = other.array(c),
	       *ime = array(c)+nbS(c)/*, imee = m1.nbS(1-c)*/; im1 != ime; im1++, im2++)
	  res[c] += (*im1 == *im2);
      }
      const double nbV = nbS1()+other.nbS2();
      return static_cast<double>(res[0])/nbV+static_cast<double>(res[1])/nbV;
    }
  
  };

  // ==================================================================
  // Global operators
  template <typename IndexType>
  bool operator==(const Matching<IndexType> &m1, const Matching<IndexType> &m2)
  {
    if ((m1.nbS1() != m2.nbS1()) || (m1.nbS2() != m2.nbS2())) return false;
    for (IndexType im1 = 0; im1 < m1.nbS1(); im1++)
      if (m1(im1) != m2(im1)) return false;
    for (IndexType im2 = 0; im2 < m2.nbS2(); im2++)
      if (m1(1,im2) != m2(1,im2)) return false;
    return true;
  }

  // --------------------------------------------------------------
  template <typename IndexType>
  long unsigned int nb_ec_matchings(IndexType n1, IndexType n2)
  {
    long unsigned int res = 0, mn12 = std::min(n1,n2);
    long unsigned int f1 = factorial(n1), f2 = factorial(n2);
    long unsigned int f12 = f1*f2;
    for (IndexType i = 1; i <= mn12; i++)
      res += f12/(factorial(n1-i)*factorial(n2-i)*factorial(i));
    return res+1;
  }
  
  // ==================================================================
  // Similarity for matchings
  /*template<typename IndexType>
  IndexType sim(Matching<IndexType> &m1, Matching<IndexType> &m2, bool ec = false)
  {
    if (m1.nbS1() != m2.nbS1() || m1.nbS2() != m2.nbS2()) return 0;
    unsigned short idx_set = (m1.nbS1() <= m2.nbS2());
    IndexType res = 0, nbS = m1.nbS(idx_set);
    if (ec) {
      for (IndexType im1 = 0; im1 < m1.nbS1(); im1++)
	if (m1(im1) == m2(im1)) res += 1;
    }
    else {
      for (IndexType im1 = 0; im1 < nbS; im1++)
	if (m1(idx_set,im1) == m2(idx_set,im1)) res += 1;
    }
    return res;
    }*/
  
  // --------------------------------------------------------------
  /*template<typename IndexType>
  double nsim(Matching<IndexType> &m1, Matching<IndexType> &m2)
  {
    if (m1.nbS1() != m2.nbS1() || m1.nbS2() != m2.nbS2()) return 0;
    unsigned long int res[2];
    for (unsigned short c = 0; c < 2; c++) {
      res[c] = 0; 
      for (const IndexType *im1 = m1.array(c), *im2 = m2.array(c),
	     *ime = m1.array(c)+m1.nbS(c); im1 != ime; im1++, im2++)
	res[c] += (*im1 == *im2);
    }
    const double nbV = m1.nbS1()+m1.nbS2();
    return ((double)res[0])/nbV+((double)res[1])/nbV;
  }*/
  
  // --------------------------------------------------------------
  /*template<typename IndexType, class Container = std::list<Matching<IndexType> > >
  IndexType sim(Container &matchings)
  {
    if (matchings.size() < 2) return 0;
    typename Container::const_iterator it2;
    IndexType res = 0;
    for (typename Container::iterator it = matchings.cbegin(), end = matchings.cend();
	 it != end; ++it) {
      it2 = it;
      ++it2;
      for (; it2 != end; ++it2) res += sim<IndexType>(*it,*it2);
    }
    return res;
    }*/

  
  
  /*  template<typename IndexType, class Container = std::list<Matching<IndexType> > >
  double sim(Container &matchings)
  {
    // construct histogram
    if (matchings.size() < 2) return 0;
    typedef std::map< std::pair<IndexType,IndexType>, IndexType > MyMap;
    MyMap histo;
    std::pair<typename MyMap::iterator,bool> edg_p;
    typename Container::const_iterator it = matchings.cbegin(), end = matchings.cend();
    IndexType i, j, n1 = it->nbS1(), n2 = it->nbS2();
    IndexType resi = 0;
    double res = 0;
    for (; it != end; ++it)
    {
      for (i = 0; i < n1; i++)
      {
	edg_p = histo.insert(std::make_pair(std::make_pair(i,(*it)(i)),1));
	if (!edg_p.second) edg_p.first->second++;
      }
      for (j = 0; j < n2; j++)
      {
	if ((*it)(1,j) == n1) {
	  edg_p = histo.insert(std::make_pair(std::make_pair((*it)(1,j),j),1));
	  if (!edg_p.second) edg_p.first->second++;
	}
      }
    }
    // compute total similarity
    for (typename MyMap::const_iterator hit = histo.cbegin(),
	   hend = histo.cend(); hit != hend; ++hit)
    {
      std::cout << hit->first.first << "->" << hit->first.second << " w=" << hit->second << std::endl;
      resi += ((hit->first.first == n1 || hit->first.second == n2) ?
	       hit->second * hit->second : 2 * hit->second * hit->second);
    }
    res = ((double)resi) / ((double)(matchings.size() * (matchings.size() - 1) * (n1+n2)));
    res -= 1.0 / ((double)(matchings.size() -1));
    return res;
    }*/

  // =============================================
  template <typename IndexType = int, class NMapType = Matching<IndexType> >
  class NodeMaps : public std::vector<NMapType> {
  public:
    typedef NMapType NMap;
    typedef std::vector<NMapType> Container;
    using Container::size;
    using Container::reserve;
    using Container::begin;
    using Container::end;
    // -------------------------
    NodeMaps() { }
    // -------------------------
    ~NodeMaps() {
      /*for (typename std::vector<NMap*>::iterator it = this->begin(), end = this->end();
	   it != end; ++it)
	   if (*it) { delete *it; *it = NULL; }*/
    }
    // -------------------------
    /*const Matching<IndexType>& nodeMap(IndexType i) const { return *((*this)[i]); }*/
    NMap& nodeMap(IndexType i) { return ((*this)[i]); }
    // -------------------------
    NMap& push(const IndexType &n1, const IndexType &n2,
	       bool ec = false, bool init_values = true)
    {
      /*NMap *m = new Matching<IndexType>(n1,n2,ec,init_values);
	this->push_back(m);*/
      this->emplace_back(n1,n2,ec,init_values);
      return this->back();
    }
    // -------------------------
    /*void clear()
    {
      //for (typename std::vector<NMap*>::iterator it = this->begin(), last = this->end(); it != last; ++it)
      //if (*it) { delete *it; *it = NULL; }
      Container::clear();
      }*/
    // -------------------------
    /*void erase(typename Container::iterator first, typename Container::const_iterator last)
    {
      //for (typename std::vector<NMap*>::iterator it = first; it != last; ++it)
      //if (*it) { delete *it; *it = NULL; }
      Container::erase(first,last);
      }*/
    // -------------------------
    void print()
    {
      for (typename Container::iterator it = this->begin(), end = this->end();
	   it != end; ++it)
	//if (*it) (*it)->print();
	it->print();
    }
    // --------------------------------------------------------------
    bool allUnique()
    {
      if (this->size() < 2) return true;

#ifdef LSAPE_PARA // parallel version
    
      const IndexType iend = this->size();
      const IndexType i1end = (iend-1) / 2;
      bool stop = false;

#ifdef LSAPE_MX_THREADS
#pragma omp parallel for num_threads(LSAPE_MX_THREADS)
#else
#pragma omp parallel for
#endif
    
      for (IndexType i1 = 0; i1 < i1end; i1++)
      {
	bool stop_tmp = false;
#pragma omp critical
	{ stop_tmp = stop; }
	if (!stop_tmp) {
	  IndexType j = i1+1;
	  for (; j < iend; j++)
	    //if (*(*this)[i1] == *(*this)[j]) break;
	    if ((*this)[i1] == (*this)[j]) break;
	  if (j == iend) {
#pragma omp critical
	    { stop_tmp = stop; }
	    if (!stop_tmp) {
	      IndexType i2 = iend-i1-2;
	      for (j = i2+1; j < iend; j++)
		//if (*(*this)[i2] == *(*this)[j]) break;
		if ((*this)[i2] == (*this)[j]) break;
	      if (j < iend) {
#pragma omp critical
		{ stop = true;}
	      }
	    }
	  }
	  else {
#pragma omp critical
	    { stop = true; }
	  }
	}
      }
      if (stop) return false;
      if (iend % 2 == 0) { // the element iend/2 is not yet treated
	IndexType j = i1end+1;
	for (; j < iend; j++)
	  //if (*(*this)[i1end] == *(*this)[j]) break;
	  if ((*this)[i1end] == (*this)[j]) break;
	if (j < iend) return false;
      }
      
#else // sequential version
      
      for (typename Container::const_iterator it = this->cbegin(), end = this->cend();
	   it != end; ++it) {
	typename Container::const_iterator it2 = it;
	++it2;
	for (; it2 != end; ++it2)
	  //if (**it == **it2)
	  if (*it == *it2)
	  {
	    //	    std::cout << "isSet ..." << std::endl;
	    //(*it2)->print();
	    return false;
	  }
      }
      
#endif
      
      return true;
    }
    // --------------------------------------------------------------
    double nsim()
    {
      if (this->size() < 2) return 0;
      
      const double den = (double(this->size()*(this->size()-1)))/2.0;
      double res = 0;
      
#ifdef LSAPE_PARA // parallel version
      
      const IndexType iend = this->size();
      const IndexType i1end = (iend-1) / 2;
      
#ifdef LSAPE_MX_THREADS
#pragma omp parallel for num_threads(LSAPE_MX_THREADS)
#else
#pragma omp parallel for
#endif
      
      for (IndexType i1 = 0; i1 < i1end; i1++)
      {
	double r = 0;
	for (IndexType j = i1+1; j < iend; j++) r += (*this)[i1].nsim((*this)[j]);
	IndexType i2 = iend-i1-2;
	for (IndexType j = i2+1; j < iend; j++) r += (*this)[i2].nsim((*this)[j]);
#pragma omp critical
	{ res += r/den; }
      }
      
      if (iend % 2 == 0) { // the element iend/2 is not yet treated
	double r = 0;
	for (IndexType j = i1end+1; j < iend; j++) r += (*this)[i1end].nsim((*this)[j]);
	res += r/den;
      }
      
#else // sequential version
      
      for (typename Container::iterator it = this->begin(), end = this->end(); it != end; ++it) {
	typename Container::iterator it2 = it;
	++it2;
	double r = 0;
	for (; it2 != end; ++it2) r += it->nsim(*it2);
	res += r / den;
      }
      
#endif
    
      return res;
    }
  };

  
} // end namespace

#endif
