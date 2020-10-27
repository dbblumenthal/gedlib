// -----------------------------------------------------------
/**
 * @file enumerators.hh
 * @brief Enumerate matchings in oriented bipartite graphs
 * @author Sebastien Bougleux
 * @date September 07 2019
 * @institution UNICAEN, ENSICAEN, CNRS, GREYC, France
 */
/*
 * -----------------------------------------------------------
 * This file is part of LSAPE.
 * LSAPE is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL-C License. See README for more
 * details.
 * -----------------------------------------------------------
*/

#ifndef __ENUMERATORS_HH__
#define __ENUMERATORS_HH__

#include "scc.hh"
#include "hungarian-lsap.h"

#include <list>
#include <vector>

namespace liblsap {

  // ==============================================================
  enum ENUM_OPTION {
    ENUM_1ST_SCC_1ST_EDGE,
    ENUM_RAND_SCC_RAND_EDGE,
    ENUM_LB_NBCYCLES_SCC_MIN_DEGREE_EDGE//,
    //ENUM_MAX_NBNODES_SCC_MIN_DEGREE_EDGE
  };

  // ==============================================================
  enum ENUM_EDG_SELECT {
    ENUM_EDG_SELECT_NAIVE,
    ENUM_EDG_SELECT_RAND,
    ENUM_EDG_SELECT_BALANCED,
    ENUM_EDG_SELECT_MAX_FREQ
  };
  
  // ==============================================================
  enum ENUM_CYCLE_SELECT {
    ENUM_CYCLE_SELECT_NAIVE,
    ENUM_CYCLE_SELECT_RAND,
    ENUM_CYCLE_SELECT_MXW,
    ENUM_CYCLE_SELECT_MAX_EDG_FREQ
  };
  
  // ==============================================================
  enum ENUM_DISS_EDG {
    ENUM_DISS_EDG_MXW,
    ENUM_DISS_EDG_RAND
  };
  
  enum ENUM_DISS_ALGO {
    ENUM_DISS_DFS,
    ENUM_DISS_MXW_DFS,
    ENUM_DISS_MXW_DFS_DUMMY1_FREE,
    ENUM_DISS_MXEX
  };
  
  // ==============================================================
  template <class GraphType, class MatchingType, class PrimCont>
  class EnumMatchings {

  public:
    typedef typename GraphType::IndexType IndexType;
    typedef typename GraphType::Node Node;
    typedef typename GraphType::Arc Arc;
    typedef typename GraphType::Cycle Cycle;
    typedef typename GraphType::NodeList NodeList;
    typedef typename StronglyConnectedComponents<GraphType>::ComponentList SCCList;
    typedef typename StronglyConnectedComponents<GraphType>::Component Component;
    typedef PrimCont PrimContainer;
    
  protected:
    GraphType& _G;
    StronglyConnectedComponents<GraphType> *_scc;
    PrimContainer* _matchings;
    bool _matchings_internal;
    ENUM_EDG_SELECT _select_opt;
    IndexType _nb_max_match;
    IndexType *_visited;
    
    // -------------------------------------------------
    void (EnumMatchings<GraphType,MatchingType,PrimContainer>::*_select_arc)(Arc*&,Component*&,Component*&) = NULL;
    // -------------------------------------------------
    void (EnumMatchings<GraphType,MatchingType,PrimContainer>::*_select_arc_in_cycle)(const Cycle&,Arc*&) = NULL;
    // -------------------------------------------------
    void (EnumMatchings<GraphType,MatchingType,PrimContainer>::*_select_arc_in_scc)(Component*,Arc*&) = NULL;
    // -------------------------------------------------
    bool (EnumMatchings<GraphType,MatchingType,PrimContainer>::*_find_cycle)(Arc*,Cycle&) = NULL;
    // -------------------------------------------------
    void (EnumMatchings<GraphType,MatchingType,PrimContainer>::*_find_cycles)(Component&,std::list<Cycle>&) = NULL;
    // -------------------------------------------------
    void (EnumMatchings<GraphType,MatchingType,PrimContainer>::*_select_cycle)
    (std::list<Cycle>&,const Cycle*&) = NULL;

    // --------------------------------------------------------------
    struct PathCandidate {

      Arc *arc12;
      Arc *arc21;
      typename Arc::Weight weight;
      
      PathCandidate() : arc12(NULL), arc21(NULL), weight(0) { }
      PathCandidate(Arc *a21, Arc *a12, const  typename Arc::Weight w)
	: arc12(a12), arc21(a21), weight(w) { }
      PathCandidate(const PathCandidate &p)
	: arc12(p.arc12), arc21(p.arc21), weight(p.weight) { }
      ~PathCandidate() { }
      void update(Arc *a21, Arc *a12, const  typename Arc::Weight w)
      { arc12 = a12; arc21 = a21; weight = w; }
      void update(const PathCandidate &p) { arc12 = p.arc12; arc21 = p.arc21; weight = p.weight; }
    };
    PathCandidate **_cands;
    PathCandidate **_cands_tmp;
    size_t _nb_cands;
    size_t _nb_max_cands;
    std::pair<Arc*,Arc*> *_parent;
    typename Arc::Weight *_w;
    typename Arc::Weight _ws;
    
    // --------------------------------------------------------------
    void delete_diss()
    {
      if (_cands) {
	for (PathCandidate **pit = _cands, **pend = _cands+_nb_max_cands-1; pit <= pend; pit++)
	  if (*pit) delete *pit;
	delete[] _cands;
      }
      if (_cands_tmp) {
	for (PathCandidate **pit = _cands_tmp, **pend = _cands_tmp+_nb_max_cands-1; pit <= pend; pit++)
	  if (*pit) delete *pit;
	delete[] _cands_tmp;
      }
      _cands = NULL;
      _cands_tmp = NULL;
      _nb_max_cands = 0;
      _nb_cands = 0;
      if (_parent) delete[] _parent;
      _parent = NULL;
      if (_w) delete[] _w;
      _w = NULL;
      _ws = 0;
    }
    
    // --------------------------------------------------------------
    void set_options()
    {
      switch(_select_opt) {
      case ENUM_EDG_SELECT_NAIVE:
	_select_arc = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_arc_naive;
	break;
      case ENUM_EDG_SELECT_RAND:
	_select_arc = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_arc_rand;
	break;
      case ENUM_EDG_SELECT_MAX_FREQ:
	_select_arc = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_arc_max_freq;
	break;
      default: _select_arc = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_arc_lbc_min_deg;
      }
      _find_cycle = &EnumMatchings<GraphType,MatchingType,PrimCont>::dfs_cycle;
    }
    
    // --------------------------------------------------------------
    void set_options(const ENUM_DISS_ALGO &algo_opt, const ENUM_CYCLE_SELECT &select_cycle_opt)
    {
      switch(_select_opt) {
      case ENUM_EDG_SELECT_NAIVE:
	_select_arc_in_cycle = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_arc_naive_in_cycle;
	_select_arc_in_scc = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_arc_naive_in_scc;
	break;
      case ENUM_EDG_SELECT_RAND:
	_select_arc_in_cycle = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_arc_rand_in_cycle;
	_select_arc_in_scc = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_arc_rand_in_scc;
	break;
      case ENUM_EDG_SELECT_MAX_FREQ:
	_select_arc_in_cycle = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_arc_max_freq_in_cycle;
	_select_arc_in_scc = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_arc_max_freq_in_scc;
	break;
      default:
	_select_arc_in_cycle = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_arc_lbc_min_deg_in_cycle;
	_select_arc_in_scc = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_arc_lbc_min_deg_in_scc;
      }
      
      switch(algo_opt) {
      case ENUM_DISS_MXEX:
	_find_cycles = (_G.nbDummyNodes() > 0 ? 
			&EnumMatchings<GraphType,MatchingType,PrimCont>::max_ec_exchange :
			&EnumMatchings<GraphType,MatchingType,PrimCont>::max_exchange);
	break;
      case ENUM_DISS_DFS:
	_find_cycles = &EnumMatchings<GraphType,MatchingType,PrimCont>::dfs_cycles;
	break;
      case ENUM_DISS_MXW_DFS:
	_find_cycles = &EnumMatchings<GraphType,MatchingType,PrimCont>::mxw_dfs_cycles;
	break;
      case ENUM_DISS_MXW_DFS_DUMMY1_FREE:
	_find_cycles = &EnumMatchings<GraphType,MatchingType,PrimCont>::mxw_dfs_cycles_dummyS1_free;
	break;
      default:
	_find_cycles = &EnumMatchings<GraphType,MatchingType,PrimCont>::mxw_dfs_cycles;
      }

      switch(select_cycle_opt) {
      case ENUM_CYCLE_SELECT_RAND:
	_select_cycle = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_cycle_rand;
	break;
      case ENUM_CYCLE_SELECT_MXW:
	_select_cycle = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_cycle_mxw;
	break;
      case ENUM_CYCLE_SELECT_MAX_EDG_FREQ:
	_select_cycle = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_cycle_mxw_arc;
	break;
      case ENUM_CYCLE_SELECT_NAIVE:
	_select_cycle = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_cycle_naive;
	break;
      default:
	_select_cycle = &EnumMatchings<GraphType,MatchingType,PrimCont>::select_cycle_rand;
      }
    }
    
    // --------------------------------------------------------------
    virtual void saveMatching()
    {
      Arc *arc = NULL;
      MatchingType &m = _matchings->push(_G.nbRealNodes(0),_G.nbRealNodes(1),(_G.nbDummyNodes()!=0),false);

      // on real nodes
      if (_G.nbRealNodes(0) > 0)
      {
	for (IndexType n_idx = _G.minRealNodeIdxS1(), n_end = _G.maxRealNodeIdxS1();
	     n_idx <= n_end; n_idx++)
	{
	  arc = _G.outgoingArc(n_idx);
	  if (arc) // TODO: to remove, should always be true
	  {
	    m(0,_G.idx2setidx(n_idx)) = _G.idx2setidx(arc->to->index);
	    if (!_G.isDummyS2(arc->to->index)) {
	      m(1,_G.idx2setidx(arc->to->index)) = _G.idx2setidx(n_idx);
	      arc->weight += 2;
	    }
	    else arc->weight += 1;
	  }
	  else std::cout << "ERROR saveMatching()\n"; // TODO: to remove
	}
      }
      // on S1's dummy nodes
      if (_G.nbDummyNodes() != 0) {
	for (arc = _G.outgoingArc(_G.dummyS1idx()); arc != NULL; arc = arc->next)
	  if (!_G.isDummyS2(arc->to->index)) {
	    m(1,_G.idx2setidx(arc->to->index)) = _G.idx2setidx(_G.dummyS1idx());
	    arc->weight += 1;
	  }
      }
    }
    // --------------------------------------------------------------
    void init_weights()
    {
      Arc *arc = NULL;
      if (_G.nbRealNodes(0) > 0)
      {
	for (IndexType n_idx = _G.minRealNodeIdxS1(), n_end = _G.maxRealNodeIdxS1();
	     n_idx <= n_end; n_idx++)
	{
	  arc = _G.outgoingArc(n_idx);
	  if (arc) // TODO: to remove, should always be true
	  {
	    if (!_G.isDummyS2(arc->to->index)) arc->weight = 2;
	    else arc->weight = 1;
	  }
	  else std::cout << "ERROR init_weights()\n"; // TODO: to remove
	}
      }
      if (_G.nbDummyNodes() != 0) // on S1's dummy nodes
	for (arc = _G.outgoingArc(_G.dummyS1idx()); arc != NULL; arc = arc->next)
	  if (!_G.isDummyS2(arc->to->index)) arc->weight = 1;
    }
    // --------------------------------------------------------------
    void init_diss()
    {
      _nb_cands = 0;
      _nb_max_cands = std::max(_G.nbNodes(0),_G.nbNodes(1))+1;
      _cands = new PathCandidate*[_nb_max_cands];
      _cands_tmp = new PathCandidate*[_nb_max_cands];
      for (PathCandidate **pit = _cands, **ptit = _cands_tmp, **pend = _cands+_nb_max_cands-1;
	   pit <= pend; pit++, ptit++)
	{ *pit = NULL; *ptit = NULL; }
      _parent = new std::pair<Arc*,Arc*>[_G.nbNodes(1)];
      _ws = 0;
      _w = new typename Arc::Weight[_G.nbNodes(1)];
    }
    // --------------------------------------------------------------
    /*bool dfsCycle(const IndexType &v_idx_from, const IndexType &v_idx, Cycle &cycle)
    {
      _processedS2[v_idx] = 1;
      for (Arc* arcVU = _G.outgoingArc(v_idx); arcVU != NULL; arcVU = arcVU->next)
      {
	cycle.push_back(arcVU);
	for (Arc* arcUV = _G.outgoingArc(arcVU->to->index); arcUV != NULL; arcUV = arcUV->next) {
	  if (arcUV->to->index == v_idx) continue;  // avoid 2-cycles
	  if (arcUV->to->index == v_idx_from) {          // cycle found
	    cycle.push_front(arcUV);
	    return true;
	  }
	  if (_processedS2[arcUV->to->index] == 0) {
	    cycle.push_back(arcUV);
	    if (dfsCycle(v_idx_from,arcUV->to->index,cycle)) return true;
	    cycle.pop_back();
	  }
	}
	cycle.pop_back();
      }
      return false;
      }*/
    // --------------------------------------------------------------
    /*bool dfsCycle(const IndexType &v_idx_from, const IndexType &v_idx, Cycle &cycle)
    {
      if (!_G.isDummy(v_idx))
	_visited[_G.idx2key(v_idx)] = 1; // do not mark dummy nodes
      for (Arc* arc = _G.outgoingArc(v_idx); arc != NULL; arc = arc->next)
      {
	// avoid 2-cycles (occur only for dummy nodes)
	if (!cycle.empty() && arc->to->index == cycle.back()->from->index) continue;
	// cycle found
	if (arc->to->index == v_idx_from) {
	  (_G.inS1(v_idx_from) ? cycle.push_back(arc) : cycle.push_front(arc));
	  return true;
	}
	if (_visited[_G.idx2key(arc->to->index)] == 0) {
	  cycle.push_back(arc);
	  if (dfsCycle(v_idx_from,arc->to->index,cycle)) return true;
	  cycle.pop_back();
	}
      }
      return false;
      }*/
    // --------------------------------------------------------------
    /*void cycle_from_rand_arc(const Component &cc, Cycle &cycle)
    {
      IndexType nid = randInt(1,cc.nbNodes), cpt = 1, rid = cc.lNodes->front();
      typename Component::NodeList::const_iterator nit = cc.lNodes->cbegin(), nend = cc.lNodes->cend();
      
      for (; nit != nend; ++nit, cpt++) {
	if (cpt <= nid) {
	  if (!_G.isDummy(*nit)) rid = *nit;
	  else cpt++;
	}
	_visited[_G.idx2key(*nit)] = 0;
      }
      dfsCycle(rid,rid,cycle);
      
      for (nit = cc.lNodes->cbegin(); nit != nend; ++nit) _visited[_G.idx2key(*nit)] = 1;
      }*/
    // --------------------------------------------------------------
    /*std::pair<Component*,Component*> _LB_NBCYCLES_SCC_MIN_DEGREE_EDGE(Cycle &cycle)
    {
      Arc const *arc = NULL;
      Arc const *arc_min = NULL;
      IndexType din, dout, d_min = std::numeric_limits<IndexType>::max(), dg;
      Component *cc = NULL, *cc_prev = NULL;
      IndexType wmn = std::numeric_limits<IndexType>::max(), wmn_tmp = 0, wmx = 0;
      cycle.clear();
      
      // find SCC with min LB |E|-|V|+1
      // if several, take one with highest number of nodes
      for (Component *ccit = _scc->head(), *ccit_prev = NULL; ccit != NULL; ccit = ccit->next) {
	wmn_tmp = 0;
	for (typename Component::NodeList::const_iterator nit = ccit->lNodes->cbegin(),
	       nend = ccit->lNodes->cend(); nit != nend; ++nit)
	  if (_G.inS2(*nit)) wmn_tmp += (_G.node(*nit)->degOut()-1);
	  else if (_G.isDummyS1(*nit)) wmn_tmp += _G.node(*nit)->degOut()-1;
	if (wmn_tmp < wmn) { cc = ccit; wmn = wmn_tmp; cc_prev = ccit_prev; wmx = ccit->nbNodes; }
	else if (wmn_tmp == wmn && wmx < ccit->nbNodes)
	  { cc = ccit; wmx = ccit->nbNodes; cc_prev = ccit_prev; }
	ccit_prev = ccit;
      }

      // find an arc not incident to a dummy node
      // s.t. 
      for (typename Component::NodeList::const_iterator nit = cc->lNodes->cbegin(),
	     nend = cc->lNodes->cend(); nit != nend; ++nit)
      {
	_visited[_G.idx2key(*nit)] = 0;
	if (_G.inS1(*nit) && !_G.isDummyS1(*nit)) {  
	  arc = _G.outgoingArc(*nit);
	  if (!_G.isDummyS2(arc->to->index)) {
	    din = _G.node(*nit)->degIn();
	    dout = arc->to->degOut();
	    dg = din + dout;
	    if (dg  < d_min || (dg == d_min && din > 1 && dout > 1))
	      { d_min = dg; arc_min = arc; }
	  }
	}
      }
      
      // no arc not incident to a dummy node (in one case only)
      // there is at least one arc from S2 to S1 exluding dummy nodes
      if (arc_min == NULL)
      {
	arc = _G.outgoingArc(_G.dummyS1idx());
	if (arc != NULL && _scc->scc(arc->from->index) == cc) {
	  for (; arc != NULL; arc = arc->next) {
	    if (!_G.isDummyS2(arc->to->index) && _scc->scc(arc->to->index) == cc) {
	      dout = arc->to->degOut();
	      if (dout < d_min) { d_min = dout; arc_min = arc; }
	    }
	  }
	}
	if (arc_min == NULL) std::cout << "Error no arc_min in findCycle" << std::endl;
      }
      // find a cycle from the arc
      if (dfsCycle(arc_min->to->index,arc_min->to->index,cycle))
	{
	  for (typename Component::NodeList::const_iterator nit = cc->lNodes->cbegin(),
		 nend = cc->lNodes->cend(); nit != nend; ++nit) _visited[_G.idx2key(*nit)] = 1;
	  return std::make_pair(cc,cc_prev);
	}
      else std::cout << "Bad finCycle from " << arc_min->to->index << std::endl;
      return std::make_pair(cc,cc_prev);
      }*/
    
    // --------------------------------------------------------------
    /*std::pair<Component*,Component*> _MAX_NBNODES_SCC_MIN_DEGREE_EDGE(Cycle &cycle)
    {
      Component *cc = _scc->head(), *cc_prev = NULL;
      IndexType wmx = cc->nbNodes;
      Component *ccit = cc->next, *ccit_prev = cc;
      Arc const *arc = NULL;
      Arc const *arc_min = NULL;
      IndexType din, dout, d_min = std::numeric_limits<IndexType>::max();
      typename Component::NodeList::const_iterator nit, nend = cc->lNodes->cend();
      
      cycle.clear();
      
      for (; ccit != NULL; ccit = ccit->next) {
	if (ccit->nbNodes > wmx) { cc = ccit; wmx = ccit->nbNodes; cc_prev = ccit_prev; }
	ccit_prev = ccit;
      }
      
      // find an arc not incident to a dummy node
      for (nit = cc->lNodes->cbegin(); nit != nend; ++nit)
      {
	_visited[_G.idx2key(*nit)] = 0;
	if (_G.inS1(*nit) && !_G.isDummyS1(*nit)) {
	  arc = _G.outgoingArc(*nit);
	  if (!_G.isDummyS2(arc->to->index)) {
	    din = _G.node(*nit)->degIn();
	    dout = arc->to->degOut();
	    if (din > 1 && dout > 1)
	      if (d_min < 0) { if (-(din+dout) > d_min) { d_min = -(din+dout); arc_min = arc; }}
	      else { d_min = -(din+dout); arc_min = arc; }
	    else if (din+dout < d_min) { d_min = din+dout; arc_min = arc; }
	  }
	}
      }
      
      // no arc not incident to a dummy node (in one case only)
      // there is at least one arc from S2 to S1 exluding dummy nodes
      if (arc_min == NULL)
      {
	arc = _G.outgoingArc(_G.dummyS1idx());
	if (arc != NULL && _scc->scc(arc->from->index) == cc) {
	  for (; arc != NULL; arc = arc->next) {
	    if (!_G.isDummyS2(arc->to->index) && _scc->scc(arc->to->index) == cc) {
	      dout = arc->to->degOut();
	      if (dout < d_min) { d_min = dout; arc_min = arc; }
	    }
	  }
	}
	if (arc_min == NULL) std::cout << "Error no arc_min in findCycle" << std::endl;
      }
      // find a cycle from the arc
      if (dfsCycle(arc_min->to->index,arc_min->to->index,cycle))
      {
	for (nit = cc->lNodes->cbegin(); nit != nend; ++nit) _visited[_G.idx2key(*nit)] = 1;
	return std::make_pair(cc,cc_prev);
      }
      else std::cout << "Bad finCycle from " << arc_min->to->index << std::endl;
      return std::make_pair((Component*)NULL,(Component*)NULL);
      }*/
    // -------------------------------------------------
    void select_arc_naive(Arc* &arc, Component* &cc, Component* &cc_prev)
    {
      cc_prev = NULL; cc = NULL;
      for (typename Component::NodeList::const_iterator it = _scc->head()->lNodes->cbegin(),
	     end = _scc->head()->lNodes->cend(); it != end; ++it)
      {
	if (_G.inS1(*it) && !_G.isDummyS1(*it)) {
	  arc = _G.outgoingArc((*it));
	  if (arc != NULL) { cc = _scc->head(); return; }
	}
      }
      arc = NULL;
    }
    // --------------------------------------------------------------
    void select_arc_naive_in_cycle(const Cycle &cycle, Arc* &arc)
    {
      for (typename Cycle::const_iterator ci = cycle.cbegin(), last = cycle.cend();
	   ci != last; ++ci)
	if (_G.inS1((*ci)->from->index) && !_G.isDummyS1((*ci)->from->index)) { arc = *ci; return; }
      arc = NULL;
    }
    // --------------------------------------------------------------
    void select_arc_naive_in_scc(Component *cc, Arc* &arc)
    {
      for (typename Component::NodeList::const_iterator ci = cc->lNodes->cbegin(), last = cc->lNodes->cend();
	   ci != last; ++ci) {
	if (_G.inS1(*ci) && !_G.isDummyS1(*ci)) { arc = _G.outgoingArc(*ci); return; }
      }
      arc = NULL;
    }
    // -------------------------------------------------
    void select_arc_lbc_min_deg(Arc* &arc, Component* &cc, Component* &cc_prev)
    {
      Arc *atmp = NULL;
      IndexType din, dout, d_min = std::numeric_limits<IndexType>::max(), dg;
      IndexType wmn = std::numeric_limits<IndexType>::max(), wmn_tmp = 0, wmx = 0;
      
      // find SCC with min LB |E|-|V|+1
      // if several, take one with highest number of nodes
      for (Component *ccit = _scc->head(), *ccit_prev = NULL; ccit != NULL; ccit = ccit->next) {
	wmn_tmp = 0;
	for (typename Component::NodeList::const_iterator nit = ccit->lNodes->cbegin(),
	       nend = ccit->lNodes->cend(); nit != nend; ++nit)
	  if (_G.inS2(*nit)) wmn_tmp += (_G.node(*nit)->degOut()-1);
	  else if (_G.isDummyS1(*nit)) wmn_tmp += _G.node(*nit)->degOut()-1;
	if (wmn_tmp < wmn) { cc = ccit; wmn = wmn_tmp; cc_prev = ccit_prev; wmx = ccit->nbNodes; }
	else if (wmn_tmp == wmn && wmx < ccit->nbNodes)
	  { cc = ccit; wmx = ccit->nbNodes; cc_prev = ccit_prev; }
	ccit_prev = ccit;
      }

      // find an arc
      arc = NULL;
      for (typename Component::NodeList::const_iterator nit = cc->lNodes->cbegin(),
	     nend = cc->lNodes->cend(); nit != nend; ++nit)
      {
	if (_G.inS1(*nit) && !_G.isDummyS1(*nit)) {  
	  atmp = _G.outgoingArc(*nit);
	  //if (!_G.isDummyS2(atmp->to->index)) {
	    din = _G.node(*nit)->degIn();
	    dout = atmp->to->degOut();
	    dg = din + dout;
	    if (dg  < d_min || (dg == d_min && din > 1 && dout > 1))
	      { d_min = dg; arc = atmp; }
	    //}
	}
      }
      
      if (arc != NULL) return;
      cc = NULL; cc_prev = NULL;
    }
    // --------------------------------------------------------------
    void select_arc_lbc_min_deg_in_cycle(const Cycle &cycle, Arc* &arc) // TODO
    {
      for (typename Cycle::const_iterator ci = cycle.cbegin(), last = cycle.cend();
	   ci != last; ++ci) {
	if (_G.inS1((*ci)->from->index) && !_G.isDummyS1((*ci)->from->index)) { arc = *ci; return; }
      }
      arc = NULL;
    }
    // --------------------------------------------------------------
    void select_arc_lbc_min_deg_in_scc(Component *cc, Arc* &arc) // TODO
    {
      for (typename Component::NodeList::const_iterator ci = cc->lNodes->cbegin(), last = cc->lNodes->cend();
	   ci != last; ++ci) {
	if (_G.inS1(*ci) && !_G.isDummyS1(*ci)) { arc = _G.outgoingArc(*ci); return; }
      }
      arc = NULL;
    }
    // -------------------------------------------------
    void select_arc_rand(Arc* &arc, Component* &cc, Component* &cc_prev)
    {
      // random scc
      std::pair<Component*,Component*> cmp = _scc->randComponent();
      cc = cmp.first; cc_prev = cmp.second;
      // random node in scc
      IndexType rd = randInt(1,cmp.first->nbNodes), n_idx = std::numeric_limits<IndexType>::max();
      
      typename Component::NodeList::const_iterator it = cmp.first->lNodes->cbegin(),
	end = cmp.first->lNodes->cend();
      
      for (; rd > 0 && it != end; ++it, rd--) {
	if (_G.inS1(*it) && !_G.isDummyS1(*it)) n_idx = *it;
      }
      
      if (n_idx != std::numeric_limits<IndexType>::max()) { arc = _G.outgoingArc(n_idx); return; }

      for (; it != end; ++it) {
	if (_G.inS1(*it) && !_G.isDummyS1(*it)) { arc = _G.outgoingArc(*it); return; }
      }
      cc = NULL; cc_prev = NULL; arc = NULL;
    }
    // --------------------------------------------------------------
    void select_arc_rand_in_cycle(const Cycle &cycle, Arc* &arc)
    {
      IndexType rd = randInt(1,cycle.size());
      typename Cycle::const_iterator ci = cycle.cbegin(), last = cycle.cend();
      arc = NULL;
      for (; rd > 0 && ci != last; ++ci, rd--) {
	if (_G.inS1((*ci)->from->index) && !_G.isDummyS1((*ci)->from->index)) { arc = *ci; }
      }
      if (arc != NULL) return;
      
      for (; ci != last; ++ci)
	if (_G.inS1((*ci)->from->index) && !_G.isDummyS1((*ci)->from->index)) { arc = *ci; return; }
      
      select_arc_naive_in_cycle(cycle,arc);
    }
    // --------------------------------------------------------------
    void select_arc_rand_in_scc(Component *cc, Arc* &arc)
    {
      IndexType rd = randInt(1,cc->nbNodes);
      typename Component::NodeList::const_iterator ci = cc->lNodes->cbegin(), last = cc->lNodes->cend();
      arc = NULL;
      for (; rd > 0 && ci != last; ++ci, rd--) {
	if (_G.inS1(*ci) && !_G.isDummyS1(*ci)) { arc = _G.outgoingArc(*ci); }
      }
      if (arc != NULL) return;
      
      for (; ci != last; ++ci)
	if (_G.inS1(*ci) && !_G.isDummyS1(*ci)) { arc = _G.outgoingArc(*ci); return; }
      select_arc_naive_in_scc(cc,arc);
    }
    // -------------------------------------------------
    Arc* select_arc_max_freq(Component* cc)
    {
      typename Arc::Weight wmx = -1;
      Arc *res = NULL, *arc = NULL;
      for (typename Component::NodeList::const_iterator nit = cc->lNodes->cbegin(),
	     nend = cc->lNodes->cend(); nit != nend; ++nit) {
	if (_G.inS1(*nit) && !_G.isDummyS1(*nit)) {
	  arc = _G.outgoingArc(*nit);
	  if (arc->weight > wmx) { res = arc; wmx = arc->weight; }
	}
      }
      return res;
    }
    // -------------------------------------------------
    void select_arc_max_freq(Arc* &arc, Component* &cc, Component* &cc_prev)
    {
      Arc *atmp = NULL;
      IndexType wmx = 0;
      typename Arc::Weight awmx = -1;
      cc = NULL; cc_prev = NULL; arc = NULL;
      for (Component *ccit = _scc->head(), *ccit_prev = NULL; ccit != NULL; ccit = ccit->next) {
	atmp = select_arc_max_freq(ccit);
	if (atmp) {
	  if (atmp->weight > awmx)
	    { awmx = atmp->weight; cc = ccit; wmx = ccit->nbNodes; cc_prev = ccit_prev; arc = atmp; }
	  else if (atmp->weight == awmx && ccit->nbNodes > wmx)
	    { cc = ccit; wmx = ccit->nbNodes; cc_prev = ccit_prev; arc = atmp; }
	}
	ccit_prev = ccit;
      }
    }
    // --------------------------------------------------------------
    void select_arc_max_freq_in_cycle(const Cycle &cycle, Arc* &arc)
    {
      typename Arc::Weight wmx = -1;
      Arc *atmp = NULL;
      arc = NULL;
      for (typename Cycle::const_iterator ci = cycle.cbegin(), last = cycle.cend();
	   ci != last; ++ci) {
	if (_G.inS1((*ci)->from->index) && !_G.isDummyS1((*ci)->from->index) && (*ci)->weight > wmx)
	  { atmp = *ci; arc = atmp; }
      }
    }
    // --------------------------------------------------------------
    void select_arc_max_freq_in_scc(Component *cc, Arc* &arc)
    {
      typename Arc::Weight wmx = -1;
      Arc *atmp = NULL;
      arc = NULL;
      for (typename Component::NodeList::const_iterator ci = cc->lNodes->cbegin(), last = cc->lNodes->cend();
	   ci != last; ++ci) {
	if (_G.inS1(*ci) && !_G.isDummyS1(*ci)) {
	  atmp = _G.outgoingArc(*ci);
	  if (atmp && atmp->weight > wmx) { wmx = atmp->weight; arc = atmp; }
	}
      }
    }
    // --------------------------------------------------------------
    bool dfs_cycle(Arc* root_arc, Cycle &cycle)
    { return _G.dfsCycle(root_arc->from->index,cycle); }
    // --------------------------------------------------------------
    void enum_balanced(MatchingType &cur_m)
    {
      /*if (_G.findArc(-1,0)==NULL || _G.findArc(0,-1)==NULL) {
	std::cout << "where is dummy ?\n";
	_G.print();
	_scc->print();
	exit(1);
	}*/
      
      // check stopping conditions
      if (_nb_max_match <= 0 || _scc->nbSCC() == 0) return;

      // find arc
      Component *cmp = NULL, *cmp_prev = NULL;
      Arc *arc = NULL;
      (*this.*_select_arc)(arc,cmp,cmp_prev);
      if (arc == NULL) {
	std::cout << "no arc !!!" << std::endl;
	_scc->print();
	_G.print();
	exit(1);
      }
      
      // find a cycle containing the arc (starting from arc's begin node)
      Cycle cycle;
      (*this.*_find_cycle)(arc,cycle);
      if (cycle.empty()) {
	std::cout << "no cycle from arc " << arc->from->index << "->" 
		  << arc->to->index << std::endl;
	_scc->print();
	_G.print();
	std::cerr << "EXIT" << std::endl;
	exit(1);
	return;
      }
      
      const IndexType scc_idx = cmp->index;//_scc->sccIndex(arc->from->index);
      if (scc_idx == 0) {
	std::cout << "Error selection scc in enumeration\n";
	exit(1);
      }
      //std::pair<Component*,Component*> scc = _scc->sccFromIndex(scc_idx);
	// Construct new matching from the current one and the cycle, and save it
      _G.swap(cycle);
      saveMatching();
      _nb_max_match--;

      // Construct new graph not containing the selected arc from V1 to V2 in the SCC
      // represents all the matchings not containing the arc
      // the arc is not removed from memory, it is only removed from the graph structure
      _G.unconnect(arc);

      // Update the decomposition of the current scc
      // and trim the graph, arcs are only removed from the graph structure
      std::forward_list<Arc*> Larcs;
      std::pair<IndexType,IndexType> param = _scc->updateDecompAndTrim(cmp,cmp_prev,Larcs);
      
      // enum all matchings not containing arc
      enum_balanced(_matchings->nodeMap(_matchings->size()-1));

      // reconstruct G: reconnect arcs from V2 to V1
      _scc->recompose(scc_idx,cmp,cmp_prev,param.first,param.second);
      while (!Larcs.empty()) { _G.connect(Larcs.front()); Larcs.pop_front(); }
      _G.connect(arc);
      _G.swap(cycle);
      cycle.clear();
      
      // Update the decomposition of the current scc
      // and trim the graph, arcs are only removed from the graph structure
      if (!_G.isDummyS1(arc->from->index)) {
	for (Arc *tarc = _G.ingoingArc(arc->from->index), *a_next = NULL; tarc != NULL; tarc = a_next)
	  { a_next = tarc->next_to; _G.unconnect(tarc);  Larcs.push_front(tarc); }	
      }
      if (!_G.isDummyS2(arc->to->index)) {
	for (Arc *tarc = _G.outgoingArc(arc->to->index), *a_next = NULL; tarc != NULL; tarc = a_next)
	  { a_next = tarc->next; _G.unconnect(tarc);  Larcs.push_front(tarc); }	
      }
      param = _scc->updateDecompAndTrim(cmp,cmp_prev,Larcs);
            
      // enum all matchings containing arc
      enum_balanced(cur_m);
      
      // reconstruct _G: reconnect arcs 
      _scc->recompose(scc_idx,cmp,cmp_prev,param.first,param.second);
      while (!Larcs.empty()) { _G.connect(Larcs.front()); Larcs.pop_front(); }
    }
    // --------------------------------------------------------------
    void dfs_cycles(Component &cc, std::list<Cycle> &cycles)
    {
      // select an arc to begin the cycle
      Arc *arc = NULL;
      (*this.*_select_arc_in_scc)(&cc,arc);
      if (arc != NULL) {
	cycles.emplace_back();
	dfs_cycle(arc,cycles.back());
      }
    }
    
    // --------------------------------------------------------------
    // Counting sort of candidate arcs w.r.t. weight (used by findCycle at each step)
    // wmin: minimum weight in candidate list
    // wmax: maximum weight in candidate list
    void candsCountingSort(const typename Arc::Weight &wmin, const typename Arc::Weight &wmax)
    {
      const size_t nbEl = std::abs(wmin) + wmax + 1;
      size_t *count = new size_t[nbEl];
      size_t *cit = count, *citm = count;
      size_t const * const cend = count+nbEl-1;
      for (; cit <= cend; cit++) *cit = 0;
      for (size_t i = 0; i < _nb_cands; i++) count[_cands[i]->weight-wmin]++;
      for (cit = count+1; cit <= cend; cit++, citm++) *cit += *citm;
      for (size_t i = 0; i < _nb_cands; i++) {
	_cands_tmp[count[_cands[i]->weight-wmin]-1]->update(*_cands[i]);
	count[_cands[i]->weight-wmin]--;
      }
      PathCandidate** ptmp = _cands;
      _cands = _cands_tmp;
      _cands_tmp = ptmp;
      ptmp = NULL;
      delete[] count;
    }
    // --------------------------------------------------------------
    void findCycleUpdate(const IndexType &cc_nidx, Arc *arc21, Arc *arc12, const Arc *root_arc,
			 Cycle &cycle,
			 PathCandidate** &pcur, PathCandidate** &ptmp,
			 PathCandidate &mx_cand, typename Arc::Weight &mx_wv,
			 typename Arc::Weight &mnd2, typename Arc::Weight &mxd2)
    {
      const typename Arc::Weight wdiff = arc12->weight - arc21->weight; // wdiff may be negative
      const typename Arc::Weight wr = _w[arc21->from->index] + wdiff; // sum of weights from root node to arc12->to
      std::pair<Arc*,Arc*> parcs;

      // cycle found
      if (arc12 == root_arc) { 	// cannot be dummy S2 by definition. TODO why not ?
	if (wr > _ws) {         // longer cycle, save it
	  _ws = wr;
	  //std::cout << "cycle found " << _ws << std::endl;
	  cycle.clear();
	  cycle.push_front(arc21);
	  for (parcs = _parent[arc21->from->index]; parcs.first != NULL;
	       parcs = _parent[parcs.second->from->index])
	    { cycle.push_front(parcs.first); cycle.push_front(parcs.second); }
	  cycle.push_front(arc12);
	}
      }
      
      // try update
      else {
	const IndexType vp_idx = arc12->to->index;
	if (_visited[_G.idx2key(vp_idx)] == cc_nidx) {  // so never or already in stack
	  if (!_G.isDummyS2(vp_idx)) {                  // so arc12->to is unique for v
	    //std::cout << "vp_idx=" << vp_idx << " " << _visited[_G.idx2key(vp_idx)] << " "
	    //	      << cc_nidx << std::endl;

	    if (wr > _w[vp_idx]) {                      // longer path then update
	      _w[vp_idx] = wr;
	      _parent[vp_idx] = std::make_pair(arc12,arc21);
	      if (*pcur) (*pcur)->update(arc21,arc12,wdiff);
	      else { *pcur = new PathCandidate(arc21,arc12,wdiff); *ptmp = new PathCandidate(); }
	      pcur++; _nb_cands++; ptmp++;
	      if (wdiff > mxd2) mxd2 = wdiff;  // update min and max weight diff (for counting sort)
	      if (wdiff < mnd2) mnd2 = wdiff;  // ! no else due to the 1rst time wv takes a value!
	    }
	  }
	  else // vp is dummmy, so there is possibly several u for which u->vp, decide in findCycle
	    if (wdiff > mx_wv) { mx_wv = wdiff; mx_cand.update(arc21,arc12,wdiff); }
	}
      }
    }
    // --------------------------------------------------------------
    void toPositive(IndexType &a) { a = ~a+1; }
    // --------------------------------------------------------------
    typename Arc::Weight mxw_dfs_cycle_dummyS1_free(const IndexType &cc_nidx,
						    const Arc *root_arc, Cycle &cycle)
    {
      PathCandidate *mx_cand = new PathCandidate();     // candidates to dummy S2
      const Arc *cur_arc = root_arc;
      Arc *arc21 = NULL, *arc12 = NULL;
      IndexType v_key, u_idx, v_idx = root_arc->to->index;
      PathCandidate **pcur = NULL, **ptmp = NULL;
      typename Arc::Weight mxd2 = 0, mnd2 = 0, mx_wv;
      
      std::forward_list<Arc const *> nstack;            // stack
      nstack.push_front(root_arc);                      // push root in stack
      _ws = std::numeric_limits<typename Arc::Weight>::min();     // current cycle weight
      _w[v_idx] = 0;
      _parent[v_idx] = std::make_pair((Arc*)NULL,(Arc*)NULL);

      //std::cout << "SCC idx=" << cc_nidx << "\n Start with arc " << root_arc->from->index
      //	<< " -> " << root_arc->to->index << "\n";
      
      while (!nstack.empty())
      {
	/*std::cout << "bgin from " << nstack.front()->from->index
	  << " to " << nstack.front()->to->index << std::endl;*/
	//take the 1rst non-processed arc from stack
	do {
	  cur_arc = nstack.front();
	  v_idx = cur_arc->to->index;
	  v_key = _G.idx2key(v_idx);
	  nstack.pop_front();
	} while (!nstack.empty() && _visited[v_key] > 0);
	
	if (_visited[v_key] > 0) { break; }
	
	toPositive(_visited[v_key]);
	
	_nb_cands = 0; pcur = _cands; ptmp = _cands_tmp;
	mx_wv = mxd2 = std::numeric_limits<typename Arc::Weight>::min();
	mnd2 = std::numeric_limits<typename Arc::Weight>::max();
	arc21 = _G.outgoingArc(v_idx);
	if (arc21 == NULL) {
	  std::cout << "node=" << v_idx << std::endl;
	  std::cout << "cur arc=" << cur_arc->from->index << "->" << cur_arc->to->index << std::endl;
	  _scc->print();
	  _G.print();
	  exit(1);
	}
	// find candidate paths of length 2
	for (arc21 = _G.outgoingArc(v_idx); arc21 != NULL; arc21 = arc21->next)
	{
	  u_idx = arc21->to->index;
	  //std::cout << "next arc21=" << arc21->from->index << "->" << u_idx
	  //	    << " " << _visited[_G.idx2key(u_idx)] << std::endl;
	  // already visited (in another retained cycle) or dummy cycle
	  if (_visited[_G.idx2key(u_idx)] != cc_nidx || cur_arc->from->index == u_idx) continue;
	  if (_G.isDummyS1(u_idx)) { // find next arcs from dummy node in S1 to a node in S2
	    for (arc12 = _G.outgoingArc(u_idx); arc12 != NULL; arc12 = arc12->next)
	      if (arc12->to->index != v_idx) { // not a dymmy cycle
		//std::cout << "next arc12=" << arc12->from->index << "->" << arc12->to->index << std::endl;
		findCycleUpdate(cc_nidx,arc21,arc12,root_arc,cycle,pcur,ptmp,*mx_cand,mx_wv,mnd2,mxd2);
	      }
	  }
	  else {// u is not Dummy, so there is only one outgoing arc from u to S2
	    //std::cout << "next arc12=" << _G.outgoingArc(u_idx)->from->index << "->"
	    //	      << _G.outgoingArc(u_idx)->to->index << std::endl;
	    findCycleUpdate(cc_nidx,arc21,_G.outgoingArc(u_idx),root_arc,cycle,
			    pcur,ptmp,*mx_cand,mx_wv,mnd2,mxd2);
	  }
	}
	//std::cout << "next step mx_wv=" << mx_wv << std::endl;
	// check if there is an extension to dummy node in S2
	// and add arc with max weight to list of candidates
	if (mx_wv > std::numeric_limits<typename Arc::Weight>::min()) {
	  v_idx = mx_cand->arc12->to->index; // must be dummy S2
	  _w[v_idx] = mx_cand->weight + _w[mx_cand->arc21->from->index];
	  _parent[v_idx] = std::make_pair(mx_cand->arc12,mx_cand->arc21);
	  if (*pcur) (*pcur)->update(*mx_cand);
	  else { *pcur = new PathCandidate(*mx_cand); *ptmp = new PathCandidate(); }
	  pcur++; _nb_cands++; ptmp++;
	  if (mx_cand->weight > mxd2) mxd2 = mx_cand->weight;
	  if (mx_cand->weight < mnd2) mnd2 = mx_cand->weight;
	}
	
	//std::cout << "nbCands=" << _nb_cands << std::flush << std::endl;
	/*	for (size_t k = 0; k < _nb_cands; k++) {
	  std::cout << "(" << _cands[k]->arc21->from->index << " " 
		    << _cands[k]->arc21->to->index << " " << _cands[k]->arc12->from->index << " " 
		    << _cands[k]->arc12->to->index << " " << _cands[k]->weight << ")";
	}
	std::cout << std::endl;*/
	
	if (mnd2 > 0) mnd2 = 0;

	// add candidates to stack
	if (_nb_cands > 0) {
	  if (_nb_cands > 1) candsCountingSort(mnd2,mxd2); // TODO if _nb_cands==2
	  for (size_t k = 0; k < _nb_cands; k++) nstack.push_front(_cands[k]->arc12);
	}
      }
      
      delete mx_cand; mx_cand = NULL;
      return _ws;
    }
    // --------------------------------------------------------------
    typename Arc::Weight mxw_dfs_cycle(const IndexType &cc_nidx, const Arc *root_arc, Cycle &cycle)
    {
      PathCandidate *mx_cand = new PathCandidate();     // candidates to dummy S2
      const Arc *cur_arc = root_arc;
      Arc *arc21 = NULL, *arc12 = NULL;
      IndexType v_key, u_idx, v_idx = root_arc->to->index;
      PathCandidate **pcur = NULL, **ptmp = NULL;
      typename Arc::Weight mxd2 = 0, mnd2 = 0, mx_wv;
      IndexType state_u = 0;
      
      std::forward_list<Arc const *> nstack;            // stack
      nstack.push_front(root_arc);                      // push root in stack
      _ws = std::numeric_limits<typename Arc::Weight>::min();     // current cycle weight
      _w[v_idx] = 0;
      _parent[v_idx] = std::make_pair((Arc*)NULL,(Arc*)NULL);

      //std::cout << "SCC idx=" << cc_nidx << "\n Start with arc " << root_arc->from->index
      //	<< " -> " << root_arc->to->index << "\n";
      
      while (!nstack.empty())
      {
	/*std::cout << "bgin from " << nstack.front()->from->index
	  << " to " << nstack.front()->to->index << std::endl;*/
	//take the 1rst non-processed arc from stack
	do {
	  cur_arc = nstack.front();
	  v_idx = cur_arc->to->index;
	  v_key = _G.idx2key(v_idx);
	  nstack.pop_front();
	} while (!nstack.empty() && _visited[v_key] > 0);
	
	if (_visited[v_key] > 0) { break; }
	
	toPositive(_visited[v_key]);
	if (_G.isDummyS1(cur_arc->from->index)) state_u = 2;
	
	_nb_cands = 0; pcur = _cands; ptmp = _cands_tmp;
	mx_wv = mxd2 = std::numeric_limits<typename Arc::Weight>::min();
	mnd2 = std::numeric_limits<typename Arc::Weight>::max();
	arc21 = _G.outgoingArc(v_idx);
	if (arc21 == NULL) {
	  std::cout << "node=" << v_idx << std::endl;
	  std::cout << "cur arc=" << cur_arc->from->index << "->" << cur_arc->to->index << std::endl;
	  _scc->print();
	  _G.print();
	  exit(1);
	}
	// find candidate paths of length 2
	for (arc21 = _G.outgoingArc(v_idx); arc21 != NULL; arc21 = arc21->next)
	{
	  u_idx = arc21->to->index;
	  //std::cout << "next arc21=" << arc21->from->index << "->" << u_idx
	  //	    << " " << _visited[_G.idx2key(u_idx)] << std::endl;
	  // already visited (in another retained cycle) or dummy cycle
	  if (_visited[_G.idx2key(u_idx)] != cc_nidx || cur_arc->from->index == u_idx) continue;
	  if (_G.isDummyS1(u_idx)) { // find next arcs from dummy node in S1 to a node in S2
	    if (state_u != 2) {
	      for (arc12 = _G.outgoingArc(u_idx); arc12 != NULL; arc12 = arc12->next)
		if (arc12->to->index != v_idx) { // not a dymmy cycle
		  //std::cout << "next arc12=" << arc12->from->index << "->" << arc12->to->index << std::endl;
		  findCycleUpdate(cc_nidx,arc21,arc12,root_arc,cycle,pcur,ptmp,*mx_cand,mx_wv,mnd2,mxd2);
		}
	    }
	  }
	  else {// u is not Dummy, so there is only one outgoing arc from u to S2
	    //std::cout << "next arc12=" << _G.outgoingArc(u_idx)->from->index << "->"
	    //	      << _G.outgoingArc(u_idx)->to->index << std::endl;
	    findCycleUpdate(cc_nidx,arc21,_G.outgoingArc(u_idx),root_arc,cycle,
			    pcur,ptmp,*mx_cand,mx_wv,mnd2,mxd2);
	  }
	}
	//std::cout << "next step mx_wv=" << mx_wv << std::endl;
	// check if there is an extension to dummy node in S2
	// and add arc with max weight to list of candidates
	if (mx_wv > std::numeric_limits<typename Arc::Weight>::min()) {
	  v_idx = mx_cand->arc12->to->index; // must be dummy S2
	  _w[v_idx] = mx_cand->weight + _w[mx_cand->arc21->from->index];
	  _parent[v_idx] = std::make_pair(mx_cand->arc12,mx_cand->arc21);
	  if (*pcur) (*pcur)->update(*mx_cand);
	  else { *pcur = new PathCandidate(*mx_cand); *ptmp = new PathCandidate(); }
	  pcur++; _nb_cands++; ptmp++;
	  if (mx_cand->weight > mxd2) mxd2 = mx_cand->weight;
	  if (mx_cand->weight < mnd2) mnd2 = mx_cand->weight;
	}
	
	//std::cout << "nbCands=" << _nb_cands << std::flush << std::endl;
	/*	for (size_t k = 0; k < _nb_cands; k++) {
	  std::cout << "(" << _cands[k]->arc21->from->index << " " 
		    << _cands[k]->arc21->to->index << " " << _cands[k]->arc12->from->index << " " 
		    << _cands[k]->arc12->to->index << " " << _cands[k]->weight << ")";
	}
	std::cout << std::endl;*/
	
	if (mnd2 > 0) mnd2 = 0;

	// add candidates to stack
	if (_nb_cands > 0) {
	  if (_nb_cands > 1) candsCountingSort(mnd2,mxd2); // TODO if _nb_cands==2
	  for (size_t k = 0; k < _nb_cands; k++) nstack.push_front(_cands[k]->arc12);
	}
      }
      
      delete mx_cand; mx_cand = NULL;
      return _ws;
    }
    // --------------------------------------------------------------
    void mxw_dfs_cycles(Component &cc, std::list<Cycle> &cycles)
    {
      // select an arc to begin the cycle
      Arc *root_arc = NULL;
      typename Arc::Weight cycle_weight = 0;
      (*this.*_select_arc_in_scc)(&cc,root_arc);
      if (root_arc != NULL)
      {
	// init dfs
	const IndexType cc_idx = -cc.index;
	for (typename Component::NodeList::const_iterator nit = cc.lNodes->cbegin(),
	       nend = cc.lNodes->cend(); nit != nend; ++nit) {
	  _visited[_G.idx2key(*nit)] = cc_idx;
	  if (_G.inS2(*nit)) _w[*nit] = std::numeric_limits<typename Arc::Weight>::min();
	}
	// find cycle
	cycles.emplace_back();
	cycle_weight = mxw_dfs_cycle(-cc.index,root_arc,cycles.back());
	if (cycles.back().empty() || cycle_weight <= 0) cycles.pop_back();
      }
    }
    // --------------------------------------------------------------
    void mxw_dfs_cycles_dummyS1_free(Component &cc, std::list<Cycle> &cycles)
    {
      // select an arc to begin the cycle
      Arc *root_arc = NULL;
      typename Arc::Weight cycle_weight = 0;
      (*this.*_select_arc_in_scc)(&cc,root_arc);
      if (root_arc != NULL)
      {
	// init dfs
	const IndexType cc_idx = -cc.index;
	for (typename Component::NodeList::const_iterator nit = cc.lNodes->cbegin(),
	       nend = cc.lNodes->cend(); nit != nend; ++nit) {
	  _visited[_G.idx2key(*nit)] = cc_idx;
	  if (_G.inS2(*nit)) _w[*nit] = std::numeric_limits<typename Arc::Weight>::min();
	}
	// find cycle
	cycles.emplace_back();
	cycle_weight = mxw_dfs_cycle_dummyS1_free(-cc.index,root_arc,cycles.back());
	if (cycles.back().empty() || cycle_weight <= 0) cycles.pop_back();
      }
    }
    // --------------------------------------------------------------
    void max_exchange(Component &cc, std::list<Cycle> &cycles)
    {
      // MIN-PM instance
      IndexType n = cc.nbNodes;
      typename Arc::Weight *C = new typename Arc::Weight[n*n], c_mx = 2*_matchings->size()*n, wmx = 0;
      const Arc* arc = NULL;
      IndexType i1 = 0, i2 = 0;
      IndexType *g = new IndexType[n], *ginv = new IndexType[_G.nbNodes()];
      
      for (typename Component::NodeList::const_iterator nit = cc.lNodes->cbegin(),
	     nend = cc.lNodes->cend(); nit != nend; ++nit)
      {
	if (_G.inS1(*nit)) {
	  ginv[_G.idx2key(*nit)] = i1;
	  g[i1] = *nit;
	  i1++;
	}
	else {
	  ginv[_G.idx2key(*nit)] = i2;
	  i2++;
	}
      }
      
      // construct cost matrix
      for (typename Arc::Weight *cp = C, *cp_end = C+n*n; cp != cp_end; cp++) *cp = c_mx;
      for (IndexType i = 0; i < n; i++) C[i+i*n] = 0;
      for (typename Component::NodeList::const_iterator nit = cc.lNodes->cbegin(),
	     nend = cc.lNodes->cend(); nit != nend; ++nit)
      {
	if (_G.inS1(*nit)) {
	  arc = _G.outgoingArc(*nit);
	  C[(ginv[_G.idx2key(arc->to->index)]+i1)*n+ginv[_G.idx2key(*nit)]] = -arc->weight;
	  if (arc->weight > wmx) wmx = arc->weight;
	}
	else {
	  g[ginv[_G.idx2key(*nit)]+i1] = *nit;
	  for (arc = _G.outgoingArc(*nit); arc != NULL; arc = arc->next)
	    C[ginv[_G.idx2key(arc->to->index)]*n+ginv[_G.idx2key(*nit)]+i1] = arc->weight;
	}
      }
      
      // rescale to positive costs
      for (typename Arc::Weight *cp = C, *cp_end = C+n*n; cp != cp_end; cp++) *cp += wmx;
      
      IndexType *m12 = new IndexType[n];
      typename Arc::Weight *u1 = new typename Arc::Weight[n], *u2 = new typename Arc::Weight[n];
      
      
      hungarianSquareLSAP<typename Arc::Weight,IndexType>(C,n,m12,u1,u2);

      typename Arc::Weight cycle_w = 0;
      Arc *a = NULL;
      for (IndexType i = 0; i < n; i++) _visited[i] = 0;

      for (IndexType i = 0, ii = 0, j = -1; ii < i1; ii++)
      {
	if (_visited[ii] == 1) continue;
	i = ii; j = -1;
	cycle_w = 0;
	cycles.emplace_back();
	while (m12[i] != i)
	{
	  _visited[i] = 1;
	  j = m12[i];
	  a = _G.findArc(g[i],g[j]);
	  cycles.back().push_back(a);
	  cycle_w += (_G.inS1(a->from->index) ? a->weight : -a->weight);
	  if (j == ii) break;
	  i = j;
	}
	if (j != ii || cycle_w <= 0) cycles.pop_back();
      }
      
      delete[] m12; delete[] u1; delete[] u2; delete[] g; delete[] ginv; delete[] C;
    }
    // --------------------------------------------------------------
    void max_ec_exchange(Component &cc, std::list<Cycle> &cycles)
    {
      // init functions g and ginv
      IndexType *g = new IndexType[2*cc.nbNodes], *ginv = new IndexType[_G.nbNodes()];
      IndexType nbS1 = 0, nbS2 = 0, nbDumS1 = 0, nbDumS2 = 0;
      for (typename Component::NodeList::const_iterator nit = cc.lNodes->cbegin(),
	     nend = cc.lNodes->cend(); nit != nend; ++nit)
      {
	if (_G.inS1(*nit)) {
	  if (!_G.isDummyS1(*nit)) {
	    ginv[_G.idx2key(*nit)] = nbS1;
	    g[nbS1] = *nit;
	    nbS1++;
	  }
	  else nbDumS1++;
	}
	else {
	  if (!_G.isDummyS2(*nit)) {
	    ginv[_G.idx2key(*nit)] = nbS2;
	    nbS2++;
	  }
	  else nbDumS2++;
	}
      }
      IndexType R = std::min(nbS1,nbS2), nbS = nbS1+nbS2;
      IndexType n = nbS+2*R;
      if (nbDumS1 > 0) for (IndexType r = nbS, rend = nbS+R; r < rend; r++) g[r] = _G.dummyS1idx();
      if (nbDumS2 > 0) for (IndexType r = nbS+R; r < n; r++) g[r] = _G.dummyS2idx();

      // construct max ec-exchange instance
      typename Arc::Weight *C = new typename Arc::Weight[n*n], c_mx = 2*_matchings->size()*n, wmx = 0;
      const Arc* arc = NULL;
      for (typename Arc::Weight *cp = C, *cp_end = C+n*n; cp != cp_end; cp++) *cp = c_mx;
      for (IndexType i = 0; i < n; i++) C[i+i*n] = 0;
      for (IndexType i = nbS, iend = nbS+R; i < iend; i++) { C[(R+i)*n+i] = 0; C[i*n+R+i] = 0; }
      for (typename Component::NodeList::const_iterator nit = cc.lNodes->cbegin(),
	     nend = cc.lNodes->cend(); nit != nend; ++nit)
      {
	if (_G.inS1(*nit)) {
	  if (!_G.isDummyS1(*nit)) {
	    arc = _G.outgoingArc(*nit);
	    if (!_G.isDummyS2(arc->to->index))
	      C[(ginv[_G.idx2key(arc->to->index)]+nbS1)*n+ginv[_G.idx2key(*nit)]] = -arc->weight;
	    else
	      for (IndexType r = nbS+R; r < n; r++) C[r*n+ginv[_G.idx2key(*nit)]] = -arc->weight;
	    if (arc->weight > wmx) wmx = arc->weight;
	  }
	  else { // dummy S1
	    for (arc = _G.outgoingArc(*nit); arc != NULL; arc = arc->next) {
	      if (_G.isDummyS2(arc->to->index) || _scc->sccIndex(arc->to->index) != cc.index) continue;
	      for (IndexType r = nbS, rend = nbS+R; r < rend; r++) {
		C[(ginv[_G.idx2key(arc->to->index)]+nbS1)*n+r] = -arc->weight;
		if (arc->weight > wmx) wmx = arc->weight;
	      }
	    }
	  }
	}
	else { // in S2
	  if (!_G.isDummyS2(*nit)) {
	    g[ginv[_G.idx2key(*nit)]+nbS1] = *nit;
	    for (arc = _G.outgoingArc(*nit); arc != NULL; arc = arc->next) {
	      if (!_G.isDummyS1(arc->to->index))
		C[ginv[_G.idx2key(arc->to->index)]*n+ginv[_G.idx2key(*nit)]+nbS1] = arc->weight;
	      else
		for (IndexType r = nbS, rend = nbS+R; r < rend; r++)
		  C[r*n+ginv[_G.idx2key(*nit)]+nbS1] = arc->weight;
	    }
	  }
	  else { // dummy S2
	    for (arc = _G.outgoingArc(*nit); arc != NULL; arc = arc->next) {
	      if (_G.isDummyS1(arc->to->index)) continue;
	      for (IndexType r = nbS+R; r < n; r++)
		C[ginv[_G.idx2key(arc->to->index)]*n+r] = arc->weight;
	    }
	  }
	}
      }
      /*for (IndexType b1 = 0; b1 < n; b1++){
	for (IndexType b2 = 0; b2 < n; b2++)
	  std::cout << C[sub2idx(b1,b2,n)] << "\t";
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
      // rescale to positive costs
      for (typename Arc::Weight *cp = C, *cp_end = C+n*n; cp != cp_end; cp++) *cp += wmx;
      
      IndexType *m12 = new IndexType[n];
      typename Arc::Weight *u1 = new typename Arc::Weight[n], *u2 = new typename Arc::Weight[n];
      
      hungarianSquareLSAP<typename Arc::Weight,IndexType>(C,n,m12,u1,u2);

      typename Arc::Weight cycle_w = 0;
      Arc *a = NULL;
      IndexType *visited = new IndexType[n]; // TODO REMOVE ?
      for (IndexType i = 0; i < n; i++) visited[i] = 0;
      for (IndexType i = 0, ii = 0, j = -1; ii < nbS1; ii++)
      {
	if (visited[ii] == 1) continue;
	i = ii; j = -1;
	cycle_w = 0;
	cycles.emplace_back();
	while (m12[i] != i)
	{
	  visited[i] = 1;
	  j = m12[i];
	  if (!(i > nbS && j > nbS)) {
	    a = _G.findArc(g[i],g[j]);
	    cycles.back().push_back(a);
	    cycle_w += (_G.inS1(a->from->index) ? a->weight : -a->weight);
	  }
	  if (j == ii) break;
	  i = j;
	}
	if (j != ii || cycle_w <= 0) cycles.pop_back();
      }
      delete[] m12; delete[] u1; delete[] u2; delete[] g; delete[] ginv; delete[] C;
      delete[] visited;
    }
    // --------------------------------------------------------------
    void select_cycle_naive(std::list<Cycle> &cycles, const Cycle* &cycle)
    { cycle = &cycles.front(); }
    // --------------------------------------------------------------
    void select_cycle_rand(std::list<Cycle> &cycles, const Cycle* &cycle)
    {
      IndexType rd = randInt(1,cycles.size());
      for (typename std::list<Cycle>::const_iterator cit = cycles.cbegin(), last = cycles.cend(),
	     cres = cycles.cend(); cit != last; ++cit, rd--)
	if (rd == 1) { cycle = &(*cit); break; }
    }
    // --------------------------------------------------------------
    void select_cycle_mxw_arc(std::list<Cycle> &cycles, const Cycle* &cycle)
    {
      typename Arc::Weight mxw = -1;
      Arc *arc = NULL;
      for (typename std::list<Cycle>::const_iterator cit = cycles.begin(), last = cycles.cend(),
	     cres = cycles.cend(); cit != last; ++cit) {
	select_arc_max_freq_in_cycle(*cit,arc);
	if (arc && arc->weight > mxw) { mxw = arc->weight; cycle = &(*cit); }
      }
    }
    // --------------------------------------------------------------
    typename Arc::Weight cycle_weight(const Cycle &cycle)
    {
      typename Arc::Weight cw = 0;
      for (typename Cycle::const_iterator ci = cycle.cbegin(), last = cycle.cend();
	   ci != last; ++ci) {
	if (_G.inS1((*ci)->from->index)) cw += (*ci)->weight;
	else cw -= (*ci)->weight;
      }
      return cw;
    }
    // --------------------------------------------------------------
    void select_cycle_mxw(std::list<Cycle> &cycles, const Cycle* &cycle)
    {
      typename Arc::Weight mx_cycle_w = std::numeric_limits<typename Arc::Weight>::min(), wtmp = 0;
      size_t mx_nb_nodes = 0;
      for (typename std::list<Cycle>::const_iterator cit = cycles.begin(), last = cycles.cend(),
	     cres = cycles.cend(); cit != last; ++cit) {
	wtmp = cycle_weight(*cit);
	if (wtmp > mx_cycle_w) { mx_cycle_w = wtmp; cycle = &(*cit); mx_nb_nodes = cit->size(); }
	else if (wtmp == mx_cycle_w && cit->size() > mx_nb_nodes)
	  { cycle = &(*cit); mx_nb_nodes = cit->size(); }
      }
    }
    // --------------------------------------------------------------
    void find_cycles_and_arc(std::list<Cycle> &cycles, Arc* &arc)
    {
      for (Component *ccit = _scc->head(); ccit != NULL; ccit = ccit->next)
      { (*this.*_find_cycles)(*ccit,cycles); }
      
      if (!cycles.empty()) {
	const Cycle *cycle = NULL;
	(*this.*_select_cycle)(cycles,cycle);
	if (cycle != NULL) (*this.*_select_arc_in_cycle)(*cycle,arc);
      }
      else { arc = NULL; }
    }
    // --------------------------------------------------------------
    void enum_balanced_scc(const MatchingType &cur_m)
    {
      // try to find a cycle
      if (_nb_max_match <= 0 || _scc->nbSCC() == 0) return;
      std::list<Cycle> cycles;
      Arc *arc = NULL;
      //std::cout << "nbprim=" << _nb_max_match << " nbscc=" << _scc->nbSCC() << std::endl;
      find_cycles_and_arc(cycles, arc);
      //std::cout << "csize=" << cycles.size() << std::endl;
      if (cycles.empty()) { return; }
      
      if (arc == NULL) {
	std::cout << "no arc !!!" << std::endl;
	std::cout << "nb_cycles=" << cycles.size() << std::endl;
	_scc->print();
	_G.print();
	exit(1);
      }

      const IndexType scc_idx = _scc->sccIndex(arc->from->index);
      std::pair<Component*,Component*> sccp = _scc->sccFromIndex(scc_idx);

      //std::cout << std::endl << "arc=" << arc->from->index << " " << arc->to->index << std::endl;
      //_G.print();
      
      // select arc and swap along cycles and save new matching
      for (typename std::list<Cycle>::iterator cit = cycles.begin(), cend = cycles.end();
	   cit != cend; ++cit)
	_G.swap(*cit);
      saveMatching();
      _nb_max_match--;
      
      _G.unconnect(arc);
      std::forward_list<Arc*> Larcs;
      std::pair<IndexType,IndexType> param = _scc->updateDecompAndTrim(sccp.first,sccp.second,Larcs);

      enum_balanced_scc(_matchings->nodeMap(_matchings->size()-1));
      
      _scc->recompose(scc_idx,sccp.first,sccp.second,param.first,param.second);
      while (!Larcs.empty()) { _G.connect(Larcs.front()); Larcs.pop_front(); }
      _G.connect(arc);
      for (typename std::list<Cycle>::iterator cit = cycles.begin(), cend = cycles.end();
	   cit != cend; ++cit)
	_G.swap(*cit);
      cycles.clear();

      if (!_G.isDummyS1(arc->from->index)) {
	for (Arc *tarc = _G.ingoingArc(arc->from->index), *a_next = NULL; tarc != NULL; tarc = a_next)
	  { a_next = tarc->next_to; _G.unconnect(tarc);  Larcs.push_front(tarc); }	
      }
      if (!_G.isDummyS2(arc->to->index)) {
	for (Arc *tarc = _G.outgoingArc(arc->to->index), *a_next = NULL; tarc != NULL; tarc = a_next)
	  { a_next = tarc->next; _G.unconnect(tarc);  Larcs.push_front(tarc); }	
      }
      param = _scc->updateDecompAndTrim(sccp.first,sccp.second,Larcs);
      
      enum_balanced_scc(cur_m);
      
      // reconstruct _G: reconnect arcs 
      _scc->recompose(scc_idx,sccp.first,sccp.second,param.first,param.second);
      while (!Larcs.empty()) { _G.connect(Larcs.front()); Larcs.pop_front(); }
    }
    
  public:
    EnumMatchings(GraphType &G)
      : _G(G), _scc(NULL), _matchings(NULL), _matchings_internal(true),
	_select_opt(ENUM_EDG_SELECT_BALANCED), _nb_max_match(-1),
	_visited(new IndexType[G.nbNodes()]), _cands(NULL), _cands_tmp(NULL),
	_nb_cands(0), _nb_max_cands(0), _parent(NULL), _w(NULL), _ws(0)
    {
      // init _visited (DFS for cycles)
      for (IndexType *vit = _visited, *vend = _visited+_G.nbNodes(); vit != vend; vit++) *vit = 1;
    }
    // --------------------------------------------------------------
    virtual ~EnumMatchings()
    {
      delete_diss();
      if (_scc) delete _scc;
      if (_matchings_internal && _matchings) delete _matchings;
      if (_visited) delete[] _visited;
    }
    // --------------------------------------------------------------
    /*template <class Path>
    bool dfsEvenAlternatingPath(const IndexType &n_index_from, IndexType n_index, bool *visited, Path &path)
    {
      visited[n_index] = true;
      Arc *arc = _G.outgoingArc(n_index);
      if (arc == NULL) return ((path.size() > 1) && (path.size()%2 == 0));
      for (; arc != NULL; arc = arc->next)
      {
	if (visited[arc->to->index] == false)
	{
	  path.push_back(arc);
	  if (dfsPath(n_index_from,arc->to->index,visited,path)) return true;
	  path.pop_back();
	}
      }
      return ((path.size() > 1) && (path.size()%2 == 0));
    }
    // --------------------------------------------------------------
    template <class Path>
    bool dfsEvenAlternatingPath(const IndexType &n_index, Path &path)
    {
      bool *visited = new bool[nbNodes()], path_found = false;
      for (bool *v = visited, *v_end = visited+nbNodes(); v != v_end; ++v) *v = false;
      path.clear();
      path_found = dfsPath(n_index,n_index,visited,path);
      delete[] visited;
      return path_found;
    }
    // --------------------------------------------------------------
    template <class Path>
    bool findEvenAlternatingPath(Path &path, Matching<IndexType> *cur_m)
    {
      for (typename std::vector<IndexType>::const_iterator it = _nodes_S12[1].begin(),
	     end = _nodes_S12[1].end();
	   it != end; ++it)
	if (cur_m->covered(it->index) && dfsEvenAlternatingPath(it->index,path)) return true;
      return false;
      }*/
    // --------------------------------------------------------------
    void enumerate(Matching<IndexType> &matching, IndexType nbMatch,
		   PrimContainer &Lout, ENUM_EDG_SELECT option)
    {
      _select_opt = option;
      set_options();
      std::forward_list<Arc*> L;
      if (_matchings != NULL && _matchings_internal == true) delete _matchings;
      _matchings_internal = false;
      _matchings = &Lout;
      _nb_max_match = nbMatch;
      _matchings->reserve(nbMatch+1);
      if (_matchings->capacity() < static_cast<size_t>(nbMatch+1)) {
	std::cerr << "EnumMatchings::enumerate(...): enumeration may be compromized\n"
		  << "number of primal solution reduced to " << _matchings->capacity() << "\n";
	_nb_max_match = _matchings->capacity()-1;
      }

      init_weights();
      
      if (_scc != NULL) delete _scc;
      _scc = new StronglyConnectedComponents<GraphType>(_G);
      _scc->decompose();
      _scc->trimArcs(L);
      
      // enumerate the other matchings 
      if (_G.isBalanced()) enum_balanced(matching);
      
      // reconstruct initial graph
      for (typename std::forward_list<Arc*>::const_iterator lit = L.cbegin(),
	     lend = L.cend(); lit != lend; ++lit)
	_G.connect(*lit);
      
      delete _scc; _scc = NULL;
    }
    // --------------------------------------------------------------
    void enumerateDissimilar(Matching<IndexType> &matching, IndexType nbMatch, PrimContainer &Lout,
			     ENUM_DISS_ALGO algo_opt, ENUM_EDG_SELECT arc_select_opt,
			     ENUM_CYCLE_SELECT cycle_select_opt, bool several_scc = true)
    {
      _select_opt = arc_select_opt;
      set_options(algo_opt,cycle_select_opt);
      if (algo_opt == ENUM_DISS_MXW_DFS || algo_opt == ENUM_DISS_MXW_DFS_DUMMY1_FREE) init_diss();
      
      std::forward_list<Arc*> L;
      if (_matchings != NULL && _matchings_internal == true) delete _matchings;
      _matchings_internal = false;
      _matchings = &Lout;
      _nb_max_match = nbMatch;
      _matchings->reserve(nbMatch+1);
      if (_matchings->capacity() < static_cast<size_t>(nbMatch+1)) {
	std::cerr << "EnumMatchings::enumerate(...): enumeration may be compromized\n"
		  << "number of primal solution reduced to " << _matchings->capacity() << "\n";
	_nb_max_match = _matchings->capacity()-1;
      }
      
      init_weights();
      
      
      if (_scc != NULL) delete _scc;
      _scc = new StronglyConnectedComponents<GraphType>(_G);
      _scc->decompose();
      _scc->trimArcs(L);
	
      // enumerate the other matchings
      if (_G.isBalanced()) enum_balanced_scc(matching);

      // reconstruct initial graph
      for (typename std::forward_list<Arc*>::const_iterator lit = L.cbegin(),
	     lend = L.cend(); lit != lend; ++lit)
	_G.connect(*lit);
      
      delete _scc; _scc = NULL;
      if (algo_opt == ENUM_DISS_MXW_DFS) delete_diss();
    }
  };
  
}

#endif
