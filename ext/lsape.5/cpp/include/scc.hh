// -----------------------------------------------------------
/**
 * @file scc.hh
 * @brief Strongly connected component decomposition
 * @author Sebastien Bougleux
 * @date September 07 2019
 * @institution Normandie Univ, UNICAEN, ENSICAEN, CNRS, GREYC, France
 */
/*
 * -----------------------------------------------------------
 * This file is part of LIBLSAP.
 * LIBGP is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL-C License. See README for more
 * details.
 * -----------------------------------------------------------
*/

#ifndef __SCC_HH__
#define __SCC_HH__

#include <forward_list>
#include <stack>
#include <iostream>

namespace liblsap {

  template <typename IndexType>
  struct ConnectedComponent {
    typedef typename std::forward_list<IndexType> NodeList;
    IndexType index;
    NodeList *lNodes;
    IndexType nbNodes;
    IndexType weight;
    ConnectedComponent *next;
    ConnectedComponent()
      : index(0), lNodes(new NodeList()), nbNodes(0),
	weight(std::numeric_limits<IndexType>::min()), next(NULL) { }
    ConnectedComponent(const IndexType &id)
      : index(id), lNodes(new NodeList()), nbNodes(0),
	weight(std::numeric_limits<IndexType>::min()), next(NULL) { }
    ~ConnectedComponent() { if (lNodes) delete lNodes; }
    // void clear() { lNodes->clear(); nbNodes = 0; weight = std::numeric_limits<IndexType>::min(); }
  };
  
  // ==============================================================
  template <typename IndexType>
  class GraphComponent : public std::vector<IndexType> {
  public:
    typedef typename std::vector<IndexType> NodeList;
    IndexType index;    
    GraphComponent()
      : std::vector<IndexType>(), index(0) { }
    GraphComponent(const IndexType &id)
      : std::vector<IndexType>(), index(id) { }
    ~GraphComponent() { }
    void print() {
      std::cout << "idx=" << index << ": ";
      for (typename NodeList::const_iterator it = this->cbegin(), cend = this->cend(); it != cend; ++it)
	std::cout << *it << " ";
      std::cout << std::endl;
    }
  };

  // ==============================================================
  template <typename IndexType, class ComponentType = GraphComponent<IndexType>,
	    class ComponentsContainer = std::vector<ComponentType*> >
  class SccDecomp : public ComponentsContainer {
  public:
    typedef ComponentType Component;
    IndexType num;
    IndexType *vnum;
    IndexType *vaccess;
    bool *in_S;
    std::stack<IndexType> S;
    Component **node_comp;
    Component *cur_scc;
    IndexType nb_scc;
    // --------------------------------------------------------------
    SccDecomp() : num(0), vnum(NULL), vaccess(NULL), in_S(NULL), node_comp(NULL),
		  cur_scc(NULL), nb_scc(0) { }
    // --------------------------------------------------------------
    SccDecomp(IndexType nb_nodes)
      : num(0), vnum(new IndexType[nb_nodes]), vaccess(new IndexType[nb_nodes]),
	in_S(new bool[nb_nodes]), node_comp(new Component*[nb_nodes]),
	cur_scc(NULL), nb_scc(0) {
      for (ComponentType **c = node_comp, **cend = node_comp+nb_nodes; c != cend; c++) *c = NULL;
    }
    // --------------------------------------------------------------
    ~SccDecomp() { clear_containers(); }
    // --------------------------------------------------------------
    void clear_containers()
    {
      for (typename ComponentsContainer::iterator it = this->begin(), end = this->end();
	   it != end; ++it) { if (*it) { delete *it; *it = NULL; } }
      if (vnum) delete[] vnum;
      if (vaccess) delete[] vaccess;
      if (in_S) delete[] in_S;
      if (node_comp) delete[] node_comp;
      vnum = NULL; vaccess = NULL; in_S = NULL; node_comp = NULL;
      cur_scc = NULL; nb_scc = 0;
    }
    // --------------------------------------------------------------
    void init(IndexType nb_nodes)
    {
      clear_containers();
      vnum = new IndexType[nb_nodes];
      vaccess = new IndexType[nb_nodes];
      in_S = new bool[nb_nodes];
      node_comp = new Component*[nb_nodes];
      for (ComponentType **c = node_comp, **cend = node_comp+nb_nodes; c != cend; c++) *c = NULL;
      num = 0;	
    }
    // --------------------------------------------------------------
    Component* randComponent()
    {
      IndexType rd = randInt(1,this->nb_scc);
      for (typename ComponentsContainer::iterator it = this->begin(), end = this->end();
	   it != end; ++it) {
	if (*it) {
	  if (rd == 1) return *it;
	  rd--;
	}
      }
      return NULL;
    }
    // --------------------------------------------------------------
    void print() {
      for (typename ComponentsContainer::iterator it = this->begin(), end = this->end();
	   it != end; ++it) {
	if (*it) {
	  std::cout << "scc ";
	  (*it)->print();
	}
      }
      
    }
  };

  // ==============================================================
  /** @class StronglyConnectedComponents
   *  @brief Decomposition of a graph into SCC
   */
  template <class GraphType>
  class StronglyConnectedComponents {

  public:
    typedef typename GraphType::IndexType IndexType;
    typedef ConnectedComponent<IndexType> Component;
    typedef std::forward_list<Component*> ComponentList;
    typedef typename Component::NodeList NodeList;
    
  protected:
    typedef typename GraphType::Arc Arc;
    // for algo
    IndexType _num;
    IndexType *_vnum;
    IndexType *_vaccess;
    std::stack<IndexType> _stack;
    bool *_in_stack;
    // the graph
    GraphType &_G;
    // SCCs
    IndexType _nbSCC;
    //ComponentList _scc_list;
    Component *_scc_list_head;
    Component **_node_scc;
    IndexType _cur_scc_idx;
    bool _cur_scc_used;
    Component *_cur_scc;
    bool _already_decomposed;
    NodeList _unused_nodes;
    IndexType _nb_unused_nodes;
    IndexType _mx_idx;
    
    // -------------------------------------------------------------
    // Tarjan decomposition for bipartite graph from a node v (its key)
    void decomposeRec(IndexType v, IndexType v_last)
    {
      _vaccess[v] = _vnum[v] = _num;
      _stack.push(v); _in_stack[v] = true;
      _num++;
      
      for (const Arc *arc = _G.outgoingArc(_G.key2idx(v)); arc != NULL; arc = arc->next)
      {
	if (_vnum[_G.idx2key(arc->to->index)] == -1) // 1st access
	{
	  decomposeRec(_G.idx2key(arc->to->index), v);
	  _vaccess[v] = std::min(_vaccess[v],_vaccess[_G.idx2key(arc->to->index)]);
	}
	else
	  if (/*(_G.idx2key(arc->to->index) != v_last) &&*/ _in_stack[_G.idx2key(arc->to->index)])
	    _vaccess[v] = std::min(_vaccess[v],_vnum[_G.idx2key(arc->to->index)]);
      }
      
      if (_vaccess[v] == _vnum[v]) // v is a root, extract scc
      {
	IndexType w = 0;
	_nbSCC++;
	Component *scc = new Component(_nbSCC);

	do {
	  w = _stack.top();
	  _stack.pop();
	  _in_stack[w] = false;
	  scc->lNodes->push_front(_G.key2idx(w));
	  scc->nbNodes++;
	  _node_scc[w] = scc;
	} while (w != v);

	if (scc->nbNodes < 3) 
	{
	  for (typename NodeList::const_iterator nit = scc->lNodes->cbegin(),
		 nend = scc->lNodes->cend(); nit != nend; ++nit) {
	    _node_scc[_G.idx2key(*nit)] = NULL;
	    _unused_nodes.push_front(*nit);
	    _nb_unused_nodes++;
	  }
	  delete scc; scc = NULL; _nbSCC--;
	}
	else pushSCC(scc);
      }
    }
    // --------------------------------------------------------------
    void pushSCC(Component *scc)
    {
      if (_scc_list_head) {
	scc->next = _scc_list_head;
	_scc_list_head = scc;
      }
      else _scc_list_head = scc;
    }
    // --------------------------------------------------------------
    void popSCC()
    {
      if (_scc_list_head) _scc_list_head = _scc_list_head->next;
    }
    // --------------------------------------------------------------
    void connectSCC(Component *scc, Component *prev)
    {
      if (prev) { scc->next = prev->next; prev->next = scc; }
      else { scc->next = _scc_list_head; _scc_list_head = scc; }
    }
    // --------------------------------------------------------------
    void unconnectSCC(Component *scc, Component *prev)
    {
      if (prev) prev->next = scc->next;
      else _scc_list_head = scc->next;
      scc->next = NULL;
    }
    // --------------------------------------------------------------
    // decompose only current scc index (_cur_scc) from node v (its key)
    void updateDecompositionRec(IndexType v, IndexType v_last)
    {
      _vaccess[v] = _vnum[v] = _num;
      _stack.push(v); _in_stack[v] = true;
      _num++;
      
      for (Arc *arc = _G.outgoingArc(_G.key2idx(v)); arc != NULL; arc = arc->next)
      {
	if (_node_scc[_G.idx2key(arc->to->index)] == _cur_scc) { // restrict to same scc
	  if (_vnum[_G.idx2key(arc->to->index)] == -1) {
	    updateDecompositionRec(_G.idx2key(arc->to->index),v);
	    _vaccess[v] = std::min(_vaccess[v],_vaccess[_G.idx2key(arc->to->index)]);
	  }
	  else
	    if (/*(_G.idx2key(arc->to->index) != v_last) &&*/ _in_stack[_G.idx2key(arc->to->index)])
	      _vaccess[v] = std::min(_vaccess[v],_vnum[_G.idx2key(arc->to->index)]);
	}
      }

      if (_vaccess[v] == _vnum[v]) // v is a root
      {
	IndexType w = 0;
	Component *scc = NULL;
	
	if (_cur_scc_used) { _nbSCC++; scc = new Component(_nbSCC); }
	else { scc = _cur_scc; }

	do {
	  w = _stack.top();
	  _stack.pop();
	  _in_stack[w] = false;
	  scc->lNodes->push_front(_G.key2idx(w));
	  scc->nbNodes++;
	  _node_scc[w] = scc;
	} while (w != v);
	
	if (scc->nbNodes < 3) // not an SCC
	{
	  for (typename NodeList::const_iterator nit = scc->lNodes->cbegin(),
		 nend = scc->lNodes->cend(); nit != nend; ++nit) {
	    _node_scc[_G.idx2key(*nit)] = NULL;
	    _unused_nodes.push_front(*nit);
	    _nb_unused_nodes++;
	  }
	  if (_cur_scc_used) { delete scc; scc = NULL; _nbSCC--; }
	  else { scc->lNodes->clear(); scc->nbNodes = 0; }
	}
	else { // an SCC
	  if (_cur_scc_used) pushSCC(scc);
	  else { _cur_scc_used = true; }
	  _mx_idx++; scc->index = _mx_idx;
	}
      }
    }
    
  public:
    // -----------------------------------------------------------------
    StronglyConnectedComponents(GraphType &G)
      : _num(0), _vnum(new IndexType[G.nbNodes()]), _vaccess(new IndexType[G.nbNodes()]),
	_in_stack(new bool[G.nbNodes()]), _G(G), _nbSCC(0), _scc_list_head(NULL),
	_node_scc(new Component*[G.nbNodes()]),
	_cur_scc_idx(0), _cur_scc_used(false), _already_decomposed(false),
	_nb_unused_nodes(0), _mx_idx(0)
    { for (int i = 0; i < _G.nbNodes(); i++) _node_scc[i] = NULL; }
    // -----------------------------------------------------------------
    ~StronglyConnectedComponents()
    {
      if (_in_stack) delete[] _in_stack;
      if (_vnum) delete[] _vnum;
      if (_vaccess) delete[] _vaccess;
      if (_node_scc) delete[] _node_scc;
      for (Component *c = _scc_list_head, *cnext = NULL; c != NULL; c = cnext) {
	cnext = c->next;
	delete c;
      }
      _scc_list_head = NULL;
    }
    // -----------------------------------------------------------------
    IndexType nbSCC() { return _nbSCC; }
    Component* head() { return _scc_list_head; }
    typename NodeList::const_iterator unused_nodes_cbegin() { return _unused_nodes.cbegin(); }
    typename NodeList::const_iterator unused_nodes_cend() { return _unused_nodes.cend(); }
    // -----------------------------------------------------------------
    IndexType sccIndex(const IndexType &n_index)
    {
      const Component* cc = _node_scc[_G.idx2key(n_index)];
      return (cc == NULL ? 0 : cc->index);
    }
    // -----------------------------------------------------------------
    std::pair<Component*,Component*> sccFromIndex(IndexType scc_index)
    {
      Component* cc = head(), *cc_prev = NULL;
      for (; cc != NULL; cc = cc->next) {
	if (cc->index == scc_index) return std::make_pair(cc,cc_prev);
	cc_prev = cc;
      }
      return std::make_pair((Component*)NULL,(Component*)NULL);
    }
    // -----------------------------------------------------------------
    std::pair<Component*,Component*> randComponent()
    {
      if (_nbSCC == 1) return std::make_pair(head(),(Component*)NULL);
      const IndexType r = randInt(1,_nbSCC);
      IndexType i = 1;
      for (Component* cc = head(), *cc_prev = NULL; cc != NULL; cc = cc->next) {
	if (r == i) return std::make_pair(cc,cc_prev);
	cc_prev = cc;
	i++;
      }
      return std::make_pair((Component*)NULL,(Component*)NULL);
    }
    // -----------------------------------------------------------------
    std::pair<Component*,Component*> sccFromIndex(IndexType scc_index, Component** acomp)
    {
      Component* cc = head(), *cc_prev = NULL;
      std::pair<Component*,Component*> res = std::make_pair((Component*)NULL,(Component*)NULL);
      
      for (; cc != NULL; cc = cc->next) {
	acomp[cc->index] = cc;
	if (scc_index == cc->index) {
	  res = std::make_pair(cc,cc_prev);
	  cc = cc->next;
	  break;
	}
	cc_prev = cc;
      }
      for (; cc != NULL; cc = cc->next) acomp[cc->index] = cc;
      return res;
    }
    // -----------------------------------------------------------------
    Component* scc(const IndexType &n_index) { return _node_scc[_G.idx2key(n_index)]; }
    const Component* scc(const IndexType &n_index) const { return _node_scc[_G.idx2key(n_index)]; }
    // -----------------------------------------------------------------
    void clear()
    {
      _already_decomposed = false;
      for (IndexType i = 0; i < _G.nbNodes(); i++) _node_scc[i] = -1;
    }
    // ------------------------------------------------------------------------
    void decompose()
    {
      if (!_already_decomposed) {
	_num = 0; _nbSCC = 0;
	IndexType *it = _vnum;
	bool *its = _in_stack;
	for (IndexType *end = _vnum + _G.nbNodes(); it != end; it++, its++)
	  { *it = -1; *its = false; }
	it = _vnum;
	for (IndexType n_key = 0; n_key < _G.nbNodes(); n_key++, it++)
	  if (*it == -1) decomposeRec(n_key,n_key);
	_already_decomposed = true;
	_mx_idx = _nbSCC;
      }
    }
    // ------------------------------------------------------------------------
    /*std::pair<IndexType,IndexType> updateDecomposition(Component* cur_scc,
						       Component* prev_scc)
    {
      if (_already_decomposed)
      {
	IndexType n_key, nb_scc = _nbSCC;
	NodeList *nlist = cur_scc->lNodes;
	_num = 0;
	_cur_scc = cur_scc;
	_cur_scc_used = false;
	cur_scc->lNodes = new NodeList();
	cur_scc->nbNodes = 0;
	_nb_unused_nodes = 0;
	
	typename NodeList::const_iterator it = nlist->cbegin(), end = nlist->cend();
	
	for (; it != end; ++it) _vnum[_G.idx2key(*it)] = -1;
	
	for (it = nlist->cbegin(); it != end; ++it) {
	  n_key = _G.idx2key(*it);
	  if (_vnum[n_key] == -1) updateDecompositionRec(n_key,n_key);
	}
	
	delete nlist;

	if (!_cur_scc_used) // no SCC
	{
	  _nbSCC--;
	  unconnectSCC(cur_scc,prev_scc);
	}
	return std::make_pair(_nbSCC-nb_scc,_nb_unused_nodes);
      }
      return std::make_pair((IndexType)0,(IndexType)0);
      }*/
    // ------------------------------------------------------------------------
    std::pair<IndexType,IndexType> updateDecompAndTrim(Component* cur_scc,
						       Component* prev_scc,
						       std::forward_list<Arc*> &Larcs)
    {
      if (_already_decomposed)
      {
	IndexType n_key, nb_scc = _nbSCC;
	NodeList *nlist = cur_scc->lNodes;
	_num = 0;
	_cur_scc = cur_scc;
	_cur_scc_used = false;
	cur_scc->lNodes = new NodeList;
	cur_scc->nbNodes = 0;
	_nb_unused_nodes = 0;
	
	typename NodeList::const_iterator it = nlist->cbegin(), end = nlist->cend();
	
	for (; it != end; ++it) _vnum[_G.idx2key(*it)] = -1;
	
	for (it = nlist->cbegin(); it != end; ++it) {
	  n_key = _G.idx2key(*it);
	  if (_vnum[n_key] == -1) updateDecompositionRec(n_key,n_key);
	}

	if (!_cur_scc_used) // no SCC
	{
	  _nbSCC--;
	  unconnectSCC(cur_scc,prev_scc);
	}

	trimArcs(*nlist,Larcs);
	
	delete nlist; nlist = NULL;
	
	return std::make_pair(_nbSCC-nb_scc,_nb_unused_nodes);
      }
      return std::make_pair<IndexType>(0,0);
    }
    // ------------------------------------------------------------------------
    // recompose a decomposition from a given scc (no exploration of the arcs is needed)
    void recompose(IndexType scc_idx, Component* cur_scc, Component* prev_scc,
		   IndexType nb_other_scc, IndexType nb_unused_nodes)
    {
      Component *other_scc = NULL;
      if (cur_scc->nbNodes == 0) // cur_scc not in the list
	{ _nbSCC++; connectSCC(cur_scc,prev_scc); }
      cur_scc->index = scc_idx;
      for (; nb_other_scc > 0; nb_other_scc--)
      {
	other_scc = _scc_list_head;
	for (typename NodeList::const_iterator it = other_scc->lNodes->cbegin(),
	       end = other_scc->lNodes->cend(); it != end; ++it)
	{
	  cur_scc->lNodes->push_front(*it);
	  cur_scc->nbNodes++;
	  _node_scc[_G.idx2key(*it)] = cur_scc;
	}
	popSCC();
	delete other_scc; other_scc = NULL; _mx_idx--;
	_nbSCC--;
      }
      for (; nb_unused_nodes > 0; nb_unused_nodes--)
      {
	cur_scc->lNodes->push_front(_unused_nodes.front());
	cur_scc->nbNodes++;
	_node_scc[_G.idx2key(_unused_nodes.front())] = cur_scc;
	_unused_nodes.pop_front();
      }
    }
    // --------------------------------------------------------------
    // Unconnect arcs (still in memory) not in a strongly connected component
    // and starting at a node of V2, and return them as a list
    // in O(|E|)
    void trimArcs(std::forward_list<Arc*> &L)
    {
      for (IndexType n = 0, mx = _G.maxNodeIdx()+1; n != mx; n++) // nodes of V2
      {
	Component *cc = _node_scc[_G.idx2key(n)];
	if (cc == NULL) {  // not in an SCC, unconnect arcs
	  for (Arc *arc = _G.outgoingArc(n), *a_next = NULL; arc != NULL; arc = a_next)
	    { a_next = arc->next; _G.unconnect(arc); L.push_front(arc); }
	}
	else { // in an SCC, unconnect arcs to a node in a different SCC
	  for (Arc *arc = _G.outgoingArc(n), *a_next = NULL; arc != NULL; arc = a_next) {
	    a_next = arc->next;
	    if (_node_scc[_G.idx2key(arc->to->index)] != cc)
	      { _G.unconnect(arc); L.push_front(arc); }
	  }
	}
      }
    }
    // --------------------------------------------------------------
    void trimArcs(Component *cc, IndexType nb_other, IndexType nb_unused_nodes,
		  std::forward_list<Arc*> &Larcs)
    {
      Arc* arc = NULL, *a_next = NULL;
      for (typename NodeList::const_iterator nit = cc->lNodes->cbegin(),
	     nend = cc->lNodes->cend(); nit != nend; ++nit)
      {
	if (_G.inS2(*nit))
	  for (arc = _G.outgoingArc(*nit); arc != NULL; arc = a_next)
	  {
	    a_next = arc->next;
	    if (_node_scc[_G.idx2key(arc->to->index)] != cc) {
	      _G.unconnect(arc); Larcs.push_front(arc);
	    }
	  }
      }

      for (Component *cur_scc = head(); nb_other > 0 && cur_scc;
	   cur_scc = cur_scc->next, nb_other--)
      {
	for (typename Component::NodeList::const_iterator nit = cur_scc->lNodes->cbegin(),
	       nend = cur_scc->lNodes->cend(); nit != nend; ++nit)
	{
	  if (_G.inS2(*nit))
	    for (arc = _G.outgoingArc(*nit); arc != NULL; arc = a_next)
	    {
	      a_next = arc->next;
	      if (scc(arc->to->index) != cur_scc)
	      {
		_G.unconnect(arc);
		Larcs.push_front(arc);
	      }
	    }
	}
      }
      
      // arcs in NULL scc and in parent SCC
      if (nb_unused_nodes > 0) {
	for (typename Component::NodeList::const_iterator nit = unused_nodes_cbegin(),
	       nend = unused_nodes_cend();
	     nit != nend && nb_unused_nodes > 0; ++nit, nb_unused_nodes--)
	{
	  if (_G.inS2(*nit)) {
	    for (arc = _G.outgoingArc(*nit); arc != NULL; arc = a_next)
	    {
	      a_next = arc->next;
	      _G.unconnect(arc);
	      Larcs.push_front(arc);
	    }
	  }
	}
      }
    }
    // --------------------------------------------------------------
    void trimArcs(NodeList &L, std::forward_list<Arc*> &Larcs)
    {
      Arc* arc = NULL, *a_next = NULL;
      for (typename NodeList::const_iterator nit = L.cbegin(), nend = L.cend();
	   nit != nend; ++nit)
      {
	Component* cmp = _node_scc[_G.idx2key(*nit)];
	if (_G.inS2(*nit)) {
	  if (cmp == NULL) {  // isolated node
	    if (!_G.isDummyS2(*nit)) {
	      for (arc = _G.outgoingArc(*nit); arc != NULL; arc = a_next)
		{ a_next = arc->next; _G.unconnect(arc); Larcs.push_front(arc); }
	    }
	    else {
	      for (arc = _G.outgoingArc(*nit); arc != NULL; arc = a_next) {
		a_next = arc->next;
		if (!_G.isDummyS1(arc->to->index)) { _G.unconnect(arc); Larcs.push_front(arc); }
	      }
	    }
	  }
	  else { // not an isolated node
	    if (!_G.isDummyS2(*nit)) {
	      for (arc = _G.outgoingArc(*nit); arc != NULL; arc = a_next) {
		a_next = arc->next;
		if (_node_scc[_G.idx2key(arc->to->index)] != cmp)
		  { _G.unconnect(arc); Larcs.push_front(arc); }
	      }
	    }
	    else {
	      for (arc = _G.outgoingArc(*nit); arc != NULL; arc = a_next) {
		a_next = arc->next;
		if (!_G.isDummyS1(arc->to->index) && _node_scc[_G.idx2key(arc->to->index)] != cmp)
		  { _G.unconnect(arc); Larcs.push_front(arc); }
	      }
	    }
	  }
	}
	else { // in S1
	  if (cmp == NULL) { // meaning not anymore in an SCC (isolated node)
	    if (!_G.isDummyS1(*nit)) {
	      for (arc = _G.ingoingArc(*nit); arc != NULL; arc = a_next) {
		a_next = arc->next_to;
		if (_node_scc[_G.idx2key(arc->from->index)] != cmp)
		  { _G.unconnect(arc); Larcs.push_front(arc); }
	      }
	    }
	    else {
	      for (arc = _G.ingoingArc(*nit); arc != NULL; arc = a_next) {
		a_next = arc->next_to;
		if (!_G.isDummyS2(arc->from->index) && _node_scc[_G.idx2key(arc->from->index)] != cmp)
		  { _G.unconnect(arc); Larcs.push_front(arc); }
	      }
	    }
	  }
	}
      }
    }
    // ------------------------------------------------------------------------
    void print()
    {
      std::cout << std::endl << _nbSCC << " SCCs\n";
      for (Component *c = _scc_list_head; c != NULL; c = c->next)
      {
	std::cout << "scc index=" << " " << c->index << " w=" << c->weight << " nodes=[";
	for (typename NodeList::const_iterator itn = c->lNodes->cbegin(),
	       endn = c->lNodes->cend(); itn != endn; ++itn) std:: cout << *itn << ",";
	std::cout << "]" << std::endl;
      }
      
      for (IndexType i = _G.minNodeIdx(); i <= _G.maxNodeIdx(); i++)
	std::cout << i << ':' << sccIndex(i) << ", ";
      std::cout << std::endl << std::endl;
    }
  };
}

#endif
