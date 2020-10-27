// -----------------------------------------------------------
/**
 * @file graphs.hh
 * @brief Data structures for graphs, bipartite graphs
 * @author Sebastien Bougleux
 * @date September 07 2019
 * @institution Normandie Univ, UNICAEN, ENSICAEN, CNRS, GREYC, France
 */
/*
 * -----------------------------------------------------------
 * This file is part of LSAPE.
 * LSAPE is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL-C License. See README for more
 * details.
 * -----------------------------------------------------------
*/

#ifndef __GP_GRAPHS_H__
#define __GP_GRAPHS_H__

#include <typeinfo>
#include <iomanip>
#include <limits>
#include <forward_list>
#include <vector>
#include <queue>
#include <stack>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "matchings.hh"
#include "enumerators.hh"

namespace liblsap {
  
  // ==============================================================
  enum GRAPH_TYPE { UNDIRECTED = 0, DIRECTED, BIDIRECTED, ORIENTED };
  enum GRAPH_DUMMY_TYPE { DUMMY_0 = 0, DUMMY_12, DUMMY_2, DUMMY_1 };

  // ==============================================================
  /** @class IndexNode
   *  @brief Basic node data structure with index
   */
  template <typename IdxType = int>
  class IndexNode {

  public:
    typedef IdxType IndexType;
    IndexType index;  // index in graph array

    IndexNode() : index(0) { }
    IndexNode(IndexType id) : index(id) { }
    ~IndexNode() { }
    const IndexType& idx() const { return index; }
    IndexType deg() const { return 0; }
    IndexType degIn() const { return 0; }
    IndexType degOut() const { return 0; }
    void degOutInc() { }
    void degOutDec() { }
    void degInInc() { }
    void degInDec() { }
  };

  // ==============================================================
  /** @class IndexNode
   *  @brief Basic node data structure with index
   */
  template <typename IdxType = int>
  class DegreeNode : public IndexNode<IdxType> {
  public:
    typedef IdxType IndexType;
    IndexType deg_out; // degree in graph
    IndexType deg_in;

    DegreeNode() : IndexNode<IdxType>(), deg_out(0), deg_in(0) { }
    DegreeNode(IndexType id) : IndexNode<IdxType>(id), deg_out(0), deg_in(0) { }
    ~DegreeNode() { }
    IndexType deg() const { return (deg_in + deg_out); }
    IndexType degIn() const { return deg_in; }
    IndexType degOut() const { return deg_out; }
    void degOutInc() { deg_out++; }
    void degOutDec() { deg_out--; }
    void degInInc() { deg_in++; }
    void degInDec() { deg_in--; }
  };

  // ==============================================================
  /** @class NeibArc
   *  @brief Arc data structure
   */
  template <class NodeType = IndexNode<> >
  struct NeibArc {
    typedef NodeType Node;
    
    Node *from;               // start node
    Node *to;                 // end node
    NeibArc<Node> *next;      // next arc incident to from node
    NeibArc<Node> *prev;      // previous arc around start node
    NeibArc<Node> *next_to;   // next arc around end node
    NeibArc<Node> *prev_to;   // previous arc around end node
    
    NeibArc()
      : from(NULL), to(NULL), next(NULL), prev(NULL), next_to(NULL), prev_to(NULL) { }
    NeibArc(Node* const &n_from, Node* const &n_to)
      : from(n_from), to(n_to), next(NULL), prev(NULL), next_to(NULL), prev_to(NULL) { }
    ~NeibArc() { }
    void print() { std::cout << from->index << "->" << to->index; }
  };

  // ==============================================================
  /** @class WeightedArc
   *  @brief Weighted arc data structure, similar to NeibArc
   */
  template <class WeightType = double, class NodeType = IndexNode<> >
  struct WeightedArc {
    typedef NodeType Node;
    typedef WeightType Weight;
    
    Node *from;
    Node *to;
    WeightedArc<WeightType,Node> *next;
    WeightedArc<WeightType,Node> *prev;
    WeightedArc<WeightType,Node> *next_to;
    WeightedArc<WeightType,Node> *prev_to;
    WeightType weight;
    
    WeightedArc()
      : from(NULL), to(NULL), next(NULL), prev(NULL), next_to(NULL), prev_to(NULL), weight(0) { }
    WeightedArc(Node* const &n_from, Node * const &n_to)
      : from(n_from), to(n_to), next(NULL), prev(NULL), next_to(NULL), prev_to(NULL), weight(0) { }
    ~WeightedArc() { }
    void print() { std::cout << from->index << "->" << to->index << " w=" << weight; }
  };

  // ==============================================================
  /** @class NeibGraph
   *  @brief Basic graph data structure based on a list of neighbors for each node
   */
  template <class ArcType = NeibArc<>,
	    class NodeContainer = std::vector<typename ArcType::Node*>,
	    class ArcContainer = std::vector<ArcType*> >
  class NeibGraph {

  public:
    typedef NeibGraph<ArcType,NodeContainer,ArcContainer> GraphType;
    typedef ArcType Arc;
    typedef typename Arc::Node Node;
    typedef typename Node::IndexType IndexType;
    typedef NodeContainer NodeList;
    typedef ArcContainer ArcList;
    typedef std::list<Arc*> Cycle;
    
  protected:
    IndexType _nb_nodes;
    IndexType _nb_arcs;
    GRAPH_TYPE _gtype;
    GRAPH_DUMMY_TYPE _dummy_type;
    NodeList _nodes;
    ArcList _outgoing;
    ArcList _ingoing;

    // --------------------------------------------------------------
    void connectOutgoing(Arc &arc)
    {
      Arc* &inc = outgoingArc(arc.from->index);
      if (inc == NULL) { arc.prev = NULL; arc.next = NULL; }
      else { arc.next = inc; inc->prev = &arc; }
      inc = &arc;
      arc.from->degOutInc();
      arc.to->degInInc();
      _nb_arcs++;
    }
    // --------------------------------------------------------------
    void connectIngoing(Arc &arc)
    {
      Arc* &inc = ingoingArc(arc.to->index);
      arc.prev_to = NULL;
      if (inc == NULL) arc.next_to = NULL;
      else { arc.next_to = inc; inc->prev_to = &arc; }
      inc = &arc;
    }
    // --------------------------------------------------------------
    void unconnectOutgoing(Arc *arc)
    {
      Arc* &inc = outgoingArc(arc->from->index);
      if (inc == arc) {
	inc = arc->next;
	if (arc->next) arc->next->prev = NULL;
      }
      else {
	arc->prev->next = arc->next;
	if (arc->next) arc->next->prev = arc->prev;
      }
      arc->from->degOutDec();
      arc->to->degInDec();
      _nb_arcs--;
    }
    // --------------------------------------------------------------
    void unconnectIngoing(Arc *arc)
    {
      Arc* &inc = ingoingArc(arc->to->index);
      if (inc == arc) {
	inc = arc->next_to;
	if (arc->next_to) arc->next_to->prev_to = NULL;
      }
      else {
	arc->prev_to->next_to = arc->next_to;
	if (arc->next_to) arc->next_to->prev_to = arc->prev_to;
      }
    }
    

  public:
    NeibGraph(GRAPH_TYPE gtype = UNDIRECTED)
      : _nb_nodes(0), _nb_arcs(0), _gtype(gtype), _dummy_type(DUMMY_0) { }
    // --------------------------------------------------------------
    virtual ~NeibGraph() { clear(); }
    // --------------------------------------------------------------
    virtual NeibGraph<Arc,NodeContainer,ArcContainer>& clear() { deleteNodes(); return *this; }
    // --------------------------------------------------------------
    const IndexType& nbNodes() const { return _nb_nodes; }
    // --------------------------------------------------------------
    virtual IndexType nbRealNodes() { return _nb_nodes; }
    // --------------------------------------------------------------
    virtual IndexType nbDummyNodes() const { return 0; }
    // --------------------------------------------------------------
    const IndexType& nbArcs() const { return _nb_arcs; }
    // --------------------------------------------------------------
    virtual IndexType key2idx(const IndexType &n_key) const { return n_key; }
    virtual IndexType idx2key(const IndexType &n_index) const { return n_index; }
    virtual IndexType idx2key(const IndexType &n_index) { return n_index; }
    // --------------------------------------------------------------
    virtual bool validNodeIndex(const IndexType &n_index)
    { return (n_index >= 0 && n_index < nbNodes()); }
    // --------------------------------------------------------------
    virtual IndexType minNodeIdx() const { return 0; }
    // --------------------------------------------------------------
    virtual IndexType maxNodeIdx() const { return nbNodes()-1; }
    // --------------------------------------------------------------
    virtual Node* node(const IndexType& n_index) { return _nodes[n_index]; }
    // --------------------------------------------------------------
    bool isolatedNode(const IndexType &n_index) { return (outgoingArc(n_index) == NULL); }
    // --------------------------------------------------------------
    virtual bool isDummy(const IndexType &n_index) const { return false; }
    virtual bool isDummy(const IndexType &n_index) { return false; }
    // --------------------------------------------------------------
    virtual Node* insertNode(unsigned short always_0_not_used = 0)
    {
      Node *n = new Node(nbNodes());
      _nodes.push_back(n);
      _outgoing.push_back(NULL);
      _ingoing.push_back(NULL);
      _nb_nodes++;
      return n;
    }
    // --------------------------------------------------------------
    void insertNodes(IndexType nb_nodes, unsigned short always_0_not_used = 0)
    { for (; nb_nodes > 0; nb_nodes--) insertNode(always_0_not_used); }
    // --------------------------------------------------------------
    /*virtual void deleteNode(const IndexType &n_index)
    {
      deleteArcsFrom(n_index);
      delete _nodes[n_index]; _nodes[n_index] = NULL;
      _nb_nodes--;  // TODO: decay indices
      }*/
    // --------------------------------------------------------------
    void deleteNodes()
    {
      deleteArcs();
      _outgoing.clear();
      _ingoing.clear();
      for (typename NodeList::iterator it = _nodes.begin(), end = _nodes.end(); it != end; ++it)
	{ delete *it; *it = NULL; }
      _nodes.clear();
      _nb_nodes = 0;
    }
    // --------------------------------------------------------------
    virtual Arc* arc()
    {
      Arc *res = NULL;
      for (typename NodeList::const_iterator nit = _nodes.cbegin(), nend = _nodes.cend();
	   nit != nend; ++nit)
      {
	res = outgoingArc((*nit)->index);
	if (res != NULL) return res;
      }
      return NULL;
    }    
    // --------------------------------------------------------------
    Arc*& outgoingArc(const IndexType &n_from) { return _outgoing[idx2key(n_from)]; }
    // --------------------------------------------------------------
    const Arc* outgoingArc(const IndexType &n_from) const { return _outgoing[idx2key(n_from)]; }
    // --------------------------------------------------------------
    Arc*& ingoingArc(const IndexType &n_index) { return _ingoing[idx2key(n_index)]; }
    // --------------------------------------------------------------
    const Arc* ingoingArc(const IndexType &n_index) const { return _ingoing[idx2key(n_index)]; }
    // --------------------------------------------------------------
    Arc* findArc(const IndexType &n_from, const IndexType &n_to)
    {
      for (Arc *arc = outgoingArc(n_from); arc != NULL; arc = arc->next)
	if (arc->to->index == n_to) return arc;
      return NULL;
    }
    // --------------------------------------------------------------
    Arc* findArc(const IndexType &n_from, const IndexType &n_to, Arc* &prev_arc)
    {
      for (Arc *arc = outgoingArc(n_from); arc != NULL; arc = arc->next)
      {
	if (arc->to->index == n_to) return arc;
	prev_arc = arc;
      }
      return NULL;
    }
    // --------------------------------------------------------------
    // insert an arc as the 1st outgoing arc of n_from
    Arc* insertArc(const IndexType &n_from, const IndexType &n_to)
    {
      Arc* arc = NULL;
      if (validNodeIndex(n_from) && validNodeIndex(n_to)) {
	arc = new Arc(node(n_from),node(n_to));
	Arc* &a_inc = outgoingArc(n_from);
	if (a_inc) { arc->next = a_inc; a_inc->prev = arc; }
	a_inc = arc;
	arc->from->degOutInc();
	arc->to->degInInc();
	_nb_arcs++;
	Arc* &a_inc2 = ingoingArc(n_to);
	if (a_inc2) { arc->next_to = a_inc2; a_inc2->prev_to = arc; }
	a_inc2 = arc;
      }
      return arc;
    }
    // --------------------------------------------------------------
    void connect(Arc *arc) { connectOutgoing(*arc); connectIngoing(*arc); }
    // --------------------------------------------------------------
    void unconnect(Arc *arc) { unconnectOutgoing(arc); unconnectIngoing(arc); }
    // --------------------------------------------------------------
    template <class ArcsContainer>
    void unconnect_and_trim(Arc *arc, ArcsContainer &L)
    {
      std::queue<Arc*> q;
      Arc *ca = NULL;
      IndexType nfrom, nto, nbeg = arc->from->index;
      q.push(arc);
      while (!q.empty())
      {
	ca = q.front(); q.pop();
	nfrom = ca->from->index; nto = ca->to->index;
	if (!isDummy(nfrom) || !isDummy(nto)) { unconnect(ca); L.push_front(ca); }
	ca = ingoingArc(nto); 
	if (ca == NULL)
	  for (Arc *na = outgoingArc(nto); na != NULL; na = na->next) q.push(na);
      }
      if (outgoingArc(nbeg) == NULL) {
	for (Arc *na = ingoingArc(nbeg); na != NULL; na = na->next_to) q.push(na);
	while (!q.empty())
	{
	  ca = q.front(); q.pop();
	  nfrom = ca->from->index; nto = ca->to->index;
	  if (!isDummy(nfrom) || !isDummy(nto)) { unconnect(ca); L.push_front(ca); }
	  ca = outgoingArc(nfrom); 
	  if (ca == NULL)
	    for (Arc *na = ingoingArc(nfrom); na != NULL; na = na->next_to) q.push(na);
	}
      }
      //std::cout << "not NULL beg node " << nbeg << std::endl;
      // dummy nodes may be still connected
      // remove if
      /*if (_dummy_type == DUMMY_12) {
	nto = -1;
	ca = ingoingArc(nto);
	if (ca && ca->next_to == NULL && isDummy(ca->from->index)) {
	  nto = 0;
	  ca = ingoingArc(nto);
	  if (ca && ca->next_to == NULL && isDummy(ca->from->index)) {
	    for (Arc *na = outgoingArc(nto); na != NULL; na = na->next) q.push_back(na);
	  }
	}	
      }
      */
    }
    // --------------------------------------------------------------
    template <class ArcsContainer>
    void unconnect_outarcs_and_trim(Arc *arc, ArcsContainer &L)
    {
      std::queue<Arc*> q;
      Arc *ca = NULL;
      IndexType nfrom, nto;
      q.push(arc);
      while (!q.empty())
      {
	ca = q.front(); q.pop();
	nfrom = ca->from->index; nto = ca->to->index;
	if (!isDummy(nfrom) || !isDummy(nto)) { unconnect(ca); L.push_front(ca); }
	ca = ingoingArc(nto); 
	if (ca == NULL)
	  for (Arc *na = outgoingArc(nto); na != NULL; na = na->next) q.push(na);
      }
      // dummy nodes may be still connected
      // remove if
      /*if (_dummy_type == DUMMY_12) {
	nto = -1;
	ca = ingoingArc(nto);
	if (ca && ca->next_to == NULL && isDummy(ca->from->index)) {
	  nto = 0;
	  ca = ingoingArc(nto);
	  if (ca && ca->next_to == NULL && isDummy(ca->from->index)) {
	    for (Arc *na = outgoingArc(nto); na != NULL; na = na->next) q.push_back(na);
	  }
	}	
      }
      */
    }
    // --------------------------------------------------------------
    void connectOutgoingArcs(Arc *arc)
    {
      outgoingArc(arc->from->index) = arc;
      for (; arc != NULL; arc = arc->next) { connectIngoing(*arc); _nb_arcs++; }
    }
    // --------------------------------------------------------------
    void connectIngoingArcs(Arc *in_arc)
    {
      ingoingArc(in_arc->to->index) = in_arc;
      for (; in_arc != NULL; in_arc = in_arc->next_to) connectOutgoing(*in_arc);
    }
    // --------------------------------------------------------------
    Arc* unconnectOutgoingArcs(Node &n)
    {
      Arc* rec = outgoingArc(n.index);
      for (Arc *arc = rec; arc != NULL; arc = arc->next) { unconnectIngoing(arc); _nb_arcs--; }
      outgoingArc(n.index) = NULL;
      return rec;
    }
    // --------------------------------------------------------------
    Arc* unconnectIngoingArcs(Node &n)
    {
      Arc* rec = ingoingArc(n.index);
      for (Arc *arc = rec; arc != NULL; arc = arc->next_to) unconnectOutgoing(arc);
      ingoingArc(n.index) = NULL;
      return rec;
    }
    // --------------------------------------------------------------
    // delete arc from node index n_from to node index n_to
    void deleteArc(const IndexType &n_from, const IndexType &n_to)
    {
      if (n_from < nbNodes() && n_to < nbNodes()) {
	Arc *prev_arc = NULL;
	Arc *arc = findArc(n_from,n_to,prev_arc);
	if (arc)
	{
	  if (prev_arc) {
	    prev_arc->next = arc->next;
	    if (arc->next) arc->next->prev = prev_arc;
	  }
	  else outgoingArc(n_from) = NULL;
	  if (arc->prev_to) {
	    arc->prev_to->next_to = arc->next_to;
	    if (arc->next_to) arc->next_to->prev_to = arc->prev_to;
	  }
	  else ingoingArc(n_to) = NULL;
	  arc->from->degOutDec();
	  arc->to->degInDec();
	  delete arc;
	  _nb_arcs--;
	}
      }
    }
    // --------------------------------------------------------------
    // delete all arcs incident to the node of index n_from
    void deleteArcs(const IndexType &n_from)
    {
      if (n_from >= minNodeIdx() && n_from <= maxNodeIdx()) {
	for (Arc *arc = outgoingArc(n_from), *a_next = NULL; arc != NULL; arc = a_next)
	{
	  a_next = arc->next;
	  unconnect(arc);
	  delete arc; arc = NULL;
	}
	outgoingArc(n_from) = NULL;
      }
    }
    // --------------------------------------------------------------
    // delete all arcs
    void deleteArcs()
    {
      if (_nb_arcs > 0)
	for (typename NodeList::const_iterator it = _nodes.cbegin(), end = _nodes.cend();
	     it != end; ++it) deleteArcs((*it)->index);
    }
    // --------------------------------------------------------------
    void swap(Arc *arc)
    {
      unconnect(arc);
      Node *n_tmp = arc->from;
      arc->from = arc->to;
      arc->to = n_tmp;
      connect(arc);
    }
    // --------------------------------------------------------------
    void swap(Cycle &cycle)
    {
      for (typename Cycle::iterator it = cycle.begin(), end = cycle.end();
	   it != end; ++it)
	if (!isDummy((*it)->from->index) || !isDummy((*it)->to->index)) swap(*it);
    }
    // --------------------------------------------------------------
    template <class GT>
    void trimArcs(StronglyConnectedComponents<GT>& scc)
    {
      IndexType scc_id = 0;
      for (IndexType n = 0; n < nbNodes(); n++)
      {
	scc_id = scc[key2idx(n)];
	for (Arc *arc = outgoingArc(key2idx(n)), *a_next = NULL; arc != NULL; arc = a_next)
	{
	  a_next = arc->next;
	  if (scc[arc->to->index] != scc_id) // not in same SCC
	  {
	    unconnect(arc);
	    delete arc;
	  }
	}
      }
    }
    // --------------------------------------------------------------
    void unconnect(SccDecomp<IndexType> &s_scc, std::forward_list<Arc*> &L)
    {
      IndexType vid = 0;
      typename SccDecomp<IndexType>::Component *scc = NULL;
      for (IndexType n = 0; n < nbNodes(); n++)
      {
	scc = s_scc.node_comp[n];
	vid = key2idx(n);
	if (scc == NULL) { // isolated node
	  if (!isDummy(vid)) {
	    for (Arc *arc = outgoingArc(vid), *a_next = NULL; arc != NULL; arc = a_next)
	      { a_next = arc->next; unconnect(arc); L.push_front(arc); }
	  }
	  else {
	    for (Arc *arc = outgoingArc(vid), *a_next = NULL; arc != NULL; arc = a_next) {
	      a_next = arc->next;
	      if (isDummy(arc->to->index)) continue;
	      unconnect(arc); L.push_front(arc);
	    }
	  }
	}
	else { // not an isolated node
	  if (!isDummy(vid)) {
	    for (Arc *arc = outgoingArc(vid), *a_next = NULL; arc != NULL; arc = a_next)
	    {
	      a_next = arc->next;
	      if (s_scc.node_comp[idx2key(arc->to->index)] != scc) // not in same SCC
		{ unconnect(arc); L.push_front(arc); }
	    }
	  }
	  else {
	    for (Arc *arc = outgoingArc(vid), *a_next = NULL; arc != NULL; arc = a_next)
	    {
	      a_next = arc->next;
	      if (isDummy(arc->to->index)) continue;
	      if (s_scc.node_comp[idx2key(arc->to->index)] != scc) // not in same SCC
		{ unconnect(arc); L.push_front(arc); }
	    }
	  }
	}
      }
    }
    // --------------------------------------------------------------
    bool dfsCycle(const IndexType &n_index_from, const IndexType &n_index, const IndexType &n_prev,
		  bool *visited, Cycle &cycle)
    {
      if (!isDummy(n_index)) visited[idx2key(n_index)] = true;
      for (Arc *arc = outgoingArc(n_index); arc != NULL; arc = arc->next)
      {
	// avoid 2-cycles (occur only for dummy nodes)
	if (!cycle.empty() && arc->to->index == cycle.back()->from->index) continue;
	if (arc->to->index == n_index_from) { cycle.push_back(arc); return true; }
	if (visited[idx2key(arc->to->index)] == false)
	{
	  cycle.push_back(arc);
	  if (dfsCycle(n_index_from,arc->to->index,n_index,visited,cycle)) return true;
	  cycle.pop_back();
	}
      }
      return false;
    }
    // --------------------------------------------------------------
    bool dfsCycle(const IndexType n_index, Cycle &cycle)
    {
      bool *visited = new bool[nbNodes()], cycle_found = false;
      for (bool *v = visited, *v_end = visited+nbNodes(); v != v_end; ++v) *v = false;
      cycle_found = dfsCycle(n_index,n_index,n_index,visited,cycle);
      delete[] visited;
      return cycle_found;
    }
    // --------------------------------------------------------------
    bool loadFromTextFile(const char *file_name)
    {
      std::ifstream infile(file_name);
      if (!infile) {
	// TODO: error msg
	return false;
      }
      IndexType nb_nodes, nb_arcs, n_from, n_to;
      clear();
      std::cout << "read file " << file_name << std::endl;
      infile >> nb_nodes >> nb_arcs;
      insertNodes(nb_nodes);
      std::cout << "#nodes=" << nbNodes() << std::endl;
      for (IndexType i = 0; i < nb_arcs; i++) {
	infile >> n_from >> n_to;;
	insertArc(n_from,n_to);
      }
      infile.close();
      std::cout << "#arcs=" << nbArcs() << std::endl;
      std::cout << "file closed" << std::endl;
      return true;
    }
    // --------------------------------------------------------------
    void print()
    {
      std::cout << "print graph\n" << "#nodes=" << nbNodes() << " #arcs=" << nbArcs() << std::endl;
      IndexType n = 0;
      Arc *arc = NULL;
      for (typename ArcList::iterator it = _outgoing.begin(), end = _outgoing.end();
	   it != end; ++it, n++)
      {
	arc = *it;
	if (arc == NULL) std::cout << "no outgoing arc for node " << n << std::endl;
	else
	  for (; arc != NULL; arc = arc->next)
	    std::cout << arc->from->index << " -> " << arc->to->index << std::endl;
      }
      std::cout << "ingoing structure\n";
      for (typename ArcList::iterator it = _ingoing.begin(), end = _ingoing.end();
	   it != end; ++it, n++)
      {
	arc = *it;
	if (arc == NULL) std::cout << "no ingoing arc for node " << n << std::endl;
	else
	  for (; arc != NULL; arc = arc->next_to)
	    std::cout << arc->from->index << " -> " << arc->to->index << std::endl;
      }
    }
  };

  // ==============================================================
  template <class ArcType = NeibArc<> >
  class BipartiteGraph : public NeibGraph<ArcType,
					  std::deque<typename ArcType::Node*>,
					  std::deque<ArcType*> > {

  public:
    typedef ArcType Arc;
    typedef typename Arc::Node Node;
    typedef NeibGraph<ArcType,
		      std::deque<Node*>,
		      std::deque<Arc*> > ParentGraph;
    typedef BipartiteGraph<ArcType> GraphType;
    typedef typename ParentGraph::IndexType IndexType;
    typedef typename ParentGraph::NodeList NodeList;
    typedef typename ParentGraph::ArcList ArcList;
    typedef typename ParentGraph::Cycle Cycle;    
    using ParentGraph::nbNodes;
    using ParentGraph::outgoingArc;
    using ParentGraph::ingoingArc;
    typedef Matching<IndexType> NodeMapType;
    typedef std::vector<NodeMapType> NodeMapsContainer;
    
  protected:
    IndexType _nb_nodes_S12[2];
    
  public:
    BipartiteGraph(GRAPH_TYPE gtype = DIRECTED) : ParentGraph(gtype)
    { _nb_nodes_S12[0] = _nb_nodes_S12[1] = 0; }
    BipartiteGraph(const IndexType &n1, const IndexType &n2, GRAPH_TYPE gtype = DIRECTED)
      : ParentGraph(gtype)
    {
      _nb_nodes_S12[0] = _nb_nodes_S12[1] = 0;
      ParentGraph::insertNodes(n1,0); ParentGraph::insertNodes(n2,1);
    }
    // --------------------------------------------------------------
    ~BipartiteGraph() { deleteNodes(); }
    // --------------------------------------------------------------
    virtual bool isBalanced() { return (_nb_nodes_S12[0] == _nb_nodes_S12[1]); }
    // --------------------------------------------------------------
    const IndexType& nbNodes(const unsigned short &set_index) { return _nb_nodes_S12[set_index]; }
    // --------------------------------------------------------------
    const IndexType& nbNodes(const unsigned short &set_index) const
    { return _nb_nodes_S12[set_index]; }
    // --------------------------------------------------------------
    virtual IndexType nbRealNodes(const unsigned short &set_index) 
    { return _nb_nodes_S12[set_index]; }
    // --------------------------------------------------------------
    Node* node(const IndexType &n_index) { return ParentGraph::node(idx2key(n_index)); }
    // --------------------------------------------------------------
    //virtual bool isDummyS1(const IndexType &n_index) { return false; }
    virtual bool isDummyS1(const IndexType &n_index) const { return false; }
    //virtual bool isDummyS2(const IndexType &n_index) { return false; }
    virtual bool isDummyS2(const IndexType &n_index) const { return false; }
    virtual IndexType dummyS1idx() const { return std::numeric_limits<IndexType>::min(); }
    virtual IndexType dummyS2idx() const { return std::numeric_limits<IndexType>::max(); }
    // --------------------------------------------------------------
    bool inS1(const IndexType &n_index) { return (n_index < 0); }
    bool inS2(const IndexType &n_index) { return (n_index >= 0); }
    // --------------------------------------------------------------
    typename NodeList::const_iterator nodes_S1_cbegin()
    { return ParentGraph::_nodes.cbegin(); }
    typename NodeList::const_iterator nodes_S1_cend()
    { return ParentGraph::_nodes.cbegin()+_nb_nodes_S12[0]; }
    // --------------------------------------------------------------
    bool validNodeIndex(const IndexType &n_index)
    { return (n_index >= -_nb_nodes_S12[0] && n_index < _nb_nodes_S12[1]); }
    // --------------------------------------------------------------
    IndexType minNodeIdx() const { return -_nb_nodes_S12[0]; }
    // --------------------------------------------------------------
    IndexType maxNodeIdx() const { return _nb_nodes_S12[1]-1; }
    // --------------------------------------------------------------
    IndexType minRealNodeIdxS1() const { return -_nb_nodes_S12[0]; }
    virtual IndexType minRealNodeIdxS2() const { return 0; }
    // --------------------------------------------------------------
    virtual IndexType maxRealNodeIdxS1() const { return -1; }
    IndexType maxRealNodeIdxS2() const { return _nb_nodes_S12[1]-1; }
    // --------------------------------------------------------------
    virtual IndexType idx2key(const IndexType &n_index) { return (n_index + _nb_nodes_S12[0]); }
    virtual IndexType idx2key(const IndexType &n_index) const
    { return (n_index + _nb_nodes_S12[0]); }
    // --------------------------------------------------------------
    IndexType key2idx(const IndexType &n_key) const
    { return (n_key - _nb_nodes_S12[0]); }
    // --------------------------------------------------------------
    virtual IndexType idx2setidx(const IndexType &n_idx) const
    { return (n_idx < 0 ? -n_idx-1 : n_idx); }
    // --------------------------------------------------------------
    virtual IndexType setidx2idx(unsigned short set_idx, const IndexType &n_setidx) const
    { return (set_idx == 0 ? -n_setidx-1 : n_setidx); }
    // --------------------------------------------------------------
    virtual IndexType idxg2idx(const IndexType &n_idxg) const
    { return (n_idxg < nbNodes(0) ? -n_idxg-1 : n_idxg-nbNodes(0)); }
    // --------------------------------------------------------------
    Node* insertNode(unsigned short set_index = 0)
    {
      Node *n = NULL;
      if (set_index == 0) {
	_nb_nodes_S12[set_index]++;
	n = new Node(-_nb_nodes_S12[set_index]);
	ParentGraph::_nodes.push_front(n);
	ParentGraph::_outgoing.push_front(NULL);
	ParentGraph::_ingoing.push_front(NULL);
      }
      else {
	n = new Node(_nb_nodes_S12[set_index]);
	_nb_nodes_S12[set_index]++;
	ParentGraph::_nodes.push_back(n);
	ParentGraph::_outgoing.push_back(NULL);
	ParentGraph::_ingoing.push_back(NULL);
      }
      ParentGraph::_nb_nodes++;
      return n;
    }
    // --------------------------------------------------------------
    void deleteNodes() { ParentGraph::deleteNodes(); _nb_nodes_S12[0] = _nb_nodes_S12[1] = 0; }
    // --------------------------------------------------------------
    virtual Arc* arc()
    {
      Arc *res = NULL;
      for (IndexType i = minRealNodeIdxS1(), iend = maxRealNodeIdxS1(); i <= iend; i++)
      {
	res = outgoingArc(i);
	if (res != NULL) return res;
      }
      return NULL;
    }    
    // --------------------------------------------------------------
    void scc_decomp(IndexType v, IndexType v_last, SccDecomp<IndexType> &s_scc)
    {
      s_scc.vaccess[v] = s_scc.vnum[v] = s_scc.num;
      s_scc.S.push(v); s_scc.in_S[v] = true;
      s_scc.num++;
      
      for (const Arc *arc = outgoingArc(key2idx(v)); arc != NULL; arc = arc->next)
      {
	if (s_scc.vnum[idx2key(arc->to->index)] == -1) // 1st access
	{
	  scc_decomp(idx2key(arc->to->index),v,s_scc);
	  s_scc.vaccess[v] = std::min(s_scc.vaccess[v],s_scc.vaccess[idx2key(arc->to->index)]);
	}
	else
	  if ((idx2key(arc->to->index) != v_last) && s_scc.in_S[idx2key(arc->to->index)])
	    s_scc.vaccess[v] = std::min(s_scc.vaccess[v],s_scc.vnum[idx2key(arc->to->index)]);
      }
      
      if (s_scc.vaccess[v] == s_scc.vnum[v]) // v is a root, extract scc
      {
	IndexType w = 0;
	typename SccDecomp<IndexType>::Component *scc =
	  new typename SccDecomp<IndexType>::Component(s_scc.size());
	
	do {
	  w = s_scc.S.top();
	  s_scc.S.pop();
	  s_scc.in_S[w] = false;
	  scc->push_back(key2idx(w));
	  s_scc.node_comp[w] = scc;
	} while (w != v);

	// avoid 2-cycles (in case of dummy nodes or unoriented graphs)
	if (scc->size() > 2) { s_scc.push_back(scc); s_scc.nb_scc++; }
	else {
	  for (auto nit = scc->cbegin(), nend = scc->cend(); nit != nend; ++nit)
	    s_scc.node_comp[idx2key(*nit)] = NULL;
	  delete scc; scc = NULL;
	}
      }
    }
    // --------------------------------------------------------------
    void sccDecomp(SccDecomp<IndexType> &s_scc)
    {
      // init
      s_scc.init(nbNodes());
      IndexType *it = s_scc.vnum;
      bool *its = s_scc.in_S;
      for (IndexType *end = s_scc.vnum + nbNodes(); it != end; it++, its++) { *it = -1; *its = false; }

      // decompose
      it = s_scc.vnum;
      for (IndexType n_key = 0; n_key < nbNodes(); n_key++, it++)
	if (*it == -1) scc_decomp(n_key,n_key,s_scc);
    }
    // --------------------------------------------------------------
    // decompose only current scc index (_cur_scc) from node v (its key)
    void scc_update_decomp(IndexType v, IndexType v_last, SccDecomp<IndexType> &s_scc)
    {
      s_scc.vaccess[v] = s_scc.vnum[v] = s_scc.num;
      s_scc.S.push(v); s_scc.in_S[v] = true;
      s_scc.num++;
      
      for (Arc *arc = outgoingArc(key2idx(v)); arc != NULL; arc = arc->next)
      {
	if (s_scc.node_comp[idx2key(arc->to->index)] == s_scc.cur_scc) { // restrict to same scc
	  if (s_scc.vnum[idx2key(arc->to->index)] == -1) {
	    scc_update_decomp(idx2key(arc->to->index),v,s_scc);
	    s_scc.vaccess[v] = std::min(s_scc.vaccess[v],s_scc.vaccess[idx2key(arc->to->index)]);
	  }
	  else
	    if ((idx2key(arc->to->index) != v_last) && s_scc.in_S[idx2key(arc->to->index)])
	      s_scc.vaccess[v] = std::min(s_scc.vaccess[v],s_scc.vnum[idx2key(arc->to->index)]);
	}
      }
      
      if (s_scc.vaccess[v] == s_scc.vnum[v]) // v is a root
      {
	IndexType w = 0;
	typename SccDecomp<IndexType>::Component *scc =
	  new typename SccDecomp<IndexType>::Component(s_scc.size());
	
	do {
	  w = s_scc.S.top();
	  s_scc.S.pop();
	  s_scc.in_S[w] = false;
	  scc->push_back(key2idx(w));
	  s_scc.node_comp[w] = scc;
	} while (w != v);
	
	if (scc->size() > 2) { s_scc.push_back(scc); s_scc.nb_scc++; }
	else { // not an SCC
	  for (auto nit = scc->cbegin(), nend = scc->cend(); nit != nend; ++nit)
	    s_scc.node_comp[idx2key(*nit)] = NULL;
	  delete scc; scc = NULL;
	}
      }
    }
    // --------------------------------------------------------------
    void sccUpdateDecompAndTrim(IndexType scc_id, SccDecomp<IndexType> &s_scc, IndexType &nb_scc,
				std::forward_list<Arc*> &L, std::forward_list<IndexType> &Lnodes)
    {
      // init
      typename SccDecomp<IndexType>::Component *cur_scc = s_scc[scc_id], *new_scc = NULL;
      Arc *arc = NULL, *a_next = NULL;
      IndexType nb_scc_init = s_scc.size(), n_key;
      s_scc.num = 0;
      s_scc.cur_scc = cur_scc;

      typename SccDecomp<IndexType>::Component::const_iterator it = cur_scc->cbegin(),
	end = cur_scc->cend();
      for (; it != end; ++it) s_scc.vnum[idx2key(*it)] = -1;
      
      // update decomposition
      for (it = cur_scc->cbegin(); it != end; ++it) {
	n_key = idx2key(*it);
	if (s_scc.vnum[n_key] == -1) scc_update_decomp(n_key,n_key,s_scc);
      }

      // trim arcs
      for (it = cur_scc->cbegin(); it != end; ++it) {
	new_scc = s_scc.node_comp[idx2key(*it)];
	if (new_scc == NULL) { // isolated node
	  Lnodes.push_front(*it);
	  if (!this->isDummy(*it)) {
	    for (arc = outgoingArc(*it), a_next = NULL; arc != NULL; arc = a_next)
	      { a_next = arc->next; this->unconnect(arc); L.push_front(arc); }
	    for (arc = ingoingArc(*it), a_next = NULL; arc != NULL; arc = a_next)
	      { a_next = arc->next_to; this->unconnect(arc); L.push_front(arc); }
	  }
	  else {
	    for (arc = outgoingArc(*it), a_next = NULL; arc != NULL; arc = a_next) {
	      a_next = arc->next;
	      if (this->isDummy(arc->to->index)) continue;
	      this->unconnect(arc); L.push_front(arc);
	    }
	    for (arc = ingoingArc(*it), a_next = NULL; arc != NULL; arc = a_next) {
	      a_next = arc->next_to;
	      if (this->isDummy(arc->from->index)) continue;
	      this->unconnect(arc); L.push_front(arc);
	    }
	  }
	}
	else {
	  if (!this->isDummy(*it)) {
	    for (arc = outgoingArc(*it), a_next = NULL; arc != NULL; arc = a_next) {
	      a_next = arc->next;
	      if (s_scc.node_comp[idx2key(arc->to->index)] != new_scc) // not in same SCC
		{ this->unconnect(arc); L.push_front(arc); }
	    }
	    for (arc = ingoingArc(*it), a_next = NULL; arc != NULL; arc = a_next) { //TODO check utility
	      a_next = arc->next_to;
	      if (s_scc.node_comp[idx2key(arc->from->index)] != new_scc) // not in same SCC
		{ this->unconnect(arc); L.push_front(arc); }
	    }
	  }
	  else {
	    for (arc = outgoingArc(*it), a_next = NULL; arc != NULL; arc = a_next) {
	      a_next = arc->next;
	      if (this->isDummy(arc->to->index)) continue;
	      if (s_scc.node_comp[idx2key(arc->to->index)] != new_scc) // not in same SCC
		{ this->unconnect(arc); L.push_front(arc); }
	    }
	    for (arc = ingoingArc(*it), a_next = NULL; arc != NULL; arc = a_next) { //TODO check utility
	      a_next = arc->next_to;
	      if (this->isDummy(arc->from->index)) continue;
	      if (s_scc.node_comp[idx2key(arc->from->index)] != new_scc) // not in same SCC
		{ this->unconnect(arc); L.push_front(arc); }
	    }
	  }
	}
      }
      delete cur_scc; cur_scc = NULL; s_scc[scc_id] = NULL; // to be reconstructed later if needed
      nb_scc = s_scc.size()-nb_scc_init;
      s_scc.cur_scc = NULL;
      s_scc.nb_scc--;
    }
    // ------------------------------------------------------------------------
    // recompose a decomposition from a given scc (no exploration of the arcs is needed)
    void sccRecomp(IndexType scc_id, SccDecomp<IndexType> &s_scc, 
		   IndexType nb_scc, std::forward_list<IndexType> &Lnodes)
    {
      typename SccDecomp<IndexType>::Component *other_scc = NULL, *cur_scc = NULL;
      if (nb_scc == 0) cur_scc = new typename SccDecomp<IndexType>::Component(scc_id);
      else {
	cur_scc = s_scc.back(); s_scc.pop_back(); cur_scc->index = scc_id; nb_scc--; s_scc.nb_scc--;
	for (typename SccDecomp<IndexType>::Component::const_iterator it = cur_scc->cbegin(),
	       end = cur_scc->cend(); it != end; ++it)
	  s_scc.node_comp[idx2key(*it)] = cur_scc;
	for (; nb_scc > 0; nb_scc--)
	{
	  other_scc = s_scc.back(); s_scc.pop_back(); s_scc.nb_scc--;
	  for (typename SccDecomp<IndexType>::Component::const_iterator it = other_scc->cbegin(),
		 end = other_scc->cend(); it != end; ++it)
	  {
	    cur_scc->push_back(*it);
	    s_scc.node_comp[idx2key(*it)] = cur_scc;
	  }
	  delete other_scc; other_scc = NULL;
	}
      }
      while (!Lnodes.empty())
      {
	cur_scc->push_back(Lnodes.front());
	s_scc.node_comp[idx2key(Lnodes.front())] = cur_scc;
	Lnodes.pop_front();
      }
      s_scc[scc_id] = cur_scc; s_scc.nb_scc++;
    }
    // --------------------------------------------------------------
    /*Arc* randArc()
    {
      const IndexType rd = -randInt(1,nbNodes(0));
      Arc *res = outgoingArc(rd);
      if (res && !(ParentGraph::isDummy(res->from->index) || ParentGraph::isDummy(res->to->index)))
	return res;
      for (IndexType i = rd-1; i >= minNodeIdx(); i--) {
	res = outgoingArc(i);
	if (res && !(ParentGraph::isDummy(res->from->index) || ParentGraph::isDummy(res->to->index)))
	  return res;
      }
      for (IndexType i = rd+1; i <= maxNodeIdx(); i++) {
	res = outgoingArc(i);
	if (res && !(ParentGraph::isDummy(res->from->index) || ParentGraph::isDummy(res->to->index))) {
	  if (i >= 0)
	    std::cout << res->from->index << "->" << res->to->index << std::endl;
	  return res;
	}
      }
      return NULL;
      }*/
    // --------------------------------------------------------------
    void sym_diff(NodeMapType &M, Cycle &cycle)
    {
      typename Cycle::iterator it = cycle.begin(), end = cycle.end();
      if (inS1((*it)->from->index)) ++it;
      for (; it != end; ++it) {
	if (!(isDummyS2((*it)->from->index) && isDummyS1((*it)->to->index))) {
	  if (!isDummyS1((*it)->to->index))
	  {
	    M(0,idx2setidx((*it)->to->index)) = idx2setidx((*it)->from->index);
	    if (!isDummyS2((*it)->from->index)) {
	      M(1,idx2setidx((*it)->from->index)) = idx2setidx((*it)->to->index);
	    }
	  }
	  else M(1,idx2setidx((*it)->from->index)) = idx2setidx((*it)->to->index);
	}
	++it;
	if (it == end) break;
      }
    }
    // --------------------------------------------------------------
    /* void enum_matchings(IndexType &nb_mx_match, Matching<IndexType> &M,
			NodeMaps<IndexType> &matchs)
    {
      // only isolated nodes => M is a leaf for recursion tree
      if (ParentGraph::nbArcs() == 0 || nb_mx_match == 0) return;

      // find cycle and arc (1st arc of the cycle)
      Arc *arc = em_select_arc_naive();
      if (arc == NULL) return;
      
      Cycle cycle;
      enum_matchings_find_cycle(arc,cycle);

      if (cycle.empty()) return;
      
      // save matching
      matchs.push_back(new Matching<IndexType>(M));
      sym_diff(matchs.nodeMap(matchs.size()-1),cycle);
      nb_mx_match--;

      // enumerate other matchings
      std::forward_list<Arc*> save_arcs;
	
      if (inS1(cycle.front()->from->index))
      {
	// swap arcs along cycle and unconnect unecessary arcs and save them
	ParentGraph::swap(cycle);
	ParentGraph::unconnect_and_trim(cycle.front(),save_arcs);

	// enumerate all matchings without arc cycle.front()
	enum_matchings(nb_mx_match,matchs.nodeMap(matchs.size()-1),matchs);
	
	// reconstruct graph
	while (!save_arcs.empty()) { ParentGraph::connect(save_arcs.front()); save_arcs.pop_front(); }
	ParentGraph::swap(cycle);
	
	// unconnect unecessary arcs and save them
	ParentGraph::unconnect_and_trim(cycle.front(),save_arcs);
	// enumerate all matchings with arc cycle.front()
	enum_matchings(nb_mx_match,M,matchs);

	// reconnect arcs and free space
	while (!save_arcs.empty()) { ParentGraph::connect(save_arcs.front()); save_arcs.pop_front(); }
      }
      else { // cycle's 1st node is in S2
	std::cout << "other side !!!" << std::endl;
	// unconnect unecessary arcs and save them
	unconnect_and_trim<std::forward_list<Arc*> >(cycle.front(),save_arcs);

	// enumerate all matchings without arc cycle.front()
	enum_matchings(nb_mx_match,M,matchs);
	
	// reconstruct graph
	while (!save_arcs.empty()) { connect(save_arcs.front()); save_arcs.pop_front(); }
	
	// swap arcs along cycle and unconnect unecessary arcs and save them
	swap(cycle);
	unconnect(cycle.front());
	save_arcs.push_front(cycle.front());
	if (!isDummyS1(cycle.front()->from->index))
	  unconnect_inarcs_and_trim(*(cycle.front()->from),save_arcs);
	if (!isDummyS2(cycle.front()->to->index))
	  unconnect_outarcs_and_trim(*(cycle.front()->to),save_arcs);
	
	// enumerate all matchings with arc cycle.front()
	enum_matchings(nb_mx_match,Mp,matchs);

	// reconstruct graph and free space
	while (!save_arcs.empty()) { connect(save_arcs.front()); save_arcs.pop_front(); }
	swap(cycle);
      }
    }*/
    // --------------------------------------------------------------
    /*void enumMatchings(Matching<IndexType> &cur_match, IndexType nb_mx_match,
		       NodeMaps<IndexType> &L, ENUM_EDG_SELECT select_edg = ENUM_EDG_SELECT_RAND)
    {
      EnumMatch<GraphType,Matching<IndexType>,NodeMaps<IndexType> > E(*this);
      E.enumerate(cur_match,nb_mx_match,L,select_edg);
      }*/
    // --------------------------------------------------------------
    template <class MatchContainer>
    void enumMaximumMatchings(Matching<IndexType> &match_init, unsigned int nb_match,
			      MatchContainer &L, ENUM_EDG_SELECT option = ENUM_EDG_SELECT_BALANCED)
    {
      EnumMatchings<BipartiteGraph<Arc>,Matching<IndexType>,MatchContainer> E(*this);
      E.enumerate(match_init,nb_match,L,option);
    }
    // --------------------------------------------------------------
    template <class MatchContainer>
    void enumDissimilarMaximumMatchings(Matching<IndexType> &match_init, unsigned int nb_match,
					MatchContainer &L,
					ENUM_DISS_ALGO algo = ENUM_DISS_MXW_DFS,
					ENUM_EDG_SELECT edge_select = ENUM_EDG_SELECT_RAND,
					ENUM_CYCLE_SELECT cycle_select_opt = ENUM_CYCLE_SELECT_MXW, 
					bool several_scc = true)
    {
      EnumMatchings<BipartiteGraph<Arc>,Matching<IndexType>,MatchContainer> E(*this);
      E.enumerateDissimilar(match_init,nb_match,L,algo,edge_select,cycle_select_opt,several_scc);
    }
    // --------------------------------------------------------------
    bool loadFromTextFile(const char *file_name)
    {
      std::ifstream infile(file_name);
      if (!infile) {
	// TODO: error msg
	return false;
      }
      IndexType nb_nodes1, nb_nodes2, nb_arcs, n_from, n_to, n_from_index, n_to_index;
      this->clear();
      std::cout << "read file " << file_name << std::endl;
      infile >> nb_nodes1 >> nb_nodes2 >> nb_arcs;
      ParentGraph::insertNodes(nb_nodes1,0);
      ParentGraph::insertNodes(nb_nodes2,1);
      std::cout << "#nodes1=" << nbNodes(0) << " #nodes2=" << nbNodes(1) << std::endl;
      for (IndexType i = 0; i < nb_arcs; i++) {
	infile >> n_from >> n_to;
	n_from_index = idxg2idx(n_from);
	n_to_index = idxg2idx(n_to);
	ParentGraph::insertArc(n_from_index,n_to_index);
      }
      infile.close();
      std::cout << "#arcs=" << ParentGraph::nbArcs() << std::endl;
      std::cout << "file closed" << std::endl;
      return true;
    }
    // --------------------------------------------------------------
    void print(bool print_ingoing = false)
    {
      std::cout << "print bipartite graph\n" << "#nodes1=" << nbNodes(0)
		<< " #nodes2=" << nbNodes(1) << " #arcs=" << this->nbArcs() << std::endl;
      IndexType n = 0;
      Arc *arc = NULL;
      /*for (typename NodeList::const_iterator it = ParentGraph::_nodes.cbegin(),
	     end = ParentGraph::_nodes.cend(); it != end; ++it)
	     std::cout << (*it)->index << " " << (*it)->deg() << std::endl;*/
      for (typename ArcList::iterator it = this->_outgoing.begin(),
	     end = this->_outgoing.end();
	   it != end; ++it, n++)
      {
	arc = *it;
	if (arc == NULL) std::cout << "no outgoing arc for node " << key2idx(n) << std::endl;
	else
	  for (; arc != NULL; arc = arc->next)
	    {
	      arc->print();
	      if (this->isDummy(arc->from->index) || this->isDummy(arc->to->index))
		std::cout << " dummy"; std::cout << std::endl;
	    }
      }
      if (print_ingoing) {
	std::cout << "ingoing structure\n";
	n = 0;
	for (typename ArcList::iterator it = this->_ingoing.begin(),
	       end = this->_ingoing.end();
	     it != end; ++it, n++)
	{
	  arc = *it;
	  if (arc == NULL) std::cout << "no ingoing arc for node " << key2idx(n) << std::endl;
	  else
	    for (; arc != NULL; arc = arc->next_to) { arc->print(); std::cout << std::endl; }
	}
      }
    }
  };

  // ==============================================================
  // always two nodes, the dummy ones at index -1 for S1 and 0 for S2
  template <class ArcType = NeibArc<> >
  class BipartiteGraphEC : public BipartiteGraph<ArcType> {

  public:
    typedef ArcType Arc;
    typedef typename Arc::Node Node;
    typedef BipartiteGraph<Arc> ParentGraph;
    typedef typename ParentGraph::IndexType IndexType;
    typedef typename ParentGraph::NodeList NodeList;
    typedef typename ParentGraph::Cycle Cycle;    

  
  protected:
    
    // --------------------------------------------------------------
    void insertDummy()
    {
      if (ParentGraph::nbNodes() == 0) { ParentGraph::insertNode(0); ParentGraph::insertNode(1); }
      ParentGraph::insertArc(-1,0);
      ParentGraph::insertArc(0,-1);
    }
        
  public:  
    BipartiteGraphEC(GRAPH_TYPE gtype = DIRECTED)
      : BipartiteGraph<Arc>(gtype)
    { insertDummy(); ParentGraph::_dummy_type = DUMMY_12; }
    // --------------------------------------------------------------
    BipartiteGraphEC(const IndexType &n1, const IndexType &n2, GRAPH_TYPE gtype = DIRECTED)
      : BipartiteGraph<Arc>(gtype)
    {
      insertDummy();
      ParentGraph::insertNodes(n1,0);
      ParentGraph::insertNodes(n2,1);
      ParentGraph::_dummy_type = DUMMY_12;
    }
    // --------------------------------------------------------------
    ~BipartiteGraphEC() { this->deleteNodes(); }
    // --------------------------------------------------------------
    BipartiteGraphEC<Arc>& clear()
    {
      if (ParentGraph::nbNodes(0) > 1 || ParentGraph::nbNodes(1) > 1)
	{ParentGraph::clear(); insertDummy(); }
      return *this;
    }
    // --------------------------------------------------------------
    IndexType nbRealNodes() { return ParentGraph::_nb_nodes-2; }
    // --------------------------------------------------------------
    IndexType nbRealNodes(const unsigned short &set_index)
    { return ParentGraph::_nb_nodes_S12[set_index]-1; }
    // --------------------------------------------------------------
    IndexType nbRealNodes(const unsigned short &set_index) const
    { return ParentGraph::_nb_nodes_S12[set_index]-1; }
    // --------------------------------------------------------------
    IndexType nbDummyNodes() const { return 2; }
    // --------------------------------------------------------------
    bool isDummy(const IndexType &n_index) { return (isDummyS1(n_index) || isDummyS2(n_index)); }
    bool isDummy(const IndexType &n_index) const
    { return (isDummyS1(n_index) || isDummyS2(n_index)); }
    bool isDummyS1(const IndexType &n_index) { return (n_index == -1); }
    bool isDummyS2(const IndexType &n_index) { return (n_index == 0); }
    bool isDummyS1(const IndexType &n_index) const { return (n_index == -1); }
    bool isDummyS2(const IndexType &n_index) const { return (n_index == 0); }
    IndexType dummyS1idx() const { return -1; }
    IndexType dummyS2idx() const { return 0; }
    // --------------------------------------------------------------
    IndexType maxRealNodeIdxS1() const { return -2; }
    IndexType minRealNodeIdxS2() const { return 1; }
    // --------------------------------------------------------------
    bool isBalanced() { return true; }
    // --------------------------------------------------------------
    IndexType idx2setidx(const IndexType &n_idx) const
    { return (n_idx < -1 ? -n_idx-2 :
	      (n_idx == -1 ? ParentGraph::nbNodes(0)-1 :
	       (n_idx > 0 ? n_idx-1 : ParentGraph::nbNodes(1)-1))); }
    // --------------------------------------------------------------
    IndexType setidx2idx(unsigned short set_idx, const IndexType &n_setidx) const
    { return (set_idx == 0 ?
	      (n_setidx < nbRealNodes(set_idx) ? -(n_setidx+2) : -1) :
	      (n_setidx < nbRealNodes(set_idx) ? n_setidx+1 : 0)); }
    // --------------------------------------------------------------
    IndexType setidx2idx(unsigned short set_idx, const IndexType &n_setidx)
    { return (set_idx == 0 ?
	      (n_setidx < nbRealNodes(set_idx) ? -(n_setidx+2) : -1) :
	      (n_setidx < nbRealNodes(set_idx) ? n_setidx+1 : 0)); }
    // --------------------------------------------------------------
    IndexType idxg2idx(const IndexType &n_idxg) const
    { return (n_idxg+1 < ParentGraph::nbNodes(0) ? -n_idxg-2 : n_idxg+2-ParentGraph::nbNodes(0)); }
    // --------------------------------------------------------------
    Arc* arc()
    {
      Arc *res = NULL;
      for (IndexType i = this->minRealNodeIdxS1(), iend = dummyS1idx(); i != iend; i++)
      {
	res = ParentGraph::outgoingArc(i);
	if (res != NULL) return res;
      }
      return NULL;
    }    
    // --------------------------------------------------------------
    template <class MatchContainer = std::list<Matching<IndexType> > >
    void enumMaximumMatchings(Matching<IndexType> &match_init, unsigned int nb_match,
			      MatchContainer &L, ENUM_EDG_SELECT option = ENUM_EDG_SELECT_RAND)
    {
      EnumMatchings<BipartiteGraphEC<Arc>,Matching<IndexType>,MatchContainer> E(*this);
      E.enumerate(match_init,nb_match,L,option);
    }
    // --------------------------------------------------------------
    void generateRandom(IndexType n1, IndexType n2, IndexType dmax, IndexType insdel_percent = 20)
    {
      std::vector<IndexType> v;
      
      this->clear();
      ParentGraph::insertNodes(n1-1,0);
      ParentGraph::insertNodes(n2-1,1);

      std::srand(unsigned(std::time(0)));
	
      if (n1 < n2) {
      	IndexType n_idx = minRealNodeIdxS2(), n_end = maxRealNodeIdxS1(),
	  deg = 1, n_idx1, n_end2 = this->maxRealNodeIdxS2(), dcpt;

	// permutation of V2 (which is larger than V1)
	for (; n_idx <= n_end2; n_idx++) v.push_back(n_idx);
	std::random_shuffle(v.begin(),v.end());
	
	// insert arcs from V1 to V2 (all nodes of V1 are covered by a matching)
	typename std::vector<IndexType>::const_iterator it = v.cbegin(), end = v.cend();
	IndexType cent = 100;
	IndexType nbsub = nbRealNodes(0) * (100-insdel_percent)/cent, nbs = 0;
	// substitutions
	for (n_idx = this->minRealNodeIdxS1(); nbs < nbsub; n_idx++, ++it, nbs++)
	  ParentGraph::insertArc(n_idx,*it);
	// insertions + deletions
	for (; n_idx <= n_end; n_idx++, ++it) {
	  this->insertArc(dummyS1idx(),*it);
	  ParentGraph::insertArc(n_idx,dummyS2idx());
	}
	// other arcs
	v.clear();
	dmax++;
	if (dmax > nbRealNodes(0)) dmax = nbRealNodes(0);
	// save indicies of V1
	for (n_idx = this->minRealNodeIdxS1(); n_idx <= n_end; n_idx++) v.push_back(n_idx);
	std::random_shuffle(v.begin(),v.end());
	
	for (n_idx = this->minRealNodeIdxS2(), n_idx1 = 0; n_idx <= n_end2; n_idx++) {
	  // insert arcs from dummy node in V1 to uncovered nodes in V2
	  if (ParentGraph::ingoingArc(n_idx) == NULL) ParentGraph::insertArc(dummyS1idx(),n_idx);
	  // all nodes of V1 and V2 are covered
	  // insert arcs from V2 to V1 by iteratively permuting indicies of V1
	  // so that output degree of nodes in V2 is random between 0 and dmax
	  deg = std::rand()%dmax;
	  if (n_idx1+deg > n_end) { std::random_shuffle(v.begin(),v.end()); n_idx1 = 0; }
	  for (dcpt = 1; dcpt <= deg; dcpt++)
	    { ParentGraph::insertArc(n_idx,v[n_idx1]); n_idx1++; }
	}
      }
    }
    // --------------------------------------------------------------
    bool loadFromTextFile(const char *file_name)
    {
      std::ifstream infile(file_name);
      if (!infile) {
	// TODO: error msg
	return false;
      }
      IndexType nb_nodes1, nb_nodes2, nb_arcs, n_from, n_to, n_from_index, n_to_index;
      this->clear();
      std::cout << "read file " << file_name << std::endl;
      infile >> nb_nodes1 >> nb_nodes2 >> nb_arcs;
      ParentGraph::insertNodes(nb_nodes1-1,0);
      ParentGraph::insertNodes(nb_nodes2-1,1);
      std::cout << "#nodes1=" << ParentGraph::nbNodes(0) << " #nodes2=" << ParentGraph::nbNodes(1) << std::endl;
      for (IndexType i = 0; i < nb_arcs; i++) {
	infile >> n_from >> n_to;
	ParentGraph::insertArc(n_from,n_to);
      }
      infile.close();
      std::cout << "#arcs=" << ParentGraph::nbArcs() << std::endl;
      std::cout << "file closed" << std::endl;
      return true;
    }
    
    
  };
  
}

#endif
