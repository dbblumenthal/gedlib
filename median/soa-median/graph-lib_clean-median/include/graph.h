
#ifndef __PGRAPHH__
#define __PGRAPHH__

#include <vector>
#include <map>
#include <fstream>
#include <cstdlib>
#include <algorithm>    // for std::random_shuffle
#include <iostream>

#include <tinyxml.h>

#include "utils.h"

/** @brief An oriented edge of a graph.
 *
 * The class <code>GEdge</code> defines an oriented graph.
 * An edge points towards a Node identified by its id. Each edge is associated to a templated attribute.
 */

// GEdge is a class encoding edges within a graph. The label associated to each edge is templeted by EdgeAttribute
template<class EdgeAttribute>
class GEdge{
private:
  /** The next neighbourhood node. */
  GEdge *next;
  /** The rank of the node in the graph array. */
  int incident_node; //XXX : How is it difficult to link a GNode ?
  /** The number identifiing a given edge. */
   int edge_id;

public:
   /** Attribute of the edge */
  EdgeAttribute attr; //TODO : remove from public

   /**
    * Creates a new edge towards node with id n, attributed by attr and connecting to the GEdge list adj. Note that an edge is always oriented. Non oriented edges are encoded by two symmetric edges
    * @param n	the incident node.
    * @param adj the next edge (in a list).
    * @param attr the label.
    */
   GEdge( int n, GEdge *adj, EdgeAttribute attr ): next(adj), incident_node(n), edge_id(-1), attr(attr) {};


   /**
    * Creates a new edge with id i towards node with id n, with attributed by attr and connecting to the GEdge list adj. Note that an edge is always oriented. Non oriented edges are encoded by two symmetric edges
    * @param n	the incident node.
    * @param i the identifier of edge.
    * @param adj the next edge (in a list).
    * @param attr the label.
    */
   GEdge( int n, GEdge *adj, int i, EdgeAttribute attr): next(adj), incident_node(n), edge_id(i), attr(attr) {};

   /**
    * Deletes the edge
    */
   ~GEdge(){};

   /**
    * Returns the number of the incident incident_node.
    * @return	the number of the connected incident_node.
    */
   int IncidentNode() const { return incident_node; };

  /**
   * Sets the incident node for a given edge.
   * @param new_node	the number of the connected incident_node.
   * XXX: take care of symmetric edges
   */
  void  setIncidentNode(int new_node)  {  this->incident_node = new_node; };

   /**
    * Returns the next edge in the list of edge. Useful to traverse all incident edges from a node (const version).
    * @return	the next edge.
    */
  GEdge<EdgeAttribute>* Next() const { return next; };

   /**
    * Returns the next edge in the list of edge. Useful to traverse all incident edges from a node.
    * @return	the next edge.
    */
  GEdge<EdgeAttribute>* Next( GEdge<EdgeAttribute>* n ) { return next=n; };

   /**
    * Returns the index of the referenced edge.
    * @return	the index.
    */
   int EdgeId() const { return edge_id; }

   /**
    * Sets the new index of the referenced edge.
    * @param i	the new index.
    * @return	the new index of the object.
    */
   int EdgeId( int i ) { return edge_id=i; }
};

/** @brief A node of a graph.
 *
 * The class <code>GNode</code> defines a node.
 * A node corresponds to an element of a graph. Each Node is associated to an templated attribute and linked to other nodes through <code><GEdge/code>.
 */
template <class NodeAttribute, class EdgeAttribute>
class GNode{
private:
  /** The list of incident edges. */
  GEdge<EdgeAttribute> *adjacents;
  /** The id of the node in the graph. */
  int item;

public :
  /** The attribute of the node.
   * TODO: remove from public
   */
  NodeAttribute attr;

  /**
   * Creates a new node with the specified id,
   * and the specified attribute.
   * @param i	the identifier of the node in the graph.
   * @param attr  the label associated to the node.
   */
  GNode<NodeAttribute, EdgeAttribute>( int i, NodeAttribute attr ): adjacents(0), item(i), attr(attr) { };

  /*
   * Node destructor.
   * XXX: -> Destroy the list of adjacent node, without worrying about linked nodes.  What about attr ?
   */
  ~GNode(){
    GEdge<EdgeAttribute> *q,*p=adjacents;

    while ((q=p)) {
      p=p->Next();
      delete q;
    }
  };

  /**
   * Returns the list of all incident edges.
   * @return	the list of incident edges.
   */
  GEdge<EdgeAttribute> * getIncidentEdges() const { return adjacents; };

  /**
   * Connect the current node to another node identified by incidentNode. Specify the label of corresponding edge
   * @param incidentNode The identified to the node to be connected
   * @param label  The label of new edge
   * XXX: Check for directed graphs.
   * XXX: Only manage GNode ?
   */
  GEdge<EdgeAttribute> * Connect( int incidentNode, EdgeAttribute label ) { //XXX: Check for link already here
    return ( adjacents=new GEdge<EdgeAttribute>( incidentNode, adjacents, label ) );
  };


  /**
   * Connect the current node to another node identified by incidentNode. Specify the label of corresponding edge.
   * @param incidentNode The identified to the node to be connected
   * @param edge_id  The identifier of new edge
   * @param label  The label of new edge
   * XXX: Check for directed graphs.
   * XXX: Only manage GNode ?
   * XXX: Check for integrity of edge_id
   */
  GEdge<EdgeAttribute> * Connect( int incidentNode, int edge_id, EdgeAttribute attr){
    return ( adjacents=new GEdge<EdgeAttribute>( incidentNode, adjacents, edge_id, attr ) );
  };

  /**
   * Returns the degree of a node
   * @return the degree
   */
  int Degree(){
    int degree = 0;
    GEdge<EdgeAttribute>* p = adjacents;
    while(p) {degree ++;p=p->Next();}
    return degree;
  };

  /**
   * Deletes the specified node from the list of connected nodes by deleting corresponding edge.
   * @param incidentNode  the specified node to unlink.
   * @return the updated list of edges.
   */
  GEdge<EdgeAttribute>* UnConnect( int incidentNode ){
    GEdge<EdgeAttribute> *p = getIncidentEdges();
    GEdge<EdgeAttribute> *q;

    if (!p) return NULL;
    if(p->IncidentNode() == incidentNode){ //First edge is the one; we delete it (special case).
      adjacents = p->Next();
      delete p;
      return adjacents;
    }


    while (p->Next())
      {
	if(p->Next()->IncidentNode() == incidentNode){//We found corresponding edge; we delete it.
	  q = p->Next();
	  p->Next(q->Next());
	  delete q;
	}else
	  p = p->Next();
      }
    return adjacents;
  };

  /**
   * Returns the index of the referenced object.
   * @return	the index.
   */
  int Item() const { return item; }

  /**
   * Sets the new index of the referenced object.
   * @param i	the new index.
   * @return	the new index of the object.
   */
  int Item( int i ) { return item=i; }
};

/** @brief A 2D graph.
 *
 * A graph is a set of nodes connected together. A graph can be directed and undirected;
 * The graph is associated to two templated parameters : Edge and Node attributes.
 */


template< class NodeAttribute, class EdgeAttribute>
class Graph {
protected :
  std::vector<GNode<NodeAttribute, EdgeAttribute> *> tnode; //List of nodes. tnode[i] corresponds to node with id i
  int nbNodes;
  int nbEdges;
  bool _directed;

  friend class GEdge<EdgeAttribute>;

public :
  /**
   * Constructor from a gxl file. GXL file corresponds to a generic file format for graphs. See ?? for a complete description
   * @param GXL filename
   * @param readNodeLabel function to read a Node attribute node in gxl file. Must be specific to NodeAttribute type
   * @param readEdgeLabel function to read an Edge attribute node in gxl file. Must be specific to EdgeAttribute type
   * XXX: find reference for gxl
   */
  Graph(const char * filename,
	NodeAttribute (*readNodeLabel)(TiXmlElement *elem),
	EdgeAttribute (*readEdgeLabel)(TiXmlElement *elem));


  Graph(const Graph<NodeAttribute,EdgeAttribute>& g ):
    nbNodes(0),
    nbEdges(0)
    //_directed(g._directed)
  {
    _directed = 0;
    //std::cout << "graph size = " << g.Size() << std::endl;
    for (int i=0; i<g.Size(); i++){
      //std::cout << "copying vertex " << i << std::endl;
      this->Add(new GNode<NodeAttribute,EdgeAttribute> (i, g[i]->attr));
    }


    for (int i=0; i<g.Size(); i++){
      GEdge<EdgeAttribute> *p = g[i]->getIncidentEdges();
      while(p){

	if(this->_directed || i<p->IncidentNode())
        this->Link(i, p->IncidentNode() , p->attr);
        // std::cout << "linked  " << i <<  " to " << p->IncidentNode() << std::endl;}
        p = p->Next();

      }
    }


  }
  // alternate constructor that copies graph g and adds an extra disconnecte vertex new_node
  Graph(const Graph<NodeAttribute,EdgeAttribute>  & g, const GNode<NodeAttribute,EdgeAttribute> * new_node ):
    nbNodes(0),
    nbEdges(0)
    //_directed(g._directed)
  {
    _directed = 0;

    for (int i=0; i<g.Size(); i++){
      this->Add(new GNode<NodeAttribute,EdgeAttribute> (i, g[i]->attr));
      //std::cout << "copying node "<< i << " with attribute " << g[i]->attr << std::endl;
    }
    this->Add(new GNode<NodeAttribute,EdgeAttribute> (g.Size(), new_node->attr));
     //std::cout << "adding node "<< g.Size() << " with attribute " << new_node->attr << std::endl;

    for (int i=0; i<g.Size(); i++){
      GEdge<EdgeAttribute> *p = g[i]->getIncidentEdges();
      while(p){

	if(this->_directed || i<p->IncidentNode())
        this->Link(i, p->IncidentNode() , p->attr);
        // std::cout << "linked  " << i <<  " to " << p->IncidentNode() << std::endl;}
        p = p->Next();

      }
    }


  }

    // alternate constructor that copies a graph excluding the vertex of a certain index
   Graph(const Graph<NodeAttribute,EdgeAttribute>& g , const int excluded_node_ind):
    nbNodes(0),
    nbEdges(0),
    _directed(g._directed)
  {
    for (int i=0; i<g.Size(); i++){
          if (i<excluded_node_ind)
          this->Add(new GNode<NodeAttribute,EdgeAttribute> (i, g[i]->attr));
          else if (i>excluded_node_ind)
           this->Add(new GNode<NodeAttribute,EdgeAttribute> (i-1, g[i]->attr));
    }
    for (int i=0; i<excluded_node_ind; i++){
        if (i != excluded_node_ind){
            GEdge<EdgeAttribute> *p = g[i]->getIncidentEdges();
            while(p){
                if((this->_directed || i<p->IncidentNode()) && p->IncidentNode() < excluded_node_ind )
                this->Link(i, p->IncidentNode() , p->attr);
                else if((this->_directed || i<p->IncidentNode()) && p->IncidentNode() > excluded_node_ind )
                this->Link(i, p->IncidentNode()-1, p->attr);
                // std::cout << "linked  " << i <<  " to " << p->IncidentNode() << std::endl;}
                p = p->Next();
            }
        }
    }
    for (int i=excluded_node_ind+1; i<g.Size(); i++){
        if (i != excluded_node_ind){
            GEdge<EdgeAttribute> *p = g[i]->getIncidentEdges();
            while(p){
                if((this->_directed || i<p->IncidentNode()) && p->IncidentNode() < excluded_node_ind )
                this->Link(i-1, p->IncidentNode() , p->attr);
                else if((this->_directed || i<p->IncidentNode()) && p->IncidentNode() > excluded_node_ind )
                this->Link(i-1, p->IncidentNode()-1, p->attr);
                // std::cout << "linked  " << i <<  " to " << p->IncidentNode() << std::endl;}
                p = p->Next();
            }
        }
    }
  }


  /**
   * Fills current graph with contents read from a gxl file
   * @param GXL filename
   * @param readNodeLabel function to read a Node attribute node in gxl file. Must be specific to NodeAttribute type
   * @param readEdgeLabel function to read an Edge attribute node in gxl file. Must be specific to EdgeAttribute type
   */
  virtual
  void GraphLoadGXL(const char * filename,
		    NodeAttribute (*readNodeLabel)(TiXmlElement *elem),
		    EdgeAttribute (*readEdgeLabel)(TiXmlElement *elem));

  /**
   * Deletes the graph.
   */
  virtual ~Graph(){
    for (int i=0;i<nbNodes;i++) {
      if (tnode[i])
	delete tnode[i];
    }
  };

  /**
   *
   * @return true if the graph is directed.
   */
  bool isDirected() const { return _directed; }
  /**
   * Returns the number of nodes.
   * @return	the size.
   */
  int Size() const { return nbNodes; };
  /**
   * Returns the number of edges.
   * @return	the number of edges.
   * XXX: Accord it with graph theory terms order and size
   */

  int getNbEdges() const { return _directed?nbEdges:nbEdges/2; }
  /**
   * Creates a new graph with no data.
   * @param directed true for creating a directed graph.
   */
  Graph( bool directed =false): tnode(0), nbNodes(0), nbEdges(0), _directed(directed) { }


  /**
   * Returns the node with the  specified identifier.
   * @param id	the identifier.
   * @return	the corresponding node.
   */
  GNode<NodeAttribute, EdgeAttribute> *operator[]( int id ){ return(tnode[id]); }

  /**
   * Returns the node with the  specified identifier (const version).
   * @param id	the identifier.
   * @return	the corresponding node.
   */
  const GNode<NodeAttribute, EdgeAttribute> *operator[]( int pos ) const { return(tnode[pos]); }

  /**
   * Adds a new node to the graph.
   * @param node the new Node
   * @return  id of the node.
   */
  int Add( GNode<NodeAttribute, EdgeAttribute>* node ){
    tnode.push_back(node);
    nbNodes ++;
    return tnode.size();
  };

  /**
   * Deletes the specified node from the graph. Unlinks it
   * from connected nodes.
   * @param s	the node to be deleted.
   * @return	SUCCESS or FAILURE.
   */
  GNode<NodeAttribute, EdgeAttribute> * Del( int s ){
    for (int i=0;i<nbNodes;i++){
        Unlink(s,i);
    }
    GNode<NodeAttribute, EdgeAttribute> *oldNode;
    oldNode = tnode[s];
    tnode[s] = NULL;
    nbNodes --;
    return oldNode;

  };

  /**
   * Connects two nodes in the graph. Unlike <code>Connect</code> method, if the graph is undirected, creates the symmetric edge
   * @param firstNode first node to be connected
   * @param secondNode second node to be connected
   * @param label the attribute associated to created edge(s)
   * @return created <code>GEdge</code>.
   */

  GEdge<EdgeAttribute> * Link(int firstNode, int secondNode , EdgeAttribute label){
    GEdge<EdgeAttribute> * e = NULL;
    if (tnode[firstNode] != NULL && tnode[secondNode] != NULL){
      e = tnode[firstNode]->Connect(secondNode,nbEdges,label);
      nbEdges ++;
      if(!_directed){
	tnode[secondNode]->Connect(firstNode,nbEdges,label);
	nbEdges ++;
     }
    }
    return e;
  };
  //TODO : UnLink !

    GEdge<EdgeAttribute> * Unlink(int firstNode, int secondNode){
    GEdge<EdgeAttribute> * e = NULL;
    if (tnode[firstNode] != NULL && tnode[secondNode] != NULL){
      e = tnode[firstNode]->UnConnect(secondNode);
      nbEdges --;
      if(!_directed){
	tnode[secondNode]->UnConnect(firstNode);
	nbEdges --;

     }
    }
    return e;
  };


  //checks symmetry of a graph

  bool isSymmetric(){


	bool isSymmetric = true;
	int n = this->Size();
	  for (int i=0; i<n; i++){
	    GEdge<EdgeAttribute> *p1 = (*this)[i]->getIncidentEdges();
		while (p1){
			int j = p1->IncidentNode();
			GEdge<EdgeAttribute> *p2 = this->getEdge(j,i);
			if ((!p2) || (p2->attr != p1->attr)){
                               // std::cout << "Symmetry broken between vertices " <<i<< " and " << j << std::endl;
				isSymmetric = false;
                                }
			p1= p1->Next();
		 }

	}

  	return isSymmetric;
  };

  /*
   * Returns true if both nodes are linked.
   * @return TRUE if edge exists, FALSE otherwise.
   */
  bool isLinked(int firstNode, int secondNode) const {
    // Re-implemented from getEdge to keep const qualifier
   GEdge<EdgeAttribute> *p = tnode[firstNode]->getIncidentEdges();
    while(p){
      if (p->IncidentNode() == secondNode)
	return true;
      else
	p=p->Next();
    }
    return false;
  };

  /*
   * Retrieves <code>GEdge</code> between two nodes.
   * @param firstNode first incident node connecting the desired edge
   * @param secondNode sedond incident node connecting the desired edge
   * @return TRUE if edge exists, FALSE otherwise.
   */
  GEdge<EdgeAttribute> * getEdge(int firstNode, int secondNode){
    GEdge<EdgeAttribute> *p = tnode[firstNode]->getIncidentEdges();
    while(p){
      if (p->IncidentNode() == secondNode)
	return p;
      else
	p=p->Next();
    }
    return NULL;
  };

  /*
   * Retrieves the symmetric <code>GEdge</code> of a given edge emaning from a node.
   * @param nodeId Node incident to desired Edge.
   * @param p the edge
   * @return the symmetric <code>GEdge</code> of p. NULL if not found.
   */
  GEdge<EdgeAttribute> * getSymmetricEdge(int nodeId,const GEdge<EdgeAttribute> * p){
   if(!_directed){
      return getEdge(p->IncidentNode(), nodeId);
   }
   return NULL;
  };
  /* Compute a random permutation of the list of nodes  contained in the graph. Modify the underlying structure according to this permutatation
   *
   */
  void shuffleize(){
    // set some values:
    std::vector<int> perm;
    for (int i=0; i<nbNodes; ++i) perm.push_back(i);
    // using built-in random generator:

    std::random_shuffle ( perm.begin(), perm.end() ); //seed ?

    //for(int i=0;i<nbNodes;i++){
    //  std::cout << perm[i] << " ";
    //}
    // std::cout << std::endl;
    std::vector<GNode<NodeAttribute,EdgeAttribute>*> new_tnode;
    std::vector<int> inv_tnode(Size());
    //Modification of list nodes
    for(int i=0;i<nbNodes;i++){
      new_tnode.push_back(tnode[perm[i]]);
      inv_tnode[perm[i]] = i;
    }
    //Application of node permuations to graph'edges
    for(int i=0;i<nbNodes;i++){
      GEdge<EdgeAttribute> *p = tnode[i]->getIncidentEdges();
      while(p){
	p->setIncidentNode(inv_tnode[p->IncidentNode()]);
	p = p->Next();
      }
    }
    //Update the list of nodes
    for(int i=0;i<nbNodes;i++)tnode[i] = new_tnode[i];
  }
}; //End of class graph

template < class NodeAttribute, class EdgeAttribute>
void Graph<NodeAttribute,EdgeAttribute>::GraphLoadGXL(const char * filename,
						       NodeAttribute (*readNodeLabel)(TiXmlElement *elem),
						       EdgeAttribute (*readEdgeLabel)(TiXmlElement *elem)){
  std::ifstream file(filename,std::ios::in);
  std::vector<char*> v;
  //XXX: find gxl property for this point
  _directed = false;
  TiXmlDocument doc(filename );
  if(!doc.LoadFile()){
    std::cerr << "Error while loading file" << std::endl;
    std::cerr << "error #" << doc.ErrorId() << " : " << doc.ErrorDesc() << std::endl;
  }


  TiXmlHandle hdl(&doc);
  std::map<int,int> id_to_index;
  //TiXmlElement *elem = hdl.FirstChildElement().FirstChildElement().FirstChildElement().Element();
  TiXmlElement *elem = hdl.FirstChildElement().FirstChildElement("graph").FirstChildElement().Element();
  while (elem){
    if(strcmp(elem->Value(),"node") == 0){
      int id;
      if (elem->Attribute("id")[0] == 'n')
        id = std::stoi((elem->Attribute("id")) + 1);
      else
	try{
          id= std::stoi( elem->Attribute("id"));
	}catch(std::exception& e){
          std::hash<std::string> hash_fn;
	  id= hash_fn(std::string(elem->Attribute("id")));
	}

      NodeAttribute label = readNodeLabel(elem);
      id_to_index[id] = nbNodes;
      Add(new GNode<NodeAttribute, EdgeAttribute>(id,label));
    }else if (strcmp(elem->Value(),"edge") == 0){
        int from=-1;
        int to=-1;
        const char* s_from = elem->Attribute("from");
        const char* s_to = elem->Attribute("to");
        if (s_from == NULL || s_to == NULL){
          s_from = elem->Attribute("source") + 1;
	  s_to = elem->Attribute("target") + 1;
	}

	try{from = std::stoi(s_from);
	to = std::stoi(s_to);}
        catch(std::exception& e){
          std::hash<std::string> hash_fn;
	  from= hash_fn(std::string(s_from));
          to= hash_fn(std::string(s_to));

	}

	EdgeAttribute label = readEdgeLabel(elem);
	Link(id_to_index[from], id_to_index[to],label);

    }
    elem = elem->NextSiblingElement(); // iteration
  }
}

template<class NodeAttribute, class EdgeAttribute>
Graph<NodeAttribute,EdgeAttribute>::Graph(const char * filename, NodeAttribute (*readNodeLabel)(TiXmlElement *elem),EdgeAttribute (*readEdgeLabel)(TiXmlElement *elem)){
  GraphLoadGXL(filename, readNodeLabel,readEdgeLabel);
}

#endif // __PGRAPHH__
