/**
 * @file GraphEditDistance.h
 * @author Benoit Gaüzère <<benoit.gauzere@insa-rouen.fr>>
 * @version     0.0.1 - Tue Jan 26 2016
 *
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __GRAPHEDITDISTANCE_H__
#define __GRAPHEDITDISTANCE_H__

#include "graph.h"

// A TRANSFORMER EN CLASSE ABSTRAITE

template<class NodeAttribute, class EdgeAttribute>
class EditDistanceCost{

public :

  virtual double NodeSubstitutionCost(GNode<NodeAttribute,EdgeAttribute> * n1,
				      GNode<NodeAttribute,EdgeAttribute> * n2,
				      Graph<NodeAttribute,EdgeAttribute> * g1,
				      Graph<NodeAttribute,EdgeAttribute> * g2)=0;

  virtual double NodeDeletionCost(GNode<NodeAttribute,EdgeAttribute> * n1,
				    Graph<NodeAttribute,EdgeAttribute> * g1)=0;

  virtual double NodeInsertionCost(GNode<NodeAttribute,EdgeAttribute> * n2,
				      Graph<NodeAttribute,EdgeAttribute> * g2)=0;

  virtual double EdgeSubstitutionCost(GEdge<EdgeAttribute> * e1,
					GEdge<EdgeAttribute> * e2,
					Graph<NodeAttribute,EdgeAttribute> * g1,
					Graph<NodeAttribute,EdgeAttribute> * g2)=0;

  virtual double EdgeDeletionCost(GEdge<EdgeAttribute> * e1,
				  Graph<NodeAttribute,EdgeAttribute> * g1)=0;

  virtual double EdgeInsertionCost(GEdge<EdgeAttribute> * e2,
				   Graph<NodeAttribute,EdgeAttribute> * g2)=0;



  virtual EditDistanceCost * clone() const = 0;
};




template<class NodeAttribute, class EdgeAttribute>
class GraphEditDistance
{
protected:
  EditDistanceCost<NodeAttribute,EdgeAttribute> * cf;

public:

  enum Operation {node_add, node_del,node_sub, edge_add, edge_del,edge_sub};

  struct EditOperation{
    Operation op;
    int v1_i = 0; // index noeud 1 de G1
    int v1_j = 0; // index noeud 2 de G1 (pour les edges)
     int v2_i = 0; // index noeud 1 de G2
    int v2_j = 0; // index noeud 2 de G2 (pour les edges)
    double cost;
  };

  //Mapping is an array encoding the mapping of each node
  GraphEditDistance(  EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction):cf(costFunction){};

  double GedFromMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
			Graph<NodeAttribute,EdgeAttribute> * g2,
			int * G1toG2, int n,
			int * G2toG1, int m);


  std::vector<std::vector<struct EditOperation *>> EditPathFromMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
							Graph<NodeAttribute,EdgeAttribute> * g2,
							int * G1toG2,  int n,
							 int * G2toG1, int m);

  virtual double operator()(Graph<NodeAttribute,EdgeAttribute> * g1,
			    Graph<NodeAttribute,EdgeAttribute> * g2);

  virtual void getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
				 Graph<NodeAttribute,EdgeAttribute> * g2,
				 int * G1_to_G2, int * G2_to_G2)=0;

  EditDistanceCost<NodeAttribute,EdgeAttribute> * getCostFunction(){ return cf; }
  void setCostFunction(EditDistanceCost<NodeAttribute, EdgeAttribute> * ncf){ cf = ncf; }

  virtual ~GraphEditDistance(){};

  virtual GraphEditDistance<NodeAttribute,EdgeAttribute> * clone() const = 0;
};



//TODO mapping a éclaicir. Voir la methode de seb sur hungarian LSAPE
template<class NodeAttribute, class EdgeAttribute>
double GraphEditDistance<NodeAttribute, EdgeAttribute>::GedFromMapping(Graph<NodeAttribute, EdgeAttribute> * g1,
								       Graph<NodeAttribute, EdgeAttribute> * g2,
								       int * G1toG2, int n,
								       int * G2toG1, int m){
  int node_ins =0, node_sub=0, node_del = 0,
    edge_ins=0, edge_sub = 0, edge_del = 0;

  /*Edit distance computation*/
  double cost = 0.0;
  for (int i=0; i<n; ++i)//We process each G1 node's appariemment
    if(G1toG2[i] >= m){ //node onto esp, Deletion
      cost += cf->NodeDeletionCost((*g1)[i],g1);
      node_del ++;
    }else{
      //Substitution
      cost += cf->NodeSubstitutionCost((*g1)[i],(*g2)[G1toG2[i]],g1,g2);
      node_sub ++;
    }

  for (int i=0; i<m; ++i)//We process each G2 node's appariemment
    //We only care about G2 node's insertions
    if (G2toG1[i]>= n){
      cost += cf->NodeInsertionCost((*g2)[i],g2);
      node_ins ++;
    }

  //Edges

  // bool * g2_processed_edges = new bool[g2->getNbEdges()*2]; // Fois deux ??
  // memset(g2_processed_edges,0,sizeof(bool)*g2->getNbEdges()*2);
  double cost_edges = 0.0;
  for (long i =0; i<n; ++i){ // G1's edges traversal
    GEdge<EdgeAttribute> *p = (*g1)[i]->getIncidentEdges();
    while (p){
      /*2 possibilities : edge is deleted, or subtitued
	Subtitution condition : e = (start, end) with start and end mapped onto G2 in f_start and f_end,
	and there is an edge between f_start and f_end*/
	int start = i;
	int end = p->IncidentNode();
	int f_start = G1toG2[start];
	int f_end = G1toG2[end];

	if(( f_start < m) && (f_end < m)){
	  //Mapping of (start, end) onto G2 exists, check if an edge exists
	  GEdge<EdgeAttribute>  *mappedEdge = g2->getEdge(f_start,f_end);
	  if( mappedEdge != NULL){
	    //Edge exists !!
	    cost_edges += cf->EdgeSubstitutionCost( p,mappedEdge,g1,g2);
	    edge_sub ++;
	    // g2_processed_edges[mappedEdge->EdgeId()] = true;
	    // g2_processed_edges[g2->getEdge(f_end,f_start)->EdgeId()] = true;
	  }else{
	    //Edge does not exist in G2 -> Deletion
	    cost_edges += cf->EdgeDeletionCost(p,g1);
	    edge_del ++;
	  }
	}else{
	  //start or end has been deleted => we delete the edge.
	  cost_edges += cf->EdgeDeletionCost(p,g1);
	  edge_del ++;
	}
	p= p->Next();
    }
  }

  for (int i=0; i<m; ++i){
    GEdge<EdgeAttribute> *p = (*g2)[i]->getIncidentEdges();
    while(p){
      int start = i;
      int end = p->IncidentNode();
      int f_start = G2toG1[start];
      int f_end = G2toG1[end];
      if(( f_start < n) && (f_end < n)){
	//Mapping of (start, end) onto G2 exists, check if an edge exists
	GEdge<EdgeAttribute>  *mappedEdge = g1->getEdge(f_start,f_end);
	if( mappedEdge != NULL){
	  //Edge exists !! But it has already been processed
	  // cost_edges += cf->EdgeSubstitutionCost( p,mappedEdge,g1,g2);
	  // edge_sub ++;
	  ;
	}else{
	  //Edge does not exist in G1 -> Insertion
	  cost_edges += cf->EdgeInsertionCost(p,g2);
	  edge_ins ++;
	}
      }else{
	//start or end have been inserted => we insert the edge.
	cost_edges += cf->EdgeInsertionCost(p,g2);
	edge_ins ++;
      }
      p = p->Next();
    }
  }


  if(! g1->isDirected()){
    cost_edges = cost_edges/2;
    }
#if DEBUG

  std::cerr << "Node substitutions : " << node_sub << std::endl;
  std::cerr << "Node addition : " << node_ins << std::endl;
  std::cerr << "Node deletion : " << node_del << std::endl;
  std::cerr << "Edge substitutions : " << edge_sub << std::endl;
  std::cerr << "Edge insertion  : " << edge_ins << std::endl;
  std::cerr << "Edge deletion  : " << edge_del << std::endl;
#endif
  return cost + cost_edges;


}
template<class NodeAttribute, class EdgeAttribute>
double GraphEditDistance<NodeAttribute, EdgeAttribute>::
operator()(Graph<NodeAttribute,EdgeAttribute> * g1,
	   Graph<NodeAttribute,EdgeAttribute> * g2){
  int n=g1->Size();
  int m=g2->Size();
  int * G1_to_G2 = new int[n];
  int * G2_to_G1 = new int[m];
  this->getOptimalMapping(g1,g2,G1_to_G2,G2_to_G1);
  double ged = this->GedFromMapping(g1,g2,G1_to_G2,n,G2_to_G1,m);
  delete [] G1_to_G2;
  delete [] G2_to_G1;
  return ged;
}

template<class NodeAttribute, class EdgeAttribute>
std::vector<std::vector<struct GraphEditDistance<NodeAttribute, EdgeAttribute>::EditOperation*>>
GraphEditDistance<NodeAttribute, EdgeAttribute>::EditPathFromMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
								     Graph<NodeAttribute,EdgeAttribute> * g2,
								      int * G1toG2,  int n,
								      int * G2toG1, int m){

  int node_ins =0, node_sub=0, node_del = 0,
    edge_ins=0, edge_sub = 0, edge_del = 0;

  /*Edit distance computation*/
  std::vector<std::vector<struct GraphEditDistance<NodeAttribute, EdgeAttribute>::EditOperation*>> operations;
  for (int i = 0;i<6;i++){
    std::vector<struct GraphEditDistance<NodeAttribute, EdgeAttribute>::EditOperation*> sub_vector;
    operations.push_back(sub_vector);
  }
  //[1]: node del
  //[3]: node sub
  //[4]: node add
  //[2]: edge sub
  //[0]: edge del
  //[5]: edge add


  double cost = 0.0;
  for (unsigned int i=0; i<n; ++i)//We process each G1 node's appariemment
    if((G1toG2[i]>=m)){ //node onto esp, Deletion
      cost = cf->NodeDeletionCost((*g1)[i],g1);
      struct EditOperation * cur_node_del = new struct EditOperation();
    //  cur_node_del->op = node_del;
      cur_node_del->v1_i = i;
      cur_node_del->v2_i = 0; //Reste à 0 car non utilisé
      cur_node_del->cost = cost;
      operations[1].push_back(cur_node_del);
      node_del ++;

    }else{
      //Substitution
      cost = cf->NodeSubstitutionCost((*g1)[i],(*g2)[G1toG2[i]],g1,g2);
       if(cost>0){
          struct EditOperation * cur_node_sub = new struct EditOperation();
        //  cur_node_sub->op = node_sub;
          cur_node_sub->v1_i = i;
          cur_node_sub->v2_i = G1toG2[i];
          cur_node_sub->cost = cost;

          operations[3].push_back(cur_node_sub);
      node_sub ++;
      }
    }

  for (unsigned int i=0; i<m; ++i)//We process each G2 node's appariemment
    //We only care about G2 node's insertions
    if (G2toG1[i]>= n){
      cost = cf->NodeInsertionCost((*g2)[i],g2);
      node_ins ++;
      struct EditOperation * cur_node_ins = new struct EditOperation();
     // cur_node_ins->op = node_add;
      cur_node_ins->v1_i = 0;
      cur_node_ins->v2_i = i;
      cur_node_ins->cost = cost;
      operations[4].push_back(cur_node_ins);
    }

  //Edges

  // bool * g2_processed_edges = new bool[g2->getNbEdges()*2]; // Fois deux ??
  // memset(g2_processed_edges,0,sizeof(bool)*g2->getNbEdges()*2);
  double cost_edges = 0.0;
  for (unsigned int i =0; i<n; ++i){ // G1's edges traversal
    GEdge<EdgeAttribute> *p = (*g1)[i]->getIncidentEdges();
    while (p){
      /*2 possibilities : edge is deleted, or subtitued
	Subtitution condition : e = (start, end) with start and end mapped onto G2 in f_start and f_end,
	and there is an edge between f_start and f_end*/
	unsigned int start = i;
	unsigned int end = p->IncidentNode();
	unsigned int f_start = G1toG2[start];
	unsigned int f_end = G1toG2[end];

	if((f_start<m) && (f_end<m)){
	  //Mapping of (start, end) onto G2 exists, check if an edge exists
	  GEdge<EdgeAttribute>  *mappedEdge = g2->getEdge(f_start,f_end);
	  if( mappedEdge != NULL){
	    //Edge exists !!
	    cost_edges = cf->EdgeSubstitutionCost( p,mappedEdge,g1,g2);
	    if(cost_edges > 0){
            edge_sub ++;
            struct EditOperation * cur_edge_sub = new struct EditOperation();
           // cur_edge_sub->op = edge_sub;
            cur_edge_sub->v1_i = start;
            cur_edge_sub->v1_j = end;

            cur_edge_sub->v2_i = f_start;
            cur_edge_sub->v2_j = f_end;
            cur_edge_sub->cost = cost_edges;
            operations[2].push_back(cur_edge_sub);
	    }

	    // g2_processed_edges[mappedEdge->EdgeId()] = true;
	    // g2_processed_edges[g2->getEdge(f_end,f_start)->EdgeId()] = true;
	  }else{
	    //Edge does not exist in G2 -> Deletion
	    cost_edges = cf->EdgeDeletionCost(p,g1);
	    edge_del ++;
	    struct EditOperation * cur_edge_del = new struct EditOperation();
	    //cur_edge_del->op = edge_del;
	    cur_edge_del->v1_i = start;
	    cur_edge_del->v1_j = end;

	    cur_edge_del->v2_i = 0;
	    cur_edge_del->v2_j = 0;
	    cur_edge_del->cost = cost_edges;
	    operations[0].push_back(cur_edge_del);

	  }
	}else{
	  //start or end has been deleted => we delete the edge.
	  cost_edges = cf->EdgeDeletionCost(p,g1);
	  edge_del ++;
	  struct EditOperation *  cur_edge_del = new struct EditOperation();
	 // cur_edge_del->op = edge_del;
	  cur_edge_del->v1_i = start;
	  cur_edge_del->v1_j = end;

	  cur_edge_del->v2_i = 0;
	  cur_edge_del->v2_j = 0;
	  cur_edge_del->cost = cost_edges;
	  operations[0].push_back(cur_edge_del);

	}
	p= p->Next();
    }
  }

  for (unsigned int i=0; i<m; ++i){
    GEdge<EdgeAttribute> *p = (*g2)[i]->getIncidentEdges();
    while(p){
      unsigned int start = i;
      unsigned int end = p->IncidentNode();
      unsigned int f_start = G2toG1[start];
      unsigned int f_end = G2toG1[end];
      if((f_start<n) && (f_end<n)){
	//Mapping of (start, end) onto G2 exists, check if an edge exists
	GEdge<EdgeAttribute>  *mappedEdge = g1->getEdge(f_start,f_end);
	if( mappedEdge != NULL){
	  //Edge exists !! But it has already been processed
	  // cost_edges += cf->EdgeSubstitutionCost( p,mappedEdge,g1,g2);
	  // edge_sub ++;
	  ;
	}else{
	  //Edge does not exist in G1 -> Insertion
	  cost_edges = cf->EdgeInsertionCost(p,g2);
	  edge_ins ++;

	  struct EditOperation * cur_edge_add = new struct EditOperation();
	//  cur_edge_add->op = edge_add;
	  cur_edge_add->v1_i = f_start;
	  cur_edge_add->v1_j = f_end;


	  cur_edge_add->v2_i = start;
	  cur_edge_add->v2_j = end;
	  cur_edge_add->cost = cost_edges;
	  operations[5].push_back(cur_edge_add);

	}
      }else{
	//start or end have been inserted => we insert the edge.
	cost_edges = cf->EdgeInsertionCost(p,g2);
	edge_ins ++;

	struct EditOperation * cur_edge_add = new struct EditOperation();
//	cur_edge_add->op = edge_add;

	cur_edge_add->v1_i = 0;
	cur_edge_add->v1_j = 0;
	if(f_start<n) cur_edge_add->v1_i = f_start;
	if(f_end<n) cur_edge_add->v1_j = f_end;

	cur_edge_add->v2_i = start;
	cur_edge_add->v2_j = end;
	cur_edge_add->cost = cost_edges;
	operations[5].push_back(cur_edge_add);
      }
      p = p->Next();
    }
  }
  if(! g1->isDirected()){
  //delete the double edges operations
  for (int i=0;i<operations[0].size();i++){
    for (int j =i+1;j<operations[0].size();j++){
        if ((operations[0][i]->v1_i==operations[0][j]->v1_j) && (operations[0][i]->v1_j==operations[0][j]->v1_i)){
            operations[0].erase(operations[0].begin()+j);}
        }
    }

  for (int i=0;i<operations[2].size();i++){
    for (int j =i+1;j<operations[2].size();j++){
        if ((operations[2][i]->v1_i==operations[2][j]->v1_j) && (operations[2][i]->v1_j==operations[2][j]->v1_i)){
            operations[2].erase(operations[2].begin()+j);
        }
    }
  }
  for (int i=0;i<operations[5].size();i++){
    for (int j =i+1;j<operations[5].size();j++){
        if ((operations[5][i]->v1_i==operations[5][j]->v1_j) && (operations[5][i]->v1_j==operations[5][j]->v1_i)){
            operations[5].erase(operations[5].begin()+j);
        }
    }

  }



  }


  // if(! g1->isDirected())
  //   cost_edges = cost_edges/2;

#if DEBUG

  std::cerr << "Node substitutions : " << node_sub << std::endl;
  std::cerr << "Node addition : " << node_ins << std::endl;
  std::cerr << "Node deletion : " << node_del << std::endl;
  std::cerr << "Edge substitutions : " << edge_sub << std::endl;
  std::cerr << "Edge insertion  : " << edge_ins << std::endl;
  std::cerr << "Edge deletion  : " << edge_del << std::endl;
#endif
  return operations;

}

#endif // __GRAPHEDITDISTANCE_H__
