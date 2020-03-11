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
class EditOperation{
protected:

double OperationCost;

public:

virtual void Apply<NodeAttribute, EdgeAttribute>(Graph<NodeAttribute,EdgeAttribute> * g1) const =0;
double getCost<NodeAttribute, EdgeAttribute>() {return OperationCost;};
virtual ~EditOperation(){};
}


template<class NodeAttribute, class EdgeAttribute>
class NodeEditOperation : {

}



template<class NodeAttribute, class EdgeAttribute>
class GraphEditDistance
{
protected:
  EditDistanceCost<NodeAttribute,EdgeAttribute> * cf;

public:

  //Mapping is an array encoding the mapping of each node
  GraphEditDistance(  EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction):cf(costFunction){};

  double GedFromMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
			Graph<NodeAttribute,EdgeAttribute> * g2,
			int * G1toG2, int n,
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


#endif // __GRAPHEDITDISTANCE_H__
