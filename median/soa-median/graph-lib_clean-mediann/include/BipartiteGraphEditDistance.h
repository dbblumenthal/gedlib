/**
 * @file BipartiteGraphEditDistance.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>>
 * @version     0.0.1 - Mon Feb  6 2017
 *
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *
 * Bipartite Graph edit distance algorithm as described in Riesen's
 * book.  Riesen, K. (2015). Structural pattern recognition with graph
 * edit distance. Advances in Computer Vision and Pattern Recognition,
 * Cham.
 */

#ifndef __BIPARTITEGRAPHEDITDISTANCE_H__
#define __BIPARTITEGRAPHEDITDISTANCE_H__

#include "hungarian-lsape.hh"
#include "GraphEditDistance.h"
#include "utils.h"
//TODO : donner la possibilité de récupérer le mapping ?
template<class NodeAttribute, class EdgeAttribute>
class BipartiteGraphEditDistance:
  public GraphEditDistance<NodeAttribute, EdgeAttribute>
{
protected:
  double * C;

protected:
  virtual void computeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
				 Graph<NodeAttribute,EdgeAttribute> * g2);

  double SubstitutionCost(GNode<NodeAttribute,EdgeAttribute> * v1,
			  GNode<NodeAttribute,EdgeAttribute> * v2,
			  Graph<NodeAttribute,EdgeAttribute> * g1,
			  Graph<NodeAttribute,EdgeAttribute> * g2);

  double DeletionCost(GNode<NodeAttribute,EdgeAttribute> * v1,
		      Graph<NodeAttribute,EdgeAttribute> * g1);

  double InsertionCost(GNode<NodeAttribute,EdgeAttribute> * v2,
		       Graph<NodeAttribute,EdgeAttribute> * g2);

public:
  BipartiteGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction):
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
    C(NULL)
  {};

  // virtual double operator()(Graph<NodeAttribute,EdgeAttribute> * g1,
  // 			    Graph<NodeAttribute,EdgeAttribute> * g2);

  virtual void getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
				 Graph<NodeAttribute,EdgeAttribute> * g2,
				 int * G1_to_G2, int * G2_to_G2);

  double getLowerBound(Graph<NodeAttribute,EdgeAttribute> * g1,
		       Graph<NodeAttribute,EdgeAttribute> * g2,
		       int * G1_to_G2,int * G2_to_G1);

    virtual ~BipartiteGraphEditDistance(){
    if (this->C != NULL) delete [] this->C;
  }

  virtual BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute> * clone() const
  {
    return new BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute> (*this);
  }

};


//A remonter a GraphEditDistance.h ?

template<class NodeAttribute, class EdgeAttribute>
void BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
		  Graph<NodeAttribute,EdgeAttribute> * g2,
		  int * G1_to_G2,int * G2_to_G1){
  int n=g1->Size();
  int m=g2->Size();
  // Compute C
  delete [] this->C;
  computeCostMatrix(g1,g2);
  // for (int i=0;i<n+1;i++){
  //   for (int j=0;j<m+1;j++)
  //     std::cout << C[sub2ind(i,j,n+1)]<< " ";
  //   std::cout << std::endl;
  // }
  //Compute optimal assignement
  double *u = new double[n+1];
  double *v = new double[m+1];
  hungarianLSAPE(C,n+1,m+1, G1_to_G2,G2_to_G1, u,v,false);
  delete [] u;
  delete [] v;

}


template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
getLowerBound(Graph<NodeAttribute,EdgeAttribute> * g1,
	      Graph<NodeAttribute,EdgeAttribute> * g2,
	      int * G1_to_G2,int * G2_to_G1){
  int n = g1->Size();
  int m = g2->Size();
  double lower_bound = 0.0;
  for (int i =0;i<n;i++){
    int phi_i = G1_to_G2[i];
    if (phi_i < m){//substitution
      double total_substitution_cost =  this->SubstitutionCost((*g1)[i],(*g2)[phi_i],g1,g2);
      double node_substitution_cost =  this->cf->NodeSubstitutionCost((*g1)[i],(*g2)[phi_i],g1,g2);
      total_substitution_cost -= node_substitution_cost;
      total_substitution_cost = 0.5*total_substitution_cost + node_substitution_cost;
      lower_bound += total_substitution_cost;
    }else{ //deletion
      double total_deletion_cost =  this->DeletionCost((*g1)[i],g1);
      double node_deletion_cost =  this->cf->NodeDeletionCost((*g1)[i],g1);
      total_deletion_cost -= node_deletion_cost;
      total_deletion_cost = 0.5*total_deletion_cost + node_deletion_cost;
      lower_bound += total_deletion_cost;
    }
  }
  for (int j=0;j<m;j++)
    if(G2_to_G1[j] >= n){
      double total_insertion_cost =  this->InsertionCost((*g2)[j],g2);
      double node_insertion_cost =  this->cf->NodeInsertionCost((*g2)[j],g2);
      total_insertion_cost -= node_insertion_cost;
      total_insertion_cost = 0.5*total_insertion_cost + node_insertion_cost;
      lower_bound += total_insertion_cost;
    }
  return lower_bound;
}






// template<class NodeAttribute, class EdgeAttribute>
// double BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
// operator()(Graph<NodeAttribute,EdgeAttribute> * g1,
// 	   Graph<NodeAttribute,EdgeAttribute> * g2){
//   int n=g1->Size();
//   int m=g2->Size();
//   int * G1_to_G2 = new int[n];
//   int * G2_to_G1 = new int[m];
//   this->getOptimalMapping(g1,g2,G1_to_G2,G2_to_G1);
//   return this->GedFromMapping(g1,g2,G1_to_G2,n,G2_to_G1,m);
// }


template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
SubstitutionCost(GNode<NodeAttribute,EdgeAttribute> * v1,
		 GNode<NodeAttribute,EdgeAttribute> * v2,
		 Graph<NodeAttribute,EdgeAttribute> * g1,
		 Graph<NodeAttribute,EdgeAttribute> * g2){
  int n=v1->Degree();
  int m=v2->Degree();

  GEdge<EdgeAttribute> * e1 = v1->getIncidentEdges(); //edge from v1 in G1
  GEdge<EdgeAttribute> * _e2 = v2->getIncidentEdges(); //edge from v2 in G2

  GEdge<EdgeAttribute> * e2 = _e2; //We keep a copy of e2 to start again an iteration

  double * local_C = new double[(n+1) * (m+1)];
  memset(local_C,0,sizeof(double)*(n+1) * (m+1));
  for (int i=0;e1;i++){
    e2 = _e2;
    for (int j=0;e2;j++){
      local_C[sub2ind(i,j,n+1)] = this->cf->EdgeSubstitutionCost(e1,e2,g1,g2);
      // this->cf->NodeSubstitutionCost((*g1)[e1->IncidentNode()],
      // 								 (*g2)[e2->IncidentNode()],
      // 								 g1,g2) +
      e2 = e2->Next();
    }
    e1 = e1->Next();
  }

  e1 = v1->getIncidentEdges();
  for (int i=0;e1;i++){
    local_C[sub2ind(i,m,n+1)] = this->cf->EdgeDeletionCost(e1,g1); //this->cf->NodeDeletionCost((*g1)[e1->IncidentNode()],g1) +
    e1 = e1->Next();
  }
  e2 = v2->getIncidentEdges();
  for (int j=0;e2;j++){
    local_C[sub2ind(n, j, n+1)] = this->cf->EdgeInsertionCost(e2,g2); //this->cf->NodeInsertionCost((*g2)[e2->IncidentNode()],g2) +
    e2 = e2->Next();
  }
  local_C[sub2ind(n,m,n+1)] = 0;
  int *rho = new int[n];
  int *varrho = new int[m];
  double *u = new double[n+1];
  double *v = new double[m+1];
  hungarianLSAPE(local_C,n+1,m+1, rho,varrho, u,v,false);
  double cost=0.0;
  for (int i =0;i<n+1;i++)
    cost += u[i];
  for (int j =0;j<m+1;j++)
    cost += v[j];
  delete [] u;delete [] v;
  delete [] rho;delete [] varrho;
  delete [] local_C;
  return cost + this->cf->NodeSubstitutionCost(v1,v2,g1,g2);
}


template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
DeletionCost(GNode<NodeAttribute,EdgeAttribute> * v1,Graph<NodeAttribute,EdgeAttribute> * g1){
  int n=v1->Degree();

  GEdge<EdgeAttribute> * e1 = v1->getIncidentEdges();
  double cost = 0.0;
  for (int i=0;i<n;i++, e1 = e1->Next())
    cost += this->cf->EdgeDeletionCost(e1,g1); //this->cf->NodeDeletionCost((*g1)[e1->IncidentNode()],g1) +

  return cost + this->cf->NodeDeletionCost(v1,g1);
}

template<class NodeAttribute, class EdgeAttribute>
double BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
InsertionCost(GNode<NodeAttribute,EdgeAttribute> * v2,Graph<NodeAttribute,EdgeAttribute> * g2){
  int n=v2->Degree();

  GEdge<EdgeAttribute> * e2 = v2->getIncidentEdges();
  double cost = 0.0;
  for (int i=0;i<n;i++, e2 = e2->Next())
    cost  += this->cf->EdgeInsertionCost(e2,g2);// this->cf->NodeInsertionCost((*g2)[e2->IncidentNode()],g2) +
  return cost + this->cf->NodeInsertionCost(v2,g2);
}

template<class NodeAttribute, class EdgeAttribute>
void BipartiteGraphEditDistance<NodeAttribute, EdgeAttribute>::
computeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
		  Graph<NodeAttribute,EdgeAttribute> * g2){

  delete[] C;

  int n=g1->Size();
  int m=g2->Size();
  //std::cout << "graph 1 size " << n << " graph 2 size " << m <<std::endl;
  double * D = new double[(n+1) * (m+1)];
  //std::cout << "memory for cost matrix assigned" << std::endl;
  C=D;
  for (int i =0;i<n;i++)
    for(int j = 0;j<m;j++){
       //std::cout << " sub " << i << " " << j << " = ";
      C[sub2ind(i,j,n+1)] = this->SubstitutionCost((*g1)[i],(*g2)[j],g1,g2);
       //std::cout << C[sub2ind(i,j,n+1)] <<std::endl;
       }
  for (int i =0;i<n;i++)
    C[sub2ind(i,m,n+1)] = this->DeletionCost((*g1)[i],g1);

  for (int j =0;j<m;j++)
    C[sub2ind(n,j,n+1)] = this->InsertionCost((*g2)[j],g2);

  C[sub2ind(n,m,n+1)] = 0;
}
#endif // __BIPARTITEGRAPHEDITDISTANCE_H__
