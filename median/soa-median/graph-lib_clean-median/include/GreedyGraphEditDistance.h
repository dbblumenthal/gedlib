/**
 * @file GreedyBipartiteGraphEditDistance.h
 * @author Evariste <<evariste.daller@unicaen.fr>>
 *
 *
 */


#ifndef __GREEDYBIPARTITEGED_H__
#define __GREEDYBIPARTITEGED_H__


#include "BipartiteGraphEditDistanceMulti.h"


template <class NodeAttribute, class EdgeAttribute>
class GreedyGraphEditDistance :
  public BipartiteGraphEditDistanceMulti<NodeAttribute, EdgeAttribute>
{

protected:

  //int* _costs;

public:

  GreedyGraphEditDistance( EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
			   int nep ):
    BipartiteGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),
    BipartiteGraphEditDistanceMulti<NodeAttribute,EdgeAttribute>(costFunction, nep)
  {};


public:

  /**
   * @brief Allocate and retruns $k$ optimal mappings between <code>g1</code> and <code>g2</code>
   * @param k  The number of mappings to compute, -1 to get all perfect matchings
   * @param C  The cost matrix
   * @return  A list of mappings given as arrays of int. For each mapping M, <code>M[i]</code> is the mapping, in g2, of node i in g1
   * @note  Each array is allocated here and have to be deleted manually
   */
  virtual std::list<int*> getMappings(Graph<NodeAttribute,EdgeAttribute> * g1,
				      Graph<NodeAttribute,EdgeAttribute> * g2,
				      int k );

};






template<class NodeAttribute, class EdgeAttribute>
std::list<int*> GreedyGraphEditDistance<NodeAttribute, EdgeAttribute>::
getMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
	     Graph<NodeAttribute,EdgeAttribute> * g2,
	     int k)
{
  int n=g1->Size();
  int m=g2->Size();
  
  this->computeCostMatrix(g1, g2);
  
  // Compute integer cost matrix
  int* Ci = new int[(n+1)*(m+1)];

  for (int j=0; j<=m; j++){
    for (int i=0; i<=n; i++){
      Ci[sub2ind(i,j,n+1)] = (int)(this->C[sub2ind(i,j,n+1)]);
    }
  }
  // the returned mappings
  std::list<int*> mappings;

  cDigraph<int> dg = greedySortDigraph<int, int>(Ci, n+1, m+1);
  AllPerfectMatchingsEC<int> apm(dg,n,m);
  apm.enumPerfectMatchings(dg, this->_nep);
  mappings = apm.getPerfectMatchings();
  
  delete [] Ci;

  return mappings;
}


#endif
