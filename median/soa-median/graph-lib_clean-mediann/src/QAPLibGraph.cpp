#include "QAPLibGraph.h"
#include "utils.h"


QAPLibGraph::QAPLibGraph(const int* adj_matrix, const int& n):
  Graph<int,int>(true)
{
  this->_directed = true;
  
  for (int i=0; i<n; i++){
    Add(new GNode<int, int> (i, adj_matrix[i*n+i]));
  }

  for (int j=0; j<n; j++){
    for (int i=0; i<n; i++){
      if (i != j && adj_matrix[sub2ind(i,j,n)] != 0){
        Link(i,j,adj_matrix[sub2ind(i,j,n)]);
      }
    }
  }
}
