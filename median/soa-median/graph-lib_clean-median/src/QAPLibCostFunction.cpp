#include "QAPLibCostFunction.h"



double QAPLibCost::NodeSubstitutionCost(GNode<int,int> * n1, GNode<int,int> * n2,
					Graph<int,int> * g1, Graph<int,int> * g2)
{
  return (n1->attr * n2->attr);
}

double QAPLibCost::NodeDeletionCost(GNode<int,int> * n1, Graph<int,int> * g1)
{
  return _INFINITY_;
}

double QAPLibCost::NodeInsertionCost(GNode<int,int> * n2, Graph<int,int> * g2)
{
  return _INFINITY_;
}

double QAPLibCost::EdgeSubstitutionCost(GEdge<int> * e1, GEdge<int> * e2,
					Graph<int,int> * g1, Graph<int,int> * g2)
{
  return e1->attr * e2->attr;
}

double QAPLibCost::EdgeDeletionCost(GEdge<int> * e1,Graph<int,int> * g1)
{
  return _INFINITY_;
}

double QAPLibCost::EdgeInsertionCost(GEdge<int> * e2,Graph<int,int> * g2)
{
  return _INFINITY_;
}

