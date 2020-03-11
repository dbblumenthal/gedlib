#ifndef __QAPLIB_COSTFUNCTION_H__
#define __QAPLIB_COSTFUNCTION_H__

#include "GraphEditDistance.h"
#include "QAPLibGraph.h"


class QAPLibCost : public EditDistanceCost <int,int>
{

private:

  double _INFINITY_;

public:

  virtual double NodeSubstitutionCost(GNode<int,int> * n1,GNode<int,int> * n2,Graph<int,int> * g1,Graph<int,int> * g2);
  virtual double NodeDeletionCost(GNode<int,int> * n1,Graph<int,int> * g1);
  virtual double NodeInsertionCost(GNode<int,int> * n2,Graph<int,int> * g2);
  virtual double EdgeSubstitutionCost(GEdge<int> * e1,GEdge<int> * e2,Graph<int,int> * g1,Graph<int,int> * g2);
  virtual double EdgeDeletionCost(GEdge<int> * e1,Graph<int,int> * g1);
  virtual double EdgeInsertionCost(GEdge<int> * e2,Graph<int,int> * g2);

  virtual QAPLibCost * clone() const { return new QAPLibCost(*this);}

  QAPLibCost () :
    _INFINITY_(100000.0)
  {}
  
};


#endif
