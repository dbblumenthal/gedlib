
#include "WebCostFunction.h"
#include <math.h>

double WebDistanceCost::NodeSubstitutionCost(GNode<WebNAtt,WebEAtt> * n1,GNode<WebNAtt,WebEAtt> * n2,Graph<WebNAtt,WebEAtt> * g1,Graph<WebNAtt,WebEAtt> * g2)
{
  return (n1->attr.id.compare(n2->attr.id) == 0 ? 
          fabs(n1->attr.freq-n2->attr.freq)*this->_cneqs : this->_cns);
}


double WebDistanceCost::NodeDeletionCost(GNode<WebNAtt,WebEAtt> * n1,Graph<WebNAtt,WebEAtt> * g1)
{
  return this->_cnd;
}


double WebDistanceCost::NodeInsertionCost(GNode<WebNAtt,WebEAtt> * n2,Graph<WebNAtt,WebEAtt> * g2)
{
  return this->_cni;
}



double WebDistanceCost::EdgeSubstitutionCost(GEdge<WebEAtt> * e1,GEdge<WebEAtt> * e2,Graph<WebNAtt,WebEAtt> * g1,Graph<WebNAtt,WebEAtt> * g2)
{  
  return (e1->attr.id.compare(e2->attr.id) == 0 ? 
          fabs(e1->attr.val-e2->attr.val)*this->_ceeqs : this->_ces);
}


double WebDistanceCost::EdgeDeletionCost(GEdge<WebEAtt> * e1,Graph<WebNAtt,WebEAtt> * g1)
{
  return this->_ced;
}


double WebDistanceCost::EdgeInsertionCost(GEdge<WebEAtt> * e2,Graph<WebNAtt,WebEAtt> * g2)
{
  return this->_cei;
}

