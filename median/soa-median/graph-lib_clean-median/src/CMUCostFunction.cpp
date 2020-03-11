#include <math.h>
#include "CMUCostFunction.h"


double CMUDistanceCost::NodeSubstitutionCost(GNode<CMUPoint,double> * n1,GNode<CMUPoint,double> * n2,Graph<CMUPoint,double> * g1,Graph<CMUPoint,double> * g2)
{
  return .0;
}


double CMUDistanceCost::NodeDeletionCost(GNode<CMUPoint,double> * n1,Graph<CMUPoint,double> * g1)
{
  return _INFINITY_;
}


double CMUDistanceCost::NodeInsertionCost(GNode<CMUPoint,double> * n2,Graph<CMUPoint,double> * g2)
{
  return _INFINITY_;
}



double CMUDistanceCost::EdgeSubstitutionCost(GEdge<double> * e1,GEdge<double> * e2,Graph<CMUPoint,double> * g1,Graph<CMUPoint,double> * g2)
{
  return fabs((1-_alpha)*(e1->attr - e2->attr));
}


double CMUDistanceCost::EdgeDeletionCost(GEdge<double> * e1,Graph<CMUPoint,double> * g1)
{
  return (1-_alpha) * e1->attr;
}


double CMUDistanceCost::EdgeInsertionCost(GEdge<double> * e2,Graph<CMUPoint,double> * g2)
{
  return (1-_alpha) * e2->attr;
}


