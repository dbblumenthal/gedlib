

#ifndef __CMUMEDIANLABEL_H__
#define __CMUMEDIANLABEL_H__
#include "GraphEditDistance.h"
#include "CMUGraph.h"
#include "MedianGraph.h"
#include "CMUCostFunction.h"


class CMUMedianLabel:public MedianLabel<CMUPoint,double,int>
{
protected :

CMUDistanceCost * cf;


public:


  virtual CMUPoint MedianNodeLabel(int * * mappingsFromMedian, int node1, Dataset<CMUPoint,double,int> * ds);
  virtual double MedianEdgeLabel(Graph<CMUPoint,double> * medianGraph, int * * mappingsFromMedian, int node1, int node2, Dataset<CMUPoint,double,int> * ds);

  CMUMedianLabel(CMUDistanceCost* costFunction):cf(costFunction){}; 

};

#endif // __CMUMEDIANLABEL_H__
