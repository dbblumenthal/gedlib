

#ifndef __LETTERMEDIANLABEL_H__
#define __LETTERMEDIANLABEL_H__
#include "GraphEditDistance.h"
#include "LetterGraph.h"
#include "MedianGraph.h"
#include "LetterCostFunction.h"


class LetterMedianLabel:public MedianLabel<CMUPoint,double,int>
{
protected :

LetterDistanceCost * cf;


public:


  virtual CMUPoint MedianNodeLabel(int * * mappingsFromMedian, int node1, Dataset<CMUPoint,double,int> * ds);
  virtual double MedianEdgeLabel(Graph<CMUPoint,double> * medianGraph, int * * mappingsFromMedian, int node1, int node2, Dataset<CMUPoint,double,int> * ds);
  virtual double NodeDelta(Graph<CMUPoint,double> * Gbar,int * * mappingsFromMedian, int nodeIndex, Dataset<CMUPoint,double,int> * dataset);
  virtual double NodeLabelDelta(int median_size, int * * mappingsToMedian,int * * mappingsFromMedian, Graph<CMUPoint,double> *& Gbar, CMUPoint* NodeLabel, Dataset<CMUPoint,double,int> * dataset);
  virtual double WeightedEdgeMeanLabel(double label1, double label2, double alpha);
  virtual CMUPoint WeightedVertexMeanLabel(CMUPoint label1, CMUPoint label2, double alpha);

  LetterMedianLabel(LetterDistanceCost* costFunction):cf(costFunction){};

};

#endif // __LETTERMEDIANLABEL_H__
