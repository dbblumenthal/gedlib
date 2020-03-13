#ifndef __IBD_MEDIANLABEL_H__
#define __IBD_MEDIANLABEL_H__
#include "GraphEditDistance.h"
#include "IBDGraph.h"
#include "MedianGraph.h"
#include "IBDCostFunction.h"


class IBDMedianLabel : public MedianLabel<int,double,int>
{
 protected:
  
  IBDDistanceCost *cf;
  
 public:
  
  IBDMedianLabel(IBDDistanceCost* costFunction) : cf(costFunction) {}

  int MedianNodeLabel(int * * mappingsFromMedian, int node1, Dataset<int,double,int> * ds) { return 0; }
  double MedianEdgeLabel(Graph<int,double> * medianGraph, int * * mappingsFromMedian, int node1, int node2, Dataset<int,double,int> * ds) { return 0; }
  double NodeDelta(Graph<int,double> * Gbar,int * * mappingsFromMedian, int nodeIndex, Dataset<int,double,int> * dataset) { return 0; }
  double NodeLabelDelta(int median_size, int * * mappingsToMedian,int * * mappingsFromMedian, Graph<int,double> *& Gbar, int* NodeLabel, Dataset<int,double,int> * dataset) { return 0; }
  
  double WeightedEdgeMeanLabel(double label1, double label2, double alpha)
  {
    alpha /= cf->EdgeSubstitutionCost(label1,label2);
    return alpha * label2 + (1.0 - alpha) * label1;
  }
  
  int WeightedVertexMeanLabel(int label1, int label2, double alpha)
  {
    if (alpha < cf->SubstitutionCost(label1,label2)/2.0) return label1;
    return label2;
  }

};

#endif
