

#ifndef __CONSTANTMEDIANLABEL_H__
#define __CONSTANTMEDIANLABEL_H__
#include "GraphEditDistance.h"
#include "SymbolicGraph.h"
#include "MedianGraph.h"
#include "ConstantGraphEditDistance.h"


class ConstantMedianLabel:public MedianLabel<int,int,double>
{
protected :

ConstantEditDistanceCost * cf;


public:


  virtual int MedianNodeLabel(int * * mappingsFromMedian, int node1, Dataset<int,int,double> * ds);
  virtual int MedianEdgeLabel(Graph<int,int> * medianGraph , int * * mappingsFromMedian, int node1, int node2, Dataset<int,int,double> * ds);
  virtual double NodeDelta(Graph<int,int> * Gbar,int * * mappingsFromMedian, int nodeIndex, Dataset<int,int,double> * dataset);
  virtual double NodeLabelDelta(int median_size, int * * mappingsToMedian, int * * mappingsFromMedian, Graph<int,int>*& Gbar, int * NodeLabel, Dataset<int,int,double> * dataset);
  virtual int WeightedEdgeMeanLabel(int label1, int label2, double alpha);
  virtual int WeightedVertexMeanLabel(int label1, int label2, double alpha);

  ConstantMedianLabel(ConstantEditDistanceCost* costFunction):cf(costFunction){};



};

#endif // __CONSTANTMEDIANLABEL_H__
