

#ifndef __WEBMEDIANLABEL_H__
#define __WEBMEDIANLABEL_H__
#include "GraphEditDistance.h"
#include "WebGraph.h"
#include "MedianGraph.h"
#include "WebCostFunction.h"


class WebMedianLabel:public MedianLabel<WebNAtt, WebEAtt, int>
{
protected :

WebDistanceCost * cf;


public:


  virtual WebNAtt MedianNodeLabel(int * * mappingsFromMedian, int node1, Dataset<WebNAtt, WebEAtt, int> * ds);
  virtual WebEAtt MedianEdgeLabel(Graph<WebNAtt, WebEAtt> * medianGraph, int * * mappingsFromMedian, int node1, int node2, Dataset<WebNAtt, WebEAtt, int> * ds);
  virtual double NodeDelta(Graph<WebNAtt,WebEAtt> * Gbar,int * * mappingsFromMedian, int nodeIndex, Dataset<WebNAtt, WebEAtt, int> * dataset);

  WebMedianLabel(WebDistanceCost* costFunction):cf(costFunction){};

};

#endif // __WEBMEDIANLABEL_H__
