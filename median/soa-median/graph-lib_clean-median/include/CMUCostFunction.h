
#ifndef __CMUCOSTFUNCTION_H__
#define __CMUCOSTFUNCTION_H__

//#include <limits>
#include "GraphEditDistance.h"
#include "CMUGraph.h"



class CMUDistanceCost:public EditDistanceCost<CMUPoint,double>
{

private:

  double _INFINITY_;
  double _alpha;
  double _tnodes;

public:

  virtual double NodeSubstitutionCost(GNode<CMUPoint,double> * n1,GNode<CMUPoint,double> * n2,Graph<CMUPoint,double> * g1,Graph<CMUPoint,double> * g2);
  virtual double NodeDeletionCost(GNode<CMUPoint,double> * n1,Graph<CMUPoint,double> * g1);
  virtual double NodeInsertionCost(GNode<CMUPoint,double> * n2,Graph<CMUPoint,double> * g2);
  virtual double EdgeSubstitutionCost(GEdge<double> * e1,GEdge<double> * e2,Graph<CMUPoint,double> * g1,Graph<CMUPoint,double> * g2);
  virtual double EdgeDeletionCost(GEdge<double> * e1,Graph<CMUPoint,double> * g1);
  virtual double EdgeInsertionCost(GEdge<double> * e2,Graph<CMUPoint,double> * g2);
  


  virtual CMUDistanceCost * clone() const {return new CMUDistanceCost(*this);}

  CMUDistanceCost () :
    //_INFINITY_(std::numeric_limits<double>::max())
    _INFINITY_(100000.0),
    _alpha(0.5),
    _tnodes(100000.0)
  {}
  


};


#endif
