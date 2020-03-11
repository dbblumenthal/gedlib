
#ifndef __LETTERCOSTFUNCTION_H__
#define __LETTERCOSTFUNCTION_H__

//#include <limits>
#include "GraphEditDistance.h"
#include "LetterGraph.h"



class LetterDistanceCost:public EditDistanceCost<CMUPoint,double>
{

private:

  double _tnodes;
  double _tedges;
  double _alpha;

public:

  virtual double NodeSubstitutionCost(GNode<CMUPoint,double> * n1,GNode<CMUPoint,double> * n2,Graph<CMUPoint,double> * g1,Graph<CMUPoint,double> * g2);
  virtual double NodeDeletionCost(GNode<CMUPoint,double> * n1,Graph<CMUPoint,double> * g1);
  virtual double NodeInsertionCost(GNode<CMUPoint,double> * n2,Graph<CMUPoint,double> * g2);
  virtual double EdgeSubstitutionCost(GEdge<double> * e1,GEdge<double> * e2,Graph<CMUPoint,double> * g1,Graph<CMUPoint,double> * g2);
  virtual double EdgeDeletionCost(GEdge<double> * e1,Graph<CMUPoint,double> * g1);
  virtual double EdgeInsertionCost(GEdge<double> * e2,Graph<CMUPoint,double> * g2);
   double SubstitutionCost(CMUPoint label_1, CMUPoint label_2);


  virtual LetterDistanceCost * clone() const {return new LetterDistanceCost(*this);}

  LetterDistanceCost (double tn, double te, double a) :
    _tnodes(tn),
    _tedges(te),
    _alpha(a)
  {}


};


#endif
