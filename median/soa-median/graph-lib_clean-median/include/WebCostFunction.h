
#ifndef __WebCOSTFUNCTION_H__
#define __WebCOSTFUNCTION_H__

//#include <limits>
#include "GraphEditDistance.h"
#include "WebGraph.h"



class WebDistanceCost:public EditDistanceCost<WebNAtt,WebEAtt>
{

private:

  double _cns;
  double _cni;
  double _cnd;
  double _cneqs;
  double _ces;
  double _cei;
  double _ced;
  double _ceeqs;

public:

  virtual double NodeSubstitutionCost(GNode<WebNAtt,WebEAtt> * n1,GNode<WebNAtt,WebEAtt> * n2,Graph<WebNAtt,WebEAtt> * g1,Graph<WebNAtt,WebEAtt> * g2);
  virtual double NodeDeletionCost(GNode<WebNAtt,WebEAtt> * n1,Graph<WebNAtt,WebEAtt> * g1);
  virtual double NodeInsertionCost(GNode<WebNAtt,WebEAtt> * n2,Graph<WebNAtt,WebEAtt> * g2);
  virtual double EdgeSubstitutionCost(GEdge<WebEAtt> * e1,GEdge<WebEAtt> * e2,Graph<WebNAtt,WebEAtt> * g1,Graph<WebNAtt,WebEAtt> * g2);
  virtual double EdgeDeletionCost(GEdge<WebEAtt> * e1,Graph<WebNAtt,WebEAtt> * g1);
  virtual double EdgeInsertionCost(GEdge<WebEAtt> * e2,Graph<WebNAtt,WebEAtt> * g2);


  virtual WebDistanceCost * clone() const {return new WebDistanceCost(*this);}

  WebDistanceCost () : _cns(1000),_cni(1), _cnd(1), _cneqs(0.0001),_ces(1), _cei(1), _ced(1), _ceeqs(0.0001)
  {}

  WebDistanceCost (double cns,double cni, double cnd, double cneqs,
			   double ces,double cei, double ced, double ceeqs)
	: _cns(cns),_cni(cni), _cnd(cnd),_cneqs(cneqs),_ces(ces), _cei(cei), _ced(ced),_ceeqs(ceeqs){};

  double cns(){return _cns;};
  double cni(){return _cni;};
  double cnd(){return _cnd;};
  double cneqs(){return _cneqs;};
  double ces(){return _ces;};
  double cei(){return _cei;};
  double ced(){return _ced;};
  double ceeqs(){return _ceeqs;};



};


#endif
