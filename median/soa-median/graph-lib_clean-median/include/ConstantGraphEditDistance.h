/**
 * @file ConstantGraphEditDistance.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Sun Feb  5 2017
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __CONSTANTGRAPHEDITDISTANCE_H__
#define __CONSTANTGRAPHEDITDISTANCE_H__
#include "GraphEditDistance.h"
#include "SymbolicGraph.h"


class ConstantEditDistanceCost:public EditDistanceCost<int,int>
{
private:
  double _cns;
  double _cni;
  double _cnd;
  double _ces;
  double _cei;
  double _ced;
public:

  virtual double NodeSubstitutionCost(GNode<int,int> * n1,GNode<int,int> * n2,Graph<int,int> * g1,Graph<int,int> * g2);  
  virtual double NodeDeletionCost(GNode<int,int> * n1,Graph<int,int> * g1);
  virtual double NodeInsertionCost(GNode<int,int> * n2,Graph<int,int> * g2);
  virtual double EdgeSubstitutionCost(GEdge<int> * e1,GEdge<int> * e2,Graph<int,int> * g1,Graph<int,int> * g2);
  virtual double EdgeDeletionCost(GEdge<int> * e1,Graph<int,int> * g1);
  virtual double EdgeInsertionCost(GEdge<int> * e2,Graph<int,int> * g2);
  
  virtual ConstantEditDistanceCost * clone() const { return new ConstantEditDistanceCost(*this); }
  
  ConstantEditDistanceCost(const ConstantEditDistanceCost& other) :
    _cns(other._cns), _cni(other._cni), _cnd(other._cnd),
    _ces(other._ces), _cei(other._cei), _ced(other._ced)
  {}
  
  ConstantEditDistanceCost(double cns,double cni, double cnd,
			   double ces,double cei, double ced):_cns(cns),_cni(cni), _cnd(cnd),
								 _ces(ces), _cei(cei), _ced(ced){};

  double cns(){return _cns;};
  double cni(){return _cni;};
  double cnd(){return _cnd;};

  double ces(){return _ces;};
  double cei(){return _cei;};
  double ced(){return _ced;};


};

#endif // __CONSTANTGRAPHEDITDISTANCE_H__
