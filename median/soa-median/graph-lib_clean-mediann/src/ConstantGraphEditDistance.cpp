/*
 * @file ConstantGraphEditDistance.cpp
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Sun Feb  5 2017
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */
#include "ConstantGraphEditDistance.h"

double ConstantEditDistanceCost::NodeSubstitutionCost(GNode<int,int> * n1,
						      GNode<int,int> * n2,
						      Graph<int,int> * g1,
						      Graph<int,int> * g2){
  return (n1->attr != n2->attr) * _cns;};
  
double  ConstantEditDistanceCost::NodeDeletionCost(GNode<int,int> * n1,
						   Graph<int,int> * g1){return _cnd;};

  
double  ConstantEditDistanceCost::NodeInsertionCost(GNode<int,int> * n2,
						    Graph<int,int> * g2){return _cni;};

  
double  ConstantEditDistanceCost::EdgeSubstitutionCost(GEdge<int> * e1,
						       GEdge<int> * e2,
						       Graph<int,int> * g1,
						       Graph<int,int> * g2){
  return (e1->attr != e2->attr) * _ces;};
  
double  ConstantEditDistanceCost::EdgeDeletionCost(GEdge<int> * e1,
						   Graph<int,int> * g1){return _ced;};

  
double  ConstantEditDistanceCost::EdgeInsertionCost(GEdge<int> * e2,
						    Graph<int,int> * g2){return _cei;};
