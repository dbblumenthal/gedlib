/**
 * @file RandomWalksGraphEditDistance.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Wed Feb  8 2017
 *
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __RANDOMWALKSGRAPHEDITDISTANCE_H__
#define __RANDOMWALKSGRAPHEDITDISTANCE_H__

#include <Eigen/Dense>
using namespace Eigen;

#include "BipartiteGraphEditDistance.h"
#include "ConstantGraphEditDistance.h"

#include "SymbolicGraph.h"

class RandomWalksGraphEditDistance:
  public virtual BipartiteGraphEditDistance<int, int>
{
protected:
  ConstantEditDistanceCost * cf;
  int _k;

  static int * labeledKron(int *m1, int nb_rows_m1,int nb_cols_m1,
			   int * m2, int nb_rows_m2, int nb_cols_m2,
			   int sizeWx[2]);
  static MatrixXi histoLab(int nbLab,  RowVectorXi IL, MatrixXi W);

  void computeCostMatrix(Graph<int,int> * g1,
			 Graph<int,int> * g2);

public:
  RandomWalksGraphEditDistance(ConstantEditDistanceCost * costFunction, int k):
    BipartiteGraphEditDistance<int,int>(costFunction),cf(costFunction),_k(k){};
};

#endif // __RANDOMWALKSGRAPHEDITDISTANCE_H__
