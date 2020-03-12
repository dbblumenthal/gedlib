#include "LetterCostFunction.h"
#include <cmath>

double LetterDistanceCost::SubstitutionCost(CMUPoint label_1, CMUPoint label_2)
{
  double dx = label_1.x - label_2.x;
  double dy = label_1.y - label_2.y;
  if (_power == 1) return _alpha * std::sqrt(dx*dx + dy*dy);
  if (_power == 2) return _alpha * (dx*dx + dy*dy);
  return _alpha *std::pow(std::sqrt(dx*dx + dy*dy),_power);
}

double LetterDistanceCost::NodeSubstitutionCost(GNode<CMUPoint,double> * n1,
						GNode<CMUPoint,double> * n2,
						Graph<CMUPoint,double> * g1,
						Graph<CMUPoint,double> * g2)
{
  double dx = n2->attr.x - n1->attr.x;
  double dy = n2->attr.y - n1->attr.y;
  if (_power == 1) return _alpha * std::sqrt(dx*dx + dy*dy);
  if (_power == 2) return _alpha * (dx*dx + dy*dy);
  return _alpha *std::pow(std::sqrt(dx*dx + dy*dy),_power);
}

double LetterDistanceCost::NodeDeletionCost(GNode<CMUPoint,double> * n1,Graph<CMUPoint,double> * g1)
{
  if (_alpha == 0) return _tnodes;
  return _alpha * _tnodes;
}

double LetterDistanceCost::NodeInsertionCost(GNode<CMUPoint,double> * n2,Graph<CMUPoint,double> * g2)
{
  if (_alpha == 0) return _tnodes;
  return _alpha * _tnodes;
}

double LetterDistanceCost::EdgeSubstitutionCost(GEdge<double> * e1,GEdge<double> * e2,Graph<CMUPoint,double> * g1,Graph<CMUPoint,double> * g2)
{
  double d = std::abs(e1->attr - e2->attr);
  return (1-_alpha) * d;
}

double LetterDistanceCost::EdgeDeletionCost(GEdge<double> * e1,Graph<CMUPoint,double> * g1)
{
  if (_alpha == 1) return _tedges;
  return (1-_alpha) * _tedges;
}

double LetterDistanceCost::EdgeInsertionCost(GEdge<double> * e2,Graph<CMUPoint,double> * g2)
{
  if (_alpha == 1) return _tedges;
  return (1-_alpha) * _tedges;
}
