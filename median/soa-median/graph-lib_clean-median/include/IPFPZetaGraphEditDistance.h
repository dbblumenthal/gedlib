/**
 * @file IPFPGraphEditDistance.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Wed Mar  8 2017
 * 
 * @todo allow a random init, or alternative initializations
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __IPFPZETAGRAPHEDITDISTANCE_H__
#define __IPFPZETAGRAPHEDITDISTANCE_H__
#include <Eigen/Dense>
using namespace Eigen;
#include "hungarian-lsape.hh"
#include "IPFPGraphEditDistance.h"
#include "RandomWalksGraphEditDistance.h"
#include "utils.h"

template<class NodeAttribute, class EdgeAttribute>
class IPFPZetaGraphEditDistance:
  public virtual IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>{
protected:
  
private:
  double _zeta;
protected:
  virtual void LinearSubProblem();
  virtual   double getCost(int * G1_to_G2,int * G2_to_G1, int n, int m);
  virtual   double getCost(double * , int n, int m);
  virtual   double getAlpha();
  virtual   double getBeta();
  
public:
  IPFPZetaGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
			    GraphEditDistance<NodeAttribute,EdgeAttribute> * ed_init,
			    double zeta):
    IPFPGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction,ed_init),_zeta(zeta){};
  IPFPZetaGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
			    double zeta):
    IPFPGraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),_zeta(zeta){};
  
  void setZeta(double zeta){
    this->_zeta = zeta;
  };
  double * getCurrentMatrix(){
    return this->Xk;
  }

  void setCurrentMatrix(double * Matrix, int n, int m ){//N and m are matrix sizes
    if(this->Xk == NULL)
      this->Xk = new double[n*m];
    memcpy(this->Xk,Matrix,n*m*sizeof(double));  
  }

  void setCurrentMatrix(int * G1_to_G2, int * G2_to_G1, int n, int m ){//N and m are graph sizes
    if(this->Xk == NULL)
      this->Xk = new double[(n+1)*(m+1)];
    this->mappingsToMatrix(G1_to_G2,G2_to_G1,n,m, this->Xk);
  }

  
};


template<class NodeAttribute, class EdgeAttribute>
void IPFPZetaGraphEditDistance<NodeAttribute, EdgeAttribute>::
LinearSubProblem(){
  Map<MatrixXd> m_linearSubProblem(this->linearSubProblem,this->_n+1,this->_m+1);

  Map<MatrixXd> m_XkD(this->XkD,this->_n+1,this->_m+1);
  Map<MatrixXd> m_Xk(this->Xk,this->_n+1,this->_m+1);
  Map<MatrixXd> m_C(this->C,this->_n+1,this->_m+1);
  
  m_linearSubProblem = ((m_XkD + m_C) * (1.-fabs(this->_zeta)) + m_Xk*this->_zeta*2) ;
  
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPZetaGraphEditDistance<NodeAttribute, EdgeAttribute>::
getCost(double * Matrix, int n, int m){
  double S_k = IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::getCost(Matrix, n, m);
  return S_k*(1-fabs(this->_zeta)) + this->_zeta*this->linearCost(this->Xk,this->Xk,n+1,m+1);
}


template<class NodeAttribute, class EdgeAttribute>
double IPFPZetaGraphEditDistance<NodeAttribute, EdgeAttribute>::
getCost(int * G1_to_G2,int * G2_to_G1, int n, int m){
  double S_k = IPFPGraphEditDistance<NodeAttribute, EdgeAttribute>::getCost(G1_to_G2,G2_to_G1, n, m);
  return S_k*(1-fabs(this->_zeta)) + this->_zeta*this->linearCost(this->bkp1,this->bkp1,n+1,m+1);
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPZetaGraphEditDistance<NodeAttribute, EdgeAttribute>::
getAlpha(){
  return this->R.back() - 2 * this->S[this->k] + (1.-fabs(this->_zeta))*this->oldLterm;
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPZetaGraphEditDistance<NodeAttribute, EdgeAttribute>::
getBeta(){
  return this->S.back() + this->S[this->k] - this->R.back() - (1.-fabs(this->_zeta))*this->oldLterm;
}





#endif // __IPFPZETAGRAPHEDITDISTANCE_H__
