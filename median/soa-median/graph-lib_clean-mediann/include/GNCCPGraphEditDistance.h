/**
 * @file GNCCPGraphEditDistance.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>> 
 * @version     0.0.1 - Wed Mar  8 2017
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __GNCCPGRAPHEDITDISTANCE_H__
#define __GNCCPGRAPHEDITDISTANCE_H__
#include <Eigen/Dense>
using namespace Eigen;

#include "GraphEditDistance.h"
#include "IPFPZetaGraphEditDistance.h"
#include "utils.h"

template<class NodeAttribute, class EdgeAttribute>
class GNCCPGraphEditDistance:
  public GraphEditDistance<NodeAttribute, EdgeAttribute>{
protected:
  
  double _d = 0.1;
  double _zeta;
  IPFPZetaGraphEditDistance<NodeAttribute,EdgeAttribute> * sub_algo;
  GraphEditDistance<NodeAttribute,EdgeAttribute> * _ed_init;
public:
  GNCCPGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction):
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),_ed_init(0){};
  GNCCPGraphEditDistance(EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction,
			 GraphEditDistance<NodeAttribute,EdgeAttribute> * ed_init):
    GraphEditDistance<NodeAttribute,EdgeAttribute>(costFunction),_ed_init(ed_init){};

  virtual void getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
				 Graph<NodeAttribute,EdgeAttribute> * g2,
				 int * G1_to_G2, int * G2_to_G1);
  
  ~GNCCPGraphEditDistance(){}

  virtual GNCCPGraphEditDistance<NodeAttribute, EdgeAttribute>* clone() const {
    return new GNCCPGraphEditDistance<NodeAttribute, EdgeAttribute>(*this);
  }

};




template<class NodeAttribute, class EdgeAttribute>
void GNCCPGraphEditDistance<NodeAttribute, EdgeAttribute>::
getOptimalMapping(Graph<NodeAttribute,EdgeAttribute> * g1,
			       Graph<NodeAttribute,EdgeAttribute> * g2,
			       int * G1_to_G2, int * G2_to_G1){
  int n = g1->Size();
  int m = g2->Size();
  
  if(this->_ed_init){
    this->_ed_init->getOptimalMapping(g1,g2,G1_to_G2, G2_to_G1);
  }else{
    for(int i =0;i<g1->Size();i++)
    G1_to_G2[i] = (i>g2->Size())?g2->Size():i;
  for(int j=0;j<g2->Size();j++)
    G2_to_G1[j] = (j>g1->Size())?g1->Size():j;
  }
  
  
  this->_zeta = 1;
  this->sub_algo = new IPFPZetaGraphEditDistance<NodeAttribute, EdgeAttribute>(this->cf, this->_zeta);

  this->sub_algo->setCurrentMatrix(G1_to_G2,G2_to_G1,n,m);
  double * Xk = this->sub_algo->getCurrentMatrix();
  Map<MatrixXd> m_Xk(Xk,n+1,m+1);
#if DEBUG
  IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
  std::cout << m_Xk.format(OctaveFmt) << std::endl;
#endif
  this->sub_algo->setMaxIter(50);
  this->sub_algo->setEpsilon(0.005);
  bool flag = true;
  while((this->_zeta > -1) && flag){
    //this->sub_algo->setMaxIter(30+ 70*(1-fabs(this->_zeta)));
    this->sub_algo->setZeta(this->_zeta);
    this->sub_algo->IPFPalgorithm(g1,g2);
#if DEBUG
    std::cout << "zeta : " << this->_zeta << std::endl;
    std::cout << m_Xk.format(OctaveFmt) << std::endl;
#endif
    this->_zeta -= this->_d;
    flag = ((m_Xk.array().round() - m_Xk.array()).abs().sum() !=0);
  }
  

#if DEBUG
  std::cout << "zeta final : " << this->_zeta << std::endl;
  std::cout << m_Xk.format(OctaveFmt) << std::endl;
#endif
  m_Xk= MatrixXd::Ones(n+1, m+1)-m_Xk;
  double *u = new double[n+1];
  double *v = new double[m+1];
  hungarianLSAPE(Xk,  n+1,  m+1, G1_to_G2,G2_to_G1, u,v,false);

#if DEBUG
  std::cout << " G1_to_G2 :" << std::endl;
  for (int i = 0; i < n ;i ++)
    std::cout << i << " -> " << G1_to_G2[i] << std::endl;
  std::cout << " G2_to_G1 :" << std::endl;
  for (int j = 0; j < m ;j ++)
    std::cout << j << " -> " << G2_to_G1[j] << std::endl;
#endif	 
  delete [] u;
  delete [] v;   
}





#endif // __GNCCPGRAPHEDITDISTANCE_H__
