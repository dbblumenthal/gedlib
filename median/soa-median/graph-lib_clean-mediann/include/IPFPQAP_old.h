#ifndef __IPGPQUAP_H__
#define __IPGPQUAP_H__


#include <Eigen/Dense>
using namespace Eigen;
#include "hungarian-lsap.hh"
#include "GraphEditDistance.h"
#include "MappingRefinement.h"
#include "utils.h"

template<class NodeAttribute, class EdgeAttribute>
class IPFPQAP: 
  public MappingRefinement<NodeAttribute, EdgeAttribute>
{

protected:

  EditDistanceCost<NodeAttribute,EdgeAttribute> * costFunction; 
        //!< Cost function used to compute node and edge substitution costs

  int maxIter = 100;
  double epsilon = 0.001;


  //Data inherent to *one* computation
  double * C = 0;
  double * linearSubProblem = 0;
  double * XkD = 0;
  double * Xk = 0;
  double Lterm = 0;
  double oldLterm = 0;
  double * Xkp1tD = 0;
  double * bkp1 = 0;
  int _n = -1;
  int _m = -1;
  int k = -1;
  bool _directed = false;
  std::vector<double> S;
  std::vector<double> R;
  
  double* J = NULL;
  bool recenter=false;

  void (*_MappingInit)(Graph<NodeAttribute,EdgeAttribute> * g1,
                       Graph<NodeAttribute,EdgeAttribute> * g2,
                       int * G1_to_G2, int * G2_to_G2);

  virtual
  void NodeCostMatrix(Graph<NodeAttribute,EdgeAttribute> * g1,
                      Graph<NodeAttribute,EdgeAttribute> * g2);

  virtual
  double * QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
                         Graph<NodeAttribute,EdgeAttribute> * g2,
                         int * G1_to_G2, double * XkD);

  virtual
  double * QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
                         Graph<NodeAttribute,EdgeAttribute> * g2,
                         double * Matrix, double * XkD);
  
  virtual
  double * QuadraticTerm(Graph<NodeAttribute,EdgeAttribute> * g1,
                         Graph<NodeAttribute,EdgeAttribute> * g2,
                         std::vector<std::pair<std::pair<int,int>,double> > mappings, double * XkD);

  //This linearCost is efficient for sparse Xk matrices
  // n is nb rows of matrices, m is nb columns
  virtual double linearCost(double * CostMatrix, int * G1_to_G2, int n, int m);
  // n is nb rows of matrices, m is nb columns
  virtual double linearCost(double * CostMatrix, double * Xk, int n, int m);

  // Fill this->linearSubProblem with appropriatelinear problem
  virtual void LinearSubProblem();
  virtual double getCost(int * G1_to_G2, int n, int m);
  virtual double getCost(double * Matrix , int n, int m);
  virtual double getAlpha();
  virtual double getBeta();

  virtual double * mappingsToMatrix(int * G1_to_G2, int n, int m, double * Matrix);



public:
  IPFPQAP( EditDistanceCost<NodeAttribute,EdgeAttribute> * cf ):
    costFunction(cf),
    J(NULL)
  {};


  /**
   * Returns a refined mapping using IPFP algorithm
   * The initialisation is a flat continuous solution
   * i.e. $x_{ik} = 1/n$ if G1_to_G2 == G2_to_G1 == NULL
   */
  virtual void getBetterMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                                 Graph<NodeAttribute,EdgeAttribute> * g2,
                                 int * G1_to_G2, int* G2_to_G1=NULL, bool fromInit=true);

  /**
   * Get an optimal mapping from the continunous initialization given by Xk0
   */
  virtual void getBetterMappingFromInit( Graph<NodeAttribute,EdgeAttribute> * g1,
                                         Graph<NodeAttribute,EdgeAttribute> * g2,
                                         int* G1_to_G2, double * X0 = NULL);

  virtual void IPFPalgorithm(Graph<NodeAttribute,EdgeAttribute> * g1,
                     Graph<NodeAttribute,EdgeAttribute> * g2);


  virtual double mappingCost( Graph<NodeAttribute,EdgeAttribute> * g1,
                      Graph<NodeAttribute,EdgeAttribute> * g2,
                      int* G1_to_G2, int* G2_to_G1 = NULL );

  
  /**
   * @brief  At the beginning of IPFP algorithm, the initialization will be recentered in the direction of the current <code>this->J</code> matrix
   *  
   *  If J is NULL, the geometrical barycenter of doubly stochastic matrices J = (1/n) will be used
   */
  virtual void recenterInit(bool yes=true){
    recenter = yes;
  }
  
  
  /**
   * @brief  Set the centering matrix to J of size \f$n\times n\f$ and activate the centering
   *
   *  On the next IPFP algorithm call, the initialization \f$X_0\f$ will be translated to \f$\frac{1}{2}\times(X_0+J)\f$.
   *  If <code>nJ==NULL</code>, then the geometrical barycenter of doubly stochastic matrices \f$J = (1/n)\f$ will be used
   */
  virtual void recenterInit(double* nJ, int n);

  IPFPQAP * clone() const { return new IPFPQAP(*this); }

  void setMaxIter(int mi){ this->maxIter=mi; }
  void setEpsilon(double eps){ this->epsilon = eps; }

  virtual ~IPFPQAP(){
    if (this->C != NULL) delete [] this->C;
    if (this->linearSubProblem != NULL) delete [] this->linearSubProblem;
    if (this->XkD != NULL) delete [] this->XkD;
    if (this->Xk != NULL) delete [] this->Xk;
    if (this->Xkp1tD != NULL) delete [] this->Xkp1tD;
    if (this->bkp1 != NULL) delete [] this->bkp1;
    if (this->J != NULL) delete [] this->J;
  }

};


template<class NodeAttribute, class EdgeAttribute>
void IPFPQAP<NodeAttribute, EdgeAttribute>::
recenterInit(double * nJ, int n)
{
  if ( this->J != NULL) delete[] this->J;
  this->J = new double[n*n];
  
  this->recenter = true;
  if (nJ == NULL){
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
        this->J[sub2ind(i,j,n)] = 1.0/n;
  }
  else{
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
        this->J[sub2ind(i,j,n)] = nJ[sub2ind(i,j,n)];
  }
}


template<class NodeAttribute, class EdgeAttribute>
void IPFPQAP<NodeAttribute, EdgeAttribute>::
NodeCostMatrix( Graph<NodeAttribute,EdgeAttribute> * g1,
                Graph<NodeAttribute,EdgeAttribute> * g2)
{
  int n=g1->Size();
  int m=g2->Size();

  this->C = new double[n*m];
  for(int i=0;i<n;i++)
    for(int j=0;j<m;j++)
      C[sub2ind(i,j,n)] = this->costFunction->NodeSubstitutionCost((*g1)[i],(*g2)[j],g1,g2);
}


template<class NodeAttribute, class EdgeAttribute>
double * IPFPQAP<NodeAttribute,EdgeAttribute>::
QuadraticTerm( Graph<NodeAttribute,EdgeAttribute> * g1,
               Graph<NodeAttribute,EdgeAttribute> * g2,
               double * Matrix, double * XkD)
{
  int n = g1->Size();
  int m = g2->Size();

  std::vector<std::pair<std::pair<int,int>, double>> mappings;
  for(int i=0;i<n;i++)
    for(int j=0;j<m;j++){
      double value = Matrix[sub2ind(i,j,n)];
      if(value > 0.){
        std::pair<int,int> tmp = std::pair<int,int>(i,j);
        mappings.push_back(std::pair<std::pair<int,int>,double>(tmp,value));
      }
    }
  if (! XkD)
    XkD=new double[n*m];

  return this->QuadraticTerm(g1,g2,mappings, XkD);
}



template<class NodeAttribute, class EdgeAttribute>
double * IPFPQAP<NodeAttribute, EdgeAttribute>::
QuadraticTerm( Graph<NodeAttribute,EdgeAttribute> * g1,
               Graph<NodeAttribute,EdgeAttribute> * g2,
               int * G1_to_G2,double * XkD )
{
  
  int n = g1->Size();
  int m = g2->Size();
  
  //Reconstruction d'un mapping
  std::vector<std::pair<std::pair<int,int>, double>> mappings;
  for (int i =0;i<n;i++){
    if (G1_to_G2[i] >= 0){
      std::pair<int,int> tmp = std::pair<int,int>(i,G1_to_G2[i]);
      mappings.push_back(std::pair<std::pair<int,int>,double>(tmp,1.));
    }
  }
  if (! XkD)
    XkD=new double[n*m];

  return this->QuadraticTerm(g1,g2,mappings,XkD);

}



template<class NodeAttribute, class EdgeAttribute>
double * IPFPQAP<NodeAttribute, EdgeAttribute>::
QuadraticTerm( Graph<NodeAttribute,EdgeAttribute> * g1,
               Graph<NodeAttribute,EdgeAttribute> * g2,
               std::vector<std::pair<std::pair<int,int>,double> > mappings,
               double * quadraticTerm )
{
  int n = g1->Size();
  int m = g2->Size();

  if (! quadraticTerm)
    quadraticTerm=new double[n*m];

  memset(quadraticTerm,0,sizeof(double)*n*m);

  for(int j = 0; j < n; j++){
    for(int l = 0; l < m;l++){

      std::vector<std::pair<std::pair<int,int>,double> >::iterator it = mappings.begin();
      int __i=0;
      for(;it != mappings.end();it++){
        int i = it->first.first;
        int k = it->first.second;

        if(i != j  &&  k != l){
          GEdge<EdgeAttribute> * e1 = NULL;
          bool delta_e1 = false;
          e1 = g1->getEdge(i,j);
          delta_e1 = (e1 !=NULL);


          GEdge<EdgeAttribute> * e2 = NULL;
          bool delta_e2 = false;
          e2=g2->getEdge(k,l);
          delta_e2 = (e2 != NULL);

          //TODO : Optimize if sequence
          //If (i,j) and (k,l) are both same nodes,
          //no edges between them, so delta_e1 and delta_e2 are both 0, and so the cost

          double cost = 0.0;
          if (delta_e1 && delta_e2) // sub
          cost = this->costFunction->EdgeSubstitutionCost(e1,e2,g1,g2);

          quadraticTerm[sub2ind(j,l,n)] +=  cost*it->second;
        }
        __i++;
      }
      if(! this->_directed)
        quadraticTerm[sub2ind(j,l,n)] *= 0.5;
    }
  }

  return quadraticTerm;
}

/*
template<class NodeAttribute, class EdgeAttribute>
void IPFPQAP<NodeAttribute, EdgeAttribute>::
getOptimalMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                   Graph<NodeAttribute,EdgeAttribute> * g2,
                   int * G1_to_G2, int* G2_to_G1 )
{
  //Compute Mapping init
  // => shuffleize
  int n=g1->Size();
  int m=g2->Size();
  //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  unsigned seed = 22463276112;

  std::vector<int> perm(n);
  for (int i=0; i<m; i++) perm[i] = i;
  for (int i=m; i<n; i++) perm[i] = -1;
  std::shuffle(perm.begin(), perm.end(), std::default_random_engine(seed));

  for (int i=0; i<n; i++) {
    G1_to_G2[i] = perm[i];
    std::cerr << G1_to_G2[i] << " ";
  }
  std::cerr << std::endl;

  getOptimalMappingFromInit(g1, g2, G1_to_G2);

  for (int i=0; i<n; i++)
    if (G1_to_G2[i] >= 0)
      G2_to_G1[G1_to_G2[i]] = i;

  for (int j=n; j<m; j++)
    G2_to_G1[j] = -1;
}
*/



template<class NodeAttribute, class EdgeAttribute>
void IPFPQAP<NodeAttribute, EdgeAttribute>::
getBetterMapping( Graph<NodeAttribute,EdgeAttribute> * g1,
                  Graph<NodeAttribute,EdgeAttribute> * g2,
                  int * G1_to_G2, int* G2_to_G1, bool fromInit )
{
  this->_n = g1->Size();
  this->_m = g2->Size();
  
  this->Xk = new double[(this->_n) * (this->_m)];
  
  if (fromInit){
    this->Xk = this->mappingsToMatrix(G1_to_G2,this->_n,this->_m,this->Xk);
  }
  else{
    for (int j=0; j<this->_m; j++)
      for (int i=0; i<this->_n; i++)
        this->Xk[sub2ind(i,j,this->_n)] = 1.0/this->_n;
  }

  getBetterMappingFromInit(g1, g2, G1_to_G2, this->Xk);

}



template<class NodeAttribute, class EdgeAttribute>
void IPFPQAP<NodeAttribute, EdgeAttribute>::
getBetterMappingFromInit( Graph<NodeAttribute,EdgeAttribute> * g1,
                           Graph<NodeAttribute,EdgeAttribute> * g2,
                           int* G1_to_G2, double * X0)
{
  this->_n = g1->Size();
  this->_m = g2->Size();

  if (this->Xk == NULL)
    this->Xk = new double[(this->_n) * (this->_m)];

  if (this->Xk != X0){
    for (int j=0; j<this->_m; j++)
      for (int i=0; i<this->_n; i++)
	this->Xk[sub2ind(i,j,this->_n)] = X0[sub2ind(i,j,this->_n)];
  }

  Map<MatrixXd> m_Xk(Xk, this->_n, this->_m);

  this->IPFPalgorithm(g1,g2);


  m_Xk= m_Xk *-1;
  double *u = new double[this->_n];
  double *v = new double[this->_m];
  hungarianLSAP<double, int>(this->Xk,  this->_n,  this->_m, G1_to_G2, u,v);
  delete [] this->Xk; this->Xk=NULL;
  delete [] u;
  delete [] v;


}


template<class NodeAttribute, class EdgeAttribute>
void IPFPQAP<NodeAttribute, EdgeAttribute>::
IPFPalgorithm( Graph<NodeAttribute,EdgeAttribute> * g1,
               Graph<NodeAttribute,EdgeAttribute> * g2 )
{
  
  // Recenter the init mapping ?
  if (this->recenter){
    if (this->J == NULL) this->recenterInit(NULL, this->_n); // @see recenterInit(double*, int)
    for (int j=0; j<this->_n;  j++)
      for (int i=0; i<this->_n;  i++)
        this->Xk[sub2ind(i,j,this->_n)] = (this->Xk[sub2ind(i,j,this->_n)] + this->J[sub2ind(i,j,this->_n)]) / 2;
  }

  this->_directed = (g1->isDirected() || g2->isDirected());
  //We assume that Xk is filled with a matrix, binary or not

  this->_n = g1->Size();
  this->_m = g2->Size();

  S.clear();
  R.clear();

  NodeCostMatrix(g1,g2);//REdondant for GNCCP
  Map<MatrixXd> m_C(this->C, this->_n, this->_m); //REdondant for GNCCP

  this->bkp1 = new double [(this->_n) * (this->_m)];
  Map<MatrixXd> m_bkp1 (this->bkp1,  this->_n,  this->_m);
  //this->bkp1 = mappingsToMatrix(G1_to_G2,  this->_n,  this->_m,this->bkp1);

  this->XkD = this->QuadraticTerm(g1,g2,this->Xk,NULL); //REdondant for GNCCP
  Map<MatrixXd> m_XkD(XkD,  this->_n,  this->_m);

  Map<MatrixXd> m_Xk(this->Xk,  this->_n,  this->_m);

  Lterm = linearCost(this->C,this->Xk,  this->_n,  this->_m);
  S.push_back(this->getCost(Xk,  this->_n,  this->_m));
#if DEBUG
  std::cout << "S(0) = " << S.back() << std::endl;
#endif
  k=0;
  this->linearSubProblem = new double [(  this->_n) * (  this->_m)];
  Map<MatrixXd> m_linearSubProblem(linearSubProblem,  this->_n,  this->_m);


  this->Xkp1tD = new double [(  this->_n) * (  this->_m)];
  Map<MatrixXd> m_Xkp1tD (this->Xkp1tD,  this->_n,  this->_m);


  double *u = new double[this->_n];
  double *v = new double[this->_m];  int * G1_to_G2 = new int[this->_n];
  bool flag_continue = true;

  //BipartiteGraphEditDistanceMulti<int,int> ed_multi(this->cf, 30); // To know how many solutions to lsap per iteration
  while((k < this->maxIter) && flag_continue){ //TODO : fixer un epsilon, param ?
    this->XkD = QuadraticTerm(g1, g2, Xk, this->XkD);
    this->LinearSubProblem();//    should call it gradient direction

    hungarianLSAP<double,int>(linearSubProblem,  this->_n,  this->_m, G1_to_G2, u,v);
    //bkp1 is the matrix version of mapping G1_to_G2 so a binary matrix
    this->bkp1 = mappingsToMatrix(G1_to_G2,  this->_n,  this->_m,this->bkp1);
    R.push_back(linearCost(linearSubProblem,G1_to_G2,  this->_n,  this->_m));
    //std::list<int*> mappings = ed_multi.getKOptimalMappings(g1, g2, linearSubProblem, 30);
    //std::cout << mappings.size() << ", ";

#if DEBUG
    IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    std::cout << "XkD" << std::endl;
    std::cout << m_XkD.format(OctaveFmt) << std::endl;
    std::cout << "linearSubProblem" << std::endl;
    std::cout << m_linearSubProblem.format(OctaveFmt) << std::endl;
    std::cout << "bkp1" << std::endl;
    std::cout << m_bkp1.format(OctaveFmt) << std::endl;
    std::cout << "R : " << R.back() << std::endl;
#endif

    this->oldLterm = Lterm;
    this->Lterm = linearCost(this->C,G1_to_G2,  this->_n,  this->_m);
    XkD = QuadraticTerm(g1,g2,G1_to_G2,XkD);
    S.push_back(this->getCost(G1_to_G2,  this->_n,  this->_m));

#if DEBUG
    std::cout << "S : " << S.back() << std::endl;

#endif

    double alpha = getAlpha();
    double beta= getBeta();
    double t0 =0.0;
    if(beta > 0.000001)
      t0 = -alpha / (2.*beta);
    //Built a new Xk matrix (possibly not permutation)
#if DEBUG

      std::cout << "t0 : " << t0 << std::endl;
      std::cout << "Alpha : " << alpha << std::endl;
      std::cout << "Beta : " << beta << std::endl;
#endif

    /*
    if (R.back() < 0.0001 || -alpha / R.back() < this->epsilon){
      flag_continue = false;
    }
    //*/

    //*
    if (R.back() < 0.0001)
      flag_continue = (fabs(alpha) > this->epsilon);
    else
      flag_continue = (fabs(alpha / R.back()) > this->epsilon);
    //*/

    if ((beta < 0.00001) || (t0 >= 1))
      //if(flag_continue)
        memcpy(this->Xk,bkp1,sizeof(double)*(  this->_n)*(  this->_m));
        
      //Lterm = Lterm_new;
    else{
      //Line search
      MatrixXd maj_matrix(  this->_n,  this->_m);
      maj_matrix = t0*(m_bkp1 - m_Xk);
#if DEBUG
      std::cout << "line search" << std::endl;
      std::cout << "Norm de la maj : " << maj_matrix.norm() << std::endl;
#endif
      m_Xk = m_Xk + t0*(m_bkp1 - m_Xk);
      S[k+1] = S[k] - ((pow(alpha,2))/(4*beta));
      this->Lterm = linearCost(this->C, Xk,   this->_n,  this->_m);
    }
#if DEBUG
    std::cout << "Xk à l'itération " << k << std::endl;
    std::cout << m_Xk.format(OctaveFmt) << std::endl;
    std::cout << "------------------------------------------------------------"  << std::endl << std::endl;
#endif

    this->k++;
  }


  delete [] this->Xkp1tD;this->Xkp1tD=0;
  delete [] this->linearSubProblem;this->linearSubProblem=0;
  delete [] this->XkD;this->XkD = 0;
  delete [] this->C;this->C = 0;
  delete [] this->bkp1;this->bkp1 = 0;
  delete [] u;
  delete [] v;
  delete [] G1_to_G2;
#if DEBUG
  std::cout << S.back() << std::endl;
  std::cout << "Fin d'IPFP : "<< k << "iterations " << std::endl;
#endif
  //Xk contains the optimal bistochastic matrix, binary or not.

}



template<class NodeAttribute, class EdgeAttribute>
double IPFPQAP<NodeAttribute, EdgeAttribute>::
linearCost(double * CostMatrix, int * G1_to_G2, int n, int m){
  double sum = 0.0;
  for(int i=0;i<n;i++)
    if (G1_to_G2[i]>=0) sum += CostMatrix[sub2ind(i,G1_to_G2[i],n)];
  return sum;
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPQAP<NodeAttribute, EdgeAttribute>::
linearCost(double * CostMatrix, double * X, int n, int m){
  //Todo : optimiser avec dot product ?
  double sum = 0.0;
  for(int i=0;i<n;i++)
    for(int j=0;j<m;j++)
      sum += CostMatrix[sub2ind(i,j,n)] * X[sub2ind(i,j,n)];
  return sum;
}

template<class NodeAttribute, class EdgeAttribute>
void IPFPQAP<NodeAttribute, EdgeAttribute>::
LinearSubProblem(){
  Map<MatrixXd> m_linearSubProblem(this->linearSubProblem,this->_n,this->_m);

  Map<MatrixXd> m_XkD(this->XkD,this->_n,this->_m);
  Map<MatrixXd> m_C(this->C,this->_n,this->_m);

  m_linearSubProblem = 2*m_XkD + m_C;
}



template<class NodeAttribute, class EdgeAttribute>
double IPFPQAP<NodeAttribute, EdgeAttribute>::
getCost(double * Matrix, int n, int m){
  return linearCost(this->XkD,Matrix, n,m )+ this->Lterm;
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPQAP<NodeAttribute, EdgeAttribute>::
getCost(int * G1_to_G2, int n, int m){

  return linearCost(this->XkD,G1_to_G2, n, m)+ this->Lterm;
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPQAP<NodeAttribute, EdgeAttribute>::
getAlpha(){
  return this->R.back() - 2 * this->S[k] + this->oldLterm;
}

template<class NodeAttribute, class EdgeAttribute>
double IPFPQAP<NodeAttribute, EdgeAttribute>::
getBeta(){
  return S.back() + S[k] -R.back() - this->oldLterm;
}


template<class NodeAttribute, class EdgeAttribute>
double * IPFPQAP<NodeAttribute, EdgeAttribute>::
mappingsToMatrix(int * G1_to_G2, int n, int m, double * Matrix){
  memset(Matrix,0,sizeof(double)*n*m);
  for (int i =0;i<n;i++){
    //if (G1_to_G2[i] >= m) std::cout << G1_to_G2[i] << std::endl;
    if (G1_to_G2[i] >= 0)  Matrix[sub2ind(i,G1_to_G2[i],n)] = 1;
  }
  return Matrix;
}


template<class NodeAttribute, class EdgeAttribute>
double IPFPQAP<NodeAttribute, EdgeAttribute>::
mappingCost( Graph<NodeAttribute,EdgeAttribute> * g1,
             Graph<NodeAttribute,EdgeAttribute> * g2,
             int* G1_to_G2 , int* G2_to_G1)
{
  //*
  int n = g1->Size(); this->_n = n;
  int m = g2->Size(); this->_m = m;

  this->_directed = true;

  this->Xk = new double[n*m];
  this->Xk = this->mappingsToMatrix(G1_to_G2, n, m, this->Xk);
  this->XkD = QuadraticTerm(g1, g2, Xk, this->XkD);

  NodeCostMatrix(g1, g2);
  this->linearSubProblem = new double [n*m];
  this->LinearSubProblem();
  this->Lterm = linearCost(this->C, this->Xk, n, m);

  double cost = getCost(G1_to_G2, n, m);

  delete[] this->linearSubProblem; this->linearSubProblem = NULL;
  delete[] this->C;   this->C = NULL;
  delete[] this->XkD; this->XkD = NULL;
  delete[] this->Xk;  this->Xk = NULL;

  return cost;
  //*/

  /*
  double cost = 0.0;
  for (int i=0; i<g1->Size(); i++){
    for (int j=0; j<g1->Size(); j++){
      GEdge<EdgeAttribute> * e1 = g1->getEdge(i,j);
      GEdge<EdgeAttribute> * e2 = g2->getEdge(G1_to_G2[i], G1_to_G2[j]);
      if (e1 && e2)
        cost += e1->attr * e2->attr;
    }
  }
  return cost;
  //*/
}


#endif
