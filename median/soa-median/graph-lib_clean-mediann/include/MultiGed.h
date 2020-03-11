/**
 * @file MultiGed.h
 * @author Évariste <<evariste.daller@unicaen.fr>>
 * @version     0.0.1 -  Apr 6 2017
 *
 */

#ifndef __MULTIGED_H__
#define __MULTIGED_H__

#include "GraphEditDistance.h"
#include "AllPerfectMatchingsEC.h"
#include "hungarian-lsap.hh"
#include "hungarian-lsape.hh"
#include "MappingGenerator.h"


/**
 * @brief Base class for graph edit distance methods that consider multiple linear approximations via a modified Uno algorithm
 * 
 *   This class implements some methods to compute several assignments through the LSAPE library.
 *   It may be considered as a MappingGenerator using Uno's algorithm to compute bipartite optimal mappings 
 *   from an arbitrary one, solution of the hungarian algorithm.
 * 
 * @see RandomMappings
 * @see GreedyGraphEditDistance
 */
template<class NodeAttribute, class EdgeAttribute>
  class MultiGed : public MappingGenerator<NodeAttribute, EdgeAttribute>
{

protected: /* MEMBERS */

  int _nep; //!< number of edit paths to compute GED
  double* _Clsap; //! lsap version of cost matrix
  double _ged;


protected: /* PROTECTED MEMBER FUNCTIONS */

  /**
   * @brief Compute the cost matrix corresponding to the lsap given C
   */
  virtual void computeCostMatrixLSAP(double* C, int n, int m);


public: /* CONSTRUCTORS AND ACCESSORS */

  MultiGed(int _k):
    _nep(_k), _Clsap(NULL), _ged(-1)
  {}


  virtual ~MultiGed(){
    if (_Clsap != NULL) delete[] _Clsap;
    _Clsap = NULL;
  }


  virtual void setK(int newk) { _nep = newk; } //!< Set the number of assignments to be calculated to newk
  virtual int getK() { return _nep; }  //!< Returns the number of assignments

  virtual double getGED() {return _ged; } //!< if the returned value is -1, the ged has not been computed

public: /* PUBLIC MEMBER FUNCTIONS */

  /**
   * @brief The method to compute <code>k</code> mappings is delegated to derivated classes
   *
   *  This method can include a way to compute a biartite or quadratic cost matrix, to choose 
   *  the returned assignments and so on.
   */
  
/* nico virtual std::list<int*> getMappings( Graph<NodeAttribute, EdgeAttribute>* g1,
				       Graph<NodeAttribute, EdgeAttribute>* g2,
				       int k ) = 0;
*/  
  /**
   * @brief Allocate and retruns $k$ optimal mappings between <code>g1</code> and <code>g2</code>
   * @param k  The number of mappings to compute, -1 to get all perfect matchings
   * @param C  The cost matrix
   * @return  A list of mappings given as arrays of int. For each mapping M, <code>M[i]</code> is the mapping, in g2, of node i in g1
   * @note  Each array is allocated here and have to be deleted manually
   */
  virtual std::list<int*> getKOptimalMappings(Graph<NodeAttribute,EdgeAttribute> * g1,
                                              Graph<NodeAttribute,EdgeAttribute> * g2,
                                              double* C,
                                              const int& k);

  /**
   * @brief call to getKOptimalMappings(g1, g2, C, this->_nep)
   */
  virtual std::list<int*> getKOptimalMappings(Graph<NodeAttribute,EdgeAttribute> * g1,
                                              Graph<NodeAttribute,EdgeAttribute> * g2,
                                              double* C);

 /**
  * @brief compute an optimal mapping between <code>g1</code> and <code>g2</code>
  *        from k different optimal mappings by minimizing the ged optained, according
  *        to the ged obtained from <code>graphdistance</code>
  * @note The GED is computed and set in <code>ged</code>
  */
 virtual double computeOptimalMapping( GraphEditDistance<NodeAttribute,EdgeAttribute> * graphdistance,
                                       Graph<NodeAttribute,EdgeAttribute> * g1,
                                       Graph<NodeAttribute,EdgeAttribute> * g2,
                                       double* C,
                                       int * G1_to_G2, int * G2_to_G1 );


 virtual MultiGed<NodeAttribute, EdgeAttribute> * clone() const = 0;

};


//----


template<class NodeAttribute, class EdgeAttribute>
void MultiGed<NodeAttribute, EdgeAttribute>::
computeCostMatrixLSAP(double* C, int n, int m)
{
  _Clsap = new  double[(n+m)*(n+m)];
  //memset(_Clsap, -1.0, sizeof(double)*(n+m)*(n+m)); // inf costs to all non-possible mappings
  for (int j=0; j<m+n; j++)
    for (int i=0; i<m+n; i++)
      if (i>=n && j>=m) _Clsap[sub2ind(i,j,n+m)] = 0;
      else _Clsap[sub2ind(i,j,n+m)] = -1.0;

  //XXX changer column first
  for(int i=0;i<n;i++)
    for(int j=0;j<m;j++)
      _Clsap[sub2ind(i,j,n+m)] = C[sub2ind(i,j,n+1)];

   for(int i=0;i<n;i++)
     _Clsap[sub2ind(i,m+i,n+m)] = C[sub2ind(i,m,n+1)];
   for(int j=0;j<m;j++)
     _Clsap[sub2ind(n+j,j,n+m)] = C[sub2ind(n,j,n+1)];
     
     
     //nico C représente une matrice n+1 x m+1 ? si oui, on peut donner directement la matrice _exploredMappingsSum à manger à cette procédure, et ça devrait marcher (?)
     // donc au final, lancer un genre de getKOptimalMappings( g1,g2,_psi, k)
                         
}



template<class NodeAttribute, class EdgeAttribute>
std::list<int*> MultiGed<NodeAttribute, EdgeAttribute>::
getKOptimalMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
                     Graph<NodeAttribute,EdgeAttribute> * g2,
                     double* C )
{
  return getKOptimalMappings(g1, g2, C, _nep);
}



template<class NodeAttribute, class EdgeAttribute>
std::list<int*> MultiGed<NodeAttribute, EdgeAttribute>::
getKOptimalMappings( Graph<NodeAttribute,EdgeAttribute> * g1,
                     Graph<NodeAttribute,EdgeAttribute> * g2,
                     double* C,     const int& k)
{
 
  int n=g1->Size();
  int m=g2->Size();

  // compute _Clsap
  this->computeCostMatrixLSAP(C, n, m);

  // the returned mappings
  std::list<int*> mappings;


  // Compute an optimal assignement
  double *u = new double[n+1];
  double *v = new double[m+1];
  int* G1_to_G2 = new int[n+1];
  int* G2_to_G1 = new int[m+1];

  //hungarianLSAP<double,int>(this->_Clsap, n+m, m+n, G1_to_G2, u, v, true);
  hungarianLSAPE(C, n+1, m+1, G1_to_G2, G2_to_G1, u, v, false);

  // Compute LSAP solution from LSAPE
  int* rhoperm = new int[n+m];
  bool* epsAssign = new bool[n]; // is eps[i] assigned
  for (int j=0; j<n; j++) epsAssign[j] = false;
  for (int i=0; i<n; i++){
    if (G1_to_G2[i] < m)
      rhoperm[i] = G1_to_G2[i];
    else{
      rhoperm[i] = i+m;
      epsAssign[i] = true;
    }
  }
  int firstEpsNonAssign = 0;
  for (int j=0; j<m; j++){
    if (G2_to_G1[j] == n)
      rhoperm[j+n] = j;

    else{ // find the first epsilon not assigned
      while (firstEpsNonAssign < n && epsAssign[firstEpsNonAssign]) firstEpsNonAssign++;
      rhoperm[j+n] = firstEpsNonAssign + m;
      epsAssign[firstEpsNonAssign] = true;
    }
  }

  // Compute LSAP u and v
  double *lu = new double[n+m];
  double *lv = new double[n+m];
  for (int i=0; i<n; i++) lu[i] = u[i];
  for (int i=n; i<n+m; i++) lu[i] = 0;
  for (int j=0; j<m; j++) lv[j] = v[j];
  for (int j=m; j<n+m; j++) lv[j] = 0;


  // Compute the k optimal mappings
  cDigraph<int> edg = equalityDigraph<double,int> (this->_Clsap, n+m, n+m, rhoperm, lu, lv);
  AllPerfectMatchingsEC<int> apm(edg, n, m);
  apm.enumPerfectMatchings(edg,k);
  mappings = apm.getPerfectMatchings();

  // Add the first one to the list
  mappings.push_front(rhoperm);

  delete [] this->_Clsap; this->_Clsap = NULL;
  delete [] epsAssign;
  delete [] u;
  delete [] v;
  delete [] lu;
  delete [] lv;
  delete [] G2_to_G1;
  delete [] G1_to_G2;


  return mappings;
}


template<class NodeAttribute, class EdgeAttribute>
double MultiGed<NodeAttribute, EdgeAttribute>::
computeOptimalMapping ( GraphEditDistance<NodeAttribute,EdgeAttribute> * graphdistance,
                        Graph<NodeAttribute,EdgeAttribute> * g1,
                        Graph<NodeAttribute,EdgeAttribute> * g2,
                        double* C,
                        int * G1_to_G2, int * G2_to_G1 )
{
  int n=g1->Size();
  int m=g2->Size();

  std::list<int*> mappings = getKOptimalMappings(g1, g2, C, _nep);
  //std::cerr << mappings.size() << std::endl;

  typename std::list<int*>::const_iterator it;


  // Get the min of ged;
  _ged = -1;
  double nged;

  int* local_G1_to_G2 = new int[n+1];
  int* local_G2_to_G1 = new int[m+1];

  for (it=mappings.begin(); it!=mappings.end(); it++){
    int* lsapMapping = *it;

    // computation of G1_to_G2 and G2_to_G1
    for (int j=0; j<m; j++){ // connect all to epsilon by default
      local_G2_to_G1[j] = n;
    }

    for (int i=0; i<n; i++){
      if (lsapMapping[i] >= m)
        local_G1_to_G2[i] = m; // i -> epsilon
      else{
        local_G1_to_G2[i] = lsapMapping[i];
        local_G2_to_G1[lsapMapping[i]] = i;
      }
    }

    for (int j=0; j<m; j++){
      if (lsapMapping[n+j] < m){
        local_G2_to_G1[j] = n; // epsilon -> j
      }
    }

    nged = graphdistance->GedFromMapping(g1, g2, local_G1_to_G2,n, local_G2_to_G1,m);

    // if nged is better : save the mapping and the ged
    if (_ged > nged || _ged == -1){
      _ged = nged;
      for (int i=0; i<n; i++) G1_to_G2[i] = local_G1_to_G2[i];
      for (int j=0; j<m; j++) G2_to_G1[j] = local_G2_to_G1[j];
    }
  }

  delete [] local_G1_to_G2;
  delete [] local_G2_to_G1;

  typename std::list<int*>::iterator it_del;
  for (it_del=mappings.begin(); it_del!=mappings.end(); it_del++)
    delete[] *it_del;

  return _ged;
}


#endif
