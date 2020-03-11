/**
 * @file AuctionLSAP.h
 * @author Évariste <<evariste.daller@unicaen.fr>>
 * @version Apr 12 2017
 */

#ifndef __AUCTIONLSAP_H__
#define __AUCTIONLSAP_H__


#include <list>
#include "hungarian-lsap.hh"
#include "utils.hh"


/**
 * @class AuctionLSAP
 * @author evariste
 * @date 12/04/17
 * @file AuctionLSAP.h
 * @brief Base class of sequential auction algorithms
 */
template <typename DT, typename IT>
class AuctionLSAP {

protected: /* MEMBERS */

  DT _epsilon;        //!< Current minimum bid
  DT _finalEpsilon;   //!< Lower bound of epsilon to stop eps-scaling
  DT _scalingFactor;  //!< Scaling factor for epsilon scaling phase

  bool*  _unassigned;     //!< <code>_unassigned[i]</code> is true if vertex i in $G_1$ is not mapped
  IT     _nbUnassigned;   //!< Number of $i$ such that <code>_unassigned[i] == true</code>

  DT* _bids;   //! Array of current bids corresponding to bidders, index corresponds to targets @see _bidders
  IT*     _bidders; //! Array of current bidders corresponding to bids, index corresponds to targets @see _bids

  IT     _n;   //!< Size of $G_1$
  IT     _m;   //!< Size of $G_2$


protected: /* VIRTUAL MEMBER FUNCTIONS */

  /**
   * @brief The initialization phase
   *
   *    Basically allocate the arrays and give initial values to members.
   *
   * @note All parameters can be seen as initial values to begin with
   */
  virtual void initialization( const DT* C, const IT& n, const IT& m, DT* u, DT* v,
                               IT* G1_to_G2, IT* G2_to_G1, const unsigned short init_type )
  {
    if (_unassigned) delete [] _unassigned;
    if (_bids)       delete [] _bids;
    if (_bidders)    delete [] _bidders;
    _nbUnassigned = n;
    _n = n;
    _m = m;

    _unassigned = new bool[n];
    _bids = new DT[m];
    _bidders = new IT[m];

    for (IT i=0; i<n; i++) _unassigned[i] = true;
    //for (IT i=0; i<m; i++) _bids[i] = -1;
    //for (IT i=0; i<m; i++) _bidders[i] = 0;

    if(init_type == 1){
      if (n == m) basicPreprocLSAP<DT,IT> (C, n, G1_to_G2, G2_to_G1, u, v, 0);
      else if(n < m)  basicPreprocRectRowsLSAP<DT,IT> (C, n, m, G1_to_G2, G2_to_G1, u, v, 0);
      else  basicPreprocRectColsLSAP<DT,IT> (C, n, m, G1_to_G2, G2_to_G1, u, v, 0);
    }
    else{
      for (IT i=0; i<n; i++) { G1_to_G2[i] = -1; u[i] = 0; }
      for (IT i=0; i<m; i++) { G2_to_G1[i] = -1; v[i] = 0; }
    }
  }


  /**
   * @brief The method to choose the set of bidders at each iteration
   *
   *    This can be one of the three possible approaches :
   *    * The <em>Gauss-Seidel</em> approach, where only one bidder is selected at the time
   *    * The <em>Jacobi</em> approach, where all vertex in $G_1$ are selected as bidders
   *    * A combination of both, with a number of bidder between 1 and $n$
   *
   * @return A list of vertex chosen as bidders
   */
  virtual std::list<IT> chooseBidders() const = 0;


  /**
   * @brief Implements the bidding round of the auction algorithm
   *
   *    This step updates the members <code>_bids</code> and <code>_bidder</code>
   *    regarding the behavior of each bidder in <code>bidders</code>
   *
   * @param C [in]   The cost matrix
   * @param v [in]  The dual variable associated to $G_2$, $i.e.$ the price of each vertex
   * @param bidders [in]     The list of the current bidders
   */
  virtual void biddingRound( const DT* C, DT* v, const std::list<IT>& bidders )
  {
    // find the best candidate ji for each i
    IT ji;   DT pi;   //< best profit
    DT ps;   //< second best profit
    for(typename std::list<IT>::const_iterator it=bidders.begin(); it != bidders.end(); it++){
      IT i = *it;
      pi = C[sub2ind(i,0,_n)] - v[0];
      ji = 0;
      // Best profit search
      for (IT j=1; j<_m; j++){
        if (C[sub2ind(i,j,_n)] - v[j] > pi){
          pi = C[sub2ind(i,j,_n)] - v[j];
          ji = j;
        }
      }
      // Second best profit search
      if (ji == 0 && _m > 1)
        ps = C[sub2ind(i,1,_n)] - v[1];
      else
        ps = C[sub2ind(i,0,_n)] - v[0];
      for (IT j=0; j<_m; j++){
        if (j != ji && C[sub2ind(i,j,_n)] - v[j] > ps){
          ps = C[sub2ind(i,j,_n)] - v[j];
        }
      }
      bid(*it, ji, pi - ps + _epsilon);
    }
  }


  /**
   * @brief Implements the matching round of the auction algorithm
   * @param C [in]    The cost matrix
   * @param u [inout] Benefits
   * @param v [inout] Prices
   * @param G1_to_G2  [inout]  Matching
   */
  virtual void matchingRound( const DT* C, DT* u, DT* v, IT* G1_to_G2 , IT* G2_to_G1)
  {
    // for all bid > 0 :
    for (int j=0; j<_m; j++){
      if (_bids[j] > 0){
        // unmatch j if matched
        if (G2_to_G1[j] != -1){
          int i=G2_to_G1[j];
          _unassigned[i] = true;
          G1_to_G2[i] = -1;
          // swap : same nbUnassigned
        }
        else{
          _nbUnassigned--; // A new element is assigned
        }
        int i=_bidders[j];
        // match the best bidder to j
        G1_to_G2[i] = j;
        G2_to_G1[j] = i;
        _unassigned[i] = false;

        // update the dual variables
        v[j] = v[j] + _bids[j];
        u[i] = C[sub2ind(i,j,_n)] - v[j];
      }
    }
  }



protected: /* PROTECTED MEMBER FUNCTIONS */

  /**
   * @brief The given bidder tries to bid on <code>target</code> the bid <code>increment</code>
   *
   *    The procedure check if the given bid (<code>increment</code>) is greater than the
   *    current bid on target (if any). That way, only the biggest bid and the associate
   *    bidder is saved for each vertex in $G_2$.
   *
   * @param bidder
   * @param target
   * @param increment
   */
  virtual void bid(IT bidder, IT target, DT increment){
    if (_bids[target] < increment){
      _bids[target] = increment;
      _bidders[target] = bidder;
    }
  }


  /**
   * @brief The main auction algorithm, which calls to the different pure virtual steps
   * @param C   The cost matrix of size $n \times n$
   * @param n   Size of $G_1$
   * @param m   Size of $G_2$
   * @param u   Dual variables associated to $G_1$
   * @param v   Dual variables associated to $G_2$
   * @param G1_to_G2    The assignment
   */
  virtual void auctionAlgorithm( const DT* C, const int& n, const int& m, DT* u, DT* v, IT* G1_to_G2, IT* G2_to_G1 ){
    do{
      // Reset all bids
      for (int j=0; j<m; j++) _bids[j] = -1;

      std::list<IT> bidders = chooseBidders();
      biddingRound(C, v, bidders);
      matchingRound(C, u, v, G1_to_G2, G2_to_G1);
    } while (_nbUnassigned > 0);
  }


public: /* PUBLIC MEMBER FUNCTIONS */

  AuctionLSAP( const DT& firstEpsilon,
               const DT& finalEpsilon,
               const DT& scalingFactor ):
    _epsilon(firstEpsilon), _finalEpsilon(finalEpsilon), _scalingFactor(scalingFactor),
    _unassigned(NULL), _nbUnassigned(0),
    _bids(NULL), _bidders(NULL)
  {}


  /**
   * finalEpsilon will be lower than $\frac{1}{n}$
   */
  AuctionLSAP( const DT& firstEpsilon,
               const int& n ):
    _epsilon(firstEpsilon), _finalEpsilon(1.0/(n+1)), _scalingFactor(5.0),
    _unassigned(NULL), _nbUnassigned(0),
    _bids(NULL), _bidders(NULL)
  {}


  /**
   * Initialization as in [1] :
   * * firstEpsilon = max(C) / 5
   * * finalEpsilon = 1/(n+1)
   * * scaling factor = 5
   * [1] Castanon, Reverse auction algorithms for assignment problems, 1993
   */
  AuctionLSAP( const DT* C,
               const int& n,
	       const int& m):
    _epsilon(0.0), _finalEpsilon(1.0/(n+1)), _scalingFactor(5.0),
    _unassigned(NULL), _nbUnassigned(0),
    _bids(NULL), _bidders(NULL)
  {
    DT max = C[0];
    for (int j=0; j<m; j++){
      for (int i=0; i<n; i++){
        if (C[sub2ind(i,j,n)] > max) max = C[sub2ind(i,j,n)];
      }
    }
    _epsilon = max / 5;
  }

  ~AuctionLSAP(){
    if (_unassigned) delete [] _unassigned;
    if (_bids)       delete [] _bids;
    if (_bidders)    delete [] _bidders;
  }


  /**
   * @brief Apply the sequential auction algorithm several times with epsilon scaling.
   *
   *    Returns the dual variables $u$ and $v$ as well as and assignment in G1_to_G2
   *
   * @param C   The cost matrix of size $n \times n$
   * @param n   Size of $G_1$ (rows)
   * @param m   Size of $G_2$ (cols)
   * @param u   Dual variables associated to $G_1$
   * @param v   Dual variables associated to $G_2$
   * @param G1_to_G2    The assignment
   * @param G2_to_G1    The assignment from G2
   * @param init_type   0: no initialization, 1: classical (default)
   */
  void operator() ( const DT* C, const int& n, const int& m, DT* u, DT* v,
                    IT* G1_to_G2, IT* G2_to_G1, const unsigned short& init_type=1 )
  {
    initialization(C, n, m, u, v, G1_to_G2, G2_to_G1, init_type);

    while(_epsilon > _finalEpsilon){
      for (int i=0; i<n; i++){
        G1_to_G2[i] = -1;
        _unassigned[i] = true;
        _nbUnassigned = n;
      }
      for (int j=0; j<m; j++){
        G2_to_G1[j] = -1;
      }

      auctionAlgorithm(C, n, m, u, v, G1_to_G2, G2_to_G1);

      _epsilon = _epsilon / _scalingFactor;
    }
  }


  void operator() ( const DT* C, const int& n, const int& m, DT* u, DT* v,
                    IT* G1_to_G2, const unsigned short& init_type=1 )
  {
    IT* local_G2_to_G1 = new IT[m];
    (*this)(C, n, m, u, v, G1_to_G2, local_G2_to_G1, init_type);
    delete[] local_G2_to_G1;
  }


  void operator() ( const DT* C, const int& n, DT* u, DT* v,
                    IT* G1_to_G2, const unsigned short& init_type=1 )
  {
    IT* local_G2_to_G1 = new IT[n];
    (*this)(C, n, n, u, v, G1_to_G2, local_G2_to_G1, init_type );
    delete[] local_G2_to_G1;
  }


  void operator() ( const DT* C, const int& n, IT* G1_to_G2, const unsigned short& init_type=1 )
  {
    IT* local_G2_to_G1 = new IT[n];
    DT* local_u = new DT[n];
    DT* local_v = new DT[n];

    (*this)(C, n, n, local_u, local_v, G1_to_G2, local_G2_to_G1, init_type);

    delete[] local_v;
    delete[] local_u;
    delete[] local_G2_to_G1;
  }


};


#endif
