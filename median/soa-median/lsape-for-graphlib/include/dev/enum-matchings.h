// -----------------------------------------------------------
/**
 * @file enum-matchings.h
 * @brief Enumeration of matchings in a bipartite graph
 * @author Evariste Daller
 * @author Sebastien Bougleux
 * @date 14/03/2017
 * @institution Normandie Univ, UNICAEN, ENSICAEN, CNRS, GREYC, France
 */
/*
 * -----------------------------------------------------------
 * This file is part of LSAPE.
 * LSAPE is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL-C License. See README for more
 * details.
 * -----------------------------------------------------------
*/

#ifndef __ENUM_MATCHINGS_H__
#define __ENUM_MATCHINGS_H__

#include <vector>
#include <stack>
#include <list>

namespace lsape {

  /**
   * @class BipartiteSCC
   * @author evariste
   * @date 14/03/17
   * @brief Denote the presence of nodes in a bipartite-graph
   *      related strongly connected component (SCC)
   *
   * @long Describe a SCC in a bipartite graph $G=(X \cup Y, E)$.
   *      * \f$\forall 0 < i < card(X)\f$, `u[i]` is true iff \f$x_i\f$ is in the SCC
   *      * \f$\forall 0 < j < card(Y)\f$, `v[j]` is true iff \f$y_j\f$ is in the SCC
   */
  class BipartiteSCC{

  public:

    std::vector<bool> u;
    std::vector<bool> v;

  public:

    //BipartiteSCC(unsigned int size_u, unsigned int size_v);

  };

  // --------------------------------------------------------------
  /**
   * An edge
   */
  template <typename IT>
    class _Edge{

  public:
    IT x;
    IT y;

    _Edge(IT, IT);
  };

  // --------------------------------------------------------------
  /**
   * A bipartite graph as a matrix of edges
   */
  template <typename DT, typename IT>
    class Digraph{

  private:
    IT _rows;
    IT _cols;
    DT** _data;

  public:

  Digraph() :
    _rows(0),
      _cols(0),
      _data(nullptr)
      {}

  Digraph(IT rows, IT cols):
    _rows(rows),
      _cols(cols),
      _data(nullptr)
      {
	_data = new DT*[rows];
	for (IT i=0; i<rows; i++)  _data[i] = new DT[cols];
	for (IT i=0; i<rows; i++){
	  for(IT j=0; j<cols; j++){
	    _data[i][j] = 0;
	  }
	}
      }

    ~Digraph(){
      if (_data){
	for (IT i=0; i<_rows; i++)  delete[] _data[i];
	delete[] _data;
	_data = nullptr;
      }
    }

    void re_alloc(IT rows, IT cols){
      _rows = rows;
      _cols = cols;
      if (_data != nullptr){
	for (IT i=0; i<_rows; i++)  delete[] _data[i];
	delete[] _data;
	_data = nullptr;
      }

      _data = new DT*[rows];
      for (IT i=0; i<rows; i++)  _data[i] = new DT[cols];
      for (IT i=0; i<rows; i++){
	for(IT j=0; j<cols; j++){
	  _data[i][j] = 0;
	}
      }
    }

    const IT& rows() const {return _rows;}
    const IT& cols() const {return _cols;}

    DT& operator() (IT r, IT c){
      return _data[r][c];
    }

    const DT& operator() (IT r, IT c) const {
      return _data[r][c];
    }

    Digraph<DT,IT> operator+ (const Digraph<DT,IT>& other) {
      Digraph<DT,IT> res(_rows, _cols);
      for (IT i=0; i<_rows; i++){
	for (IT j=0; j<_cols; j++){
	  res(i,j) = _data[i][j] + other(i,j);
	}
      }
      return res;
    }
  };

  // --------------------------------------------------------------
  template<typename IT>
    using cDigraph = Digraph<char,IT>;

  // --------------------------------------------------------------
  /**
   * @brief returns the equality digraph as a matrix which rows represents
   *      nodes in X and columns, nodes in Y.
   * @long  For each potential edge `(x,y)` in the matrix, 3 values are possible :
   *      * 0 if there is no edge `(x,y)` in the equality digraph
   *      * 1 if there is an edge from `x` to `y`
   *      * -1 it there is an edge from `y` to `x`
   */
  template <typename DT, typename IT>
    cDigraph<IT> equalityDigraph(const DT *C, const IT &n, const IT &m,
				 const IT *rho12, const DT *u, const DT *v);
  
  // --------------------------------------------------------------
  /**
   * @class AllPerfectMatchings
   * @author evariste
   * @date 14/03/17
   * @brief
   */
  template <typename IT>
    class AllPerfectMatchings {

  protected:

    // Ressources for Tarjan (bipartite)

    int num;       //!< num of the current node
    int access;    //!< numAccess of the current node

    std::vector<int> vnum;     //!< list of number of nodes (first X then Y) @see offset
    std::vector<int> vaccess;
    std::vector<bool> instack;
    std::stack<int> tarjanStack;
    std::stack<bool> setstack; //!< in which subset is each elem of tarjanStack (true <=> X)
    std::list<BipartiteSCC> scc;

    // Storage of the perfect matchings (idx of the matched nodes in Y)
    std::list<IT*> &perfectMatchings;
    std::list<IT*> _perfectMatchings;
    bool _pMatchInClass;
    
    int offset; //!< offset to access the right subset (0 if X, n if Y)
    int offsize; //!< size of the offset

    // Storage of the deleted edges through iterations
    cDigraph<IT> gc;

    int nbmatch;  //!< number of found matchings
    int maxmatch; //!< maximum number of matching wanted

    int depth;    //!< current depth in the research tree
    int maxdepth; //!< maximum depth

  protected:
    // Iter of Tarjan
    void
      strongConnect( const cDigraph<IT>& gm, int v );


    /**
     * Iteration of Uno's algorithm
     */
    virtual void
      enumPerfectMatchings_iter( cDigraph<IT>& gm );


    /**
     * Euristic for the choice of an edge in Uno
     * @param x [inout]  the first node of the edge : \f$x \in X\f$
     * @param y [inout]  the second node of the edge : \f$y \in Y\f$
     * @return true if there is an edge, false else;
     */
    virtual bool
      chooseEdge( const cDigraph<IT>& gm, IT* x, IT* y );


    /**
     * Reverse edges in gm
     */
    void
      reverseEdges( cDigraph<IT>& gm, const std::list<_Edge<IT> >& edges ) const;


    /**
     * Save the matching corresponding to the cycle
     */
    virtual void
      saveMatching( const cDigraph<IT>& gm );


    /**
     * @brief Perform a depth first search on gm to find a cycle begining from x
     * @param gm  bipartite graph
     * @param x   first node of the cycle
     * @param cycle   where to store the cycle
     */
    bool
      depthFirstSearchCycle( const cDigraph<IT>& gm, const IT& x, std::list<_Edge<IT> >& cycle, std::vector<bool>& visited ) const;

    // Access
  public:

    /**
     * @brief returns a reference to the matchings found
     */
    const std::list<IT*> &getPerfectMatchings() const;

    /**
     * @brief Allow the user to have a quick wat to free the memory when wanted
     */
    void
      deleteMatching();

    /**
     * Initialize with the right dimensions
     */
    AllPerfectMatchings(cDigraph<IT>& gm);

    /**
     * Initialize with the right dimensions
     */
    AllPerfectMatchings(cDigraph<IT>& gm, std::list<IT*> &perfMatchings);

    /**
     * @brief free the allocated memory
     */
    ~AllPerfectMatchings();

  public:

    /**
     * @brief Find all the strongly connected components of the bipartite
     *      graph denoted by gm (\f$G=(X \cup Y, E)\f$)
     * @param gm  Adjacency matrix with $\card(X)$ rows and \f$\card(Y)\f$
     *      columns with :
     *      * \f$gm_{ij} =  1\f$ iff an edge exists from \f$x_i\f$ to \f$y_i\f$
     *      * \f$gm_{ij} = -1\f$ iff an edge exists from \f$y_j\f$ to \f$x_i\f$
     *      * \f$gm_{ij} =  0\f$ iff there is no edge between $x_i$ and \f$x_j\f$
     * @return  The SCCs in the graph
     * @see BipartiteSCC
     */
    const std::list<BipartiteSCC>&
      findSCC( const cDigraph<IT>& gm );


    /**
     * @brief Removes all edges in the graph denoted by gm that are
     *      out a scc (edges between SCC).
     * @long Applies for all potential edge in gm :
     *     \f[ gm[i,j] := gm[i,j] * \sum_{(u,v) \in scc} u_i v_j \f]
     * @param gm  [inout] The graph
     * @param scc_in List of the SCCs in the graph
     * @see findSCC
     * @return  the edges deleted
     */
    std::list<_Edge<IT> >
      rmUnnecessaryEdges( cDigraph<IT>& gm, const std::list<BipartiteSCC>& scc_in );

    /**
     * @brief Alias of rmUnnecessaryEdges(cDigraph<IT>&)
     * @long The current state of scc is used
     * @see scc
     */
    std::list<_Edge<IT> >
      rmUnnecessaryEdges( cDigraph<IT>& gm);


    /**
     * @brief enumerate all perfect matchings according to Uno's algorithm
     * @long matchings found are stored in perfectMatchings
     * @see getPerfectMatchings()
     * @param gm  Adjacency matrix with \f$\card(X)\f$ rows and \f$\card(Y)\f$
     *      columns with :
     *      * \f$gm_{ij} =  1\f$ iff an edge exists from \f$x_i\f$ to \f$y_i\f$
     *      * \f$gm_{ij} = -1\f$ iff an edge exists from \f$y_j\f$ to \f$x_i\f$
     *      * \f$gm_{ij} =  0\f$ iff there is no edge between $x_i$ and \f$x_j\f$
     *      This matrix represent the M matching in the G graph (1's)
     *
     * @param k maximum number of matching wanted before the algorithm stops (-1 for all matchings)
     * @param maxDepth  Maximum depth to reach in the research tree
     * @note The matrix will be altered (edges that are not in a SCC will be destroyed)
     */
    virtual void
      enumPerfectMatchings( cDigraph<IT>& gm , int k=-1, int maxDp=-1);

  };

} // end namespace lsape
 
#include "enum-matchings.tpp"

#endif
