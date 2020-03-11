/**
 * @file QAPlibGraph.h
 * @author Évariste DALLER <<evariste.daller@unicaen.fr>>
 *
 */

#ifndef __QAPLIBGRAPH_H__
#define __QAPLIBGRAPH_H__

#include "graph.h"

class QAPLibGraph : public Graph <int, int>
{

private:

public:

  /**
   * Load Graph from an adjacency matrix
   */
  QAPLibGraph (const int* adj_matrix, const int& n );


};


#endif
