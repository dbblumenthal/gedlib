/**
 * @file CMUGraph.h
 * @author Évariste DALLER <<evariste.daller@unicaen.fr>>
 *
 */

#ifndef __CMUGRAPH_H__
#define __CMUGRAPH_H__

#include "graph.h"


struct CMUPoint{
  double x;
  double y;

  bool operator==(const CMUPoint& a) const
{
    return (x == a.x && y == a.y);
}

   bool operator!=(const CMUPoint& a) const
{
    return (x != a.x || y != a.y);
}


} ;


class CMUGraph : public Graph< CMUPoint, double >
{

private:

  static double   readCMUEdgeLabel(TiXmlElement *elem);
  static CMUPoint readCMUNodeLabel(TiXmlElement *elem);


public:

  /**
   * Load Graph from a GXL file of the CMU Dataset
   * \cite{CMUDataset}
   */
  CMUGraph(const char* fileName);

};


#endif
