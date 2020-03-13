#ifndef __IBD__GRAPH_H__
#define __IBD__GRAPH_H__

#include "graph.h"

class IBDGraph : public Graph< int, double >
{

private:
  
  static double readIBDEdgeLabel(TiXmlElement *elem);
  static int    readIBDNodeLabel(TiXmlElement *elem);


public:
  
  /**
   * Load Graph from a GXL file of IBD Dataset 
   * \cite{IBDDataset}
   */
  IBDGraph(const char* fileName);
  
  void GraphLoadGXL( const char * filename,
		     int (*readNodeLabel)(TiXmlElement *elem),
		     double (*readEdgeLabel)(TiXmlElement *elem) );
  
};


#endif
