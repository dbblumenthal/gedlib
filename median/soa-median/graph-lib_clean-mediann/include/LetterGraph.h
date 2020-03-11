/**
 * @file LetterGraph.h
 * @author Évariste DALLER <<evariste.daller@unicaen.fr>> 
 *
 */
 
#ifndef __LETTERGRAPH_H__
#define __LETTERGRAPH_H__

#include "graph.h"
#include "CMUGraph.h"


class LetterGraph : public Graph< CMUPoint, double >
{

private:
  
  static double   readLetterEdgeLabel(TiXmlElement *elem);
  static CMUPoint readLetterNodeLabel(TiXmlElement *elem);


public:
  
  /**
   * Load Graph from a GXL file of the Letter Dataset 
   * \cite{LetterDataset}
   */
  LetterGraph(const char* fileName);
  
  virtual void
  GraphLoadGXL( const char * filename,
		            CMUPoint (*readNodeLabel)(TiXmlElement *elem),
        		    double (*readEdgeLabel)(TiXmlElement *elem) 
        		  );
  
};


#endif
