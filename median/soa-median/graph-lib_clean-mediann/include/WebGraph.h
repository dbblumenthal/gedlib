/**
 * @file WebGraph.h
 * @author Évariste DALLER <<evariste.daller@unicaen.fr>>
 *
 */

#ifndef __WebGRAPH_H__
#define __WebGRAPH_H__

#include "graph.h"


struct WebNAtt{
  double freq;
  std::string id;

  bool operator==(const WebNAtt & a) const
{
    return (freq == a.freq && id == a.id);
}

   bool operator!=(const WebNAtt & a) const
{
    return (freq != a.freq || id != a.id);
}

} ;

struct  WebEAtt{
  double val;
  std::string id;
    bool operator==(const WebEAtt & a) const
{
    return (val == a.val && id == a.id);
}

   bool operator!=(const WebEAtt & a) const
{
    return (val != a.val || id != a.id);
}

};




class WebGraph : public Graph< WebNAtt, WebEAtt >
{

private:

  static WebEAtt   readWebEdgeLabel(TiXmlElement *elem);
  static WebNAtt readWebNodeLabel(TiXmlElement *elem);


public:

  /**
   * Load Graph from a GXL file of the Web Dataset
   * \cite{WebDataset}
   */
  WebGraph(const char* fileName);

};


#endif
