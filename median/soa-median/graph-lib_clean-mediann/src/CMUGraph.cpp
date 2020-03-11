#include "CMUGraph.h"


double  CMUGraph::readCMUEdgeLabel(TiXmlElement *elem)
{
  double dist = -1.0;
  TiXmlElement* child = elem->FirstChildElement();
  
  while ( child ){
    std::string childName =  child->Attribute("name");
    TiXmlElement * child2 = child->FirstChildElement();
    if ( child2 ){
      if(childName.compare("dist")==0){
	      dist = std::stod(child2->GetText());
      }
    }
    child = child->NextSiblingElement(); // iteration
  }
  return dist;
}




CMUPoint CMUGraph::readCMUNodeLabel(TiXmlElement *elem)
{
  CMUPoint point;
  point.x = -1.0;
  point.y = -1.0;
  
  TiXmlElement* child1 = elem->FirstChildElement();
  TiXmlElement* child2 = elem->LastChild()->ToElement();
  
  if ( child1 && child2 ){
    std::string child1Name =  child1->Attribute("name");
    std::string child2Name =  child2->Attribute("name");
    TiXmlElement * child12 = child1->FirstChildElement();
    TiXmlElement * child22 = child2->FirstChildElement();
    if ( child12 && child22 ){
      if(child1Name.compare("x")==0 && child2Name.compare("y")==0){
	      point.x = std::stod(child12->GetText());
	      point.y = std::stod(child22->GetText());
      }
    }
  }
  return point;
}



CMUGraph::CMUGraph(const char* filename) : Graph< CMUPoint, double > (false)
{
  GraphLoadGXL(filename,readCMUNodeLabel,readCMUEdgeLabel);
}

