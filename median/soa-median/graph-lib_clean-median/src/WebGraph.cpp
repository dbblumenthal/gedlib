#include "WebGraph.h"


WebEAtt  WebGraph::readWebEdgeLabel(TiXmlElement *elem) // edge
{

  WebEAtt ne;
  std::string s1 = elem->Attribute("from");
  std::string s2 = elem->Attribute("to");
  ne.id =  s1+s2;
  TiXmlElement* child = elem->FirstChildElement(); // attribute
  ne.val = 0;
  while ( child ){
    std::string childName =  child->Attribute("name"); // text
    if(childName.compare("TEXT")==0){
	 ne.val =  std::stod(child->Attribute("value"));
    }

    child = child->NextSiblingElement(); // iteration
  }
  //std::cout << " EDGE   id =" << ne.id << " value =" << ne.val << std::endl;
  return ne;
}




WebNAtt WebGraph::readWebNodeLabel(TiXmlElement *elem)  // node
{

  WebNAtt na;

  TiXmlElement* child1 = elem->FirstChildElement(); // attribute
  na.id = elem->Attribute("id");
  if ( child1){
    na.freq =  std::stod(child1->Attribute("value"));

  }
  //std::cout << " NODE   id =" << na.id << " value =" << na.freq << std::endl;
  return na;
}



WebGraph::WebGraph(const char* filename) : Graph< WebNAtt, WebEAtt > (false)
{
  GraphLoadGXL(filename,readWebNodeLabel,readWebEdgeLabel);
}

