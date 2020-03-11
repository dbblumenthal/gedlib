#include "LetterGraph.h"


double  LetterGraph::readLetterEdgeLabel(TiXmlElement *elem)
{
  return 1;
}




CMUPoint LetterGraph::readLetterNodeLabel(TiXmlElement *elem)
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



LetterGraph::LetterGraph(const char* filename) : Graph< CMUPoint, double > (false)
{
  GraphLoadGXL(filename,readLetterNodeLabel,readLetterEdgeLabel);
}



void LetterGraph::GraphLoadGXL( const char * filename,
						                     CMUPoint (*readNodeLabel)(TiXmlElement *elem),
						                     double (*readEdgeLabel)(TiXmlElement *elem))
{
  std::ifstream file(filename,std::ios::in);
  std::vector<char*> v;
  //XXX: find gxl property for this point
  _directed = false;
  TiXmlDocument doc(filename );
  if(!doc.LoadFile()){
    std::cerr << "Error while loading file" << std::endl;
    std::cerr << "error #" << doc.ErrorId() << " : " << doc.ErrorDesc() << std::endl;
  }
   
  TiXmlHandle hdl(&doc);
  std::map<int,int> id_to_index;
  TiXmlElement *elem = hdl.FirstChildElement().FirstChildElement().FirstChildElement().Element();
  while (elem){
    if(strcmp(elem->Value(),"node") == 0){
      int id = std::stoi( elem->Attribute("id")+1);
      
      CMUPoint label = readNodeLabel(elem);
      id_to_index[id] = nbNodes;
      Add(new GNode<CMUPoint, double>(id,label));
    }else if (strcmp(elem->Value(),"edge") == 0){
        int from=-1;
        int to=-1;
	from = std::stoi(elem->Attribute("from")+1);
	to = std::stoi(elem->Attribute("to")+1);
	double label = readEdgeLabel(elem);
	Link(id_to_index[from], id_to_index[to],label);	      
    }
    elem = elem->NextSiblingElement(); // iteration
  }
}
