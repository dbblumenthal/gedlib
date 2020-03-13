#include "IBDGraph.h"

// -------------------------------------------------------------------------
double IBDGraph::readIBDEdgeLabel(TiXmlElement *elem)
{
  double lab = 0;
  TiXmlElement* child1 = elem->FirstChildElement();
  if ( child1 ) {
    std::string child1Name =  child1->Attribute("name");
    if(child1Name.compare("nlogratio")==0) {
      TiXmlElement *child12 = child1->FirstChildElement();
      if ( child12 ) lab = std::stod(child12->GetText());
    }
  }
  return lab;
}
// -------------------------------------------------------------------------
int IBDGraph::readIBDNodeLabel(TiXmlElement *elem)
{
  int lab = 0;
  TiXmlElement* child1 = elem->FirstChildElement();
  if ( child1 ) {
    std::string child1Name =  child1->Attribute("name");
    if(child1Name.compare("OTU")==0) {
      TiXmlElement *child12 = child1->FirstChildElement();
      if ( child12 ) lab = std::stoi(child12->GetText());
    }
  }
  return lab;
}
// -------------------------------------------------------------------------
IBDGraph::IBDGraph(const char* filename) : Graph< int, double > (false)
{
  GraphLoadGXL(filename,readIBDNodeLabel,readIBDEdgeLabel);
}
// -------------------------------------------------------------------------
void IBDGraph::GraphLoadGXL( const char * filename,
			     int (*readNodeLabel)(TiXmlElement *elem),
			     double (*readEdgeLabel)(TiXmlElement *elem))
{
  std::ifstream file(filename,std::ios::in);
  std::vector<char*> v;
  //XXX: find gxl property for this point
  _directed = false;
  TiXmlDocument doc(filename);
  if(!doc.LoadFile()){
    std::cerr << "Error while loading file" << std::endl;
    std::cerr << "error #" << doc.ErrorId() << " : " << doc.ErrorDesc() << std::endl;
  }
   
  TiXmlHandle hdl(&doc);
  std::map<int,int> id_to_index;
  TiXmlElement *elem = hdl.FirstChildElement().FirstChildElement("graph").FirstChildElement().Element();
  while (elem){
    if(strcmp(elem->Value(),"node") == 0){
      int id = std::stoi( elem->Attribute("id")+1);
      int label = readNodeLabel(elem);
      id_to_index[id] = nbNodes;
      Add(new GNode<int, double>(id,label));
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
