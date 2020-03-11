
#include "CMUMedianLabel.h"
#include <cmath>



CMUPoint CMUMedianLabel::MedianNodeLabel(int * * mappingsFromMedian, int node1, Dataset<CMUPoint,double,int> * dataset){
//std::cout << " Computation of median label for vertex " <<  node1  << std::endl;
CMUPoint MedianLabel;
double sumX = 0;
double sumY = 0;
int N= dataset->size();

for (int i=0;i<N;i++){
  int mappedIndex = mappingsFromMedian[i][node1];
  if (mappedIndex < (*dataset)[i]-> Size()){
  	CMUPoint localPoint =  (*(*dataset)[i])[mappedIndex]->attr  ;
  	//std::cout << "label sommet " << node1 << " dans graphe " << i << " est " << label << std::endl;
  	sumX+= localPoint.x;
        sumY+= localPoint.y;
    }
}

MedianLabel.x = sumX/N;
MedianLabel.y = sumY/N;

//std::cout << "label sommet " << node1 << " dans graphe  median est " << Medianlabel << std::endl;

return MedianLabel;

};

double CMUMedianLabel::MedianEdgeLabel(Graph<CMUPoint,double> * Gbar, int * * mappingsFromMedian, int node1, int node2, Dataset<CMUPoint,double,int> * dataset){

//std::cout << " Computation of median edge label for edge " <<  node1 << node2 << std::endl;

Graph<CMUPoint,double> *MedianGraph = new Graph<CMUPoint,double>(*Gbar);


bool isConnected;
double absentEdge = std::numeric_limits<double>::max();
double labelIfPresent = std::pow(std::pow((*MedianGraph)[node1]->attr.x - (*MedianGraph)[node2]->attr.x, 2) + std::pow((*MedianGraph)[node1]->attr.y - (*MedianGraph)[node2]->attr.y, 2), 1/2);
double costIfLinked = 0;
double costIfNotLinked = 0;
int N= dataset->size();
/*
int MedianlabelNotAbsent = 0;
int N= dataset->size();
std::map<int,int> m;
std::map<int,int>::iterator it;
double label;
*/

if (!(MedianGraph->isLinked(node1,node2)))
MedianGraph->Link(node1,node2,labelIfPresent);

for (int i=0;i<N;i++){
  if (!(*dataset)[i]->isLinked(mappingsFromMedian[i][node1],mappingsFromMedian[i][node2]))
  costIfLinked+= this->cf->EdgeDeletionCost(MedianGraph->getEdge(node1,node2),MedianGraph);
  else
  costIfLinked+= this->cf->EdgeSubstitutionCost(MedianGraph->getEdge(node1,node2),(*dataset)[i]->getEdge(mappingsFromMedian[i][node1],mappingsFromMedian[i][node2]),MedianGraph,(*dataset)[i]);
}

MedianGraph->Unlink(node1,node2);

for (int i=0;i<N;i++){
  if ((*dataset)[i]->isLinked(mappingsFromMedian[i][node1],mappingsFromMedian[i][node2]))
  costIfNotLinked+= this->cf->EdgeInsertionCost((*dataset)[i]->getEdge(mappingsFromMedian[i][node1],mappingsFromMedian[i][node2]),(*dataset)[i]);
}

delete MedianGraph;
if (costIfLinked < costIfNotLinked)
return labelIfPresent;
else
return absentEdge;
};
  

