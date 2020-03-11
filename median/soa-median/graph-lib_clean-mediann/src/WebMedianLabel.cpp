
#include "WebMedianLabel.h"
#include <cmath>
#include <string>



WebNAtt WebMedianLabel::MedianNodeLabel(int * * mappingsFromMedian, int node1, Dataset<WebNAtt, WebEAtt, int> * dataset){
//std::cout << " Computation of median label for vertex " <<  node1  << std::endl;
//std::cout << "\n --- calcul de label moyen pour sommet " << node1 << " du médian \n";
//std::cout<<"MedianNode computation\n";
int maxFrequence = 0;
int N= dataset->size();
std::map<std::string,int> m;
std::map<std::string,int>::iterator it;
std::map<std::string,double> m2;
std::map<std::string,double>::iterator it2;
std::string label;
for (int i=0;i<N;i++){
  int mappedIndex = mappingsFromMedian[i][node1];
  if (mappedIndex < (*dataset)[i]-> Size()){
  	label =  (*(*dataset)[i])[mappedIndex]->attr.id  ;
  	//std::cout << "label sommet " << node1 << " dans graphe " << i << " est " << label << std::endl;
  	it = m.find(label);
  	if (it == m.end()) {// label non présent
     		m[label] = 1;
     		m2[label] = (*(*dataset)[i])[mappedIndex]->attr.freq  ;
     		}
     		       // ajout d'un nouveau label avec une fréquence de 1
  	else{
        it->second++;
        it2=m2.find(label);
        it2->second += (*(*dataset)[i])[mappedIndex]->attr.freq  ;
        }
  }       // incrémentation de la fréquence du label
}

std::string MedianId;
for (it=m.begin();it!=m.end();it++){
if (it->second > maxFrequence){
   maxFrequence = it->second;
   MedianId = it->first;
}
}
WebNAtt medianLabel;

medianLabel.id = MedianId;
medianLabel.freq = m2[MedianId]/static_cast<double>(maxFrequence);

//std::cout<<"MedianNode computed\n";
return medianLabel;

};

WebEAtt WebMedianLabel::MedianEdgeLabel(Graph<WebNAtt,WebEAtt> * Gbar, int * * mappingsFromMedian, int node1, int node2, Dataset<WebNAtt, WebEAtt, int> * dataset){
//std::cout<<"MedianEdge computation\n";
int maxFrequence = 0;
int N= dataset->size();
WebEAtt absentEdge = std::numeric_limits<WebEAtt>::max();
std::map<std::string,int> m;
std::map<std::string,int>::iterator it;
std::map<std::string,double> m2;
std::map<std::string,double>::iterator it2;
WebEAtt label;
//std::cout << "initialiazed \n" ;
for (int i=0;i<N;i++){
  int mappedIndex1 = mappingsFromMedian[i][node1];
   int mappedIndex2 = mappingsFromMedian[i][node2];
    //std::cout << "("<< node1 <<","<<node2<<") in graph with size " << Gbar->Size() << " mapped to ("<< mappedIndex1 << "," << mappedIndex2 <<"). in graph with size " << (*dataset)[i]-> Size() << std::endl;
  if (mappedIndex1 < (*dataset)[i]-> Size() && mappedIndex2 < (*dataset)[i]-> Size() && (*dataset)[i]->isLinked(mappedIndex1,mappedIndex2)){
    //std::cout << "linked \n";
  	label =  ((*dataset)[i]->getEdge(mappingsFromMedian[i][node1],mappingsFromMedian[i][node2]))->attr;
  	}
  	else
  	{
  //	std::cout << "not linked \n";
  	label = absentEdge;
  	}
  	//std::cout << "label sommet " << node1 << " dans graphe " << i << " est " << label << std::endl;
  	it = m.find(label.id);
  	if (it == m.end()) {// label non présent
     		m[label.id] = 1;
     		m2[label.id] =   label.val  ;
     		}
     		       // ajout d'un nouveau label avec une fréquence de 1
  	else{
        it->second++;
        it2=m2.find(label.id);
        it2->second +=   label.val  ;
        }
  }       // incrémentation de la fréquence du label
//std::cout << "labels separated ! \n";
std::string MedianId;
for (it=m.begin();it!=m.end();it++){
if (it->second > maxFrequence){
   maxFrequence = it->second;
   MedianId = it->first;
}
}
WebEAtt medianLabel;

medianLabel.id = MedianId;
medianLabel.val = m2[MedianId]/static_cast<double>(maxFrequence);

//std::cout<<"MedianEdge computed\n";
return medianLabel;
};

double WebMedianLabel::NodeDelta(Graph<WebNAtt,WebEAtt> * Gbar,int * * mappingsFromMedian, int nodeIndex, Dataset<WebNAtt, WebEAtt, int> * dataset){

int N= dataset->size();
double actual_cost =0;
double potential_cost =0;

for (int i=0;i<N;i++){
  int mappedIndex = mappingsFromMedian[i][nodeIndex];
  if (mappedIndex < (*dataset)[i]-> Size()){
  	actual_cost += this->cf->NodeSubstitutionCost((* Gbar)[nodeIndex],(*(*dataset)[i])[mappedIndex],Gbar,(*dataset)[i]);
  	potential_cost += this->cf->NodeInsertionCost((*(*dataset)[i])[mappedIndex],(*dataset)[i]);
  	for (int j=0;j<N;j++){
        GEdge<WebEAtt> * EdgeInGbar=Gbar->getEdge(nodeIndex,j);
        if (EdgeInGbar){
            GEdge<WebEAtt> * EdgeInDatasetGraph=(*dataset)[i]->getEdge(mappedIndex,mappingsFromMedian[i][j]);
            if (EdgeInDatasetGraph){//edge exists in both graphs
                actual_cost += this->cf->EdgeSubstitutionCost(EdgeInGbar,EdgeInDatasetGraph,Gbar,(*dataset)[i]);
                potential_cost += this->cf->EdgeInsertionCost(EdgeInDatasetGraph,(*dataset)[i]);
            }
            else{
                actual_cost += this->cf->EdgeDeletionCost(EdgeInGbar,Gbar);
            }
        }
  	}
  }
  else actual_cost += this->cf->NodeDeletionCost((* Gbar)[nodeIndex],Gbar);

}
 //std::cout << "--- médian a les coordonnées ("<< MedianLabel.x << "," << MedianLabel.y <<")\n";
//std::cout << "label sommet " << node1 << " dans graphe  median est " << Medianlabel << std::endl;

return potential_cost - actual_cost;

};


