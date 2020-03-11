
#include "ConstantMedianLabel.h"
 #include <map>



int ConstantMedianLabel::MedianNodeLabel(int * * mappingsFromMedian, int node1, Dataset<int,int,double> * dataset){
//std::cout << " Computation of median label for vertex " <<  node1  << std::endl;
int maxFrequence = 0;
int Medianlabel = 0;
int N= dataset->size();
std::map<int,int> m;
std::map<int,int>::iterator it;
int label;
for (int i=0;i<N;i++){
  int mappedIndex = mappingsFromMedian[i][node1];
  if (mappedIndex < (*dataset)[i]-> Size()){
  	label =  (*(*dataset)[i])[mappedIndex]->attr  ;
  	//std::cout << "label sommet " << node1 << " dans graphe " << i << " est " << label << std::endl;
  	it = m.find(label);
  	if (it == m.end()) // label non présent
     		m[label] = 1;       // ajout d'un nouveau label avec une fréquence de 1
  	else it->second++;
  }       // incrémentation de la fréquence du label
}

for (it=m.begin();it!=m.end();it++){
if (it->second > maxFrequence){
   maxFrequence = it->second;
   Medianlabel = it->first;
}
}

//std::cout << "label sommet " << node1 << " dans graphe  median est " << Medianlabel << std::endl;

return Medianlabel;

};
int ConstantMedianLabel::MedianEdgeLabel(Graph<int,int> * medianGraph, int * * mappingsFromMedian, int node1, int node2, Dataset<int,int,double> * dataset){

//std::cout << " Computation of median edge label for edge " <<  node1 << node2 << std::endl;

  double ces=this->cf->ces();
  double cei=this->cf->cei();
  double ced=this->cf->ced();


bool isConnected;
int absentEdge = std::numeric_limits<int>::max();
int maxFrequenceNotAbsent = 0;
int MedianlabelNotAbsent = 0;
int N= dataset->size();
std::map<int,int> m;
std::map<int,int>::iterator it;
int label;

for (int i=0;i<N;i++){
  if (mappingsFromMedian[i][node1] < (*dataset)[i]->Size() && mappingsFromMedian[i][node2] < (*dataset)[i]->Size()){
  isConnected = (*dataset)[i]->isLinked(mappingsFromMedian[i][node1],mappingsFromMedian[i][node2]);
  if (isConnected)
  label =  ((*dataset)[i]->getEdge(mappingsFromMedian[i][node1],mappingsFromMedian[i][node2]))->attr;
  else
  label = absentEdge;
}
  else
  label = absentEdge;
  //std::cout << "label arete " << node1 <<","<< node2 << " dans graphe " << i << " est " << label << std::endl;
  it = m.find(label);
  if (it == m.end()) // label non présent
     m[label] = 1;       // ajout d'un nouveau label avec une fréquence de 1
  else it->second++;          // incrémentation de la fréquence du label
}

for (it=m.begin();it!=m.end();it++){
if (it->second > maxFrequenceNotAbsent && it->first!= absentEdge){
   maxFrequenceNotAbsent = it->second;
   MedianlabelNotAbsent = it->first;
}
}


if (maxFrequenceNotAbsent < N*(1-(cei/ces)) - m[absentEdge]*(1 - (ced + cei) / ces) )
{
  //std::cout << "arete " <<  node1 <<","<< node2 <<  " est absente dans median graphe " << std::endl;
  return absentEdge ;
}
else
{
  //std::cout << "arete " <<  node1 <<","<< node2 <<  " dans median graphe a le label" << MedianlabelNotAbsent << std::endl;
  return MedianlabelNotAbsent ;

}
};

double ConstantMedianLabel::NodeDelta(Graph<int,int> * Gbar,int * * mappingsFromMedian, int nodeIndex, Dataset<int,int,double> * dataset){

int N= dataset->size();
double actual_cost =0;
double potential_cost =0;
//std::cout << "computing NodeDelta node " << nodeIndex << std::endl;
for (int i=0;i<N;i++){
  int mappedIndex = mappingsFromMedian[i][nodeIndex];
  if (mappedIndex < (*dataset)[i]-> Size()){
  	actual_cost += this->cf->NodeSubstitutionCost((* Gbar)[nodeIndex],(*(*dataset)[i])[mappedIndex],Gbar,(*dataset)[i]);
  	potential_cost += this->cf->NodeInsertionCost((*(*dataset)[i])[mappedIndex],(*dataset)[i]);
  	for (int j=0;j<Gbar->Size();j++){
        GEdge<int> * EdgeInGbar=Gbar->getEdge(nodeIndex,j);
        if (EdgeInGbar){
            GEdge<int>* EdgeInDatasetGraph=(*dataset)[i]->getEdge(mappedIndex,mappingsFromMedian[i][j]);
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
//std::cout << "potential = " << potential_cost<< " , actual = "<< actual_cost << std::endl;
return potential_cost - actual_cost;

};


double ConstantMedianLabel::NodeLabelDelta(int median_size, int * * mappingsToMedian, int ** mappingsFromMedian, Graph<int,int>* &Gbar, int * nodeLabel, Dataset<int,int,double> * dataset){
    double cvs=this->cf->cns();
    double cvi=this->cf->cni();
    double cvd=this->cf->cnd();
    double DeltaCost = 0;
    int label;
    int N= dataset->size();
    std::map<int,int> num_ins_by_label;
    std::map<int,int>::iterator it;
    std::map<int,bool> ins_found_by_label;
    std::map<int,bool>::iterator it2;
    int num_edit_path_with_insertions = 0;
    for (int i=0;i<N;i++){
        for (it2=ins_found_by_label.begin();it2!=ins_found_by_label.end();it2++){
            it2->second = false;
        }
        bool has_insertions = false;
        int graph_size=dataset->getGraph(i)->Size();
        for (int j=0;j<graph_size;j++){
            if (mappingsToMedian[i][j]>= median_size){
                has_insertions = true;
                label =  (*(*dataset)[i])[j]->attr;
                it = num_ins_by_label.find(label);
                if (it == num_ins_by_label.end()){ // label non présent
                    num_ins_by_label[label] = 1;       // ajout d'un nouveau label avec une fréquence de 1
                    ins_found_by_label[label] = true;
                }
                else if (!ins_found_by_label[label]){
                    it->second++;
                    ins_found_by_label[label] = true;
                }
            }
        }
        num_edit_path_with_insertions += has_insertions;
    }
    int maxFrequence=0;
    for (it=num_ins_by_label.begin();it!=num_ins_by_label.end();it++){
        if (it->second > maxFrequence){
            maxFrequence = it->second;
            *nodeLabel = it->first;
        }
    }

    //std::cout << "taille du dataset = " << N << " , nb graphs with insertions = " << num_edit_path_with_insertions<< " , maxFrequence = " << maxFrequence <<  " of label " << (*nodeLabel)<<std::endl;
    DeltaCost -= cvi * maxFrequence; // vertices had to be inserted before, and now must only be substituted with same label.
    DeltaCost -= (cvi - cvs) * (num_edit_path_with_insertions - maxFrequence); // vertices had to be inserted before, and now must only be substituted with different label.
    DeltaCost += cvd * (N - num_edit_path_with_insertions);// vertex in median must now be deleted.
     if (DeltaCost < 0){
        int nb_same_label =0;
        int nb_insertions =0;
        GNode<int,int> * Newnode = new GNode<int,int>(median_size,(*nodeLabel));
        // à remettre
            //  delete Gbar;
        Gbar->Add(Newnode);
        //std::cout<<"adding node " << Gbar->Size()-1 << " in MedianGraph, with DeltaCost = " << DeltaCost << " , and Label = " << (*nodeLabel) << std::endl;
        int M = dataset->size();
        for(int i=0;i<M;i++){
            int size_graph_i = dataset->getGraph(i)->Size();
            //updating mappingsToMedian so that vertices mapped to nMedian (supposed to be deleted) are still deleted w.r.t. a bigger median graph.
            int best_assigned_vertex = size_graph_i;//initialize as deletion
            bool same_label_found = false;
            bool insertion_found = false;
            for (int j=0;j<size_graph_i;j++){
                if (mappingsToMedian[i][j] == median_size)
                    mappingsToMedian[i][j] = median_size+1;
                if  (!same_label_found && mappingsToMedian[i][j] >= median_size){
                    insertion_found = true; // node i in graph j is inserted
                    int label = (*(dataset->getGraph(i)))[j]->attr;
                    if(label == (*nodeLabel)){
                        best_assigned_vertex = j;
                        same_label_found = true;
                    }
                    else if (cvs <= cvd)
                        best_assigned_vertex = j;
                }
            }
            //if (insertion_found) nb_insertions++;
            //if (same_label_found) nb_same_label++;
            //reallocating memory for new mappingsFromMedian in order to take in account the new size of Gmed.

            int * UpdatedMappingsFromMedian_to_i = new int [Gbar->Size()];
            for(int j=0;j<Gbar->Size()-1;j++){
                UpdatedMappingsFromMedian_to_i[j]=mappingsFromMedian[i][j];
            }
            //newly added vertex is assigned to inserted vertex with same label, or inserted vertex with different label, or epsilon.
            UpdatedMappingsFromMedian_to_i[Gbar->Size()-1]=best_assigned_vertex;
            //std::cout << "pre last update of mappings to graph " << i << std::endl;
            //PrintPointerContent(mappingsToMedian[i],this->ds->getGraph(i)->Size());
            if (best_assigned_vertex <size_graph_i){
              //  std::cout << "best_assigned_vertex = " << best_assigned_vertex << std::endl;
                mappingsToMedian[i][best_assigned_vertex]=Gbar->Size()-1;}
            // à remettre
            //delete[] mappingsFromMedian[i];
            mappingsFromMedian[i]=UpdatedMappingsFromMedian_to_i;
            //std::cout << "post last update of mappings to graph " << i << std::endl;
            //PrintPointerContent(mappingsToMedian[i],this->ds->getGraph(i)->Size());
            }

        //std::cout<<"mappings updated " << std::endl;
        //std::cout<<"updated median created " << std::endl;
        //std::cout<<"old median -->"<< std::endl;
        //printGraph(*(Gbar));
         //std::cout<<"node deleted"<< std::endl;
        //std::cout << "taille du dataset = " << M << " , nb graphs with insertions = " << nb_insertions << " , maxFrequence = " << nb_same_label << std::endl;
    }


    return DeltaCost;
};

 int ConstantMedianLabel::WeightedEdgeMeanLabel(int label1, int label2, double alpha){
    if (alpha < this->cf->ces()/2)
        return label1;
    else
        return label2;

};
int ConstantMedianLabel::WeightedVertexMeanLabel(int label1, int label2, double alpha){
   if (alpha < this->cf->cns()/2)
        return label1;
    else
        return label2;
};
