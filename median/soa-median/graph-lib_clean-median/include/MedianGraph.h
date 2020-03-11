#ifndef __MEDIANGRAPH_H__
#define __MEDIANGRAPH_H__

#include "graph.h"
#include "Dataset.h"
#include "GraphEditDistance.h"

template<typename T1, typename T2>
void PrintPointerContent(T1 * pointer, T2 length){
    std::cout<<"[";
    for (int j = 0; j < length;j++){
        std::cout << pointer[j] << " ";
    }
    std::cout<< "]"<<std::endl;
}

template<class NodeAttribute, class EdgeAttribute,class PropertyType>
class MedianLabel{


public :


  virtual NodeAttribute MedianNodeLabel(int * * mappingsFromMedian, int nodeIndex, Dataset<NodeAttribute, EdgeAttribute, PropertyType> * ds)=0;

  virtual EdgeAttribute MedianEdgeLabel(Graph<NodeAttribute, EdgeAttribute> * medianGraph, int * * mappingsFromMedian, int node1Index, int node2Index, Dataset<NodeAttribute, EdgeAttribute, PropertyType> * ds)=0;

  virtual double NodeDelta(Graph<NodeAttribute, EdgeAttribute> * medianGraph,int * * mappingsFromMedian, int nodeIndex, Dataset<NodeAttribute, EdgeAttribute, PropertyType> * ds)=0;

  // NodeDelta computes the cost of not having the vertex in the median w.r.t to having it.

  virtual double NodeLabelDelta(int median_size, int * * mappingsToMedian, int * * mappingsFromMedian, Graph<NodeAttribute, EdgeAttribute> *& medianGraph,NodeAttribute* NodeLabel, Dataset<NodeAttribute, EdgeAttribute, PropertyType> * ds)=0;

  // NodeDelta computes the cost of having a new vertex with given label w.r.t to not having it,
  // it also updates the median graph passed as a pointer, and the assignments, if Delta is negative.

  virtual EdgeAttribute WeightedEdgeMeanLabel(EdgeAttribute label1, EdgeAttribute label2, double alpha)=0;

  virtual NodeAttribute WeightedVertexMeanLabel(NodeAttribute label1, NodeAttribute label2, double alpha)=0;


};

template<class NodeAttribute, class EdgeAttribute, class PropertyType>
class MedianGraph
{
protected:
  MedianLabel<NodeAttribute,EdgeAttribute,PropertyType> * ml;
  Graph<NodeAttribute,EdgeAttribute> * Gbar;
  Dataset<NodeAttribute, EdgeAttribute, PropertyType> * ds;

public:

  //Mapping is an array encoding the mapping of each node
  MedianGraph(Graph<NodeAttribute,EdgeAttribute> * g, Dataset< NodeAttribute,EdgeAttribute, PropertyType> * dataset,MedianLabel<NodeAttribute,EdgeAttribute,PropertyType> * medianlabel ):ml(medianlabel),ds(dataset),Gbar(g){};

  bool updateMedianGraph(int * * mappingsFromMedian);
  bool checkNodeDeletion(int * * mappingsFromMedian, int ** mappingsToMedian,int nMedian);
  bool checkNodeInsertion(int * * mappingsFromMedian, int ** mappingsToMedian,int nMedian);
  Graph<NodeAttribute,EdgeAttribute>* getGraph(){return Gbar;};
  Graph<NodeAttribute, EdgeAttribute>* ComputeWeightedMeanGraph(Graph<NodeAttribute,EdgeAttribute> * G1, Graph<NodeAttribute,EdgeAttribute> * G2, int * G1toG2, int * G2toG1, double alpha, GraphEditDistance<NodeAttribute, EdgeAttribute> * ed);
  Graph<NodeAttribute, EdgeAttribute>* ComputeWeightedMeanGraph(int indexGraph1, int indexGraph2, int * G1toG2, int * G2toG1, double alpha, GraphEditDistance<NodeAttribute, EdgeAttribute> * ed);


  ~MedianGraph(){
   delete ml;
   delete Gbar;
   delete ds;
   };
};
template<class NodeAttribute, class EdgeAttribute, class PropertyType>
bool MedianGraph<NodeAttribute,EdgeAttribute,PropertyType>::updateMedianGraph(int * * mappingsFromMedian){
//std::cout << "call to updateMedianGraph" ;
bool isModified = false;
EdgeAttribute absentEdge = std::numeric_limits<EdgeAttribute>::max();
GEdge<EdgeAttribute> * currentEdge = NULL;
EdgeAttribute newEdgeAttribute;
NodeAttribute newNodeAttribute;
int currentNbEdges;
int N = ds->size();
int n = Gbar->Size();
//GNode<NodeAttribute,EdgeAttribute>** assignedNodes = new int*[n];
NodeAttribute medianAttribute;
for (int i=0;i<n;i++){
   newNodeAttribute = ml-> MedianNodeLabel(mappingsFromMedian,i,this->ds); // updating label of node i
   if (newNodeAttribute != (*Gbar)[i]->attr){
   //std::cout << "Label Change, vertex " << i << std::endl;
   (*Gbar)[i]->attr = newNodeAttribute;
   isModified = true;

   }
}

for (int i=0;i<n;i++){
   for (int j=i+1;j<n;j++){
      newEdgeAttribute = ml->MedianEdgeLabel(Gbar,mappingsFromMedian,i,j,this->ds);

      if (newEdgeAttribute!=absentEdge) {
          if (Gbar->isLinked(i,j)){
              currentEdge = Gbar->getEdge(i,j);
              if (currentEdge->attr != newEdgeAttribute){
              //std::cout << "Label Change, edge (" << i <<","<<j<< ") from "<< currentEdge->attr << " to " << newEdgeAttribute << std::endl;
              currentEdge->attr = newEdgeAttribute; //updating label of edge (i,j)
              Gbar->getEdge(j,i)->attr = newEdgeAttribute;
              isModified = true;

              }
          }
          else{
              Gbar->Link(i,j,newEdgeAttribute); // creating the edge with the new label
              //(*Gbar)[i]->Connect(j,newEdgeAttribute);
              //(*Gbar)[j]->Connect(i,newEdgeAttribute);
              isModified = true;
              //std::cout << "EdgeInsertion, (" << i <<","<<j<< ")" <<  std::endl;
              if (!Gbar->getEdge(j,i)){
                    Gbar->Link(j,i,newEdgeAttribute);
              }
              }
      }
      else {
         if (Gbar->isLinked(i,j)){
              Gbar->Unlink(i,j); // removing the edge if newEdgeAttribute == NULL
              //(*Gbar)[i]->UnConnect(j);
              //(*Gbar)[j]->UnConnect(i);
              isModified = true;
              //std::cout << "EdgeDeletion, (" << i <<","<<j<< ")" << std::endl;
         }
      }
    }

  }

//std::cout << " ----> "<< isModified <<std::endl;
return isModified;

};

template<class NodeAttribute, class EdgeAttribute, class PropertyType>
bool MedianGraph<NodeAttribute,EdgeAttribute,PropertyType>::checkNodeDeletion(int * * mappingsFromMedian, int ** mappingsToMedian,int nMedian){
    int M = this->ds->size();
    bool isModified = false;
    double minDelta{0};
    int indMinDelta{0};
    double * deltas = new double [nMedian];
    for (int i=0;i<nMedian;i++){
        deltas[i]=ml->NodeDelta(Gbar,mappingsFromMedian,i,this->ds);
        //std::cout<<"Delta node " << i<< " = " << deltas[i] << std::endl;
        if (i==0 || deltas[i]<minDelta){
            minDelta = deltas[i];
            indMinDelta = i;
        }
    }

    if (minDelta < 0){

        //std::cout<<"deleting node " << indMinDelta<< " in MedianGraph, with DeltaCost = " << minDelta << std::endl;

            for(int i=0;i<M;i++){
                for (int j=indMinDelta;j<nMedian-1;j++){
                    mappingsFromMedian[i][j]=mappingsFromMedian[i][j+1];
                }
                int size_graph_i = this->ds->getGraph(i)->Size();
                for (int j=0;j<size_graph_i;j++){
                    if (mappingsToMedian[i][j] == indMinDelta)
                        mappingsToMedian[i][j] = nMedian+2;
                    if (mappingsToMedian[i][j] > indMinDelta)
                        mappingsToMedian[i][j] = mappingsToMedian[i][j]-1;
                }
            }
        //std::cout<<"mappings updated " << std::endl;
        Graph<NodeAttribute,EdgeAttribute>* new_median = new Graph<NodeAttribute,EdgeAttribute>(*(this->Gbar),indMinDelta);
        //std::cout<<"updated median created " << std::endl;
        //std::cout<<"old median -->"<< std::endl;
        //printGraph(*(Gbar));
         //std::cout<<"new median -->"<< std::endl;
         //printGraph(*(new_median));
        delete Gbar;
        Gbar = new_median;
        isModified = true;
         //std::cout<<"node deleted"<< std::endl;
    }
    //à remettre
    //delete[] deltas;
    return isModified;
}

template<class NodeAttribute, class EdgeAttribute, class PropertyType>
bool MedianGraph<NodeAttribute,EdgeAttribute,PropertyType>::checkNodeInsertion(int * * mappingsFromMedian, int ** mappingsToMedian,int nMedian){
    /*
    int N = this->ds->size();
    std::cout << "mappingsfrom =" << std::endl;
    for (int i = 0; i<N;i++){
            PrintPointerContent(mappingsFromMedian[i],nMedian);
    }
  std::cout << "mappingsTo =" << std::endl;
    for (int i = 0; i<N;i++){
            PrintPointerContent(mappingsToMedian[i],this->ds->getGraph(i)->Size());
    }
*/

    bool isModified = false;
    NodeAttribute * nodeLabel = new NodeAttribute;
    double Delta = ml->NodeLabelDelta(nMedian,mappingsToMedian, mappingsFromMedian,this->Gbar,nodeLabel,this->ds);
    //std::cout << "best Delta for insertion of vertex in Median = " << Delta << std::endl;
    nMedian=this->Gbar->Size();
    /*
      std::cout << "mappingsfrom =" << std::endl;
    for (int i = 0; i<N;i++){
            PrintPointerContent(mappingsFromMedian[i],nMedian);
    }
  std::cout << "mappingsTo =" << std::endl;
    for (int i = 0; i<N;i++){
            PrintPointerContent(mappingsToMedian[i],this->ds->getGraph(i)->Size());
    }
    */
    delete nodeLabel;
    if (Delta < 0) return true;
    return false;
}
template<class NodeAttribute, class EdgeAttribute, class PropertyType>
Graph<NodeAttribute, EdgeAttribute>* MedianGraph<NodeAttribute,EdgeAttribute,PropertyType>::ComputeWeightedMeanGraph(Graph<NodeAttribute,EdgeAttribute> * G1, Graph<NodeAttribute,EdgeAttribute> * G2, int * G1toG2, int * G2toG1, double alpha, GraphEditDistance<NodeAttribute, EdgeAttribute> * ed){
    int n = G1->Size();
    int m = G2->Size();
    std::vector<std::vector<struct GraphEditDistance<NodeAttribute, EdgeAttribute>::EditOperation*>> operations;
    operations = ed->EditPathFromMapping(G1,G2,G1toG2,n,G2toG1,m);
    //std::cout<<"Complete Edit Path : " << std::endl;
    for (int i=0;i<operations.size();i++){
       // std::cout<<"       Operation " <<  i<< std::endl;
        for (int j = 0;j<operations[i].size();j++){
         //   std::cout << "       v1_i = "<< operations[i][j]->v1_i << " , v1_j = "<< operations[i][j]->v1_j << " , v2_i = "<< operations[i][j]->v2_i << " , v2_j = "<< operations[i][j]->v2_j << " , cost = "<< operations[i][j]->cost<< std::endl;

        }
    }
    //std::cout<<"Effective Edit Path : " << std::endl;
    Graph<NodeAttribute, EdgeAttribute>* Gmean = new Graph<NodeAttribute, EdgeAttribute>(*G1);
    //std::cout << "ok 1" << std::endl;
    double total_edit_cost=0;
    int ind_operation = 0;
    //-----------//
    // edge deletions
    //-----------//
    while ((ind_operation<operations[0].size()) && (total_edit_cost + operations[0][ind_operation]->cost < alpha)){
        // application operation
        Gmean->Unlink(operations[0][ind_operation]->v1_i,operations[0][ind_operation]->v1_j);
        //std::cout << "       v1_i = "<< operations[0][ind_operation]->v1_i << " , v1_j = "<< operations[0][ind_operation]->v1_j << " , v2_i = "<< operations[0][ind_operation]->v2_i << " , v2_j = "<< operations[0][ind_operation]->v2_j << std::endl;
        //incrementation cost and index of edit-operation
        total_edit_cost += operations[0][ind_operation]->cost;
        ind_operation ++;
        //std::cout << "       current cost  = " << total_edit_cost << " , objective = " << alpha << std::endl;
    }
    if (ind_operation<operations[0].size()){
        // necessary decision whether the current operation must be applied or not
        //std::cout<< "ok 1.0";
        if (alpha - total_edit_cost > operations[0][ind_operation]->cost/2.0){
        //std::cout<< "ok 1.1";
            Gmean->Unlink(operations[0][ind_operation]->v1_i,operations[0][ind_operation]->v1_j);
            //std::cout << "       v1_i = "<< operations[0][ind_operation]->v1_i << " , v1_j = "<< operations[0][ind_operation]->v1_j << " , v2_i = "<< operations[0][ind_operation]->v2_i << " , v2_j = "<< operations[0][ind_operation]->v2_j << std::endl;
            //std::cout << "       final operation " << std::endl;
        }
         //std::cout<< "ok 1.2";
        return Gmean;
    }
    //std::cout << "ok 2" << std::endl;
    //-----------//
    //node deletions
    //-----------//
    ind_operation = 0;
    while (ind_operation<operations[1].size() && total_edit_cost + operations[1][ind_operation]->cost < alpha){
        //std::cout << "ok 2.00"<< std::endl;
        // application operation
        int v1= operations[1][ind_operation]->v1_i;
        Graph<NodeAttribute, EdgeAttribute>* Gmean_temp = new Graph<NodeAttribute, EdgeAttribute>(*Gmean,v1);  //constructeur de copie qui exlue le sommet v1_i;
        delete Gmean;
        Gmean=Gmean_temp;
        //incrementation cost and index of edit-operation
        total_edit_cost += operations[1][ind_operation]->cost;
        //std::cout << "       v1_i = "<< operations[1][ind_operation]->v1_i << " , v1_j = "<< operations[1][ind_operation]->v1_j << " , v2_i = "<< operations[1][ind_operation]->v2_i << " , v2_j = "<< operations[1][ind_operation]->v2_j << std::endl;
        //std::cout << "       current cost  = " << total_edit_cost << " , objective = " << alpha << std::endl;
        //update of edit operations with new indices
          for (int type_op=0;type_op<6;type_op++){
            for (int ind_op=0;ind_op<operations[type_op].size();ind_op++){
                if (operations[type_op][ind_op]->v1_i>v1) operations[type_op][ind_op]->v1_i --;
                if (operations[type_op][ind_op]->v1_j>v1) operations[type_op][ind_op]->v1_j --;
            }
        }
        //std::cout << "       vertex deleted, indices updated" << std::endl;
        ind_operation ++;
    }
    //std::cout << "ok 2.01"<< std::endl;
    if (ind_operation<operations[1].size()){
    //std::cout<< "ok 2.0"<< std::endl;
        // necessary decision whether the current operation must be applied or not
        if (alpha - total_edit_cost > operations[1][ind_operation]->cost/2){
        //std::cout<< "ok 2.1"<< std::endl;
             int v1= operations[1][ind_operation]->v1_i;
            Graph<NodeAttribute, EdgeAttribute>* Gmean_temp = new Graph<NodeAttribute, EdgeAttribute>(*Gmean,v1);    //constructeur de copie qui exlue le sommet v1_i;
            delete Gmean;
            Gmean=Gmean_temp;
            //std::cout << "       v1_i = "<< operations[1][ind_operation]->v1_i << " , v1_j = "<< operations[1][ind_operation]->v1_j << " , v2_i = "<< operations[1][ind_operation]->v2_i << " , v2_j = "<< operations[1][ind_operation]->v2_j << std::endl;
            //std::cout << "       final operation   " << std::endl;
        }
        //std::cout<< "ok 2.2"<< std::endl;
        return Gmean;
    }
    //std::cout << "ok 3" << std::endl;
    //-----------//
    //edge substitutions
    //-----------//
    ind_operation = 0;
    while (ind_operation<operations[2].size() && total_edit_cost + operations[2][ind_operation]->cost < alpha){
        // application operation
        //std::cout << "       v1_i = "<< operations[2][ind_operation]->v1_i << " , v1_j = "<< operations[2][ind_operation]->v1_j << " , v2_i = "<< operations[2][ind_operation]->v2_i << " , v2_j = "<< operations[2][ind_operation]->v2_j << std::endl;
        GEdge<EdgeAttribute>* e1 = Gmean->getEdge(operations[2][ind_operation]->v1_i,operations[2][ind_operation]->v1_j);
        GEdge<EdgeAttribute>* e2 = G2->getEdge(operations[2][ind_operation]->v2_i,operations[2][ind_operation]->v2_j);
        e1->attr = e2->attr;
        //incrementation cost and index of edit-operation
        total_edit_cost += operations[2][ind_operation]->cost;
        //std::cout << "       current cost  = " << total_edit_cost << " , objective = " << alpha << std::endl;
        //update of edit operations with new indices
        ind_operation ++;

    }
    if (ind_operation<operations[2].size()){
        // necessary decision whether the current operation must be applied or not
        GEdge<EdgeAttribute>* e1 = Gmean->getEdge(operations[2][ind_operation]->v1_i,operations[2][ind_operation]->v1_j);
        GEdge<EdgeAttribute>* e2 = G2->getEdge(operations[2][ind_operation]->v2_i,operations[2][ind_operation]->v2_j);
        e1->attr = ml->WeightedEdgeMeanLabel(e1->attr,e2->attr,alpha - total_edit_cost);
        //std::cout << "       v1_i = "<< operations[2][ind_operation]->v1_i << " , v1_j = "<< operations[2][ind_operation]->v1_j << " , v2_i = "<< operations[2][ind_operation]->v2_i << " , v2_j = "<< operations[2][ind_operation]->v2_j << std::endl;
        //std::cout << "       final operation " << std::endl;
        return Gmean;
    }
    //std::cout << "ok 4" << std::endl;
    //-----------//
    //vertex substitution
    //-----------//
    ind_operation = 0;
    while (ind_operation<operations[3].size() && total_edit_cost + operations[3][ind_operation]->cost < alpha){
        // application operation
        GNode<NodeAttribute,EdgeAttribute>* v1 = (*Gmean)[operations[3][ind_operation]->v1_i];
        GNode<NodeAttribute,EdgeAttribute>* v2 = (*G2)[operations[3][ind_operation]->v2_i];
        v1->attr = v2->attr;
        //incrementation cost and index of edit-operation
        total_edit_cost += operations[3][ind_operation]->cost;
        //std::cout << "       v1_i = "<< operations[3][ind_operation]->v1_i << " , v1_j = "<< operations[3][ind_operation]->v1_j << " , v2_i = "<< operations[3][ind_operation]->v2_i << " , v2_j = "<< operations[3][ind_operation]->v2_j << std::endl;
        //std::cout << "       current cost  = " << total_edit_cost << " , objective = " << alpha << std::endl;
        ind_operation ++;
    }
    if (ind_operation<operations[3].size()){
        // necessary decision whether the current operation must be applied or not
        GNode<NodeAttribute,EdgeAttribute>* v1 = (*Gmean)[operations[3][ind_operation]->v1_i];
        GNode<NodeAttribute,EdgeAttribute>* v2 = (*G2)[operations[3][ind_operation]->v2_i];
        v1->attr = ml->WeightedVertexMeanLabel(v1->attr,v2->attr,alpha - total_edit_cost);
        //std::cout << "       v1_i = "<< operations[3][ind_operation]->v1_i << " , v1_j = "<< operations[3][ind_operation]->v1_j << " , v2_i = "<< operations[3][ind_operation]->v2_i << " , v2_j = "<< operations[3][ind_operation]->v2_j << std::endl;
        //std::cout << "       final operation" << alpha << std::endl;
        return Gmean;
    }
    //std::cout << "ok 5" << std::endl;
    //-----------//
    //node insertion
    //-----------//
    ind_operation = 0;
    while (ind_operation<operations[4].size() && total_edit_cost + operations[4][ind_operation]->cost < alpha){
        // application operation
        int ind_inserted = operations[4][ind_operation]->v2_i;
        GNode<NodeAttribute,EdgeAttribute>* v2 = (*G2)[ind_inserted];
        int new_index_in_G1 = Gmean->Size();
        Gmean->Add(new GNode<NodeAttribute,EdgeAttribute> (new_index_in_G1, v2->attr));
        //update of edit operations with new indices, only for edge insertions
        for (int ind_op=0;ind_op<operations[5].size();ind_op++){
             if (operations[5][ind_op]->v2_i== ind_inserted) operations[5][ind_op]->v1_i= new_index_in_G1;
             if (operations[5][ind_op]->v2_j== ind_inserted) operations[5][ind_op]->v1_j= new_index_in_G1;
        }
        //incrementation cost and index of edit-operation
        total_edit_cost += operations[4][ind_operation]->cost;
        //std::cout << "       v1_i = "<< operations[4][ind_operation]->v1_i << " , v1_j = "<< operations[4][ind_operation]->v1_j << " , v2_i = "<< operations[4][ind_operation]->v2_i << " , v2_j = "<< operations[4][ind_operation]->v2_j << std::endl;
        //std::cout << "       current cost  = " << total_edit_cost << " , objective = " << alpha << std::endl;
        ind_operation ++;
    }
    if (ind_operation<operations[4].size()){
        // necessary decision whether the current operation must be applied or not
        if (alpha - total_edit_cost > operations[4][ind_operation]->cost/2){
            int ind_inserted = operations[4][ind_operation]->v2_i;
            GNode<NodeAttribute,EdgeAttribute>* v2 = (*G2)[ind_inserted];
            int new_index_in_G1 = Gmean->Size();
            Gmean->Add(new GNode<NodeAttribute,EdgeAttribute> (new_index_in_G1, v2->attr));
           // std::cout << "       v1_i = "<< operations[4][ind_operation]->v1_i << " , v1_j = "<< operations[4][ind_operation]->v1_j << " , v2_i = "<< operations[4][ind_operation]->v2_i << " , v2_j = "<< operations[4][ind_operation]->v2_j << std::endl;
           // std::cout << "       final operation" << std::endl;
        }
        return Gmean;
    }
    //std::cout << "ok 6" << std::endl;
    //-----------//
    //Edge insertion
    //-----------//
    ind_operation = 0;
    while (ind_operation<operations[5].size() && total_edit_cost + operations[5][ind_operation]->cost < alpha){
        // application operation
        GEdge<EdgeAttribute>* e2 = G2->getEdge(operations[5][ind_operation]->v2_i,operations[5][ind_operation]->v2_j);
        Gmean->Link(operations[5][ind_operation]->v1_i,operations[5][ind_operation]->v1_j,e2->attr);
        //incrementation cost and index of edit-operation
        total_edit_cost += operations[5][ind_operation]->cost;
        //std::cout << "       v1_i = "<< operations[5][ind_operation]->v1_i << " , v1_j = "<< operations[5][ind_operation]->v1_j << " , v2_i = "<< operations[5][ind_operation]->v2_i << " , v2_j = "<< operations[5][ind_operation]->v2_j << std::endl;
        //std::cout << "       current cost  = " << total_edit_cost << " , objective = " << alpha << std::endl;
        ind_operation ++;
    }
    if (ind_operation<operations[5].size()){
        // necessary decision whether the current operation must be applied or not
        if (alpha - total_edit_cost > operations[5][ind_operation]->cost/2){
             GEdge<EdgeAttribute>* e2 = G2->getEdge(operations[5][ind_operation]->v2_i,operations[5][ind_operation]->v2_j);
             Gmean->Link(operations[5][ind_operation]->v1_i,operations[5][ind_operation]->v1_j,e2->attr);
             //std::cout << "       v1_i = "<< operations[5][ind_operation]->v1_i << " , v1_j = "<< operations[5][ind_operation]->v1_j << " , v2_i = "<< operations[5][ind_operation]->v2_i << " , v2_j = "<< operations[5][ind_operation]->v2_j << std::endl;
             //std::cout << "       final operation" << std::endl;
        }
        return Gmean;
    }
    //std::cout << "ok 7" << std::endl;
    return Gmean;
};

template<class NodeAttribute, class EdgeAttribute, class PropertyType>
Graph<NodeAttribute, EdgeAttribute>* MedianGraph<NodeAttribute,EdgeAttribute,PropertyType>::ComputeWeightedMeanGraph(int indexGraph1, int indexGraph2, int * G1toG2, int * G2toG1, double alpha, GraphEditDistance<NodeAttribute, EdgeAttribute> * ed){
    Graph<NodeAttribute,EdgeAttribute> * G1 = this->ds->getGraph(indexGraph1);
    Graph<NodeAttribute,EdgeAttribute> * G2 = this->ds->getGraph(indexGraph2);
    return this->ComputeWeightedMeanGraph(G1,G2,G1toG2,G2toG1,alpha,ed);
};


#endif // __MEDIANGRAPH_H__

