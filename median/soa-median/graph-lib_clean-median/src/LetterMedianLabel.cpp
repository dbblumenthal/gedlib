#include "LetterMedianLabel.h"
#include <cmath>

double* computeGeometricMedian_bis(double* vectors, int lengthVector, int numVectors){
/*
    for (int i=0;i< lengthVector;i++){
        std::cout << "[" ;
        for (int j=0;j< numVectors;j++){
            std::cout << vectors[sub2ind(i,j,lengthVector)] << " ";
            }
        std::cout << "]" << std::endl;
    }
    */
    double SOD = 0;
    double old_SOD = 1000000.0;
    double* current_sol= new double[lengthVector];
    double* numerator= new double[lengthVector];
    double* denominator= new double[lengthVector];

    for (int i=0;i< lengthVector;i++){
        current_sol[i]=0;
        for (int j=0;j< numVectors;j++){
            current_sol[i]+=vectors[sub2ind(i,j,lengthVector)];
        }
        current_sol[i]/=numVectors;
        numerator[i]=current_sol[i];
        denominator[i]=current_sol[i];
    }
    double* norme = new double [numVectors];
    double delta =10;
    int num_it=0;
    while(delta > 0.00001){
    SOD = 0;
    //while(delta > 0.0001 && num_it<100){
        num_it ++;
        double* next_sol= new double[lengthVector];
        delta = 0;

        for (int j=0;j< numVectors;j++){
            norme[j]=0;
            for (int i=0;i< lengthVector;i++){
                norme[j] += (vectors[sub2ind(i,j,lengthVector)]-current_sol[i])* (vectors[sub2ind(i,j,lengthVector)]-current_sol[i]);
            }
            norme[j]=std::sqrt(norme[j]);
            SOD+= norme[j];
        }
        //std::cout << "SOD at iteration " << num_it << " = " << SOD<<  std::endl;

        for (int i=0;i< lengthVector;i++){
              numerator[i]=0;
             denominator[i]=0;
             for (int j=0;j< numVectors;j++){
                 if (norme[j] != 0) {
                    numerator[i]+= vectors[sub2ind(i,j,lengthVector)]/norme[j];
                    denominator[i]+= 1.0/norme[j];
                }
                // ??
                else{
                //numerator[i]+= 0.00001;
                //denominator[i]+= 0.00001;
                //std::cout << "!!!!! norme = 0 !!!!!" << std::endl;
                }
             }
             next_sol[i]=numerator[i]/denominator[i];
             //delta += std::abs(next_sol[i]-current_sol[i]);
        }
        delta = old_SOD - SOD;
        old_SOD = SOD;
        delete[] current_sol;
        current_sol = next_sol;
    }
    delete[] norme;
    delete[] numerator;
    delete[] denominator;
    return current_sol;
}

CMUPoint LetterMedianLabel::MedianNodeLabel(int * * mappingsFromMedian, int node1, Dataset<CMUPoint,double,int> * dataset){
//std::cout << " Computation of median label for vertex " <<  node1  << std::endl;
//std::cout << "\n --- calcul de label moyen pour sommet " << node1 << " du médian \n";
CMUPoint MedianLabel;
double sumX = 0;
double sumY = 0;
int num_assignments=0;
int N= dataset->size();
double * vectors = new double[N*2];

for (int i=0;i<N;i++){
  int mappedIndex = mappingsFromMedian[i][node1];
  if (mappedIndex < (*dataset)[i]-> Size()){
  	CMUPoint localPoint =  (*(*dataset)[i])[mappedIndex]->attr  ;
  	//std::cout << "label sommet " << node1 << " dans graphe " << i << " est " << label << std::endl;
  	sumX+= localPoint.x;
    sumY+= localPoint.y;
    vectors[sub2ind(0,num_assignments,2)]=localPoint.x;
    vectors[sub2ind(1,num_assignments,2)]=localPoint.y;
    num_assignments++;


    //std::cout << "------- sommet assigné dans graphe "<< i << " a les coordonnées ("<< localPoint.x << "," << localPoint.y <<")\n";
    }
}
if (num_assignments>0){
double* median_vector = new double[2];
median_vector = computeGeometricMedian_bis(vectors,2,num_assignments);
MedianLabel.x = median_vector[0];
MedianLabel.y = median_vector[1];
//MedianLabel.x = sumX/static_cast<double>(num_assignments);
//MedianLabel.y = sumY/static_cast<double>(num_assignments);
 //std::cout << "--- médian a les coordonnées ("<< MedianLabel.x << "," << MedianLabel.y <<")\n";
//std::cout << "label sommet " << node1 << " dans graphe  median est " << Medianlabel << std::endl;
}
else{
MedianLabel = (*(*dataset)[0])[0]->attr;
}
return MedianLabel;

};

double LetterMedianLabel::MedianEdgeLabel(Graph<CMUPoint,double> * Gbar, int * * mappingsFromMedian, int node1, int node2, Dataset<CMUPoint,double,int> * dataset){


    Graph<CMUPoint,double> *MedianGraph = new Graph<CMUPoint,double>(*Gbar);
    double absentEdge = std::numeric_limits<double>::max();
    int num_graphs_with_edge = 0;
    int N = dataset->size();
    for (int i=0;i<N;i++ ){
        //std::cout << "isEdge in graph " << i ;
        if ((mappingsFromMedian[i][node1]<dataset->getGraph(i)->Size())&&(mappingsFromMedian[i][node2]<dataset->getGraph(i)->Size())&&(dataset->getGraph(i)->getEdge(mappingsFromMedian[i][node1],mappingsFromMedian[i][node2]))){
        num_graphs_with_edge ++;
        //std::cout << " yes" << std::endl;
        }
        //std::cout << " no" << std::endl;
    }
    if(num_graphs_with_edge > N/2)
    return 1;
    else
    return absentEdge;
};

double LetterMedianLabel::NodeDelta(Graph<CMUPoint,double> * Gbar,int * * mappingsFromMedian, int nodeIndex, Dataset<CMUPoint,double,int> * dataset){

double best_Delta = 10000.0;
int N= dataset->size();
double actual_cost =0;
double potential_cost =0;
 //std::cout << "--- médian a les coordonnées ("<< MedianLabel.x << "," << MedianLabel.y <<")\n";
//std::cout << "label sommet " << node1 << " dans graphe  median est " << Medianlabel << std::endl;
for (int i=0;i<N;i++){
  int mappedIndex = mappingsFromMedian[i][nodeIndex];
  if (mappedIndex < (*dataset)[i]-> Size()){
  	actual_cost += this->cf->NodeSubstitutionCost((* Gbar)[nodeIndex],(*(*dataset)[i])[mappedIndex],Gbar,(*dataset)[i]);
  	potential_cost += this->cf->NodeInsertionCost((*(*dataset)[i])[mappedIndex],(*dataset)[i]);
  	for (int j=0;j<Gbar->Size();j++){
        GEdge<double> * EdgeInGbar=Gbar->getEdge(nodeIndex,j);
        if (EdgeInGbar){
            GEdge<double>* EdgeInDatasetGraph=(*dataset)[i]->getEdge(mappedIndex,mappingsFromMedian[i][j]);
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


return potential_cost - actual_cost;

};

double LetterMedianLabel::NodeLabelDelta(int median_size, int * * mappingsToMedian, int * * mappingsFromMedian, Graph<CMUPoint,double> * &Gbar, CMUPoint* NodeLabel,Dataset<CMUPoint,double,int> * dataset){
//return 100.0;

//Initialization
double deletion_cost = this->cf->NodeDeletionCost((*(dataset->getGraph(0)))[0],dataset->getGraph(0));
double insertion_cost = this->cf->NodeInsertionCost((*(dataset->getGraph(0)))[0],dataset->getGraph(0));
double best_Delta = 10000.0;
int ind_best_candidate;
int max_iter = 100;
int max_k_means = 3;
int num_candidates = max_k_means*(max_k_means +1)/2;
int  ** mappings_from_new_vertices = new int*[num_candidates];
int  * mappings_from_best_candidate = new int[dataset->size()];

//collecting all inserted labels of inserted vertices

std::vector<CMUPoint> inserted_labels;
double ** init_clusters_vectors = new double * [num_candidates];
std::vector<CMUPoint>* init_clusters_labels = new std::vector<CMUPoint>[num_candidates];
std::vector<CMUPoint>* updated_clusters_labels = new std::vector<CMUPoint>[num_candidates];
CMUPoint * candidate_labels = new CMUPoint[num_candidates];
bool * graph_has_insertions = new bool [dataset->size()];
std::vector<CMUPoint>* inserted_labels_by_graph = new std::vector<CMUPoint> [dataset->size()];
std::vector<int>* inserted_ids_by_graph = new std::vector<int> [dataset->size()];
for (int i = 0; i<dataset->size();i++){
    graph_has_insertions[i]=false;
    for (int j =0;j < dataset->getGraph(i)->Size();j++){
        if (mappingsToMedian[i][j]>= median_size){
            inserted_labels.push_back((*(dataset->getGraph(i)))[j]->attr);
            inserted_labels_by_graph[i].push_back((*(dataset->getGraph(i)))[j]->attr);
            inserted_ids_by_graph[i].push_back(j);
            graph_has_insertions[i]=true;
        }
    }
}
//k_means to generate candidate labels //initialization
for (int k_means = 1;k_means<=max_k_means;k_means++){
    std::random_shuffle(inserted_labels.begin(), inserted_labels.end());
    int i=0;
    int j=0;
    while (i+j<inserted_labels.size()){
            init_clusters_labels[(k_means*(k_means-1))/2+j].push_back(inserted_labels[i]);
            j++;
            if(j>=k_means) j=0;
            i++;
    }

}
for (int i=0;i<num_candidates;i++){
    init_clusters_vectors[i] = new double[init_clusters_labels[i].size()*2];
    for (int j = 0;j<init_clusters_labels[i].size();j++){
        //init_clusters_vectors[i][sub2ind(j,0,init_clusters_labels[i].size())]=init_clusters_labels[i][j].x;
        //init_clusters_vectors[i][sub2ind(j,1,init_clusters_labels[i].size())]=init_clusters_labels[i][j].y;
        init_clusters_vectors[i][sub2ind(0,j,2)]=init_clusters_labels[i][j].x;
        init_clusters_vectors[i][sub2ind(1,j,2)]=init_clusters_labels[i][j].y;
    }
    double * candidate_i_coordinates = new double[2];
    candidate_i_coordinates = computeGeometricMedian_bis(init_clusters_vectors[i],2,init_clusters_labels[i].size());
    candidate_labels[i].x=candidate_i_coordinates[0];
    candidate_labels[i].y=candidate_i_coordinates[1];
    delete[] candidate_i_coordinates;
}
for (int k_means = 2;k_means<=max_k_means;k_means++){
    bool convergence = false;
    int num_iter = 0;
    while (!convergence && num_iter < max_iter){
    double SOD = 0;
    double init_SOD = 0;
    //std::cout << "iteration " << num_iter << " of while loop" << std::endl;
        convergence = true;
         //update clusters
        for (int ind_cluster = (k_means*(k_means-1))/2;ind_cluster < (k_means*(k_means+1))/2;ind_cluster ++){
            for (int i = 0;i<init_clusters_labels[ind_cluster].size();i++){
                double best_sq_distance = -1.0;
                int best_index = -1;
                for (int ind_candidate = (k_means*(k_means-1))/2;ind_candidate < (k_means*(k_means+1))/2;ind_candidate ++){
                    double dx = init_clusters_labels[ind_cluster][i].x - candidate_labels[ind_candidate].x;
                    double dy = init_clusters_labels[ind_cluster][i].y - candidate_labels[ind_candidate].y;
                    double current_sq_distance = dx*dx + dy*dy;
                    if (ind_cluster == ind_candidate){
                        init_SOD+=std::sqrt(current_sq_distance);
                    }
                    if(best_sq_distance == -1.0 || current_sq_distance<best_sq_distance){
                        best_sq_distance=current_sq_distance;
                        best_index = ind_candidate;
                    }
                }
                SOD+=std::sqrt(best_sq_distance);
                if(best_index!=ind_cluster) {
                    convergence = false;
                    //std::cout << "best index = " << best_index << " , ind_cluster = " << ind_cluster << std::endl;
                }
                updated_clusters_labels[best_index].push_back(init_clusters_labels[ind_cluster][i]);
            }
        }
        /*
        std::cout << "clusters =" << std::endl;
        for (int ind_cluster = (k_means*(k_means-1))/2;ind_cluster < (k_means*(k_means+1))/2;ind_cluster ++){
            std::cout << "[ ";
            for(int i=0;i<updated_clusters_labels[ind_cluster].size();i++){
                std::cout << "("<< updated_clusters_labels[ind_cluster][i].x << ";" << updated_clusters_labels[ind_cluster][i].y <<  ") ";
            }
            std::cout << "]" << std::endl;
        */
        /*
        std::cout << "candidate Labels = ";
        for (int ind_cluster = (k_means*(k_means-1))/2;ind_cluster < (k_means*(k_means+1))/2;ind_cluster ++){
            std::cout << "("<< candidate_labels[ind_cluster].x << ";" << candidate_labels[ind_cluster].y <<  ") " ;
        }
        std::cout << "init_SOD before updating clusters = "<< init_SOD<< ", updated SOD = " << SOD <<  std::endl;
        */
        if (init_SOD < SOD){
        std::cout << "!!!!!! SOD going up on updating clusters !!!!!!!   init_SOD before updating clusters = "<< init_SOD<< ", updated SOD = " << SOD << std::endl;
        }
        init_SOD = SOD;
        SOD = 0;
        //update medians if not converged
        if (! convergence){
            for (int ind_cluster = (k_means*(k_means-1))/2;ind_cluster < (k_means*(k_means+1))/2;ind_cluster ++){
                init_clusters_labels[ind_cluster].clear();
                init_clusters_labels[ind_cluster]= updated_clusters_labels[ind_cluster];
                updated_clusters_labels[ind_cluster].clear();
                delete[] init_clusters_vectors[ind_cluster];
                init_clusters_vectors[ind_cluster] = new double[init_clusters_labels[ind_cluster].size()*2];
                for (int j = 0;j<init_clusters_labels[ind_cluster].size();j++){
                    //init_clusters_vectors[ind_cluster][sub2ind(j,0,init_clusters_labels[ind_cluster].size())]=init_clusters_labels[ind_cluster][j].x;
                    //init_clusters_vectors[ind_cluster][sub2ind(j,1,init_clusters_labels[ind_cluster].size())]=init_clusters_labels[ind_cluster][j].y;
                    init_clusters_vectors[ind_cluster][sub2ind(0,j,2)]=init_clusters_labels[ind_cluster][j].x;
                    init_clusters_vectors[ind_cluster][sub2ind(1,j,2)]=init_clusters_labels[ind_cluster][j].y;
                }
                double * candidate_i_coordinates = new double[2];
                candidate_i_coordinates = computeGeometricMedian_bis(init_clusters_vectors[ind_cluster],2,init_clusters_labels[ind_cluster].size());
                candidate_labels[ind_cluster].x=candidate_i_coordinates[0];
                candidate_labels[ind_cluster].y=candidate_i_coordinates[1];
                delete[] candidate_i_coordinates;
                double cluster_SOD = 0;
                for (int i = 0;i<init_clusters_labels[ind_cluster].size();i++){
                    double dx = init_clusters_labels[ind_cluster][i].x - candidate_labels[ind_cluster].x;
                    double dy = init_clusters_labels[ind_cluster][i].y - candidate_labels[ind_cluster].y;
                    double current_sq_distance = dx*dx + dy*dy;
                    SOD+= std::sqrt(current_sq_distance);
                }
            }
            /*
        if (init_SOD < SOD){
        std::cout << "!!!!!! SOD going up on updating centers !!!!!!!   init_SOD before updating clusters = "<< init_SOD<< ", updated SOD = " << SOD << std::endl;
        }
        */
        }
    num_iter ++;
    }
}
//cleaning memory

for (int i =0;i<num_candidates;i++){
    delete[] init_clusters_vectors[i];
}
delete[] init_clusters_vectors;
delete[] init_clusters_labels;

// candidate labels computed, run Block gradient descent for each label
for (int ind_candidate=0;ind_candidate<num_candidates;ind_candidate++){
    mappings_from_new_vertices[ind_candidate]=new int[dataset->size()];
    for (int i = 0;i < dataset->size();i++){
        mappings_from_new_vertices[ind_candidate][i] = dataset->getGraph(i)->Size();//initializing with only deletions
    }
    bool convergence =  false;
    int num_iter = 0;
    while(!convergence && num_iter < max_iter){
        convergence = true;
        //step 1: find optimal asssignement wrt label
        for (int graph_id = 0;graph_id<dataset->size();graph_id++){
            if(graph_has_insertions[graph_id]){
                int best_index = -1;
                double best_sq_distance;
                for(int i=0;i<inserted_ids_by_graph[graph_id].size();i++){
                    double dx = inserted_labels_by_graph[graph_id][i].x - candidate_labels[ind_candidate].x;
                    double dy = inserted_labels_by_graph[graph_id][i].y - candidate_labels[ind_candidate].y;
                    double current_sq_distance = (dx * dx) + (dy * dy);
                    if (best_index == -1 || current_sq_distance < best_sq_distance){
                        best_index = inserted_ids_by_graph[graph_id][i];
                        best_sq_distance = current_sq_distance;
                    }
                }
                if (std::sqrt(best_sq_distance) > deletion_cost)
                    best_index = dataset->getGraph(graph_id)->Size();
                if (best_index != mappings_from_new_vertices[ind_candidate][graph_id]){
                    mappings_from_new_vertices[ind_candidate][graph_id]=best_index;
                    convergence = false;
                }

            }
        }
        //step 2: recompute median if not converged
        if (!convergence){
            std::vector<CMUPoint> assigned_labels;
            for (int graph_id = 0;graph_id<dataset->size();graph_id++){
                if (mappings_from_new_vertices[ind_candidate][graph_id] < dataset->getGraph(graph_id)->Size()){
                    assigned_labels.push_back((*(dataset->getGraph(graph_id)))[mappings_from_new_vertices[ind_candidate][graph_id]]->attr);
                }
            }
            int num_subs = assigned_labels.size();
            double * vectors = new double [num_subs * 2];
            for (int i = 0;i<num_subs;i++){
                vectors[sub2ind(i,0,num_subs)]=assigned_labels[i].x;
                vectors[sub2ind(i,1,num_subs)]=assigned_labels[i].y;
            }
            double * candidate_coordinates = new double[2];
            candidate_coordinates = computeGeometricMedian_bis(vectors,num_subs,2);
            candidate_labels[ind_candidate].x=candidate_coordinates[0];
            candidate_labels[ind_candidate].y=candidate_coordinates[1];
            delete[] candidate_coordinates;
            delete[] vectors;
        }
        num_iter++;
    }
    // optimal assignments and candidate computed, evaluating Delta Cost
    double current_Delta = 0;
    for (int i = 0;i < dataset->size();i++){
        if (mappings_from_new_vertices[ind_candidate][i] >= dataset->getGraph(i)->Size()){
            current_Delta += deletion_cost;
        }
        else{
            CMUPoint substituted_label = (*(dataset->getGraph(i)))[mappings_from_new_vertices[ind_candidate][i]]->attr;
            current_Delta += this->cf->SubstitutionCost(substituted_label,candidate_labels[ind_candidate]) - insertion_cost;
        }
    }
    if(current_Delta<best_Delta){
        best_Delta = current_Delta;
        ind_best_candidate = ind_candidate;
        for (int i = 0;i<dataset->size();i++){
            mappings_from_best_candidate[i]=mappings_from_new_vertices[ind_candidate][i];
        }
    }
    delete[] mappings_from_new_vertices[ind_candidate];
}
//Best candidate and mappings identified, updating AllPairMappings and Graph if Delta <0
if (best_Delta < 0){
    //std::cout << "improving insertion found !!" << std::endl;
    GNode<CMUPoint,double> * Newnode = new GNode<CMUPoint,double>(median_size,candidate_labels[ind_best_candidate]);
    Gbar->Add(Newnode);
    int M = dataset->size();
    for(int i=0;i<M;i++){
        int size_graph_i = dataset->getGraph(i)->Size();
        //updating mappingsToMedian so that vertices mapped to nMedian (supposed to be deleted) are still deleted w.r.t. a bigger median graph
        //and assigning the proper nodes to the new vertex.
        int best_assigned_vertex = size_graph_i+1;
        for (int j=0;j<size_graph_i;j++){
            if (mappingsToMedian[i][j] == median_size)
                mappingsToMedian[i][j] = median_size+1;
            if (j== mappings_from_best_candidate[i])
                mappingsToMedian[i][j] = median_size;
        }
        //reallocating memory for new mappingsFromMedian in order to take in account the new size of Gmed.

        int * UpdatedMappingsFromMedian_to_i = new int [median_size+1];
        for(int j=0;j<median_size;j++){
            UpdatedMappingsFromMedian_to_i[j]=mappingsFromMedian[i][j];
        }
        //newly added vertex is assigned to propoer vertex in graph i
        UpdatedMappingsFromMedian_to_i[median_size]=mappings_from_best_candidate[i];
        //std::cout << "pre last update of mappings to graph " << i << std::endl;
        //PrintPointerContent(mappingsToMedian[i],this->ds->getGraph(i)->Size());
        delete[] mappingsFromMedian[i];
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
*NodeLabel= candidate_labels[ind_best_candidate];
return best_Delta;

}

double LetterMedianLabel::WeightedEdgeMeanLabel(double label1, double label2, double alpha) { return 1; }


CMUPoint LetterMedianLabel::WeightedVertexMeanLabel(CMUPoint label1, CMUPoint label2, double alpha)
{
  alpha /= cf->SubstitutionCost(label1,label2);
  double a1 = 1.0 - alpha;
  CMUPoint nodeLabel;
  nodeLabel.x = alpha * label2.x + a1 * label1.x;
  nodeLabel.y = alpha * label2.y + a1 * label1.y;
  return nodeLabel;
}
