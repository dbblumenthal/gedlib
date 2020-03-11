/*
 * @file test_graph.cpp
 * @author Évariste <<evariste.daller@unicaen.fr>>
 *
 *
 * Calcule toutes les distances entre toutes les molécules d'un datatset
 *
 */


#include <unistd.h>
#include <iostream>
#include <string>
#include <sstream>

#include <limits>

#include "graph.h"
#include "tinyxml.h"
#include "Dataset.h"
#include "GraphEditDistance.h"
#include "SymbolicGraph.h"
#include "ConstantGraphEditDistance.h"
#include "MedianGraph.h"
#include "ConstantMedianLabel.h"
#include "CMUGraph.h"
#include "CMUDataset.h"
#include "CMUMedianLabel.h"
#include "CMUCostFunction.h"
#include "LetterGraph.h"
#include "LetterDataset.h"
#include "LetterMedianLabel.h"
#include "LetterCostFunction.h"
#include "BipartiteGraphEditDistance.h"
#include "BipartiteGraphEditDistanceMulti.h"
#include "GreedyGraphEditDistance.h"
#include "RandomWalksGraphEditDistance.h"
#include "RandomWalksGraphEditDistanceMulti.h"
#include "IPFPGraphEditDistance.h"
#include "RandomMappings.h"
#include "MultistartRefinementGraphEditDistance.h"
#include "GNCCPGraphEditDistance.h"
#include "utils.h"
using namespace std;



void usage (char * s)
{
  cerr << "Usage : "<< s << " dataset " << " options"<<endl;
  cerr << "options:" << endl;
  cerr << "\t -m method " << endl;
  cerr << "\t \t Specify the algorithm used to compute edit distance" << endl;
  cerr << "\t -o output_file " << endl;
  cerr << "\t \t Specify a filename for outputing gram matrix" << endl;
  cerr << "\t -c cns,cni,cnd,ces,cei,ced " << endl;
  cerr << "\t \t Specify edit operation costs" << endl;
  cerr << "\t -p n_edit_paths " << endl;
  cerr << "\t \t Specify the number of edit paths to compute GED (lsape_multi)" << endl;
}

struct Options{
  string dataset_file = "";
  string train_file ="";
  string method = "";
  string initial_method = "";
  string output_file = "";
  double cns = 1;
  double cni = 3;
  double cnd = 3;
  double ces = 1;
  double cei = 3;
  double ced = 3;
  bool all_train_set_sizes = false;
  bool shuffle = false;
  bool cmu = false;
  bool letter = false;
  int k = 3;
  int nep = 100; // number of edit paths for allsolution
  int nrep = 100; // default number of returned solutions from the launch of nep parallel IPFPs
  bool diag = false; // if true, instances are tested against themselves only -> opt is always 0.
  double adap = 0; // adaptative factor, algorithm gets more exploratory when it starts to converge. 0 = not adaptative (up to v15 option -a does not exist and is considered 0) 1 = strongly adaptative
  int numshuffle = 0;
  bool median = 0;
};

struct Options * parseOptions(int argc, char** argv){
  struct Options * options = new struct Options();
  options->train_file = string(argv[1]);
  options->dataset_file = string(argv[2]);
  int opt;
  stringstream sstream;
  while ((opt = getopt(argc, argv, "m:i:o:c:sp:r:zyda:h:et")) != -1) {
    switch (opt) {
    case 'm':
      options->method = string(optarg);
      break;
     case 'i':
      options->initial_method = string(optarg);
      break;
    case 'o':
      options->output_file = string(optarg);
      break;
    case 'c':
      sstream << optarg;
      sstream >> options->cns;
      sstream >> options->cni;
      sstream >> options->cnd;
      sstream >> options->ces;
      sstream >> options->cei;
      sstream >> options->ced;
      break;
    case 's':
      options->shuffle=true;
      break;
     case 't':
      options->all_train_set_sizes=true;
      break;
    case 'p':
      sstream << optarg;
      sstream >> options->nep;
      break;
    case 'z':
      options->cmu = true;
      break;
     case 'y':
      options->letter=true;
      break;
    case 'd':
      options->diag = true;
      break;
    case 'r':
	options->nrep = atoi(optarg);
      //sstream << optarg;
      //sstream >> options->nrep;
      break;
     case 'a':
	options->adap = atoi(optarg);
      //sstream << optarg;
      //sstream >> options->nrep;
      break;
      case 'h':
	options->numshuffle = atoi(optarg);
      //sstream << optarg;
      //sstream >> options->nrep;
      break;
      case 'e':
        options->median = true;
      break;

    default: /* '?' */
      cerr << "Options parsing failed."  << endl;
      usage(argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  if (options->nrep > options->nep) options->nrep=options->nep;
  if (options->initial_method == "") options->initial_method=options->method;
  return options;
}

template <class NodeAttribute, class EdgeAttribute>
void printGraph(const Graph< NodeAttribute, EdgeAttribute>  &g) {
 for (int i=0; i<g.Size(); i++){
      std::cout << i << " " << g[i]->attr << std::endl;
    }

    for (int i=0; i<g.Size(); i++){
      GEdge<EdgeAttribute> *p = g[i]->getIncidentEdges();
      while(p){
        std::cout << i << " " << p->IncidentNode() << " " <<  p->attr << std::endl;
        p = p->Next();
      }
    }
}

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
double * computeGraphEditDistance(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
				  GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
				  bool shuffle, int nep=0, int numshuffle=0){
  if(shuffle)
    dataset->shuffleize();



  //double * distances =  dataset->computeGraphEditDistance(ed, true);
  int N = dataset->size();
  int numtimes = N*N;
double sumtimes = 0;
  double* distances = new double[N*N];
  struct timeval  tv1, tv2;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      #ifdef PRINT_TIMES
        gettimeofday(&tv1, NULL);
      #endif
      //std::cout << "before copy\n\n";
      //printGraph<NodeAttribute, EdgeAttribute>((*(*dataset)[j]));
      Graph<NodeAttribute, EdgeAttribute>  *Ginit1 = dataset->getGraph(i);  // pointer on the graph i
      Graph<NodeAttribute, EdgeAttribute> *Gperm1 = new Graph<NodeAttribute, EdgeAttribute>(*Ginit1);  // copy to a new graph
      Graph<NodeAttribute, EdgeAttribute>  *Ginit2 = dataset->getGraph(j);  // pointer on the graph j
      Graph<NodeAttribute, EdgeAttribute> *Gperm2 = new Graph<NodeAttribute, EdgeAttribute>(*Ginit2);  // copy to a new graph

      //Gperm->shuffleize();
      //Graph<NodeAttribute, EdgeAttribute>  G(*dataset->getGraph(j));
      //Graph<NodeAttribute, EdgeAttribute>  G2(*dataset->getGraph(i));
     // std::cout << "before perm\n\n";
      //printGraph<NodeAttribute, EdgeAttribute>(*Gperm2);

      for (int h=0; h<numshuffle; h++) {
      	        Gperm1->shuffleize();
                Gperm2->shuffleize();
               // std::cout << "shuffle " << h << std::endl << std::endl;
               // printGraph<NodeAttribute, EdgeAttribute >(*Gperm2);
	}

   //  distances[sub2ind(i,j,N)] = (*ed)((*dataset)[i], Gperm);
      distances[sub2ind(i,j,N)] = (*ed)(Gperm1, Gperm2);
  //    distances[sub2ind(i,j,N)] = (*ed)((*dataset)[i], (*dataset)[j]); //commenter cette ligne et décommenter les trois lignes d'avant pour un vrai shufflize...

      cout << (*dataset)[i]->Size() << " ";
      cout << (*dataset)[j]->Size() << " ";
      #ifdef PRINT_TIMES
        gettimeofday(&tv2, NULL);
        cout << ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec)) << " " ;
        sumtimes+= ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
      #endif

      cout << (int)distances[sub2ind(i,j,N)];
      delete Gperm1;
      delete Gperm2;
      cout << endl;
    }
  }
   std::cout << "temps moyen = " << (sumtimes/numtimes) << std::endl;
  return distances;
}




template <class NodeAttribute, class EdgeAttribute, class PropertyType>
Graph<NodeAttribute, EdgeAttribute> * computeMedianGraph(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
				  GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
				  MedianLabel< NodeAttribute, EdgeAttribute, PropertyType> * mcf, int SetMedianIndex, int * * mappingsFromMedian, int * * mappingsToMedian, double * distancesToMedian){

//std::cout <<"\n ----------\n starting generalized computation" << endl;

  double SOD =0;

  int N = dataset->size();


  Graph<NodeAttribute, EdgeAttribute>  *Gmed = dataset->getGraph(SetMedianIndex);  // pointer on the graph of index bestSODIndex
  Graph<NodeAttribute, EdgeAttribute> *medianG = new Graph<NodeAttribute, EdgeAttribute>(*Gmed);  // copy to a new graph
  MedianGraph<NodeAttribute, EdgeAttribute,PropertyType>* mg = new MedianGraph<NodeAttribute,EdgeAttribute,PropertyType>(medianG,dataset,mcf);

 int nMedian = medianG->Size();


// update vertex-label and edge-set of the median graph
 bool isModifiedGraph = mg->updateMedianGraph( mappingsFromMedian); // isModified is false when updating the median does not modifiy it.

 bool ismodifiedEP = true;

 while (isModifiedGraph || ismodifiedEP ){
   ismodifiedEP = false;
   // updating the SOD
   SOD = 0;
   for (int i=0;i<N;i++){
        Graph<NodeAttribute, EdgeAttribute>  *gi = dataset->getGraph(i);
      	int m = gi->Size();
        distancesToMedian[i] = ed->GedFromMapping(medianG,gi,mappingsFromMedian[i],nMedian,mappingsToMedian[i],m);
 	SOD+= distancesToMedian[i];
	 }


   //std::cout <<endl<< "SOD = " << SOD << endl;
   // updating the edit paths
   for (int i=0;i<N;i++){
        Graph<NodeAttribute, EdgeAttribute>  *gi = dataset->getGraph(i);
        int m = gi->Size();
        int * localMappingToMedian = new int[m];
        int * localMappingFromMedian = new int[nMedian];
        ed->getOptimalMapping(medianG,gi,localMappingFromMedian,localMappingToMedian);
        double localDistance = ed->GedFromMapping(medianG,gi,localMappingFromMedian,nMedian,localMappingToMedian,m);
        if (localDistance< distancesToMedian[i]){
            delete mappingsFromMedian[i];
            delete mappingsToMedian[i];
            mappingsFromMedian[i] = localMappingFromMedian;
            mappingsToMedian[i] = localMappingToMedian;
            ismodifiedEP = true;
            }
        else {
            delete [] localMappingFromMedian;
            delete [] localMappingToMedian;
           }


        }



   //updating the median graph
   isModifiedGraph = mg->updateMedianGraph( mappingsFromMedian);
   //if  (ismodifiedEP) std::cout <<endl<< " EP modified " << endl;
   //if  (isModifiedGraph) std::cout << " Median Graph modified " << endl;


 }
   return medianG;


}

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
int computeSetMedianGraph(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
				  GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
				  int * * mappingsFromMedian, int * * mappingsToMedian, double * distancesToMedian){
//std::cout << "starting set_median computation" << endl;

  double SOD;
  double bestSOD;
  int bestSODIndex = 0;

  int N = dataset->size();
  double* distances = new double[N*N];
   int * * AllPairMappings = new int*[N*N];

  struct timeval  tv1, tv2;


  for (int i=0; i<N-1; i++){
     distances[sub2ind(i,i,N)] = 0;
    SOD = 0;
    for (int j=i+1; j<N; j++){
      Graph<NodeAttribute, EdgeAttribute>  *g1 = dataset->getGraph(i);  // pointer on the graph i
      Graph<NodeAttribute, EdgeAttribute>  *g2 = dataset->getGraph(j);  // pointer on the graph j
      int n=g1->Size();
      int m=g2->Size();
      AllPairMappings[sub2ind(i,j,N)] = new int[n];
      AllPairMappings[sub2ind(j,i,N)] = new int[m];
      ed->getOptimalMapping(g1,g2,AllPairMappings[sub2ind(i,j,N)],AllPairMappings[sub2ind(j,i,N)]);
      distances[sub2ind(i,j,N)] = ed->GedFromMapping(g1,g2,AllPairMappings[sub2ind(i,j,N)],n,AllPairMappings[sub2ind(j,i,N)],m);
      //distances[sub2ind(i,j,N)] = (*ed)((*dataset)[i], (*dataset)[j]); //commenter cette ligne et décommenter les trois lignes d'avant pour un vrai shufflize...
      distances[sub2ind(j,i,N)] = distances[sub2ind(i,j,N)];


      //cout << (int)distances[sub2ind(i,j,N)];
      //cout << endl;
    }
    for (int j=0; j<N; j++) SOD+=distances[sub2ind(i,j,N)];
    if (i==0 || SOD < bestSOD){
       bestSOD = SOD;
       bestSODIndex = i;
    }


  }
  Graph<NodeAttribute, EdgeAttribute>  *Gmed = dataset->getGraph(bestSODIndex);  // pointer on the graph of index bestSODIndex
  //std::cout <<endl<< "index of set median = " << bestSODIndex << ", SOD = " << bestSOD << endl;


// creating and adding the self mapping into AllpairMappings only for set-median graph
 int nMedian = Gmed->Size();
 int * selfAssignementTo = new int[nMedian];
 int * selfAssignementFrom = new int[nMedian];
 for (int i=0; i<nMedian; i++) {
    selfAssignementTo[i]=i;
    selfAssignementFrom[i]=i;
    }
 for (int i=0; i<N; i++) {
    if (i==bestSODIndex){
        mappingsFromMedian[i]=selfAssignementFrom;
        mappingsToMedian[i]=selfAssignementTo;
        distancesToMedian[i]=0.0;
    }
    else{
        mappingsFromMedian[i]=AllPairMappings[sub2ind(bestSODIndex,i,N)];
        mappingsToMedian[i]=AllPairMappings[sub2ind(i,bestSODIndex,N)];
        distancesToMedian[i]=distances[sub2ind(i,bestSODIndex,N)];
        }
 }

 return bestSODIndex;

}


template <class NodeAttribute, class EdgeAttribute, class PropertyType>
double * computeGraphEditDistanceDiag(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
				  GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
				  bool shuffle, int nep=0){
  if(shuffle)
    dataset->shuffleize();

  //double * distances =  dataset->computeGraphEditDistance(ed, true);
  int N = dataset->size();
  double* distances = new double[N*N];
  struct timeval  tv1, tv2;
  for (int i=0; i<N; i++){

      #ifdef PRINT_TIMES
        gettimeofday(&tv1, NULL);
      #endif

      Graph<NodeAttribute, EdgeAttribute>  G(*dataset->getGraph(i));
      G.shuffleize();
      distances[sub2ind(i,i,N)] = (*ed)((*dataset)[i], &G);

      cout << (*dataset)[i]->Size() << " ";
      cout << (*dataset)[i]->Size() << " ";
      #ifdef PRINT_TIMES
        gettimeofday(&tv2, NULL);
        cout << ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec)) << " " ;
      #endif

      cout << (int)distances[sub2ind(i,i,N)];
      cout << endl;


  }
  return distances;
}
template <class NodeAttribute, class EdgeAttribute, class PropertyType>
double run_classification_by_1nn(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * representative_set,
                  Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
				  GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
				  int nep =0){

  struct timeval  tv1, tv2;
 #ifdef PRINT_TIMES
        gettimeofday(&tv1, NULL);
      #endif


std::vector<int> num_correct;
std::vector<int> num_total;
std::vector<PropertyType> properties;
int M = dataset->size();
//std::cout<< "\n---------\n dataset size = " << M << std::endl;
int N = representative_set->size();
//std::cout<< " representative set size = " << N << std::endl;
double best_distance = 0;
int closest_neighbor_index = 0;
for (int i=0; i< M; i++){
    //find the closest neighbor
    for (int j=0; j<N; j++){
         int n = (*representative_set)[j]->Size();
         int m = (*dataset)[i]->Size();
         int * localMappingToRepresentative = new int[m];
         int * localMappingFromRepresentative = new int[n];
         ed->getOptimalMapping((*representative_set)[j],(*dataset)[i],localMappingFromRepresentative,localMappingToRepresentative);
         double local_distance = ed->GedFromMapping((*representative_set)[j],(*dataset)[i],localMappingFromRepresentative,n,localMappingToRepresentative,m);
         if (j==0 || local_distance < best_distance){
            best_distance= local_distance;
            closest_neighbor_index = j;
            }
         delete [] localMappingFromRepresentative;
         delete [] localMappingToRepresentative;
        }
        //std::cout<< "best distance between representative graph and datagraph is " << best_distance << std::endl;
    // check if classification is correct
    //std::cout<< "representative Property = " << (*representative_set)(closest_neighbor_index) << " , Graph Property = " << (*dataset)(i) << std::endl;
    int p=0;
    while (p<=properties.size()){
        if (p<properties.size() && (*dataset)(i)==properties[p]){
            num_total[p]++;
            if  ((*dataset)(i)==(*representative_set)(closest_neighbor_index)){
                num_correct[p]++;
            }
        break;
        }
        p++;
    }
    if (p==properties.size()+1){
        num_total.push_back(1);
        num_correct.push_back(0);
        properties.push_back((*dataset)(i));
        if  ((*dataset)(i)==(*representative_set)(closest_neighbor_index)){
                num_correct[p]++;
            }

        }
    }
    int num_total_final=0;
    int num_correct_final=0;
    for (int p=0;p<properties.size();p++){
        //std::cout << "\n - results for class " <<  properties[p] << " - \n" << "number of graphs in the class = " <<  num_total[p] << " , number of correct classifications = " << num_correct[p] ;
        //std::cout << "\ntaux de precision pour la classe = " <<  100*static_cast<double>(num_correct[p])/static_cast<double>(num_total[p]) << "%";
        num_total_final += num_total[p];
        num_correct_final += num_correct[p];
    }


double taux_precision = static_cast<double>(num_correct_final)/static_cast<double>(num_total_final);
//std::cout << "\n\n - taux de precision total = " << taux_precision*100 << "%"<<std::endl;
 #ifdef PRINT_TIMES
        gettimeofday(&tv2, NULL);
        cout << " - computing time = "<<((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec)) << " seconds. \n \n " ;
      #endif


return taux_precision;
}

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
double run_classification_by_1nn_Letter(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * representative_set,
                  Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
				  GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
				  int nep =0){
struct timeval  tv1, tv2;
 #ifdef PRINT_TIMES
        gettimeofday(&tv1, NULL);
      #endif
std::vector<int> num_correct;
std::vector<int> num_total;
std::vector<PropertyType> properties;
int M = dataset->size();
//std::cout<< "\n---------\n dataset size = " << M << std::endl;
int N = representative_set->size();
//std::cout<< " representative set size = " << N << std::endl;
double best_distance = 0;
int closest_neighbor_index = 0;
for (int i=0; i< M; i++){
    //find the closest neighbor
    for (int j=0; j<N; j++){
         int n = (*representative_set)[j]->Size();
         int m = (*dataset)[i]->Size();
         int * localMappingToRepresentative = new int[m];
         int * localMappingFromRepresentative = new int[n];
         ed->getOptimalMapping((*representative_set)[j],(*dataset)[i],localMappingFromRepresentative,localMappingToRepresentative);
         double local_distance = ed->GedFromMapping((*representative_set)[j],(*dataset)[i],localMappingFromRepresentative,n,localMappingToRepresentative,m);
         if (j==0 || local_distance < best_distance){
            best_distance= local_distance;
            closest_neighbor_index = j;
            }
         delete [] localMappingFromRepresentative;
         delete [] localMappingToRepresentative;
        }
        //std::cout<< "best distance between representative graph and datagraph is " << best_distance << std::endl;
    // check if classification is correct
    //std::cout<< "representative Property = " << (*representative_set)(closest_neighbor_index) << " , Graph Property = " << (*dataset)(i) << std::endl;
    int p=0;
    while (p<=properties.size()){
        if (p<properties.size() && (*dataset)(i)==properties[p]){
            num_total[p]++;
            if  ((*dataset)(i)==(*representative_set)(closest_neighbor_index)){
                num_correct[p]++;
            }
        break;
        }
        p++;
    }
    if (p==properties.size()+1){
        num_total.push_back(1);
        num_correct.push_back(0);
        properties.push_back((*dataset)(i));
        if  ((*dataset)(i)==(*representative_set)(closest_neighbor_index)){
                num_correct[p-1]++;
            }

        }
    }
    int num_total_final=0;
    int num_correct_final=0;
    for (int p=0;p<properties.size();p++){
        //std::cout << "\n - results for class " <<  (char)properties[p] << " - \n" << "number of graphs in the class = " <<  num_total[p] << " , number of correct classifications = " << num_correct[p] ;
        //std::cout << "\ntaux de precision pour la classe = " <<  100*static_cast<double>(num_correct[p])/static_cast<double>(num_total[p]) << "%";
        num_total_final += num_total[p];
        num_correct_final += num_correct[p];
    }

double taux_precision = static_cast<double>(num_correct_final)/static_cast<double>(num_total_final);
//std::cout << "\n\n - taux de precision total = " << taux_precision*100 << "%"<<std::endl;
std::cout << ";"<< taux_precision*100;
#ifdef PRINT_TIMES
        gettimeofday(&tv2, NULL);
        cout << ";"<< ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
      #endif

return taux_precision;
}






int main (int argc, char** argv)
{

  struct Options * options =   parseOptions(argc,argv);


  options->k = 3;







 if(options->letter){
  LetterDistanceCost *cf = new LetterDistanceCost(0.9,1.7,0.75);
  LetterMedianLabel * mcf =  new LetterMedianLabel(cf);


  // IPFP used as a refinement method
  IPFPGraphEditDistance<CMUPoint,double> * algoIPFP = new IPFPGraphEditDistance<CMUPoint,double>(cf);
  algoIPFP->recenterInit();

  // Sinkhorn balanced random init
  if (options->method == string("ipfpe_random_sh")){
    algoIPFP->continuousRandomInit(true);
    options->method = string("ipfpe_multi_random");
  }

  GraphEditDistance<CMUPoint,double>* ed;
  if(options->method == string("lsape_bunke"))
    ed = new BipartiteGraphEditDistance<CMUPoint,double>(cf);
  else if( options->method == string("lsape_multi_bunke") )
    ed = new BipartiteGraphEditDistanceMulti<CMUPoint,double>(cf, options->nep);
  else if( options->method == string("lsape_multi_greedy") )
    ed = new GreedyGraphEditDistance<CMUPoint,double>(cf, options->nep);
  //else if( options->method == string("multi_random") )
  //  ed = new RandomInitForIPFP<CMUPoint,double>(cf, options->nep);

  else if(options->method == string("ipfpe_flat")){
    algoIPFP->continuousFlatInit(true);
    ed = algoIPFP->clone();
  }

  else if(options->method == string("ipfpe_bunke")){
    BipartiteGraphEditDistance<CMUPoint,double> *ed_init = new BipartiteGraphEditDistance<CMUPoint,double>(cf);
    ed =new IPFPGraphEditDistance<CMUPoint,double>(cf,ed_init);
  } else if(options->method == string("ipfpe_multi_bunke")){
    BipartiteGraphEditDistanceMulti<CMUPoint,double> *ed_init = new BipartiteGraphEditDistanceMulti<CMUPoint,double>(cf, options->nrep);
    ed = new MultistartRefinementGraphEditDistance<CMUPoint,double>(cf, ed_init, options->nep, algoIPFP,options->nrep,options->adap);
  } else if(options->method == string("ipfpe_multi_random")){
    RandomMappingsGED<CMUPoint,double> *init = new RandomMappingsGED<CMUPoint,double>();
    ed = new MultistartRefinementGraphEditDistance<CMUPoint,double>(cf, init, options->nep, algoIPFP,options->nrep,options->adap);

  } else if(options->method == string("ipfpe_multi_greedy")){
    GreedyGraphEditDistance<CMUPoint,double> *ed_init = new GreedyGraphEditDistance<CMUPoint,double>(cf, options->nep);
    ed = new MultistartRefinementGraphEditDistance<CMUPoint,double>(cf, ed_init, options->nep, algoIPFP);

  } else if(options->method == string("gnccp")){
    //RandomWalksGraphEditDistance *ed_init = new RandomWalksGraphEditDistance(cf,3 );
    ed = new GNCCPGraphEditDistance<CMUPoint,double>(cf);//,ed_init);

  } else{
    cerr << "Undefined graph edit distance algorithm "<< endl;
    usage(argv[0]);
    return EXIT_FAILURE;
  }
     GraphEditDistance<CMUPoint,double>* i_ed;
  if(options->initial_method == string("lsape_bunke"))
    i_ed = new BipartiteGraphEditDistance<CMUPoint,double>(cf);
  else if( options->initial_method == string("lsape_multi_bunke") )
    i_ed = new BipartiteGraphEditDistanceMulti<CMUPoint,double>(cf, options->nep);
  else if( options->initial_method == string("lsape_multi_greedy") )
    i_ed = new GreedyGraphEditDistance<CMUPoint,double>(cf, options->nep);
  //else if( options->method == string("multi_random") )
  //  ed = new RandomInitForIPFP<CMUPoint,double>(cf, options->nep);

  else if(options->method == string("ipfpe_flat")){
    algoIPFP->continuousFlatInit(true);
    i_ed = algoIPFP->clone();
  }

  else if(options->initial_method == string("ipfpe_bunke")){
    BipartiteGraphEditDistance<CMUPoint,double> *ed_init = new BipartiteGraphEditDistance<CMUPoint,double>(cf);
    i_ed =new IPFPGraphEditDistance<CMUPoint,double>(cf,ed_init);
  } else if(options->initial_method == string("ipfpe_multi_bunke")){
    BipartiteGraphEditDistanceMulti<CMUPoint,double> *ed_init = new BipartiteGraphEditDistanceMulti<CMUPoint,double>(cf, options->nrep);
    i_ed = new MultistartRefinementGraphEditDistance<CMUPoint,double>(cf, ed_init, options->nep, algoIPFP,options->nrep,options->adap);
  } else if(options->initial_method == string("ipfpe_multi_random")){
    RandomMappingsGED<CMUPoint,double> *init = new RandomMappingsGED<CMUPoint,double>();
    i_ed = new MultistartRefinementGraphEditDistance<CMUPoint,double>(cf, init, options->nep, algoIPFP,options->nrep,options->adap);

  } else if(options->initial_method == string("ipfpe_multi_greedy")){
    GreedyGraphEditDistance<CMUPoint,double> *ed_init = new GreedyGraphEditDistance<CMUPoint,double>(cf, options->nep);
    i_ed = new MultistartRefinementGraphEditDistance<CMUPoint,double>(cf, ed_init, options->nep, algoIPFP);

  } else if(options->initial_method == string("gnccp")){
    //RandomWalksGraphEditDistance *ed_init = new RandomWalksGraphEditDistance(cf,3 );
    i_ed = new GNCCPGraphEditDistance<CMUPoint,double>(cf);//,ed_init);

  } else{
    cerr << "Undefined graph edit distance algorithm "<< endl;
    usage(argv[0]);
    return EXIT_FAILURE;
  }
  int dataset_reduction =0;
  int dataset_max_reduction=std::numeric_limits<int>::max();

  while(dataset_reduction < dataset_max_reduction){
      LetterDataset * dataset = new LetterDataset(options->dataset_file.c_str());
      LetterDataset * trainset = new LetterDataset(options->train_file.c_str());
      bool next_train_test = true;
        struct timeval  tv1, tv2;
        #ifdef PRINT_TIMES
        //gettimeofday(&tv1, NULL);
        #endif

      //Separation des différentes classes en différents datasets

       //initialisation
       int M = trainset->size();
       std::vector<LetterDataset*> datasets_by_class;
       LetterDataset* initial_dataset = new LetterDataset();
       initial_dataset->add((*trainset)[0],(*trainset)(0));
       datasets_by_class.push_back(initial_dataset);
       //separation
       for (int i=1;i<M;i++){
            int j=0;
            while (j<datasets_by_class.size() && (*trainset)(i)!= (*datasets_by_class[j])(0)){
                j++;
                }
            if (j<datasets_by_class.size()){
                datasets_by_class[j]->add((*trainset)[i],(*trainset)(i));
                }
            else{
                LetterDataset* new_dataset = new LetterDataset();
                new_dataset->add((*trainset)[i],(*trainset)(i));
                datasets_by_class.push_back(new_dataset);
                }
            }
        //std::cout << M ;
        // Medians computation as a new dataset

        LetterDataset * medians_dataset = new LetterDataset();
        LetterDataset * Smedians_dataset = new LetterDataset();

        int P = datasets_by_class.size();

         for(int i = 0; i<P;i++){
            if (dataset_max_reduction > datasets_by_class[i]->size())
            dataset_max_reduction = datasets_by_class[i]->size();
            }
        //std::cout << P << " classes computed \n";
        // for each class
        //train_set_reduction
          for(int i = 0; i<P;i++){
                for (int j=0;j<dataset_reduction;j++){
                datasets_by_class[i]->pop_back();
                }
            }

          LetterDataset * reduced_trainset = new LetterDataset();

          for(int i = 0; i<P;i++){
                for (int j=0;j<datasets_by_class[i]->size();j++){
                    reduced_trainset->add((*datasets_by_class[i])[j],(*datasets_by_class[i])(j));
                }
            }
            std::cout<<reduced_trainset->size();

        for(int i = 0; i<P;i++){
            //computation of set median
            int N = datasets_by_class[i]->size();
            double * distancesToMedian = new double [N];
            int * * mappingsFromMedian = new int*[N];
            int * * mappingsToMedian = new int*[N];
            int  SetMedianIndex = computeSetMedianGraph(datasets_by_class[i],i_ed,mappingsFromMedian,mappingsToMedian, distancesToMedian);
            Smedians_dataset->add((*datasets_by_class[i])[SetMedianIndex],(*datasets_by_class[i])(SetMedianIndex));

            //computation of generalized median
            Graph<CMUPoint,double>* local_median = computeMedianGraph(datasets_by_class[i],ed,mcf,SetMedianIndex,mappingsFromMedian,mappingsToMedian, distancesToMedian);
            //std::cout << "median computed for class with Property " << (char)(*datasets_by_class[i])(0) << " that has " << N << " graphs.\n";
            medians_dataset->add(local_median,(*datasets_by_class[i])(0));
            //std::cout << "median added to the dataset \n";
            //cleaning the memory
            delete [] distancesToMedian;
            for (int j = 0;j<N;j++){
                delete [] mappingsFromMedian[j];
                delete [] mappingsToMedian[j];
                }
            delete [] mappingsFromMedian;
            delete [] mappingsToMedian;

            }
            //std::cout << "set of medians computed , dataset size = " << dataset->size();
            #ifdef PRINT_TIMES
            //gettimeofday(&tv2, NULL);
            //cout << ";"<< ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
            #endif
          //std::cout<< "--- 1nn classification with generalized medians ---\n";
          run_classification_by_1nn_Letter(medians_dataset,dataset,ed);
          //std::cout<< "--- 1nn classification with set medians --- \n";
          run_classification_by_1nn_Letter(Smedians_dataset,dataset,ed);
          //std::cout<< "--- 1nn classification with train set --- \n";
          run_classification_by_1nn_Letter(reduced_trainset,dataset,ed);

      delete trainset;
      delete medians_dataset;
      delete dataset;

      if (options->all_train_set_sizes){
        dataset_reduction++;
        }
      else{
        dataset_reduction=dataset_max_reduction;
        }
        std::cout <<"\n";

  }

  delete cf;
  delete mcf;
  delete options;
  delete algoIPFP;
  delete i_ed;
  delete ed;
  return 0;

 }

else{
 ConstantEditDistanceCost* cf = new ConstantEditDistanceCost(options->cns,options->cni, options->cnd,
							       options->ces,options->cei, options->ced);
  ConstantMedianLabel* mcf = new ConstantMedianLabel(cf);


  // IPFP used as a refinement method
  IPFPGraphEditDistance<int,int> * algoIPFP = new IPFPGraphEditDistance<int,int>(cf);
  algoIPFP->recenterInit();

  // Sinkhorn balanced random init
  if (options->method == string("ipfpe_random_sh")){
    algoIPFP->continuousRandomInit(true);
    options->method = string("ipfpe_multi_random");
  }

  GraphEditDistance<int,int>* ed;
  if(options->method == string("lsape_bunke"))
    ed = new BipartiteGraphEditDistance<int,int>(cf);
  else if( options->method == string("lsape_multi_bunke") )
    ed = new BipartiteGraphEditDistanceMulti<int,int>(cf, options->nep);
  else if( options->method == string("lsape_rw"))
    ed = new RandomWalksGraphEditDistance(cf,options->k);
  else if( options->method == string("lsape_multi_rw") )
    ed = new RandomWalksGraphEditDistanceMulti(cf, options->k, options->nep);
  else if( options->method == string("lsape_multi_greedy") )
    ed = new GreedyGraphEditDistance<int,int>(cf, options->nep);
  //else if( options->method == string("multi_random") )
  //  ed = new RandomInitForIPFP<int,int>(cf, options->nep);

  else if(options->method == string("ipfpe_flat")){
    algoIPFP->continuousFlatInit(true);
    ed = algoIPFP->clone();
  }

  else if(options->method == string("ipfpe_bunke")){
    BipartiteGraphEditDistance<int,int> *ed_init = new BipartiteGraphEditDistance<int,int>(cf);
    ed =new IPFPGraphEditDistance<int,int>(cf,ed_init);
  } else if(options->method == string("ipfpe_multi_bunke")){
    BipartiteGraphEditDistanceMulti<int,int> *ed_init = new BipartiteGraphEditDistanceMulti<int,int>(cf, options->nrep);
    ed = new MultistartRefinementGraphEditDistance<int,int>(cf, ed_init, options->nep, algoIPFP,options->nrep,options->adap);
  } else if(options->method == string("ipfpe_multi_rw")){
    RandomWalksGraphEditDistanceMulti *ed_init = new RandomWalksGraphEditDistanceMulti(cf, options->k, options->nep);
    ed = new MultistartRefinementGraphEditDistance<int,int>(cf, ed_init, options->nep, algoIPFP);
  } else if(options->method == string("ipfpe_rw")){
    RandomWalksGraphEditDistance *ed_init = new RandomWalksGraphEditDistance(cf,options->k);
    ed =new IPFPGraphEditDistance<int,int>(cf,ed_init);

  } else if(options->method == string("ipfpe_multi_random")){
    RandomMappingsGED<int,int> *init = new RandomMappingsGED<int,int>();
    ed = new MultistartRefinementGraphEditDistance<int,int>(cf, init, options->nep, algoIPFP,options->nrep,options->adap);

  } else if(options->method == string("ipfpe_multi_greedy")){
    GreedyGraphEditDistance<int,int> *ed_init = new GreedyGraphEditDistance<int,int>(cf, options->nep);
    ed = new MultistartRefinementGraphEditDistance<int,int>(cf, ed_init, options->nep, algoIPFP);

  } else if(options->method == string("gnccp")){
    //RandomWalksGraphEditDistance *ed_init = new RandomWalksGraphEditDistance(cf,3 );
    ed = new GNCCPGraphEditDistance<int,int>(cf);//,ed_init);

  } else{
    cerr << "Undefined graph edit distance algorithm "<< endl;
    usage(argv[0]);
    return EXIT_FAILURE;
  }

   GraphEditDistance<int,int>* i_ed;
  if(options->initial_method == string("lsape_bunke"))
    i_ed = new BipartiteGraphEditDistance<int,int>(cf);
  else if( options->initial_method == string("lsape_multi_bunke") )
    i_ed = new BipartiteGraphEditDistanceMulti<int,int>(cf, options->nep);
  else if( options->initial_method == string("lsape_multi_greedy") )
    i_ed = new GreedyGraphEditDistance<int,int>(cf, options->nep);
  //else if( options->method == string("multi_random") )
  //  ed = new RandomInitForIPFP<CMUPoint,double>(cf, options->nep);

  else if(options->method == string("ipfpe_flat")){
    algoIPFP->continuousFlatInit(true);
    i_ed = algoIPFP->clone();
  }

  else if(options->initial_method == string("ipfpe_bunke")){
    BipartiteGraphEditDistance<int,int> *ed_init = new BipartiteGraphEditDistance<int,int>(cf);
    i_ed =new IPFPGraphEditDistance<int,int>(cf,ed_init);
  } else if(options->initial_method == string("ipfpe_multi_bunke")){
    BipartiteGraphEditDistanceMulti<int,int> *ed_init = new BipartiteGraphEditDistanceMulti<int,int>(cf, options->nrep);
    i_ed = new MultistartRefinementGraphEditDistance<int,int>(cf, ed_init, options->nep, algoIPFP,options->nrep,options->adap);
  } else if(options->initial_method == string("ipfpe_multi_random")){
    RandomMappingsGED<int,int> *init = new RandomMappingsGED<int,int>();
    i_ed = new MultistartRefinementGraphEditDistance<int,int>(cf, init, options->nep, algoIPFP,options->nrep,options->adap);

  } else if(options->initial_method == string("ipfpe_multi_greedy")){
    GreedyGraphEditDistance<int,int> *ed_init = new GreedyGraphEditDistance<int,int>(cf, options->nep);
    i_ed = new MultistartRefinementGraphEditDistance<int,int>(cf, ed_init, options->nep, algoIPFP);

  } else if(options->initial_method == string("gnccp")){
    //RandomWalksGraphEditDistance *ed_init = new RandomWalksGraphEditDistance(cf,3 );
    i_ed = new GNCCPGraphEditDistance<int,int>(cf);//,ed_init);

  } else{
    cerr << "Undefined graph edit distance algorithm "<< endl;
    usage(argv[0]);
    return EXIT_FAILURE;
  }
 int dataset_reduction =0;
  int dataset_max_reduction=std::numeric_limits<int>::max();

  while(dataset_reduction < dataset_max_reduction){
      ChemicalDataset<double> * dataset = new ChemicalDataset<double>(options->dataset_file.c_str());
      ChemicalDataset<double> * trainset = new ChemicalDataset<double>(options->train_file.c_str());
      bool next_train_test = true;
        struct timeval  tv1, tv2;
        #ifdef PRINT_TIMES
        //gettimeofday(&tv1, NULL);
        #endif

      //Separation des différentes classes en différents datasets

       //initialisation
       int M = trainset->size();
       std::vector<ChemicalDataset<double>*> datasets_by_class;
       ChemicalDataset<double>* initial_dataset = new ChemicalDataset<double>();
       initial_dataset->add((*trainset)[0],(*trainset)(0));
       datasets_by_class.push_back(initial_dataset);
       //separation
       for (int i=1;i<M;i++){
            int j=0;
            while (j<datasets_by_class.size() && (*trainset)(i)!= (*datasets_by_class[j])(0)){
                j++;
                }
            if (j<datasets_by_class.size()){
                datasets_by_class[j]->add((*trainset)[i],(*trainset)(i));
                }
            else{
                ChemicalDataset<double>* new_dataset = new ChemicalDataset<double>();
                new_dataset->add((*trainset)[i],(*trainset)(i));
                datasets_by_class.push_back(new_dataset);
                }
            }
        //std::cout << M ;
        // Medians computation as a new dataset

        ChemicalDataset<double> * medians_dataset = new ChemicalDataset<double>();
        ChemicalDataset<double> * Smedians_dataset = new ChemicalDataset<double>();

        int P = datasets_by_class.size();

         for(int i = 0; i<P;i++){
            if (dataset_max_reduction > datasets_by_class[i]->size())
            dataset_max_reduction = datasets_by_class[i]->size();
            }
        //std::cout << P << " classes computed \n";
        // for each class
        //train_set_reduction
          for(int i = 0; i<P;i++){
                for (int j=0;j<dataset_reduction;j++){
                datasets_by_class[i]->pop_back();
                }
            }

          ChemicalDataset<double> * reduced_trainset = new ChemicalDataset<double>();

          for(int i = 0; i<P;i++){
                for (int j=0;j<datasets_by_class[i]->size();j++){
                    reduced_trainset->add((*datasets_by_class[i])[j],(*datasets_by_class[i])(j));
                }
            }
            std::cout<<reduced_trainset->size();

        for(int i = 0; i<P;i++){
            //computation of set median
            int N = datasets_by_class[i]->size();
            double * distancesToMedian = new double [N];
            int * * mappingsFromMedian = new int*[N];
            int * * mappingsToMedian = new int*[N];
            int  SetMedianIndex = computeSetMedianGraph(datasets_by_class[i],i_ed,mappingsFromMedian,mappingsToMedian, distancesToMedian);
            Smedians_dataset->add((*datasets_by_class[i])[SetMedianIndex],(*datasets_by_class[i])(SetMedianIndex));

            //computation of generalized median
            Graph<int,int>* local_median = computeMedianGraph(datasets_by_class[i],ed,mcf,SetMedianIndex,mappingsFromMedian,mappingsToMedian, distancesToMedian);
            //std::cout << "median computed for class with Property " << (char)(*datasets_by_class[i])(0) << " that has " << N << " graphs.\n";
            medians_dataset->add(local_median,(*datasets_by_class[i])(0));
            //std::cout << "median added to the dataset \n";
            //cleaning the memory
            delete [] distancesToMedian;
            for (int j = 0;j<N;j++){
                delete [] mappingsFromMedian[j];
                delete [] mappingsToMedian[j];
                }
            delete [] mappingsFromMedian;
            delete [] mappingsToMedian;

            }
            //std::cout << "set of medians computed , dataset size = " << dataset->size();
            #ifdef PRINT_TIMES
            //gettimeofday(&tv2, NULL);
            //cout << ";"<< ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
            #endif
          //std::cout<< "--- 1nn classification with generalized medians ---\n";
          run_classification_by_1nn_Letter(medians_dataset,dataset,ed);
          //std::cout<< "--- 1nn classification with set medians --- \n";
          run_classification_by_1nn_Letter(Smedians_dataset,dataset,ed);
          //std::cout<< "--- 1nn classification with train set --- \n";
          run_classification_by_1nn_Letter(reduced_trainset,dataset,ed);

      delete trainset;
      delete medians_dataset;
      delete dataset;

      if (options->all_train_set_sizes){
        dataset_reduction++;
        }
      else{
        dataset_reduction=dataset_max_reduction;
        }
        std::cout <<"\n";

  }

  delete cf;
  delete mcf;
  delete options;
  delete algoIPFP;
  delete i_ed;
  delete ed;
  return 0;

 }


}
