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
#include <cstdlib>
#include <cmath>

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
#include "WebGraph.h"
#include "WebDataset.h"
#include "WebMedianLabel.h"
#include "WebCostFunction.h"
#include "LetterGraph.h"
#include "LetterDataset.h"
#include "LetterMedianLabel.h"
#include "LetterCostFunction.h"
#include "IBDGraph.h"
#include "IBDDataset.h"
#include "IBDMedianLabel.h"
#include "IBDCostFunction.h"
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
	double cns = 2;
	double cni = 4;
	double cnd = 4;
	double ces = 1;
	double cei = 1;
	double ced = 1;
	bool all_train_set_sizes = false;
	double size_random_trainset = 0;
	double proportion_random_trainset = 100;
	double time_limit = 600;
	bool shuffle = false;
	bool cmu = false;
	bool letter = false;
	unsigned int power_cost = 1;
	bool web = false;
	int k = 3;
	int nep = 100; // number of edit paths for allsolution
	int nrep = 100; // default number of returned solutions from the launch of nep parallel IPFPs
	bool diag = false; // if true, instances are tested against themselves only -> opt is always 0.
	double adap = 0; // adaptative factor, algorithm gets more exploratory when it starts to converge. 0 = not adaptative (up to v15 option -a does not exist and is considered 0) 1 = strongly adaptative
	int numshuffle = 0;
	bool median = 0;
	int num_repetitions_xp = 1;
	int proportion_random_init_median = 100;
	int dataset_seed = 123;
	string collection_id ="0";
	string collection_percentage="10";
	string ibd_costs = "";
};

struct Options * parseOptions(int argc, char* argv[]){
	struct Options * options = new struct Options();
	options->dataset_file = string(argv[1]);
	int opt;
	stringstream sstream;
	optind = 2;
	while ((opt = getopt(argc, argv, "m:i:o:c:sp:r:zywdaq:l:g:h:etu:v:x:b:f:I:P:")) != -1) {
		switch (opt) {
		case 'I':
			options->collection_id = string(optarg);
			break;
		case 'P':
			options->collection_percentage = string(optarg);
			break;
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
		case 'g':
			options->power_cost = (unsigned int)atoi(optarg);
			break;
		case 'w':
			options->web=true;
			break;
		case 'd':
			options->diag = true;
			break;
		case 'l':
			options->time_limit = atof(optarg);
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
		case 'q':
			options->ibd_costs = std::string(optarg);
			break;
		case 'u':
			options->size_random_trainset = atoi(optarg);
			//sstream << optarg;
			//sstream >> options->nrep;
			break;
		case 'x':
			options->num_repetitions_xp = atoi(optarg);
			//sstream << optarg;
			//sstream >> options->nrep;
			break;
		case 'v':
			options->proportion_random_trainset = atoi(optarg);
			//sstream << optarg;
			//sstream >> options->nrep;
			break;
		case 'h':
			options->numshuffle = atoi(optarg);
			//sstream << optarg;
			//sstream >> options->nrep;
			break;
		case 'f':
			options->dataset_seed = atoi(optarg);
			//sstream << optarg;
			//sstream >> options->nrep;
			break;
		case 'b':
			options->proportion_random_init_median = atoi(optarg);
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
	if (options->train_file=="") options->train_file = options->dataset_file;
	return options;
}


void printGraph(const Graph< int, int>  &g) {
	for (int i=0; i<g.Size(); i++){
		std::cout << i << " " << g[i]->attr << std::endl;
	}

	for (int i=0; i<g.Size(); i++){
		GEdge<int> *p = g[i]->getIncidentEdges();
		while(p){
			std::cout << i << " " << p->IncidentNode() << " " <<  p->attr << std::endl;
			p = p->Next();
		}
	}
}
void printLetterGraph(Graph<CMUPoint,double> * set_median){
	std::cout << "Graph " << std::endl;
	std::cout << "  --  vertices -- " << std::endl;
	for (int j =0;j<set_median->Size();j++){
		std::cout<<"[" << j << ",(" << (*set_median)[j]->attr.x << " , " << (*set_median)[j]->attr.y << " )],";
	}
	std::cout<<std::endl;
	std::cout << "  --  edges -- " << std::endl;
	for (int j =0;j<set_median->Size();j++){
		for (int k = j+1; k <set_median->Size();k++){
			if (set_median->getEdge(j, k))
				std::cout << "(" << j << "," << k << "),";
		}
	}
	std::cout<<std::endl;
}

double* computeGeometricMedian(double* vectors, int lengthVector, int numVectors){
	/*
    for (int i=0;i< lengthVector;i++){
        std::cout << "[" ;
        for (int j=0;j< numVectors;j++){
            std::cout << vectors[sub2ind(i,j,lengthVector)] << " ";
            }
        std::cout << "]" << std::endl;
    }
	 */
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
		}

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
					numerator[i]+= 0.0001;
					denominator[i]+= 0.0001;
				}
			}
			next_sol[i]=numerator[i]/denominator[i];
			delta += std::abs(next_sol[i]-current_sol[i]);
		}
		delete[] current_sol;
		current_sol = next_sol;
	}
	delete[] norme;
	delete[] numerator;
	delete[] denominator;
	return current_sol;
}

double computeAlpha_and_ProjectedVector(double* p, double* A, double* B, int n, double * projection){
	/*
    std::cout<<"Alpha computation :" << std::endl;
    std::cout<<"A = [";
        for (int j = 0; j < n;j++){
        std::cout << A[j] << " ";
        }
        std::cout<< "]"<<std::endl;
    std::cout<<"B = [";
        for (int j = 0; j < n;j++){
        std::cout << B[j] << " ";
        }
        std::cout<< "]"<<std::endl;
	 */
	double alpha;
	double* AB = new double[n];
	double sq_dist=0;
	for (int i=0;i<n;i++){
		AB[i]=B[i]-A[i];
		sq_dist+= AB[i]*AB[i];
		//std::cout << "A[i] = " << A[i] <<" , B[i] = " <<B[i]  <<" , sq_dist = " << sq_dist << std::endl;
	}
	if (sq_dist == 0){
		//projection = A;
		for (int i=0;i<n;i++){
			projection[i]=A[i];
		}
		delete[] AB;
		return 0;
	}
	else{
		double* Ap = new double[n];
		for (int i=0;i<n;i++){
			Ap[i]=p[i]-A[i];
		}
		double dot_Ap_AB = 0;
		for (int i=0;i<n;i++){
			dot_Ap_AB+= Ap[i]*AB[i];
		}
		alpha = dot_Ap_AB / sq_dist;
		if (alpha <= 0) {
			//projection = A;
			for (int i=0;i<n;i++){
				projection[i]=A[i];
			}
			delete[] AB; delete[] Ap;
			return 0;
		}
		else if (alpha >= 1){
			//projection = B;
			for (int i=0;i<n;i++){
				projection[i]=B[i];
			}
			delete[] AB; delete[] Ap;
			return 1;
		}
		else{
			//double* q = new double[n];
			for (int i=0;i<n;i++){
				projection[i]=A[i] + alpha * AB[i];
			}

			delete[] AB; delete[] Ap;
			return alpha;
		}

	}

}

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
double computeSOD(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		Graph<NodeAttribute,EdgeAttribute> * medianG, double time_limit=600.0
)
{
	struct timeval  tv1, tv2;
	gettimeofday(&tv1, NULL);
	int N= dataset->size();
	double SOD = 0;
	for (int i =0;i<N;i++){
		SOD+=(*ed)(dataset->getGraph(i), medianG);
		gettimeofday(&tv2, NULL);
		double running_time= (double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec) ;
		if (running_time > time_limit){
			return -1;
		}
	}
	return SOD;
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
bool computeMappingsAndDistances(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		int** AllPairMappings, double * distances, double time_limit = 1800.00){
	//double * distances =  dataset->computeGraphEditDistance(ed, true);
	int N = dataset->size();
	struct timeval  tv1, tv2;
	gettimeofday(&tv1, NULL);

	for (int i=0; i<N; i++){
		Graph<NodeAttribute,EdgeAttribute>  *g1 = dataset->getGraph(i);
		int n=g1->Size();

		for (int j=i+1; j<N; j++){
			// pointer on the graph i
			Graph<NodeAttribute,EdgeAttribute>  *g2 = dataset->getGraph(j);  // pointer on the graph j
			int m=g2->Size();
			AllPairMappings[sub2ind(i,j,N)] = new int[n];
			AllPairMappings[sub2ind(j,i,N)] = new int[m];
		}
	}

	for (int i=0; i<N; i++){
		Graph<NodeAttribute,EdgeAttribute>  *g1 = dataset->getGraph(i);
		int n=g1->Size();
		distances[sub2ind(i,i,N)] = 0;
		int * selfAssignement = new int[n];
		for (int j=0; j<n; j++) {
			selfAssignement[j]=j;
		}
		AllPairMappings[sub2ind(i,i,N)]=selfAssignement;

		for (int j=i+1; j<N; j++){
			// pointer on the graph i
			Graph<NodeAttribute,EdgeAttribute>  *g2 = dataset->getGraph(j);  // pointer on the graph j
			int m=g2->Size();
			ed->getOptimalMapping(g1,g2,AllPairMappings[sub2ind(i,j,N)],AllPairMappings[sub2ind(j,i,N)]);
			distances[sub2ind(i,j,N)] = ed->GedFromMapping(g1,g2,AllPairMappings[sub2ind(i,j,N)],n,AllPairMappings[sub2ind(j,i,N)],m);
			//distances[sub2ind(i,j,N)] = (*ed)((*dataset)[i], (*dataset)[j]); //commenter cette ligne et décommenter les trois lignes d'avant pour un vrai shufflize...
			distances[sub2ind(j,i,N)] = distances[sub2ind(i,j,N)];
			//cout << (int)distances[sub2ind(i,j,N)];
			//cout << endl;
			gettimeofday(&tv2, NULL);
			double time=((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
			if (time > time_limit){
				return false;
			}
		}
	}
	/*
     for (int i=0; i<N; i++){
      std::cout<< "[" << std::endl;
        for (int j=0; j<N; j++) {
            std::cout << distances[sub2ind(i,j,N)] << " ";
        }
     std::cout<< "]" << std::endl;
     }
	 */
	return true;

}



template <class NodeAttribute, class EdgeAttribute, class PropertyType>
void printDatasetStatistics(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset){
	//std::cout << "number of graphs = "<<dataset->size() << std::endl;
	double sumnodes = 0;
	double sumedges =0;
	for (int i = 0; i<dataset->size();i++){
		sumnodes += (*dataset)[i]->Size();
		sumedges += (*dataset)[i]->getNbEdges();
	}
	//std::cout << "mean_number of nodes = " <<  sumnodes / dataset->size()<< std::endl;
	// std::cout << "mean_number of edges = " <<  sumedges / dataset->size()<< std::endl;
	std::cout <<  sumnodes / dataset->size()<< ";"<<sumedges / dataset->size()<<";";
	// std::cout << "mean_number of edges = " <<  sumedges / dataset->size()<< std::endl;


}

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
void printGraphStatistics(Graph< NodeAttribute, EdgeAttribute> * graph){
	//std::cout << "number of graphs = "<<dataset->size() << std::endl;


	//std::cout << "mean_number of nodes = " <<  sumnodes / dataset->size()<< std::endl;
	// std::cout << "mean_number of edges = " <<  sumedges / dataset->size()<< std::endl;
	std::cout <<  graph->Size()<< ";"<<graph->getNbEdges()<<";";
	// std::cout << "mean_number of edges = " <<  sumedges / dataset->size()<< std::endl;


}
template <class NodeAttribute, class EdgeAttribute, class PropertyType>
Graph<NodeAttribute, EdgeAttribute> * computeRecursiveWeightedMins(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		MedianLabel< NodeAttribute, EdgeAttribute, PropertyType> * mcf, int max_recursion)
		{      int N = dataset->size();
		double* distances =new double[N*N];
		int * * AllPairMappings = new int*[N*N];
		for (int i=0; i<N; i++){
			Graph<int,int>  *g1 = dataset->getGraph(i);
			int n=g1->Size();
			distances[sub2ind(i,i,N)] = 0;
			int * selfAssignement = new int[n];
			for (int j=0; j<n; j++) {
				selfAssignement[j]=j;
			}
			AllPairMappings[sub2ind(i,i,N)]=selfAssignement;
			for (int j=i+1; j<N; j++){
				// pointer on the graph i
				Graph<int,int>  *g2 = dataset->getGraph(j);  // pointer on the graph j
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
		}
		Graph<NodeAttribute, EdgeAttribute> * medianG = new Graph<NodeAttribute, EdgeAttribute> ;
		medianG=computeRecursiveWeightedMins(dataset,ed,mcf,max_recursion, AllPairMappings,distances);
		return medianG;
		}

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
Graph<NodeAttribute, EdgeAttribute> * computeSODBasedRecursiveWeightedMins(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		MedianLabel< NodeAttribute, EdgeAttribute, PropertyType> * mcf, int max_recursion, int ** AllPairMappings, double * distances, bool refine = 0,bool keep_best =0)
		{      int N = dataset->size();

		double* SODs = new double[N];
		int * sorted_indices = new int[max_recursion+1];
		for (int i=0;i<max_recursion+1;i++){
			sorted_indices[i]=-1;
		}
		for (int i = 0;i<N;i++){
			//compute distance between graph i and median in vector space
			/*
        std::cout<<"[";
        for (int j = 0; j < max_recursion+1;j++){
        std::cout << sorted_indices[j] << " ";
        }
        std::cout<< "]"<<std::endl;
			 */
			SODs[i] = 0;
			for (int j = 0;j<N;j++){
				SODs[i]+= distances[sub2ind(i,j,N)];
			}
			//std::cout << "Squared distance from Graph " << i << " = " << vec_distances[i] << std::endl;
			//std::cout << "distance from Graph " << i << " = " << vec_distances[i] << std::endl;
			//checks if new_graph belongs in the max_recursion+1 closest to median, and updates the sorted list.
			for (int j = 0; j < max_recursion+1;j++){
				if (sorted_indices[j] == -1){
					sorted_indices[j]=i;
					j=max_recursion+1;
				}
				else if (SODs[i] < SODs[sorted_indices[j]]){
					for (int k= max_recursion;k>j;k--){
						sorted_indices[k]=sorted_indices[k-1];
					}
					sorted_indices[j]=i;
					j=max_recursion+1;
				}
			}
		}
		std::cout<<"graphs with minimal SODs = [";
		for (int j = 0; j < max_recursion+1;j++){
			std::cout << sorted_indices[j] << " ";
		}
		std::cout<< "]"<<std::endl;
		//initialization

		double SOD = 0;
		double* distancesFromMedian = new double [N];
		for (int i = 0;i<N;i++){
			SOD+=distances[sub2ind(sorted_indices[0],i,N)];
			distancesFromMedian[i]=distances[sub2ind(sorted_indices[0],i,N)];
		}
		//std::cout<<"initial SOD of Graph " << sorted_indices[0] << " = " << SOD << std::endl;
		std::cout<<"initial SOD of Graph " << sorted_indices[0] << " = " << SOD << std::endl;
		double * reduced_distances = new double[(max_recursion+1)*(max_recursion+1)];

		// Copying relevant mappings to new pointers
		int ** mappingsFromMedian = new int*[N];
		int ** mappingsToMedian = new int*[N];
		bool self_assigned = false;
		for (int i = 0;i<N;i++){
			//std::cout << "copying mapping from " << sorted_indices[0] << " to " << i << std::endl;
			mappingsFromMedian[i] = new int[dataset->getGraph(sorted_indices[0])->Size()];
			for (int j = 0;j<dataset->getGraph(sorted_indices[0])->Size();j++){
				mappingsFromMedian[i][j]= AllPairMappings[sub2ind(sorted_indices[0],i,N)][j];
			}

			//std::cout << "copying mapping from " <<  i << " to " <<  sorted_indices[0]<< std::endl;
			mappingsToMedian[i] = new int[dataset->getGraph(i)->Size()];
			for (int j = 0; j<dataset->getGraph(i)->Size();j++){
				mappingsToMedian[i][j]= AllPairMappings[sub2ind(i,sorted_indices[0],N)][j];
			}
		}
		//std::cout << "copying distances"<<std::endl;
		for (int i = 0;i<max_recursion+1;i++){
			//std::cout<<"[";
			for (int j = 0; j<max_recursion+1;j++){
				reduced_distances[sub2ind(i,j,max_recursion+1)]= distances[sub2ind(sorted_indices[i],sorted_indices[j],N)];
				//std::cout << " " << reduced_distances[sub2ind(i,j,max_recursion+1)];
			}
			//std::cout<<"]"<<std::endl;
		}

		double* median_vector = computeGeometricMedian(reduced_distances,max_recursion+1,max_recursion+1);



		double* current_projection = new double[max_recursion+1];
		double* current_vector = new double[max_recursion+1];
		double * new_projection=new double[max_recursion+1];

		for (int i = 0;i<max_recursion+1;i++){
			current_projection[i]=reduced_distances[sub2ind(0,i,max_recursion+1)];
		}
		//Graph<int,int> * medianG = new Graph<int,int>(*(test_set->getGraph(sorted_indices[0])));
		Graph<NodeAttribute, EdgeAttribute> * medianG = new Graph<NodeAttribute, EdgeAttribute>(*(dataset->getGraph(sorted_indices[0])));
		Graph<NodeAttribute, EdgeAttribute> * mediantmp = NULL;
		MedianGraph<NodeAttribute, EdgeAttribute, PropertyType>* mg = new MedianGraph<NodeAttribute, EdgeAttribute, PropertyType>(mediantmp,dataset,mcf);
		double alpha=0;
		//recursion
		for (int i = 0; i< max_recursion;i++){
			// compute new projected vector and corresponding alpha
			for (int j = 0;j<max_recursion+1;j++){
				current_vector[j]=reduced_distances[sub2ind(i+1,j,max_recursion+1)];
			}
			alpha = computeAlpha_and_ProjectedVector(median_vector,current_projection,current_vector,max_recursion+1,new_projection);
			//alpha = 0;

			//std::cout << "alpha between current projection and graph " << sorted_indices[i+1] << " = " << alpha << std::endl;
			//create new weighted mean Graph
			//medianG=mg->ComputeWeightedMeanGraph(sorted_indices[i],sorted_indices[i+1],mappingsFromMedian[sorted_indices[i+1]],mappingsToMedian[sorted_indices[i+1]],alpha*distancesFromMedian[sorted_indices[i+1]],ed);
			mediantmp=mg->ComputeWeightedMeanGraph(medianG,(*dataset)[sorted_indices[i+1]],mappingsFromMedian[sorted_indices[i+1]],mappingsToMedian[sorted_indices[i+1]],alpha*distancesFromMedian[sorted_indices[i+1]],ed);
			//delete medianG;

			//compute assignments and SOD

			double SOD_temp=0;
			int ** mapfrom_temp = new int*[N];
			int ** mapto_temp = new int*[N];
			for (int j=0; j<N; j++){
				Graph<NodeAttribute, EdgeAttribute>  *g2 = dataset->getGraph(j);  // pointer on the graph j
				int n=mediantmp->Size();
				int m=g2->Size();
				mapfrom_temp[j] = new int[n];
				mapto_temp[j] = new int[m];
				//ed->getOptimalMapping(medianG,g2,mappingsFromMedian[j],mappingsToMedian[j]);
				//std::cout << "ok A graph" << j << " , median Size = " << n << " , current graph Size = " << m << std::endl;
				//printGraph(*medianG);
				ed->getOptimalMapping(mediantmp,g2,mapfrom_temp[j],mapto_temp[j]);
				//std::cout << "ok B" << std::endl;
				//distancesFromMedian[j] = ed->GedFromMapping(medianG,g2,mappingsFromMedian[j],n,mappingsToMedian[j],m);
				distancesFromMedian[j] = ed->GedFromMapping(mediantmp,g2,mapfrom_temp[j],n,mapto_temp[j],m);

				//std::cout << "ok C" << std::endl;
				SOD_temp+= distancesFromMedian[j] ;
			}
			std::cout << "new SOD at recursion "<< i << " = " << SOD_temp << " ";
			if (SOD_temp < SOD || !keep_best){
				for (int j=0; j<N; j++){
					delete[] mappingsFromMedian[j];
					delete[] mappingsToMedian[j];
					mappingsFromMedian[j] = mapfrom_temp[j];
					mappingsToMedian[j] = mapto_temp[j];
				}
				medianG=mediantmp;
				mediantmp=NULL;
				SOD=SOD_temp;
				if (current_projection != new_projection){
					delete[] current_projection;
					current_projection = new_projection;
				}
				std::cout << "---> median updated" << std::endl;
			}
			else{
				for (int j=0; j<N; j++){
					delete[] mapfrom_temp[j];
					delete[] mapto_temp[j];
				}
				mediantmp=NULL;
				std::cout << "---> median not updated" << std::endl;
			}
			delete[] mapfrom_temp;
			delete[] mapto_temp;
			if (refine){
				Graph<NodeAttribute,EdgeAttribute>* local_median = computeMedianGraph(dataset,ed,mcf,medianG,mappingsFromMedian,mappingsToMedian, distancesFromMedian);

			}


			//break;
		}

		// Cleaning memory
		//delete dataset;
		//for (int i = 0;i<N;i++){
		//    for (int j = 0;j<N;j++){
		/*
            std::cout << "cleaning mapping between " << i << " and " << j << " = [";
            for (int jj = 0; jj < N;jj++){
                std::cout << AllPairMappings[sub2ind(i,j,N)][jj] << " ";
            }
            std::cout<< "]"<<std::endl;
		 */
		//       delete[] AllPairMappings[sub2ind(i,j,N)];
		//   }
		// }
		//delete[] distances;
		return medianG;



		}

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
std::vector<Graph<NodeAttribute, EdgeAttribute> *> ComputeProjectedGraphs(int* indices_uplet, int size_uplet, Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed, MedianLabel< NodeAttribute, EdgeAttribute, PropertyType> * mcf, int ** AllPairMappings, double * distances, double * median_vector, double * median_vector_out ){
	std::vector<Graph<NodeAttribute, EdgeAttribute> *> projected_graphs;
	int N = dataset->size();
	double SOD = 0;
	double* distancesFromMedian = new double [N];
	for (int i = 0;i<N;i++){
		SOD+=distances[sub2ind(indices_uplet[0],i,N)];
		distancesFromMedian[i]=distances[sub2ind(indices_uplet[0],i,N)];
	}
	//std::cout<<"initial SOD of Graph " << indices_uplet[0] << " = " << SOD << std::endl;
	// Copying relevant mappings to new pointers
	int * current_mappingFromMedian = new int[dataset->getGraph(indices_uplet[0])->Size()];
	int * current_mappingToMedian = new int[dataset->getGraph(indices_uplet[1])->Size()];
	for (int i = 0 ; i < dataset->getGraph(indices_uplet[0])->Size() ; i++ ){
		current_mappingFromMedian[i]=AllPairMappings[sub2ind(indices_uplet[0],indices_uplet[1],N)][i];
	}
	for (int i = 0 ; i < dataset->getGraph(indices_uplet[1])->Size() ; i++ ){
		current_mappingToMedian[i]=AllPairMappings[sub2ind(indices_uplet[1],indices_uplet[0],N)][i];
	}
	double* current_projection = new double[N];
	double* current_vector = new double[N];
	double * new_projection=new double[N];
	for (int i = 0;i<N;i++){
		current_projection[i]=distances[sub2ind(indices_uplet[0],i,N)];
	}
	//Graph<int,int> * medianG = new Graph<int,int>(*(test_set->getGraph(sorted_indices[0])));
	Graph<NodeAttribute, EdgeAttribute> * medianG = new Graph<NodeAttribute, EdgeAttribute>(*(dataset->getGraph(indices_uplet[0])));
	Graph<NodeAttribute, EdgeAttribute> * mediantmp = NULL;
	MedianGraph<NodeAttribute, EdgeAttribute, PropertyType>* mg = new MedianGraph<NodeAttribute, EdgeAttribute, PropertyType>(mediantmp,dataset,mcf);
	double alpha=0;
	//recursion
	for (int i = 0; i< size_uplet-1;i++){
		// compute new projected vector and corresponding alpha
		for (int j = 0;j<N;j++){
			current_vector[j]=distances[sub2ind(indices_uplet[i+1],j,N)];
		}
		alpha = computeAlpha_and_ProjectedVector(median_vector,current_projection,current_vector,N,new_projection);
		//alpha = 0;
		if (current_projection != new_projection){
			delete[] current_projection;
			current_projection = new_projection;
		}
		//std::cout << "alpha between current projection and graph " << indices_uplet[i+1] << " = " << alpha << std::endl;
		//create new weighted mean Graph
		//medianG=mg->ComputeWeightedMeanGraph(sorted_indices[i],sorted_indices[i+1],mappingsFromMedian[sorted_indices[i+1]],mappingsToMedian[sorted_indices[i+1]],alpha*distancesFromMedian[sorted_indices[i+1]],ed);
		mediantmp=mg->ComputeWeightedMeanGraph(medianG,(*dataset)[indices_uplet[i+1]],current_mappingFromMedian,current_mappingToMedian,alpha*distancesFromMedian[indices_uplet[i+1]],ed);
		//std::cout << "Weighted mean between current projection and "<< indices_uplet[i+1] << " computed." << std::endl;
		//delete medianG;
		medianG=mediantmp;
		projected_graphs.push_back(mediantmp);
		mediantmp=NULL;
		//compute assignments and SOD
		if (i<size_uplet-2){
			//ed->getOptimalMapping(medianG,g2,mappingsFromMedian[j],mappingsToMedian[j]);
			//std::cout << "ok A graph" << j << " , median Size = " << n << " , current graph Size = " << m << std::endl;
			//printGraph(*medianG);
			delete[] current_mappingFromMedian;
			delete[] current_mappingToMedian;
			current_mappingFromMedian = new int[medianG->Size()];
			current_mappingToMedian = new int[dataset->getGraph(indices_uplet[i+2])->Size()];
			ed->getOptimalMapping(medianG,dataset->getGraph(indices_uplet[i+2]),current_mappingFromMedian,current_mappingToMedian);
			//std::cout << "ok B" << std::endl;
			//distancesFromMedian[j] = ed->GedFromMapping(medianG,g2,mappingsFromMedian[j],n,mappingsToMedian[j],m);
			distancesFromMedian[i] = ed->GedFromMapping(medianG,dataset->getGraph(indices_uplet[i+2]),current_mappingFromMedian,medianG->Size(),current_mappingToMedian,dataset->getGraph(indices_uplet[i+2])->Size());
			//std::cout << "ok C" << std::endl;

		}
	}
	for (int i = 0; i < N; i++){
		median_vector_out[i]=  new_projection[i];
	}
	delete[] current_mappingFromMedian;
	delete[] current_mappingToMedian;
	return projected_graphs;
}



template <class NodeAttribute, class EdgeAttribute, class PropertyType>
Graph<NodeAttribute, EdgeAttribute> * computeRecursiveLinearMedian(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		MedianLabel< NodeAttribute, EdgeAttribute, PropertyType> * mcf, int ** AllPairMappings, double * distances, int **&mappingsFromMedian,int **&mappingsToMedian,double *& distancesToMedian,int size_uplets, bool best = 0, bool refine = 0, double time_limit = 600.0 ,int max_uplets = 0)
		{      //std::cout << "time limit = " << time_limit << std::endl;
	struct timeval  tv1, tv2;
	gettimeofday(&tv1, NULL);

	MedianGraph<NodeAttribute, EdgeAttribute, PropertyType>* mg = new MedianGraph<NodeAttribute, EdgeAttribute, PropertyType>(nullptr,dataset,mcf);
	int N = dataset->size();
	double* median_vector = computeGeometricMedian(distances,N,N);
	double* vec_distances = new double[N];
	int * sorted_indices = new int[N];
	double SOD_refined_out = -1.0;
	double SOD_out = -1.0;
	int order_out = 0;
	int order_refined_out = 0;
	int i_out;
	int j_out;
	Graph<NodeAttribute, EdgeAttribute> * median_graph_out = nullptr;

	for (int i=0;i<N;i++){
		sorted_indices[i]=-1;
	}
	for (int i = 0;i<N;i++){
		//compute distance between graph i and median in vector space
		vec_distances[i]=0;
		for (int j = 0;j<N;j++){
			vec_distances[i]+= (median_vector[j]-distances[sub2ind(i,j,N)])*(median_vector[j]-distances[sub2ind(i,j,N)]);
		}
		// std::cout << "Squared distance from Graph " << i << " = " << vec_distances[i] << std::endl;
		vec_distances[i]=std::sqrt(vec_distances[i]);
		//std::cout << "distance from Graph " << i << " = " << vec_distances[i] << std::endl;
		//checks if new_graph belongs in the max_recursion+1 closest to median, and updates the sorted list.
		for (int j = 0; j < N;j++){
			if (sorted_indices[j] == -1){
				sorted_indices[j]=i;
				j=N;
			}
			else if (vec_distances[i] < vec_distances[sorted_indices[j]]){
				for (int k= N-1;k>j;k--){
					sorted_indices[k]=sorted_indices[k-1];
				}
				sorted_indices[j]=i;
				j=N;
			}
		}
	}
	// Generation of uplets
	int index_uplet = 0;
	int num_uplets;
	if (max_uplets ==0)
		num_uplets = N/size_uplets;
	else
		num_uplets =std::min(max_uplets,  N/size_uplets);
	int * * uplets = new int * [num_uplets];
	while ((index_uplet+1)*size_uplets <= N && index_uplet < num_uplets){
		uplets[index_uplet] = new int[size_uplets];
		for (int i = 0;i<size_uplets;i++){
			uplets[index_uplet][i] = sorted_indices[index_uplet*size_uplets + i];
		}
		index_uplet ++;
	}
	/*
      std::cout << "uplets generated = " << std::endl;
      for (int i = 0;i<num_uplets;i++){
        std::cout << "[" ;
        for (int j =0;j<size_uplets;j++){
            std::cout << uplets[i][j] << " ";
        }
        std::cout << "]" <<std::endl;
      }
	 */
	//initialization$
	double alpha;
	Graph<NodeAttribute, EdgeAttribute> * v_i_Graph = nullptr;
	double * v_i = new double[N];
	std::vector<Graph<NodeAttribute, EdgeAttribute> *>* candidate_Graphs = new std::vector<Graph<NodeAttribute, EdgeAttribute> *>[num_uplets];
	// eventually contains num_uplet vectors having each size_uplets graph pointers except for the first which has size_uplets-1 graph pointers.
	//

	//main loop


	for (int i = 0;i<num_uplets;i++){
		double * v_i_minus_1= new double [N];
		double * uplet_median_vector = new double [N];
		for (int j = 0;j<N;j++){
			v_i_minus_1[j]=v_i[j];
		}
		//computing candidates
		if (i==0){
			candidate_Graphs[i]=ComputeProjectedGraphs(uplets[i],size_uplets,dataset,ed,mcf,AllPairMappings,distances,median_vector,v_i);
		}
		if (i>0){
			candidate_Graphs[i]=ComputeProjectedGraphs(uplets[i],size_uplets,dataset,ed,mcf,AllPairMappings,distances,median_vector,uplet_median_vector);
			alpha = computeAlpha_and_ProjectedVector(median_vector,v_i_minus_1,uplet_median_vector,N,v_i);
			Graph<NodeAttribute, EdgeAttribute> * v_i_minus_1_Graph = candidate_Graphs[i-1].back();
			Graph<NodeAttribute, EdgeAttribute> * uplet_median_vector_Graph = candidate_Graphs[i].back();
			int n = v_i_minus_1_Graph->Size();
			int m = uplet_median_vector_Graph->Size();
			int *mapfrom = new int[n];
			int *mapto = new int[m];
			ed->getOptimalMapping(v_i_minus_1_Graph,uplet_median_vector_Graph,mapfrom,mapto);
			double current_distance = ed->GedFromMapping(v_i_minus_1_Graph,uplet_median_vector_Graph,mapfrom,n,mapto,m);
			candidate_Graphs[i].push_back(mg->ComputeWeightedMeanGraph(v_i_minus_1_Graph,uplet_median_vector_Graph,mapfrom,mapto,alpha*current_distance,ed));
			delete[] mapfrom;
			delete[] mapto;
		}
		delete[] v_i_minus_1;
		delete[] uplet_median_vector;

		//computing SODs and keeping best
		int current_candidate_set_size = candidate_Graphs[i].size();
		for (int j = 0;j<current_candidate_set_size;j++){
			Graph<NodeAttribute, EdgeAttribute> * g1 = candidate_Graphs[i][j];
			//std::cout<< "current candidate =" << std::endl;
			//printGraph(*g1);
			int **mapsfrom = new int*[N];
			int **mapsto = new int*[N];
			int n = g1->Size();
			double SOD_refined = 0;
			double SOD = 0;
			double * local_distances = new double[N];
			for (int k =0;k<N;k++){
				Graph<NodeAttribute, EdgeAttribute> * g2 = dataset->getGraph(k);
				int m = g2->Size();
				mapsfrom[k] = new int[n];
				mapsto[k] = new int[m];
				//std::cout << "computing distance from candidate "<< j << " with graph "<< k << std::endl;
				ed->getOptimalMapping(g1,g2,mapsfrom[k],mapsto[k]);
				local_distances[k] = ed->GedFromMapping(g1,g2,mapsfrom[k],n,mapsto[k],m);
				SOD+= local_distances[k];
			}
			//std::cout << "SOD of candidate "<< j << " of "<< i << "th uplet = " << SOD << std::endl;
			if (SOD_out == -1.0 || SOD_out > SOD){
				if (SOD_out != -1.0){
					delete[] distancesToMedian;
					for (int k =0;k<N;k++){
						delete[] mappingsFromMedian[k];
						delete[] mappingsToMedian[k];
					}
					delete[] mappingsFromMedian;
					delete[] mappingsToMedian;
				}
				distancesToMedian = local_distances;
				mappingsFromMedian = mapsfrom;
				mappingsToMedian = mapsto;
				SOD_out = SOD;
				i_out = i;
				j_out = j;
				order_out= g1->Size();
				median_graph_out = new Graph<NodeAttribute, EdgeAttribute>(*g1);
				//std::cout<< "current best candidate with order " << std::endl;
				//printGraph(*median_graph_out);
			}
			else{
				for (int k =0;k<N;k++){
					delete[] mapsfrom[k];
					delete[] mapsto[k];
				}
				delete[] local_distances;
				delete[] mapsfrom;
				delete[] mapsto;
			}
			// checking time limit
			gettimeofday(&tv2, NULL);
			double time=((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));

			if (time > time_limit){
				j=current_candidate_set_size;
				i=num_uplets;// get out of the loop
			}

		}
	}
	delete[] candidate_Graphs;
	std::cout << SOD_out <<",";
	//std::cout << "final solution for current projections computed =" << std::endl;
	//printGraph(*median_graph_out);
	return median_graph_out;
		}


template <class NodeAttribute, class EdgeAttribute, class PropertyType>
Graph<NodeAttribute, EdgeAttribute> * computeRecursiveWeightedMins(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		MedianLabel< NodeAttribute, EdgeAttribute, PropertyType> * mcf, int max_recursion, int ** AllPairMappings, double * distances, bool refine = 0)
		{      int N = dataset->size();

		double* median_vector = computeGeometricMedian(distances,N,N);

		double* vec_distances = new double[N];
		int * sorted_indices = new int[max_recursion+1];
		for (int i=0;i<max_recursion+1;i++){
			sorted_indices[i]=-1;
		}
		for (int i = 0;i<N;i++){


			vec_distances[i]=0;
			for (int j = 0;j<N;j++){
				vec_distances[i]+= (median_vector[j]-distances[sub2ind(i,j,N)])*(median_vector[j]-distances[sub2ind(i,j,N)]);
			}
			//std::cout << "Squared distance from Graph " << i << " = " << vec_distances[i] << std::endl;
			vec_distances[i]=std::sqrt(vec_distances[i]);
			//std::cout << "distance from Graph " << i << " = " << vec_distances[i] << std::endl;
			//checks if new_graph belongs in the max_recursion+1 closest to median, and updates the sorted list.
			for (int j = 0; j < max_recursion+1;j++){
				if (sorted_indices[j] == -1){
					sorted_indices[j]=i;
					j=max_recursion+1;
				}
				else if (vec_distances[i] < vec_distances[sorted_indices[j]]){
					for (int k= max_recursion;k>j;k--){
						sorted_indices[k]=sorted_indices[k-1];
					}
					sorted_indices[j]=i;
					j=max_recursion+1;
				}
			}
		}
		//computation of projections and corresponding graphs
		double * median_vector_out = new double[N];
		std::vector<Graph<NodeAttribute, EdgeAttribute> *> candidate_Graphs = ComputeProjectedGraphs(sorted_indices,max_recursion+1,dataset,ed,mcf,AllPairMappings,distances,median_vector,median_vector_out);

		//SOD computations and keeping best
		Graph<NodeAttribute, EdgeAttribute> * median_graph_out = nullptr;
		double SOD_out = -1.0;
		for (int j = 0;j<candidate_Graphs.size();j++){
			Graph<NodeAttribute, EdgeAttribute> * g1 = candidate_Graphs[j];
			int n = g1->Size();
			double SOD = 0;
			for (int k =0;k<N;k++){
				Graph<NodeAttribute, EdgeAttribute> * g2 = dataset->getGraph(k);
				int m = g2->Size();
				int *mapfrom = new int[n];
				int *mapto = new int[m];
				ed->getOptimalMapping(g1,g2,mapfrom,mapto);
				SOD+= ed->GedFromMapping(g1,g2,mapfrom,n,mapto,m);
				delete[] mapfrom;
				delete[] mapto;
			}
			if (SOD_out == -1.0 || SOD_out > SOD){
				SOD_out = SOD;
				median_graph_out = new Graph<NodeAttribute, EdgeAttribute>(g1);
			}
			delete g1;
		}

		//break;
		/*
           if (refine){
            Graph<NodeAttribute,EdgeAttribute>* local_median = computeMedianGraph(dataset,ed,mcf,medianG,mappingsFromMedian,mappingsToMedian, distancesFromMedian);

           }
		 */
		std::cout << SOD_out ;
		return median_graph_out;
		}





/*
template<typename T1, typename T2>
void PrintPointerContent(T1 * pointer, T2 length){
    std::cout<<"[";
    for (int j = 0; j < length;j++){
        std::cout << pointer[j] << " ";
    }
    std::cout<< "]"<<std::endl;
}


/*void PrintPointerContent(double * pointer, int length){
    std::cout<<"[";
    for (int j = 0; j < length;j++){
        std::cout << pointer[j] << " ";
    }
    std::cout<< "]"<<std::endl;
}*/
/*
template <class NodeAttribute, class EdgeAttribute, class PropertyType>
Graph<NodeAttribute, EdgeAttribute> * computeStochasticInitModel(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset, double size_correction = 1){
  int N = dataset->size();
  for (int i =0;i<N;i++){
    Graph<NodeAttribute,EdgeAttribute> * g = dataset->getGraph(i);
    int n = g->Size();


  }

}
 */



template <class NodeAttribute, class EdgeAttribute, class PropertyType>
Graph<NodeAttribute, EdgeAttribute> * computeMedianGraph(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		MedianLabel< NodeAttribute, EdgeAttribute, PropertyType> * mcf, int SetMedianIndex, int * * mappingsFromMedian, int * * mappingsToMedian, double * distancesToMedian, double * SOD = nullptr, int * num_it = nullptr){

	Graph<NodeAttribute, EdgeAttribute>  *Gmed = dataset->getGraph(SetMedianIndex);  // pointer on the graph of index bestSODIndex
	return computeMedianGraph(dataset,ed,mcf,Gmed,mappingsFromMedian,mappingsToMedian,distancesToMedian,SOD,num_it);
}

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
Graph<NodeAttribute, EdgeAttribute> * computeMedianGraph(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		MedianLabel< NodeAttribute, EdgeAttribute, PropertyType> * mcf, int num_graphs, GraphEditDistance<NodeAttribute, EdgeAttribute> * final_ed){
	int numgraphs;
	struct timeval  tv1, tv2;
	int N = dataset->size();
	int ** AllPairMappings = new int*[N*N];
	double  * distances = new double[N*N];
	bool * map_computed = new bool[N*N];
	int ** mappingsFromMedian = new int*[N];
	int ** mappingsToMedian = new int*[N];

	for (int i=0; i<N; i++){
		for (int j=i+1; j<N; j++){
			map_computed[sub2ind(i,j,N)] = false;
			map_computed[sub2ind(j,i,N)] = false;

		}
	}


	std::default_random_engine randGen;
	randGen.seed(123);
	//std::cout << "proportion = " << proportion_set << " ;     taille dataset = " << N <<  " ;   num_graphs = " << numgraphs << std::endl;
	int * indices = new int[N];
	for (int i=0;i<N;i++){
		indices[i]=i;
	}

	std::shuffle(&indices[0], &indices[N-1], randGen);
	int num_it_out;
	int ind_out;
	double SOD_out = -1.0;
	double SOD_unrefined_out = -1.0;
	int refined_order_out=0;
	Graph<NodeAttribute, EdgeAttribute>  * graph_out = nullptr;
	int * init_indices ;
	//set initial graphs
	if (num_graphs < 1){//pick 3 graphs, the biggest, the smallest and one with mean order
		init_indices = new int[3];
		int max_index= indices[0];
		int max_order = dataset->getGraph(indices[0])->Size();
		int min_index = indices[0];
		int min_order = max_order;
		double mean_order = dataset->getGraph(indices[0])->Size();
		for (int i = 1; i <N;i ++){
			int index = indices[i];
			int current_order = dataset->getGraph(index)->Size();
			if (current_order < min_order){
				min_order = current_order;
				min_index = index;
			}
			else  if (current_order > max_order){
				max_order = current_order;
				max_index = index;
			}
			mean_order+= current_order;
		}
		mean_order/=N;
		int mean_order_int = (int)std::round(mean_order);
		int min_sq_distance_to_mean = (mean_order - dataset->getGraph(0)->Size()) * (mean_order - dataset->getGraph(0)->Size()) ;
		int mean_index = 0;
		for (int i = 1; i <N;i ++){
			int index = indices[i];
			int current_order = dataset->getGraph(index)->Size();
			double sq_distance_to_mean = (mean_order - current_order)* (mean_order - current_order);
			if (sq_distance_to_mean < min_sq_distance_to_mean){
				mean_index = index;
				min_sq_distance_to_mean = sq_distance_to_mean;
			}
			if (sq_distance_to_mean == 0){
				i = N;
			}
		}
		init_indices[0] = min_index;
		init_indices[1] = max_index;
		init_indices[2] = mean_index;
		numgraphs = 3;
	}
	else{
		numgraphs=std::min(num_graphs,N);
		init_indices = new int[numgraphs];
		for (int i=0;i<numgraphs;i++){
			init_indices[i] = indices[i];
		}
	}

	for (int i=0;i<numgraphs;i++){
		//std::cout << "graph " << i << std::endl;
		gettimeofday(&tv1, NULL);
		double* distancesToMedian = new double[N];
		Graph<NodeAttribute,EdgeAttribute>  *g1 = new Graph<NodeAttribute,EdgeAttribute> (*dataset->getGraph(init_indices[i]));
		int n=g1->Size();
		distances[sub2ind(init_indices[i],init_indices[i],N)] = 0;
		int * selfAssignement = new int[n];
		for (int j=0; j<n; j++) {
			selfAssignement[j]=j;
		}
		AllPairMappings[sub2ind(init_indices[i],init_indices[i],N)]=selfAssignement;
		map_computed[sub2ind(init_indices[i],init_indices[i],N)]= true;
		double SOD = 0;
		double SOD_unrefined=0;
		int num_it =0;
		for (int j=0; j<N; j++){
			Graph<NodeAttribute,EdgeAttribute>  *g2 = dataset->getGraph(j);  // pointer on the graph j
			int m=g2->Size();
			if (!map_computed[sub2ind(init_indices[i],j,N)]){
				AllPairMappings[sub2ind(init_indices[i],j,N)] = new int[n];
				AllPairMappings[sub2ind(j,init_indices[i],N)] = new int[m];
				ed->getOptimalMapping(g1,g2,AllPairMappings[sub2ind(init_indices[i],j,N)],AllPairMappings[sub2ind(j,init_indices[i],N)]);
				distances[sub2ind(init_indices[i],j,N)] = ed->GedFromMapping(g1,g2,AllPairMappings[sub2ind(init_indices[i],j,N)],n,AllPairMappings[sub2ind(j,init_indices[i],N)],m);
				distances[sub2ind(j,init_indices[i],N)] = distances[sub2ind(init_indices[i],j,N)];
				map_computed[sub2ind(init_indices[i],j,N)]= true;
				map_computed[sub2ind(j,init_indices[i],N)]= true;
			}
			mappingsFromMedian[j] = new int[n];
			mappingsToMedian[j] = new int[m];
			memcpy(mappingsFromMedian[j],AllPairMappings[sub2ind(init_indices[i],j,N)],sizeof(int)*n);
			memcpy(mappingsToMedian[j],AllPairMappings[sub2ind(j,init_indices[i],N)],sizeof(int)*m);
			distancesToMedian[j]=distances[sub2ind(init_indices[i],j,N)];
			SOD_unrefined += distancesToMedian[j];
		}
		//std::cout << " maps computed " << std::endl;
		std::cout << (*dataset)(0)<< " " << N << " ";
		Graph<NodeAttribute, EdgeAttribute>  *Gmed = computeMedianGraph(dataset,ed,mcf,g1,mappingsFromMedian,mappingsToMedian,distancesToMedian,&SOD,&num_it);
		//for(int j=0;j<N;j++){
		//  SOD+=distancesToMedian[j];
		//}
		//std::cout << " median computed " << std::endl;

		if(SOD_out==-1.0 || SOD_out > SOD){
			SOD_out = SOD;
			num_it_out = num_it;
			SOD_unrefined_out=SOD_unrefined;
			graph_out = new  Graph<NodeAttribute, EdgeAttribute> (Gmed);
			ind_out = init_indices[i];
			refined_order_out = Gmed->Size();
		}

		gettimeofday(&tv2, NULL);
		double time=((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
		double finalSOD = computeSOD(dataset,final_ed,Gmed);
		//std::cout << " SOD computed " << std::endl;
		std::cout << num_it<< " "<< finalSOD << " "<< time <<" " <<  0 << " " << SOD << " 1"<< std::endl;
		delete Gmed;

	}
	//clean memory !!!
	//std::cout<<  SOD_unrefined_out << " ; " << num_it_out << " ; " << SOD_out << " ; ";
	return graph_out;

}

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
Graph<NodeAttribute, EdgeAttribute> * computeMedianGraph(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		MedianLabel< NodeAttribute, EdgeAttribute, PropertyType> * mcf, Graph<NodeAttribute, EdgeAttribute> * Gmed, int * * mappingsFromMedian, int * * mappingsToMedian, double * distancesToMedian, double * p_SOD = nullptr, int* p_num_it = nullptr){

	//std::cout <<"\n ----------\n starting generalized computation" << endl;
	struct timeval  tv1, tv2;
#ifdef PRINT_TIMES
	gettimeofday(&tv1, NULL);
#endif
	double SOD =0;

	int N = dataset->size();
	for (int i = 0;i<N;i++){
		SOD+=distancesToMedian[i];

	}
	std::cout << SOD << " ";
	//std::cout << "initial distances";
	//PrintPointerContent(distancesToMedian,N);




	//Graph<NodeAttribute, EdgeAttribute>  *Gmed = dataset->getGraph(SetMedianIndex);  // pointer on the graph of index bestSODIndex
	Graph<NodeAttribute, EdgeAttribute> *medianG = new Graph<NodeAttribute, EdgeAttribute>(*Gmed);  // copy to a new graph
	MedianGraph<NodeAttribute, EdgeAttribute,PropertyType>* mg = new MedianGraph<NodeAttribute,EdgeAttribute,PropertyType>(medianG,dataset,mcf);
	//printGraph(*medianG);

	int nMedian = medianG->Size();
	/*
 std::cout << "mappingsfrom =" << std::endl;
    for (int i = 0; i<N;i++){
            PrintPointerContent(mappingsFromMedian[i],nMedian);
    }
  std::cout << "mappingsTo =" << std::endl;
    for (int i = 0; i<N;i++){
            PrintPointerContent(mappingsToMedian[i],dataset->getGraph(i)->Size());
    }
	 */

	// update vertex-label and edge-set of the median graph
	bool isModifiedGraph = mg->updateMedianGraph( mappingsFromMedian); // isModified is false when updating the median does not modifiy it.
	medianG=mg->getGraph();
	//printGraph(*medianG);

	bool ismodifiedEP = true;
	int num_it=0;
	int num_it_wo_modification = 0;
	while (num_it_wo_modification<=3 && (ismodifiedEP || isModifiedGraph)){
		num_it ++;
		/*
   std::cout << "mappingsfrom =" << std::endl;
    for (int i = 0; i<N;i++){
            PrintPointerContent(mappingsFromMedian[i],nMedian);
    }
  std::cout << "mappingsTo =" << std::endl;
    for (int i = 0; i<N;i++){
            PrintPointerContent(mappingsToMedian[i],dataset->getGraph(i)->Size());
    }
		 */

		ismodifiedEP = false;
		medianG = mg->getGraph();
		nMedian = medianG->Size();
		// updating the SOD
		SOD = 0;
		for (int i=0;i<N;i++){
			Graph<NodeAttribute, EdgeAttribute>  *gi = dataset->getGraph(i);
			int m = gi->Size();


			distancesToMedian[i] = ed->GedFromMapping(medianG,gi,mappingsFromMedian[i],nMedian,mappingsToMedian[i],m);

			SOD+= distancesToMedian[i];
		}
		//std::cout << "updated distances";
		//PrintPointerContent(distancesToMedian,N);
		// updating the edit paths
		for (int i=0;i<N;i++){
			Graph<NodeAttribute, EdgeAttribute>  *gi = dataset->getGraph(i);
			int m = gi->Size();
			int * localMappingToMedian = new int[m];
			int * localMappingFromMedian = new int[medianG->Size()];
			ed->getOptimalMapping(medianG,gi,localMappingFromMedian,localMappingToMedian);
			double localDistance = ed->GedFromMapping(medianG,gi,localMappingFromMedian,nMedian,localMappingToMedian,m);
			if (localDistance< distancesToMedian[i]){
				//std::cout << "better mapping found " << std::endl;
				//à remettre
				delete[] mappingsFromMedian[i];
				delete[] mappingsToMedian[i];
				mappingsFromMedian[i] = localMappingFromMedian;
				mappingsToMedian[i] = localMappingToMedian;
				ismodifiedEP = true;
			}
			else {
				//std::cout << "no better mapping found " << std::endl;
				// à remetre
				delete [] localMappingFromMedian;
				delete [] localMappingToMedian;
			}

		}


		//updating the median graph
		//std::cout << "updating median" << std::endl;
		isModifiedGraph = mg->updateMedianGraph( mappingsFromMedian);
		//std::cout << "updated median" << std::endl;
		//if  (ismodifiedEP) std::cout <<endl<< " EP modified " << endl;
		//if  (isModifiedGraph) std::cout << " Median Graph modified " << endl;

		if(isModifiedGraph == false){
			//std::cout << " Checking Node Deletion " << endl;
			bool isModifiedGraphSize =true;
			while (isModifiedGraphSize){
				//std::cout<< "nodedeletioncheck" << std::endl;
				isModifiedGraphSize = mg->checkNodeDeletion(mappingsFromMedian,mappingsToMedian,nMedian);
				//std::cout << "deletion checked ------>" << isModifiedGraphSize << std::endl;
				nMedian=mg->getGraph()->Size();
				if (isModifiedGraphSize)
					isModifiedGraph=true;
			}
			if (!isModifiedGraph){
				//isModifiedGraphSize = 0;
				isModifiedGraphSize =true;
				while (isModifiedGraphSize){
					//std::cout<< "nodeinsertioncheck" << std::endl;
					isModifiedGraphSize = mg->checkNodeInsertion(mappingsFromMedian,mappingsToMedian,nMedian);
					medianG = mg->getGraph();
					nMedian = medianG->Size();
					if (isModifiedGraphSize)
						isModifiedGraph=true;
				}

				medianG = mg->getGraph();
				nMedian = medianG->Size();
				//std::cout << "insertion checked ------>" << isModifiedGraphSize << std::endl;

			}
			if (isModifiedGraph){
				medianG = mg->getGraph();
				nMedian = medianG->Size();
				//printGraph(*medianG);

			}
			//std::cout << " Node Deletion Checked new graph in while loop --->" << endl;
			//printGraph(*medianG);
		}

		if(isModifiedGraph == false){
			num_it_wo_modification++;
		}
		else{
			num_it_wo_modification=0;
		}


	}
	//std::cout <<endl<< "---->refined SOD = " << SOD <<" , graph order = "<< nMedian << endl;
	//std::cout << SOD << ";";
	if (p_SOD != nullptr)
		*p_SOD = SOD;
	//std::cout << num_it << " ; ";
	if (p_num_it != nullptr)
		*p_num_it = num_it;
	return medianG;
}


template <class NodeAttribute, class EdgeAttribute, class PropertyType>
int computeSetMedianGraphInitializedDistances(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		int * * mappingsFromMedian, int * * mappingsToMedian, double * distancesToMedian, int** AllPairMappings, double * distances){
	struct timeval  tv1, tv2;
	//std::cout << "starting set_median computation" << endl;
	double SOD;
	double bestSOD;
	int bestSODIndex = 0;
	int N = dataset->size();
	for (int i=0; i<N; i++){
		SOD = 0;
		for (int j=0; j<N; j++) SOD+=distances[sub2ind(i,j,N)];
		if (i==0 || SOD < bestSOD){
			bestSOD = SOD;
			bestSODIndex = i;
		}
	}
	//std::cout << "set_medianIndex = "<< bestSODIndex << " , SOD = " << bestSOD << std::endl ;
	//std::cout <<  bestSOD << " ; ";
	Graph<NodeAttribute, EdgeAttribute>  *Gmed = dataset->getGraph(bestSODIndex);  // pointer on the graph of index bestSODIndex
	//std::cout <<endl<< "index of set median = " << bestSODIndex << ", SOD = " << bestSOD << endl;
	// creating and adding the self mapping into AllpairMappings only for set-median graph

	for (int i=0; i<N; i++) {
		mappingsFromMedian[i]=new int[Gmed->Size()];
		for (int j = 0;j<Gmed->Size();j++){
			mappingsFromMedian[i][j]=AllPairMappings[sub2ind(bestSODIndex,i,N)][j];
		}
		mappingsToMedian[i]= new int [dataset->getGraph(i)->Size()];
		for (int j = 0; j <dataset->getGraph(i)->Size();j++){
			mappingsToMedian[i][j]=AllPairMappings[sub2ind(i,bestSODIndex,N)][j];
		}
		distancesToMedian[i]=distances[sub2ind(i,bestSODIndex,N)];
	}
	return bestSODIndex;
}


template <class NodeAttribute, class EdgeAttribute, class PropertyType>
int computeSetMedianGraph(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		int * * mappingsFromMedian, int * * mappingsToMedian, double * distancesToMedian){
	struct timeval  tv1, tv2;
	//std::cout << "starting set_median computation" << endl;
#ifdef PRINT_TIMES
	gettimeofday(&tv1, NULL);
#endif
	double SOD;
	double bestSOD;
	int bestSODIndex = 0;

	int N = dataset->size();
	double* distances = new double[N*N];
	int * * AllPairMappings = new int*[N*N];

	computeMappingsAndDistances(dataset,ed,AllPairMappings,distances);
	int SetmedianIndex = computeSetMedianGraphInitializedDistances(dataset,ed,mappingsFromMedian,mappingsToMedian,AllPairMappings,distances);
	delete[] distances;
	for (int i = 0; i < N*N;i++){
		delete[] AllPairMappings[i];
	}
	delete[] AllPairMappings;
	return SetmedianIndex;
}




template <class NodeAttribute, class EdgeAttribute, class PropertyType>
int computeBiggestGraph(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		int * * mappingsFromMedian, int * * mappingsToMedian, double * distancesToMedian){
	int biggestOrder = (*dataset)[0]->Size();
	int biggestOrderIndex = 0;
	int N = dataset->size();
	for (int i=1;i<N;i++){
		if((*dataset)[i]->Size()>biggestOrder){
			biggestOrder=(*dataset)[i]->Size();
			biggestOrderIndex=i;
		}
	}
	//std::cout << "biggest order index found = " << biggestOrderIndex << std::endl;
	Graph<NodeAttribute, EdgeAttribute> * g1 = (*dataset)[biggestOrderIndex];
	for (int i=0;i<N;i++){
		//std::cout << "i = " << i << std::endl;
		Graph<NodeAttribute, EdgeAttribute> * g2 = (*dataset)[i];
		int m = g2->Size();
		mappingsFromMedian[i]=new int[biggestOrder];
		mappingsToMedian[i]=new int[m];
		ed->getOptimalMapping(g1,g2,mappingsFromMedian[i],mappingsToMedian[i]);
		//std::cout << "mapping computed " << i << std::endl;
		distancesToMedian[i]= ed->GedFromMapping(g1,g2,mappingsFromMedian[i],biggestOrder,mappingsToMedian[i],m);
	}
	return biggestOrderIndex;
}

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
int computeSmallestGraph(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		int * * mappingsFromMedian, int * * mappingsToMedian, double * distancesToMedian){
	int smallestOrder = (*dataset)[0]->Size();
	int smallestOrderIndex = 0;
	int N = dataset->size();
	for (int i=1;i<N;i++){
		if((*dataset)[i]->Size()<smallestOrder){
			smallestOrder=(*dataset)[i]->Size();
			smallestOrderIndex=i;
		}
	}
	//std::cout << "biggest order index found = " << biggestOrderIndex << std::endl;
	Graph<NodeAttribute, EdgeAttribute> * g1 = (*dataset)[smallestOrderIndex];
	for (int i=0;i<N;i++){
		//std::cout << "i = " << i << std::endl;
		Graph<NodeAttribute, EdgeAttribute> * g2 = (*dataset)[i];
		int m = g2->Size();
		mappingsFromMedian[i]=new int[smallestOrder];
		mappingsToMedian[i]=new int[m];
		ed->getOptimalMapping(g1,g2,mappingsFromMedian[i],mappingsToMedian[i]);
		//std::cout << "mapping computed " << i << std::endl;
		distancesToMedian[i]= ed->GedFromMapping(g1,g2,mappingsFromMedian[i],smallestOrder,mappingsToMedian[i],m);
	}
	return smallestOrderIndex;
}

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
int computeConsensualGraph(Dataset< NodeAttribute, EdgeAttribute, PropertyType> * dataset,
		GraphEditDistance<NodeAttribute, EdgeAttribute> * ed,
		int * * mappingsFromMedian, int * * mappingsToMedian, double * distancesToMedian){
	int smallestOrder = (*dataset)[0]->Size();
	int smallestOrderIndex = 0;
	int N = dataset->size();
	for (int i=1;i<N;i++){
		if((*dataset)[i]->Size()<smallestOrder){
			smallestOrder=(*dataset)[i]->Size();
			smallestOrderIndex=i;
		}
	}
	//std::cout << "biggest order index found = " << biggestOrderIndex << std::endl;
	Graph<NodeAttribute, EdgeAttribute> * g1 = (*dataset)[smallestOrderIndex];
	for (int i=0;i<N;i++){
		//std::cout << "i = " << i << std::endl;
		Graph<NodeAttribute, EdgeAttribute> * g2 = (*dataset)[i];
		int m = g2->Size();
		mappingsFromMedian[i]=new int[smallestOrder];
		mappingsToMedian[i]=new int[m];
		ed->getOptimalMapping(g1,g2,mappingsFromMedian[i],mappingsToMedian[i]);
		//std::cout << "mapping computed " << i << std::endl;
		distancesToMedian[i]= ed->GedFromMapping(g1,g2,mappingsFromMedian[i],smallestOrder,mappingsToMedian[i],m);
	}
	return smallestOrderIndex;
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
		//  std::cout<< "best distance between representative graph and datagraph is " << best_distance << std::endl;
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
		std::cout << "\n - results for class " <<  (char)properties[p] << " - \n" << "number of graphs in the class = " <<  num_total[p] << " , number of correct classifications = " << num_correct[p] ;
		std::cout << "\ntaux de precision pour la classe = " <<  100*static_cast<double>(num_correct[p])/static_cast<double>(num_total[p]) << "%";
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


// =================================================================

int main (int argc, char* argv[])
{
	struct Options * options =   parseOptions(argc,argv);

	options->k = 3;
	// IBD ------------------------------------------------------------------------------------
	if (options->ibd_costs != "") {

		IBDDistanceCost *cf = new IBDDistanceCost(options->ibd_costs);
		IBDMedianLabel *mcf = new IBDMedianLabel(cf);

		// IPFP used as a refinement method
		IPFPGraphEditDistance<int,double> * algoIPFP = new IPFPGraphEditDistance<int,double>(cf);
		algoIPFP->recenterInit();
		RandomMappingsGED<int,double> *final_init = new RandomMappingsGED<int,double>();
		GraphEditDistance<int,double>* final_ed = new MultistartRefinementGraphEditDistance<int,double>(cf, final_init, 10, algoIPFP,5,options->adap);
		GraphEditDistance<int,double>* ed;
		if(options->method == "lsape_multi_bunke")
			ed = new BipartiteGraphEditDistanceMulti<int,double>(cf, options->nep);
		else {
			RandomMappingsGED<int,double> *init = new RandomMappingsGED<int,double>();
			ed = new MultistartRefinementGraphEditDistance<int,double>(cf, init, options->nep, algoIPFP,options->nrep,options->adap);
		}

		for (int k=0; k<options->num_repetitions_xp; k++) {

			double init_time = 0.0;
			IBDDataset *completeset = new IBDDataset((options->dataset_file).c_str());
			int M = completeset->size();
			int prop_init_random_median = 10;
			struct timeval  tv1, tv2;
#ifdef PRINT_TIMES
			gettimeofday(&tv1, NULL);
#endif
			GraphEditDistance<int,double>* i_ed = nullptr;
			i_ed = new BipartiteGraphEditDistanceMulti<int,double>(cf, options->nep);
			double rec_init_time;
			Graph<int,double>* local_median = nullptr;
			Graph<int,double>* init_median = nullptr;
			double * distancesToMedian = new double [M];
			int * * mappingsFromMedian = new int*[M];
			int * * mappingsToMedian = new int*[M];
			double * p_SOD= new double;
			double finalSOD;
			int * p_num_it = new int;
			int status;
			int ** AllPairMappings =  new int*[M*M];
			double * distances = new double[M*M];
#ifdef PRINT_TIMES
			gettimeofday(&tv1, NULL);
#endif
			bool init_computed = computeMappingsAndDistances(completeset,i_ed,AllPairMappings,distances,options->time_limit);
#ifdef PRINT_TIMES
			gettimeofday(&tv2, NULL);
			init_time=((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
#endif
			std::cout << "Id_collection,algo,SOD,time,refined_SOD,refinement_time,status" << std::endl;
			//Linear
			std::cout << options->collection_id <<"," << options->collection_percentage << ",linear,";
			if(!init_computed){
				std::cout << "-1,"<<init_time <<",-1,"<< init_time << ",0" << std::endl;
			}
			else{
				gettimeofday(&tv1, NULL);
				local_median = computeRecursiveLinearMedian(completeset,ed,mcf,AllPairMappings,distances,mappingsFromMedian,mappingsToMedian, distancesToMedian,2,1,0,options->time_limit-init_time,1);
				gettimeofday(&tv2, NULL);
				//std::cout<<"median_computed"<< std::endl;
				double rec_init_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
				gettimeofday(&tv1, NULL);
				finalSOD=computeSOD(completeset,final_ed,local_median,options->time_limit -(rec_init_time+init_time));
				gettimeofday(&tv2, NULL);
				//std::cout<<"median_refined" << std::endl;
				double refinement_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
				status=3;
				if(rec_init_time+init_time>options->time_limit){
					status=1;
				}
				else if(rec_init_time+init_time+refinement_time>options->time_limit){
					status=2;
					finalSOD=-1;
				}
				std::cout << init_time+rec_init_time << "," << finalSOD <<","<<refinement_time << "," << status << std::endl;
			}
			//std::cout << "descent converged !! :)" << std::endl;
			delete init_median;
			delete local_median;
			//Triangular
			std::cout << options->collection_id <<"," << options->collection_percentage << ",Triangular,";
			if(!init_computed){
				std::cout << "-1,"<<init_time <<",-1,"<< init_time << ",0" << std::endl;
			}
			else{
				gettimeofday(&tv1, NULL);
				local_median = computeRecursiveLinearMedian(completeset,ed,mcf,AllPairMappings,distances,mappingsFromMedian,mappingsToMedian, distancesToMedian,3,1,0,options->time_limit-init_time,1);
				gettimeofday(&tv2, NULL);
				//std::cout<<"median_computed"<< std::endl;
				double rec_init_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
				gettimeofday(&tv1, NULL);
				finalSOD=computeSOD(completeset,final_ed,local_median,options->time_limit -(rec_init_time+init_time));
				gettimeofday(&tv2, NULL);
				//std::cout<<"median_refined" << std::endl;
				double refinement_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
				status=3;
				if(rec_init_time+init_time>options->time_limit){
					status=1;
				}
				else if(rec_init_time+init_time+refinement_time>options->time_limit){
					status=2;
				}
				std::cout << init_time+rec_init_time << "," << finalSOD <<","<<refinement_time << "," << status << std::endl;
			}
			//std::cout << "descent converged !! :)" << std::endl;
			delete init_median;
			delete local_median;
			//std::cout << "init BestLinear" << std::endl;
			std::cout << options->collection_id <<"," << options->collection_percentage << ",BestLinear,";
			if(!init_computed){
				std::cout << "-1,"<<init_time <<",-1,"<< init_time << ",0" << std::endl;
			}
			else{
				gettimeofday(&tv1, NULL);
				local_median = computeRecursiveLinearMedian(completeset,ed,mcf,AllPairMappings,distances,mappingsFromMedian,mappingsToMedian, distancesToMedian,2,1,0,options->time_limit-init_time,0);
				gettimeofday(&tv2, NULL);
				//std::cout<<"median_computed"<< std::endl;
				double rec_init_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
				gettimeofday(&tv1, NULL);
				finalSOD=computeSOD(completeset,final_ed,local_median,options->time_limit -(rec_init_time+init_time));
				gettimeofday(&tv2, NULL);
				//std::cout<<"median_refined" << std::endl;
				double refinement_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
				status=3;
				if(rec_init_time+init_time>options->time_limit){
					status=1;
				}
				else if(rec_init_time+init_time+refinement_time>options->time_limit){
					status=2;
				}
				std::cout << init_time+rec_init_time << "," << finalSOD <<","<<refinement_time << "," << status << std::endl;
			}

			//BestLinear
			delete init_median;
			delete local_median;
			//std::cout << "init Besttriangular" << std::endl;
			std::cout << options->collection_id <<"," << options->collection_percentage << ",BestTriangular,";
			if(!init_computed){
				std::cout << "-1,"<<init_time <<",-1,"<< init_time << ",0" << std::endl;
			}
			else{
				gettimeofday(&tv1, NULL);
				local_median = computeRecursiveLinearMedian(completeset,ed,mcf,AllPairMappings,distances,mappingsFromMedian,mappingsToMedian, distancesToMedian,3,1,0,options->time_limit-init_time,0);
				gettimeofday(&tv2, NULL);
				//std::cout<<"median_computed"<< std::endl;
				double rec_init_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
				gettimeofday(&tv1, NULL);
				finalSOD=computeSOD(completeset,final_ed,local_median,options->time_limit -(rec_init_time+init_time));
				gettimeofday(&tv2, NULL);
				//std::cout<<"median_refined" << std::endl;
				double refinement_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
				status=3;
				if(rec_init_time+init_time>options->time_limit){
					status=1;
				}
				else if(rec_init_time+init_time+refinement_time>options->time_limit){
					status=2;
				}
				std::cout << init_time+rec_init_time << "," << finalSOD <<","<<refinement_time << "," << status << std::endl;
			}
			delete i_ed;
			for (int j = 0;j<M;j++){
				delete [] mappingsFromMedian[j];
				delete [] mappingsToMedian[j];
			}
			delete [] mappingsFromMedian;
			delete [] mappingsToMedian;

		}

		delete ed;
		delete cf;
		delete mcf;
		delete options;
		delete algoIPFP;

	} // end IBD
	else {
		// LETTER ------------------------------------------------------------------------------
		if (options->letter) {
			LetterDistanceCost *cf = new LetterDistanceCost(0.9,1.7,0.75,options->power_cost);
			LetterMedianLabel * mcf =  new LetterMedianLabel(cf);
			//std::cout << "median labels set\n";

			// IPFP used as a refinement method

			IPFPGraphEditDistance<CMUPoint,double> * algoIPFP = new IPFPGraphEditDistance<CMUPoint,double>(cf);
			algoIPFP->recenterInit();
			RandomMappingsGED<CMUPoint,double> *final_init = new RandomMappingsGED<CMUPoint,double>();
			GraphEditDistance<CMUPoint,double>* final_ed = new MultistartRefinementGraphEditDistance<CMUPoint,double>(cf, final_init, 10, algoIPFP,5,options->adap);
			GraphEditDistance<CMUPoint,double>* ed;
			if(options->method == "lsape_multi_bunke")
				ed = new BipartiteGraphEditDistanceMulti<CMUPoint,double>(cf, options->nep);
			else{
				RandomMappingsGED<CMUPoint,double> *init = new RandomMappingsGED<CMUPoint,double>();
				ed = new MultistartRefinementGraphEditDistance<CMUPoint,double>(cf, init, options->nep, algoIPFP,options->nrep,options->adap);
			}


			int dataset_reduction =0;
			int dataset_max_reduction=std::numeric_limits<int>::max();

			for(int k=0; k<options->num_repetitions_xp;k++){
				double init_time = 0.0;
				LetterDataset * completeset = new LetterDataset((options->dataset_file+"-"+options->collection_percentage+"-"+options->collection_id+".ds").c_str());
				int M = completeset->size();
				std::vector<LetterDataset*> complete_datasets_by_class;
				std::vector<LetterDataset*> datasets_by_class;
				LetterDataset* initial_dataset = new LetterDataset();
				initial_dataset->add((*completeset)[0],(*completeset)(0));
				complete_datasets_by_class.push_back(initial_dataset);
				//separation
				datasets_by_class.push_back(completeset);
				int P=datasets_by_class.size(); //P now equal to the number of classes in dataset_by_class
				for (int class_index = 0;class_index<P;class_index++){
					LetterDataset* test_set = datasets_by_class[class_index];
					int N=test_set->size();
					int prop_init_random_median = 10;
					bool init_distances = true;
					struct timeval  tv1, tv2;
#ifdef PRINT_TIMES
					gettimeofday(&tv1, NULL);
#endif

					for (int init_method = 0; init_method<1; init_method++){
						init_distances = true;
						GraphEditDistance<CMUPoint,double>* i_ed = nullptr;
						if(init_method == 0)
							i_ed = new BipartiteGraphEditDistanceMulti<CMUPoint,double>(cf, options->nep);
						else if (init_method == 1){
							RandomMappingsGED<CMUPoint,double> *init = new RandomMappingsGED<CMUPoint,double>();
							i_ed = new MultistartRefinementGraphEditDistance<CMUPoint,double>(cf, init, options->nep, algoIPFP,options->nrep,options->adap);
						}
						else{
							init_distances=false;
						}
						double rec_init_time;
						Graph<CMUPoint,double>* local_median = nullptr;
						Graph<CMUPoint,double>* init_median = nullptr;
						double * distancesToMedian = new double [N];
						int * * mappingsFromMedian = new int*[N];
						int * * mappingsToMedian = new int*[N];
						double * p_SOD= new double;
						double finalSOD;
						int * p_num_it = new int;
						//std::cout << std::endl <<"----- Method Based on Block gradient descent with set median init------"<< std::endl << std::endl;
						if (init_distances){
							// init_distances
							int status;
							int ** AllPairMappings =  new int*[N*N];
							double * distances = new double[N*N];
#ifdef PRINT_TIMES
							gettimeofday(&tv1, NULL);
#endif
							bool init_computed = computeMappingsAndDistances(test_set,i_ed,AllPairMappings,distances,options->time_limit);
#ifdef PRINT_TIMES
							gettimeofday(&tv2, NULL);
							init_time=((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
#endif
							//std::cout << "distances initialized in " << init_time << "seconds" << std::end

							std::cout << "Id_collection,percentage,algo,SOD,time,refined_SOD,refinement_time,status" << std::endl;
							//Linear
							std::cout << options->collection_id <<"," << options->collection_percentage << ",linear,";
							if(!init_computed){
								std::cout << "-1,"<<init_time <<",-1,"<< init_time << ",0" << std::endl;
							}
							else{
								gettimeofday(&tv1, NULL);
								local_median = computeRecursiveLinearMedian(test_set,ed,mcf,AllPairMappings,distances,mappingsFromMedian,mappingsToMedian, distancesToMedian,2,1,0,options->time_limit-init_time,1);
								gettimeofday(&tv2, NULL);
								//std::cout<<"median_computed"<< std::endl;
								double rec_init_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
								gettimeofday(&tv1, NULL);
								finalSOD=computeSOD(test_set,final_ed,local_median,options->time_limit -(rec_init_time+init_time));
								gettimeofday(&tv2, NULL);
								//std::cout<<"median_refined" << std::endl;
								double refinement_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
								status=3;
								if(rec_init_time+init_time>options->time_limit){
									status=1;
								}
								else if(rec_init_time+init_time+refinement_time>options->time_limit){
									status=2;
									finalSOD=-1;
								}
								std::cout << init_time+rec_init_time << "," << finalSOD <<","<<refinement_time << "," << status << std::endl;
							}
							//std::cout << "descent converged !! :)" << std::endl;
							delete init_median;
							delete local_median;
							//Triangular
							std::cout << options->collection_id <<"," << options->collection_percentage << ",Triangular,";
							if(!init_computed){
								std::cout << "-1,"<<init_time <<",-1,"<< init_time << ",0" << std::endl;
							}
							else{
								gettimeofday(&tv1, NULL);
								local_median = computeRecursiveLinearMedian(test_set,ed,mcf,AllPairMappings,distances,mappingsFromMedian,mappingsToMedian, distancesToMedian,3,1,0,options->time_limit-init_time,1);
								gettimeofday(&tv2, NULL);
								//std::cout<<"median_computed"<< std::endl;
								double rec_init_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
								gettimeofday(&tv1, NULL);
								finalSOD=computeSOD(test_set,final_ed,local_median,options->time_limit -(rec_init_time+init_time));
								gettimeofday(&tv2, NULL);
								//std::cout<<"median_refined" << std::endl;
								double refinement_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
								status=3;
								if(rec_init_time+init_time>options->time_limit){
									status=1;
								}
								else if(rec_init_time+init_time+refinement_time>options->time_limit){
									status=2;
								}
								std::cout << init_time+rec_init_time << "," << finalSOD <<","<<refinement_time << "," << status << std::endl;
							}
							//std::cout << "descent converged !! :)" << std::endl;
							delete init_median;
							delete local_median;
							//std::cout << "init BestLinear" << std::endl;
							std::cout << options->collection_id <<"," << options->collection_percentage << ",BestLinear,";
							if(!init_computed){
								std::cout << "-1,"<<init_time <<",-1,"<< init_time << ",0" << std::endl;
							}
							else{
								gettimeofday(&tv1, NULL);
								local_median = computeRecursiveLinearMedian(test_set,ed,mcf,AllPairMappings,distances,mappingsFromMedian,mappingsToMedian, distancesToMedian,2,1,0,options->time_limit-init_time,0);
								gettimeofday(&tv2, NULL);
								//std::cout<<"median_computed"<< std::endl;
								double rec_init_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
								gettimeofday(&tv1, NULL);
								finalSOD=computeSOD(test_set,final_ed,local_median,options->time_limit -(rec_init_time+init_time));
								gettimeofday(&tv2, NULL);
								//std::cout<<"median_refined" << std::endl;
								double refinement_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
								status=3;
								if(rec_init_time+init_time>options->time_limit){
									status=1;
								}
								else if(rec_init_time+init_time+refinement_time>options->time_limit){
									status=2;
								}
								std::cout << init_time+rec_init_time << "," << finalSOD <<","<<refinement_time << "," << status << std::endl;
							}

							//BestLinear
							delete init_median;
							delete local_median;
							//std::cout << "init Besttriangular" << std::endl;
							std::cout << options->collection_id <<"," << options->collection_percentage << ",BestTriangular,";
							if(!init_computed){
								std::cout << "-1,"<<init_time <<",-1,"<< init_time << ",0" << std::endl;
							}
							else{
								gettimeofday(&tv1, NULL);
								local_median = computeRecursiveLinearMedian(test_set,ed,mcf,AllPairMappings,distances,mappingsFromMedian,mappingsToMedian, distancesToMedian,3,1,0,options->time_limit-init_time,0);
								gettimeofday(&tv2, NULL);
								//std::cout<<"median_computed"<< std::endl;
								double rec_init_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
								gettimeofday(&tv1, NULL);
								finalSOD=computeSOD(test_set,final_ed,local_median,options->time_limit -(rec_init_time+init_time));
								gettimeofday(&tv2, NULL);
								//std::cout<<"median_refined" << std::endl;
								double refinement_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
								status=3;
								if(rec_init_time+init_time>options->time_limit){
									status=1;
								}
								else if(rec_init_time+init_time+refinement_time>options->time_limit){
									status=2;
								}
								std::cout << init_time+rec_init_time << "," << finalSOD <<","<<refinement_time << "," << status << std::endl;
							}
							delete i_ed;
							for (int j = 0;j<N;j++){
								delete [] mappingsFromMedian[j];
								delete [] mappingsToMedian[j];
							}
							delete [] mappingsFromMedian;
							delete [] mappingsToMedian;
						}
						else{
							local_median = computeMedianGraph(test_set,ed,mcf,0,final_ed); // tests min, mean, and max order
							local_median = computeMedianGraph(test_set,ed,mcf,40,final_ed);
						}
					}

				}

			}
			delete ed;
			delete cf;
			delete mcf;
			delete options;
			delete algoIPFP;
		} // end LETTER
		// SYMBOLIC ------------------------------------------------------------------------------
		else {

			ConstantEditDistanceCost* cf = new ConstantEditDistanceCost(options->cns,options->cni, options->cnd,
					options->ces,options->cei, options->ced);
			ConstantMedianLabel* mcf = new ConstantMedianLabel(cf);
			//std::cout << "median labels set\n";

			// IPFP used as a refinement method

			IPFPGraphEditDistance<int,int> * algoIPFP = new IPFPGraphEditDistance<int,int>(cf);
			algoIPFP->recenterInit();
			RandomMappingsGED<int,int> *final_init = new RandomMappingsGED<int,int>();
			//GraphEditDistance<int,int>* final_ed = new MultistartRefinementGraphEditDistance<int,int>(cf, final_init, 40, 10,options->adap);
			GraphEditDistance<int,int>* final_ed = new MultistartRefinementGraphEditDistance<int,int>(cf, final_init, 10, algoIPFP,5,options->adap); // fast final GED fpr testing
			GraphEditDistance<int,int>* ed;
			if(options->method == "lsape_multi_bunke")
				ed = new BipartiteGraphEditDistanceMulti<int,int>(cf, options->nep);
			else{
				RandomMappingsGED<int,int> *init = new RandomMappingsGED<int,int>();
				ed = new MultistartRefinementGraphEditDistance<int,int>(cf, init, options->nep, algoIPFP,options->nrep,options->adap);
			}


			int dataset_reduction =0;
			int dataset_max_reduction=std::numeric_limits<int>::max();

			for(int k=0; k<options->num_repetitions_xp;k++){
				double init_time = 0.0;
				ChemicalDataset<double> * completeset = new ChemicalDataset<double>((options->dataset_file+"-"+options->collection_percentage+"-"+options->collection_id+".ds").c_str());
				int M = completeset->size();
				std::vector<ChemicalDataset<double>*> complete_datasets_by_class;
				std::vector<ChemicalDataset<double>*> datasets_by_class;
				ChemicalDataset<double>* initial_dataset = new ChemicalDataset<double>();
				initial_dataset->add((*completeset)[0],(*completeset)(0));
				complete_datasets_by_class.push_back(initial_dataset);
				//separation
				datasets_by_class.push_back(completeset);

				int P=datasets_by_class.size(); //P now equal to the number of classes in dataset_by_class
				for (int class_index = 0;class_index<P;class_index++){
					std::vector<int> set_sizes;
					if(options->all_train_set_sizes){
						for (int i = 0;i<19;i++){
							if (i<10){
								set_sizes.push_back((i+1)*10);
							}
							else
								set_sizes.push_back((i-8)*100);

						}

					}
					else
						set_sizes.push_back(datasets_by_class[class_index]->size());

					for (int size_index = 0; size_index < set_sizes.size();size_index++){
						int current_size;
						ChemicalDataset<double>* test_set = nullptr;
						if (set_sizes[size_index]>datasets_by_class[class_index]->size()){
							current_size = datasets_by_class[class_index]->size();
						}
						else
							current_size = set_sizes[size_index];

						if (current_size == datasets_by_class[class_index]->size()){
							test_set = datasets_by_class[class_index];
						}
						else{
							test_set = new ChemicalDataset<double>;
							for (int i =0;i<current_size;i++){
								test_set->add(datasets_by_class[class_index]->getGraph(i),datasets_by_class[class_index]->getProperty(i));
							}
							//std::cout << " taille dataset " << test_set->size() << std::endl;
						}
						int N=test_set->size();
						int prop_init_random_median = 10;
						bool init_distances = true;
						struct timeval  tv1, tv2;
#ifdef PRINT_TIMES
						gettimeofday(&tv1, NULL);
#endif

						for (int init_method = 0; init_method<1; init_method++){ // only bipartite initialization
							init_distances = true;
							GraphEditDistance<int,int>* i_ed = nullptr;
							if(init_method == 0)
								i_ed = new BipartiteGraphEditDistanceMulti<int,int>(cf, options->nep);
							else if (init_method == 1){
								RandomMappingsGED<int,int> *init = new RandomMappingsGED<int,int>();
								i_ed = new MultistartRefinementGraphEditDistance<int,int>(cf, init, options->nep, algoIPFP,options->nrep,options->adap);
							}
							else{
								init_distances=false;
							}
							double rec_init_time;
							Graph<int,int>* local_median = nullptr;
							Graph<int,int>* init_median = nullptr;
							double * distancesToMedian = new double [N];
							int * * mappingsFromMedian = new int*[N];
							int * * mappingsToMedian = new int*[N];
							double finalSOD;
							double * p_SOD= new double;
							int * p_num_it = new int;
							//std::cout << std::endl <<"----- Method Based on Block gradient descent with set median init------"<< std::endl << std::endl;
							if (init_distances){
								// init_distances
								int status;
								int ** AllPairMappings =  new int*[N*N];
								double * distances = new double[N*N];
#ifdef PRINT_TIMES
								gettimeofday(&tv1, NULL);
#endif
								bool init_computed = computeMappingsAndDistances(test_set,i_ed,AllPairMappings,distances,options->time_limit);
#ifdef PRINT_TIMES
								gettimeofday(&tv2, NULL);
								init_time=((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
#endif
								//std::cout << "distances initialized in " << init_time << "seconds" << std::end

								std::cout << "Id_collection,percentage,algo,SOD,time,refined_SOD,refinement_time,status" << std::endl;
								//Linear
								std::cout << options->collection_id <<"," << options->collection_percentage << ",linear,";
								if(!init_computed){
									std::cout << "-1,"<<init_time <<",-1,"<< init_time << ",0" << std::endl;
								}
								else{
									gettimeofday(&tv1, NULL);
									local_median = computeRecursiveLinearMedian(test_set,ed,mcf,AllPairMappings,distances,mappingsFromMedian,mappingsToMedian, distancesToMedian,2,1,0,options->time_limit-init_time,1);
									gettimeofday(&tv2, NULL);
									//std::cout<<"median_computed"<< std::endl;
									double rec_init_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
									gettimeofday(&tv1, NULL);
									finalSOD=computeSOD(test_set,final_ed,local_median,options->time_limit -(rec_init_time+init_time));
									gettimeofday(&tv2, NULL);
									//std::cout<<"median_refined" << std::endl;
									double refinement_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
									status=3;
									if(rec_init_time+init_time>options->time_limit){
										status=1;
									}
									else if(rec_init_time+init_time+refinement_time>options->time_limit){
										status=2;
									}
									std::cout << init_time+rec_init_time << "," << finalSOD <<","<<refinement_time << "," << status << std::endl;
								}
								//std::cout << "descent converged !! :)" << std::endl;
								delete init_median;
								delete local_median;
								//Triangular
								std::cout << options->collection_id <<"," << options->collection_percentage << ",Triangular,";
								if(!init_computed){
									std::cout << "-1,"<<init_time <<",-1,"<< init_time << ",0" << std::endl;
								}
								else{
									gettimeofday(&tv1, NULL);
									local_median = computeRecursiveLinearMedian(test_set,ed,mcf,AllPairMappings,distances,mappingsFromMedian,mappingsToMedian, distancesToMedian,3,1,0,options->time_limit-init_time,1);
									gettimeofday(&tv2, NULL);
									//std::cout<<"median_computed"<< std::endl;
									double rec_init_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
									gettimeofday(&tv1, NULL);
									finalSOD=computeSOD(test_set,final_ed,local_median,options->time_limit -(rec_init_time+init_time));
									gettimeofday(&tv2, NULL);
									//std::cout<<"median_refined" << std::endl;
									double refinement_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
									status=3;
									if(rec_init_time+init_time>options->time_limit){
										status=1;
									}
									else if(rec_init_time+init_time+refinement_time>options->time_limit){
										status=2;
									}
									std::cout << init_time+rec_init_time << "," << finalSOD <<","<<refinement_time << "," << status << std::endl;
								}
								//std::cout << "descent converged !! :)" << std::endl;
								delete init_median;
								delete local_median;
								//std::cout << "init BestLinear" << std::endl;
								std::cout << options->collection_id <<"," << options->collection_percentage << ",BestLinear,";
								if(!init_computed){
									std::cout << "-1,"<<init_time <<",-1,"<< init_time << ",0" << std::endl;
								}
								else{
									gettimeofday(&tv1, NULL);
									local_median = computeRecursiveLinearMedian(test_set,ed,mcf,AllPairMappings,distances,mappingsFromMedian,mappingsToMedian, distancesToMedian,2,1,0,options->time_limit-init_time,0);
									gettimeofday(&tv2, NULL);
									//std::cout<<"median_computed"<< std::endl;
									double rec_init_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
									gettimeofday(&tv1, NULL);
									finalSOD=computeSOD(test_set,final_ed,local_median,options->time_limit -(rec_init_time+init_time));
									gettimeofday(&tv2, NULL);
									//std::cout<<"median_refined" << std::endl;
									double refinement_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
									status=3;
									if(rec_init_time+init_time>options->time_limit){
										status=1;
									}
									else if(rec_init_time+init_time+refinement_time>options->time_limit){
										status=2;
									}
									std::cout << init_time+rec_init_time << "," << finalSOD <<","<<refinement_time << "," << status << std::endl;
								}

								//BestLinear
								delete init_median;
								delete local_median;
								//std::cout << "init Besttriangular" << std::endl;
								std::cout << options->collection_id <<"," << options->collection_percentage << ",BestTriangular,";
								if(!init_computed){
									std::cout << "-1,"<<init_time <<",-1,"<< init_time << ",0" << std::endl;
								}
								else{
									gettimeofday(&tv1, NULL);
									local_median = computeRecursiveLinearMedian(test_set,ed,mcf,AllPairMappings,distances,mappingsFromMedian,mappingsToMedian, distancesToMedian,3,1,0,options->time_limit-init_time,0);
									gettimeofday(&tv2, NULL);
									//std::cout<<"median_computed"<< std::endl;
									double rec_init_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
									gettimeofday(&tv1, NULL);
									finalSOD=computeSOD(test_set,final_ed,local_median,options->time_limit -(rec_init_time+init_time));
									gettimeofday(&tv2, NULL);
									//std::cout<<"median_refined" << std::endl;
									double refinement_time = ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));
									status=3;
									if(rec_init_time+init_time>options->time_limit){
										status=1;
									}
									else if(rec_init_time+init_time+refinement_time>options->time_limit){
										status=2;
									}
									std::cout << init_time+rec_init_time << "," << finalSOD <<","<<refinement_time << "," << status << std::endl;
								}
								delete i_ed;
								for (int j = 0;j<N;j++){
									delete [] mappingsFromMedian[j];
									delete [] mappingsToMedian[j];
								}
								delete [] mappingsFromMedian;
								delete [] mappingsToMedian;
							}
							else{
								local_median = computeMedianGraph(test_set,ed,mcf,0,final_ed); // tests min, mean, and max order
								local_median = computeMedianGraph(test_set,ed,mcf,40,final_ed);
							}
						}

					}
				}
			}
			delete ed;
			delete cf;
			delete mcf;
			delete options;
			delete algoIPFP;
		} // end SYMBOLIC
	} // end else for IBD
	return 0;
}
