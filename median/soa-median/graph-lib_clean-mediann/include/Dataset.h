/**
 * @file Dataset.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>>
 * @version     0.0.1 - Mon Feb  6 2017
 *
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *
 * @brief The Dataset class implements a collection of templated graphs together with templated properties.
 */

#ifndef __DATASET_H__
#define __DATASET_H__
#define TI_USE_STL
#include <vector>
#include <string>
#include <libgen.h>
#include <map>
#include <sstream>
#include "graph.h"
#include "SymbolicGraph.h"
#include "GraphEditDistance.h"

template <class NodeAttribute, class EdgeAttribute, class PropertyType>
class Dataset
{
private:
  /* List of properties associated to graphs
   */
  std::vector<PropertyType> properties;
  /* List of graphs composing the dataset
   */
  std::vector<Graph<NodeAttribute, EdgeAttribute> * > graphs;

  // TODO: void loadCXL(const char * filename);
public:
  // TODO :  Dataset(const char * filename){};
  /**
   * Accessors to graphs
   * @param id the identifier of graph
   * @return the desired graph
   * XXX: control id param.
   */
  Graph<NodeAttribute, EdgeAttribute> * getGraph(unsigned int id) const {return graphs[id];};
  Graph<NodeAttribute, EdgeAttribute> * operator[](unsigned int id) const {return graphs[id];};

  /**
   * Accessors to properties
   * @param id the identifier of graph associated to property
   * @return the property associated to graph id
   * XXX: control id param.
   */
  PropertyType getProperty(unsigned int id) const {return properties[id];};
  PropertyType operator()(unsigned int id) const {return properties[id];};

  PropertyType getMostFrequentProperty() const {

    int maxFrequence = 0;
    int Medianlabel = 0;
    int N= size();
    std::map<PropertyType,int> m;
    typename std::map<PropertyType,int>::iterator it;
    PropertyType label;
    for (int i=0;i<N;i++){

        label =  getProperty(i)  ;
  	//std::cout << "label sommet " << node1 << " dans graphe " << i << " est " << label << std::endl;
        it = m.find(label);
        if (it == m.end()) // label non présent
                m[label] = 1;       // ajout d'un nouveau label avec une fréquence de 1
        else it->second++;
                 // incrémentation de la fréquence du label
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




  /**
   * Add a pair of graph and its attribute to the dataset
   * @param g the graph to add
   * @param y the property associated to g
   * @return the id of added graph
  */
  int add(Graph<NodeAttribute, EdgeAttribute> * g, PropertyType y){graphs.push_back(g);properties.push_back(y);return graphs.size();}

  /**
   * Returns and erase the last graph of the dataset
   * @return the last graph of the dataset, NULL if it is empty
   */
  Graph<NodeAttribute, EdgeAttribute> * pop_back(){
    if (graphs.empty()) return NULL;
    Graph<NodeAttribute, EdgeAttribute> * g = graphs.back();
    graphs.pop_back();
    return g;
  }

   /**
   * empty the dataset but does not delete the graphs
   */

   void clear_not_delete(){
    graphs.clear();
    properties.clear();
  }


  /**
   * Clear all graphs in the dataset
   */

  void clear(){
    for (unsigned int i=0; i<graphs.size(); i++)
      delete graphs[i];
    graphs.clear();
  }

  /**
   * Return the number of graphs included within the dataset
   * @return the size of dataset
   */
  int size() const {return graphs.size();};

  /**
   * Compute the graph edit distance according to a given algorithm between each pair of graphs included within dataset.
   * @param ed the <code>GraphEditDistance</code> used to compute the ged
   * @param quiet if TRUE, prints the pair of graphs currently processed while execution.
   * @return the distance matrix of size N*N computed according to ed
   */
  double * computeGraphEditDistance(GraphEditDistance<NodeAttribute,EdgeAttribute> * ed, bool quiet = true) const;

  /**
   * Apply shuffleization procedure on all graphs composing the dataset.
   */
  void shuffleize();




  /**
   * An empty dataset
   */
  Dataset(){}

  ~Dataset(){
    for (unsigned int i=0; i<graphs.size(); i++)
      delete graphs[i];
  }
};


template<class NodeAttribute,class EdgeAttribute, class PropertyType>
void Dataset<NodeAttribute,EdgeAttribute,PropertyType>::shuffleize(){
//std::mt19937 generator (seed);
int N = size();
  for (int i=0;i<N;i++)
    (*this)[i]->shuffleize();
}



template<class NodeAttribute,class EdgeAttribute, class PropertyType>
double * Dataset<NodeAttribute,EdgeAttribute,PropertyType>::computeGraphEditDistance(GraphEditDistance<NodeAttribute,EdgeAttribute> * ed, bool quiet) const{
  int N = size();
  double * distances = new double[size()*size()];
  for (int i=0;i<N;i++)
    for (int j=0;j<N;j++)
      {
      distances[sub2ind(i,j,N)] = (*ed)(graphs[i],graphs[j]);
      if(!quiet)
	std::cout << i << "," << j << '\r';
      }
  if(!quiet)
    std::cout << std::endl;
  return distances;
}

/*@brief The class ChemicalDataset implements a collection of <code>SymbolicGraph</code>
 *
 */
template<class PropertyType>
class ChemicalDataset: public Dataset<int,int,PropertyType>
{
private:
  void loadDS(const char * filename);
public:
  ChemicalDataset(const char * filename);
  ChemicalDataset(): Dataset<int,int,PropertyType>(){}
};

template<class PropertyType>
void ChemicalDataset<PropertyType>::loadDS(const char* filename){

  std::ifstream f_tmp (filename);
  char * unconst_filename = new char[strlen(filename)+1];
  unconst_filename = strcpy(unconst_filename, filename);
  char * path = dirname(unconst_filename);
  if (f_tmp.is_open()){
    std::string s;
    while (getline(f_tmp, s))
      if (s[0] != '#'){
	std::string path_ctfile(path);
	path_ctfile += std::string("/");
	std::istringstream liness(s);
	std::string ctfile;
	liness >> ctfile;
	PropertyType y;
	liness >> y;
	std::string full_ctfile = path_ctfile;
	full_ctfile += ctfile;
	SymbolicGraph * g = new SymbolicGraph(full_ctfile.c_str());
	this->add(g,y);
      }
  }
  f_tmp.close();

  delete[] unconst_filename;

}
template<class PropertyType>
ChemicalDataset<PropertyType>::ChemicalDataset(const char * filename){
  const char * ext = strrchr(filename,'.');
  if (strcmp(ext,".ds") == 0){
    loadDS(filename);
  }
}


#endif // __DATASET_H__
