/**
 * @file MultipleIPFPGraphEditDistance.h
 * @author Evariste <<evariste.daller@unicaen.fr>>
 * @version     Jun  9 2017
 *
 */

#ifndef __MULTISTARTMAPPINGREFINEMENT_H__
#define __MULTISTARTMAPPINGREFINEMENT_H__

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <sys/time.h>
#include <list>
#include "MappingRefinement.h"
#include "MappingGenerator.h"


/**
 * @brief A MappingRefinement method which uses a multistart approach
 * 
 *   The multistart approach consists in refine several initial mappings
 *   and keep the best one. To this end, the initializations are here lists of
 *   mappings
 */
template<class NodeAttribute, class EdgeAttribute>
class MultistartMappingRefinement
{

protected:

  MappingGenerator<NodeAttribute, EdgeAttribute> * initGen; //!< Generator of initializations 
  int k; //!< Number of initial mapping to generate from \ref initGen
  std::list<int*> refinedMappings; //!< The last set of refined mappings

public:

  /**
   * @brief Outputs the best mapping refined by \ref algorithm from initializations given by \ref initGen
   * 
   *  During the process, only the best mapping is kept in memory
   *  allowing the procedure to require less memory than \ref getBetterMappings.
   *
   * @param  algorithm   the refinement method
   * @param  g1          First graph
   * @param  g2          Second graph
   * @param  G1_to_G2    forward output mapping 
   * @param  G2_to_G1    reverse output mapping, useful for the graph edit distance
   * @see getBestMappingFromSet
   */
  virtual void getBestMapping( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                               Graph<NodeAttribute,EdgeAttribute> * g1,
                               Graph<NodeAttribute,EdgeAttribute> * g2,
                               int * G1_to_G2, int * G2_to_G1 );

  /**
   * @brief Returns the list of refined mappings from initializations given by \reg initGen
   *
   * @param  algorithm   the refinement method
   * @param  g1          First graph
   * @param  g2          Second graph
   * @see getBetterMappingsFromSet
   */
  virtual const std::list<int*>&
  getBetterMappings( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                     Graph<NodeAttribute,EdgeAttribute> * g1,
                     Graph<NodeAttribute,EdgeAttribute> * g2 );

  /**
   * @brief Returns the best mappings refined by \ref algorithm from the given \ref mappings
   */
  virtual void getBestMappingFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                                      Graph<NodeAttribute,EdgeAttribute> * g1,
                                      Graph<NodeAttribute,EdgeAttribute> * g2,
                                      int * G1_to_G2, int * G2_to_G1,
                                      std::list<int*>& mappings );


  /**
   * @brief Returns the list of refined mappings from the given \ref mappings
   *
   * @param  algorithm   the refinement method
   * @param  g1          First graph
   * @param  g2          Second graph
   * @param  mapping     a list of arrays representing initial mappings
   * @note   mappings are allocated on the heap and memory management is left to the user
   * @see getBetterMappings
   */
  virtual const std::list<int*>&
  getBetterMappingsFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                            Graph<NodeAttribute,EdgeAttribute> * g1,
                            Graph<NodeAttribute,EdgeAttribute> * g2,
                            std::list<int*>& mappings );


  MultistartMappingRefinement( MappingGenerator<NodeAttribute, EdgeAttribute> * gen,
                               int nSol
                             ):
    initGen(gen),
    k(nSol)
  {}

  virtual ~MultistartMappingRefinement(){}


};

//---


template<class NodeAttribute, class EdgeAttribute>
void MultistartMappingRefinement<NodeAttribute, EdgeAttribute>::
getBestMapping( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                Graph<NodeAttribute,EdgeAttribute> * g1,
                Graph<NodeAttribute,EdgeAttribute> * g2,
                int * G1_to_G2, int * G2_to_G1 )
{
  //Compute Mapping init
  struct timeval  tv1, tv2;
  gettimeofday(&tv1, NULL);

  std::list<int*> mappings = initGen->getMappings(g1, g2, k);
  gettimeofday(&tv2, NULL);

  this->getBestMappingFromSet(algorithm, g1, g2, G1_to_G2, G2_to_G1, mappings);
  
  // Memoy cleaning
  for (std::list<int*>::iterator it = mappings.begin(); it!=mappings.end(); it++)
    delete[] *it;
}



template<class NodeAttribute, class EdgeAttribute>
void MultistartMappingRefinement<NodeAttribute, EdgeAttribute>::
getBestMappingFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                Graph<NodeAttribute,EdgeAttribute> * g1,
                Graph<NodeAttribute,EdgeAttribute> * g2,
                int * G1_to_G2, int * G2_to_G1,
                std::list<int*>& mappings )
{
  struct timeval  tv1, tv2;
  int n = g1->Size();
  int m = g2->Size();

  typename std::list<int*>::const_iterator it;
  double cost = -1;
  double ncost;


  // Multithread
  #ifdef _OPENMP
    gettimeofday(&tv1, NULL);
    int** arrayMappings = new int*[mappings.size()];
    int* arrayCosts = new int[mappings.size()];
    int* arrayLocal_G1_to_G2 = new int[n * mappings.size()];
    int* arrayLocal_G2_to_G1 = NULL;
    if (G2_to_G1 != NULL)
      arrayLocal_G2_to_G1 = new int[m * mappings.size()];

    int i=0; for (it=mappings.begin(); it!=mappings.end(); it++){
      arrayMappings[i] = *it;
      i++;
    }

    //omp_set_dynamic(0);
    //omp_set_num_threads(4);
    #pragma omp parallel for schedule(dynamic) //private(tid, i, j, ncost, ipfpGed )
    for (unsigned int tid=0; tid<mappings.size(); tid++){
      int* lsapMapping = arrayMappings[tid];
      int* local_G1_to_G2 = &(arrayLocal_G1_to_G2[tid*n]);

      int* local_G2_to_G1 = NULL;
      if (G2_to_G1 != NULL)
        local_G2_to_G1 = &(arrayLocal_G2_to_G1[tid*m]);

  // Sequential
  #else
    int* local_G1_to_G2 = new int[n];
    int* local_G2_to_G1 = NULL;
    if (G2_to_G1 != NULL)
      local_G2_to_G1 = new int[m];

    double t_acc = 0; // accumulated time
    for (it=mappings.begin(); it!=mappings.end(); it++){
      gettimeofday(&tv1, NULL);
      int* lsapMapping = *it;
  #endif

    // Copy the mapping into the local array
    for (int i=0; i<n; i++)
      local_G1_to_G2[i] = lsapMapping[i];

    if (local_G2_to_G1 != NULL){
      for (int j=0; j<m; j++) local_G2_to_G1[j] = -1;
      for (int i=0; i<n; i++) 
        if (local_G1_to_G2[i] >= 0)
          local_G2_to_G1[local_G1_to_G2[i]] = i;
    }
    
    MappingRefinement<NodeAttribute, EdgeAttribute> * local_method;
    
    #ifdef _OPENMP
      local_method = algorithm->clone();
    #else
      local_method = algorithm;
    #endif
    
    local_method->getBetterMapping(g1, g2, local_G1_to_G2, local_G2_to_G1, true);
    ncost = local_method->mappingCost(g1, g2, local_G1_to_G2, local_G2_to_G1);


    // Multithread
    #ifdef _OPENMP
      // save the approx cost
      arrayCosts[tid] = ncost;
     // _distances_[tid] = ncost;

    // Sequential
    #else
      // if ncost is better : save the mapping and the cost
      if (cost > ncost || cost == -1){
        cost = ncost;
        for (int i=0; i<n; i++) G1_to_G2[i] = local_G1_to_G2[i];
        if (G2_to_G1 != NULL)
          for (int j=0; j<m; j++) G2_to_G1[j] = local_G2_to_G1[j];
      }
      gettimeofday(&tv2, NULL);
      t_acc += ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec));

    #endif

  } //end for

  // Multithread : Reduction
  #ifdef _OPENMP
    gettimeofday(&tv2, NULL);
    
    gettimeofday(&tv1, NULL);

    int i_optim;
    for (unsigned int i=0; i<mappings.size(); i++){
      if (cost > arrayCosts[i] || cost == -1){
         cost = arrayCosts[i];
         i_optim = i;
      }
    }
    for (int i=0; i<n; i++) G1_to_G2[i] = arrayLocal_G1_to_G2[i_optim*n + i];
    
    if (G2_to_G1 != NULL)
      for (int j=0; j<m; j++) G2_to_G1[j] = arrayLocal_G2_to_G1[i_optim*m + j];

    // To match the output format size in XPs
    //for (int i=mappings.size(); i<k; i++) _distances_[i] = 9999;

    gettimeofday(&tv2, NULL);
    //_xp_out_ <<  ((double)(tv2.tv_usec - tv1.tv_usec)/1000000 + (double)(tv2.tv_sec - tv1.tv_sec)) << ", ";
    //_xp_out_ << ((float)t) / CLOCKS_PER_SEC << ", ";

    delete[] arrayLocal_G1_to_G2;
    delete[] arrayCosts;
    delete[] arrayMappings;
    if (arrayLocal_G2_to_G1 != NULL)
      delete[] arrayLocal_G2_to_G1;

  // Sequential : deletes
  #else

    delete [] local_G1_to_G2;
    if (G2_to_G1 != NULL)
      delete [] local_G2_to_G1;

  #endif


}


template<class NodeAttribute, class EdgeAttribute>
const std::list<int*>& MultistartMappingRefinement<NodeAttribute, EdgeAttribute>::
getBetterMappings( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                   Graph<NodeAttribute,EdgeAttribute> * g1,
                   Graph<NodeAttribute,EdgeAttribute> * g2 )
{
  //Compute Mapping init
  struct timeval  tv1, tv2;
  gettimeofday(&tv1, NULL);

  std::list<int*> mappings = initGen->getMappings(g1, g2, k);
  gettimeofday(&tv2, NULL);

  return this->getBetterMappingsFromSet(algorithm, g1, g2, mappings);
}



template<class NodeAttribute, class EdgeAttribute>
const std::list<int*>& MultistartMappingRefinement<NodeAttribute, EdgeAttribute>::
getBetterMappingsFromSet( MappingRefinement<NodeAttribute, EdgeAttribute> * algorithm,
                          Graph<NodeAttribute,EdgeAttribute> * g1,
                          Graph<NodeAttribute,EdgeAttribute> * g2,
                          std::list<int*>& mappings )
{
  int n = g1->Size();

  typename std::list<int*>::const_iterator it;
  
  int* arrayLocal_G1_to_G2 = new int[n * mappings.size()];

  // Multithread
  #ifdef _OPENMP
    int** arrayMappings = new int*[mappings.size()];

    int i=0; for (it=mappings.begin(); it!=mappings.end(); it++){
      arrayMappings[i] = *it;
      i++;
    }

    //omp_set_dynamic(0);
    //omp_set_num_threads(4);
    // Set a dynamic scheduler : algorithms can have different execution time given different initial mappings
    #pragma omp parallel for schedule(dynamic) //private(tid, i, j, ncost, ipfpGed )
    for (unsigned int tid=0; tid<mappings.size(); tid++){
      int* lsapMapping = arrayMappings[tid];

  // Sequential
  #else
    unsigned int tid=0;
    for (it=mappings.begin(); it!=mappings.end(); it++){
      int* lsapMapping = *it;
  #endif
  
  
    int* local_G1_to_G2 = &(arrayLocal_G1_to_G2[tid*n]);

    // Copy the mapping into the local array
    for (int i=0; i<n; i++)
      local_G1_to_G2[i] = lsapMapping[i];
    
    MappingRefinement<NodeAttribute, EdgeAttribute> * local_method;
    
    #ifdef _OPENMP
      local_method = algorithm->clone();
    #else
      local_method = algorithm;
      tid++;
    #endif
    
    local_method->getBetterMapping(g1, g2, local_G1_to_G2, NULL, true);

  } //end for

  // Reduction - indexation of the list
    refinedMappings.clear();
    int _i_=0;
    for (it=mappings.begin(); it!=mappings.end(); it++){
      refinedMappings.push_back(arrayLocal_G1_to_G2 + _i_*n);
      _i_++;
    }

   #ifdef _OPENMP
    delete[] arrayMappings;
   #endif
   
   return refinedMappings;
}


#endif // __MULTIPLEIPFPGRAPHEDITDISTANCE_H__
