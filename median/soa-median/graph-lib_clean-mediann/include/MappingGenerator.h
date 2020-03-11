/**
 * @file MappingGenerator.h
 * @author Evariste <<evariste.daller@unicaen.fr>>
 * @version  Sep  25  2017 
 */

#ifndef __MAPPING_GENERATOR_H
#define __MAPPING_GENERATOR_H

#include <list>


/**
 * @brief A mappings (or assignements) generator between two graphs
 */
template<class NodeAttribute, class EdgeAttribute>
  class MappingGenerator {

 public:

  /**
   * @brief How to generate the mappings
   */
  virtual std::list<int*> getMappings( Graph<NodeAttribute, EdgeAttribute>* g1, Graph<NodeAttribute, EdgeAttribute>* g2,
				       int k ) = 0;



  virtual MappingGenerator<NodeAttribute, EdgeAttribute>* clone() const = 0;
};


#endif
