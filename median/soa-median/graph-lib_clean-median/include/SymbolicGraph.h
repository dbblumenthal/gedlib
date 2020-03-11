/**
 * @file SymbolicGraph.h
 * @author Benoit <<benoit.gauzere@insa-rouen.fr>>
 * @version     0.0.1 - Sat Feb  4 2017
 *
 * @todo see TODO marks
 * @bug see XXX marks
 *
 * @brief The class Symbolic Graph implements a special case of graph where both edge and node attributes are symbolic.
 */

#ifndef __SYMBOLICGRAPH_H__
#define __SYMBOLICGRAPH_H__
#include "graph.h"
class SymbolicGraph: public Graph<int,int>
{
private:
  /* Function to read Edge Attribute from a GXL file, according to K. Riesen's datasets
   * @return the label read from file from elem
   */
  static int readChemicalEdgeLabel(TiXmlElement *elem);
  /* Function to read Node Attribute from a GXL file, according to K. Riesen's datasets
   * @return the label read from file from elem
   */
  static int readChemicalNodeLabel(TiXmlElement *elem);

  /* Read Edge Attributes from a Graphml file, for datasets ENZYME, and D&D
   * @return 1 everywhere
   */
  static int readGraphmlEdgeLabel(TiXmlElement *elem);

  /* Read Node Attributes from a Graphml file, for datasets ENZYME and D&D
   * @return the label read
   */
  static int readGraphmlNodeLabel(TiXmlElement *elem);

public:
  /* Constructor to fill a Symbolic graph from a ct file (ChemDraw Connection Table format)
   * @param filename path to a ct file.
   */
  SymbolicGraph(const char * filename);

  /* Constructor to fill a Symbolic graph from an adjacency matrix encoded as a n * n int array. Diagonals elements embed node's labels
   * @param am adjacency matrix encoded as a int array
   * @param nb_nodes specify the graph size
   * @param directed TRUE if graph encoded by am is directed, FALSE otherwise
   */
  SymbolicGraph(int * am, int nb_nodes, bool directed);


  /* Return a n*n int array encoding the adjacency matrix corresponding to current graph.
   * @return a pointer to adjacency matrix
   */
  int * getLabeledAdjacencyMatrix();

};


/* Write the graph into the given output stream, in the CT format
 */
bool writeCTfile(Graph<int,int>& graph, std::ofstream& output);


bool writeCTfile_std(Graph<int,int>& graph);

#endif // __SYMBOLICGRAPH_H__

