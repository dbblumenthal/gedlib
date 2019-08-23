CONTENTS OF THIS FILE
---------------------
   
 * Reference
 * Authors
 * Database Content

Reference
---------------------
A Graph Database Repository and Performance Evaluation Metrics for Graph Edit Distance, submitted to GBR2015.


Authors
---------------------
Zeina Abu-Aisheh, 3rd year PhD student at François Rabelais University, Tours, France.
Romain Raveaux, Assistant Professor at François Rabelais University, Tours, France.
Jean-Yves Ramel, Professor at François Rabelais University, Tours, France.  


Database Content
---------------------
 - CMU-cost-function

   Two meta parameters are included tau_{vertex} and tau_{edge} where tau_{vertex} denotes a vertex deletion or insertion whereas tau_{edge}
   denotes an edge deletion or insertion. Both tau_{vertex} and tau_{edge} are non-negative parameters. A third meta parameter alpha is integrated to control whether
   the edit operation cost on the vertices or on the edges is more important. tau_{node}, tau_{edge} and alpha are set to 100000, 100000 and 0.5 respectively.
   The cost function is written in Java and is found in the CMU-GED folder.
 
 - CMU-low-level-info
   
  The pairwise comparisons of each of the aforementioned subsets is performed. For each comparison, the following information is added:
  1- The name of each pair of graphs d(g_1,g_2).
  2- The number of vertices of each pair of graphs.
  3- The number of edges of each pair of graphs.
  4- The name of a graph edit distance (GED) method P that succeeds at finding the best distance for d(g_1,g_2).
  5- A parameter used for the method BeamSearch GED which represents the size of the stack (1, 10 or 100 in the experiments).
  6- The best distance that was found by method P.
  7- A boolean value that tells whether the found solution is optimal or not.
  8- The classes of graphs g_1 and g_2.
  9- The matching sequence found by method P.
  10- For instance Node:3->4=0.0 means that substituting vertex 3 on graph g_1 with vertex 4 on graph g_2 costs 0.0. Moreover, Edge:1_<>2->1_<>5=4.0 means that substituting
      1_<>2 on graph g_1 with edge 1_<>5 on graph g_2 costs 4.0. On the other hand, Node:3->eps=4.0 means that deleting vertex 3 costs 4.0 while vertex:eps->4=6.0 means that
      inserting vertex 4 costs 6.0.	 
 
 - CMU-subsets

   In CMU there are no subsets. All graphs have the same number of vertices (i.e., 30 vertices).

 - Information about CMU

   The CMU model house sequence is a database made up of a series of images of a toy house that has been captured from different viewpoints. 111 images in total 
   are publicly available. 660 comparisons are carried out. A manual identication of corner features, or points, is done to represent vertices on each of the
   rotated images. Then, the Delaunay triangulation is applied on the corner-features in order to identify edges and finally transform the images into graphs. 
   Each graph has 30 vertices. Vertices are labelled with (x,y) coordinates while edges are labelled with the distance between vertices. Moreover, the ground 
   truth is attached with each pair of graphs.
   