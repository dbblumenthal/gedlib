======= TERMS OF USAGE ======= 

The IAM-Graph DB is publicly accessible and freely available for non-commercial research purposes. If you are publishing scientific work based on the IAM-Graph DB, we request you to include the following reference to our database:

@unpublished{riesen08iam,
	Author = {Riesen, K. and Bunke, H.,
	Note = {accepted for publication in SSPR 2008,
	Title = {{IAM Graph Database Repository for Graph Based Pattern Recognition and Machine Learning


=======  DATA SET ======= 

Each of the data sets available on this repository is divided into three disjoint subsets, which can be used for training, validation, and testing novel learning algorithms (train.cxl, valid.cxl, test.cxl). The protein data set consists of graphs representing proteins originally used in[borgwardt05protein]. The graphs are constructed from the Protein Data Bank  and labeled with their corresponding enzyme class labels from the BRENDA enzyme database. The proteins database consists of six classes (EC 1, EC 2, EC 3, EC 4, EC 5, EC 6), which represent proteins out of the six enzyme commission top level hierarchy (EC classes). The proteins are converted into graphs by representing the secondary structure elements of a protein with nodes and edges of an attributed graph. Nodes are labeled with their type (helix, sheet, or loop) and their amino acid sequence (e.g. TFKEVVRLT). Every node is connected with an edge to its three nearest neighbors in space. Edges are labeled with their type and the distance they represent in angstroms. In sample.eps six images of proteins of all six classes are given. There are 600 proteins totally, 100 per class. We use a training, validation and test set of equal size (200). The classification task on this data set consists in predicting the enzyme class membership. We achieve a classification rate of 65.5% on the test set. 

=======  REFERENCES ======= 

This data set is employed in the following publications (this list does not claim to be exhaustive, of course):

@phdthesis{borgwardt07graph,
	Author = {Borgwardt, K.},
	School = {Ludwig-Maximilians-University Munich},
	Title = {Graph Kernels},
	Year = {2007}}

@article{borgwardt05protein,
	Author = {Borgwardt, K. and Ong, C. and Sch{\"o}nauer, S. and Vishwanathan, S. and Smola, A. and Kriegel, H.-P.},
	Journal = {Bioinformatics},
	Keywords = {Graph Kernel, Protein Prediction},
	Number = {1},
	Pages = {47--56},
	Title = {Protein Function Prediction via Graph Kernels},
	Volume = {21},
	Year = {2005}}


=======  CONTACT INFORMATION ======= 

If you have any question concerning this data set, do not hesitate to contact me: riesen@iam.unibe.ch

