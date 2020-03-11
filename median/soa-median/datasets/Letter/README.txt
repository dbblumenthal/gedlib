======= TERMS OF USAGE ======= 

The IAM-Graph DB is publicly accessible and freely available for non-commercial research purposes. If you are publishing scientific work based on the IAM-Graph DB, we request you to include the following reference to our database:

@unpublished{riesen08iam,
	Author = {Riesen, K. and Bunke, H.},
	Note = {accepted for publication in SSPR 2008},
	Title = {{IAM} Graph Database Repository for Graph Based Pattern Recognition and Machine Learning}}


=======  DATA SET ======= 

Each of the data sets available on this repository is divided into three disjoint subsets, which can be used for training, validation, and testing novel learning algorithms (train.cxl, valid.cxl, test.cxl). 
This graph data set involves graphs that represent distorted letter drawings. We consider the 15 capital letters of the Roman alphabet that consist of straight lines only (A, E, F, H, I, K, L, M, N, T, V, W, X, Y, Z). For each class, a prototype line drawing is manually constructed. These prototype drawings are then converted into prototype graphs by representing lines by undirected edges and ending points of lines by nodes. Each node is labeled with a two-dimensional attribute giving its position relative to a reference coordinate system. Edges are unlabeled. The graph database consists of a training set, a validation set, and a test set of size 750 each. The graphs are uniformly distributed over the 15 classes. In order to test classifiers under different conditions, distortions are applied on the prototype graphs with three different levels of strength, viz. low, medium and high. Hence, our experimental data set comprises 6,750 graphs altogether. The classification rates achieved on this data set are 99.6% (low), 94.0% (medium), and 90.0% (high). In sample.eps the prototype graph and a graph instance for each distortion level representing the letter A are illustrated.


=======  REFERENCES ======= 

This data set is employed in the following publications (this list does not claim to be exhaustive, of course):

@book{neuhaus07bridging,
	Author = {Neuhaus, M. and Bunke, H.},
	Publisher = {World Scientific},
	Title = {Bridging the Gap Between Graph Edit Distance and Kernel Machines},
	Year = {2007}}

@inproceedings{bunke07family,
	Author = {Bunke, H. and Riesen, K.},
	Booktitle = {Proc. 12th Iberoamerican Congress on Pattern Recognition},
	Editor = {Rueda, L. and Mery, D. and Kittler, J.},
	Pages = {20--31},
	Series = {LNCS 4756},
	Title = {A Family of Novel Graph Kernels for Structural Pattern Recognition},
	Year = {2007}}

@inproceedings{riesen06embedding,
	Author = {Riesen, K. and Neuhaus, M. and Bunke, H.},
	Booktitle = {Proc.\ 6th Int.\ Workshop on Graph Based Representations in Pattern Recognition},
	Editor = {Escolano, F. and Vento, M.},
	Pages = {383--393},
	Series = {LNCS 4538},
	Title = {Graph Embedding in Vector Spaces by Means of Prototype Selection},
	Year = {2007}}


=======  CONTACT INFORMATION ======= 

If you have any question concerning this data set, do not hesitate to contact me: riesen@iam.unibe.ch

