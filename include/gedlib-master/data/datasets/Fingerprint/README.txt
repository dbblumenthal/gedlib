======= TERMS OF USAGE ======= 

The IAM-Graph DB is publicly accessible and freely available for non-commercial research purposes. If you are publishing scientific work based on the IAM-Graph DB, we request you to include the following reference to our database:

@unpublished{riesen08iam,
	Author = {Riesen, K. and Bunke, H.,
	Note = {accepted for publication in SSPR 2008,
	Title = {{IAM Graph Database Repository for Graph Based Pattern Recognition and Machine Learning


=======  DATA SET ======= 

Each of the data sets available on this repository is divided into three disjoint subsets, which can be used for training, validation, and testing novel learning algorithms (train.cxl, valid.cxl, test.cxl).
Fingerprints are converted into graphs by filtering the images and extracting regions that are relevant [neuhaus05graph]. In order to obtain graphs from fingerprint images, the relevant regions are binarized and a noise removal and thinning procedure is applied. This results in a skeletonized representation of the extracted regions. Ending points and bifurcation points of the skeletonized regions are represented by nodes. Additional nodes are inserted in regular intervals between ending points and bifurcation points. Finally, undirected edges are inserted to link nodes that are directly connected through a ridge in the skeleton. Each node is labeled with a two-dimensional attribute giving its position. The edges are attributed with an angle denoting the orientation of the edge with respect to the horizontal direction. 

The fingerprint database used in our experiments is based on the NIST-4 reference database of fingerprints [nist4]. It consists of a training set of size 500, a validation set of size 300, and a test set of size 2,000. Thus, there are 2,800 fingerprint images totally out of the four classes arch, left, right, and whorl from the Galton-Henry classification system. Note that in our benchmark test only the four-class problem of fingerprint classification is considered, i.e. the fifth class tented arch is merged with the class arch. Therefore, the first class (arch) consists of about twice as many graphs as the other three classes (left, right, whorl). For examples of these fingerprint classes, see Fig. sample.eps. The classification rate achieved on this data set is 76.6%.


=======  REFERENCES ======= 

This data set is employed in the following publications (this list does not claim to be exhaustive, of course):

@inproceedings{neuhaus05graph,
	Author = {Neuhaus, M. and Bunke, H.},
	Booktitle = {Proc.\ 5th Int.\ Conf.\ on Audio- and Video-Based Biometric Person Authentication},
	Editor = {Kanade, T. and Jain, A. and Ratha, N.K.},
	Pages = {191--200},
	Publisher = {Springer},
	Series = {LNCS 3546},
	Title = {A Graph Matching Based Approach to Fingerprint Classification Using Directional Variance},
	Year = {2005}}

@manual{nist4,
	Author = {Watson, C.I. and Wilson, C.L.},
	Organization = {National Institute of Standards and Technology},
	Title = {NIST {S}pecial {D}atabase 4, {F}ingerprint {D}atabase},
	Year = {1992}}


@book{neuhaus07bridging,
	Author = {Neuhaus, M. and Bunke, H.,
	Publisher = {World Scientific,
	Title = {Bridging the Gap Between Graph Edit Distance and Kernel Machines,
	Year = {2007}}

@inproceedings{bunke07family,
	Author = {Bunke, H. and Riesen, K.,
	Booktitle = {Proc. 12th Iberoamerican Congress on Pattern Recognition,
	Editor = {Rueda, L. and Mery, D. and Kittler, J.,
	Pages = {20--31,
	Series = {LNCS 4756,
	Title = {A Family of Novel Graph Kernels for Structural Pattern Recognition,
	Year = {2007}}




=======  CONTACT INFORMATION ======= 

If you have any question concerning this data set, do not hesitate to contact me: riesen@iam.unibe.ch

