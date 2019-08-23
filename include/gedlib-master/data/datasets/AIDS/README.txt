======= TERMS OF USAGE ======= 

The IAM-Graph DB is publicly accessible and freely available for non-commercial research purposes. If you are publishing scientific work based on the IAM-Graph DB, we request you to include the following reference to our database:

@unpublished{riesen08iam,
	Author = {Riesen, K. and Bunke, H.,
	Note = {accepted for publication in SSPR 2008,
	Title = {{IAM Graph Database Repository for Graph Based Pattern Recognition and Machine Learning


=======  DATA SET ======= 

Each of the data sets available on this repository is divided into three disjoint subsets, which can be used for training, validation, and testing novel learning algorithms (train.cxl, valid.cxl, test.cxl). 
The AIDS data set consists of graphs representing molecular compounds. We construct graphs from the AIDS Antiviral Screen Database of Active Compounds [molecules]. This data set consists of two classes (active, inactive), which represent molecules with activity against HIV or not. The molecules are converted into graphs in a straightforward manner by representing atoms as nodes and the covalent bonds as edges. Nodes are labeled with the number of the corresponding chemical symbol and edges by the valence of the linkage. In sample.eps one molecular compound of both classes is illustrated. Note that different shades of grey represent different chemical symbols, i.e. node labels. We use a training set and a validation set of size 250 each, and a test set of size 1,500. Thus, there are 2,000 elements totally (1,600 inactive elements and 400 active elements). The classification result achieved on this data set is 97.3%.

=======  REFERENCES ======= 

This data set is employed in the following publications (this list does not claim to be exhaustive, of course):


@url{molecules,
	Author = {Development Therapeutics Program {DTP}},
	Note = {{\tt http://dtp.nci.nih.gov/docs/aids/aids\_data.html}},
	Title = {{AIDS} Antiviral Screen},
	Year = {2004}}

@book{neuhaus07bridging,
	Author = {Neuhaus, M. and Bunke, H.,
	Publisher = {World Scientific,
	Title = {Bridging the Gap Between Graph Edit Distance and Kernel Machines,
	Year = {2007

@inproceedings{bunke07family,
	Author = {Bunke, H. and Riesen, K.,
	Booktitle = {Proc. 12th Iberoamerican Congress on Pattern Recognition,
	Editor = {Rueda, L. and Mery, D. and Kittler, J.,
	Pages = {20--31,
	Series = {LNCS 4756,
	Title = {A Family of Novel Graph Kernels for Structural Pattern Recognition,
	Year = {2007




=======  CONTACT INFORMATION ======= 

If you have any question concerning this data set, do not hesitate to contact me: riesen@iam.unibe.ch

