======= TERMS OF USAGE ======= 

The IAM-Graph DB is publicly accessible and freely available for non-commercial research purposes. If you are publishing scientific work based on the IAM-Graph DB, we request you to include the following reference to our database:

@unpublished{riesen08iam,
	Author = {Riesen, K. and Bunke, H.,
	Note = {accepted for publication in SSPR 2008,
	Title = {{IAM Graph Database Repository for Graph Based Pattern Recognition and Machine Learning


=======  DATA SET ======= 

Each of the data sets available on this repository is divided into three disjoint subsets, which can be used for training, validation, and testing novel learning algorithms (train.cxl, valid.cxl, test.cxl). 
The GREC data set consists of graphs representing symbols from architectural and electronic drawings. The images occur at five different distortion levels. In Fig. sample.eps for each distortion level one example of a drawing is given. Depending on the distortion level, either erosion, dilation, or other morphological operations are applied. The result is thinned to obtain lines of one pixel width. Finally, graphs are extracted from the resulting denoised images by tracing the lines from end to end and detecting intersections as well as corners. Ending points, corners, intersections and circles are represented by nodes and labeled with a two-dimensional attribute giving their position. The nodes are connected by undirected edges which are labeled as line or arc. An additional attribute specifies the angle with respect to the horizontal direction or the diameter in case of arcs. From the original GREC database [dosch06report], 22 classes are considered. For an adequately sized set, all graphs are distorted nine times to obtain a data set containing 1,100 graphs uniformely distributed over the 22 classes. The resulting set is split into a traininig and a validation set of size 286 each, and a test set of size 528. The classification rate achieved on this data set is 95.5%.

=======  REFERENCES ======= 

This data set is employed in the following publications (this list does not claim to be exhaustive, of course):

@inproceedings{dosch06report,
	Author = {Dosch, Ph. and Valveny, E.},
	Booktitle = {Graphics Recognition. Ten years review and future perspectives. Proc.\ 6th Int.\ Workshop on Graphics Recognition (GREC'05)},
	Editor = {Liu Wenyin and Llad{\'o}s, J.},
	Pages = {381--397},
	Publisher = {Springer},
	Series = {LNCS 3926},
	Title = {Report on the Second Symbol Recognition Contest},
	Year = {2005}}


=======  CONTACT INFORMATION ======= 

If you have any question concerning this data set, do not hesitate to contact me: riesen@iam.unibe.ch

