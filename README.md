# BKPC R package. 


Bayesian kernel projection classifier (BKPC) is a nonlinear multicategory classifier which performs the classification of the projections of the data to the principal axes of the feature space. A Gibbs sampler is implemented to find the posterior distributions of the parameters. 

The main function is bkpc. 

The data can be passed to the bkpc function as:

- a matrix of features,

- a kernel matrix of either:

* class 'kern' (a Gaussian kernel computed using the gaussKern{BKPC} function) or 
* class 'kernelMatrix' from library kernlab. This allows for a wider selection of inbuilt kernel generating functions as well as user defined functions.

The package contains a microarray dataset and a function to extract the marginal relevance of each feature for classification.

See ?bkpc and ?marginalRelevance





