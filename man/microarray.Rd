\name{microarray}
\alias{microarray}


\docType{data}

\title{Microarray gene expression data}

\description{
Microarray gene expression data published by Khan et al. (2001). There are 2308 gene expression profiles recorded over 88 arrays. 
}

\usage{data("microarray")}
\format{
  A data frame with 88 observations on the following 2309 variables.
  
The first 2308 variables are the gene expression values for 88 arrays. The first 63 arrays correspond to the training set and the remaining 25 are from the testing set of the filtered data made available in the supplementary files. The last variable is the tumour class, a factor with levels \code{BL}, \code{EWS}, \code{NB}, \code{NORM}, \code{RMS}. 
  }


\source{
These data have been taken from the supplementary web site to the publication \url{http://research.nhgri.nih.gov/microarray/Supplement/}.}


\references{
Khan, J., Wei, J.S., Ringner, M., Saal, L.H., Ladanyi, M., Westermann, F., Berthold, F., Schwab, M., Antonescu, C.R., Peterson, C. et al.: Classification and diagnostic prediction of cancers using gene expression profiling and artificial neural networks. 2001. Nat. Med., 7, 673-679. 
}

\keyword{datasets}
