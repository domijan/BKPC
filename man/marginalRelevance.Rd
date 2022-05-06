\name{marginalRelevance}
\alias{marginalRelevance}

\title{Marginal Relevance}


\description{
Calculates Marginal Relevance of each feature (variable) useful for class (group) separation. The marginal relevance score is a ratio of the between-group to within-group sum of squares. 

}
\usage{
marginalRelevance(x, y)
}

\arguments{
  \item{x}{a data matrix.}
  \item{y}{a response vector. Should be a factor.}
}

\value{
An object of class \code{"marginalRelevance"} including:
  
  \item{score}{Marginal relevance score of each feature.}
  \item{rank}{The ranking in order of highest marginal relevance for each feature.}
  \item{orderedData}{Data matrix with columns reordered by the marginal relevance of the features.}
  \item{bestVars}{Features ordered by the marginal relevance.}

}
\references{
Dudoit S., J. Fridlyand, T. P. Speed: 
Comparison of discrimination methods for the classification of tumors using gene expression data.
Journal of the American Statistical Association, 2002, Volume 97
No 457, pp 77-87.
}
\author{K. Domijan}



\seealso{
  \code{\link{plot.marginalRelevance}} \code{\link{microarray}}
}
\examples{


data(microarray)

profiles <- as.matrix(microarray[, -2309])
tumourType <-  microarray[, 2309]
 
 margRelv <- marginalRelevance(profiles, tumourType)

# plot 30 gene profiles with highest marginal relevance score
 plot(margRelv, type = "parallelcoord", n.feat = 50, col = tumourType )

 
\dontrun{
# another example: wine data from gclus
library(gclus)
data(wine)
dt <- as.matrix(wine[, -1])
colnames(dt) <- names(wine[, -1])

label <- as.factor(wine[, 1])


margRelv <- marginalRelevance(dt, label)

#variables in order of their MR score
colnames(dt)[ margRelv$bestVars]


cparcoord(dt, order = margRelv$bestVars, col = label)
cpairs(dt, order = margRelv$bestVars, col = label)

}
}

