\name{plot.marginalRelevance}
\alias{plot.marginalRelevance}

\title{
Plot marginalRelevance Objects
}
\description{
Plots marginal relevance scores for features of a given data matrix. The default plot shows: the marginal relevance score (MR score) of each feature. The "pairs" and "parallelcoord" show scatterplot matrix and the parallel coordinates plot of features ordered by their MR score. 
}
\usage{
\method{plot}{marginalRelevance}(x, newdata = NULL, n.feat = NULL, type = "default", ...) 
}

\arguments{
  \item{x}{a \code{marginalRelevance} object.}
  \item{newdata}{a matrix containing the new input data.}
  \item{n.feat}{the number of features with highest MR score to plot. Default is all features. }
  \item{type}{"parallelcoord", "pairs" or default.}
  \item{\dots}{options directly passed to the plot function.}
}

\details{
If \code{newdata} is omitted the predictions are based on the data used for deriving the MR score. 
}

\author{K. Domijan}

\seealso{
  \code{\link{marginalRelevance}}
}

\examples{

data(microarray)

profiles <- as.matrix(microarray[, -2309])
tumourType <-  microarray[, 2309]
 
 margRelv <- marginalRelevance(profiles, tumourType)


# plot 30 gene profiles with highest marginal relevance score
 plot(margRelv, type = "parallelcoord", n.feat = 50, col = tumourType )


\dontrun{

library(kernlab)
data(spam)

test <- sample(1:4601,2000)


dt <- as.matrix(spam[-test ,-58])
labels <- spam[-test , 58]

 
margRelv <- marginalRelevance(dt, labels)

#plot MR scores
plot(margRelv)

plot(margRelv , col = labels, type = "pairs", n.feat = 5)
plot(margRelv , col = labels, type = "parallelcoord", n.feat = 30)


# test set
plot(margRelv , as.matrix(spam[test ,-58]), col = spam[test , 58], 
type = "pairs", n.feat = 5)

plot(margRelv , as.matrix(spam[test ,-58]), col = spam[test , 58], 
type = "parallelcoord", n.feat = 30)



 

 }
}
