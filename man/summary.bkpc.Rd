\name{summary.bkpc}
\alias{summary.bkpc}

\title{

Summary statistics for Markov Chain Monte Carlo chain from Bayesian Kernel Projection Classifier
}
\description{
\code{summary.bkpc} produces two sets of summary statistics for each variable: mean and standard deviation (ignoring autocorrelation of the chain) of the sample distribution and quantiles of the sample distribution using the quantiles argument.
}
\usage{
\method{summary}{bkpc}(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), n.burnin = 0, ...)
}

\arguments{
  \item{object}{an object of class \code{"bkpc"}.}
  \item{quantiles}{a vector of quantiles to evaluate for each variable.}
  \item{n.burnin}{number of burn-in iterations to discard from the thinned sample. }
  \item{\dots}{Currently not used.}
}


\author{
K. Domijan
}

\seealso{
\code{\link{bkpc}}
\code{\link{plot.bkpc}}
}
\examples{

set.seed(-88106935)

data(iris)
testset <- sample(1:150,50)

train <- as.matrix(iris[-testset,-5])
test <- as.matrix(iris[testset,-5])

wtr <- iris[-testset, 5]
wte <- iris[testset, 5]

result <- bkpc(train, y = wtr,  n.iter = 1000,  thin = 10, n.kpc = 2, 
intercept = FALSE, rotate = TRUE)


summary(result, n.burnin = 0)
}
