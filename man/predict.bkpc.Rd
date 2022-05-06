\name{predict.bkpc}
\alias{predict.bkpc}
\alias{predictMultinomSamples}

\title{
Predict Method for Bayesian Kernel Projection Classifier
}
\description{
This function predicts values based upon a model trained by \code{bkpc} for new input data. BKPC employs a Gibbs sampler to sample from the posterior distributions of the parameters, so sample probability distributions of prediction can be obtained for for new data points.
}
\usage{
\method{predict}{bkpc}(object, newdata  = NULL, n.burnin = 0, ...)
}

\arguments{
  \item{object}{a \code{bkpc} object.}
  \item{newdata}{a matrix containing the new input data}
  \item{n.burnin}{number of burn-in iterations to discard.}
  \item{\dots}{Currently not used.}
}
\details{
If \code{newdata} is omitted the predictions are based on the data used for the fit. 
}
\value{  
A list with the following components:
    \item{class}{estimated class for each observation in the input data.}
  \item{map}{maximum a posteriori probability estimate for belonging to each class for all observations in the input data.}
    \item{p}{a matrix of samples of estimated probabilities for belonging to each class for each observation in the input data.}

}
\references{
Domijan K. and Wilson S. P.:
Bayesian kernel projections for classification of high dimensional data.
Statistics and Computing, 2011, Volume 21, Issue 2, pp 203-216
}
\author{K. Domijan}

\seealso{
  \code{\link{bkpc}}
}
\examples{
set.seed(-88106935)

data(iris)
testset <- sample(1:150,30)

train <- as.matrix(iris[-testset,-5])
test <- as.matrix(iris[testset,-5])

wtr <- iris[-testset, 5]
wte <- iris[testset, 5]

result <- bkpc(train, y = wtr,  n.iter = 1000,  thin = 10, n.kpc = 2, 
intercept = FALSE, rotate = TRUE)

# predict
out <- predict(result, test, n.burnin = 20) 

# classification rate for the test set

sum(out$class == as.numeric(wte))/dim(test)[1]

table(out$class, as.numeric(wte))

# consider just misclassified observations:

missclassified <- out$class != as.numeric(wte)


tab <- cbind(out$map[missclassified, ], out$class[missclassified],  as.numeric(wte)[missclassified])
colnames(tab) = c("P(k = 1)", "P(k = 2)", "P(k = 3)", "predicted class", "true class")
tab


# consider, say, 28th observation in the test set:
# sample probability distributions of belonging to each of the three classes: 


ProbClass2samples <- out$p[28, ]
ProbClass3samples <- out$p[28 + dim(test)[1], ]
ProbClass1samples <- 1 - (ProbClass2samples + ProbClass3samples)
hist(ProbClass1samples)
hist(ProbClass2samples)
hist(ProbClass3samples)


}

