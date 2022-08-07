\name{bkpc}
\alias{bkpc}
\alias{bkpc.default}
\alias{bkpc.kern}
\alias{bkpc.kernelMatrix}

\title{Bayesian Kernel Projection Classifier}
\description{Function \code{bkpc} is used to train a Bayesian kernel projection classifier. This is a nonlinear multicategory classifier  which performs the classification of the projections of the data to the principal axes of the feature space. The Gibbs sampler is implemented to find the posterior distributions of the parameters, so probability distributions of prediction can be obtained for for new observations. }
\usage{


\method{bkpc}{default}(x, y, theta = NULL, n.kpc = NULL, thin = 100, n.iter = 1e+05, std = 10, 
g1 = 0.001, g2 = 0.001, g3 = 1, g4 = 1, initSigmasq = NULL, initBeta = NULL,
initTau = NULL, intercept = TRUE, rotate = TRUE, ...)


\method{bkpc}{kern}(x, y, n.kpc = NULL, thin = 100, n.iter = 1e+05, std = 10, 
g1 = 0.001, g2 = 0.001, g3 = 1, g4 = 1, initSigmasq = NULL, initBeta = NULL, 
initTau = NULL, intercept = TRUE, rotate = TRUE, ...)


\method{bkpc}{kernelMatrix}(x, y, n.kpc = NULL, thin = 100, n.iter = 1e+05, std = 10, 
g1 = 0.001, g2 = 0.001, g3 = 1, g4 = 1, initSigmasq = NULL, initBeta = NULL, 
initTau = NULL, intercept = TRUE, rotate = TRUE, ...)
}

\arguments{
  \item{x}{either: a data matrix, a kernel matrix of class \code{"kernelMatrix"} or a kernel matrix of class \code{"kern"}.}
  \item{y}{a response vector with one label for each row of \code{x}. Should be a factor.}
  \item{theta}{the inverse kernel bandwidth parameter.}
  \item{n.kpc}{number of kernel principal components to use.}
  \item{n.iter}{number of iterations for the MCMC algorithm.}
  \item{thin}{thinning interval.}
  \item{std}{standard deviation parameter for the random walk proposal.}
  \item{g1}{\eqn{\gamma_1} hyper-parameter of the prior inverse gamma distribution for the \eqn{\sigma^2}  parameter in the BKPC model.}
  \item{g2}{\eqn{\gamma_2} hyper-parameter of the prior inverse gamma distribution for the \eqn{\sigma^2} parameter of the BKPC model.}
  \item{g3}{\eqn{\gamma_3} hyper-parameter of the prior gamma distribution for the \eqn{\tau} parameter in the BKPC model.}
  \item{g4}{\eqn{\gamma_4} hyper-parameter of the prior gamma distribution for the \eqn{\tau} parameter in the BKPC model.}
    \item{initSigmasq}{optional specification of initial value for the \eqn{\sigma^2} parameter
in the BKPC model.}
  \item{initBeta}{optional specification of initial values for the \eqn{\beta} parameters
in the BKPC model.}
  \item{initTau}{optional specification of initial values for the \eqn{\tau} parameters
in the BKPC model.}
  \item{intercept}{if \code{intercept=TRUE}  (the default) then include the intercept in the model.}
  \item{rotate}{if \code{rotate=TRUE}  (the default) then run the BKPC model. Else run the BKMC model. }
  \item{\dots}{Currently not used.}
}
\details{Initial values for a BKPC model can be  supplied, otherwise they are generated using \code{runif} function. 

The data can be passed to the \code{bkpc} function in a matrix and the Gaussian kernel computed using the \code{gaussKern} function is then used in training the algorithm and predicting. The bandwidth parameter \code{theta} can be supplied to the \code{gaussKern} function, else a default value is used. 

In addition, \code{bkpc} also supports input in the form of a kernel matrix of class \code{"kern"} or \code{"kernelMatrix"}.The latter allows for a range of kernel functions as well as user specified ones.  

If \code{rotate=TRUE}  (the default) then the BKPC is trained. This algorithm performs the classification of the projections of the data to the principal axes of the feature space. Else the Bayesian kernel multicategory classifier (BKMC) is trained, where the data is mapped to the feature space via the kernel matrix, but not projected (rotated) to the principal axes. The hierarchichal prior structure for the two models is the same, but BKMC model is not sparse. 
}


\value{
  An object of class \code{"bkpc"} including:
  
  \item{beta}{realizations of the \eqn{\beta} parameters from the joint posterior distribution in the BKPC model.}
  \item{tau}{realizations of the \eqn{\tau} parameters from the joint posterior distribution in the BKPC model.}
  \item{z}{realizations of the latent variables \eqn{z} from the joint posterior distribution in the BKPC model.}
 \item{sigmasq}{realizations of the \eqn{\sigma^2} parameter from the joint posterior distribution in the BKPC model.}

  \item{n.class}{number of independent classes of the response variable i.e. number of classes - 1.}
  \item{n.kpc}{number of kernel principal components used.}
 \item{n.iter}{number of iterations of the MCMC algorithm.}
  \item{thin}{thinning interval.}
  \item{intercept}{if true, intercept was included in the model.}
  \item{rotate}{if true, the sparse BKPC model was fitted, else BKMC model.}
  \item{kPCA}{if \code{rotate=TRUE} an object of class \code{"kPCA"}, else \code{NULL}. }
  \item{x}{the supplied data matrix or kernel matrix.}
  \item{theta}{if data was supplied, as opposed to the kernel, this is the inverse kernel bandwidth parameter used in obtaining the Gaussian kernel, else \code{NULL}.}
  

}

\references{
Domijan K. and Wilson S. P.:
Bayesian kernel projections for classification of high dimensional data.
Statistics and Computing, 2011, Volume 21, Issue 2, pp 203-216
}

\author{K. Domijan}

\note{
If supplied, data are not scaled internally. If \code{rotate=TRUE}  the mapping is centered internally by the \code{kPCA} function. 
}



\seealso{
  \code{\link{kPCA}}
  \code{\link{gaussKern}}
  \code{\link{predict.bkpc}}
  \code{\link{plot.bkpc}}
  \code{\link{summary.bkpc}}
  \code{\link[kernlab]{kernelMatrix}} (in package \pkg{kernlab})
}
\examples{

set.seed(-88106935)

data(microarray)

# consider only four tumour classes (NOTE: "NORM" is not a class of tumour)
y <- microarray[, 2309]
train <- as.matrix(microarray[y != "NORM", -2309])
wtr <- factor(microarray[y != "NORM", 2309], levels = c("BL" ,  "EWS" , "NB" ,"RMS" ))

n.kpc <- 6
n.class <- length(levels(wtr)) - 1

K <- gaussKern(train)$K

# supply starting values for the parameters
# use Gaussian kernel as input

result <- bkpc(K, y = wtr, n.iter = 1000,  thin = 10, n.kpc = n.kpc,  
initSigmasq = 0.001, initBeta = matrix(10, n.kpc *n.class, 1), 
initTau =matrix(10, n.kpc * n.class, 1), intercept = FALSE, rotate = TRUE)

# predict

out <- predict(result, n.burnin = 10) 

table(out$class, as.numeric(wtr))

# plot the data projection on the kernel principal components

pairs(result$kPCA$KPCs[, 1 : n.kpc], col = as.numeric(wtr), 
main =  paste("symbol = predicted class", "\n", "color = true class" ), 
pch = out$class, upper.panel = NULL)
par(xpd=TRUE)
legend('topright', levels(wtr), pch = unique(out$class), 
text.col = as.numeric(unique(wtr)), bty = "n")




# Another example: Iris data

data(iris)
testset <- sample(1:150,50)

train <- as.matrix(iris[-testset,-5])
test <- as.matrix(iris[testset,-5])

wtr <- iris[-testset, 5]
wte <- iris[testset, 5]

# use default starting values for paramteres in the model.

result <- bkpc(train, y = wtr,  n.iter = 1000,  thin = 10, n.kpc = 2, 
intercept = FALSE, rotate = TRUE)

# predict
out <- predict(result, test, n.burnin = 10) 

# classification rate
sum(out$class == wte)/dim(test)[1]

table(out$class, wte)

\dontrun{
# Another example: synthetic data from MASS library

library(MASS)

train<- as.matrix(synth.tr[, -3])
test<- as.matrix(synth.te[, -3])

wtr <- as.factor(synth.tr[, 3])
wte <- as.factor(synth.te[, 3])


#  make training set kernel using kernelMatrix from kernlab library

library(kernlab)

kfunc <- laplacedot(sigma = 1)
Ktrain <- kernelMatrix(kfunc, train)

#  make testing set kernel using kernelMatrix {kernlab}

Ktest <- kernelMatrix(kfunc, test, train)

result <- bkpc(Ktrain, y = wtr, n.iter = 1000,  thin = 10,  n.kpc = 3, 
intercept = FALSE, rotate = TRUE)

# predict

out <- predict(result, Ktest, n.burnin = 10) 

# classification rate

sum(out$class == wte)/dim(test)[1]
table(out$class, wte)


# embed data from the testing set on the new space:

KPCtest <- predict(result$kPCA, Ktest)

# new data is linearly separable in the new feature space where classification takes place
library(rgl)
plot3d(KPCtest[ , 1 :  3], col = as.numeric(wte))


# another model:  do not project the data to the principal axes of the feature space. 
# NOTE: Slow
# use Gaussian kernel with the default bandwidth parameter

Ktrain <- gaussKern(train)$K

Ktest <- gaussKern(train, test, theta = gaussKern(train)$theta)$K

resultBKMC <- bkpc(Ktrain, y = wtr, n.iter = 1000,  thin = 10,  
intercept = FALSE, rotate = FALSE)

# predict
outBKMC <- predict(resultBKMC, Ktest, n.burnin = 10)

# to compare with previous model
table(outBKMC$class, wte)


# another example: wine data from gclus library

library(gclus)
data(wine)

testset <- sample(1 : 178, 90)
train <- as.matrix(wine[-testset, -1])
test <- as.matrix(wine[testset, -1])

wtr <- as.factor(wine[-testset, 1])
wte <- as.factor(wine[testset, 1])

#  make training set kernel using kernelMatrix from kernlab library

kfunc <- anovadot(sigma = 1, degree = 1)
Ktrain <- kernelMatrix(kfunc, train)

#  make testing set kernel using kernelMatrix {kernlab}
Ktest <- kernelMatrix(kfunc, test, train)

result <- bkpc(Ktrain, y = wtr, n.iter = 1000,  thin = 10,  n.kpc = 3, 
intercept = FALSE, rotate = TRUE)

out <- predict(result, Ktest, n.burnin = 10) 

# classification rate in the test set
sum(out$class == wte)/dim(test)[1]


# embed data from the testing set on the new space:
KPCtest <- predict(result$kPCA, Ktest)

# new data is linearly separable in the new feature space where classification takes place


pairs(KPCtest[ , 1 :  3], col = as.numeric(wte), 
main =  paste("symbol = predicted class", "\n", "color = true class" ), 
pch = out$class, upper.panel = NULL)

par(xpd=TRUE)

legend('topright', levels(wte), pch = unique(out$class), 
text.col = as.numeric(unique(wte)), bty = "n")


}
}

\keyword{neural}
\keyword{nonlinear}
\keyword{classif}
