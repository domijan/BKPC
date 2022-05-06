\name{kPCA}
\alias{kPCA}
\alias{kPCA.default}
\alias{kPCA.kern}
\alias{kPCA.kernelMatrix}
\alias{predict.kPCA}
\alias{getPrincipalComponents}

\title{
Kernel Principal Components Analysis
}
\description{
Kernel PCA is a nonlinear generalization of principal component analysis.
}
\usage{

\method{kPCA}{default}(x, theta = NULL, ...)

\method{kPCA}{kern}(x, ...)

\method{kPCA}{kernelMatrix}(x, ...)
}

\arguments{
  \item{x}{either: a data matrix, a kernel matrix of class \code{"kernelMatrix"} or a kernel matrix of class \code{"kern"}.}
   \item{theta}{the inverse kernel bandwidth parameter for the Gaussian kernel.}
  \item{\dots}{Currently not used.}
}
\details{ 
The data can be passed to the \code{kPCA} function in a \code{matrix} and the Gaussian kernel (via the \code{gaussKern} function) is used to map the data to the high-dimensional feature space where the principal components are computed. The bandwidth parameter \code{theta} can be supplied to the \code{gaussKern} function, else a default value is used. Alternatively, the Gaussian kernel matrix of class \code{"kern"} can be supplied to the \code{kPCA} function directly. In addition, \code{kPCA} also supports input in the form of a kernel matrix of class \code{"kernelMatrix"} (in package \pkg{kernlab}) thus allowing for other kernel functions.
}
\value{
  An object of class \code{"kPCA"} including:
  
\item{KPCs}{The original data projected  on the principal components.}
\item{Es}{the corresponding eigenvalues.}
\item{Vecs}{a matrix containing principal component vectors (column wise).}
\item{K}{a kernel matrix of class \code{"kernelMatrix"} or of class \code{"kern"}. }
\item{theta}{if Gaussian kernel was calculated, this is the bandwidth parameter used in its calculation.}
\item{x}{if supplied, the original data matrix.}
}
\references{
Schoelkopf B., A. Smola, K.-R. Mueller :
Nonlinear component analysis as a kernel eigenvalue problem. 
Neural Computation 10, 1299-1319. 
}
\author{K. Domijan}

\note{The predict function can be used to embed new data on the new space}



\seealso{
\code{\link{gaussKern}}
 \code{\link[kernlab]{kpca}} (in package \pkg{kernlab})
 \code{\link[kernlab]{kernelMatrix}} (in package \pkg{kernlab})
}
\examples{

data(iris)
testset <- sample(1:150,20)

train <- as.matrix(iris[-testset,-5])
test <- as.matrix(iris[testset,-5])


# make training set kernel

gk <- gaussKern(train)
Ktrain <- gk$K
image(Ktrain)

# make testing set kernel

gk2 <- gaussKern(train, test, gk$theta) 
Ktest <- gk2$K


#  make training set kernel using kernelMatrix from library kernlab.

library(kernlab)

kfunc <- laplacedot(sigma = 1)
Ktrain2 <- kernelMatrix(kfunc, train)
image(Ktrain2)

# make testing set kernel using kernelMatrix {kernlab}

Ktest2 <- kernelMatrix(kfunc, test, train)



# Do KPCA:

kpcData <- kPCA(train)
kpcKern <- kPCA(Ktrain)
kpcKern2 <- kPCA(Ktrain2)


# plot the data projection on the principal components

pairs(kpcData$KPCs[ , 1 : 3], col = iris[-testset, 5])

# proportion of variance explained by each PC

plot(kpcData$Es/sum(kpcData$Es), xlab = "PC", ylab = "Proportion of variance")


# embed data from the testing set on the new space:


KPCpred1 <- predict(kpcData, test)

KPCpred2 <- predict(kpcKern, Ktest)

KPCpred3 <- predict(kpcKern2, Ktest2)

#plot the test data projection on the principal components

pairs(KPCpred3[ , 1 : 3], col = iris[testset, 5])




}

\keyword{nonlinear}
\keyword{cluster}
