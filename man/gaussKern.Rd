\name{gaussKern}
\alias{gaussKern}

\title{Gaussian kernel}

\description{
Calculates Gaussian kernel: \eqn{k(x,x') = \exp(-\theta \|x - x'\|^2)} \cr
}
\usage{
gaussKern(x, newdata = x, theta = NULL)
}

\arguments{
  \item{x}{a data matrix.}
  \item{newdata}{optional second data matrix.}
  \item{theta}{the inverse kernel bandwidth parameter. If \code{NULL} a default value is used \eqn{\theta = 1/max(\|x - x'\|^2)}.}
}
\details{Also known as the radial basis kernel function, see \code{\link[kernlab]{rbfdot}} (in package \pkg{kernlab})}

\value{
 Returns a list containing the following components:
   \item{K}{a Gaussian kernel matrix of class \code{"kern"}.}
      \item{theta}{the inverse kernel bandwidth parameter.}
}
\author{K. Domijan}


\seealso{
  \code{\link{kPCA}}
  \code{\link{bkpc}}
   \code{\link[kernlab]{kernelMatrix}} (in package \pkg{kernlab})
   }
  
  
\examples{

data(iris)

testset <- sample(1:150,20)
train <- as.matrix(iris[-testset ,-5])
test <- as.matrix(iris[testset ,-5])


# make training set kernel
gk <- gaussKern(train)
Ktrain <- gk$K

image(Ktrain)

# make testing set kernel
gk2 <- gaussKern(train, test, gk$theta) 
Kest <- gk2$K

}

