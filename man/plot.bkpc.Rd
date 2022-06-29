\name{plot.bkpc}
\alias{plot.bkpc}

\title{
Plot bkpc Objects
}
\description{Plots realizations of the parameters from the joint posterior distribution in the BKPC model. The default plots show: medians, 10th and 90th percentiles. The "tracePlot" and "boxPlot" show the traceplots and boxplots of the samples. 
}
\usage{
\method{plot}{bkpc}(x, type = "default", n.burnin = 0, ...)
}

\arguments{
  \item{x}{a \code{bkpc} object.}
  \item{type}{"tracePlot", "boxPlot" or default.}
   \item{n.burnin}{number of burn-in iterations from the thinned sample to discard. }
  \item{\dots}{options directly passed to the plot function.}
}



\author{K. Domijan}

\seealso{
\code{\link{bkpc}}
\code{\link{summary.bkpc}}
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

result <- bkpc(K, y = wtr, n.iter = 10000,  thin = 100, n.kpc = n.kpc,  
initSigmasq = 0.001, initBeta = matrix(10, n.kpc *n.class, 1), 
initTau =matrix(10, n.kpc * n.class, 1), intercept = FALSE, rotate = TRUE)



plot(result, type = "tracePlot")
plot(result, type = "boxPlot", n.burnin = 20)
plot(result, n.burnin = 20)
}


