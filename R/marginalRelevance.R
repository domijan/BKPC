marginalRelevance <- function (x, y)
{
  e1 <- new.env()
  if (!is.factor(y))
    stop("error: y is not a factor")
  if (anyNA(y)) stop("error: NA in response vector")
  NF <- dim(x)[2]
  NK <- length(unique(y))
  u <- matrix(0, 1, NF)
  l <- matrix(0, 1, NF)
  xOrdered <- matrix(0, dim(x)[1], NF)
  bestVars <- matrix(0, 1, NF)
  x_dot_f <- apply(x, 2, mean, na.rm = TRUE)
  for (k in 1:NK) {
    u <- u + sum(y == levels(y)[k]) * ((apply(x[y == levels(y)[k],
                                                ], 2, mean, na.rm = TRUE) - x_dot_f)^2)
    l <- l + apply((x[y == levels(y)[k], ] - kronecker(matrix(1,
                                                              sum(y == levels(y)[k]), 1), t(as.matrix(apply(x[y ==
                                                                                                                levels(y)[k], ], 2, mean, na.rm = TRUE)))))^2, 2, sum, na.rm = TRUE)
  }
  e1$score <- u/l

  e1 <- as.list(e1)
  return(e1)
}
