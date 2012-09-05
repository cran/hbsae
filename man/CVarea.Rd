\name{CVarea}
\alias{CVarea}
\title{Compute area-level cross-validation measure for sae objects.}
\usage{
  CVarea(sae, weight = TRUE,
    cost = function(y, yhat, w) sum(w * (y - yhat)^2)/sum(w),
    K = 10, method = "hybrid", seed)
}
\arguments{
  \item{sae}{object of class sae, resulting from a call to
  \code{fSAE}, \code{fSAE.Area}, or \code{fSAE.Unit}.}

  \item{weight}{if \code{TRUE}, use weights inversely
  proportional to the MSEs of y - yhat in the cost
  function.}

  \item{cost}{cost function to be used. Defaults to a
  quadratic cost function.}

  \item{K}{K in K-fold cross-validation. Specifies in how
  many parts the dataset should be divided.}

  \item{method}{method used to refit the model. One of
  "HB", "hybrid" (default) or "REML", in the order of slow
  to fast.}

  \item{seed}{random seed used in selecting groups of areas
  to leave out in K-fold cross-validation.}
}
\value{
  The computed area-level cross-validation measure.
}
\description{
  This function computes a cross-validation measure defined
  at the area level. It can be used, for example, to
  compare the predictive ability of area and unit-level
  models. The code is based in part on that of
  \code{cv.glm} from package \pkg{boot}.
}
\examples{
d <- generateFakeData()

# compute small area estimates based on area-level and unit-level models
saeArea <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, type="area",
  silent=TRUE, keep.data=TRUE)
saeUnit <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, type="unit",
  silent=TRUE, keep.data=TRUE)

# compare area and unit-level models based on area-level cross-validation criterion
CVarea(saeArea, K=10, seed=1)  # 10-fold CV for area-level model
CVarea(saeUnit, K=10, seed=1)  # 10-fold CV for unit-level model
}

