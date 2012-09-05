\name{fSAE.Area}
\alias{fSAE.Area}
\title{Compute small area estimates based on the basic area-level model.}
\usage{
  fSAE.Area(est.init, var.init, X, x = X, ...)
}
\arguments{
  \item{est.init}{m-vector of initial estimates, where m is
  the number of sampled areas.}

  \item{var.init}{m-vector of corresponding variance
  estimates.}

  \item{X}{M x p matrix of area population means, where M
  is the number of areas for which estimates are required.}

  \item{x}{defaults to X, but in case auxiliary information
  is available at the unit level, it may be set to the
  corresponding matrix of covariate sample means.}

  \item{...}{additional arguments passed to
  \code{fSAE.Unit}.}
}
\value{
  An object of class \code{sae} containing the small area
  estimates and MSEs, the model fit, and model selection
  measures.
}
\description{
  This function returns small area estimates based on the
  basic area-level model, also known as the Fay-Herriot
  model. It calls \code{\link{fSAE.Unit}} to carry out the
  computations.
}
\examples{
d <- generateFakeData()

# first compute input estimates without "borrowing strength" over areas
sae0 <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, type="direct", keep.data=TRUE)

# compute small area estimates based on the basic area-level model
#   using the above survey regression estimates as input
sae <- fSAE.Area(EST(sae0), MSE(sae0), X=sae0$Xp)
EST(sae)  # estimates
SE(sae)  # standard errors
}
\references{
  R.E. Fay and R.A. Herriot (1979). Estimates of Income for
  Small Places: An Application of James-Stein Procedures to
  Census Data.  Journal of the American Statistical
  Association 74(366), 269-277.

  J.N.K. Rao (2003). Small Area Estimation. Wiley.
}
\seealso{
  \code{\link{sae-class}}
}

