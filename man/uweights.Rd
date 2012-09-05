\name{uweights}
\alias{uweights}
\title{Compute unit weights underlying the small area estimates or their aggregate.}
\usage{
  uweights(x, areaID = NULL, forTotal = FALSE)
}
\arguments{
  \item{x}{sae object.}

  \item{areaID}{if left unspecified (\code{NULL}), weights
  corresponding to the overall (aggregated) estimate are
  returned.  Otherwise weights that reproduce the estimate
  for a specific area are returned.}

  \item{forTotal}{if \code{FALSE} weights will be divided
  by the population size.}
}
\value{
  An object of class \code{weights}.
}
\description{
  The small area estimates can be interpreted as weighted
  sums of the response variable. This function computes the
  weights corresponding to the aggregated small area
  estimates or the weights corresponding to a specific
  small area estimate. The weights applied to the response
  variable need not exactly reproduce the Hierarchical
  Bayes estimate since the latter is averaged over the
  posterior distribution for the variance ratio whereas the
  weights are evaluated at the posterior mean. Under the
  default prior for the fixed effects, the weights applied
  to the design matrix reproduce the corresponding
  population numbers.
}
\examples{
d <- generateFakeData()

# compute small area estimates
sae <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, method="hybrid", keep.data=TRUE)

# compute unit weights
w <- uweights(sae, forTotal=TRUE)
summary(w)  # summary statistics
plot(w)  # histogram of weights
# checks
all.equal(sum(w * sae$y), sum(EST(sae) * sae$Narea))
all.equal(colSums(w * as.matrix(sae$X)), colSums(sae$Xp * sae$Narea))
}
\seealso{
  \code{\link{summary.weights}}, \code{\link{plot.weights}}
}

