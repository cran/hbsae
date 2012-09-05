\name{bench}
\alias{bench}
\title{Benchmark small area estimates.}
\usage{
  bench(x, R, rhs, mseMethod = "no")
}
\arguments{
  \item{x}{sae object to be benchmarked.}

  \item{R}{restriction matrix, r x M matrix where r is the
  number of restrictions and M the number of areas; default
  is a single constraint on the population total.  Note
  that R acts on the vector of area population totals, not
  the vector of means.}

  \item{rhs}{r-vector of benchmark totals corresponding to
  the restrictions in the rows of R.}

  \item{mseMethod}{if "no", MSEs are not updated, if
  "exact", constraints are treated as exact identities, and
  if "model", the squared differences between original and
  benchmarked estimates are added to the MSEs.}
}
\value{
  An object of class \code{sae} with adjusted estimates.
}
\description{
  Benchmark small area estimates to agree with given totals
  at aggregate levels.
}
\examples{
d <- generateFakeData()

# compute small area estimates
sae <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop)

# calibrate to overall population total
sae.c <- bench(sae, rhs=sum(d$mY0*sae$Narea))
plot(sae, sae.c)
}
\references{
  Y. You, J.N.K. Rao and P. Dick (2004). Benchmarking
  Hierarchical Bayes Small Area Estimators in the Canadian
  Census Undercoverage Estimation. Statistics in Transition
  6(5), 631-640.
}
\seealso{
  \code{\link{sae-class}}
}

