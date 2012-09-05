\name{sae-class}
\alias{cAIC}
\alias{coef.sae}
\alias{COR}
\alias{COV}
\alias{CV}
\alias{EST}
\alias{fitted.sae}
\alias{MSE}
\alias{R2}
\alias{raneff}
\alias{raneff.se}
\alias{relSE}
\alias{residuals.sae}
\alias{sae-class}
\alias{SE}
\alias{se2}
\alias{sv2}
\alias{synthetic}
\alias{vcov.sae}
\alias{wDirect}
\title{S3 class for the fitted model and SAE outcomes.}
\description{
  Functions \code{\link{fSAE}}, \code{\link{fSurvReg}},
  \code{\link{fSAE.Area}} and \code{\link{fSAE.Unit}}
  return an object of class \code{sae}. It contains
  information on the model fit as well as the small area
  estimates, error estimates and a few model selection
  measures. The functions listed below extract the main
  components from an object of class \code{sae}. \describe{
  \item{\code{EST(x, type="sae")}}{return the vector of
  small area estimates of \code{sae} object x.
  Alternatively, with \code{type} "coef" or "raneff" fixed
  or random effect estimates are returned.}
  \item{\code{MSE(x, type="sae")}}{return the vector of
  mean squared errors of \code{sae} object x.
  Alternatively, with \code{type} "coef" or "raneff" MSEs
  of fixed or random effects are returned.}
  \item{\code{SE(x, type="sae")}}{extract standard errors,
  i.e. square roots of the MSEs.} \item{\code{relSE(x,
  type="sae")}}{extract relative standard errors.}
  \item{\code{COV(x)}}{extract the covariance matrix for
  the small area estimates.} \item{\code{COR(x)}}{extract
  the correlation matrix for the small area estimates.}
  \item{\code{coef(x)}}{\code{coef} method for \code{sae}
  objects; returns vector of fixed effects.}
  \item{\code{vcov(x)}}{\code{vcov} method for \code{sae}
  objects; returns covariance matrix for fixed effects.}
  \item{\code{raneff(x, pop)}}{return vector of random
  effects. If \code{pop=TRUE} returns a vector for
  predicted areas (zero for out-of-sample areas), otherwise
  a vector for sampled areas.} \item{\code{raneff.se(x,
  pop)}}{return vector of standard errors for random
  effects.} \item{\code{residuals(x)}}{\code{residuals}
  method for \code{sae} objects; returns a vector of
  residuals.} \item{\code{fitted(x)}}{\code{fitted} method
  for \code{sae} objects; returns a vector of fitted
  values.} \item{\code{se2(x)}}{extracts within-area
  variance estimate.} \item{\code{sv2(x)}}{extracts
  between-area variance estimate.} \item{\code{wDirect(x,
  pop)}}{extract vector of weights of the survey regression
  components in the small area estimates. If
  \code{pop=TRUE} returns a vector for predicted areas
  (zero for out-of-sample areas), otherwise a vector for
  sampled areas.} \item{\code{synthetic(x)}}{extract vector
  of synthetic estimates.} \item{\code{CV(x)}}{extract
  leave-one-out cross-validation measure.}
  \item{\code{cAIC(x)}}{extract conditional AIC measure.}
  \item{\code{R2(x)}}{extract unit-level R-squared
  goodness-of-fit measure.} } Other components include
  \describe{ \item{\code{relErrM,relErrV}}{relative
  numerical integration errors in estimates and MSEs, for
  \code{method} "HB".} }
}
\examples{
d <- generateFakeData()

# compute small area estimates
sae <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop)

coef(sae)  # fixed effects
raneff(sae)  # random effects
sv2(sae)  # between-area variance
se2(sae)  # within-area variance
cAIC(sae)  # conditional AIC
}

