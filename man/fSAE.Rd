\name{fSAE}
\alias{fSAE}
\title{Fit a linear model with random area effects and compute small area estimates.}
\usage{
  fSAE(formula, data, area = NULL, popdata = NULL,
    type = "unit", model.direct = NULL,
    formula.area = NULL, contrasts.arg = NULL,
    remove.redundant = TRUE, redundancy.tol = 1e-07,
    sparse = FALSE, ...)
}
\arguments{
  \item{formula}{model formula, indicating response
  variable and covariates.}

  \item{data}{unit-level data frame containing all
  variables used in \code{formula}, \code{area} and
  \code{formula.area} arguments.  These variables should
  not contain missing values.}

  \item{area}{name of area indicator variable in
  \code{data}; if \code{NULL}, no random effects are used
  in the model.}

  \item{popdata}{data frame or matrix containing area
  population totals for all covariates. The rows should
  correspond to areas for which estimates are required.
  Column names should include those produced by
  \code{model.matrix(formula, data, contrasts.arg)}, up to
  permutations of the names in interactions.  A column
  named '(Intercept)' is required and should contain the
  area population sizes.  If \code{popdata} is \code{NULL},
  only the model fit will be returned.}

  \item{type}{type of small area estimates: "direct" for
  survey regression, "area" for area-level model, "unit"
  for unit-level model estimates.  If \code{type} is "data"
  then only the data including the model matrix and
  population means are returned.}

  \item{model.direct}{if type="area", this argument can be
  used to specify by means of a formula the covariates to
  use for the computation of the initial survey regression
  estimates.  If unspecified, the covariates specified by
  \code{formula} are used both at the unit-level (for the
  initial estimates) and at the area-level (for the
  area-level model estimates).}

  \item{formula.area}{if type="unit", this is an optional
  formula specifying covariates that should be used at the
  area level.  These covariates should be available in
  \code{popdata}.}

  \item{contrasts.arg}{list for specification of contrasts
  for factor variables. Passed to \code{model.matrix}.}

  \item{remove.redundant}{if \code{TRUE} redundant columns
  in the design matrix are removed. A warning is issued if
  the same redundancy does not show also in the
  corresponding population totals. In the case of the area
  level model there may still be redundancy at the area
  level.}

  \item{redundancy.tol}{tolerance for detecting linear
  dependencies among the columns of the design matrix. Also
  used as tolerance in the check whether the design matrix
  reduncancy is shared by the population totals.}

  \item{sparse}{if \code{TRUE} \code{sparse.model.matrix}
  (package \code{Matrix}) is used to compute the covariate
  design matrix. This can be efficient for large datasets
  and a model containing categorical variables with many
  categories.}

  \item{...}{additional arguments passed to
  \code{\link{fSAE.Unit}} or \code{\link{fSurvReg}}.}
}
\value{
  An object of class \code{sae} containing the small area
  estimates, their MSEs, and the model fit. If \code{type}
  is "data" a list containing the model matrix, response
  vector, area indicator, area population sizes and matrix
  of population means is returned.
}
\description{
  This function prepares the (unit-level) input data and
  calls one of the lower level functions
  \code{\link{fSurvReg}}, \code{\link{fSAE.Area}} or
  \code{\link{fSAE.Unit}} to compute survey regression,
  area-level model or unit-level model small area
  estimates. Area-level model estimates are computed by
  first computing survey regression estimates and using
  these as input for \code{\link{fSAE.Area}}.
}
\examples{
d <- generateFakeData()

# model fitting only
(fit <- fSAE(y0 ~ x + area2, data=d$sam, area="area"))

# model fitting and small area estimation, unit-level model
saeHB <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, silent=TRUE)
saeHB  # print a summary
EST(saeHB)  # small area estimates
SE(saeHB)  # error estimates
str(saeHB)
plot(saeHB, list(est=d$mY0), CI=2)  # compare to true population means

# unit-level model with REML model-fit instead of Bayesian approach
saeREML <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, method="REML", silent=TRUE)
plot(saeHB, saeREML)  # compare

# basic area-level model
saeA <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, type="area")
plot(saeHB, saeA)

# SAE estimates based on a linear unit-level model without area effects
saeL <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, method="synthetic")
plot(saeHB, saeL)

# model-based estimation of overall population mean without area effects
est.global <- fSAE(y0 ~ x + area2, data=d$sam, area=NULL, popdata=colSums(d$Xpop), method="synthetic")
EST(est.global); MSE(est.global)

# no model fitting or estimation, but return design matrix, variable of interest,
#   area indicator, area population sizes and matrix of population means
dat <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, type="data")
str(dat)
}
\seealso{
  \code{\link{sae-class}}
}

