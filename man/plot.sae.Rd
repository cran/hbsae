\name{plot.sae}
\alias{plot.sae}
\title{Plot method for objects of class sae.}
\usage{
  \method{plot}{sae} (..., sort.by = NULL,
    decreasing = FALSE, index = NULL, maxrows = 50,
    maxcols = 6, labels = NULL, type = "sae",
    ylab = if (type == "coef") "coef" else "area",
    main = switch(type, coef = "coefficients", raneff = "random effects", "SAE estimates"),
    offset = 0.2, cex.var = 0.9, mar = c(0, 1, 5, 0),
    CI = 1)
}
\arguments{
  \item{...}{list of \code{sae} objects (or lists with at
  least a component "est" (and optionally "mse")).  Other
  components not matching other arguments are passed to the
  lower-level plotting functions.}

  \item{sort.by}{vector by which to sort the vectors to be
  plotted; defaults to population size vector of the first
  sae object.}

  \item{decreasing}{if \code{TRUE}, sort in decreasing
  order (default).}

  \item{index}{vector of indices of the selected areas to
  be plotted.}

  \item{maxrows}{maximum number of rows in a column.}

  \item{maxcols}{maximum number of columns of estimates on
  a page.}

  \item{labels}{optional vector of lables; default is the
  vector of area names of the first sae object.}

  \item{type}{"sae" for small area estimates (default),
  "coef" for coefficients, "raneff" for random effects.}

  \item{ylab}{label for y-axis.}

  \item{main}{main title.}

  \item{offset}{space used between plots of multiple
  estimates for the same area.}

  \item{cex.var}{the fontsize of the variable names,
  default=0.9.}

  \item{mar}{a numerical vector of the form c(bottom, left,
  top, right) which gives the number of lines of margin to
  be specified on the four sides of the plot. The default
  is c(0,1,5,0).}

  \item{CI}{confidence interval, default is 1, which will
  plot plus and minus 1 standard errors. If CI=2, plots
  plus and minus 2 standard errors instead.}
}
\description{
  This function calls \code{coefplot} from package
  \pkg{arm} to plot small area estimates and standard
  errors. Multiple sets of estimates can be compared using
  this plot. The default ordering of the estimates is by
  their area population sizes.
}

