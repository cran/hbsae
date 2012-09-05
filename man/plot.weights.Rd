\name{plot.weights}
\alias{plot.weights}
\title{Plot method for objects of class \code{weights}.}
\usage{
  \method{plot}{weights} (x, log = FALSE,
    main = "Distribution of weights",
    xlab = if (log) "log(weight)" else "weight",
    ylab = "frequency", col = "cyan", ...)
}
\arguments{
  \item{x}{object of class \code{weights} as returned by
  function \code{\link{uweights}}.}

  \item{log}{whether to log-tranform the weights.}

  \item{main}{main title of plot.}

  \item{xlab}{x-axis label.}

  \item{ylab}{y-axis label.}

  \item{col}{colour.}

  \item{...}{additional arguments passed to
  \code{\link{hist}}.}
}
\description{
  Plot method for objects of class \code{weights}.
}
\seealso{
  \code{\link{uweights}}, \code{\link{summary.weights}}
}

