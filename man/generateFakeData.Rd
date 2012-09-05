\name{generateFakeData}
\alias{generateFakeData}
\title{Generate artificial dataset for demonstration and testing purposes.}
\usage{
  generateFakeData(M = 50, meanNarea = 1000, sW = 100,
    sB = 50, sBx = 0.5, samplingFraction = 0.1)
}
\arguments{
  \item{M}{number of areas.}

  \item{meanNarea}{mean number of population units per
  area.}

  \item{sW}{within area standard deviation.}

  \item{sB}{between area standard deviation.}

  \item{sBx}{random slope standard deviation.}

  \item{samplingFraction}{sampling fraction used to draw a
  random sample from the population units.}
}
\value{
  List containing sample (sam), population totals (Xpop),
  and true population means for four target variables (mY0,
  mY1, mY2, mY3).
}
\description{
  Generate artificial dataset for demonstration and testing
  purposes.
}

