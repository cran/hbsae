
#' Print method for objects of class sae.
#'
#' @export
#' @method print sae
#' @param x object of class \code{sae}.
#' @param digits number of digits to display.
#' @param ... additional arguments passed to \code{print.default}.
print.sae <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  varcomp <- c("posterior mean sv2" = x$sv2hat, "REML estimate sv2" = x$sv2REML, "posterior mean sv2/se2" = x$lambdahat, "REML estimate sv2/se2" = x$lambdaREML)
  if (!is.null(x$call)) cat("Call: ", deparse(x$call), "\n\n", sep = "")
  cat("Summaries:\n")
  summary_noNA <- function(x) summary.default(x)[seq_len(6L)]  # remove NA column, if any
  summ <- rbind("gamma weights" = summary_noNA(wDirect(x)), "random effects" = summary_noNA(raneff(x)))
  aggrs <- NULL
  if (!is.null(x$est)) {
    summ <- rbind(estimates = summary_noNA(x$est), "standard errors" = summary_noNA(sqrt(x$mse)), "rel.std.err." = summary_noNA(relSE(x)), summ)
    if (x$M > 1L && !is.null(x$Narea))
      aggrs <- list(mean = sum(x$Narea * x$est)/sum(x$Narea), total = sum(x$Narea * x$est))
  }
  print.default(format(summ, digits = digits), print.gap = 2, quote = FALSE, ...)
  if (x$benchmarked) cat("\nEstimates have been benchmarked.\n")
  if (x$aggregated) cat("\nEstimates have been aggregated.\n")
  if (!is.null(aggrs)) {
    cat("\nPopulation aggregates:\n")
    print.default(format(aggrs, digits = digits), print.gap = 2, quote = FALSE, ...)
  }
  varcomp <- c(varcomp, "se2 estimate" = se2(x))
  if (length(varcomp)) {
    cat("\n")
    print.default(format(varcomp, digits = digits), print.gap = 2, quote = FALSE, ...)
  }
  if (length(coef(x))) {
    cat("\nCoefficients:\n")
    coefs <- coef(x)
    if (length(vcov(x))) {
      coefs <- cbind(coefs, sqrt(diag(vcov(x))), coefs/sqrt(diag(vcov(x))))
      dimnames(coefs)[[2L]] <- c("Estimate", "Std. Error", "t value")
    }
    printCoefmat(coefs, digits = digits - 1L)
  }
  
  modsel <- list(p.eff = x$p.eff, cAIC = cAIC(x))
  if (length(x$CV))
    modsel <- c(list(CV = x$CV), modsel)
  cat("\nModel selection measures:\n")
  print.default(format(modsel, digits = digits), print.gap = 2, quote = FALSE, ...)
  if (length(x$R2))
    cat("R-squared:", round(x$R2, 5L), "\n")
  if (!is.null(x$low.fitted) && x$low.fitted > 0)
    cat("Fraction of fitted values below lowest data value:", round(x$low.fitted, 5L), "\n")
  if (!is.null(x$high.fitted) && x$high.fitted > 0)
    cat("Fraction of fitted values above highest data value:", round(x$high.fitted, 5L), "\n")
}


#' Plot method for objects of class sae.
#'
#' This function plots small area estimates with error bars.
#' Multiple sets of estimates can be compared. The default ordering of the estimates
#' is by their area population sizes.
#' This method uses a plot function that is adapted from function
#' \code{coefplot.default} of package \pkg{arm}.
#'
#' @export
#' @method plot sae
#' @param ... \code{sae} objects, \code{dc_summary} objects (output by the \code{summary} method for
#'  simulation objects of package \pkg{mcmcsae}), or lists.
#'  The first object must be a \code{sae} object.
#'  In case of a list the components
#'  used are those with name \code{est} for point estimates, \code{se}
#'  for standard error based intervals or \code{lower} and \code{upper} for
#'  custom intervals. Instead of \code{dc_summary} objects matrix objects are
#'  also supported as long as they contain columns named "Mean" and "SD" as do
#'  \code{dc_summary} objects. Named parameters of other types that do not match any
#'  other argument names are passed to lower-level plot functions.
#' @param n.se number of standard errors below and above the point estimates
#'  to use for error bars. By default equal to 1. This only refers to the
#'  objects of class \code{dc_summary} and \code{sae}.
#' @param est.names labels to use in the legend for the components of the \code{...} argument
#' @param sort.by vector by which to sort the coefficients, referring to the first object passed.
#' @param decreasing if \code{TRUE}, sort in decreasing order (default).
#' @param index vector of names or indices of the selected areas to be plotted.
#' @param maxrows maximum number of rows in a column.
#' @param maxcols maximum number of columns of estimates on a page.
#' @param type "sae" for small area estimates (default), "coef" for
#'  coefficients, "raneff" for random effects.
#' @param offset space used between plots of multiple estimates for the same
#'  area.
#' @param cex.var the fontsize of the variable names, default=0.8.
#' @param mar a numerical vector of the form c(bottom, left, top, right) which
#'  gives the number of lines of margin to be specified on the four sides of
#'  the plot.
plot.sae <- function(..., n.se=1, est.names, sort.by=NULL, decreasing=FALSE,
                      index=NULL, maxrows=50L, maxcols=6L, type="sae",
                      offset=0.1, cex.var=0.8, mar=c(0.1,2.1,5.1,0.1)) {

  dotargs <- list(...)
  toplot <- sapply(dotargs, function(obj) class(obj)[1L] %in% c("sae", "dc_summary", "matrix", "list"))
  grpar <- dotargs[!toplot]  # other graphical parameters
  if (!length(grpar)) grpar <- NULL

  x <- dotargs[toplot]
  if (!length(x)) stop("nothing to plot")

  # extract point estimates to plot
  ests <- lapply(x, function(obj) {
    switch(class(obj)[1L],
           sae = EST(obj, type=type),
           dc_summary=, matrix = obj[, "Mean"],
           list = obj$est,
           stop("unsupported class '", class(obj)[1L], "'")
    )
  })

  ints <- lapply(x, function(obj) {
    # for dc_summary and sae objects by default use mean +/- se intervals
    switch(class(obj)[1L],
           sae = cbind(EST(obj, type=type) - n.se * SE(obj, type=type), EST(obj, type=type) + n.se * SE(obj, type=type)),
           dc_summary=, matrix = cbind(obj[, "Mean"] - n.se * obj[, "SD"], obj[, "Mean"] + n.se * obj[, "SD"]),
           list = {
             if (!is.null(obj$se)) {
               if (is.null(obj$est)) stop("list with 'se' but no 'est' component")
               cbind(obj$est - n.se * obj$se, obj$est + n.se * obj$se)
             } else if (is.null(obj$lower)) {
               if (!is.null(obj$upper)) stop("'upper' specified without 'lower'")
               if (is.null(obj$est)) stop("object without 'est' and 'lower'")
               NULL
             } else if (is.null(obj$upper)) {
               stop("'lower' specified without 'upper'")
             } else {
               cbind(obj$lower, obj$upper)
             }
           }
    )
  })

  M <- length(ests[[1L]])
  lab <- names(ests[[1L]])
  if (is.null(lab)) lab <- as.character(seq_len(M))

  # make all objects compatible with the first, i.e. of the same length, if possible
  if (length(ests) > 1L) {
    for (i in 2:length(ests)) {
      obj.est <- ests[[i]]
      obj.int <- ints[[i]]
      if (is.null(names(obj.est))) {
        if (length(obj.est) != M) stop("unable to match components of input objects")
        names(ests[[i]]) <- lab
        if (!is.null(ints[[i]])) {
          if (nrow(obj.int) != M) stop("unable to match components of input objects")
          rownames(ints[[i]]) <- lab
        }
      } else {
        ind <- match(names(obj.est), lab)
        ests[[i]] <- NA_real_ * ests[[1L]]
        ests[[i]][ind[!is.na(ind)]] <- obj.est[!is.na(ind)]
        if (!is.null(ints[[i]])) {
          ints[[i]] <- NA_real_ * ints[[1L]]
          ints[[i]][ind[!is.na(ind)], ] <- obj.int[!is.na(ind), ]
        }
      }
    }
  }
  
  # allow plotting for a subset only
  if (is.null(index)) {
    o <- lab
  } else {
    if (is.character(index)) {
      if (anyNA(match(index, lab))) stop("some elements of index cannot be matched")
      o <- index
    } else {
      o <- lab[index]
    }
    M <- length(o)
  }
  
  if (is.null(sort.by)) {
    sort.by <- seq_along(ests[[1L]])
  }
  names(sort.by) <- lab

  o <- o[order(sort.by[o], decreasing=decreasing)]

  maxrows <- min(maxrows, M)
  cols <- ceiling(M/maxrows)
  pages <- ceiling(cols/maxcols)

  compute.xlim <- !("xlim" %in% names(grpar))

  oldpar <- par(no.readonly=TRUE)
  on.exit(par(oldpar))
  for (page in seq_len(pages)) {
    colrange <- ((page - 1L) * maxcols + 1L):min((page * maxcols), cols)
    plot.new()
    par(mfrow=c(1L, min(maxcols, cols)))
    for (co in colrange) {
      par(mfg=c(1L, co - (colrange[1L] - 1L)))  # force plot in the intended column
      rowrange <- ((co - 1L) * maxrows + 1L):min((co * maxrows), M)
      offsets <- offset * (seq_along(x) - 1L)  # - ((length(x) + 1) %/% 2))
      # first determine xlim, if not set manually
      if (compute.xlim) {
        xlim <- range(ests[[1L]][o[rowrange]][1L], na.rm=TRUE)
        for (i in seq_along(x)) {
          if (is.null(ints[[i]]))
            xlim <- range(xlim, ests[[i]][o[rowrange]], na.rm=TRUE)
          else
            xlim <- range(xlim, ints[[i]][o[rowrange], ], na.rm=TRUE)
        }
        grpar$xlim <- xlim
      }
      for (i in seq_along(x)) {
        if (is.null(ints[[i]]))
          interv <- NULL
        else
          interv <- ints[[i]][o[rowrange], ]
        if (i == 1L) {
          do.call(cplot, c(list(ests[[i]][o[rowrange]], interv, varnames=o[rowrange],
                                cex.var=cex.var, offset=offsets[i], mar=mar), grpar))
        } else {
          cplot(ests[[i]][o[rowrange]], interv, col.pts=i, add=TRUE, offset=offsets[i], mar=mar)
        }
      }
    }  # END for co
    if (!missing(est.names) || length(ests) > 1L) {  # add a legend
      par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
      plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n")
      if (missing(est.names))
        est.names <- paste0("est", seq_along(ests))
      dots.only <- sapply(x, function(obj) is.list(obj) && is.null(obj$se) && is.null(obj$lower))
      legend("topright", est.names, xpd=TRUE, horiz=TRUE, inset=c(0, 0),
             bty="n", pch=rep.int(20L, length(ests)), lty=ifelse(dots.only, 0, 1),
             col=seq_along(ests)
      )
    }
  }  # END for page
}

# Plot confidence or credible intervals for vectors of coefficients
# 
# This function is adapted from coefplot.default of package \pkg{arm}.
cplot <- function (coefs, intervals=NULL,
                   varnames=NULL, v.axis=TRUE, h.axis=TRUE,
                   cex.var=0.8, cex.pts=0.9, col.pts=1, pch.pts=20, 
                   var.las=2, xlab="", ylab="", main="",
                   mar=c(0.1,2.1,5.1,0.1), plot=TRUE,
                   add=FALSE, offset=0.1, ...) {
  
  if (is.list(coefs)) coefs <- unlist(coefs)
  m <- length(coefs)
  id <- seq_len(m)
  if (is.null(intervals)) intervals <- cbind(coefs, coefs)
  if (is.null(varnames)) {
    maxchar <- 0
  } else {
    maxchar <- max(sapply(varnames, nchar))
  }
  k <- 1/m
  oldmar <- par("mar")
  on.exit(par(mar=oldmar))
  mar[2L] <- max(oldmar[2L], trunc(mar[2L] + maxchar/10)) + 0.1
  par(mar=mar)
  if (plot) {
    if (add) {
      id <- id + offset
      points(coefs, id, pch=pch.pts, cex=cex.pts, col=col.pts)
      segments(intervals[, 1L], id, intervals[, 2L], id, lwd=1, col=col.pts)
    } else {
      plot(intervals, c(id + k, id - k), type="n",
           axes=FALSE, main=main, xlab=xlab, ylab=ylab, ...)
      if (h.axis) axis(3L)
      if (v.axis) {
        axis(2L, rev(id), varnames[rev(id)], las=var.las, 
             tck=FALSE, lty=0, cex.axis=cex.var)
      }
      abline(v=0, lty=2L)
      points(coefs, id, pch=pch.pts, cex=cex.pts, col=col.pts)
      segments(intervals[, 1L], id, intervals[, 2L], id, lwd=1, col=col.pts)
    }
  } else {  # do not plot
    plot(intervals, c(id + k, id - k), type="n",
         axes=FALSE, main="", xlab=xlab, ylab=ylab, ...)
  }
}


#' Summary method for objects of class \code{weights}.
#'
#' @export
#' @method summary weights
#' @param object object of class \code{weights} as returned by function \code{\link{uweights}}.
#' @param ... not used.
#' @seealso \code{\link{uweights}}, \code{\link{plot.weights}}
summary.weights <- function(object, ...) {
  deft <- sqrt(deff.weights(object))
  object <- unclass(object)
  cat("\n")
  print(summary(object))
  cat("\n")
  cat(sprintf("                       number:  %i \n", length(object) ))
  cat(sprintf("                          sum:  %.1f \n", sum(object) ))
  cat(sprintf("   number of negative weights:  %i \n", sum(object < 0) ))
  cat(sprintf("          number of 0 weights:  %i \n", sum(object == 0) ))
  cat(sprintf("number of weights > 0 and < 1:  %i \n", sum(object > 0 & object < 1) ))
  cat(sprintf("          number of weights 1:  %i \n", sum(object == 1) ))
  cat(sprintf("          sqrt(design-effect):  %.4f \n\n", deft ))
}

#' Plot method for objects of class \code{weights}.
#'
#' @export
#' @method plot weights
#' @param x object of class \code{weights} as returned by function \code{\link{uweights}}.
#' @param log whether to log-transform the weights.
#' @param breaks breaks argument of function \code{\link{hist}}. Default is "Scott".
#' @param main main title of plot.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param col colour.
#' @param ... additional arguments passed to \code{\link{hist}}.
#' @seealso \code{\link{uweights}}, \code{\link{summary.weights}}
plot.weights <- function(x, log=FALSE, breaks="Scott", main="Distribution of weights",
                         xlab=if (log) "log(weight)" else "weight", ylab="frequency", col="cyan", ...) {
  x <- unclass(x)
  if (log) x <- log(x)
  h <- hist(x, breaks=breaks, main=main, xlab=xlab, ylab=ylab, col=col, ...)
  boxplot(x, horizontal=TRUE, add=TRUE, at=-max(h$counts)/50, axes=TRUE)
}


