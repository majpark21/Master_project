#####################
# Analysis of trend #
#####################

#' MPFeatTrend_extrema
#'
#' Isolate local extrema and perform a linear regression.
#' Aimed at at giving an estimation of long-term trend in multi-peak signals.
#'
#' @param x a numerical vector. Trajectory with multiple peaks.
#' @param window.size integer, width of window for maxima detection.
#'  See ?detect.peak.
#' @param what character, one of "minima" or "maxima". Define whether regression
#' should be performed on maxima or minima of peaks.
#' @param robust logical. If TRUE perform robusts linear regression.
#' Median-Based Linear Models are used, see ?mblm::mblm
#'
#' @return A list of 4.
#' \itemize{
#'  \item $trend, slope of the linear model
#'  \item $model, the full linear model
#'  \item $extremes, x extrema values and their index
#'  \item $type, type of extrema, "mini" or "maxi"
#'  }
#' @export
#'
#' @examples
#' x <- sim_phase_shifted_with_fixed_trend(n = 1, noise = 0, slope = 0.1)
#' x.trend <- MPFeatTrend_extrema(x = x$value, window.size = 5, what = "maxi")
#' plot(x$value, type = "b")
#' abline(x.trend$model, col = "blue", lwd = 2, lty = "dashed")
#'
MPFeatTrend_extrema <- function(x, window.size, what, robust = FALSE){
  require(mblm)
  # Extract local extrema
  extr.pos <- which(detect.peak(x, window.size, what))
  extr.x <- x[extr.pos]
  names(extr.x) <- extr.pos

  # Perform fitting
  if(robust){
    fit <- mblm(extr.x ~ extr.pos)
    return(list(trend=fit$coefficients[2], model=fit, extremes=extr.x, type=what))
  } else {
    fit <- lm(extr.x ~ extr.pos)
    return(list(trend=fit$coefficients[2], model=fit, extremes=extr.x, type=what))
  }
}



