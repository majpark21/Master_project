#####################
# Analysis of trend #
#####################

#' MPFeatTrend_extrema
#'
#' Isolate local extrema and perform a linear regression.
#' Aimed at at giving an estimation of long-term trend in multi-peak signals.
#' @param x a numerical vector. Trajectory with multiple peaks.
#' @param window.size integer, width of window for maxima detection.
#'  See ?detect.peak.
#' @param what character, one of "minima" or "maxima". Define whether regression
#' should be performed on maxima or minima of peaks.
#' @param robust logical. If TRUE perform robusts linear regression.
#' Median-Based Linear Models are used, see ?mblm::mblm
#'
#' @return A list of 2. $trend comprises the fit object
#' @export
#'
#' @examples
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
    return(list(trend=fit$coefficients[2], trend.complete=fit, extremes=extr.x))
  } else {
    fit <- lm(extr.x ~ extr.pos)
    return(list(trend=fit$coefficients[2], trend.complete=fit, extremes=extr.x))
  }
}
