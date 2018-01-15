#####################
# Analysis of trend #
#####################

#' MPFeatTrend_extrema
#'
#' Isolate local extrema and perform a linear regression.
#' Aims at at giving an estimation of long-term trend in multi-peak signals.
#'
#' @param x a numerical vector.
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
#' @seealso MPFeatTrend_rollmean
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


#' MPFeatTrend_rollmean
#'
#' Smooth by rolling mean and perform linear regression.
#' Rolling mean is extended at extremeties by linear extrapolation, see ?rollex.
#' Aims at at giving an estimation of long-term trend in multi-peak signals.
#'
#' @param x a numerical vector.
#' @param window.size integer, width of window for rolling mean.
#' @param robust logical. If TRUE perform robusts linear regression.
#' Median-Based Linear Models are used, see ?mblm::mblm
#'
#' @return A list of 2: $model, the full linear model;
#' $trend, slope of the linear model.
#' @export
#' @seealso MPFeatTrend_extrema
#' @examples
#' width <- 30
#' x <- sim_phase_shifted_with_fixed_trend(n = 1, noise = 0, slope = 0.1)
#' # Rolling mean extended by linear extrapolation.
#' # This is also performed in MPFeatTrend_rollmean.
#' x.roll <- rollex(x$value, width)
#' plot(x$value, type = "b")
#' lines(x.roll, col ="green", lwd = 2, lty = "dashed")
#' x.trend <- MPFeatTrend_rollmean(x = x$value, window.size = width)
#' abline(x.trend$model, col = "blue", lwd = 2, lty = "dashed")
#'
MPFeatTrend_rollmean <- function(x, window.size, robust = FALSE){
  require(mblm)
  x.roll <- rollex(x, window.size)
  if(robust){
    fit <- mblm(x ~ seq_along(x))
    return(list(trend=fit$coefficients[2], model = fit))
  } else {
    fit <- lm(x ~ seq_along(x))
    return(list(trend=fit$coefficients[2], model = fit))
  }
}






