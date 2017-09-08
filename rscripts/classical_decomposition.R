#' Decompose a time series into trend, seasonal and remainder component
#' 
#' The model is assumed to be additive, i.e. the time series is of the form:
#' y(t) = trend(t) + seasonal(t) + remainder(t)
#'
#' Trends is determined using a rolling mean with window width corresponding to 
#' the frequency, extremities are padded with linear interpolation from the first and last 2 measures.
#' 
#' Season is determined by using the mean or the median of all points corresponding to a 
#' given part of a cycle, see robust.
#' 
#' @param ts: a time series or numerical vector
#' @param frequency: number of time points in one cycle
#' @param robust: if T, the median of all points of a cycle part are used to determined seasonal component;
#'                if F, the mean is used instead.
#'
#' @return a 3xn numeric matrix, with column trend, season and remainder
#' 
classical.decomposition <- function(ts, frequency, robust = T){
  if(!exists("rollex", mode="function")) source("./overlap_clipping.R")
  # 1) Get trend with extended rolling mean
  trend <- rollex(ts, frequency)
  # 2) Detrend
  detrend <- ts - trend
  # 3) Average every time point for each oscillation
  seasonal <- rep(NA, length(ts))
  if(robust){
    for(i in 1:frequency){
      index <- seq(i, length(ts), by = frequency)
      seasonal[index] <- median(detrend[index])
    } 
  } else {
    for(i in 1:frequency){
      index <- seq(i, length(ts), by = frequency)
      seasonal[index] <- mean(detrend[index])
    }
  }
  # 4) Adjust to zero
  seasonal <- scale(seasonal, center = T, scale = F)
  # 5) Remainder
  remainder <- ts - (trend + seasonal)
  out <- cbind(trend, seasonal, remainder)
  colnames(out) <- c("trend", "seasonal", "remainder")
  return(out)
}


plot.decomposition <- function(dat, dec.mat, main = "Time series decomposition"){
  xax <- seq_along(dat)
  layout(matrix(1:4, ncol = 1), widths = 1, heights = c(1.6,1,1,1.5), respect = F)
  par(mar = c(0, 5.1, 4.1, 2.1))
  plot(xax, dat, main = main, type = "l", xaxt = 'n', xlab = "", ylab = "Data", cex.lab = 1.4, cex.main = 1.5)
  par(mar = c(0, 5.1, 0, 2.1))
  plot(xax, dec.mat[,"seasonal"], type = "l", xaxt = 'n', xlab = "", ylab = "Seasonal", cex.lab = 1.4)
  par(mar = c(0, 5.1, 0, 2.1))
  plot(xax, dec.mat[,"trend"], type = "l", xaxt = 'n', xlab = "", ylab = "Trend", cex.lab = 1.4)
  par(mar = c(4.1, 5.1, 0, 2.1))
  plot(xax, dec.mat[,"remainder"], type = "h", ylab = "Remainder", xlab = "Time", cex.lab = 1.4)
  abline(h=0, lty="dashed")
}
