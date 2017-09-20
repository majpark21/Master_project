#' Find Full Width at Half maximum by spline interpolation
#'
#' @param x numeric, typically time vector
#' @param y numeric, typically measured variable
#' @param n numeric, number of points returned by interpolation. Increase improves accuracy.
#' @param method one of "minimum" or "walk". Walk method is highly recommended as it avoids many pitfalls 
#' of the minimum method. Minimum simply pick the points which values are the closest to the half of the max value.
#' Walk method walks left and right from the maximum point, and stops whenever half max value is crossed.
#'
#' @return A list of 3: left and right represent the value of x at which half of the maximum in y is reached.
#' fwhm is the difference between left and right
get.FWHM <- function(x, y, n = 30*length(x), method = "walk"){
  if(!method %in% c("walk", "minimum")) stop("Method must be one of c('walk','minimum')")
  
  # Do NOT use interpolated data as maximum estimation
  maxi <- max(y)
  tmaxi <- x[which.max(y)]
  # Spline interpolated and data
  fit <- mySpline(x, y, n)
  
  if(method == "walk"){
    thresh <- maxi/2
    start <- which(fit$x == tmaxi)[1]
    # Right walk
    right <- NA
    for(i in start:length(fit$x)){
      if(fit$y[i] < thresh){
        # Once threshold is passed, choose the value between i and i-1 that is the closest to it
        right <- ifelse(abs(thresh - fit$y[i-1]) < abs(thresh - fit$y[i]), fit$x[i-1], fit$x[i])
        break
      }
    }
    # Left walk
    left <- NA
    for(i in seq(start, 1, -1)){
      if(fit$y[i] < thresh){
        left <- ifelse(abs(thresh - fit$y[i+1]) < abs(thresh - fit$y[i]), fit$x[i+1], fit$x[i])
        break
      }
    }
    
  } else if(method == "minimum"){
    left <- fit$x[which.min(abs(fit$y[fit$x < tmaxi] - (maxi/2)))]
    right <- fit$x[which.max(which(fit$x <= tmaxi)) + which.min(abs(fit$y[fit$x > tmaxi] - (maxi/2)))]
  }
  
  return(list(fwhm = right - left, left = left, right = right))
}


# Returns the interpolated points along with the data in a single trajectory
mySpline <- function(x, y , n){
  fit <- spline(x, y, n)
  x2 <- c(x, fit$x)
  y2 <- c(y, fit$y)
  y2 <- y2[order(x2)]
  x2 <- x2[order(x2)]
  return(list(x = x2, y = y2))
}