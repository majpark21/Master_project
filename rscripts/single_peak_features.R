FeatMaxAmplitude <- function(y, basal = 1){
  # Shift the data to have basal at 0, this is important so that maximum peak is measured in amplitude instead of absolute value
  y <- y - basal
  return(list(max = max(y), time.max = which.max(y)))
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


#' Find Full Width at Half maximum by spline interpolation
#'
#' @param x numeric, typically time vector
#' @param y numeric, typically measured variable, data must be normalized!
#' @param n numeric, number of points returned by interpolation. Increase improves accuracy.
#' @param method one of "minimum" or "walk". Walk method is highly recommended as it avoids many pitfalls 
#' of the minimum method. Minimum simply pick the points which values are the closest to the half of the max value.
#' Walk method walks left and right from the maximum point, and stops whenever half max value is crossed.
#' Data must be normalized, such that "basal" represents the activity out of excitation.
#' @return A list of 3: left and right represent the value of x at which half of the maximum in y is reached.
#' fwhm is the difference between left and right
FeatFWHM <- function(y, x = seq_along(y), n = 30*length(x), method = "walk", basal = 1){
  if(!method %in% c("walk", "minimum")) stop("Method must be one of c('walk','minimum')")
  
  # Shift the data to have basal at 0, this is important so that maximum peak is measured in amplitude instead of absolute value
  y <- y - basal
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

traj <- Cora[.("P5-I100","11_9"), Ratio.norm]
plot(seq_along(traj),traj, type = "b") 
lines(mySpline(seq_along(traj), traj, 30*length(traj))$x, mySpline(seq_along(traj), traj, 30*length(traj))$y, lty = "dashed", col = "red")
traj.fwhm <- FeatFWHM(traj)
abline(v=c(traj.fwhm$left, traj.fwhm$right), col = "blue", lty = "dashed")

FeatHalfMaxDec <- function(y, ...){
  # Enforces fit to pass through maximum
  right <- FeatFWHM(y = y, ...)$right
  tmaxi <- which.max(y)
  # Use the point closest to half-max as additional point for regression (if it is interpolated then it is not an integer)
  if(right %% 1 !=0){
    splinetraj <- mySpline(x = seq_along(y), y, 30*length(y))
    newpoint <- splinetraj$y[splinetraj$x == right]
    ytrim <- c(y[tmaxi:right], newpoint) - max(y)
    t.after.max <- seq_along(ytrim) - 1
    # Replace last time by the time of interpolation
    t.after.max[length(t.after.max)] <- right - tmaxi
  } else {
    ytrim <- y[tmaxi:right] - max(y)
    t.after.max <- seq_along(ytrim) - 1
  }
  fit <- lm(ytrim ~ t.after.max + 0)
  return(coef(fit)[1])
}

FeatExpDec <- function(){}

FeatLagGrow <- function(){}

FeatHalfMaxGrow <- function(y, ...){
  left <- FeatFWHM(y = y, ...)$left
  tmaxi <- which.max(y)
  # Use the point closest to half-max as additional point for regression (if it is interpolated then it is not an integer)
  if(left %% 1 !=0){
    splinetraj <- mySpline(x = seq_along(y), y, 30*length(y))
    newpoint <- splinetraj$y[splinetraj$x == left]
    ytrim <- c(newpoint, y[(left+1):(tmaxi+1)]) - max(y)
    # Set interpolation time to 0, and shift all other times
    time.reg <- c(0,   seq(1, (length(ytrim)-1))-(1-left%%1) )
  } else {
    ytrim <- y[left:tmaxi] - max(y)
    time.reg <- seq_along(ytrim) - 1
  }
  fit <- lm(ytrim ~ time.reg)
  return(coef(fit)[2])
}


FeatAllFeat <- function(){}
