detect.maxima <- function(x, window.size){
  require(zoo)
  if(window.size %% 2 == 0) stop("window size must be odd")
  middle <- ceiling(window.size / 2)
  xz <- as.zoo(x)
  out <- rollapply(x, window.size, function(x) which.max(x) == middle)
  out <- c(rep(NA, middle-1), out, rep(NA, middle-1))
  return(out)
}

detect.minima <- function(x, window.size){
  require(zoo)
  if(window.size %% 2 == 0) stop("window size must be odd")
  middle <- ceiling(window.size / 2)
  xz <- as.zoo(x)
  out <- rollapply(x, window.size, function(x) which.min(x) == middle)
  out <- c(rep(NA, middle-1), out, rep(NA, middle-1))
  return(out)
}

x <- sin(seq(0,25,0.1))
plot(x)
abline(v = which(detect.maxima(x, 7)), col = 'red', lty = 'dashed')
abline(v = which(detect.minima(x, 7)), col = 'blue', lty = 'dashed')
