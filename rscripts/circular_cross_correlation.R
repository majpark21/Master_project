circular.cc <- function(ts1, ts2){
  # Cross-correlation with circular boundaries
  if(length(ts1) != length(ts2)) stop("Both vectors should be of the same length")
  require(wavethresh)
  out <- rep(NA, length(ts1))
  lags <- 0:length(ts1)
  for(i in 1:length(lags)){
    out[i] <- cor(x = ts1, y = guyrot(ts2, lags[i]), method = "pearson") 
  }
  return(out)
}