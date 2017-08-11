get_max_cc <- function(ts1, ts2, plot = F, ...){
  temp <- ccf(ts1, ts2, plot = plot, ...)
  max.indx <- which.max(temp$acf)
  return(list(temp$acf[max.indx], temp$lag[max.indx]))
}
