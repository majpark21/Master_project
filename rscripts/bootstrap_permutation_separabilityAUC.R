
sep.meas.along.time <- function(data1, data2, time.col, measure.col){
  timev <- unique(data1[, get(time.col)])
  if(!(identical(unique(data2[, get(time.col)]), timev))) stop("Time vectors must be identical between the two data")
  out <- separability.measures(data1[get(time.col)==timev[1], get(measure.col)], data2[get(time.col)==timev[1], get(measure.col)])
  for(t in timev[2:length(timev)]){
    out <- rbind(out, separability.measures(data1[RealTime==t, get(measure.col)], data2[RealTime==t, get(measure.col)]))
  }
  out <- cbind(timev, out)
  return(out)
}


one.permutation.auc <- function(x, y, metric){
  n <- nrow(x)
  m <- nrow(y)
  temp <- rbind(x, y)
  samp.traj <- sample(1:nrow(temp), size = n, replace = FALSE)
  x.resamp <- temp[samp.traj, ]
  y.resamp <- temp[setdiff(1:nrow(temp), samp.traj), ]
  
  seps <- sapply(1:ncol(x), function(j) separability.measures(x.resamp[, j], y.resamp[, j]))
  return(sum(unlist(seps[metric, ])))
}

permutation.auc <- function(x, y, n, metric = "jm"){
  # x,y: two matrices representing time series, row: trajectory; col: time
  # n: number of permutations
  # metric: one of "jm", "bh", "div", "tdiv", "ks"
  if(ncol(x) != ncol(y)) stop("x and y must have same number of columns")
  return(replicate(n, one.permutation.auc(x,y,metric)))
}

wrap_perm <- function(x, y, measure, n, na.fill){
  a <- as.matrix(dcast(x, Condition + Label ~ RealTime, value.var = measure)[,-c(1,2)])
  b <- as.matrix(dcast(y, Condition + Label ~ RealTime, value.var = measure)[,-c(1,2)])
  a[which(is.na(a))] <- na.fill
  b[which(is.na(b))] <- na.fill
  return(permutation.auc(a, b, n))
}

# --------------

one.bootstrap.auc.percol <- function(x, y, metric){
  samp.col <- sample(1:ncol(x), size = ncol(x), replace = TRUE)
  x.resamp <- x[, samp.col]
  y.resamp <- y[, samp.col]
  seps <- sapply(1:ncol(x), function(j) separability.measures(x.resamp[, j], y.resamp[, j]))
  return(sum(unlist(seps[metric, ])))
}

bootstrap.auc.percol <- function(x, y, B, metric = "jm"){
  # x,y: two matrices representing time series, row: trajectory; col: time
  # B: number of boostraps
  # metric: one of "jm", "bh", "div", "tdiv", "ks"
  if(ncol(x) != ncol(y)) stop("x and y must have same number of columns")
  return(replicate(B, one.bootstrap.auc.percol(x,y,metric)))
}

wrap_bootcol <- function(x, y, measure, n, na.fill){
  a <- as.matrix(dcast(x, Condition + Label ~ RealTime, value.var = measure)[,-c(1,2)])
  b <- as.matrix(dcast(y, Condition + Label ~ RealTime, value.var = measure)[,-c(1,2)])
  a[which(is.na(a))] <- na.fill
  b[which(is.na(b))] <- na.fill
  return(bootstrap.auc.percol(a, b, n))
}

# --------------

one.bootstrap.auc.perrow <- function(x, y, metric){
  samp.rowx <- sample(1:nrow(x), size = nrow(x), replace = TRUE)
  samp.rowy <- sample(1:nrow(y), size = nrow(y), replace = TRUE)
  x.resamp <- x[samp.rowx, ]
  y.resamp <- y[samp.rowy, ]
  seps <- sapply(1:ncol(x), function(j) separability.measures(x.resamp[, j], y.resamp[, j]))
  return(sum(unlist(seps[metric, ])))
}

bootstrap.auc.perrow <- function(x, y, B, metric = "jm"){
  # x,y: two matrices representing time series, row: trajectory; col: time
  # B: number of boostraps
  # metric: one of "jm", "bh", "div", "tdiv", "ks"
  if(ncol(x) != ncol(y)) stop("x and y must have same number of columns")
  return(replicate(B, one.bootstrap.auc.perrow(x,y,metric)))
}

wrap_bootrow <- function(x, y, measure, n, na.fill){
  a <- as.matrix(dcast(x, Condition + Label ~ RealTime, value.var = measure)[,-c(1,2)])
  b <- as.matrix(dcast(y, Condition + Label ~ RealTime, value.var = measure)[,-c(1,2)])
  a[which(is.na(a))] <- na.fill
  b[which(is.na(b))] <- na.fill
  return(bootstrap.auc.perrow(a, b, n))
}