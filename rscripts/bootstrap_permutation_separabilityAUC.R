#' Separability measures along time
#'
#' Takes 2 data tables in long format and computes a set of separability indices at each time point.
#' @param data1 A data.table in long format.
#' @param data2 A data.table in long format.
#' @param time.col Character, name of the column with time variable
#' @param measure.col Character, name of the colume with measured variable
#' 
#' @note The function will returns lots of warnings in case of ties, because of KS statistics.
#' @return A data frame containing the separability indices at each time point.
#' 'timev': time vector
#' 'jm': Jeffries-Matusita, a bounded version of Bhattacharyya's distance
#' 'bh': Bhattacharyya's distance
#' 'div': Kullback-Leibler's divergence
#' 'tdiv': Transformed divergence
#' 'ks': Kolmogorov-Smirnov's statistics
#' @export
#'
#' @examples
sep.meas.along.time <- function(data1, data2, time.col, measure.col){
  timev <- sort(unique(data1[, get(time.col)]))
  if(!(identical(sort(unique(data2[, get(time.col)])), timev))) stop("Time vectors must be identical between the two data")
  out <- separability.measures(data1[get(time.col)==timev[1], get(measure.col)], data2[get(time.col)==timev[1], get(measure.col)])
  for(t in timev[2:length(timev)]){
    out <- rbind(out, separability.measures(data1[get(time.col)==t, get(measure.col)], data2[get(time.col)==t, get(measure.col)]))
  }
  out <- cbind(timev, out)
  return(out)
}


# from:
# https://stats.stackexchange.com/a/78855/43610

separability.measures <- function ( Vector.1 , Vector.2 ) {
  # convert vectors to matrices in case they are not
  Matrix.1 <- as.matrix (Vector.1)
  Matrix.2 <- as.matrix (Vector.2)
  # define means
  mean.Matrix.1 <- mean ( Matrix.1 )
  mean.Matrix.2 <- mean ( Matrix.2 )
  # define difference of means
  mean.difference <- mean.Matrix.1 - mean.Matrix.2
  # define covariances for supplied matrices
  cv.Matrix.1 <- cov ( Matrix.1 )
  cv.Matrix.2 <- cov ( Matrix.2 )
  # define the halfsum of cv's as "p"
  p <- ( cv.Matrix.1 + cv.Matrix.2 ) / 2
  # --%<------------------------------------------------------------------------
  # calculate the Bhattacharryya index
  bh.distance <- 0.125 *t ( mean.difference ) * p^ ( -1 ) * mean.difference +
    0.5 * log (det ( p ) / sqrt (det ( cv.Matrix.1 ) * det ( cv.Matrix.2 )
    )
    )
  # --%<------------------------------------------------------------------------
  # calculate Jeffries-Matusita
  # following formula is bound between 0 and 2.0
  jm.distance <- 2 * ( 1 - exp ( -bh.distance ) )
  # also found in the bibliography:
  # jm.distance <- 1000 * sqrt (   2 * ( 1 - exp ( -bh.distance ) )   )
  # the latter formula is bound between 0 and 1414.0
  # --%<------------------------------------------------------------------------
  # calculate the divergence
  # trace (is the sum of the diagonal elements) of a square matrix
  trace.of.matrix <- function ( SquareMatrix ) {
    sum ( diag ( SquareMatrix ) ) }
  # term 1
  divergence.term.1 <- 1/2 * trace.of.matrix (( cv.Matrix.1 - cv.Matrix.2 ) * 
                                                ( cv.Matrix.2^ (-1) - cv.Matrix.1^ (-1) )
  )
  # term 2
  divergence.term.2 <- 1/2 * trace.of.matrix (( cv.Matrix.1^ (-1) + cv.Matrix.2^ (-1) ) *
                                                ( mean.Matrix.1 - mean.Matrix.2 ) *
                                                t ( mean.Matrix.1 - mean.Matrix.2 )
  )
  # divergence
  divergence <- divergence.term.1 + divergence.term.2
  # --%<------------------------------------------------------------------------
  # and the transformed divergence
  transformed.divergence  <- 2 * ( 1 - exp ( - ( divergence / 8 ) ) )
  
  # KS stat
  ks <- ks.test(Vector.1 , Vector.2)$statistic
  
  indices <- data.frame(
    jm=jm.distance,bh=bh.distance,div=divergence,tdiv=transformed.divergence, ks=ks)
  return(indices)
}


# -------------------------


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