perform.all.analysis <- function(data, measure.var, condition.var, time.var, color.pca.var, groups.pca = NULL, PC.to.plot = c(1,2), permutation = TRUE, nperm = 100, bootstraps = "both", nboot = 100){
  
}

# --------------

#' Report PCA
#' 
#' A function for running PCA and optionally get plots out of it.
#'
#' @param data A data table which contains time series in long format or a numeric matrix with time series per row.
#' @param what A character vector describing how the data should be handled before running pca, as well as which representations
#' should be plotted. Valid isntructions are: c("cast.and.pca", "nocast.and.pca", "pca.only", "plot.pca", "plot.variance", "plot.extremes")
#' If data are a data.table in long format, use "cast.and.pca", if data are in a matrix use "nocast.and.pca".
#' The "plot.xxx" settings call biplot.PCA and visualize.extremes.PCA.
#' @param na.fill Value to replace NA after casting data table from wide to long.
#' @param measure.var Character. Column name of the measurement that defines the time series in time.
#' @param condition.var Character vector. Column names used for casting from long to wide. Should also contain the name of 
#' the column used for coloring the PCA plot if requested.
#' @param time.var Character. Column name of the time measure.
#' @param PC Numeric length 2. Which Principal Components to plot?
#' @param label.color.pca Character or Vector used for coloring PCA. If data is long data.table (i.e. 'what' is set to "cast.and.pca")
#'  should contain the name of the column used for coloring; note that this column should also be provided in 'condition.var'.
#'  If data is a matrix (i.e. 'what' is set to "nocast.and.pca') vector of length equal to number of rows in data.
#' @param center.pca Should variables be centered before running PCA. Default is TRUE.
#' @param scale.pca Should variables be scaled before running PCA. Default is TRUE, but is susceptible to be changed.
#' @param n.extremes Numeric. How many extremes trajectories to plot.
#' @param ... additional parameters for biplot.PCA and visualize.extremes. For example var.axes=F to remove variable arrows
#' or tails = "positive" to plot only extremes trajectories on positive tail of the PCs.
#'
#' @return If 'what' is set to "pca.only", returns PCA object. Otherwise plot PCA result.
#' @export
#'
#' @examples
#' library(data.table)
#' # Create some dummy data, imagine 10 time series under three conditions A, B or C in a long data.table
#' number.measure <- 101
#' mydata <- data.table(Condition = rep(LETTERS[1:3], each = 10*number.measure),
#'  Label = rep(1:10, each = number.measure),
#'  Time=rep(seq(0,100), 30))
#' 
#' # A: oscillate around 1 for 10 time units, then shift to oscillation around 1.3
#' # B: oscillate around 1 for 10 time units, then peak to 1.3 and gets back to 1
#' # C: oscillate around 1 all along trajectory
#' 
#' mydata[Condition=="A", Measure := c(rnorm(10, 1, 0.05), rnorm(91, 1.3, 0.05)), by = "Label"]
#' mydata[Condition=="B", Measure := c(rnorm(10, 1, 0.05), rnorm(91, 1.3, 0.05) - seq(0, 0.3, length.out = 91)), by = "Label"]
#' mydata[Condition=="C", Measure := rnorm(101, 1, 0.05), by = "Label"]
#' ggplot(mydata, aes(x=Time, y=Measure)) + geom_line(aes(group=Label), alpha = 0.3) + facet_wrap("Condition") +
#'  stat_summary(fun.y = mean, geom = "line", col = "red", size = 1.25)
#' 
#' report.PCA(mydata, what = c("cast.and.pca", "plot.variance", "plot.pca", "plot.extremes"), measure.var="Measure",
#' condition.var=c("Condition", "Label"), time.var="Time", label.color.pca = "Condition", PC=c(1,2), n.extremes = 5)
#' 
report.PCA <- function(data, what, measure.var = NULL, condition.var = NULL, time.var = NULL, PC = c(1,2), na.fill = NULL, label.color.pca = NULL, center.pca = T, scale.pca = T, n.extremes = NULL,...){
  require(ggbiplot)
  # Argument check
  arg.what <- c("cast.and.pca", "nocast.and.pca", "pca.only", "plot.pca", "plot.variance", "plot.extremes")
  if(!all(what %in% arg.what)) stop(paste('what arguments must be one of', arg.what))
  if(!any(c("cast.and.pca", "nocast.and.pca") %in% what)) stop("One of c('cast.and.pca', 'nocast.and.pca') must be provided. Use cast.and.pca for data.table in long format. Use nocast.and.pca for numeric matrix ready for stats::prcomp")
  if(all(c("cast.and.pca", "nocast.and.pca") %in% what)) stop("Only one of c('cast.and.pca', 'nocast.and.pca') must be provided. Use cast.and.pca for data.table in long format. Use nocast.and.pca for numeric matrix ready for stats::prcomp")
  if(("data.frame" %in% class(data) & "nocast.and.pca" %in% what) | ("numeric" %in% class(data) & "cast.and.pca" %in% what)) warning("Use cast.and.pca for data.table in long format. Use nocast.and.pca for numeric matrix ready for stats::prcomp")
  
  # If data is a long data.table
  if("cast.and.pca" %in% what){
    cast <- cast.and.fill(data = data, condition.var = condition.var, time.var = time.var, na.fill = na.fill, measure.var = measure.var)
    pca <- prcomp(cast$mat2, center = center.pca, scale = scale.pca)
  }
  # If data is already a wide numeric matrix
  else if("nocast.and.pca" %in% what){
    pca <- prcomp(data, center = center.pca, scale. = scale.pca)
  }
  
  # Stop and return PCA object
  if("pca.only" %in% what) return(pca)
  
  # Plot explained variance
  if("plot.variance" %in% what) plot(pca)
  
  # Biplot PCA
  if("plot.pca" %in% what){
    if(is.null(label.color.pca)){
      stop("'label.color.pca' must be a name of a variable in 'data' and provided in 'condition.var' if data is a data.table. It must be a vector witht the groups if data is a matrix.")
    }
    
    if("nocast.and.pca" %in% what){
      plot(biplot.PCA(pca, label.color.pca, PC = PC, ...))
    }
    
    else if("cast.and.pca" %in% what & label.color.pca %in% colnames(cast$mat)){
      label.color.pca <- unlist(cast$mat[, label.color.pca])
      plot(biplot.PCA(pca, label.color.pca, PC = PC, ...))
    } else {
      stop("'label.color.pca' must be a name of a variable in 'data' and provided in 'condition.var' if data is a data.table. It must be a vector witht the groups if data is a matrix.")
    }
  }
  
  # Extremes plotting
  if("plot.extremes" %in% what){
    if(is.null(n.extremes)){
      n.extremes <- 3
      warning("n.extremes was not provided and automatically set to 3.")
    }
    if("cast.and.pca" %in% what) visualize.extremes.PCA(pca, cast$mat2, PC = PC, n = n.extremes, tails = "both", interact=F)
    else if("nocast.and.pca" %in% what) visualize.extremes.PCA(pca, data, PC = PC, n = n.extremes, tails = "both", interact=F)
  }
}


#' cast.and.fill
#'
#' Cast a data table from long to wide and fill missing values. Output is suitable for PCA.
#' @param data A data table in long format.
#' @param condition.var Character vector, names of variables used for rows in casting, define the conditions of the experiment.
#' @param time.var Character, name of Time variable used for columns in casting
#' @param measure.var Character, name of measurement variable, on which PCA is to be performed.
#' @param na.fill Numeric, values to replace NAs
#'
#' @return Two matrices. "mat" is a wide character matrix which contains the casted "condition.var" as first columns followed by casted measurements.
#' "mat2" is a numeric wide matrix, which is essentially the same with conditions trimmed.
#' @export
#'
#' @examples
cast.and.fill <- function(data, condition.var, time.var, measure.var, na.fill){
  formula <- as.formula(paste0(paste(condition.var, collapse = " + "), " ~ ", time.var))
  mat <- as.matrix(dcast(data, formula, value.var = measure.var))
  mat2 <- mat[,-(1:length(condition.var))]
  class(mat2) <- "numeric"
  mat2[which(is.na(mat2))] <- na.fill
  return(list(mat = mat, mat2 = mat2))
}


biplot.PCA <- function(pca.obj, group.vec, ..., PC = c(1,2), obs.scale = 1, var.scale = 1, ellipse = T, circle = T, var.axes = T){
  require(ggbiplot)
  g <- ggbiplot(pca.obj, ..., obs.scale = obs.scale, var.scale = var.scale, 
                groups = group.vec, ellipse = ellipse, 
                circle = circle, var.axes = var.axes, choices = PC)
  g <- g + scale_color_discrete(name = '')
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
  return(g)
}


visualize.extremes.PCA <- function(pca.object, matrix.data, PC, n, tails = "both", interact=T){
  # Visualize the n extreme individuals of the component PC. Tails indicates if the extremes are to be picked from left, right or both tails of PC. Matrix.data is the matrix in wide format that was used to run the pca analysis with stats::prcomp.
  
  ord <- order(pca.object$x[,PC])
  ordd <- rev(ord)
  if(tails=="both"){
    par(mfrow=c(1,2))
    for(i in 1:n){
      mini <-  min(matrix.data[ord[i],], matrix.data[ordd[i],])
      maxi <- max(matrix.data[ord[i],], matrix.data[ordd[i],])
      plot(matrix.data[ord[i],], type = "l", ylim = c(mini, maxi), xlab = "Time", ylab = "Trajectory", main = paste0("Negative Tail - Coord PC", PC, ": ", round(pca.object$x[ord[i], PC], 4)))
      plot(matrix.data[ordd[i],], type = "l", ylim = c(mini, maxi), xlab = "Time", ylab = "Trajectory", main = paste0("Positive Tail - Coord PC", PC, ": ", round(pca.object$x[ordd[i], PC], 4)))
      if(interact) readline(prompt="Press [enter] to see next plots")
    }
  }
  
  else if(tails=="positive"){
    for(i in 1:n){
      plot(matrix.data[ordd[i],], type = "l", xlab = "Time", ylab = "Trajectory", main = paste0("Positive Tail - Coord PC", PC, ": ", round(pca.object$x[ordd[i], PC], 4)))
      if(interact) readline(prompt="Press [enter] to see next plots")
    }
  }
  
  else if(tails=="positive"){
    for(i in 1:n){
      plot(matrix.data[ord[i],], type = "l", xlab = "Time", ylab = "Trajectory", main = paste0("Negative Tail - Coord PC", PC, ": ", round(pca.object$x[ord[i], PC], 4)))
      if(interact) readline(prompt="Press [enter] to see next plots")
    }
  }
}


# --------------

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