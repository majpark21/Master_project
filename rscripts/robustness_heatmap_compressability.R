#' y_gradient_heatmap
#'
#' Estimate dispersion of a set of time series by reordering them with
#' hierarchical clustering and compute time-wise difference between neighboring
#' trajectories.
#' @param data Numeric matrix in wide format (i.e. columns are time-points and rows are individuals)
#' @param distance Character. Type of distance used for hierarchical clustering. See hclust for possible values.
#' @param linkage Character. Type of linkage used for hierarchical clustering. See hclust for possible values.
#' @param n Numeric. How many resampling procedures should be performed?
#' @param size Numeric. Size of samples.
#'
#' @return A data table with various statistics of the y-derivative matrix.
#' @export
#'
#' @details The function estimates the complexity of a set of time series across
#'   trajectories. It is a measure of dispersion between trajectories.
#'
#'   Idea is similar to CV along time: in a heatmap similar trajectories are
#'   grouped together resulting in regions of low-complexity. High complexity is
#'   at the interface of clusters or within ill-defined clusters. Here we simply
#'   compute the difference between neighboring rows of the heatmap, without
#'   accounting for difference in rows (we care about difference between
#'   trajectories, not within trajectories). In a heatmap image representation,
#'   this is equivalent to linear filetring to get y-gradient. This is not using
#'   heatmap image explicetly, but does the following:
#'    1) Perform clustering on
#'   the trajectories, this clustering could be then used to reorder rows for
#'   heatmap representation
#'    2) Instead of plotting a heatmap here, simply
#'   reorder the data according to clustering order.
#'    3) Compute a “y-derivative”
#'   with difference of 1 (i.e. with closest neighbor)
#'    4) The result is in the
#'   form of a partial y-derivative matrix
#'    5) (clip small deviations)
#'    6) Extract
#'   useful stats out of that matrix
#' 
#' The procedure is repeated on subsamples of trajectories to get an estimation of the error.
#'
#' @examples
y_gradient_heatmap <- function(data, n, size, thresh.clip = NULL, distance = "euclidean", linkage = "complete"){
  if(!is.null(thresh.clip)){
    result <- data.table(matrix(0, nrow = 1, ncol = 8))
    names(result) <- c("sum", "mean", "median", "min", "max", "sum.clipped", "mean.clipped", "median.clipped")
  } else {
    result <- data.table(matrix(0, nrow = 1, ncol = 5))
    names(result) <- c("sum", "mean", "median", "min", "max")
  }
  flist <- list(sum, mean, median, min, max)
  flist2 <- list(sum, mean, median)
  for(i in 1:n){
    # 0) Subsample data
    temp <- sample(x = 1:nrow(data), size = size, replace = F)
    tempdat <- data[temp,]
    
    # 1) Hclust to get row orders
    clust <- hclust(d = dist(tempdat, method = distance), method = linkage)
    reord.data <- as.matrix(tempdat[clust$order,])
    
    # 2) y-partial derivative
    # Matrix shifted by one row, to perform differential by matrix computation
    # First row is added as last row to keep dimensions (circularization)
    diff.mat <- rbind(reord.data[-1,], reord.data[1,])
    diff.mat <- abs(diff.mat - reord.data)
    
    # 3) Extract statistics
    temp <- sapply(flist, function(f) f(diff.mat))
    # Clip diff smaller than 2%
    if(!is.null(thresh.clip)){
      diff.mat <- diff.mat[!diff.mat <= thresh.clip]
      temp2 <- sapply(flist2, function(f) f(diff.mat))
      newrow <- matrix(c(temp, temp2), nrow=1, ncol=8)
    } else {
      newrow <- matrix(c(temp), nrow=1, ncol=5)
    }
    
    result <- rbind(result, newrow, use.names=F)
  }
  # Remove initialization row and set types
  result <- result[-1]
  for(col in names(result)) set(result, j = col, value = as.numeric(result[[col]]))
  return(result)
}

result <- data.table(matrix(0, nrow = 1, ncol = 7))
names(result) <- c("sum", "mean", "median", "min", "max", "Condition", "Replicate")
for(lev in unique(Cora$Condition)){
  temp <- y_gradient_heatmap(Cora_wide[Condition==lev, -c(1,2)], n = 500, size = 25)
  temp$Condition <- lev
  temp$Replicate <- 1:nrow(temp)
  result <- rbind(result, temp)
}
result <- result[-1]

result.melt <- melt(result, id.vars = c("Condition", "Replicate"))
result.melt[, variable := factor(variable, levels = c("sum", "mean", "median", "min", "max"))]
result.melt[, freq.pulse := str_extract(Condition, "P[0-9]+(\\.[0-9]+)?")]
result.melt[, int.pulse := str_extract(Condition, "I[0-9]+")]
cond.good.order <- c("P1-I10", "P1-I25", "P1-I50", "P1-I100", "P1-I100-UO",
                     "P5-I10", "P5-I25", "P5-I50", "P5-I100", "P5-I100-UO",
                     "P10-I10", "P10-I25", "P10-I50", "P10-I100", "P10-I100-UO",
                     "P20-I10", "P20-I25", "P20-I50", "P20-I100", "P20-I100-UO")
result.melt[, Condition := factor(Condition, levels = cond.good.order)]
ggplot(result.melt, aes(x=Condition, y = value)) +
  geom_boxplot(aes(fill=freq.pulse)) +
  facet_wrap("variable", scales = "free", ncol = 1)
