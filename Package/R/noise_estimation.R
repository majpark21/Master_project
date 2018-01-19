#' noise_ydiff
#'
#' Provide an estimator of noise level within a group of time-series.
#' First step: compute hierarchical clustering and arrange data
#' according to the clustering order.
#' Second step: compute differences between neighboring rows.
#' Differences can then be raised to a power to accentuate gap
#' between large and small differences. The matrix of differences is then
#' summarized by a list of functions.
#'
#' @param x A numerical matrix. Individuals in rows, measurements in columns.
#' @param summary.fun A vector of character, give names of the functions
#' used to summarize the differences. Functions must take a numerical
#' matrix as first input. Functions are applied to the absolute values
#' of differences. Use "colMeans" to get mean difference at each time point.
#' @param power.diff numerical, to what power should differences be raised.
#' @param hcl.distance character, the method used to compute distances between
#' rows for hierarchical clustering. One of "euclidean", "maximum", "manhattan",
#' "canberra", "binary" or "minkowski". See ?stats::dist.
#' @param hcl.linkage character, the agglomeration method used in hierarchical
#' clustering. One of "ward.D", "ward.D2", "single", "complete",
#' "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or
#' "centroid" (= UPGMC). See ?stats::hclust.
#' @param plot.diff logical, whether to plot the histogram of differences.
#'
#' @return A list where each entry is the application of one of summary.fun
#' to the difference matrix.
#' @export
#'
#' @examples
#' set.seed(7)
#' x <- rnorm(n = 50)
#' x <- matrix(x, nrow=5, ncol = 10)
#' noise_ydiff(x, summary.fun= c("mean", "median", "sd", "colMeans"))
#'
noise_ydiff <- function(x, summary.fun = c("mean", "median", "sd"),
                        power.diff = 1,
                        hcl.distance = "euclidean", hcl.linkage = "complete",
                        plot.diff = T){
  # 0) Build list of summary functions
  funs <- lapply(summary.fun, function(f) eval(parse(text = f)))
  names(funs) <- summary.fun

  # 1) Get Hierarchical clustering order
  x.dist <- dist(x, method = hcl.distance)
  hcl.order <- hclust(x.dist, method = hcl.linkage)$order
  # Reorder x accordingly
  x <- x[hcl.order, ]

  # 2) Get difference between neighboring rows
  # Replicate x and shift it by one row
  x.2 <- rbind(x[-1, ], x[1, ])
  x.diff <- x - x.2

  # 3) Transform diff
  x.diff <- abs(x.diff)
  x.diff <- x.diff ^ power.diff
  if(plot.diff){
    hist(x.diff, freq = F)
    lines(density(x.diff), lwd = 2, col = "red")
  }

  # 4) Report statistics
  out <- lapply(funs, function(f) f(x.diff))
  return(out)
}
