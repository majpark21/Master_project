# 1) Get mean trajectory for a condition
# 2) Euclidean distance between each trajectory and the mean one
# 3) Distance distribution
library(data.table)
library(ggplot2)

dist_mean <- function(data, condition, tcol, measure, label, return.mean = F){
  # data: a data table containing trajectories and conditions
  # condition: column name that fully define an experimental condition, thus a set of trajectories analyzed together
  # tcol: time column name
  # measure: column name with quantity of interest
  # label: column name with label of individual objects in each condition (cell label)
  
  require(data.table)
  # 1) Get mean trajectory for each condition
  means <- data[, .(Average = mean(get(measure))), by = c(condition, tcol)]
  dat <- merge(data[, c(condition, label, tcol, measure), with = F], means, by = c(condition, tcol))
  
  # 2) Euclidean distance between each trajectory and the mean one
  distances <- dat[, .(sq_dist = (get(measure) - Average)^2) , by = c(condition, label)]
  euclid <- distances[, .(euclid_to_mean = sqrt(sum(sq_dist))), by = c(condition, label)]
  
  if(return.mean){
    return(list(means = means, euclid = euclid))
  } else {
    return(euclid)
  }
}
