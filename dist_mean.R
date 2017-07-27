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


#####

Cora <- fread("C:/Users/pixel/Dropbox/Marc-Antoine/data/set1-Coralie/tCoursesSelected.csv")
Cora[, Ratio := objCyto_Intensity_MeanIntensity_imErkCorrOrig / objNuc_Intensity_MeanIntensity_imErkCorrOrig]
setkey(Cora, Image_Metadata_Site)
#means <- Cora[, .(mean.traj = mean(Ratio)), by = .(Image_Metadata_Site, RealTime),]


euclids <- dist_mean(Cora, "Image_Metadata_Site", "RealTime", "Ratio", "objNuc_TrackObjects_Label", return.mean = F)
#setkey(euclids, Image_Metadata_Site)
euclids[, .(Mean = mean(euclid_to_mean), Variance = var(euclid_to_mean), Min = min(euclid_to_mean), Max = max(euclid_to_mean)), Image_Metadata_Site]


p <- ggplot(euclids$euclid, aes(x=Image_Metadata_Site, y=euclid_to_mean)) + geom_boxplot(aes(group=Image_Metadata_Site)) + scale_y_continuous(limits = c(0,4.5))
p
